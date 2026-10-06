"""Literature screening (SCREEN_DOCUMENTS, DEC-22): document pool → rules → LLM triage.

Layer 1 (no LLM): merge Europe PMC query results, LitVar2 variant-linked papers and
PubTator3 gene+disease papers into one pool, deduplicated by PMID, with the routes
that found each paper and the LitVar-linked variant names it contains.

Layer 2: deterministic rules record article type and access (open access, PMC,
supplementary files); an LLM reads ONLY title + abstract and answers fixed
questions, quoting the sentence its answer rests on. The host verifies the quote
is in the abstract and assigns priority by a declared rule — the model never sets
priority directly. Output: documents.csv for human review before deep reading.

Screening decides what to read, not what is true: no screening output becomes
evidence or changes any allele's rank.
"""

from __future__ import annotations

import csv
import io
import json
import random
import re
import threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from enum import StrEnum
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, Field

from tessera.llm import Budget, LLMClient, LLMError, ModelOutputStore, load_profile
from tessera.sources import literature
from tessera.schemas.base import sha256_hex
from tessera.sources.http import Fetcher, SourceError

POLICY_ID = "screening_policy_v0.1"
PROMPT_VERSION = "screen-abstract-v2"

YNU = Literal["yes", "no", "unclear"]


class VariantReporting(StrEnum):
    SPECIFIC_ALLELES_NAMED = "specific_alleles_named"
    VARIANTS_REPORTED_NOT_NAMED = "variants_reported_not_named"
    NONE = "none"
    UNCLEAR = "unclear"


class VariantClass(StrEnum):
    RARE_CODING_OR_SPLICE = "rare_coding_or_splice"
    COMMON_SNP_ASSOCIATION = "common_snp_association"
    SOMATIC_CANCER = "somatic_cancer"
    ENGINEERED_MUTATION = "engineered_mutation"
    NONE = "none"
    UNCLEAR = "unclear"


class PaperScreen(BaseModel):
    pmid: str
    about_human_gene: YNU = Field(description="Is the human gene itself a subject (not only a homolog/complex/other gene)?")
    variant_reporting: VariantReporting
    variant_class: VariantClass
    allele_functional_experiment: YNU
    systems: list[Literal["patients", "human_cells", "animal_model", "nonhuman_cells", "in_vitro", "computational",
                          "none", "unclear"]]
    case_or_cohort_study: YNU
    is_review: YNU
    evidence_quote: str = Field(description="Verbatim sentence fragment from the title/abstract supporting the answers")
    rationale: str


class ScreenBatch(BaseModel):
    papers: list[PaperScreen]


SYSTEM = """You triage scientific papers for a variant-screening project. You see only a title and abstract.
Answer strictly from that text; when the text does not say, answer "unclear".
Return only a JSON object: {"papers": [ ...one object per input paper, same pmid... ]}.
Each object has exactly these keys:
- pmid: copy from input
- about_human_gene: "yes" if the human gene named in the task is itself studied (not only a yeast/fly homolog, a protein complex in general, or a different gene with a similar name); "no"; or "unclear"
- variant_reporting: "specific_alleles_named" if the abstract names particular variants of that gene (e.g. c.4582-2delAG, p.Arg990Ter); "variants_reported_not_named" if it reports variants of that gene in people or models without naming them (e.g. "we identified rare loss-of-function variants in GENE"); "none"; or "unclear"
- variant_class: what kind of variants of that gene the paper is about: "rare_coding_or_splice" (rare germline protein-truncating, missense, splice or similar variants in patients or families); "common_snp_association" (common SNPs/rsIDs in GWAS or association studies); "somatic_cancer" (tumour mutations); "engineered_mutation" (designed mutants such as catalytic-dead or phospho-site substitutions); "none"; or "unclear"
- allele_functional_experiment: "yes" if a specific variant (not just gene knockout/knockdown) was experimentally tested; "no"; or "unclear"
- systems: list from ["patients","human_cells","animal_model","nonhuman_cells","in_vitro","computational","none","unclear"]
- case_or_cohort_study: "yes" if it analyses patients/cases/cohorts genetically; "no"; or "unclear"
- is_review: "yes" if a review, commentary or editorial; "no"; or "unclear"
- evidence_quote: a short fragment copied EXACTLY from the title or abstract that best supports your answers
- rationale: one short sentence"""


def _norm(s: str) -> str:
    s = re.sub(r"<[^>]+>", " ", s or "")
    return re.sub(r"\s+", " ", s).strip().lower()


def build_pool(req: dict, f: Fetcher, gene_symbol: str, log=print) -> tuple[dict[str, dict], dict]:
    lit = req["literature"]
    docs: dict[str, dict] = {}
    routes: dict[str, set] = defaultdict(set)
    litvar_names: dict[str, set] = defaultdict(set)
    audit = {}

    # Europe PMC query template
    q = (f'"{gene_symbol}" AND (' + " OR ".join(f'"{t}"' for t in lit["disease_terms"]) + ") AND ("
         + " OR ".join(f'"{t}"' for t in lit["variant_terms"]) + ")")
    recs, _ = literature.europepmc_search(f, q)
    for r in recs:
        if r.get("pmid"):
            docs[r["pmid"]] = r
            routes[r["pmid"]].add("europepmc")
    audit["europepmc"] = {"query": q, "records": len(recs), "with_pmid": sum(1 for r in recs if r.get("pmid"))}
    log(f"Europe PMC: {len(recs)} records")

    # LitVar2 variant-linked papers
    rows, _ = literature.litvar_gene_variants(f, gene_symbol)
    lv_failures: list[str] = []
    for row in rows:
        name = row.get("rsid") or literature.litvar_name(row)
        if name is None:
            continue  # gene-level aggregate row ("All <gene> variants")
        try:
            pubs = literature.litvar_publications(f, row["_id"])
        except SourceError as e:
            lv_failures.append(f"{row['_id']}: {e.code}")
            continue
        for pmid in pubs:
            routes[pmid].add("litvar")
            litvar_names[pmid].add(name)
    audit["litvar"] = {"variant_entries": len(rows), "papers": len(litvar_names), "failures": lv_failures}
    log(f"LitVar2: {len(rows)} variant entries → {len(litvar_names)} papers")

    # PubTator3 gene + disease entities
    pt_total = {}
    for dis in lit["pubtator_diseases"]:
        text = f"@GENE_{gene_symbol} AND @DISEASE_{dis}"
        pmids, _, total = literature.pubtator_search(f, text)
        pt_total[text] = {"count": total, "retrieved": len(pmids)}
        for p in pmids:
            routes[p].add("pubtator")
    audit["pubtator"] = pt_total
    log(f"PubTator3: {sum(v['retrieved'] for v in pt_total.values())} hits over {len(pt_total)} queries")

    missing = sorted(p for p in routes if p not in docs)
    extra, _ = literature.europepmc_by_pmids(f, missing)
    for r in extra:
        if r.get("pmid") in routes:
            docs[r["pmid"]] = r
    unresolved = sorted(p for p in routes if p not in docs)
    audit["metadata_unresolved_pmids"] = unresolved
    for p in unresolved:
        docs[p] = {"pmid": p, "title": "", "abstractText": ""}

    pool = {}
    for pmid, r in docs.items():
        types = [t for t in (r.get("pubTypeList") or {}).get("pubType", [])]
        pool[pmid] = {
            "pmid": pmid,
            "pmcid": r.get("pmcid"),
            "doi": r.get("doi"),
            "year": r.get("pubYear"),
            "journal": ((r.get("journalInfo") or {}).get("journal") or {}).get("title"),
            "title": re.sub(r"<[^>]+>", "", r.get("title") or "").strip(),
            "abstract": re.sub(r"<[^>]+>", " ", r.get("abstractText") or "").strip(),
            "pub_types": types,
            "routes": sorted(routes[pmid]),
            "litvar_variants": sorted(litvar_names.get(pmid, [])),
            "is_open_access": r.get("isOpenAccess") == "Y",
            "in_pmc": r.get("inPMC") == "Y",
            "has_suppl": r.get("hasSuppl") == "Y",
        }
    audit["pool_size"] = len(pool)
    return pool, audit


def rule_fields(d: dict) -> dict:
    types = " ".join(d["pub_types"]).lower()
    review = any(k in types for k in ("review", "editorial", "comment"))
    if d["is_open_access"] and d["in_pmc"]:
        access = "full_text_open"
    elif d["abstract"]:
        access = "abstract_only"
    else:
        access = "metadata_only"
    suppl = ("supplement_available" if d["has_suppl"] and access == "full_text_open"
             else "supplement_unavailable" if d["has_suppl"] else "no_supplement_listed")
    return {"rule_review": review, "access_status": access, "supplement_status": suppl}


def priority(d: dict, s: PaperScreen | None) -> tuple[str, str]:
    """Declared screening rule (screening_policy_v0.1). Returns (priority, reason)."""
    # LitVar-linked HGVS names (not bare rsIDs, which are mostly GWAS SNPs) mean the
    # paper names specific alleles, independently of the LLM.
    named = any(not v.startswith("rs") for v in d["litvar_variants"])
    if s is None:
        if d.get("body_alleles"):
            return "P2", "full text names variant(s) near the gene"
        if named and d["access_status"] == "full_text_open":
            # LitVar gene assignment is not trusted on its own: the full text was
            # scanned and shows no variant near the gene.
            return "P4", "LitVar link not confirmed in full text"
        if named:
            return "P3", "LitVar links named variants; full text unavailable to confirm"
        if d.get("body_gene_mentions", 0) >= 3:
            return "P3", "gene discussed in full text"
        if d.get("abstract_mentions_gene") and d.get("screen_error"):
            return "P3", "gene in abstract; LLM screen not completed"
        if d["access_status"] != "full_text_open" and not d.get("abstract_mentions_gene"):
            return "P4", "gene only outside the abstract and full text unavailable: not assessable"
        return "P4", "incidental mention"
    if s.about_human_gene == "no":
        return "P5", "not about the human gene"
    if d["rule_review"] or s.is_review == "yes":
        return "P4", "review/commentary: pointer only"
    if s.variant_class in (VariantClass.COMMON_SNP_ASSOCIATION, VariantClass.SOMATIC_CANCER):
        return "P4", f"{s.variant_class.value}: not germline rare-allele evidence"
    if s.variant_class is VariantClass.ENGINEERED_MUTATION:
        return "P3", "engineered mutants: mechanism, not patient alleles"
    rare = s.variant_class is VariantClass.RARE_CODING_OR_SPLICE
    allele = s.variant_reporting in (VariantReporting.SPECIFIC_ALLELES_NAMED,
                                     VariantReporting.VARIANTS_REPORTED_NOT_NAMED) or named
    if s.allele_functional_experiment == "yes" and rare:
        return "P1", "allele-level functional experiment on patient variants"
    if allele and rare and s.case_or_cohort_study == "yes":
        return "P1", "rare variants in patients/cohorts"
    if allele and rare:
        return "P2", "rare variants reported"
    if s.allele_functional_experiment == "yes" or allele:
        return "P3", "variants reported, class unclear"
    if s.case_or_cohort_study == "yes":
        return "P2", "gene-level patient/cohort study (supplements may list alleles)"
    if s.about_human_gene == "yes":
        return "P3", "gene-level mechanism/model study"
    return "P4", "unclear relevance"


def screen(request_path: Path, out: Path, *, llm: str = "live", refresh_pool: bool = False, log=print) -> Path:
    req = yaml.safe_load(request_path.read_text())
    cfg = req["literature"]["screening"]
    root = Path(__file__).resolve().parents[2]
    out.mkdir(parents=True, exist_ok=True)
    f = Fetcher(out / "snapshots", min_interval={"litvar2": 0.35, "pubtator3": 0.35, "europepmc": 0.1},
                cache_dir=root / "cache" / "http")
    gene = req["gene"]

    f.seed_from_snapshots(out / "snapshots", ("fullTextXML/", "publications/", "pmids/", "search/gene/"))
    pool_file = out / "pool.jsonl"
    if pool_file.exists() and not refresh_pool:  # frozen pool from an earlier run: no re-query
        pool = {d["pmid"]: d for d in map(json.loads, pool_file.read_text().splitlines())}
        audit = json.loads((out / "pool_audit.json").read_text())
    else:
        previous = set(json.loads(l)["pmid"] for l in pool_file.read_text().splitlines()) if pool_file.exists() else set()
        pool, audit = build_pool(req, f, gene, log=log)
        audit["new_since_previous_pool"] = sorted(set(pool) - previous) if previous else None
        if previous:
            log(f"Incremental: {len(set(pool) - previous)} papers new since the previous pool")
        pool_file.write_text("".join(json.dumps(d, ensure_ascii=False) + "\n" for d in pool.values()))
        (out / "pool_audit.json").write_text(json.dumps(audit, sort_keys=True))
    aliases = [gene] + list(req.get("literature", {}).get("gene_aliases", []))
    alias_re = re.compile("|".join(r"\b" + re.escape(a).replace(r"\ ", r"\s?") + r"\b" for a in aliases), re.I)
    for d in pool.values():
        d.update(rule_fields(d))
        d["abstract_mentions_gene"] = bool(alias_re.search(d["title"] + " " + d["abstract"]))
    log(f"Pool: {len(pool)} unique papers; {sum(d['abstract_mentions_gene'] for d in pool.values())} name the gene in title/abstract")

    # ---- free full-text scan for papers that mention the gene only in the body
    from tessera.supplements import C_NOTATION, P_NOTATION
    scanned = 0
    cache_file = out / "body_scan.json"
    cache = json.loads(cache_file.read_text()) if cache_file.exists() else {}
    for d in pool.values():
        d["body_gene_mentions"], d["body_alleles"] = 0, []
        if d["abstract_mentions_gene"] or d["access_status"] != "full_text_open" or not d.get("pmcid"):
            continue
        if d["pmid"] in cache:
            c = cache[d["pmid"]]
            d["body_gene_mentions"], d["body_alleles"] = c["mentions"], c["alleles"]
            if c.get("failed"):
                d["access_status"] = "full_text_fetch_failed"
            continue
        xml = literature.europepmc_fulltext(f, d["pmcid"])
        if not xml:
            d["access_status"] = "full_text_fetch_failed"
            cache[d["pmid"]] = {"mentions": 0, "alleles": [], "failed": True}
            continue
        scanned += 1
        text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", xml))
        pos = [m.start() for m in alias_re.finditer(text)]
        d["body_gene_mentions"] = len(pos)
        found = set()
        for p in pos:
            win = text[max(0, p - 300) : p + 300]
            found.update(m.group(0) for m in C_NOTATION.finditer(win))
            found.update(m.group(0).replace("(", "").replace(")", "") for m in P_NOTATION.finditer(win))
        d["body_alleles"] = sorted(found)
        cache[d["pmid"]] = {"mentions": d["body_gene_mentions"], "alleles": d["body_alleles"]}
    cache_file.write_text(json.dumps(cache, sort_keys=True))
    log(f"Full-text scan: {scanned} open-access papers; "
        f"{sum(1 for d in pool.values() if d['body_alleles'])} name variants near the gene")

    # ---- LLM triage on title + abstract
    prof = load_profile(root / "config" / "llm.yaml", cfg.get("llm_profile"), env_file=root / ".env")
    if cfg.get("max_tokens"):
        from dataclasses import replace

        prof = replace(prof, max_tokens=int(cfg["max_tokens"]))
    budget = Budget(int(cfg["max_llm_requests"]), int(cfg.get("max_input_tokens", 1_500_000)),
                    int(cfg.get("max_output_tokens", 300_000)))
    # Model outputs live in a shared cache so any run can replay them; per-paper
    # results are cached separately so a changed pool only sends NEW papers.
    client = LLMClient(prof, ModelOutputStore(root / "cache" / "model_outputs"), budget,
                       mode="replay" if llm == "off" else llm)
    item_cache = root / "cache" / "screening"
    item_cache.mkdir(parents=True, exist_ok=True)

    def item_key(p: str) -> str:
        return sha256_hex(json.dumps({"v": PROMPT_VERSION, "model": prof.model, "gene": gene,
                                      "disease": req["disease"]["input"], "pmid": p,
                                      "title": pool[p]["title"], "abstract": pool[p]["abstract"]}, sort_keys=True))
    lock = threading.Lock()
    screens: dict[str, PaperScreen] = {}
    errors: dict[str, str] = {}
    stop = threading.Event()
    wanted = sorted(p for p, d in pool.items() if d["abstract_mentions_gene"] and (d["abstract"] or d["title"]))
    todo = []
    for p in wanted:
        cf = item_cache / f"{item_key(p)}.json"
        if cf.exists():
            screens[p] = PaperScreen.model_validate_json(cf.read_text())
        else:
            todo.append(p)
    log(f"Screening cache: {len(screens)} papers reused, {len(todo)} to send")
    bs = int(cfg.get("batch_size", 10))
    batches = [todo[i : i + bs] for i in range(0, len(todo), bs)]

    def run(batch):
        if stop.is_set():
            with lock:
                for p in batch:
                    errors[p] = "not_run_after_provider_stop"
            return
        items = [{"pmid": p, "title": pool[p]["title"], "abstract": pool[p]["abstract"]} for p in batch]
        user = (f"Gene: human {gene}. Disease context: {req['disease']['input']}.\n"
                f"Papers (JSON):\n{json.dumps(items, ensure_ascii=False)}")
        try:
            res = client.structured(task_id="screen:" + ",".join(batch), prompt_version=PROMPT_VERSION,
                                    system=SYSTEM, user=user, schema=ScreenBatch)
        except LLMError as e:
            if llm == "off" and e.code == "snapshot_missing":
                with lock:
                    for p in batch:
                        errors[p] = "llm_off"
                return
            if e.code in ("provider_quota_or_auth", "budget_exhausted", "missing_credentials"):
                stop.set()
            with lock:
                for p in batch:
                    errors[p] = f"{e.code}: {e.detail[:120]}"
            return
        with lock:
            if res.value is None:
                for p in batch:
                    errors[p] = res.error or "schema_failure"
                return
            got = {s.pmid: s for s in res.value.papers}
            for p in batch:
                if p in got:
                    screens[p] = got[p]
                    (item_cache / f"{item_key(p)}.json").write_text(got[p].model_dump_json())
                else:
                    errors[p] = "missing_from_model_reply"

    with ThreadPoolExecutor(max_workers=int(cfg.get("concurrency", 4))) as ex:
        list(ex.map(run, batches))
    log(f"LLM screening: {len(screens)} screened, {len(errors)} errors; requests {budget.requests}, "
        f"tokens in {budget.input_tokens} out {budget.output_tokens}")

    # ---- host checks + priority
    rng = random.Random(int(cfg.get("audit_seed", 42)))
    rows = []
    for pmid in sorted(pool, key=lambda p: (-(int(pool[p]["year"]) if (pool[p]["year"] or "").isdigit() else 0), p)):
        d = pool[pmid]
        s = screens.get(pmid)
        quote_ok = bool(s) and _norm(s.evidence_quote) in _norm(d["title"] + " " + d["abstract"]) and len(s.evidence_quote) > 8
        d["screen_error"] = pmid in errors
        pr, why = priority(d, s)
        flags = []
        if s and not quote_ok:
            flags.append("quote_not_found_in_abstract")
        if pmid in errors:
            flags.append(f"screen_error:{errors[pmid]}")
        deep = pr in ("P1", "P2") and d["access_status"] == "full_text_open"
        if pr in ("P1", "P2") and not deep:
            flags.append("needs_manual_access")
        rows.append({**{k: d[k] for k in ("pmid", "pmcid", "year", "journal", "title")},
                     "priority": pr, "priority_reason": why,
                     "deep_read": deep, "audit_sample": False,
                     "access_status": d["access_status"], "supplement_status": d["supplement_status"],
                     "routes": ";".join(d["routes"]), "litvar_variants": ";".join(d["litvar_variants"][:15]),
                     "gene_in_abstract": d["abstract_mentions_gene"], "body_gene_mentions": d["body_gene_mentions"],
                     "body_alleles": ";".join(d["body_alleles"][:15]),
                     "about_human_gene": s.about_human_gene if s else "",
                     "variant_reporting": s.variant_reporting.value if s else "",
                     "variant_class": s.variant_class.value if s else "",
                     "allele_functional_experiment": s.allele_functional_experiment if s else "",
                     "systems": ";".join(s.systems) if s else "",
                     "case_or_cohort_study": s.case_or_cohort_study if s else "",
                     "is_review": ("yes" if d["rule_review"] else (s.is_review if s else "")),
                     "evidence_quote": s.evidence_quote if s else "", "quote_verified": quote_ok,
                     "rationale": s.rationale if s else "", "flags": ";".join(flags),
                     "pub_types": ";".join(d["pub_types"])})
    low = [r for r in rows if r["priority"] in ("P4", "P5")]
    for r in rng.sample(low, k=min(len(low), max(1, int(len(low) * float(cfg.get("audit_fraction", 0.1)))))):
        r["audit_sample"] = True  # read anyway to estimate what screening throws away

    cols = list(rows[0].keys()) if rows else []
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    (out / "documents.csv").write_text(buf.getvalue(), encoding="utf-8")
    (out / "documents.jsonl").write_text("".join(json.dumps({**pool[r["pmid"]], **r}, ensure_ascii=False) + "\n"
                                                 for r in rows), encoding="utf-8")
    f.write_index()
    audit.update({"policy": POLICY_ID, "prompt_version": PROMPT_VERSION, "llm": prof.public,
                  "llm_requests": budget.requests, "llm_input_tokens": budget.input_tokens,
                  "llm_output_tokens": budget.output_tokens, "screen_errors": len(errors),
                  "http_requests": f.logical_requests, "http_attempts": f.attempts})
    (out / "screening_manifest.json").write_text(json.dumps(audit, indent=1, sort_keys=True, ensure_ascii=False))
    return out / "documents.csv"

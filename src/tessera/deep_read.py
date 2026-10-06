"""Layer 3: document-driven deep reading (DEC-22), passage-first to save LLM calls.

For each paper selected by screening (deep_read=True, open full text):
  3a  supplementary archive → deterministic notation scan (no LLM)
  3b  full text → split into units (paragraphs, table rows, captions) → keep only units
      that name the gene or a variant notation that co-occurs with it → pack several
      papers' passages into one LLM call → structured claims with verbatim quotes.

The host verifies every quote against the passage and converts it to a text_span
locator in a derived plain-text snapshot of the article, so the pipeline's ordinary
mechanical checks apply. Alleles stay as written; identity is resolved later by the
same variant_recoder rules as every other route. Nothing here sets a rank.
"""

from __future__ import annotations

import json
import re
import threading
from concurrent.futures import ThreadPoolExecutor
from enum import StrEnum
from pathlib import Path
from typing import Literal

from pydantic import BaseModel

from tessera.llm import Budget, LLMClient, LLMError, ModelOutputStore, load_profile
from tessera.schemas.base import sha256_hex
from tessera.sources import literature
from tessera.sources.http import Fetcher
from tessera.supplements import C_NOTATION, P_NOTATION, fetch_zip, owner_is_target, scan, units

PROMPT_VERSION = "deep-read-v1"
TEXT_PARSER = "jats-units-1"


class PerturbationKind(StrEnum):
    EXACT_ALLELE_EDIT = "exact_allele_edit"
    ALLELE_CONSTRUCT = "allele_construct"
    CARRIER_COMPARISON = "carrier_comparison"
    GENE_KNOCKOUT = "gene_knockout"
    KNOCKDOWN = "knockdown"
    OVEREXPRESSION = "overexpression"
    NOT_AN_EXPERIMENT = "not_an_experiment"
    OTHER = "other"


class ExtractedClaim(BaseModel):
    passage_id: str
    quote: str
    allele_as_written: str | None
    scope: Literal["variant", "gene"]
    kind: Literal["functional_experiment", "patient_observation", "other"]
    perturbation: PerturbationKind
    species: Literal["human", "mouse", "rat", "fly", "zebrafish", "yeast", "other", "unclear"]
    cell_or_model: str | None
    comparator: str | None
    endpoint: str | None
    outcome_type: Literal["altered", "no_detected_effect", "rescue", "mixed", "inconclusive", "not_applicable"]
    direction: Literal["decreased", "increased", "none", "mixed", "unclear", "not_applicable"]


class PaperClaims(BaseModel):
    pmid: str
    claims: list[ExtractedClaim]


class ClaimBatch(BaseModel):
    papers: list[PaperClaims]


SYSTEM = """You extract evidence about a human gene from selected passages of research papers.
Use ONLY the passages given. Do not infer beyond them. Return only JSON:
{"papers": [{"pmid": "...", "claims": [ ... ]}, ...]} with one entry per input paper (empty claims allowed).
Extract a claim when a passage reports, for that gene:
 (a) a specific variant observed in patients/families/cohorts, or
 (b) an experiment that tested a specific variant (edited cells, patient-derived cells, constructs), or
 (c) a gene-level perturbation (knockout, knockdown, heterozygous loss) with a measured outcome.
Include experiments that found no effect. Skip background statements and citations of other work.
Each claim has exactly these keys:
- passage_id: the id of the passage containing the quote
- quote: a short span copied EXACTLY (character for character) from that passage
- allele_as_written: the variant exactly as written in the passage (e.g. "c.4596_4597insG", "p.Arg990Ter"), or null for gene-level claims
- scope: "variant" or "gene"
- kind: "functional_experiment", "patient_observation" or "other"
- perturbation: one of exact_allele_edit (isogenic edit installing the variant), allele_construct (expression construct carrying it), carrier_comparison (patient-derived cells/carriers vs unrelated controls), gene_knockout, knockdown, overexpression, not_an_experiment, other
- species: human, mouse, rat, fly, zebrafish, yeast, other, unclear
- cell_or_model: e.g. "iPSC-derived excitatory neurons", or null
- comparator: e.g. "isogenic wild type", "unrelated healthy donors", or null
- endpoint: what was measured, or null
- outcome_type: altered, no_detected_effect, rescue, mixed, inconclusive, not_applicable
- direction: decreased, increased, none, mixed, unclear, not_applicable"""


def article_units(xml: str) -> list[str]:
    """Paragraphs, table rows and captions of a JATS article as whitespace-normalized text."""
    out, seen = [], set()
    for tag, body in re.findall(r"<(p|tr|caption|title)\b[^>]*>(.*?)</\1>", xml, re.S):
        if tag == "tr":
            body = re.sub(r"</t[dh]>", " | ", body)
        t = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", body)).strip()
        if len(t) > 3 and t not in seen:
            seen.add(t)
            out.append(t)
    return out


def select_passages(us: list[str], alias_re: re.Pattern, max_units: int, max_chars: int) -> list[int]:
    gene_units = [i for i, u in enumerate(us) if alias_re.search(u)]
    notes = {m.group(0) for i in gene_units for m in list(C_NOTATION.finditer(us[i])) + list(P_NOTATION.finditer(us[i]))}
    extra = [i for i, u in enumerate(us) if i not in gene_units and any(n in u for n in notes)]
    both = [i for i in gene_units if C_NOTATION.search(us[i]) or P_NOTATION.search(us[i])]
    ordered = both + [i for i in gene_units if i not in both] + extra
    picked, total = [], 0
    for i in ordered:
        if len(picked) >= max_units or total + len(us[i]) > max_chars:
            continue
        picked.append(i)
        total += len(us[i])
    return sorted(picked)


def deep_read(request_path: Path, lit_dir: Path, *, llm: str = "live", log=print) -> dict:
    import csv
    import yaml

    req = yaml.safe_load(request_path.read_text())
    cfg = req["literature"]["deep_read"]
    root = Path(__file__).resolve().parents[2]
    gene = req["gene"]
    aliases = [gene] + list(req["literature"].get("gene_aliases", []))
    alias_re = re.compile("|".join(r"\b" + re.escape(a).replace(r"\ ", r"\s?") + r"\b" for a in aliases), re.I)
    f = Fetcher(lit_dir / "snapshots", min_interval={"europepmc": 0.1}, cache_dir=root / "cache" / "http",
                timeout=float(cfg.get("download_timeout_s", 120)), max_retries=1)
    docs = [r for r in csv.DictReader(open(lit_dir / "documents.csv")) if r["deep_read"] == "True"]
    log(f"Deep read: {len(docs)} papers")

    text_snaps: dict[str, dict] = {}
    supplement_hits = []
    sup_status = {}
    for n_doc, d in enumerate(docs, start=1):
        pmid, pmcid = d["pmid"], d["pmcid"]
        log(f"  [{n_doc}/{len(docs)}] PMID {pmid} {pmcid} supplements={d['supplement_status']}")
        # 3a supplements
        if d["supplement_status"] == "supplement_available":
            z, info = fetch_zip(f, pmcid)
            if z is None:
                sup_status[pmid] = info
            else:
                hs = scan(pmid, units(z), aliases)
                sup_status[pmid] = f"scanned:{len(hs)}"
                supplement_hits += [{**h.__dict__, "snapshot_id": info} for h in hs]
        # 3b full text units → derived text snapshot
        xml = literature.europepmc_fulltext(f, pmcid)
        if not xml:
            continue
        us = article_units(xml)
        text = "\n\n".join(us)
        offsets, pos = [], 0
        for u in us:
            offsets.append(pos)
            pos += len(u) + 2
        raw = text.encode("utf-8")
        snap_id = "snap_" + sha256_hex(raw)
        fname = f"europepmc_text_{sha256_hex(raw)[:16]}.txt"
        (lit_dir / "snapshots" / fname).write_bytes(raw)
        f._meta.append({"file": fname, "snapshot_id": snap_id, "source": "europepmc_text",
                        "source_record_id": f"{pmcid}:{TEXT_PARSER}", "source_version": TEXT_PARSER,
                        "source_uri": f"https://europepmc.org/article/PMC/{pmcid}",
                        "retrieved_at": "derived", "media_type": "text/plain",
                        "raw_content_digest": sha256_hex(raw), "license_status": "open_access",
                        "is_synthetic": False})
        picked = select_passages(us, alias_re, int(cfg.get("max_units_per_paper", 15)),
                                 int(cfg.get("max_chars_per_paper", 9000)))
        text_snaps[pmid] = {"snapshot_id": snap_id, "text": text, "units": us, "offsets": offsets,
                            "passages": {f"p{i}": i for i in picked}, "pmcid": pmcid}
    f.write_index()
    log(f"Supplements: {sum(1 for v in sup_status.values() if v.startswith('scanned'))} scanned, "
        f"{len(supplement_hits)} notation hits; full texts: {len(text_snaps)}; "
        f"passages: {sum(len(v['passages']) for v in text_snaps.values())}")

    # ---- pack papers into calls by character budget
    prof = load_profile(root / "config" / "llm.yaml", cfg.get("llm_profile"), env_file=root / ".env")
    budget = Budget(int(cfg["max_llm_requests"]), int(cfg.get("max_input_tokens", 2_000_000)),
                    int(cfg.get("max_output_tokens", 500_000)))
    client = LLMClient(prof, ModelOutputStore(root / "cache" / "model_outputs"), budget,
                       mode="replay" if llm == "off" else llm)
    item_cache = root / "cache" / "deep_read"
    item_cache.mkdir(parents=True, exist_ok=True)

    def key(pmid):
        t = text_snaps[pmid]
        return sha256_hex(json.dumps({"v": PROMPT_VERSION, "model": prof.model, "gene": gene,
                                      "passages": [t["units"][i] for i in t["passages"].values()]}, sort_keys=True))

    results: dict[str, PaperClaims] = {}
    pending = []
    for pmid, t in sorted(text_snaps.items()):
        cf = item_cache / f"{key(pmid)}.json"
        if not t["passages"]:
            results[pmid] = PaperClaims(pmid=pmid, claims=[])
        elif cf.exists():
            results[pmid] = PaperClaims.model_validate_json(cf.read_text())
        else:
            pending.append(pmid)
    batches, cur, size = [], [], 0
    limit = int(cfg.get("max_chars_per_call", 45000))
    for pmid in pending:
        n = sum(len(text_snaps[pmid]["units"][i]) for i in text_snaps[pmid]["passages"].values())
        if cur and size + n > limit:
            batches.append(cur)
            cur, size = [], 0
        cur.append(pmid)
        size += n
    if cur:
        batches.append(cur)
    log(f"LLM: {len(results)} papers from cache, {len(pending)} to send in {len(batches)} calls")

    errors: dict[str, str] = {}
    lock, stop = threading.Lock(), threading.Event()

    def run(batch):
        if stop.is_set():
            with lock:
                errors.update({p: "not_run_after_provider_stop" for p in batch})
            return
        payload = [{"pmid": p, "passages": [{"id": pid, "text": text_snaps[p]["units"][i]}
                                            for pid, i in text_snaps[p]["passages"].items()]} for p in batch]
        user = f"Gene: human {gene}.\nPapers:\n{json.dumps(payload, ensure_ascii=False)}"
        try:
            res = client.structured(task_id="deep:" + ",".join(batch), prompt_version=PROMPT_VERSION,
                                    system=SYSTEM, user=user, schema=ClaimBatch)
        except LLMError as e:
            if e.code in ("provider_quota_or_auth", "budget_exhausted", "missing_credentials"):
                stop.set()
            with lock:
                errors.update({p: e.code for p in batch})
            return
        with lock:
            if res.value is None:
                errors.update({p: res.error or "schema_failure" for p in batch})
                return
            got = {pc.pmid: pc for pc in res.value.papers}
            for p in batch:
                if p in got:
                    results[p] = got[p]
                    (item_cache / f"{key(p)}.json").write_text(got[p].model_dump_json())
                else:
                    errors[p] = "missing_from_model_reply"

    with ThreadPoolExecutor(max_workers=int(cfg.get("concurrency", 2))) as ex:
        list(ex.map(run, batches))

    # ---- host verification → locators
    verified, rejected = [], []
    for pmid, pc in sorted(results.items()):
        t = text_snaps[pmid]
        for c in pc.claims:
            ui = t["passages"].get(c.passage_id)
            reason = None
            if ui is None:
                reason = "unknown passage_id"
            else:
                unit = t["units"][ui]
                q = re.sub(r"\s+", " ", c.quote).strip()
                at = unit.find(q)
                if len(q) < 10 or at < 0:
                    reason = "quote not found verbatim in passage"
                elif c.allele_as_written and c.allele_as_written not in unit:
                    reason = "allele not written in the cited passage"
                elif c.allele_as_written and not owner_is_target(unit, unit.find(c.allele_as_written), aliases):
                    reason = "nearest gene symbol before the allele is not the target gene"
            rec = {"pmid": pmid, "pmcid": t["pmcid"], **c.model_dump()}
            if reason:
                rejected.append({**rec, "rejection": reason})
                continue
            start = t["offsets"][ui] + at
            rec["locator"] = {"kind": "text_span", "snapshot_id": t["snapshot_id"],
                              "text_digest": sha256_hex(t["text"]), "parser_version": TEXT_PARSER,
                              "start": start, "end": start + len(q), "excerpt": q}
            verified.append(rec)

    out = {
        "papers": len(docs), "full_texts": len(text_snaps), "llm_requests": budget.requests,
        "claims_verified": len(verified), "claims_rejected": len(rejected), "errors": errors,
        "supplement_status": sup_status, "llm": prof.public, "prompt_version": PROMPT_VERSION,
    }
    (lit_dir / "deep_claims.jsonl").write_text("".join(json.dumps(v, ensure_ascii=False) + "\n" for v in verified))
    (lit_dir / "deep_claims_rejected.jsonl").write_text("".join(json.dumps(v, ensure_ascii=False) + "\n" for v in rejected))
    (lit_dir / "supplement_hits.jsonl").write_text("".join(json.dumps(h, ensure_ascii=False) + "\n" for h in supplement_hits))
    (lit_dir / "deep_read_manifest.json").write_text(json.dumps(out, indent=1, sort_keys=True))
    log(f"Claims: {len(verified)} verified, {len(rejected)} rejected by host checks; errors {len(errors)}; "
        f"LLM requests {budget.requests}")
    return out

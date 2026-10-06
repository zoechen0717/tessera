"""M1a: fetch structured sources for one gene and assemble a frozen-ready input.

    tessera fetch --request requests/<gene>.yaml --output <dir>

writes <dir>/input.yaml plus reference, variants, annotations, snapshots,
claims, links, decisions and mechanism in the same format the M0 pipeline
freezes. Everything downstream (freeze, derive, replay) is unchanged and
offline. No LLM is involved (DEC-10).

Candidate routes: ClinVar gene-level records and SCHEMA alleles observed in at
least one case. Observation evidence (O) comes from SCHEMA case/control allele
counts. Gene mechanism comes from a declared rule over ClinGen dosage and the
SCHEMA PTV burden result; constraint metrics are recorded but never decide it.
"""

from __future__ import annotations

import json
from pathlib import Path

import yaml

from tessera.evidence.validation import resolve_pointer
from tessera.identity.normalize import Capabilities
from tessera.identity.reference import ReferenceBundle
from tessera.identity.registry import build_registry
from tessera.schemas.base import canonical_json, sha256_hex
from tessera.schemas.identity import VariantMention
from tessera.sources import curation, ensembl, literature, population
from tessera.sources.http import Fetcher, SourceError

MECHANISM_POLICY = "mechanism_rule_v0.1"


def _present(v, snap):
    return {"status": "value_present", "value": v, "source_snapshot_ids": [snap]}


def _count(v, snap):
    return _present(int(v), snap) if v is not None else {"status": "not_available", "reason_code": "not_reported"}


def fetch(request_path: Path, out: Path, *, log=print) -> Path:
    req = yaml.safe_load(request_path.read_text())
    out.mkdir(parents=True, exist_ok=True)
    cache = Path(__file__).resolve().parents[2] / "cache" / "http"
    f = Fetcher(out / "snapshots", min_interval={"ensembl": 0.1, "clinvar": 0.4}, cache_dir=cache, timeout=180)
    symbol = req["gene"]

    # ---- gene, MANE, reference slice --------------------------------------
    rel = ensembl.release(f)
    gene, gene_snap = ensembl.lookup_gene(f, symbol, rel)
    mane = ensembl.mane_select(gene)
    chrom = str(gene["seq_region_name"])
    accession, _ = ensembl.chrom_accession(f, chrom, rel)
    flank = int(req.get("flank", 2000))
    start, end = int(gene["start"]) - flank, int(gene["end"]) + flank
    seq, _ = ensembl.sequence(f, chrom, start, end, rel)
    log(f"Ensembl {rel}: {symbol} {gene['id']} {chrom}:{gene['start']}-{gene['end']} "
        f"MANE {mane['id']}.{mane['version']}; reference slice {end - start + 1} bp")

    refdir = out / "reference"
    refdir.mkdir(exist_ok=True)
    fasta = f">{chrom} GRCh38 {chrom}:{start}-{end} from Ensembl {rel}\n" + "\n".join(
        seq[i : i + 60] for i in range(0, len(seq), 60)) + "\n"
    (refdir / "reference.fa").write_text(fasta)
    (refdir / "manifest.yaml").write_text(yaml.safe_dump({
        "bundle_id": f"ensembl{rel}_GRCh38_{chrom}_{start}_{end}",
        "assembly": "GRCh38",
        "is_synthetic": False,
        "fasta": "reference.fa",
        "fasta_sha256": sha256_hex(fasta.encode()),
        "contigs": {chrom: {"accession": accession, "start": start}},
    }, sort_keys=True))
    ref = ReferenceBundle.load(refdir / "manifest.yaml")

    # ---- protein regions --------------------------------------------------
    protein_id = mane["Translation"]["id"]
    regions, regions_snap, uniprot_acc = None, None, None
    source_failures: list[str] = []
    try:
        uniprot_acc = ensembl.uniprot_accession(f, protein_id, rel)
        if uniprot_acc:
            entry, regions_snap = curation.uniprot_entry(f, uniprot_acc)
            if entry.get("sequence", {}).get("value") == ensembl.protein_sequence(f, protein_id, rel):
                regions = curation.functional_regions(entry)
            log(f"UniProt {uniprot_acc}: {'sequence identical to' if regions is not None else 'DIFFERS from'} "
                f"{protein_id}; {len(regions or [])} functional regions")
    except SourceError as e:
        # Optional source: domain membership becomes not_available, the run is partial.
        source_failures.append(f"protein_regions: {e}")
        log(f"WARNING protein regions unavailable ({e}); in_functional_domain → not_available")

    # ---- candidate routes -------------------------------------------------
    mentions: list[dict] = []
    clinvar_class: dict[str, list[dict]] = {}

    if req["routes"].get("clinvar"):
        recs = curation.clinvar_records(f, symbol)
        for rec, cv_snap in recs:
            mid = f"clinvar:{rec['accession']}"
            sets = rec.get("variation_set") or []
            spdi = sets[0].get("canonical_spdi") if len(sets) == 1 else None
            row = {"mention_id": mid, "raw_text": f"{rec['accession']} {rec.get('title', '')}".strip(),
                   "route": "clinvar", "source_refs": [f"ClinVar:{rec['accession']}"]}
            if spdi:
                acc, pos0, dele, ins = curation.parse_spdi(spdi)
                if acc == accession and ref.contains(chrom, max(pos0, 1), max(len(dele), 1)):
                    pos, r_, a_ = curation.spdi_to_vcf(pos0, dele, ins, lambda p: ref.fetch(chrom, p, 1))
                    row.update(assembly="GRCh38", chrom=chrom, pos=pos, ref=r_, alt=a_)
            mentions.append(row)
            clinvar_class.setdefault(mid, []).append({**curation.clinvar_classification(rec), "snapshot": cv_snap})
        log(f"ClinVar: {len(recs)} records for {symbol}")

    schema = None
    schema_rows: list[dict] = []
    if req["routes"].get("schema_browser"):
        schema = population.schema_fetch(f, gene["id"])
        if schema["reference_genome"] != "GRCh38":
            raise ValueError("SCHEMA browser dataset is not GRCh38")
        schema_rows = population.schema_variant_rows(schema)
        n = 0
        for r in schema_rows:
            if (r["counts"].get("ac_case") or 0) < 1:
                continue
            c, p, rf, al = r["variant_id"].split("-")
            cn = r["counts"]
            mentions.append({"mention_id": f"schema:{r['variant_id']}", "raw_text": r["variant_id"],
                             "route": "schema_browser", "assembly": "GRCh38", "chrom": c, "pos": int(p),
                             "ref": rf, "alt": al,
                             "source_refs": [f"SCHEMA:case {cn.get('ac_case')}/{cn.get('an_case')} "
                                             f"ctrl {cn.get('ac_ctrl')}/{cn.get('an_ctrl')}"]})
            n += 1
        log(f"SCHEMA: {len(schema_rows)} variants, {n} observed in ≥1 case")

    tx = f"{mane['id']}.{mane['version']}"
    tr = mane["Translation"]
    prot = f"{tr['id']}.{tr['version']}" if tr.get("version") else tr["id"]

    def resolve(items: list[tuple[dict, str, str]], label: str) -> None:
        """items: (mention, notation as written, kind rsid|protein|cdna); same rules for every
        literature-derived route (DEC-21)."""
        pending: dict[str, list[tuple[dict, str, str]]] = {}
        for m, name, kind in items:
            rin = name if kind == "rsid" else literature.recoder_input(name, tx, prot)
            if rin is None:
                m.update(resolution_status="unsupported",
                         resolution_note=f"literature notation {name!r} does not determine a genomic allele")
                mentions.append(m)
            else:
                pending.setdefault(rin, []).append((m, name, kind))
        failed_inputs: set[str] = set()
        mapped = literature.recoder(f, sorted(pending), rel, failed=failed_inputs)
        if failed_inputs:
            source_failures.append(f"{label}: variant_recoder failed for {len(failed_inputs)} inputs")
        counts = {"resolved": 0, "ambiguous": 0, "unmapped": 0, "unsupported": len(items) - sum(map(len, pending.values()))}
        for rin, group in pending.items():
            spdis = [x for x in mapped.get(rin, []) if x.startswith(accession + ":")]
            for m, name, kind in group:
                if rin in failed_inputs:
                    m.update(resolution_status="unmapped", resolution_note="variant_recoder request failed (transient)")
                    counts["unmapped"] += 1
                elif not spdis:
                    m.update(resolution_status="unmapped",
                             resolution_note=f"variant_recoder could not map {rin!r} on {tx}/{prot} (reference or notation mismatch)")
                    counts["unmapped"] += 1
                elif len(spdis) > 1 or not literature.stated_deletion_ok(name, spdis[0]):
                    why = (f"{len(spdis)} compatible genomic alleles: {', '.join(spdis)}" if len(spdis) > 1
                           else f"stated deletion in {name!r} does not match mapped allele {spdis[0]}")
                    m.update(resolution_status="ambiguous", resolution_note=why)
                    counts["ambiguous"] += 1
                else:
                    _, pos0, dele, ins = curation.parse_spdi(spdis[0])
                    if not ref.contains(chrom, max(pos0, 1), max(len(dele), 1)):
                        m.update(resolution_status="unmapped", resolution_note=f"{spdis[0]} outside the reference bundle")
                        counts["unmapped"] += 1
                    else:
                        pos, r_, a_ = curation.spdi_to_vcf(pos0, dele, ins, lambda p: ref.fetch(chrom, p, 1))
                        note = f"{rin} → {spdis[0]}" + ("; c. numbering assumed on MANE transcript" if kind == "cdna" else "")
                        m.update(assembly="GRCh38", chrom=chrom, pos=pos, ref=r_, alt=a_,
                                 resolution_method=f"variant_recoder:{kind}", resolution_note=note)
                        counts["resolved"] += 1
                mentions.append(m)
        log(f"{label}: recoder {counts}")

    if req["routes"].get("litvar"):
        rows, lv_snap = literature.litvar_gene_variants(f, symbol)
        items = []
        for row in rows:
            vid_lv = row["_id"].removeprefix("litvar@")
            name = row.get("rsid") or literature.litvar_name(row)
            if name is None:
                continue  # the gene-level aggregate row
            kind = "rsid" if row.get("rsid") else ("protein" if name.startswith("p.") else "cdna")
            items.append(({"mention_id": f"litvar:{vid_lv}", "route": "litvar", "_litvar_id": row["_id"],
                           "raw_text": f"{name} (LitVar2, {row.get('pmids_count', 0)} PMIDs)"}, name, kind))
        resolve(items, f"LitVar2 ({len(rows)} entries)")
        # PMIDs for alleles that resolved (cached for 30 days; unresolved mentions keep their count only)
        for m, _name, _kind in items:
            lid = m.pop("_litvar_id")
            if m.get("chrom"):
                try:
                    m["source_refs"] = [f"PMID:{p}" for p in literature.litvar_publications(f, lid)]
                except SourceError:
                    m["source_refs"] = ["LitVar2:publications_unavailable"]

    lit_dir = (request_path.parent / req["literature_dir"]).resolve() if req.get("literature_dir") else None
    lit_claims: list[dict] = []
    if lit_dir and req["routes"].get("literature_deep_read"):
        items, seen = [], set()
        hits_f = lit_dir / "supplement_hits.jsonl"
        for h in map(json.loads, hits_f.read_text().splitlines() if hits_f.exists() else []):
            for note, _tx in literature.clean_notations(h["notation"]):
                mid = f"supp:{h['pmid']}:{note}"
                if mid in seen:
                    continue
                seen.add(mid)
                kind = "protein" if note.startswith("p.") else "cdna"
                items.append(({"mention_id": mid, "route": "supplement",
                               "raw_text": f"{h['notation']} (PMID {h['pmid']} {h['file']} {h['where']})",
                               "source_refs": [f"PMID:{h['pmid']}"]}, note, kind))
        claims_f = lit_dir / "deep_claims.jsonl"
        mane_refseq = next((m.get("refseq_match") for m in mane.get("MANE", []) if m.get("type") == "MANE_Select"), None)
        for c in map(json.loads, claims_f.read_text().splitlines() if claims_f.exists() else []):
            lit_claims.append(c)
            c["_mention_ids"] = []
            if not c.get("allele_as_written") or c["scope"] != "variant":
                continue
            notations = literature.clean_notations(c["allele_as_written"])
            if not notations:
                m = {"mention_id": f"lit:{c['pmid']}:{c['allele_as_written']}", "route": "literature",
                     "raw_text": f"{c['allele_as_written']} (PMID {c['pmid']})", "resolution_status": "unsupported",
                     "resolution_note": "no c./p. description with a reference sequence (e.g. bare genomic coordinate without assembly)"}
                if m["mention_id"] not in seen:
                    seen.add(m["mention_id"])
                    mentions.append(m)
                continue
            for note, tx_reported in notations:
                mid = f"lit:{c['pmid']}:{note}"
                c["_mention_ids"].append(mid)
                if mid in seen:
                    continue
                seen.add(mid)
                kind = "protein" if note.startswith("p.") else "cdna"
                m = {"mention_id": mid, "route": "literature", "raw_text": f"{c['allele_as_written']} (PMID {c['pmid']})",
                     "source_refs": [f"PMID:{c['pmid']}"]}
                if tx_reported and mane_refseq and tx_reported != mane_refseq:
                    m["raw_text"] += f" [reported on {tx_reported}; mapped on MANE {mane_refseq}]"
                items.append((m, note, kind))
        resolve(items, f"Literature deep read ({len(lit_claims)} claims, {len(hits_f.read_text().splitlines()) if hits_f.exists() else 0} supplement hits)")

    for m in req.get("manual_variants") or []:
        mentions.append({**m, "route": "manual"})

    cols = ["mention_id", "raw_text", "route", "assembly", "chrom", "pos", "ref", "alt", "roles",
            "resolution_method", "resolution_status", "resolution_note", "source_refs"]
    for m in mentions:
        if isinstance(m.get("source_refs"), list):
            m["source_refs"] = ";".join(m["source_refs"])
    (out / "variants.tsv").write_text(
        "\t".join(cols) + "\n"
        + "".join("\t".join(str(m.get(c, "") if m.get(c) is not None else "") for c in cols) + "\n" for m in mentions)
    )

    # ---- identity (same function freeze will run) -------------------------
    vm = [VariantMention(mention_id=m["mention_id"], raw_text=m["raw_text"], discovery_route=m["route"],
                         reported_assembly=m.get("assembly") or None, reported_chrom=m.get("chrom") or None,
                         reported_pos=m.get("pos") or None, reported_ref=m.get("ref") or None,
                         reported_alt=m.get("alt") or None, roles=m.get("roles", []) or [],
                         resolution_method=m.get("resolution_method"), resolution_status=m.get("resolution_status"),
                         resolution_note=m.get("resolution_note"),
                         source_refs=[x for x in (m.get("source_refs") or "").split(";") if x])
          for m in mentions]
    reg = build_registry(vm, ref, gene["id"], Capabilities())
    alleles = sorted({v.vcf_key for v in reg.variants.values()}, key=lambda k: (k[1], k[2], k[3]))
    log(f"Identity: {len(mentions)} mentions → {len(alleles)} alleles; {len(reg.exclusions)} excluded")

    # ---- annotations ------------------------------------------------------
    vep = ensembl.vep(f, alleles, rel)
    gnomad, gnomad_snap = population.gnomad_gene(f, symbol)
    gnomad_by_vid: dict[str, dict] = {}
    from tessera.identity.normalize import normalize
    for gv in gnomad["variants"]:
        res = normalize(ref, chrom, gv["pos"], gv["ref"], gv["alt"])
        if res.allele:
            gnomad_by_vid[res.allele.variant_id] = gv
    vid_of = {v.vcf_key: v.variant_id for v in reg.variants.values()}
    mentions_of: dict[str, list[str]] = {}
    for mid, vid in reg.mention_to_variant.items():
        mentions_of.setdefault(vid, []).append(mid)

    ann_lines = []
    for key in alleles:
        vid = vid_of[key]
        if key not in vep:
            continue  # annotation failure stays visible as a missing record
        rec, snap = vep[key]
        a = ensembl.parse_vep(rec, snap, gene["id"], mane["id"], regions, regions_snap)
        gv = gnomad_by_vid.get(vid)
        if gv is None:
            a["population_af"] = {"status": "searched_not_found", "retrieval_event_id": gnomad_snap,
                                  "source_snapshot_ids": [gnomad_snap]}
        else:
            af, why = population.gnomad_af(gv)
            a["population_af"] = _present(af, gnomad_snap) if af is not None else \
                {"status": "not_available", "reason_code": why}
        cv = [c for mid in mentions_of.get(vid, []) for c in clinvar_class.get(mid, []) if c["description"]]
        if cv:
            a["clinvar_classification"] = {
                "status": "value_present",
                "value": "; ".join(sorted({f"{c['description']} ({c['review_status']})" for c in cv})),
                "source_snapshot_ids": sorted({c["snapshot"] for c in cv}),
            }
        a["allele"] = {"chrom": key[0], "pos": key[1], "ref": key[2], "alt": key[3]}
        ann_lines.append(a)
    (out / "annotations.jsonl").write_text("".join(json.dumps(a, sort_keys=True) + "\n" for a in ann_lines))

    # ---- SCHEMA observations → claims, links, decisions --------------------
    claims, links, decisions = [], [], []
    cohort = req.get("case_cohorts", {}).get("schema_browser", {})
    if schema:
        doc = {"variants": schema["variants"]}
        snap = schema["snapshots"]["variants"]
        for r in schema_rows:
            mid = f"schema:{r['variant_id']}"
            if mid not in reg.mention_to_variant:
                continue
            cnt = r["counts"]
            pointer = f"/variants/{r['index']}/{schema['variant_fields'].index('group_results')}/0"
            cid = f"c_{mid}"
            claims.append({
                "claim_id": cid, "revision": 1, "evidence_kind": "association", "scope": "variant",
                "subject_gene_ids": [gene["id"]], "source_reported_mention_ids": [mid],
                "claim_text": (f"SCHEMA {r['variant_id']}: {cnt['ac_case']}/{cnt['an_case']} case alleles, "
                               f"{cnt['ac_ctrl']}/{cnt['an_ctrl']} control alleles"
                               + f"; SCHEMA in_analysis={bool(cnt.get('in_analysis'))}"),
                "association": {
                    "cohort_id": "SCHEMA", "study_design": "case_control_exome_meta_analysis",
                    "count_unit": "alleles", "analysis_unit": "allele",
                    "cases_with_allele": _count(cnt.get("ac_case"), snap),
                    "case_total": _count(cnt.get("an_case"), snap),
                    "controls_with_allele": _count(cnt.get("ac_ctrl"), snap),
                    "control_total": _count(cnt.get("an_ctrl"), snap),
                },
                "locators": [{"kind": "json_field", "snapshot_id": snap, "pointer": pointer,
                              "expected_value_digest": sha256_hex(canonical_json(resolve_pointer(doc, pointer)))}],
                "source_snapshot_ids": [snap], "extraction_method": "deterministic",
            })
            links.append({"link_id": f"l_{mid}", "claim_id": cid, "claim_revision": 1, "mention_id": mid,
                          "identity_match": "exact_genomic",
                          "disease_relationship": cohort.get("disease_relationship", "unknown"),
                          "cell_context_relationship": "unknown"})
            decisions.append({"decision_id": f"vd_import_{mid}", "link_id": f"l_{mid}", "claim_id": cid,
                              "claim_revision": 1, "policy_version": "structured-import-1",
                              "reviewer_type": "importer", "reviewer_id": "schema_browser_adapter",
                              "effective_state": "accepted_auto", "sequence": len(decisions) + 1})
    # ---- literature claims (layer 3b) → claims, links, decisions -------------
    copied = set()
    for i, c in enumerate(lit_claims):
        loc = c["locator"]
        snap = loc["snapshot_id"]
        if snap not in copied:  # bring the derived article text into this input's snapshots
            for meta in map(json.loads, (lit_dir / "snapshots" / "snapshots.jsonl").read_text().splitlines()):
                if meta["snapshot_id"] == snap:
                    (out / "snapshots" / meta["file"]).write_bytes((lit_dir / "snapshots" / meta["file"]).read_bytes())
                    f._meta.append(meta)
                    copied.add(snap)
                    break
        mids = [m for m in c.get("_mention_ids", []) if m in reg.mention_to_variant]
        variant_scoped = c["scope"] == "variant" and bool(mids)
        perturb = {"not_an_experiment": "other"}.get(c["perturbation"], c["perturbation"])
        kind = "functional_experiment" if c["kind"] == "functional_experiment" else (
            "patient_observation" if c["kind"] == "patient_observation" else "review_summary")
        cid = f"c_lit_{c['pmid']}_{i}"
        claim = {
            "claim_id": cid, "revision": 1, "evidence_kind": kind,
            "scope": "variant" if variant_scoped else "gene",
            "subject_gene_ids": [gene["id"]], "source_reported_mention_ids": mids if variant_scoped else [],
            "claim_text": c["quote"][:400], "locators": [loc], "source_snapshot_ids": [snap],
            "publication_group_id": f"PMID:{c['pmid']}", "extraction_method": "llm",
        }
        if kind == "functional_experiment":
            claim["functional_assay"] = {
                "assay_name": (c.get("endpoint") or "unspecified endpoint")[:200], "assay_type": "literature_extracted",
                "tested_allele_mention_ids": mids if variant_scoped else [], "perturbation_type": perturb,
                "comparator": c.get("comparator") or "not stated", "endpoint": c.get("endpoint") or "not stated",
                "outcome_type": c["outcome_type"] if c["outcome_type"] != "not_applicable" else "inconclusive",
                "cell_context": f"{c.get('species')}; {c.get('cell_or_model') or 'unspecified'}",
                "confound_flags": ["patient_genetic_background"] if perturb == "carrier_comparison" else [],
            }
        elif kind == "patient_observation":
            na = {"status": "not_available", "reason_code": "count_not_extracted"}
            claim["association"] = {"cohort_id": f"PMID:{c['pmid']}", "study_design": "literature_report",
                                    "count_unit": "carriers", "analysis_unit": "allele" if variant_scoped else "gene",
                                    "cases_with_allele": na, "case_total": na, "controls_with_allele": na,
                                    "control_total": na}
        claims.append(claim)
        if not variant_scoped:
            continue  # gene-level literature claims are kept for audit and context only
        for j, mid in enumerate(mids):
            exact = mid.split(":", 2)[2].startswith("c.")
            lid = f"l_lit_{c['pmid']}_{i}_{j}"
            links.append({"link_id": lid, "claim_id": cid, "claim_revision": 1, "mention_id": mid,
                          "identity_match": "transcript_supported_exact" if exact else "protein_equivalent_only",
                          "disease_relationship": "unknown", "cell_context_relationship": "unknown"})
            decisions.append({"decision_id": f"vd_llm_{lid}", "link_id": lid, "claim_id": cid, "claim_revision": 1,
                              "policy_version": "llm-extract-host-quote-verified-1", "reviewer_type": "model",
                              "reviewer_id": "deep_read", "effective_state": "accepted_auto",
                              "sequence": len(decisions) + 1})

    for name, rows in (("claims", claims), ("links", links), ("decisions", decisions)):
        (out / f"{name}.jsonl").write_text("".join(json.dumps(x, sort_keys=True) + "\n" for x in rows))

    # ---- gene mechanism (declared rule) ------------------------------------
    rule = req["mechanism_rule"]
    basis, mechanism = [], "unknown"
    cg, cg_snap = curation.clingen_row(f, symbol)
    if cg:
        hi = cg.get("Haploinsufficiency Score", "")
        basis.append({"source": "ClinGen dosage", "record_id": symbol, "field": "Haploinsufficiency Score",
                      "value": f"{hi} ({cg.get('Haploinsufficiency Description', '')}; disease "
                               f"{cg.get('Haploinsufficiency Disease ID', '') or 'n/a'})", "snapshot_id": cg_snap})
        if hi in [str(x) for x in rule["clingen_hi_scores_for_lof"]]:
            mechanism = "loss_of_function"
    burden = population.schema_gene_result(schema) if schema else None
    if burden:
        p, orr = burden.get("ptv_p_value"), burden.get("ptv_odds_ratio")
        basis.append({"source": "SCHEMA gene result", "record_id": gene["id"], "field": "ptv_p_value / ptv_odds_ratio",
                      "value": f"{p} / {orr} ({burden.get('ptv_case_carrier')} case vs "
                               f"{burden.get('ptv_control_carrier')} control PTV carriers)",
                      "snapshot_id": schema["snapshots"]["gene"]})
        try:
            if mechanism == "unknown" and float(p) <= float(rule["burden_max_p"]) and float(orr) > 1:
                mechanism = "loss_of_function"
        except (TypeError, ValueError):
            pass
    cons = gnomad.get("gnomad_constraint") or {}
    if cons:
        basis.append({"source": f"gnomAD {population.GNOMAD_DATASET} constraint", "record_id": gene["id"],
                      "field": "oe_lof_upper (LOEUF), pli", "value": f"{cons.get('oe_lof_upper')}, {cons.get('pli')} "
                      "(recorded only; constraint does not decide mechanism)", "snapshot_id": gnomad_snap})
    (out / "mechanism.json").write_text(json.dumps({
        "gene_id": gene["id"], "disease_concept": req["disease"]["concept_id"], "mechanism": mechanism,
        "status": "assessed" if mechanism != "unknown" else "searched_not_found", "basis": basis,
        "decision": "rule", "policy_id": MECHANISM_POLICY,
    }, sort_keys=True))
    log(f"Mechanism: {mechanism}")

    f.write_index()
    policy = (request_path.parent / req["policy"]).resolve()
    (out / "input.yaml").write_text(yaml.safe_dump({
        "spec_version": "2.1.0-draft",
        "gene": {"gene_id": gene["id"], "symbol": symbol, "chrom": chrom,
                 "strand": "+" if gene["strand"] == 1 else "-", "is_synthetic": False},
        "disease": req["disease"],
        "screen_context": req.get("screen_context", {}),
        "reference": "reference/manifest.yaml",
        "policy": str(policy),
        "variants": "variants.tsv",
        "annotations": "annotations.jsonl",
        "snapshots": "snapshots",
        "claims": "claims.jsonl",
        "links": "links.jsonl",
        "decisions": "decisions.jsonl",
        "mechanism": "mechanism.json",
        "scope": {"allowed_consequences": req["scope"]["allowed_consequences"]},
        "observation_search": {"status": "searched" if schema else "not_searched",
                               "scope": "SCHEMA browser variant results (case/control allele counts)" if schema else None},
        "requested_panel_size": req.get("requested_panel_size"),
        "provenance": {"ensembl_release": rel, "mane_select": f"{mane['id']}.{mane['version']}",
                       "uniprot": uniprot_acc, "source_failures": source_failures, "http_requests": f.logical_requests, "http_attempts": f.attempts},
    }, sort_keys=True))
    log(f"HTTP: {f.logical_requests} requests ({f.cache_hits} from cache), {f.attempts} network attempts; "
        f"snapshots in {out / 'snapshots'}")
    return out / "input.yaml"

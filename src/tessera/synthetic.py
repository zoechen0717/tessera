"""Synthetic M0 fixture builder (EVAL B0).

Writes a complete offline input directory: a synthetic reference, allele
mentions, annotations, source snapshots, claims, links and decisions. Every
identifier is visibly synthetic (SYNTH*, synth1). Nothing here describes a real
gene, variant, cohort or publication.

Expected outcomes are defined from the fixture contents and asserted in
tests/test_m0.py; `EXPECT` names the mentions each case is built around.
"""

from __future__ import annotations

import json
import random
from pathlib import Path

import yaml

from tessera.schemas.base import canonical_json, sha256_hex

CONTIG = "synth1"
ACCESSION = "SYNTH_000001.1"
ASSEMBLY = "SYNTH-1"
GENE_ID = "SYNTH:GENE1"
SNAP_TIME = "2026-10-05T00:00:00+00:00"
HOMOPOLYMER_START0 = 300  # 0-based start of a 6×A run used for left-alignment tests

EXPECT = {
    "case_controlled": "m01",
    "nmd_escape_fs_left": "m02",
    "nmd_escape_fs_unaligned": "m03",
    "tie_a": "m04",
    "characterized_altered": "m05",
    "related_condition": "m06",
    "p_unsupported": "m07",
    "tie_b": "m08",
    "multiallelic": "m09",
    "splice_hi": "m10",
    "splice_lo": "m11",
    "user_negative_control": "m12",
    "synonymous_discovery": "m13",
    "ref_mismatch": "m14",
    "protein_only": "m15",
    "rsid_only": "m16",
    "bad_offsets": "m17",
    "human_rejected": "m18",
    "carrier_only": "m19",
    "wrong_assembly": "m20",
    "intro_mention_only": "m04",
    "nmd_trigger_ptv": "m21",
}


def _sequence() -> str:
    rng = random.Random(20261005)
    seq = [rng.choice("ACGT") for _ in range(1200)]
    for i in range(HOMOPOLYMER_START0, HOMOPOLYMER_START0 + 6):
        seq[i] = "A"
    seq[HOMOPOLYMER_START0 - 1] = "C"
    seq[HOMOPOLYMER_START0 + 6] = "G"
    return "".join(seq)


def _alt_for(ref: str) -> str:
    return {"A": "G", "C": "T", "G": "A", "T": "C"}[ref]


def build(dest: Path, mechanism: str = "loss_of_function", policy_path: Path | None = None) -> Path:
    """Write the fixture to `dest` and return the path of input.yaml."""
    dest.mkdir(parents=True, exist_ok=True)
    seq = _sequence()

    def base(pos1: int) -> str:
        return seq[pos1 - 1]

    # ---- reference bundle
    refdir = dest / "reference"
    refdir.mkdir(exist_ok=True)
    fasta = f">{CONTIG} synthetic contig, not a real sequence\n" + "\n".join(
        seq[i : i + 60] for i in range(0, len(seq), 60)
    ) + "\n"
    (refdir / "reference.fa").write_text(fasta)
    (refdir / "manifest.yaml").write_text(
        yaml.safe_dump(
            {
                "bundle_id": "synthetic_ref_v1",
                "assembly": ASSEMBLY,
                "is_synthetic": True,
                "fasta": "reference.fa",
                "fasta_sha256": sha256_hex(fasta.encode()),
                "contigs": {CONTIG: ACCESSION},
            },
            sort_keys=True,
        )
    )

    # ---- alleles: position -> (ref, alt)
    snv_pos = {f"m{n:02d}": 100 + 23 * n for n in range(1, 22)}
    alleles: dict[str, tuple[int, str, str]] = {}
    for mid, pos in snv_pos.items():
        alleles[mid] = (pos, base(pos), _alt_for(base(pos)))
    h = HOMOPOLYMER_START0  # 1-based position of the base before the run
    alleles["m02"] = (h, base(h) + "A", base(h))           # left-aligned 1-bp deletion
    alleles["m03"] = (h + 5, "AA", "A")                    # same deletion, not left-aligned
    assert seq[h + 4 : h + 6] == "AA"

    rows = []
    for mid in sorted(snv_pos):
        pos, ref, alt = alleles[mid]
        raw, chrom, assembly, roles = f"{CONTIG}:{pos}:{ref}>{alt}", CONTIG, "", ""
        if mid == "m09":
            alts = [a for a in "ACGT" if a != ref][:2]
            alt = ",".join(alts)
            raw = f"{CONTIG}:{pos}:{ref}>{alt}"
        if mid == "m14":
            true_ref, ref = ref, _alt_for(ref)  # deliberately wrong REF
            alt = next(b for b in "ACGT" if b not in (ref, true_ref))
            raw = f"{CONTIG}:{pos}:{ref}>{alt}"
        if mid == "m15":
            raw, chrom, pos, ref, alt = "p.Arg45Cys", "", "", "", ""
        if mid == "m16":
            raw, chrom, pos, ref, alt = "rs999999999", "", "", "", ""
        if mid == "m20":
            assembly = "GRCh37"
        if mid == "m12":
            roles = "negative_control"
        rows.append([mid, raw, assembly, chrom, str(pos), ref, alt, roles])
    (dest / "variants.tsv").write_text(
        "## SYNTHETIC alleles on a synthetic contig\n"
        + "\t".join(["mention_id", "raw_text", "assembly", "chrom", "pos", "ref", "alt", "roles"]) + "\n"
        + "".join("\t".join(r) + "\n" for r in rows)
    )
    m09_pos, m09_ref, _ = alleles["m09"]
    m09_alts = [a for a in "ACGT" if a != m09_ref][:2]

    # ---- snapshots
    snapdir = dest / "snapshots"
    snapdir.mkdir(exist_ok=True)
    snaps_meta = []

    def snapshot(name: str, content: bytes, media: str, source: str) -> str:
        digest = sha256_hex(content)
        (snapdir / name).write_bytes(content)
        snaps_meta.append({
            "file": name, "snapshot_id": "snap_" + digest, "source": source,
            "source_record_id": f"SYNTH-{name}", "source_version": "synthetic-1",
            "retrieved_at": SNAP_TIME, "media_type": media, "raw_content_digest": digest,
            "license_status": "synthetic", "is_synthetic": True,
        })
        return "snap_" + digest

    cohort = {"cohort": "SYNTH-COHORT-1", "rows": [
        {"allele": "m01", "cases": 3, "case_total": 1000, "controls": 0, "control_total": 2000},
        {"allele": "m17", "cases": 2, "case_total": 1000, "controls": 0, "control_total": 2000},
        {"allele": "m18", "cases": 1, "case_total": 1000, "controls": 0, "control_total": 2000},
    ]}
    related = {"cohort": "SYNTH-RELATED-COHORT", "condition": "synthetic related disorder",
               "rows": [{"allele": "m06", "cases": 1, "case_total": 500}]}
    paper = (
        "SYNTHETIC PAPER. Introduction: prior work described allele m04 in passing. "
        "Results: we installed allele m05 by exact editing in synthetic neurons; the reporter "
        "signal was decreased relative to isogenic wild type. "
        "In a separate experiment, knockout of SYNGENE1 decreased the reporter signal. "
        "Carriers of allele m19 showed decreased signal compared with unrelated controls. "
        "Allele m17 was observed in two cases."
    )
    s_cohort = snapshot("cohort.json", canonical_json(cohort).encode(), "application/json", "synthetic_cohort")
    s_related = snapshot("related.json", canonical_json(related).encode(), "application/json", "synthetic_cohort")
    s_paper = snapshot("paper.txt", paper.encode(), "text/plain", "synthetic_paper")
    s_annot = snapshot("annotation_source.json", b'{"note":"synthetic annotation source"}', "application/json", "synthetic_annotator")
    (snapdir / "snapshots.jsonl").write_text("".join(json.dumps(m, sort_keys=True) + "\n" for m in snaps_meta))

    # ---- annotations (by normalized allele)
    def present(v):
        return {"status": "value_present", "value": v, "source_snapshot_ids": [s_annot]}

    def eff(*terms):
        return [{"transcript_id": "SYNTX1.1", "consequences": list(terms), "selected_display": True}]

    def ann(mid, terms, nmd=None, domain=None, preds=None, allele=None):
        pos, ref, alt = allele or alleles[mid]
        rec = {"allele": {"chrom": CONTIG, "pos": pos, "ref": ref, "alt": alt}, "transcript_effects": eff(*terms)}
        if nmd is not None:
            rec["nmd"] = present(nmd)
        if domain is not None:
            rec["in_functional_domain"] = present(domain)
        if preds:
            rec["predictions"] = preds
        return rec

    am = lambda v: {"alphamissense": present(v)}  # noqa: E731
    sd = lambda v: {"splice_delta": present(v)}  # noqa: E731
    annotations = [
        ann("m01", ["stop_gained"], nmd="triggers"),
        ann("m21", ["stop_gained"], nmd="triggers"),
        ann("m02", ["frameshift_variant"], nmd="escapes"),
        ann("m04", ["missense_variant"], domain=True, preds=am(0.91)),
        ann("m05", ["missense_variant"], domain=False, preds=am(0.95)),
        ann("m06", ["missense_variant"], domain=True, preds=am(0.30)),
        ann("m07", ["missense_variant"], domain=True, preds={"alphamissense": {
            "status": "unsupported", "reason_code": "allele_not_in_source_release"}}),
        ann("m08", ["missense_variant"], domain=True, preds=am(0.91)),
        ann("m09", ["missense_variant"], domain=True, preds=am(0.50), allele=(m09_pos, m09_ref, m09_alts[0])),
        ann("m09", ["missense_variant"], domain=True, preds=am(0.60), allele=(m09_pos, m09_ref, m09_alts[1])),
        ann("m10", ["splice_donor_variant"], preds=sd(0.80)),
        ann("m11", ["splice_donor_variant"], preds=sd(0.20)),
        ann("m12", ["synonymous_variant"], preds=sd(0.10)),
        ann("m13", ["synonymous_variant"], preds=sd(0.05)),
        ann("m17", ["missense_variant"], domain=True, preds=am(0.70)),
        ann("m18", ["missense_variant"], domain=False, preds=am(0.40)),
        ann("m19", ["missense_variant"], domain=False, preds=am(0.20)),
    ]
    (dest / "annotations.jsonl").write_text("".join(json.dumps(a, sort_keys=True) + "\n" for a in annotations))

    # ---- claims
    def jloc(snap, doc, pointer):
        from tessera.evidence.validation import resolve_pointer
        return {"kind": "json_field", "snapshot_id": snap, "pointer": pointer,
                "expected_value_digest": sha256_hex(canonical_json(resolve_pointer(doc, pointer)))}

    def tloc(snap, text, excerpt, shift=0):
        start = text.index(excerpt) + shift
        return {"kind": "text_span", "snapshot_id": snap, "text_digest": sha256_hex(text),
                "parser_version": "plain-1", "start": start, "end": start + len(excerpt), "excerpt": excerpt}

    def count(v, snap):
        if v is None:
            return {"status": "not_available", "reason_code": "not_reported"}
        return {"status": "value_present", "value": v, "source_snapshot_ids": [snap]}

    def obs_claim(cid, mid, snap, doc, idx, loc=None):
        row = doc["rows"][idx]
        return {
            "claim_id": cid, "revision": 1, "evidence_kind": "patient_observation", "scope": "variant",
            "subject_gene_ids": [GENE_ID], "source_reported_mention_ids": [mid],
            "claim_text": f"Allele {mid} observed in {row['cases']} case(s) in {doc['cohort']}.",
            "association": {
                "cohort_id": doc["cohort"], "study_design": "case_control_counts",
                "cases_with_allele": count(row["cases"], snap), "case_total": count(row.get("case_total"), snap),
                "controls_with_allele": count(row.get("controls"), snap),
                "control_total": count(row.get("control_total"), snap), "analysis_unit": "allele",
            },
            "locators": [loc or jloc(snap, doc, f"/rows/{idx}")], "source_snapshot_ids": [snap],
            "extraction_method": "deterministic",
        }

    def func_claim(cid, scope, tested, ptype, excerpt, reported):
        return {
            "claim_id": cid, "revision": 1, "evidence_kind": "functional_experiment", "scope": scope,
            "subject_gene_ids": [GENE_ID], "source_reported_mention_ids": reported,
            "claim_text": excerpt,
            "functional_assay": {
                "assay_name": "synthetic_reporter", "assay_type": "reporter", "tested_allele_mention_ids": tested,
                "perturbation_type": ptype, "comparator": "isogenic wild type" if ptype != "carrier_comparison" else "unrelated controls",
                "endpoint": "reporter signal", "outcome_type": "altered", "cell_context": "synthetic neurons",
                "confound_flags": ["patient_genetic_background"] if ptype == "carrier_comparison" else [],
            },
            "locators": [tloc(s_paper, paper, excerpt)], "source_snapshot_ids": [s_paper],
            "extraction_method": "llm",
        }

    claims = [
        obs_claim("c_obs1", "m01", s_cohort, cohort, 0),
        obs_claim("c_obs2", "m06", s_related, related, 0),
        {**obs_claim("c_obs3", "m17", s_cohort, cohort, 1),
         "locators": [tloc(s_paper, paper, "Allele m17 was observed in two cases.", shift=7)],
         "source_snapshot_ids": [s_cohort, s_paper]},
        obs_claim("c_obs4", "m18", s_cohort, cohort, 2),
        func_claim("c_func1", "variant", ["m05"], "exact_allele_edit",
                   "we installed allele m05 by exact editing in synthetic neurons", ["m04", "m05"]),
        func_claim("c_ko", "gene", [], "gene_knockout",
                   "knockout of SYNGENE1 decreased the reporter signal", []),
        func_claim("c_carrier", "variant", ["m19"], "carrier_comparison",
                   "Carriers of allele m19 showed decreased signal", ["m19"]),
    ]
    (dest / "claims.jsonl").write_text("".join(json.dumps(c, sort_keys=True) + "\n" for c in claims))

    # ---- links (variant_id and identity decisions are filled in by the host)
    def link(lid, cid, mid, match, disease="exact"):
        return {"link_id": lid, "claim_id": cid, "claim_revision": 1, "mention_id": mid,
                "identity_match": match, "disease_relationship": disease,
                "cell_context_relationship": "approved"}

    links = [
        link("l_obs1", "c_obs1", "m01", "exact_genomic"),
        link("l_obs2", "c_obs2", "m06", "exact_genomic", disease="related"),
        link("l_obs3", "c_obs3", "m17", "exact_genomic"),
        link("l_obs4", "c_obs4", "m18", "exact_genomic"),
        link("l_func1", "c_func1", "m05", "exact_genomic"),
        link("l_intro", "c_func1", "m04", "exact_genomic"),     # mentioned, not tested
        link("l_ko_gene", "c_ko", "m09.1", "gene_only"),
        link("l_ko_bad", "c_ko", "m09.2", "exact_genomic"),     # gene claim posing as exact
        link("l_carrier", "c_carrier", "m19", "exact_genomic"),
    ]
    (dest / "links.jsonl").write_text("".join(json.dumps(l, sort_keys=True) + "\n" for l in links))

    def decision(seq, lid, cid, state, reviewer="model", human="unreviewed"):
        return {"decision_id": f"vd_{seq}_{lid}", "link_id": lid, "claim_id": cid, "claim_revision": 1,
                "policy_version": "synthetic-review-1", "reviewer_type": reviewer,
                "reviewer_id": "synthetic_reviewer" if reviewer != "model" else "synthetic-model",
                "human_review": human, "effective_state": state, "sequence": seq}

    decisions = [
        decision(1, "l_obs1", "c_obs1", "accepted_auto"),
        decision(2, "l_obs2", "c_obs2", "accepted_auto"),
        decision(3, "l_obs3", "c_obs3", "accepted_auto"),
        decision(4, "l_obs4", "c_obs4", "accepted_auto"),
        decision(5, "l_obs4", "c_obs4", "rejected", reviewer="human", human="rejected"),
        decision(6, "l_func1", "c_func1", "accepted_human", reviewer="human", human="accepted"),
        decision(7, "l_intro", "c_func1", "accepted_auto"),
        decision(8, "l_ko_gene", "c_ko", "accepted_auto"),
        decision(9, "l_ko_bad", "c_ko", "accepted_auto"),
        decision(10, "l_carrier", "c_carrier", "accepted_auto"),
    ]
    (dest / "decisions.jsonl").write_text("".join(json.dumps(d, sort_keys=True) + "\n" for d in decisions))

    # ---- gene mechanism
    (dest / "mechanism.json").write_text(json.dumps({
        "gene_id": GENE_ID, "disease_concept": "SYNTH:D1", "mechanism": mechanism,
        "status": "assessed" if mechanism != "unknown" else "searched_not_found",
        "basis": [{"source": "synthetic_dosage_curation", "record_id": "SYNTH-DOSAGE-1",
                   "field": "mechanism", "value": mechanism}],
        "decision": "imported", "policy_id": "synthetic-mechanism-1",
    }, sort_keys=True))

    if policy_path is None:
        policy_path = Path(__file__).resolve().parents[2] / "config" / "nomination_v0.1.yaml"
    input_yaml = {
        "spec_version": "2.1.0-draft",
        "gene": {"gene_id": GENE_ID, "symbol": "SYNGENE1", "chrom": CONTIG, "strand": "-", "is_synthetic": True},
        "disease": {"input": "synthetic disorder", "concept_id": "SYNTH:D1"},
        "screen_context": {"species": "synthetic", "cell_type": "synthetic neuron", "assay": "synthetic",
                           "editing_mode": "prime_editing", "host_background": None},
        "reference": "reference/manifest.yaml",
        "policy": str(policy_path),
        "variants": "variants.tsv",
        "annotations": "annotations.jsonl",
        "snapshots": "snapshots",
        "claims": "claims.jsonl",
        "links": "links.jsonl",
        "decisions": "decisions.jsonl",
        "mechanism": "mechanism.json",
        "observation_search": {"status": "searched", "scope": "synthetic cohort tables"},
        "requested_panel_size": None,
    }
    (dest / "input.yaml").write_text(yaml.safe_dump(input_yaml, sort_keys=True))
    return dest / "input.yaml"

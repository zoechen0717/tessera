"""Ensembl REST: gene/MANE lookup, chromosome accession, reference sequence, VEP.

Verified 2026-10-05 against rest.ensembl.org (release 116): /lookup/symbol with
mane=1, /info/assembly synonyms (RefSeq_genomic), /sequence/region, and POST
/vep/homo_sapiens/region returning per-transcript consequence, HGVS, an
`alphamissense` object, a `spliceai` object (DS_AG/DS_AL/DS_DG/DS_DL) and
LOFTEE `lof`/`lof_info`. VEP does NOT validate REF; tessera validates REF
against the reference bundle before calling it.
"""

from __future__ import annotations

from tessera.sources.http import Fetcher, SourceError

BASE = "https://rest.ensembl.org"
SOURCE = "ensembl"
VEP_BATCH = 200
VEP_OPTIONS = {"mane": 1, "hgvs": 1, "canonical": 1, "AlphaMissense": 1, "SpliceAI": 1, "LoF": 1}
PTV_TERMS = {"stop_gained", "frameshift_variant"}


def release(f: Fetcher) -> str:
    r = f.get(SOURCE, "info/data", "live", f"{BASE}/info/data", {"content-type": "application/json"})
    return str(max(r.json()["releases"]))


def lookup_gene(f: Fetcher, symbol: str, rel: str) -> tuple[dict, str]:
    r = f.get(SOURCE, f"lookup/symbol/{symbol}", rel, f"{BASE}/lookup/symbol/homo_sapiens/{symbol}",
              {"content-type": "application/json", "expand": 1, "mane": 1})
    return r.json(), r.snapshot_id


def mane_select(gene: dict) -> dict:
    hits = [t for t in gene.get("Transcript", []) for m in t.get("MANE", []) if m.get("type") == "MANE_Select"]
    if len(hits) != 1:
        raise SourceError(SOURCE, "mane_ambiguous", f"{len(hits)} MANE Select transcripts for {gene.get('id')}")
    return hits[0]


def chrom_accession(f: Fetcher, chrom: str, rel: str) -> tuple[str, str]:
    r = f.get(SOURCE, f"info/assembly/{chrom}", rel, f"{BASE}/info/assembly/homo_sapiens/{chrom}",
              {"content-type": "application/json", "synonyms": 1})
    refseq = [s["name"] for s in r.json().get("synonyms", []) if s.get("dbname") == "RefSeq_genomic"]
    if len(refseq) != 1:
        raise SourceError(SOURCE, "accession_missing", f"no unique RefSeq accession for {chrom}")
    return refseq[0], r.snapshot_id


def sequence(f: Fetcher, chrom: str, start: int, end: int, rel: str) -> tuple[str, str]:
    r = f.get(SOURCE, f"sequence/{chrom}:{start}-{end}", rel,
              f"{BASE}/sequence/region/human/{chrom}:{start}..{end}:1", {"content-type": "application/json"})
    seq = r.json()["seq"].upper()
    if len(seq) != end - start + 1:
        raise SourceError(SOURCE, "invalid_source_payload", "sequence length does not match region")
    return seq, r.snapshot_id


def protein_sequence(f: Fetcher, protein_id: str, rel: str) -> str:
    r = f.get(SOURCE, f"sequence/{protein_id}", rel, f"{BASE}/sequence/id/{protein_id}",
              {"content-type": "application/json", "type": "protein"})
    return r.json()["seq"]


def uniprot_accession(f: Fetcher, protein_id: str, rel: str) -> str | None:
    r = f.get(SOURCE, f"xrefs/{protein_id}", rel, f"{BASE}/xrefs/id/{protein_id}",
              {"content-type": "application/json", "external_db": "Uniprot/SWISSPROT"})
    ids = sorted({x["primary_id"] for x in r.json()})
    return ids[0] if len(ids) == 1 else None


def vep(f: Fetcher, alleles: list[tuple[str, int, str, str]], rel: str) -> dict[tuple, tuple[dict, str]]:
    """Return {(chrom, pos, ref, alt): (vep_record, snapshot_id)}."""
    out: dict[tuple, tuple[dict, str]] = {}
    for i in range(0, len(alleles), VEP_BATCH):
        batch = alleles[i : i + VEP_BATCH]
        inputs = {f"{c} {p} . {r} {a} . . .": (c, p, r, a) for c, p, r, a in batch}
        resp = f.post_json(SOURCE, f"vep/batch{i // VEP_BATCH}", rel, f"{BASE}/vep/homo_sapiens/region",
                           {"variants": list(inputs), **VEP_OPTIONS})
        for rec in resp.json():
            key = inputs.get(rec.get("input"))
            if key is not None:
                out[key] = (rec, resp.snapshot_id)
    return out


def _present(v, snap):
    return {"status": "value_present", "value": v, "source_snapshot_ids": [snap]}


def _absent(status: str, reason: str):
    return {"status": status, "reason_code": reason}


def parse_vep(rec: dict, snap: str, gene_id: str, mane_tid: str, protein_regions: list[tuple[int, int, str]] | None,
              protein_regions_snap: str | None) -> dict:
    """VEP record → VariantAnnotation fields (without `allele`/`variant_id`).

    NMD rule (declared): for PTV consequences on the MANE transcript, LOFTEE's
    50_BP_RULE PASS → `triggers`, FAIL → `escapes`; absent → not_available.
    splice_delta = max(DS_AG, DS_AL, DS_DG, DS_DL) on the MANE transcript.
    """
    tcs = [t for t in rec.get("transcript_consequences", []) if t.get("gene_id") == gene_id]
    effects = []
    mane = None
    for t in tcs:
        tid = (t.get("hgvsc") or "").split(":")[0] or t["transcript_id"]
        selected = t["transcript_id"] == mane_tid
        if selected:
            mane = t
        effects.append({
            "transcript_id": tid,
            "hgvs_c": t.get("hgvsc"),
            "hgvs_p": t.get("hgvsp"),
            "consequences": t["consequence_terms"],
            "selected_display": selected,
        })
    out: dict = {"transcript_effects": sorted(effects, key=lambda e: e["transcript_id"])}
    if mane is None:
        return out
    terms = set(mane["consequence_terms"])

    if terms & PTV_TERMS:
        info = dict(kv.split(":", 1) for kv in (mane.get("lof_info") or "").split(",") if ":" in kv)
        rule = info.get("50_BP_RULE")
        out["nmd"] = (_present({"PASS": "triggers", "FAIL": "escapes"}[rule], snap)
                      if rule in ("PASS", "FAIL") else _absent("not_available", "loftee_50bp_rule_absent"))
    else:
        out["nmd"] = _absent("not_applicable", "not_a_truncating_consequence")

    preds = {}
    am = mane.get("alphamissense")
    if "missense_variant" in terms:
        preds["alphamissense"] = (_present(float(am["am_pathogenicity"]), snap) if am
                                  else _absent("not_available", "no_alphamissense_score"))
    sai = mane.get("spliceai")
    ds = [sai.get(k) for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")] if sai else []
    preds["splice_delta"] = (_present(float(max(ds)), snap) if ds and all(v is not None for v in ds)
                             else _absent("not_available", "no_precomputed_spliceai_score"))
    out["predictions"] = preds

    ps, pe = mane.get("protein_start"), mane.get("protein_end")
    if protein_regions is None or ps is None:
        out["in_functional_domain"] = _absent("not_applicable", "no_protein_position") if ps is None \
            else _absent("not_available", "protein_isoform_not_mapped")
    else:
        hits = sorted({name for (s, e, name) in protein_regions if s <= (pe or ps) and ps <= e})
        out["in_functional_domain"] = {"status": "value_present", "value": bool(hits),
                                       "source_snapshot_ids": [protein_regions_snap]}
        out["domain_name"] = "; ".join(hits) or None
    return out

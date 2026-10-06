"""gnomAD (GraphQL) and the SCHEMA results browser.

gnomAD — verified 2026-10-05: POST https://gnomad.broadinstitute.org/api,
`gene(gene_symbol, reference_genome: GRCh38) { gnomad_constraint {...}
variants(dataset: gnomad_r4) {...} }`.

SCHEMA browser — verified 2026-10-05: GET https://schema.broadinstitute.org
/config.js (field names), /api/gene/<ENSG> (gene results incl. burden),
/api/gene/<ENSG>/variants (rows of [variant_id, pos, consequence, hgvsc, hgvsp,
info[], group_results[[ac_case, an_case, ac_ctrl, an_ctrl, n_de_novo,
in_analysis]]]). Counts are allele counts, not carriers. Field names are read
from config.js at fetch time rather than assumed.
"""

from __future__ import annotations

import json

from tessera.sources.http import Fetcher, SourceError

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET = "gnomad_r4"
SCHEMA_BASE = "https://schema.broadinstitute.org"

_GNOMAD_QUERY = """
query($symbol: String!) {
  gene(gene_symbol: $symbol, reference_genome: GRCh38) {
    gene_id
    gnomad_constraint { pli oe_lof oe_lof_lower oe_lof_upper }
    variants(dataset: %s) { variant_id pos ref alt exome { ac an filters } genome { ac an filters } }
  }
}""" % GNOMAD_DATASET


def gnomad_gene(f: Fetcher, symbol: str) -> tuple[dict, str]:
    r = f.post_json("gnomad", f"gene/{symbol}", GNOMAD_DATASET, GNOMAD_API,
                    {"query": _GNOMAD_QUERY, "variables": {"symbol": symbol}}, ttl_days=30)
    d = r.json()
    if d.get("errors"):
        raise SourceError("gnomad", "invalid_source_payload", str(d["errors"])[:300])
    return d["data"]["gene"], r.snapshot_id


def gnomad_af(variant: dict) -> tuple[float | None, str | None]:
    """Joint AF over the exome and genome calls that pass filters.

    Returns (af, None) or (None, reason). A call failing filters is not used.
    """
    ac = an = 0
    for part in ("exome", "genome"):
        call = variant.get(part)
        if call and not call.get("filters") and call.get("an"):
            ac += call["ac"]
            an += call["an"]
    if an == 0:
        return None, "no_passing_call"
    return ac / an, None


def schema_fetch(f: Fetcher, gene_id: str) -> dict:
    cfg = f.get("schema_browser", "config.js", "live", f"{SCHEMA_BASE}/config.js", media="text/plain", ttl_days=30)
    text = cfg.text()
    config = json.loads(text.split("=", 1)[1].strip().rstrip(";"))
    ds = config["datasets"][config["datasetId"]]
    gene = f.get("schema_browser", f"gene/{gene_id}", "live", f"{SCHEMA_BASE}/api/gene/{gene_id}", ttl_days=30)
    variants = f.get("schema_browser", f"gene/{gene_id}/variants", "live",
                     f"{SCHEMA_BASE}/api/gene/{gene_id}/variants", ttl_days=30)
    return {
        "dataset_id": config["datasetId"],
        "reference_genome": ds["reference_genome"],
        "variant_fields": config["variant_fields"],
        "group_fields": ds["variant_group_result_field_names"],
        "gene_fields": ds["gene_group_result_field_names"],
        "gene": gene.json()["gene"],
        "variants": variants.json()["variants"],
        "snapshots": {"config": cfg.snapshot_id, "gene": gene.snapshot_id, "variants": variants.snapshot_id},
    }


def schema_variant_rows(s: dict) -> list[dict]:
    """Rows as dicts, keeping the index for JSON Pointer locators."""
    vf, gf = s["variant_fields"], s["group_fields"]
    out = []
    for i, row in enumerate(s["variants"]):
        rec = dict(zip(vf, row))
        groups = rec.get("group_results") or []
        if len(groups) != 1:
            continue
        rec["counts"] = dict(zip(gf, groups[0]))
        rec["index"] = i
        out.append(rec)
    return out


def schema_gene_result(s: dict) -> dict | None:
    res = (s["gene"].get("gene_results") or {}).get(s["dataset_id"])
    if not res or len(res.get("group_results") or []) != 1:
        return None
    return dict(zip(s["gene_fields"], res["group_results"][0]))

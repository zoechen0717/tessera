"""M1a source parsers against verbatim excerpts of real responses (no network).

tests/fixtures/real_excerpts.json was captured by `tessera fetch` for SETD1A on
2026-10-05; expected values below are read from those excerpts, not invented.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tessera.sources import curation, ensembl, population

X = json.loads((Path(__file__).parent / "fixtures" / "real_excerpts.json").read_text())
GENE, MANE = "ENSG00000099381", "ENST00000262519"
REGIONS = [(84, 172, "Domain:RRM"), (1568, 1685, "Domain:SET")]


def _vep(kind, regions=REGIONS):
    return ensembl.parse_vep(X["vep"][kind], "snap_x", GENE, MANE, regions, "snap_u")


def test_missense_takes_alphamissense_from_mane_transcript():
    a = _vep("missense")
    sel = [e for e in a["transcript_effects"] if e["selected_display"]]
    assert len(sel) == 1 and sel[0]["transcript_id"] == "ENST00000262519.14"
    assert a["predictions"]["alphamissense"] == {
        "status": "value_present", "value": 0.0545, "source_snapshot_ids": ["snap_x"]}
    assert a["nmd"]["status"] == "not_applicable"
    assert a["in_functional_domain"]["value"] is False  # residue 294 is outside RRM/SET


def test_frameshift_nmd_from_loftee_50bp_rule():
    a = _vep("frameshift")
    assert a["nmd"]["value"] == "triggers"
    assert "alphamissense" not in a["predictions"]


def test_intronic_splice_has_no_protein_position():
    a = _vep("splice")
    assert a["in_functional_domain"]["status"] == "not_applicable"
    assert a["predictions"]["splice_delta"]["status"] in ("value_present", "not_available")


def test_unmapped_isoform_is_not_available_not_false():
    a = _vep("synonymous", regions=None)
    assert a["in_functional_domain"] == {"status": "not_available", "reason_code": "protein_isoform_not_mapped"}


@pytest.mark.parametrize("kind", ["snv", "del", "other"])
def test_spdi_to_vcf_matches_clinvar_grch38_location(kind):
    rec = next(c for c in X["clinvar"] if c["kind"] == kind)
    acc, pos0, dele, ins = curation.parse_spdi(rec["canonical_spdi"])
    assert acc == "NC_000016.10"
    pos, ref, alt = curation.spdi_to_vcf(pos0, dele, ins, lambda p: "N")
    grch38 = int(next(l for l in rec["variation_loc"] if l["assembly_name"] == "GRCh38")["start"])
    if kind == "snv":
        assert (pos, len(ref), len(alt)) == (grch38, 1, 1)
    elif kind == "del":
        assert pos == pos0 and ref == "N" + dele and alt == "N"  # anchored at the base before
    else:
        assert pos == pos0 + 1 and ref == dele and alt == ins  # expanded SPDI; normalization trims later


def test_schema_rows_bind_field_names_from_config():
    s = {"variant_fields": X["schema"]["variant_fields"], "group_fields": X["schema"]["group_fields"],
         "variants": X["schema"]["rows"]}
    rows = population.schema_variant_rows(s)
    assert rows[0]["variant_id"] == X["schema"]["rows"][0][0]
    assert set(rows[0]["counts"]) == {"ac_case", "an_case", "ac_ctrl", "an_ctrl", "n_de_novo", "in_analysis"}


def test_gnomad_filtered_calls_are_not_frequencies():
    filtered, clean = X["gnomad"]
    af, why = population.gnomad_af({**filtered, "genome": None})
    assert af is None and why == "no_passing_call"
    af, why = population.gnomad_af(clean)
    total_an = clean["exome"]["an"] + clean["genome"]["an"]
    assert why is None and af == pytest.approx((clean["exome"]["ac"] + clean["genome"]["ac"]) / total_an)

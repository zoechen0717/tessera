"""Identity: INV-01, INV-15, R-02 (ISSUE-002)."""

from __future__ import annotations

import json

import pytest
from pydantic import ValidationError

from tessera.identity.normalize import normalize
from tessera.identity.reference import ReferenceBundle, ReferenceError_
from tessera.schemas.identity import CanonicalVariant, IdentityDecision, IdentityStatus
from tessera.schemas.observed import Observed, ObservationStatus
from tessera.synthetic import CONTIG, HOMOPOLYMER_START0

from conftest import mention_label


def _ref(fixture_input):
    return ReferenceBundle.load(fixture_input.parent / "reference" / "manifest.yaml")


def test_equivalent_indel_representations_merge(run_bundle):
    a, b = mention_label(run_bundle, "m02"), mention_label(run_bundle, "m03")
    assert a == b == f"{CONTIG}:{HOMOPOLYMER_START0}:CA>C"
    variants = [json.loads(l) for l in (run_bundle / "candidate_variants.jsonl").read_text().splitlines()]
    merged = [v for v in variants if f"{v['chrom']}:{v['pos_1based']}:{v['ref']}>{v['alt']}" == a]
    assert len(merged) == 1
    assert {al["mention_id"] for al in merged[0]["aliases"]} == {"m02", "m03"}


def test_multiallelic_record_stays_two_alleles(run_bundle):
    a, b = mention_label(run_bundle, "m09.1"), mention_label(run_bundle, "m09.2")
    assert a and b and a != b


@pytest.mark.parametrize(
    "mention,status",
    [("m14", "ref_mismatch"), ("m15", "unsupported"), ("m16", "insufficient"), ("m20", "unsupported")],
)
def test_unresolvable_mentions_are_excluded_not_dropped(run_bundle, mention, status):
    excl = {e["mention_id"]: e for e in map(json.loads, (run_bundle / "exclusions.jsonl").read_text().splitlines())}
    assert excl[mention]["status"] == status
    assert mention_label(run_bundle, mention) is None
    mentions = {m["mention_id"] for m in map(json.loads, (run_bundle / "candidate_mentions.jsonl").read_text().splitlines())}
    assert mention in mentions


def test_left_alignment_and_ref_validation(fixture_input):
    ref = _ref(fixture_input)
    h = HOMOPOLYMER_START0
    r = normalize(ref, CONTIG, h + 5, "AA", "A")
    assert r.status is IdentityStatus.RESOLVED and (r.allele.pos_1based, r.allele.ref, r.allele.alt) == (h, "CA", "C")
    ins = normalize(ref, CONTIG, h + 6, "A", "AA")
    assert (ins.allele.pos_1based, ins.allele.ref, ins.allele.alt) == (h, "C", "CA")
    wrong = "A" if ref.fetch(CONTIG, 50, 1) != "A" else "C"
    assert normalize(ref, CONTIG, 50, wrong, "T" if wrong != "T" else "G").status is IdentityStatus.REF_MISMATCH


def test_reference_checksum_is_enforced(fixture_input):
    fa = fixture_input.parent / "reference" / "reference.fa"
    fa.write_text(fa.read_text().replace("A", "C", 1))
    with pytest.raises(ReferenceError_):
        _ref(fixture_input)


def test_variant_id_must_match_identity_material():
    with pytest.raises(ValidationError):
        CanonicalVariant(
            variant_id="var_forged", assembly="SYNTH-1", reference_accession="X.1", reference_bundle_id="b",
            chrom="synth1", pos_1based=1, ref="A", alt="G", allele_type="snv",
            normalization_tool="t", normalization_version="1",
        )


def test_resolved_decision_cannot_pick_among_alternatives():
    with pytest.raises(ValidationError):
        IdentityDecision(decision_id="d", mention_id="m", status="resolved", variant_ids=["var_a", "var_b"],
                         mapping_method="x", justification="x")


def test_observed_status_rules():
    with pytest.raises(ValidationError):
        Observed[float](status=ObservationStatus.NOT_SEARCHED, value=0.0)
    with pytest.raises(ValidationError):
        Observed[float](status=ObservationStatus.SEARCHED_NOT_FOUND, source_snapshot_ids=["s"])
    with pytest.raises(ValidationError):
        Observed[float](status=ObservationStatus.FAILED)
    assert Observed[float](status=ObservationStatus.VALUE_PRESENT, value=0.0, source_snapshot_ids=["s"]).value == 0.0


def test_one_ranked_row_per_allele(run_bundle):
    rows = (run_bundle / "ranked_variants.csv").read_text().splitlines()[1:]
    ids = [r.split(",")[6] for r in rows]
    assert len(ids) == len(set(ids))

"""Evidence acceptance: INV-02, INV-03, INV-04, INV-09, INV-12 (ISSUE-001/004)."""

from __future__ import annotations

import json
from collections import defaultdict

import pytest

from tessera.evidence.validation import effective_state
from tessera.pipeline import freeze
from tessera.schemas.base import ValidationState
from tessera.schemas.evidence import ValidationDecision

from conftest import load_features, mention_label


def _states(bundle):
    by = defaultdict(list)
    for d in map(json.loads, (bundle / "validation.jsonl").read_text().splitlines()):
        by[d["link_id"]].append(ValidationDecision.model_validate(d))
    return {k: effective_state(v) for k, v in by.items()}, by


def test_host_rejects_gene_claim_posing_as_exact_allele(run_bundle):
    states, by = _states(run_bundle)
    assert states["l_ko_bad"] is ValidationState.REJECTED  # INV-04
    host = [d for d in by["l_ko_bad"] if d.reviewer_type == "host"][0]
    assert any("gene-scoped" in r for r in host.reasons)


def test_paper_mentioning_allele_outside_the_experiment_is_rejected(run_bundle):
    states, by = _states(run_bundle)
    assert states["l_intro"] is ValidationState.REJECTED  # INV-02
    assert states["l_func1"] is ValidationState.ACCEPTED_HUMAN


def test_wrong_offsets_reject_the_claim(run_bundle):
    states, by = _states(run_bundle)
    assert states["l_obs3"] is ValidationState.REJECTED  # INV-03
    host = [d for d in by["l_obs3"] if d.reviewer_type == "host"][0]
    assert host.checks["citation_location"] == "fail"


def test_human_rejection_overrides_automatic_acceptance(run_bundle):
    states, _ = _states(run_bundle)
    assert states["l_obs4"] is ValidationState.REJECTED


def test_rejected_evidence_contributes_nothing(run_bundle):
    feats = load_features(run_bundle)
    m17 = feats[mention_label(run_bundle, "m17")]  # would be case_controlled if accepted
    m18 = feats[mention_label(run_bundle, "m18")]
    assert m17["o"] == "none" and m18["o"] == "none"  # INV-09
    m04 = feats[mention_label(run_bundle, "m04")]
    assert m04["characterization"] == "uncharacterized"
    m092 = feats[mention_label(run_bundle, "m09.2")]
    assert m092["characterization"] == "uncharacterized" and not m092["contributing_link_ids"]


def test_rejected_claims_stay_in_the_audit_ledger(run_bundle):
    links = {l["link_id"] for l in map(json.loads, (run_bundle / "evidence_links.jsonl").read_text().splitlines())}
    assert {"l_ko_bad", "l_intro", "l_obs3", "l_obs4"} <= links


def test_gene_level_link_is_context_only(run_bundle):
    f = load_features(run_bundle)[mention_label(run_bundle, "m09.1")]
    assert f["characterization"] == "uncharacterized"
    assert f["contributing_link_ids"] == {"context_only": ["l_ko_gene"]}


def test_carrier_comparison_does_not_characterize(run_bundle):
    f = load_features(run_bundle)[mention_label(run_bundle, "m19")]
    assert f["characterization"] == "uncharacterized"
    assert "carrier_or_non_installation_experiment_only" in f["flags"]


def test_observation_tiers(run_bundle):
    feats = load_features(run_bundle)
    assert feats[mention_label(run_bundle, "m01")]["o"] == "case_controlled"
    assert feats[mention_label(run_bundle, "m06")]["o"] == "related_condition"  # never a case tier


def test_altered_snapshot_is_detected(fixture_input, tmp_path):
    snap = fixture_input.parent / "snapshots" / "paper.txt"
    snap.write_text(snap.read_text().replace("decreased", "increased"))
    with pytest.raises(ValueError, match="does not match its digest"):
        freeze(fixture_input, tmp_path / "bad")


def test_synthetic_and_real_sources_cannot_mix(fixture_input, tmp_path):
    p = fixture_input.parent / "snapshots" / "snapshots.jsonl"
    rows = [json.loads(l) for l in p.read_text().splitlines()]
    rows[0]["is_synthetic"] = False
    p.write_text("".join(json.dumps(r) + "\n" for r in rows))
    with pytest.raises(ValueError, match="INV-12"):
        freeze(fixture_input, tmp_path / "mixed")

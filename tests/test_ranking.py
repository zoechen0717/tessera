"""nomination_v0.1: strata, ties, degeneracy, mechanism, channels.
INV-10, INV-17, INV-18, INV-19, DEC-02..07."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest
import yaml

from tessera.pipeline import run
from tessera.ranking.policy import NominationPolicy, PolicyError
from tessera.synthetic import build

from conftest import load_features, load_ranked, mention_label

CONFIG = Path(__file__).resolve().parents[1] / "config"


def _policy_dicts():
    raw = yaml.safe_load((CONFIG / "nomination_v0.1.yaml").read_text())
    table = yaml.safe_load((CONFIG / "mechanism_concordance_v0.1.yaml").read_text())
    return raw, table


def _stratum(bundle, channel, stratum):
    return [r for r in load_ranked(bundle) if r["primary_channel"] == channel and r["stratum"] == stratum]


def test_ranks_never_cross_strata(run_bundle):
    groups = {}
    for r in load_ranked(run_bundle):
        groups.setdefault((r["primary_channel"], r["stratum"]), []).append(int(r["tier_rank"]))
    for ranks in groups.values():  # INV-18: every group restarts at tier 1
        assert min(ranks) == 1 and sorted(set(ranks)) == list(range(1, max(ranks) + 1))
    header = (run_bundle / "ranked_variants.csv").read_text().splitlines()[0]
    assert "score" not in header and "probability" not in header


def test_ties_share_a_tier(run_bundle):
    a, b = mention_label(run_bundle, "m04"), mention_label(run_bundle, "m08")
    rows = {r["label"]: r for r in _stratum(run_bundle, "discovery", "missense")}
    assert rows[a]["tie_group_id"] == rows[b]["tie_group_id"]  # INV-19
    assert rows[a]["tier_rank"] == rows[b]["tier_rank"] and rows[a]["tie_group_size"] == "2"


def test_missing_prediction_sorts_after_scored_not_as_zero(run_bundle):
    rows = _stratum(run_bundle, "discovery", "missense")
    m07 = mention_label(run_bundle, "m07")
    concordant_none = [r for r in rows if r["K"] == "concordant" and r["O"] == "none"]
    assert concordant_none[-1]["label"] == m07
    assert concordant_none[-1]["P"] == "unsupported"


def test_characterized_allele_routes_to_control(run_bundle):
    f = load_features(run_bundle)[mention_label(run_bundle, "m05")]
    assert f["primary_channel"] == "control"
    assert "positive_control_candidate" in f["roles"]


def test_user_control_stays_in_control_channel(run_bundle):
    f = load_features(run_bundle)[mention_label(run_bundle, "m12")]
    assert f["primary_channel"] == "control" and "negative_control" in f["roles"]


def test_keep_in_discovery_policy(fixture_input, tmp_path):
    raw, table = _policy_dicts()
    raw["characterized_policy"] = "keep_in_discovery"
    (tmp_path / "pol").mkdir()
    (tmp_path / "pol" / "p.yaml").write_text(yaml.safe_dump(raw))
    (tmp_path / "pol" / "mechanism_concordance_v0.1.yaml").write_text(yaml.safe_dump(table))
    inp = build(tmp_path / "in2", policy_path=tmp_path / "pol" / "p.yaml")
    out = tmp_path / "run2"
    run(inp, out)
    f = load_features(out)[mention_label(out, "m05")]
    assert f["primary_channel"] == "discovery" and "positive_control_candidate" in f["roles"]


def test_mechanism_inversion_reorders_ptv_stratum(tmp_path):
    """LoF prefers NMD-triggering PTVs; dominant-negative prefers NMD-escaping ones."""
    order = {}
    for mech in ("loss_of_function", "dominant_negative"):
        out = tmp_path / mech
        run(build(tmp_path / f"in_{mech}", mechanism=mech), out)
        rows = [r for r in _stratum(out, "discovery", "ptv") if r["O"] == "none"]
        order[mech] = [r["label"] for r in sorted(rows, key=lambda r: int(r["display_index"]))]
    escaping = mention_label(tmp_path / "loss_of_function", "m02")
    assert order["loss_of_function"][-1] == escaping
    assert order["dominant_negative"][0] == escaping


def test_unknown_mechanism_makes_k_degenerate_and_says_so(tmp_path):
    out = tmp_path / "unk"
    run(build(tmp_path / "in", mechanism="unknown"), out)
    feats = load_features(out)
    assert {f["k"] for f in feats.values()} == {"not_assessable"}
    strata = json.loads((out / "strata_report.json").read_text())
    for s in strata:
        assert next(k for k in s["keys"] if k["key"] == "K")["degenerate"]
    assert "key K is constant" in (out / "report.md").read_text()


def test_mechanism_never_becomes_allele_evidence(run_bundle):
    for f in load_features(run_bundle).values():  # INV-17
        assert "mechanism" not in json.dumps(f["contributing_link_ids"])


def test_duplicate_evidence_does_not_change_tiers(fixture_input, tmp_path, run_bundle):
    base = fixture_input.parent
    links = [json.loads(l) for l in (base / "links.jsonl").read_text().splitlines()]
    decisions = [json.loads(l) for l in (base / "decisions.jsonl").read_text().splitlines()]
    dup = copy.deepcopy(next(l for l in links if l["link_id"] == "l_obs1"))
    dup["link_id"] = "l_obs1_dup"
    d = copy.deepcopy(next(x for x in decisions if x["link_id"] == "l_obs1"))
    d.update(link_id="l_obs1_dup", decision_id="vd_dup", sequence=99)
    (base / "links.jsonl").write_text("".join(json.dumps(x) + "\n" for x in links + [dup]))
    (base / "decisions.jsonl").write_text("".join(json.dumps(x) + "\n" for x in decisions + [d]))
    out = tmp_path / "dup"
    run(fixture_input, out)
    before = {r["variant_id"]: (r["tier_rank"], r["tie_group_size"]) for r in load_ranked(run_bundle)}
    after = {r["variant_id"]: (r["tier_rank"], r["tie_group_size"]) for r in load_ranked(out)}
    assert before == after  # INV-10


def test_features_do_not_depend_on_other_candidates(fixture_input, tmp_path, run_bundle):
    p = fixture_input.parent / "variants.tsv"
    p.write_text("".join(l for l in p.read_text().splitlines(keepends=True) if not l.startswith("m10\t")))
    out = tmp_path / "fewer"
    run(fixture_input, out)
    before, after = load_features(run_bundle), load_features(out)
    for label, f in after.items():
        strip = lambda x: {k: v for k, v in x.items() if k != "evidence_digest"}  # noqa: E731
        assert strip(f) == strip(before[label])


def test_policy_rejects_avi_with_alphamissense():
    raw, table = _policy_dicts()
    raw["predictors"]["splice"] = "avi"
    with pytest.raises(PolicyError):
        NominationPolicy.from_dicts(raw, table)


def test_policy_rejects_unknown_mechanism_with_concordance():
    raw, table = _policy_dicts()
    table["stratum_level"]["unknown"]["ptv"] = "concordant"
    with pytest.raises(PolicyError):
        NominationPolicy.from_dicts(raw, table)


def test_report_states_synthetic_and_provisional(run_bundle):
    text = (run_bundle / "report.md").read_text()
    assert "SYNTHETIC FIXTURE RUN" in text and "provisional" in text and "not_assessed" in text

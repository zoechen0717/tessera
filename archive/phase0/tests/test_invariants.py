"""Structural invariants (spec 49.1).

These are properties of the code, not statistics about its output. Each one
corresponds to a defect the architecture review found, and each must FAIL LOUDLY
rather than be a convention someone remembers. A change that breaks one of these
is an architectural regression, not a tuning difference, and must not be resolved
by adjusting a threshold.
"""

from __future__ import annotations

import ast
import math
from pathlib import Path

import pytest

from conftest import make_annotations, make_evidence, make_variant
from tessera.config.anchors import AnchorTable, MissingAnchorError
from tessera.config.scoring_config import ScoringConfig, ScoringConfigError
from tessera.enums import (
    ConsequenceClass,
    DiseaseMatchType,
    EvidenceLevel,
    Mechanism,
    MechanismConfidence,
    Missingness,
    RankingMode,
    ScoreComponent,
)
from tessera.ranking import components
from tessera.ranking.ranker import diagnose, rank_stratified, score_variant
from tessera.schemas.common import Quantity
from tessera.schemas.gene import GeneMechanism

SRC = Path(__file__).resolve().parents[1] / "src" / "tessera"


def _mech(mechanism=Mechanism.HAPLOINSUFFICIENCY, confidence=MechanismConfidence.HIGH):
    return GeneMechanism(gene_symbol="SETD1A", mechanism=mechanism, confidence=confidence)


def _score(cfg, *, mode=RankingMode.DISCOVERY, variant=None, annotations=None, evidence=None,
           mechanism=None):
    variant = variant or make_variant()
    annotations = annotations if annotations is not None else make_annotations()
    evidence = evidence if evidence is not None else [make_evidence()]
    return score_variant(cfg, mode, variant, annotations, evidence, mechanism or _mech())


# ---------------------------------------------------------------- spec 19.5


def test_novelty_multiplier_is_independent_of_disease_and_mechanism(cfg):
    """The exact defect in the rejected formulation.

    In the linear form, novelty was `(0.45*D + 0.35*M + ...) * (1 - saturation)`,
    so the novelty term moved when D or M moved -- which reintroduced D and M into
    the composite a second time. The multiplier must depend ONLY on saturation.
    """
    low = _score(cfg, annotations=make_annotations(alphamissense=0.05, clinvar_class="benign",
                                                   stars=0, af=0.05))
    high = _score(cfg, annotations=make_annotations(alphamissense=0.99, clinvar_class="pathogenic",
                                                    stars=4, af=1e-7))

    assert low.disease.value != pytest.approx(high.disease.value)
    assert low.mechanism.value != pytest.approx(high.mechanism.value)
    assert low.novelty_multiplier == pytest.approx(high.novelty_multiplier), (
        "the novelty multiplier moved when D and M moved -- the rejected linear "
        "formulation has crept back in (spec 19.5)"
    )


def test_configured_weights_are_effective_weights(cfg):
    """Perturbing one component must move the score by exactly its weight.

    Measured on the weighted sum, because the novelty multiplier is a common
    factor over all components and so does not distort their RATIOS.
    """
    mode = cfg.mode(RankingMode.DISCOVERY)

    def sum_with(am: float) -> tuple[float, float]:
        s = _score(cfg, annotations=make_annotations(alphamissense=am))
        return s.weighted_sum, s.mechanism.value

    s_lo, m_lo = sum_with(0.10)
    s_hi, m_hi = sum_with(0.90)

    observed = (s_hi - s_lo) / (m_hi - m_lo)
    assert observed == pytest.approx(mode.weights.mechanism, abs=1e-9), (
        f"sensitivity to M is {observed:.4f} but the configured weight is "
        f"{mode.weights.mechanism}. A component is entering the composite more "
        f"than once (spec 19.5)."
    )


def test_functional_evidence_raises_the_score_unambiguously(cfg):
    """In the rejected form, F's sign depended on that variant's D and M values.

    Adding direct supporting functional evidence must never lower the score,
    whatever D and M happen to be.
    """
    for am, cls_, stars in ((0.05, "benign", 0), (0.5, "vus", 1), (0.99, "pathogenic", 4)):
        ann = make_annotations(alphamissense=am, clinvar_class=cls_, stars=stars)
        none_found = _score(cfg, annotations=ann, evidence=[])
        with_ev = _score(cfg, annotations=ann, evidence=[make_evidence()])
        assert with_ev.functional.value >= (none_found.functional.value or 0.0), (
            f"direct functional evidence lowered F at am={am} -- F's sign is "
            f"context dependent (spec 19.5)"
        )


def test_saturation_is_not_one_minus_functional(cfg):
    """"Tested, null result" must be distinguishable from "never tested".

    A variant with strong direct evidence of NO effect has HIGH saturation and LOW
    F. Deriving saturation from `1 - F` would invert that and re-promote it as a
    discovery candidate.
    """
    contradicting = make_evidence(evidence_id="ev_null")
    contradicting.direction = contradicting.direction.__class__("contradicts_priority")
    tested_null = _score(cfg, evidence=[contradicting])
    never_tested = _score(cfg, evidence=[])

    assert tested_null.functional.value is not None
    assert tested_null.functional.value <= 0.0, "a contradicting result should not raise F"
    assert tested_null.saturation.value is not None and tested_null.saturation.value > 0.5, (
        "a variant that was directly tested must have HIGH saturation even though F "
        "is low (spec 19.5)"
    )
    assert never_tested.saturation.value is None, "never-tested saturation must be unknown"
    assert tested_null.novelty_multiplier < 1.0 + cfg.mode(RankingMode.DISCOVERY).novelty_alpha


# ---------------------------------------------------------------- spec 19.6


def test_scores_are_independent_of_candidate_set(cfg):
    """Guarantee 9. Enforced structurally: score_variant never sees the set."""
    v = make_variant()
    ann = make_annotations()
    alone = score_variant(cfg, RankingMode.DISCOVERY, v, ann, [make_evidence()], _mech())

    others = [
        score_variant(
            cfg,
            RankingMode.DISCOVERY,
            make_variant(vrs_id=f"ga4gh:VA.{i}"),
            make_annotations(vrs_id=f"ga4gh:VA.{i}", alphamissense=0.99, af=1e-9),
            [],
            _mech(),
        )
        for i in range(2, 12)
    ]
    again = score_variant(cfg, RankingMode.DISCOVERY, v, ann, [make_evidence()], _mech())

    assert alone.final_score == again.final_score == pytest.approx(alone.final_score)
    assert all(o.vrs_id != alone.vrs_id for o in others)


def test_every_raw_input_must_have_an_anchor(cfg, tmp_path):
    """A missing anchor raises at config load. It never defaults."""
    with pytest.raises(MissingAnchorError):
        cfg.anchors.require(["a_predictor_nobody_anchored"])

    bad = tmp_path / "anchors.yaml"
    bad.write_text('version: "0.1"\nanchors:\n  gnomad_af:\n    type: identity\n')
    table = AnchorTable.load(bad)
    with pytest.raises(MissingAnchorError):
        table.apply("alphamissense_score", 0.5)


def test_no_candidate_set_normalization_in_source():
    """No component may min-max, z-score, quantile, or rank-normalize.

    Checked by source inspection because this is exactly the kind of convenience
    that gets reintroduced by someone reaching for a quick normalization.
    """
    banned = ("minmax", "min_max", "zscore", "z_score", "quantile", "rankdata", "percentile")
    for path in (SRC / "ranking").rglob("*.py"):
        text = path.read_text().lower()
        for token in banned:
            assert token not in text, (
                f"{path.name} mentions {token!r}: candidate-set-dependent "
                f"normalization is prohibited (spec 19.6, Guarantee 9)"
            )


def test_anchor_mappings_are_monotone_where_declared(cfg):
    """A piecewise anchor must not accidentally invert partway through."""
    for name, anchor in cfg.anchors.anchors.items():
        if anchor.type.value != "piecewise":
            continue
        xs = anchor.breakpoints or []
        ys = [anchor.apply(x) for x in xs]
        assert ys == sorted(ys) or ys == sorted(ys, reverse=True), (
            f"anchor {name!r} is not monotone across its breakpoints: {ys}"
        )


# ---------------------------------------------------------------- spec 42.1


def test_not_searched_is_not_scored_as_zero(cfg):
    """The silent failure that fills a panel with unresearched variants.

    Because the novelty multiplier rewards low saturation, a variant that was
    merely skipped would otherwise be promoted as a discovery candidate.
    """
    searched = make_annotations()
    searched.protein_effect.alphamissense_score = Quantity.searched_not_found()
    unsearched = make_annotations()
    unsearched.protein_effect.alphamissense_score = Quantity.not_searched()

    s_searched = _score(cfg, annotations=searched)
    s_unsearched = _score(cfg, annotations=unsearched)

    assert s_searched.mechanism.value != pytest.approx(s_unsearched.mechanism.value), (
        "searched_not_found and not_searched produced the same M -- they are being "
        "collapsed (spec 42.1)"
    )
    assert s_unsearched.mechanism.value > s_searched.mechanism.value, (
        "a missing input must be excluded with weights renormalized, not "
        "substituted with 0.0"
    )
    assert s_unsearched.mechanism.weights_renormalized


def test_quantity_refuses_a_value_it_should_not_have():
    with pytest.raises(ValueError):
        Quantity(value=0.5, state=Missingness.NOT_SEARCHED)


def test_failed_retrieval_blocks_panel_inclusion(cfg):
    """A transient API failure must not quietly become a biological conclusion."""
    ann = make_annotations()
    ann.splice.spliceai_max_delta = Quantity.failed(source="spliceai")
    s = _score(cfg, annotations=ann)
    assert not s.panel_eligible
    assert any("retrieval_failure" in f for f in s.blocking_flags)


def test_never_searched_gets_no_novelty_bonus(cfg):
    s = _score(cfg, evidence=[])
    assert s.saturation.value is None
    assert s.novelty_multiplier == pytest.approx(1.0), (
        "an unresearched variant received the novelty bonus (spec 42.1)"
    )
    assert "incomplete_evidence" in s.flags


# ---------------------------------------------------------------- spec 43.1


def test_scoring_never_reads_agent_confidence():
    """Enforced by AST inspection, not by code review.

    An LLM's self-rated certainty is not calibrated and its scale drifts across
    models and prompt versions; letting it into a score would make the ranking
    silently model dependent.
    """
    forbidden = {"agent_extraction_confidence", "agent_confidence"}
    for path in (SRC / "ranking").rglob("*.py"):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute) and node.attr in forbidden:
                pytest.fail(f"{path.name} reads {node.attr} (spec 43.1)")
            if isinstance(node, ast.Constant) and node.value in forbidden:
                pytest.fail(f"{path.name} references {node.value!r} (spec 43.1)")


# ---------------------------------------------------------------- spec 19.1


def test_recurrence_is_not_a_count_of_reports(cfg):
    """Spec 19.1.1. Same proband, four sources -> one observation.

    Otherwise ClinVar submissions derived from published case reports get counted
    again by the literature route, inflating D for exactly the variants most
    likely to be already known.
    """
    one_proband_four_sources = [
        make_evidence(evidence_id="ev_pub", source_type="publication", source_id="PMID:9",
                      proband="P1"),
        make_evidence(evidence_id="ev_cv", source_type="clinvar", source_id="PMID:9",
                      proband="P1"),
        make_evidence(evidence_id="ev_rev", source_type="review", source_id="PMID:9",
                      proband="P1"),
        make_evidence(evidence_id="ev_meta", source_type="publication", source_id="PMID:9",
                      proband="P1"),
    ]
    four_probands = [
        make_evidence(evidence_id=f"ev_{i}", source_id=f"PMID:{i}", proband=f"P{i}")
        for i in range(1, 5)
    ]
    r_dup = components._cohort_normalized_recurrence(one_proband_four_sources)
    r_real = components._cohort_normalized_recurrence(four_probands)
    assert r_dup is not None and r_real is not None
    assert r_real > r_dup, "four sources for one proband scored as much as four probands"


def test_recurrence_requires_a_denominator(cfg):
    """An unnormalizable count must not become a score.

    Without the cohort-size and mutational-target correction the statistic
    collapses into a gene-length artifact.
    """
    no_cohort = [make_evidence(cohort_size=None)]
    assert components._cohort_normalized_recurrence(no_cohort) is None


def test_related_disease_does_not_score_as_exact(cfg):
    """Spec 12.1.1. The most likely way for this system to be confidently wrong.

    A gene may have abundant, high-quality clinical records whose condition is a
    DIFFERENT disease from the query.
    """
    exact = _score(cfg, evidence=[make_evidence(match=DiseaseMatchType.EXACT)])
    broader = _score(cfg, evidence=[make_evidence(match=DiseaseMatchType.BROADER)])
    unrelated = _score(cfg, evidence=[make_evidence(match=DiseaseMatchType.UNRELATED)])

    assert exact.disease.value > broader.disease.value > 0
    assert unrelated.functional.value == pytest.approx(0.0), (
        "evidence for an unrelated condition entered F (spec 12.1.1)"
    )


# ------------------------------------------------------------ spec 19.3, 9.8


def test_mechanism_gating_changes_the_ordering(cfg):
    """If flipping the gene mechanism reorders nothing, the gate is not wired in."""
    lof_like = make_variant(vrs_id="ga4gh:VA.lof")
    lof_ann = make_annotations(vrs_id="ga4gh:VA.lof", alphamissense=0.95)
    ptv = make_variant(vrs_id="ga4gh:VA.ptv", consequence_class=ConsequenceClass.PTV)
    ptv_ann = make_annotations(vrs_id="ga4gh:VA.ptv")

    def order(mechanism):
        scores = [
            score_variant(cfg, RankingMode.DISCOVERY, v, a, [], _mech(mechanism))
            for v, a in ((lof_like, lof_ann), (ptv, ptv_ann))
        ]
        return {s.vrs_id: s.mechanism.value for s in scores}

    hi = order(Mechanism.HAPLOINSUFFICIENCY)
    gof = order(Mechanism.GAIN_OF_FUNCTION)
    assert hi["ga4gh:VA.ptv"] > gof["ga4gh:VA.ptv"], (
        "a PTV scored the same under haploinsufficiency and gain-of-function -- "
        "mechanism gating is not applied (spec 9.8)"
    )


def test_unknown_mechanism_does_not_gate(cfg):
    gated = _score(cfg, mechanism=_mech(Mechanism.HAPLOINSUFFICIENCY))
    ungated = _score(cfg, mechanism=_mech(Mechanism.UNKNOWN, MechanismConfidence.LOW))
    assert ungated.mechanism.gate_multiplier == 1.0
    assert gated.mechanism.gate_multiplier < 1.0
    assert "mechanism:mechanism_unknown_no_gating" in ungated.flags


def test_low_confidence_mechanism_is_disclosed(cfg):
    s = _score(cfg, mechanism=_mech(Mechanism.HAPLOINSUFFICIENCY, MechanismConfidence.LOW))
    assert s.mechanism_is_assumed
    assert "mechanism:mechanism_assumed_low_confidence" in s.flags


def test_each_consequence_class_uses_its_own_branch(cfg):
    """M_coding and M_splice come from different predictors; no class may borrow
    another's formula."""
    seen = {}
    for cls_ in (ConsequenceClass.MISSENSE, ConsequenceClass.SPLICE, ConsequenceClass.PTV,
                 ConsequenceClass.SYNONYMOUS, ConsequenceClass.NONCODING):
        v = make_variant(vrs_id=f"ga4gh:VA.{cls_.value}", consequence_class=cls_)
        a = make_annotations(vrs_id=f"ga4gh:VA.{cls_.value}")
        s = score_variant(cfg, RankingMode.DISCOVERY, v, a, [], _mech())
        seen[cls_] = s.mechanism.branch
    assert seen[ConsequenceClass.MISSENSE] == "M_missense"
    assert seen[ConsequenceClass.SPLICE] == "M_splice"
    assert len(set(seen.values())) == len(seen)


def test_ranking_is_stratified_by_consequence_class(cfg):
    scores = []
    for i, cls_ in enumerate((ConsequenceClass.MISSENSE, ConsequenceClass.SPLICE,
                              ConsequenceClass.MISSENSE)):
        v = make_variant(vrs_id=f"ga4gh:VA.{i}", consequence_class=cls_)
        a = make_annotations(vrs_id=f"ga4gh:VA.{i}")
        scores.append(score_variant(cfg, RankingMode.DISCOVERY, v, a, [], _mech()))
    strata = rank_stratified(scores)
    assert set(strata) == {ConsequenceClass.MISSENSE, ConsequenceClass.SPLICE}
    assert len(strata[ConsequenceClass.MISSENSE]) == 2


# ---------------------------------------------------------------- spec 17.1.2


def test_host_line_conflict_excludes_from_panel(cfg):
    """A hard filter, not a score penalty. Nothing in a purely biological ranking
    would surface this, and finding it after a panel is ordered is expensive."""
    ann = make_annotations()
    ann.editability.host_line_conflict = "heterozygous"
    ann.editability.variant_status = ann.editability.resolve_status()
    s = _score(cfg, annotations=ann)
    assert not s.panel_eligible
    assert any("not_editable_in_host_line" in f for f in s.blocking_flags)


def test_ungenotyped_host_line_is_not_a_conflict(cfg):
    """editability_unknown must not be scored as not_editable."""
    ann = make_annotations()
    ann.editability.host_line_genotype_verified = False
    s = _score(cfg, annotations=ann)
    assert s.panel_eligible, "an ungenotyped position was treated as a conflict"


# ---------------------------------------------------------- spec 48.9, 53.6


def test_degenerate_component_is_reported(cfg):
    """The SETD1A + schizophrenia case.

    Disease evidence is real but gene-level, so D is constant across candidates
    and contributes nothing to the ordering -- while a configured w_D of 0.32
    would otherwise imply it did.
    """
    scores = []
    for i in range(6):
        v = make_variant(vrs_id=f"ga4gh:VA.{i}")
        a = make_annotations(vrs_id=f"ga4gh:VA.{i}", alphamissense=0.2 + 0.1 * i)
        scores.append(score_variant(cfg, RankingMode.DISCOVERY, v, a, [], _mech()))

    diag = diagnose(cfg, RankingMode.DISCOVERY, scores)
    assert ScoreComponent.DISEASE in diag.degenerate_components
    assert any("contributed nothing to the ordering" in w for w in diag.warnings)
    disease_row = next(r for r in diag.discrimination if r.component is ScoreComponent.DISEASE)
    assert disease_row.configured_weight > 0
    assert disease_row.effective_contribution == pytest.approx(0.0, abs=1e-6)
    mech_row = next(r for r in diag.discrimination if r.component is ScoreComponent.MECHANISM)
    assert not mech_row.degenerate


def test_final_score_always_reports_its_components(cfg):
    """Spec 21: no opaque single number."""
    s = _score(cfg)
    for comp in ScoreComponent:
        assert s.component(comp) is not None
        assert s.component(comp).terms, f"{comp.value} has no term breakdown"


# ---------------------------------------------------------------- config guards


def test_weights_must_sum_to_one(cfg):
    """Pydantic wraps the ScoringConfigError, so match on the message: the point is
    that the guard fires and says why, not which exception class surfaces."""
    from pydantic import ValidationError

    from tessera.config.scoring_config import ModeWeights

    with pytest.raises(ValidationError, match="must sum to 1.0"):
        ModeWeights(disease=0.5, functional=0.5, mechanism=0.5, editability=0.5)


def test_novelty_is_not_a_weighted_component(cfg):
    """Novelty must not appear among the weights at all."""
    for mode_cfg in cfg.modes.values():
        assert set(mode_cfg.weights.as_dict()) == set(ScoreComponent)
        assert "novelty" not in {c.value for c in mode_cfg.weights.as_dict()}
        assert math.isclose(sum(mode_cfg.weights.as_dict().values()), 1.0, abs_tol=1e-9)


def test_functional_aggregation_cannot_become_a_sum(cfg):
    """Summing evidence tiers is how paper-counting re-enters after being excluded
    from recurrence (spec 19.2)."""
    from pydantic import ValidationError

    from tessera.config.scoring_config import FunctionalConfig

    with pytest.raises(ValidationError, match="must be 'max'"):
        FunctionalConfig(aggregation="sum")


def test_no_branch_double_counts_alphamissense(cfg):
    """AVI is a supervised model that TAKES AlphaMissense as an input feature.

    So raw `alphagenome_avi` beside `alphamissense_score` in the same branch counts
    AlphaMissense twice. Coding branches must use `avi_regulatory_shap` -- the
    SHAP-isolated AlphaGenome-modality contribution -- instead.

    This is the same class of error the architecture review flagged elsewhere, and
    it was present in the first version of mechanism_branches.yaml. Hence a test.
    """
    for cls_, weights in cfg.branches.branches.items():
        if "alphamissense_score" in weights:
            assert "alphagenome_avi" not in weights, (
                f"branch {cls_.value!r} contains both alphamissense_score and raw "
                f"alphagenome_avi. AVI already contains AlphaMissense as an input "
                f"feature -- use avi_regulatory_shap (spec 9.5.1)."
            )


def test_avi_anchor_is_not_identity(cfg):
    """AVI is PHRED-scaled, not [0,1].

    Anchoring it as `identity` silently truncates every real AVI value to 1.0,
    which looks like a working pipeline producing a constant.
    """
    anchor = cfg.anchors["alphagenome_avi"]
    assert anchor.type.value != "identity", (
        "alphagenome_avi is anchored as identity, but AVI is PHRED-scaled "
        "(10 = top 10%, 20 = top 1%, 30 = top 0.1%)"
    )
    assert cfg.anchors.apply("alphagenome_avi", 20.0) < 1.0
    assert cfg.anchors.apply("alphagenome_avi", 5.0) < cfg.anchors.apply(
        "alphagenome_avi", 25.0
    )


def test_mechanism_branch_weights_sum_to_one(cfg):
    for cls_, weights in cfg.branches.branches.items():
        assert math.isclose(sum(weights.values()), 1.0, abs_tol=1e-9), cls_


def test_validation_mode_has_no_novelty_bonus(cfg):
    assert cfg.mode(RankingMode.VALIDATION).novelty_alpha == 0.0
    s = _score(cfg, mode=RankingMode.VALIDATION)
    assert s.novelty_multiplier == pytest.approx(1.0)

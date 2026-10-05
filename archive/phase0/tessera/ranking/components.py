"""Component scores D, F, M, E.

Two properties hold throughout this module and are tested in
tests/test_invariants.py:

  1. Every function here is PURE PER VARIANT. Nothing reads the candidate set.
     That is what makes Guarantee 9 true: adding a variant to the run cannot
     change another variant's score.

  2. Nothing here reads agent_extraction_confidence (spec 43.1).
"""

from __future__ import annotations

from tessera.config.scoring_config import ScoringConfig
from tessera.enums import (
    ConsequenceClass,
    EvidenceLevel,
    Mechanism,
    Missingness,
    ScoreComponent,
)
from tessera.schemas.annotations import VariantAnnotations
from tessera.schemas.common import Quantity
from tessera.schemas.evidence import Evidence, distinct_observations
from tessera.schemas.gene import GeneMechanism
from tessera.schemas.scoring import ComponentScore, SaturationResult, TermContribution
from tessera.schemas.variant import CanonicalVariant


def _weighted(
    cfg: ScoringConfig,
    component: ScoreComponent,
    weights: dict[str, float],
    annotations: VariantAnnotations,
    *,
    branch: str | None = None,
) -> ComponentScore:
    """Anchor each raw input, then take a weighted mean over the PRESENT ones.

    Missing inputs are dropped and the remaining weights renormalize. The
    alternative -- substituting 0.0 for a missing input -- would make an
    unmeasured variant look measured-and-bad, which is the distinction spec 42.1
    exists to protect.
    """
    terms: list[TermContribution] = []
    present: list[tuple[float, float]] = []
    missing = 0

    for raw_name, weight in weights.items():
        q: Quantity = annotations.raw_input(raw_name)
        if not q.enters_weighted_sum:
            missing += 1
            terms.append(
                TermContribution(
                    raw_input=raw_name, weight=weight, state=q.state.value, contribution=0.0
                )
            )
            continue
        if q.is_present:
            raw_value = q.anchor_input
            anchored = cfg.anchors.apply(raw_name, raw_value)
        else:
            # searched_not_found: contributes a real 0.0, without being anchored --
            # there is no raw value to map.
            raw_value, anchored = None, 0.0
        present.append((anchored, weight))
        terms.append(
            TermContribution(
                raw_input=raw_name,
                raw_value=raw_value if isinstance(raw_value, (int, float, str)) else str(raw_value),
                anchored=anchored,
                weight=weight,
                state=q.state.value,
            )
        )

    total_weight = sum(w for _, w in present)
    if total_weight <= 0.0:
        return ComponentScore(
            component=component,
            value=None,
            terms=terms,
            branch=branch,
            inputs_present=0,
            inputs_missing=missing,
            flags=["uncomputable_no_inputs"],
        )

    value = sum(a * w for a, w in present) / total_weight
    renormalized = abs(total_weight - 1.0) > 1e-9

    by_name = {t.raw_input: t for t in terms}
    for raw_name, weight in weights.items():
        t = by_name[raw_name]
        if t.anchored is not None:
            t.weight_after_renormalization = weight / total_weight
            t.contribution = t.anchored * t.weight_after_renormalization

    return ComponentScore(
        component=component,
        value=value,
        terms=terms,
        branch=branch,
        weights_renormalized=renormalized,
        inputs_present=len(present),
        inputs_missing=missing,
        flags=["weights_renormalized"] if renormalized else [],
    )


# --------------------------------------------------------------- disease (D)


def disease_score(
    cfg: ScoringConfig,
    variant: CanonicalVariant,
    annotations: VariantAnnotations,
    evidence: list[Evidence],
) -> ComponentScore:
    """Spec 19.1.

    Clinical evidence is multiplied by the disease-match anchor: a well-classified
    ClinVar record for a DIFFERENT condition is real evidence, but not evidence for
    the query disease. Recurrence comes from deduplicated observations normalized
    by cohort size and mutational target -- never from a count of reports.
    """
    sub = cfg.disease
    terms: list[TermContribution] = []
    present: list[tuple[float, float]] = []
    missing = 0

    def add(name: str, weight: float, anchored: float | None, raw: object, state: str) -> None:
        nonlocal missing
        if anchored is None:
            missing += 1
            terms.append(
                TermContribution(raw_input=name, weight=weight, state=state, contribution=0.0)
            )
            return
        present.append((anchored, weight))
        terms.append(
            TermContribution(
                raw_input=name,
                raw_value=raw if isinstance(raw, (int, float, str)) else None,
                anchored=anchored,
                weight=weight,
                state=state,
            )
        )

    # 1. fine mapping -- PIP preferred over raw p-value
    pip = annotations.gwas.posterior_inclusion_probability
    logp = annotations.gwas.log10_p
    if pip.enters_weighted_sum:
        add(
            "gwas_pip",
            sub.fine_mapping,
            cfg.anchors.apply("gwas_pip", pip.anchor_input) if pip.is_present else 0.0,
            pip.value,
            pip.state.value,
        )
    elif logp.enters_weighted_sum:
        add(
            "gwas_log10_p",
            sub.fine_mapping,
            cfg.anchors.apply("gwas_log10_p", logp.anchor_input) if logp.is_present else 0.0,
            logp.value,
            logp.state.value,
        )
    else:
        add("gwas_pip", sub.fine_mapping, None, None, pip.state.value)

    # 2. clinical evidence, scaled by how well the source condition matches
    cls_ = annotations.clinvar.germline_classification
    stars = annotations.clinvar.review_stars
    if cls_ and stars.enters_weighted_sum:
        cls_anchored = cfg.anchors.apply("clinvar_classification", cls_)
        star_anchored = cfg.anchors.apply("clinvar_review_stars", stars.anchor_input)
        match = _best_disease_match(cfg, evidence)
        add(
            "clinvar_classification",
            sub.clinical_evidence,
            cls_anchored * star_anchored * match,
            cls_,
            "present",
        )
    else:
        add("clinvar_classification", sub.clinical_evidence, None, None, stars.state.value)

    # 3. recurrence -- enrichment ratio over deduplicated observations
    recurrence = _cohort_normalized_recurrence(evidence)
    if recurrence is None:
        add("cohort_normalized_recurrence", sub.recurrence, None, None, "not_searched")
    else:
        add(
            "cohort_normalized_recurrence",
            sub.recurrence,
            cfg.anchors.apply("cohort_normalized_recurrence", recurrence),
            recurrence,
            "present",
        )

    # 4. rarity
    af = annotations.frequency.global_af
    if af.enters_weighted_sum:
        add(
            "gnomad_af",
            sub.rarity,
            cfg.anchors.apply("gnomad_af", af.anchor_input) if af.is_present else 0.0,
            af.value,
            af.state.value,
        )
    else:
        add("gnomad_af", sub.rarity, None, None, af.state.value)

    # 5. segregation / de novo
    seg = _segregation_key(evidence)
    add(
        "segregation_support",
        sub.segregation,
        cfg.anchors.apply("segregation_support", seg),
        seg,
        "present" if seg != "none" else "searched_not_found",
    )

    total_weight = sum(w for _, w in present)
    if total_weight <= 0.0:
        return ComponentScore(
            component=ScoreComponent.DISEASE,
            value=None,
            terms=terms,
            inputs_missing=missing,
            flags=["uncomputable_no_inputs"],
        )

    value = sum(a * w for a, w in present) / total_weight
    by_name = {t.raw_input: t for t in terms}
    for t in by_name.values():
        if t.anchored is not None:
            t.weight_after_renormalization = t.weight / total_weight
            t.contribution = t.anchored * t.weight_after_renormalization

    return ComponentScore(
        component=ScoreComponent.DISEASE,
        value=value,
        terms=terms,
        weights_renormalized=abs(total_weight - 1.0) > 1e-9,
        inputs_present=len(present),
        inputs_missing=missing,
    )


def _best_disease_match(cfg: ScoringConfig, evidence: list[Evidence]) -> float:
    """The most favourable disease-match multiplier among usable evidence.

    Defaults to the `related` anchor rather than to 1.0: assuming an exact match
    in the absence of evidence is precisely the error spec 12.1.1 warns about.
    """
    best = cfg.anchors.apply("disease_match_type", "related")
    for ev in evidence:
        if not ev.usable_in_scoring:
            continue
        best = max(best, cfg.anchors.apply("disease_match_type", ev.disease_match_type.value))
    return best


def _cohort_normalized_recurrence(evidence: list[Evidence]) -> float | None:
    """observed independent probands / expected under null (spec 19.1.1).

    Returns None when no cohort size is available: an unnormalizable count must
    not be turned into a score, because without the mutational-target and
    cohort-size correction the statistic collapses into a gene-length artifact.
    """
    usable = [
        ev
        for ev in evidence
        if ev.usable_in_scoring and ev.evidence_level.is_variant_specific() and ev.observation
    ]
    if not usable:
        return None
    groups = distinct_observations(usable)
    observed = sum(
        1
        for group in groups.values()
        if any(e.observation and e.observation.proband_resolution == "proband_level" for e in group)
    )
    if observed == 0:
        return None
    cohorts = {ev.observation.cohort_id: ev.cohort_size for ev in usable if ev.observation}
    total_cohort = sum(size for size in cohorts.values() if size)
    if not total_cohort:
        return None
    # Placeholder null model. Phase 2 replaces this with a per-gene, per-class
    # mutation-rate model; the interface is deliberately the ratio, not the count,
    # so that swapping the null model does not change the anchor.
    expected = max(total_cohort * 1e-4, 1e-6)
    return observed / expected


def _segregation_key(evidence: list[Evidence]) -> str:
    ranked = [
        "de_novo_confirmed",
        "segregates_in_family",
        "de_novo_unconfirmed",
        "case_control_only",
        "single_case_report",
    ]
    found = {
        ev.evidence_type
        for ev in evidence
        if ev.usable_in_scoring and ev.evidence_level.is_variant_specific()
    }
    for key in ranked:
        if key in found:
            return key
    return "none"


# ------------------------------------------------------------ functional (F)


def functional_score(
    cfg: ScoringConfig,
    annotations: VariantAnnotations,
    evidence: list[Evidence],
) -> ComponentScore:
    """Spec 19.2. MAX over tiers, not a sum.

    Three papers reporting the same assay are not three times the evidence, and
    summing tiers is how paper-counting re-enters after being excluded from
    recurrence. Replication earns a bounded additive bonus, counted on
    deduplicated observations.
    """
    terms: list[TermContribution] = []
    candidates: list[float] = []

    # A MaveDB measurement of the exact variant is variant_direct evidence.
    mave = annotations.mave.functional_class
    if mave.enters_weighted_sum and mave.is_present:
        anchored = cfg.anchors.apply("mave_functional_class", mave.anchor_input)
        candidates.append(anchored)
        terms.append(
            TermContribution(
                raw_input="mave_functional_class",
                raw_value=str(mave.anchor_input),
                anchored=anchored,
                weight=1.0,
                state="present",
            )
        )

    usable = [ev for ev in evidence if ev.usable_in_scoring]
    for ev in usable:
        if not ev.directness_key:
            continue
        directness = cfg.anchors.apply("evidence_directness", ev.directness_key)
        relevance = (
            cfg.anchors.apply("model_system_relevance", ev.model_system_key)
            if ev.model_system_key
            else 1.0
        )
        match = cfg.anchors.apply("disease_match_type", ev.disease_match_type.value)
        # Contradicting evidence is not discarded: it lowers F rather than being
        # averaged into agreement (spec 41).
        sign = -1.0 if ev.direction.value == "contradicts_priority" else 1.0
        tier = directness * relevance * match * sign
        candidates.append(tier)
        terms.append(
            TermContribution(
                raw_input=f"evidence:{ev.evidence_id}",
                raw_value=ev.directness_key,
                anchored=tier,
                weight=1.0,
                state="present",
            )
        )

    if not candidates:
        searched = any(ev for ev in evidence)
        state = Missingness.SEARCHED_NOT_FOUND if searched else Missingness.NOT_SEARCHED
        if state is Missingness.SEARCHED_NOT_FOUND:
            return ComponentScore(
                component=ScoreComponent.FUNCTIONAL,
                value=0.0,
                terms=terms,
                inputs_present=0,
                inputs_missing=0,
                flags=["searched_not_found"],
            )
        return ComponentScore(
            component=ScoreComponent.FUNCTIONAL,
            value=None,
            terms=terms,
            inputs_missing=1,
            flags=["not_searched"],
        )

    base = max(candidates)
    n_independent = len(distinct_observations([e for e in usable if e.observation]))
    bonus = 0.0
    if n_independent > 1:
        rep = cfg.anchors.apply("replication_count", n_independent)
        bonus = cfg.functional.replication_bonus * rep
        terms.append(
            TermContribution(
                raw_input="replication_count",
                raw_value=n_independent,
                anchored=rep,
                weight=cfg.functional.replication_bonus,
                contribution=bonus,
                state="present",
            )
        )

    return ComponentScore(
        component=ScoreComponent.FUNCTIONAL,
        value=min(1.0, max(0.0, base + bonus)),
        terms=terms,
        inputs_present=len(candidates),
        branch="max",
    )


# ------------------------------------------------------------- mechanism (M)


def mechanism_score(
    cfg: ScoringConfig,
    variant: CanonicalVariant,
    annotations: VariantAnnotations,
    gene_mechanism: GeneMechanism,
) -> ComponentScore:
    """Spec 19.3 + 9.8. Branches by consequence class, then gated by mechanism.

    The branch choice is not cosmetic: M_coding and M_splice are computed from
    different predictors on different scales, so one formula across all classes
    would be a category error, and a pooled cross-class ranking of the results is
    not interpretable.
    """
    cls_ = variant.consequence_class
    weights = cfg.branches.for_class(cls_)
    score = _weighted(
        cfg, ScoreComponent.MECHANISM, weights, annotations, branch=f"M_{cls_.value}"
    )
    if not score.computable:
        return score

    gate_key = _gate_key(cls_, annotations)
    multiplier, key_used = cfg.gating.multiplier(gene_mechanism.mechanism, gate_key)

    score.pre_gate_value = score.value
    score.gate_multiplier = multiplier
    score.gate_key = key_used
    score.value = (score.value or 0.0) * multiplier
    if gene_mechanism.is_assumed:
        score.flags.append("mechanism_assumed_low_confidence")
    if gene_mechanism.mechanism is Mechanism.UNKNOWN:
        score.flags.append("mechanism_unknown_no_gating")
    return score


def _gate_key(cls_: ConsequenceClass, annotations: VariantAnnotations) -> str:
    """Map a variant onto the gating vocabulary in mechanism_gating.yaml."""
    if cls_ is ConsequenceClass.PTV:
        return "ptv"
    if cls_ is ConsequenceClass.SYNONYMOUS:
        return "synonymous"
    if cls_ is ConsequenceClass.INFRAME_INDEL:
        return "inframe_indel"
    if cls_ is ConsequenceClass.SPLICE:
        return "splice_lof"
    if cls_ is ConsequenceClass.NONCODING:
        return "noncoding_other"
    if cls_ is ConsequenceClass.MISSENSE:
        ctx = annotations.protein_context
        if ctx.active_site or ctx.binding_site:
            return "missense_in_functional_site"
        am = annotations.protein_effect.alphamissense_score
        if am.is_present and (am.value or 0.0) >= 0.7:
            return "missense_predicted_lof"
        if am.is_present and (am.value or 0.0) <= 0.3:
            return "missense_predicted_tolerated"
        return "missense_predicted_lof"
    return "default"


# ----------------------------------------------------------- editability (E)


def editability_score(cfg: ScoringConfig, annotations: VariantAnnotations) -> ComponentScore:
    """Spec 19.4, 17.1.1. Derived entirely from the named design tool's output.

    No LLM assigns this. The anchor maps one specific tool's output distribution,
    so a run whose design tool is unrecorded cannot be scored reproducibly.
    """
    sub = cfg.editability
    weights = {
        "pegrna_predicted_efficiency": sub.predicted_efficiency,
        "pegrna_design_count": sub.design_count,
        "bystander_risk": sub.bystander,
        "off_target_risk": sub.off_target,
        "indel_length": sub.indel_length,
    }
    score = _weighted(cfg, ScoreComponent.EDITABILITY, weights, annotations)
    if score.computable and not annotations.editability.design_tool:
        score.flags.append("design_tool_unrecorded")
    return score


# ------------------------------------------------------------- saturation (N)


def saturation(
    cfg: ScoringConfig,
    annotations: VariantAnnotations,
    evidence: list[Evidence],
) -> SaturationResult:
    """Spec 19.5. Existence and directness of direct evidence -- never 1 - F.

    A variant with strong direct evidence of NO effect has HIGH saturation and LOW
    F. Deriving saturation from F would invert that and re-promote it as a
    discovery candidate.
    """
    usable = [ev for ev in evidence if ev.usable_in_scoring]
    direct = [ev for ev in usable if ev.effective_level is EvidenceLevel.VARIANT_DIRECT]

    if not usable:
        # Never searched: saturation is UNKNOWN, not zero. A variant that was
        # merely skipped must not be able to claim the novelty bonus.
        return SaturationResult(
            value=None,
            level_key="none",
            known=False,
            may_claim_understudied=False,
            flags=["incomplete_evidence"],
        )

    if not direct:
        if annotations.mave.functional_class.is_present:
            return SaturationResult(
                value=cfg.saturation.level("mave_score_available"),
                level_key="mave_score_available",
            )
        assoc = [
            ev for ev in usable if ev.effective_level is EvidenceLevel.VARIANT_ASSOCIATION
        ]
        key = "variant_association_only" if assoc else "none"
        return SaturationResult(value=cfg.saturation.level(key), level_key=key)

    n_independent = len(distinct_observations([e for e in direct if e.observation]))
    relevant = [
        ev
        for ev in direct
        if ev.model_system_key in ("exact_cell_type", "related_human_neuronal")
    ]
    if relevant:
        key = (
            "variant_direct_relevant_cell_replicated"
            if n_independent > 1
            else "variant_direct_relevant_cell_single"
        )
    elif any(ev.model_system_key == "in_vitro_biochemical" for ev in direct):
        key = "variant_direct_biochemical_only"
    else:
        key = "variant_direct_other_model"

    return SaturationResult(value=cfg.saturation.level(key), level_key=key)

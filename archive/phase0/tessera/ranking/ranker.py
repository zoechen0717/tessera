"""Composite score and stratified ranking.

The composite (spec 19.5, adopted form):

    weighted_sum = w_D*D + w_F*F + w_M*M + w_E*E     (renormalized over present)
    S            = weighted_sum * (1 + alpha * (1 - saturation))

Novelty is a MULTIPLIER. It introduces no second copy of D or M, which is what
makes the configured weights equal the effective weights and keeps the sign of
F's contribution well defined. The rejected linear form and the algebra showing
why it fails are in spec 19.5.
"""

from __future__ import annotations

from statistics import variance

from tessera.config.scoring_config import ScoringConfig
from tessera.enums import ConsequenceClass, RankingMode, ScoreComponent
from tessera.schemas.annotations import VariantAnnotations
from tessera.schemas.evidence import Evidence
from tessera.schemas.gene import GeneMechanism
from tessera.schemas.scoring import (
    ComponentDiscrimination,
    RunDiagnostics,
    VariantScore,
)
from tessera.schemas.variant import CanonicalVariant
from tessera.ranking import components

DEGENERATE_VARIANCE_THRESHOLD = 1e-6


def score_variant(
    cfg: ScoringConfig,
    mode: RankingMode,
    variant: CanonicalVariant,
    annotations: VariantAnnotations,
    evidence: list[Evidence],
    gene_mechanism: GeneMechanism,
) -> VariantScore:
    """Score one variant.

    Pure with respect to the candidate set: this function never sees other
    variants, which is how Guarantee 9 (scores independent of candidate-set
    composition) is enforced structurally rather than by convention.
    """
    mode_cfg = cfg.mode(mode)
    my_evidence = [ev for ev in evidence if ev.vrs_id in (None, variant.vrs_id)]

    d = components.disease_score(cfg, variant, annotations, my_evidence)
    f = components.functional_score(cfg, annotations, my_evidence)
    m = components.mechanism_score(cfg, variant, annotations, gene_mechanism)
    e = components.editability_score(cfg, annotations)

    scores = {
        ScoreComponent.DISEASE: d,
        ScoreComponent.FUNCTIONAL: f,
        ScoreComponent.MECHANISM: m,
        ScoreComponent.EDITABILITY: e,
    }
    configured = mode_cfg.weights.as_dict()

    # Renormalize over computable components only. An uncomputable component is
    # excluded, never substituted with 0.0 (spec 42.1).
    usable = {c: s for c, s in scores.items() if s.computable}
    total_weight = sum(configured[c] for c in usable)
    effective: dict[ScoreComponent, float] = {}
    flags: list[str] = []

    if total_weight <= 0.0:
        weighted_sum = 0.0
        flags.append("no_computable_components")
    else:
        for c in usable:
            effective[c] = configured[c] / total_weight
        weighted_sum = sum((scores[c].value or 0.0) * effective[c] for c in usable)
        if abs(total_weight - 1.0) > 1e-9:
            flags.append("component_weights_renormalized")

    biological_priority = sum(
        (scores[c].value or 0.0) * effective.get(c, 0.0)
        for c in (ScoreComponent.DISEASE, ScoreComponent.FUNCTIONAL, ScoreComponent.MECHANISM)
        if c in usable
    )

    sat = components.saturation(cfg, annotations, my_evidence)
    alpha = mode_cfg.novelty_alpha
    if sat.value is None:
        # Saturation unknown -> no novelty bonus. An unresearched variant must not
        # be promoted as understudied (spec 42.1).
        multiplier = 1.0
        flags.extend(sat.flags)
    else:
        multiplier = 1.0 + alpha * (1.0 - sat.value)

    final = weighted_sum * multiplier

    blocking = list(annotations.blocking_flags())
    status = annotations.editability.resolve_status()
    if status.blocks_panel_inclusion():
        blocking.append(f"not_editable:{status.value}")

    for comp in scores.values():
        flags.extend(f"{comp.component.value}:{flag}" for flag in comp.flags)

    return VariantScore(
        vrs_id=variant.vrs_id,
        consequence_class=variant.consequence_class,
        ranking_mode=mode,
        disease=d,
        functional=f,
        mechanism=m,
        editability=e,
        biological_priority=biological_priority,
        weighted_sum=weighted_sum,
        saturation=sat,
        novelty_alpha=alpha,
        novelty_multiplier=multiplier,
        final_score=final,
        effective_weights=effective,
        gene_mechanism=gene_mechanism.mechanism,
        mechanism_is_assumed=gene_mechanism.is_assumed,
        components_that_did_work=sorted(usable, key=lambda c: c.value),
        uncomputable_components=sorted(
            (c for c in scores if c not in usable), key=lambda c: c.value
        ),
        flags=flags,
        blocking_flags=blocking,
        scoring_fingerprint=cfg.fingerprint(),
    )


def rank_stratified(
    scores: list[VariantScore],
) -> dict[ConsequenceClass, list[VariantScore]]:
    """Rank WITHIN consequence class (spec 19.3).

    Cross-class pooling is not offered here: M_missense and M_splice come from
    different predictors on different scales, so a pooled ordering has no physical
    meaning. A cross-class view may be displayed, but it is presentational and
    must be labeled as such -- see `presentational_global_order`.
    """
    strata: dict[ConsequenceClass, list[VariantScore]] = {}
    for s in scores:
        strata.setdefault(s.consequence_class, []).append(s)
    for group in strata.values():
        group.sort(key=lambda s: (-s.final_score, s.vrs_id))
    return strata


def presentational_global_order(scores: list[VariantScore]) -> list[VariantScore]:
    """A cross-class ordering for display only. NOT a ranking."""
    return sorted(scores, key=lambda s: (-s.final_score, s.vrs_id))


def diagnose(cfg: ScoringConfig, mode: RankingMode, scores: list[VariantScore]) -> RunDiagnostics:
    """Spec 48.9. Report which components actually drove the ordering.

    A component whose values are identical across candidates contributes nothing
    regardless of its configured weight. This is the check that catches the
    SETD1A + schizophrenia case, where disease evidence is real but gene-level and
    therefore constant across every candidate, leaving a configured w_D = 0.32
    inert while the report would otherwise imply it was doing work.
    """
    configured = cfg.mode(mode).weights.as_dict()
    finals = [s.final_score for s in scores]
    final_var = variance(finals) if len(finals) > 1 else 0.0

    discrimination: list[ComponentDiscrimination] = []
    degenerate: list[ScoreComponent] = []
    sparse: list[ScoreComponent] = []

    for comp, weight in configured.items():
        values = [s.component(comp).value for s in scores]
        present = [v for v in values if v is not None]
        n_present = len(present)
        coverage = n_present / len(values) if values else 0.0

        # Variance is only meaningful with at least two present values. A component
        # present on ONE candidate is sparse, not degenerate: it still moves that
        # candidate's score relative to the others, so claiming it "contributed
        # nothing" would overstate the finding.
        var = variance(present) if n_present > 1 else 0.0
        at_floor = sum(1 for v in present if v <= 1e-9)
        contribution = (weight * var / final_var) if final_var > 0 else 0.0

        is_sparse = weight > 0 and n_present < 2
        is_degenerate = weight > 0 and n_present >= 2 and var < DEGENERATE_VARIANCE_THRESHOLD
        if is_sparse:
            sparse.append(comp)
        if is_degenerate:
            degenerate.append(comp)

        discrimination.append(
            ComponentDiscrimination(
                component=comp,
                configured_weight=weight,
                variance=var,
                n_present=n_present,
                n_missing=len(values) - n_present,
                coverage=coverage,
                fraction_at_floor=(at_floor / n_present) if present else 1.0,
                effective_contribution=contribution,
                degenerate=is_degenerate,
                sparse=is_sparse,
            )
        )

    warnings: list[str] = []
    for comp in degenerate:
        warnings.append(
            f"component {comp.value!r} has a configured weight of "
            f"{configured[comp]:.2f} but is CONSTANT across all candidates that have "
            f"it: it contributed nothing to the ordering. Do not describe this run's "
            f"output as {comp.value}-weighted."
        )
    for comp in sparse:
        row = next(r for r in discrimination if r.component is comp)
        if row.n_present == 0:
            warnings.append(
                f"component {comp.value!r} has a configured weight of "
                f"{configured[comp]:.2f} but could not be computed for ANY candidate. "
                f"Its weight was renormalized away; the ranking is not "
                f"{comp.value}-informed."
            )
        else:
            warnings.append(
                f"component {comp.value!r} was computable for only {row.n_present} of "
                f"{len(scores)} candidates, so its discriminative power cannot be "
                f"assessed. It does affect those candidates' scores relative to the "
                f"rest -- that asymmetry is a coverage gap, not a finding."
            )

    return RunDiagnostics(
        n_candidates=len(scores),
        discrimination=discrimination,
        degenerate_components=degenerate,
        sparse_components=sparse,
        failed_retrievals=sum(
            1 for s in scores if any(f.startswith("retrieval_failure") for f in s.blocking_flags)
        ),
        panel_ineligible=sum(1 for s in scores if not s.panel_eligible),
        warnings=warnings,
    )

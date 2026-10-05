"""Score results. Every number carries how it was produced."""

from __future__ import annotations

from pydantic import Field

from tessera.enums import ConsequenceClass, Mechanism, RankingMode, ScoreComponent
from tessera.schemas.common import Strict


class TermContribution(Strict):
    """One raw input's journey into a component score."""

    raw_input: str
    raw_value: float | str | None = None
    anchored: float | None = None
    weight: float = 0.0
    weight_after_renormalization: float = 0.0
    contribution: float = 0.0
    state: str = "present"


class ComponentScore(Strict):
    """One of D, F, M, E.

    `value is None` means the component could not be computed at all (every input
    was NOT_SEARCHED or NOT_APPLICABLE). It is then excluded from the composite and
    the remaining weights renormalize -- it is NOT treated as zero (spec 42.1).
    """

    component: ScoreComponent
    value: float | None = None
    terms: list[TermContribution] = Field(default_factory=list)
    branch: str | None = Field(default=None, description="Which M branch was used")
    gate_multiplier: float = 1.0
    gate_key: str | None = None
    pre_gate_value: float | None = None
    weights_renormalized: bool = False
    inputs_present: int = 0
    inputs_missing: int = 0
    flags: list[str] = Field(default_factory=list)

    @property
    def computable(self) -> bool:
        return self.value is not None


class SaturationResult(Strict):
    """Spec 19.5. Defined on evidence existence and directness, NOT on 1 - F.

    That is what keeps "tested, null result" (high saturation, low F) separable
    from "never tested" (unknown saturation). Using 1 - F would re-promote a
    variant already shown to do nothing as a discovery candidate.
    """

    value: float | None = None
    level_key: str = "none"
    known: bool = True
    may_claim_understudied: bool = True
    flags: list[str] = Field(default_factory=list)


class VariantScore(Strict):
    """The complete scoring record for one variant.

    Deliberately verbose. Spec 21 forbids reporting a final score without its
    components, and spec 48.9 requires knowing which components actually did work.
    """

    vrs_id: str
    consequence_class: ConsequenceClass
    ranking_mode: RankingMode

    disease: ComponentScore
    functional: ComponentScore
    mechanism: ComponentScore
    editability: ComponentScore

    biological_priority: float
    weighted_sum: float
    saturation: SaturationResult
    novelty_alpha: float
    novelty_multiplier: float
    final_score: float

    effective_weights: dict[ScoreComponent, float] = Field(default_factory=dict)
    gene_mechanism: Mechanism = Mechanism.UNKNOWN
    mechanism_is_assumed: bool = False

    components_that_did_work: list[ScoreComponent] = Field(default_factory=list)
    uncomputable_components: list[ScoreComponent] = Field(default_factory=list)
    flags: list[str] = Field(default_factory=list)
    blocking_flags: list[str] = Field(default_factory=list)

    scoring_fingerprint: dict[str, str] = Field(default_factory=dict)

    def component(self, which: ScoreComponent) -> ComponentScore:
        return {
            ScoreComponent.DISEASE: self.disease,
            ScoreComponent.FUNCTIONAL: self.functional,
            ScoreComponent.MECHANISM: self.mechanism,
            ScoreComponent.EDITABILITY: self.editability,
        }[which]

    @property
    def panel_eligible(self) -> bool:
        return not self.blocking_flags


class ComponentDiscrimination(Strict):
    """Spec 48.9. The metric that catches a zero-variance component.

    A component whose values are identical across all candidates contributes
    nothing to the ordering regardless of its configured weight. Reporting a
    "disease-weighted score" under those conditions is misleading.
    """

    component: ScoreComponent
    configured_weight: float
    variance: float
    n_present: int
    n_missing: int
    coverage: float
    fraction_at_floor: float
    effective_contribution: float
    degenerate: bool = Field(
        default=False,
        description="Present for >=2 candidates yet constant: contributes nothing to "
        "the ordering, whatever its configured weight says.",
    )
    sparse: bool = Field(
        default=False,
        description="Present for fewer than 2 candidates, so variance cannot be "
        "assessed. NOT the same as degenerate: the component still moves the score "
        "of the candidates that do have it, relative to those that do not.",
    )


class RunDiagnostics(Strict):
    n_candidates: int
    discrimination: list[ComponentDiscrimination] = Field(default_factory=list)
    degenerate_components: list[ScoreComponent] = Field(default_factory=list)
    sparse_components: list[ScoreComponent] = Field(default_factory=list)
    quarantined_variants: int = 0
    failed_retrievals: int = 0
    panel_ineligible: int = 0
    warnings: list[str] = Field(default_factory=list)

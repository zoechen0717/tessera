"""Typed evidence objects (spec 12) and critic results (spec 14)."""

from __future__ import annotations

import hashlib
from datetime import datetime

from pydantic import Field, model_validator

from tessera.enums import (
    DiseaseMatchType,
    EvidenceDirection,
    EvidenceLevel,
    VerificationStatus,
)
from tessera.schemas.common import Provenance, Strict


class ObservationKey(Strict):
    """Spec 12.4. Deduplication at the PMID level is insufficient.

    The same proband appears in a case report, a later meta-analysis, a review,
    and a ClinVar submission: four sources, one observation. Counting those as
    four inflates D for exactly the variants most likely to be already known,
    which then also distorts the novelty multiplier in the wrong direction.
    """

    first_report_source_id: str
    cohort_id: str | None = None
    proband_id: str | None = None
    proband_resolution: str = Field(
        default="cohort_only",
        description="proband_level | cohort_only -- cohort_only caps how much this "
        "observation may contribute to recurrence",
    )

    @property
    def key(self) -> str:
        return f"{self.first_report_source_id}|cohort:{self.cohort_id or '-'}|proband:{self.proband_id or '-'}"

    def __hash__(self) -> int:  # type: ignore[override]
        return hash(self.key)


class DeterministicCriticResult(Strict):
    """Spec 14.1. A hard gate: failure REJECTS, it does not warn."""

    status: str = Field(default="pending", description="pending | passed | rejected")
    source_reachable: bool | None = None
    variant_alias_matched: bool | None = None
    matched_alias: str | None = None
    span_present_in_source: bool | None = None
    schema_valid: bool | None = None
    duplicate_group: str | None = None
    rejection_reasons: list[str] = Field(default_factory=list)


class LlmCriticResult(Strict):
    """Spec 14.2. Advisory. May downgrade evidence_level; may never promote it.

    The asymmetry is deliberate: the model may express doubt about evidence, but
    the category an item belongs to is established deterministically.
    """

    status: str = Field(default="pending", description="pending | reviewed | skipped")
    assay_directness: EvidenceLevel | None = None
    model_system_relevance: str | None = None
    claim_supported_by_span: bool | None = None
    disease_match_type: DiseaseMatchType | None = None
    level_downgraded_from: EvidenceLevel | None = None
    warnings: list[str] = Field(default_factory=list)


class Evidence(Strict):
    """One typed evidence object.

    Never prose. An extractor that cannot fill `supporting_span` and `source_id`
    has not produced evidence (spec 28).
    """

    evidence_id: str
    vrs_id: str | None = None
    gene_symbol: str
    disease_id: str | None = None

    evidence_type: str
    evidence_level: EvidenceLevel
    direction: EvidenceDirection = EvidenceDirection.NEUTRAL

    assay: str | None = None
    model_system: str | None = None
    cell_type: str | None = None
    comparison: str | None = None
    phenotype: str | None = None
    effect_direction: str | None = None
    result_summary: str | None = None

    directness_key: str | None = Field(
        default=None,
        description="Key into the evidence_directness anchor, e.g. "
        "'variant_direct_relevant_human_cell'",
    )
    model_system_key: str | None = Field(
        default=None, description="Key into the model_system_relevance anchor"
    )

    disease_match_type: DiseaseMatchType = DiseaseMatchType.RELATED
    disease_match_detail: str | None = None

    observation: ObservationKey | None = None
    cohort_size: int | None = None

    source_type: str
    source_id: str
    source_url: str | None = None
    publication_year: int | None = None

    claim: str | None = None
    supporting_span: str | None = None

    agent_extraction_confidence: float | None = Field(
        default=None,
        ge=0.0,
        le=1.0,
        description="Extraction certainty, NOT biological truth and NOT calibrated. "
        "May order the review queue; must never enter any score (spec 43.1). "
        "Enforced by source inspection in tests/test_invariants.py.",
    )

    verification_status: VerificationStatus = VerificationStatus.PENDING
    critic_deterministic: DeterministicCriticResult = Field(
        default_factory=DeterministicCriticResult
    )
    critic_llm: LlmCriticResult = Field(default_factory=LlmCriticResult)

    retrieval_query: str | None = None
    retrieved_at: datetime | None = None
    extractor_model: str | None = None
    extractor_prompt_version: str | None = None
    provenance: Provenance | None = None

    @model_validator(mode="after")
    def _variant_specific_requires_variant(self) -> "Evidence":
        if self.evidence_level.is_variant_specific() and not self.vrs_id:
            raise ValueError(
                f"evidence_level={self.evidence_level.value} requires a vrs_id. "
                f"Variant-specific evidence without a resolved variant is the "
                f"gene-level-mislabeled-as-variant-level failure (spec 53.3)."
            )
        return self

    @model_validator(mode="after")
    def _citation_required(self) -> "Evidence":
        """Spec 28: no claim enters the store without source id and location."""
        if self.claim and not self.supporting_span:
            raise ValueError(
                "a claim requires supporting_span (spec 28). Evidence without a "
                "quoted location is rejected, not stored with a warning."
            )
        return self

    @property
    def observation_key(self) -> str:
        """Falls back to the source id so that dedup always has a key to group on."""
        return self.observation.key if self.observation else f"{self.source_id}|cohort:-|proband:-"

    @property
    def usable_in_scoring(self) -> bool:
        """Only verified evidence scores. Pending and rejected do not."""
        if self.critic_deterministic.status == "rejected":
            return False
        return self.verification_status.usable_in_scoring()

    @property
    def effective_level(self) -> EvidenceLevel:
        """The LLM critic may downgrade; it may never promote."""
        proposed = self.critic_llm.assay_directness
        if proposed is not None and proposed.rank > self.evidence_level.rank:
            return proposed
        return self.evidence_level


def make_evidence_id(*parts: str) -> str:
    return "ev_" + hashlib.sha256("|".join(parts).encode()).hexdigest()[:16]


def distinct_observations(items: list[Evidence]) -> dict[str, list[Evidence]]:
    """Group by observation key (spec 19.1.2).

    D is computed over deduplicated observation keys, not over evidence objects.
    Where a ClinVar submission and a publication resolve to the same key, the
    publication is retained as primary and the submission linked as secondary.
    """
    groups: dict[str, list[Evidence]] = {}
    for ev in items:
        groups.setdefault(ev.observation_key, []).append(ev)
    for key, group in groups.items():
        group.sort(key=lambda e: (e.source_type != "publication", e.source_id))
    return groups

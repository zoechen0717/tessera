"""Claims, evidence links and validation decisions (SCHEMA §6–8).

A Claim is what a source says, stored once. An EvidenceLink says which allele
the claim applies to and how directly. A ValidationDecision says whether the
link may be used. Ranking reads only links whose effective state is accepted.
"""

from __future__ import annotations

from enum import StrEnum
from typing import Literal, Self

from pydantic import Field, model_validator

from tessera.schemas.base import SCHEMA_VERSION, Strict, ValidationState
from tessera.schemas.identity import IdentityMatch
from tessera.schemas.observed import Observed
from tessera.schemas.source import SupportingLocator


class EvidenceKind(StrEnum):
    PATIENT_OBSERVATION = "patient_observation"
    ASSOCIATION = "association"
    SEGREGATION = "segregation"
    CLINICAL_CLASSIFICATION = "clinical_classification"
    FUNCTIONAL_EXPERIMENT = "functional_experiment"
    COMPUTATIONAL_PREDICTION = "computational_prediction"
    GENE_MECHANISM = "gene_mechanism"
    DATASET_DESCRIPTION = "dataset_description"
    REVIEW_SUMMARY = "review_summary"
    HYPOTHESIS = "hypothesis"


class ClaimScope(StrEnum):
    VARIANT = "variant"
    RESIDUE = "residue"
    DOMAIN = "domain"
    MUTATION_CLASS = "mutation_class"
    GENE = "gene"
    LOCUS = "locus"
    PATHWAY = "pathway"


class PerturbationType(StrEnum):
    EXACT_ALLELE_EDIT = "exact_allele_edit"
    ALLELE_CONSTRUCT = "allele_construct"
    CARRIER_COMPARISON = "carrier_comparison"
    GENE_KNOCKOUT = "gene_knockout"
    KNOCKDOWN = "knockdown"
    OVEREXPRESSION = "overexpression"
    OTHER = "other"

    @property
    def installs_allele(self) -> bool:
        """Whether the experiment tested the allele itself, not a carrier or the gene."""
        return self in (PerturbationType.EXACT_ALLELE_EDIT, PerturbationType.ALLELE_CONSTRUCT)


class OutcomeType(StrEnum):
    ALTERED = "altered"
    NO_DETECTED_EFFECT = "no_detected_effect"
    RESCUE = "rescue"
    MIXED = "mixed"
    INCONCLUSIVE = "inconclusive"


class FunctionalAssay(Strict):
    assay_name: str
    assay_type: str
    tested_allele_mention_ids: list[str] = Field(default_factory=list)
    perturbation_type: PerturbationType
    zygosity: Literal["het", "hom", "hemi", "unknown"] = "unknown"
    comparator: str
    endpoint: str
    outcome_type: OutcomeType
    cell_context: str | None = None
    confound_flags: list[str] = Field(default_factory=list)


class AssociationRecord(Strict):
    """Case/control counts for one cohort. Counts are never inferred."""

    cohort_id: str
    study_design: str
    cases_with_allele: Observed[int]
    case_total: Observed[int]
    controls_with_allele: Observed[int]
    control_total: Observed[int]
    analysis_unit: Literal["allele", "mutation_class", "gene", "locus"]

    @property
    def has_control_denominator(self) -> bool:
        return self.controls_with_allele.is_present and self.control_total.is_present


class Claim(Strict):
    schema_version: str = SCHEMA_VERSION
    claim_id: str
    revision: int = Field(ge=1)
    evidence_kind: EvidenceKind
    scope: ClaimScope
    subject_gene_ids: list[str]
    source_reported_mention_ids: list[str] = Field(default_factory=list)
    claim_text: str
    functional_assay: FunctionalAssay | None = None
    association: AssociationRecord | None = None
    locators: list[SupportingLocator] = Field(min_length=1)
    source_snapshot_ids: list[str] = Field(min_length=1)
    publication_group_id: str | None = None
    extraction_method: Literal["human", "deterministic", "llm"]

    @model_validator(mode="after")
    def _payload_matches_kind(self) -> Self:
        if self.evidence_kind is EvidenceKind.FUNCTIONAL_EXPERIMENT and self.functional_assay is None:
            raise ValueError("functional_experiment claims require functional_assay")
        if self.evidence_kind in (EvidenceKind.PATIENT_OBSERVATION, EvidenceKind.ASSOCIATION) and (
            self.association is None
        ):
            raise ValueError(f"{self.evidence_kind.value} claims require an association record")
        if self.evidence_kind is EvidenceKind.HYPOTHESIS:
            raise ValueError("hypotheses belong in interpretation outputs, never in claims")
        for loc in self.locators:
            if loc.snapshot_id not in self.source_snapshot_ids:
                raise ValueError(f"locator snapshot {loc.snapshot_id} not in source_snapshot_ids")
        return self

    @property
    def key(self) -> tuple[str, int]:
        return (self.claim_id, self.revision)


class ContextRelationship(StrEnum):
    EXACT = "exact"
    APPROVED = "approved"
    RELATED = "related"
    OTHER = "other"
    UNKNOWN = "unknown"


class DiseaseRelationship(StrEnum):
    EXACT = "exact"
    APPROVED_SYNONYM = "approved_synonym"
    RELATED = "related"
    UNRELATED = "unrelated"
    UNKNOWN = "unknown"

    @property
    def is_query_disease(self) -> bool:
        return self in (DiseaseRelationship.EXACT, DiseaseRelationship.APPROVED_SYNONYM)


class EvidenceLink(Strict):
    schema_version: str = SCHEMA_VERSION
    link_id: str
    claim_id: str
    claim_revision: int
    variant_id: str | None
    mention_id: str | None = Field(
        default=None, description="Source mention this link was resolved through"
    )
    identity_match: IdentityMatch
    identity_decision_ids: list[str] = Field(default_factory=list)
    disease_relationship: DiseaseRelationship = DiseaseRelationship.UNKNOWN
    cell_context_relationship: ContextRelationship = ContextRelationship.UNKNOWN
    inference: bool = False


class CheckResult(StrEnum):
    PASS = "pass"
    FAIL = "fail"
    UNKNOWN = "unknown"
    NOT_APPLICABLE = "not_applicable"


class ValidationDecision(Strict):
    schema_version: str = SCHEMA_VERSION
    decision_id: str
    link_id: str
    claim_id: str
    claim_revision: int
    policy_version: str
    checks: dict[str, CheckResult] = Field(default_factory=dict)
    reasons: list[str] = Field(default_factory=list)
    reviewer_type: Literal["host", "model", "human"]
    reviewer_id: str
    human_review: Literal["unreviewed", "accepted", "rejected", "needs_revision"] = "unreviewed"
    effective_state: ValidationState
    sequence: int = Field(ge=0, description="Append order; later decisions do not erase earlier ones")

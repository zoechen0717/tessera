"""Deterministic annotations and gene mechanism (SCHEMA §2.5, §4.1, §9)."""

from __future__ import annotations

from enum import StrEnum
from typing import Literal

from pydantic import Field

from tessera.schemas.base import SCHEMA_VERSION, Strict
from tessera.schemas.observed import Observed


class TranscriptEffect(Strict):
    transcript_id: str
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    consequences: list[str] = Field(min_length=1, description="Sequence Ontology terms")
    selected_display: bool = False


class VariantAnnotation(Strict):
    """All deterministic annotations for one allele.

    Each field is Observed so that "not looked up", "looked up and absent" and
    "lookup failed" stay distinct all the way into the features.
    """

    schema_version: str = SCHEMA_VERSION
    variant_id: str
    transcript_effects: list[TranscriptEffect] = Field(default_factory=list)
    nmd: Observed[Literal["triggers", "escapes"]] = Field(default_factory=Observed.not_searched)
    in_functional_domain: Observed[bool] = Field(default_factory=Observed.not_searched)
    domain_name: str | None = None
    predictions: dict[str, Observed[float]] = Field(
        default_factory=dict, description="predictor id -> raw score, source-defined direction"
    )
    population_af: Observed[float] = Field(default_factory=Observed.not_searched)
    clinvar_classification: Observed[str] = Field(default_factory=Observed.not_searched)

    def selected_effect(self) -> TranscriptEffect | None:
        chosen = [e for e in self.transcript_effects if e.selected_display]
        return chosen[0] if len(chosen) == 1 else None


class Mechanism(StrEnum):
    LOSS_OF_FUNCTION = "loss_of_function"
    GAIN_OF_FUNCTION = "gain_of_function"
    DOMINANT_NEGATIVE = "dominant_negative"
    MIXED = "mixed"
    UNKNOWN = "unknown"


class MechanismBasis(Strict):
    source: str
    record_id: str
    field: str
    value: str
    snapshot_id: str | None = None


class GeneMechanismAssessment(Strict):
    """Gene–disease mechanism. Sets K only through the concordance table; never
    copied into any allele's evidence (INV-17)."""

    schema_version: str = SCHEMA_VERSION
    gene_id: str
    disease_concept: str
    mechanism: Mechanism
    status: Literal["assessed", "searched_not_found", "not_searched"]
    basis: list[MechanismBasis] = Field(default_factory=list)
    decision: Literal["rule", "human", "imported"]
    policy_id: str


class GeneEntity(Strict):
    gene_id: str
    symbol: str
    chrom: str
    strand: Literal["+", "-"]
    is_synthetic: bool = False

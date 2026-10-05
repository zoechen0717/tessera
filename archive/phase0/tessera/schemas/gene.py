"""Resolved gene and its mechanism prior."""

from __future__ import annotations

from datetime import datetime

from pydantic import Field

from tessera.enums import Mechanism, MechanismConfidence
from tessera.schemas.common import Provenance, Quantity, Strict


class ResolvedGene(Strict):
    """Spec 7.1. Resolution happens before any search."""

    symbol: str
    hgnc_id: str | None = None
    ensembl_gene_id: str | None = None
    entrez_id: str | None = None
    uniprot_id: str | None = None
    mane_select_transcript: str | None = None
    mane_plus_clinical: list[str] = Field(default_factory=list)
    transcripts: list[str] = Field(default_factory=list)
    chrom: str | None = None
    start: int | None = None
    end: int | None = None
    strand: int | None = None
    genome_build: str = "GRCh38"
    aliases: list[str] = Field(default_factory=list)
    prev_symbols: list[str] = Field(default_factory=list)
    resolved_at: datetime | None = None
    provenance: list[Provenance] = Field(default_factory=list)


class MechanismEvidenceItem(Strict):
    source: str
    value: str


class GeneMechanism(Strict):
    """Spec 9.8. A first-class input, not a future extension.

    `UNKNOWN` with `LOW` confidence is a correct answer and must be preferred to
    a guess: the ranker then falls back to a mechanism-agnostic profile and the
    report says so, rather than silently gating on an assumption.
    """

    gene_symbol: str
    mechanism: Mechanism = Mechanism.UNKNOWN
    confidence: MechanismConfidence = MechanismConfidence.LOW
    evidence: list[MechanismEvidenceItem] = Field(default_factory=list)

    gnomad_loeuf: Quantity = Field(default_factory=Quantity.not_searched)
    gnomad_pli: Quantity = Field(default_factory=Quantity.not_searched)
    clingen_haploinsufficiency_score: Quantity = Field(default_factory=Quantity.not_searched)

    ptv_burden_significant: bool | None = None
    missense_burden_significant: bool | None = None

    resolved_at: datetime | None = None
    provenance: list[Provenance] = Field(default_factory=list)

    @property
    def is_assumed(self) -> bool:
        """True when gating is being applied on weak grounds. The report must
        disclose this rather than presenting a gated score as settled."""
        return self.mechanism is not Mechanism.UNKNOWN and (
            self.confidence is MechanismConfidence.LOW
        )

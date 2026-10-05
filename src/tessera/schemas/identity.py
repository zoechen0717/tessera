"""Allele identity: mentions, identity decisions, canonical variants (SCHEMA §4).

The canonical genomic allele is the primary key. HGVS, rsIDs and protein
shorthand are aliases; none of them is unique, so none is an identifier.
"""

from __future__ import annotations

from enum import StrEnum
from typing import Literal, Self

from pydantic import Field, model_validator

from tessera.schemas.base import SCHEMA_VERSION, Strict, content_digest

IDENTITY_VERSION = "tessera-identity-1"

_DNA = set("ACGT")


class AlleleType(StrEnum):
    SNV = "snv"
    DELETION = "deletion"
    INSERTION = "insertion"
    DELINS = "delins"
    MNV = "mnv"


class IdentityStatus(StrEnum):
    RESOLVED = "resolved"
    AMBIGUOUS = "ambiguous"
    INSUFFICIENT = "insufficient"
    REF_MISMATCH = "ref_mismatch"
    UNSUPPORTED = "unsupported"
    UNMAPPED = "unmapped"


class IdentityMatch(StrEnum):
    """SCHEMA §2.2. Only the first two establish an exact genomic allele link."""

    EXACT_GENOMIC = "exact_genomic"
    TRANSCRIPT_SUPPORTED_EXACT = "transcript_supported_exact"
    PROTEIN_EQUIVALENT_ONLY = "protein_equivalent_only"
    SAME_RESIDUE = "same_residue"
    SAME_DOMAIN = "same_domain"
    MUTATION_CLASS = "mutation_class"
    GENE_ONLY = "gene_only"
    LOCUS_ONLY = "locus_only"
    PATHWAY_ONLY = "pathway_only"
    AMBIGUOUS = "ambiguous"
    NO_MATCH = "no_match"

    @property
    def is_exact(self) -> bool:
        return self in (IdentityMatch.EXACT_GENOMIC, IdentityMatch.TRANSCRIPT_SUPPORTED_EXACT)


def compute_variant_id(
    *, identity_version: str, assembly: str, reference_accession: str, pos_1based: int, ref: str, alt: str
) -> str:
    """SCHEMA §4.1: SHA-256 of canonical identity material, prefix `var_`."""
    material = {
        "identity_version": identity_version,
        "assembly": assembly,
        "reference_accession": reference_accession,
        "pos_1based": pos_1based,
        "ref": ref,
        "alt": alt,
    }
    return "var_" + content_digest(material)


class VariantMention(Strict):
    """A variant as a source wrote it. Kept even when it cannot be resolved."""

    schema_version: str = SCHEMA_VERSION
    mention_id: str
    raw_text: str
    discovery_route: str
    source_snapshot_id: str | None = None
    reported_assembly: str | None = None
    reported_chrom: str | None = None
    reported_pos: int | None = None
    reported_ref: str | None = None
    reported_alt: str | None = None
    parent_mention_id: str | None = Field(
        default=None, description="Set when split from a multiallelic record"
    )
    roles: list[str] = Field(default_factory=list)


class IdentityDecision(Strict):
    schema_version: str = SCHEMA_VERSION
    decision_id: str
    mention_id: str
    status: IdentityStatus
    variant_ids: list[str] = Field(default_factory=list)
    mapping_method: str
    resource_versions: dict[str, str] = Field(default_factory=dict)
    justification: str
    validation_errors: list[str] = Field(default_factory=list)
    decided_by: str = "host"

    @model_validator(mode="after")
    def _resolved_means_exactly_one(self) -> Self:
        if self.status is IdentityStatus.RESOLVED and len(self.variant_ids) != 1:
            raise ValueError(
                "resolved requires exactly one concrete allele, not an arbitrary choice "
                f"among {len(self.variant_ids)}"
            )
        if self.status is not IdentityStatus.RESOLVED and len(self.variant_ids) == 1:
            raise ValueError(f"status={self.status.value} cannot name a single resolved allele")
        return self


class Alias(Strict):
    raw_value: str
    alias_type: Literal["vcf", "hgvs_g", "hgvs_c", "hgvs_p", "rsid", "shorthand", "custom"]
    mention_id: str
    decision_id: str


class CanonicalVariant(Strict):
    schema_version: str = SCHEMA_VERSION
    variant_id: str
    identity_version: str = IDENTITY_VERSION
    assembly: str
    reference_accession: str
    reference_bundle_id: str
    chrom: str
    pos_1based: int = Field(ge=1)
    ref: str
    alt: str
    allele_type: AlleleType
    normalization_status: Literal["reference_validated"] = "reference_validated"
    normalization_tool: str
    normalization_version: str
    gene_ids: list[str] = Field(default_factory=list)
    aliases: list[Alias] = Field(default_factory=list)

    @model_validator(mode="after")
    def _identity(self) -> Self:
        for name in ("ref", "alt"):
            seq = getattr(self, name)
            if not seq or not set(seq) <= _DNA:
                raise ValueError(f"{name} must be a nonempty uppercase literal ACGT allele")
        if self.ref == self.alt:
            raise ValueError("ref and alt must differ")
        expected = compute_variant_id(
            identity_version=self.identity_version,
            assembly=self.assembly,
            reference_accession=self.reference_accession,
            pos_1based=self.pos_1based,
            ref=self.ref,
            alt=self.alt,
        )
        if self.variant_id != expected:
            raise ValueError("variant_id does not match canonical identity material")
        return self

    @property
    def vcf_key(self) -> tuple[str, int, str, str]:
        return (self.chrom, self.pos_1based, self.ref, self.alt)

    @property
    def label(self) -> str:
        return f"{self.chrom}:{self.pos_1based}:{self.ref}>{self.alt}"


class Exclusion(Strict):
    """A mention that does not enter the ranked registry, with the reason."""

    mention_id: str
    decision_id: str
    status: IdentityStatus
    reason: str

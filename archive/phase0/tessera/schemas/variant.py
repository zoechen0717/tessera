"""Canonical variant identity, aliases, and the quarantine record."""

from __future__ import annotations

from datetime import datetime

from pydantic import Field, model_validator

from tessera.enums import ConsequenceClass, ControlType
from tessera.schemas.common import Provenance, Strict


class VariantAlias(Strict):
    """Every observed spelling, with where it came from.

    This set is not bookkeeping: the deterministic critic's variant-mention check
    (spec 14.1) is implemented as a string match of this set against retrieved
    full text, and it is the highest-value check in the verification layer.
    """

    alias: str
    kind: str = Field(description="rsid | hgvs_g | hgvs_c | hgvs_p | shorthand | caid | custom")
    source: str
    transcript: str | None = None


class CanonicalVariant(Strict):
    """Spec 7.3.1. VRS ID is the internal primary key; CAID resolves aliases."""

    vrs_id: str
    caid: str | None = None

    genome_build: str = "GRCh38"
    chrom: str
    pos: int
    ref: str
    alt: str

    rsid: str | None = None
    hgvs_g: str | None = None
    hgvs_c: str | None = None
    hgvs_p: str | None = None

    gene_symbol: str
    transcript: str | None = None
    variant_type: str | None = None

    consequence: str | None = None
    consequence_class: ConsequenceClass = ConsequenceClass.OTHER

    source_build: str | None = Field(
        default=None,
        description="Build the variant was reported in, before liftover. Much of the "
        "psychiatric genetics literature reports GRCh37 (spec 7.3.2).",
    )
    lifted_over: bool = False

    aliases: list[VariantAlias] = Field(default_factory=list)

    is_control: bool = False
    control_type: ControlType | None = None

    discovered_by_routes: list[str] = Field(default_factory=list)
    created_at: datetime | None = None
    provenance: list[Provenance] = Field(default_factory=list)

    @model_validator(mode="after")
    def _build_is_grch38(self) -> "CanonicalVariant":
        if self.genome_build != "GRCh38":
            raise ValueError(
                f"canonical variants are GRCh38 only, got {self.genome_build!r}. "
                f"Lift over before constructing, and record source_build."
            )
        return self

    @model_validator(mode="after")
    def _control_type_matches_flag(self) -> "CanonicalVariant":
        if self.is_control and self.control_type is None:
            raise ValueError(
                "controls must declare a control_type -- spec 8.5 requires them to be "
                "explicitly labeled so they are not mixed into discovery ranking"
            )
        if self.control_type is not None and not self.is_control:
            raise ValueError("control_type set on a non-control variant")
        return self

    def alias_strings(self) -> set[str]:
        out = {a.alias for a in self.aliases}
        for field in (self.rsid, self.hgvs_g, self.hgvs_c, self.hgvs_p, self.caid):
            if field:
                out.add(field)
        return out

    @property
    def short_label(self) -> str:
        return self.hgvs_p or self.hgvs_c or self.rsid or self.vrs_id


class QuarantinedVariant(Strict):
    """Spec 7.3.3. A variant that could not be normalized is recorded, never dropped.

    Silent dropping makes candidate_recall report a falsely high number: a variant
    that was never normalized cannot be recovered, but must not be counted as
    absent from the source either.
    """

    raw_string: str
    source: str
    source_context: str | None = None
    failure_reason: str = Field(
        description="ambiguous_transcript | liftover_failed | unparseable_hgvs "
        "| ref_mismatch | not_found | multi_allelic_ambiguous"
    )
    attempted_resolutions: list[str] = Field(default_factory=list)
    quarantined_at: datetime | None = None

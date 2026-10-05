"""Shared value types: Quantity and Provenance.

Quantity exists so that "no value" always carries a reason. Every numeric input
to the scorer is a Quantity, never a bare float, which makes it impossible to
write a scorer that treats an unsearched field as a zero.
"""

from __future__ import annotations

from datetime import date, datetime
from typing import Any, Self

from pydantic import BaseModel, ConfigDict, Field, model_validator

from tessera.enums import Missingness


class Strict(BaseModel):
    """Base model: extra fields are an error, not a silent pass-through.

    An extractor that invents a field should fail validation rather than have the
    field quietly dropped -- otherwise a schema change looks like a no-op.
    """

    model_config = ConfigDict(extra="forbid", frozen=False, use_enum_values=False)


class Quantity(Strict):
    """A number that may be absent, with the reason it is absent.

    The scorer branches on `state` (spec 42.1):

      PRESENT             -> use `value`
      SEARCHED_NOT_FOUND  -> contributes 0.0, and legitimately lowers saturation
      NOT_SEARCHED        -> excluded from the weighted sum; weights renormalize
      NOT_APPLICABLE      -> excluded; does not count as incomplete
      NOT_AVAILABLE       -> excluded
      FAILED_RETRIEVAL    -> excluded, and blocks panel inclusion
    """

    value: float | None = None
    state: Missingness = Missingness.NOT_SEARCHED
    raw: Any = Field(default=None, description="Pre-anchor raw value, for audit")
    source: str | None = None
    source_release: str | None = None
    retrieved_at: datetime | None = None

    @model_validator(mode="after")
    def _state_matches_value(self) -> Self:
        if self.state is Missingness.PRESENT and self.value is None and self.raw is None:
            raise ValueError(
                "state=present requires a value or a raw value "
                "(categorical inputs carry only `raw`)"
            )
        if self.state is not Missingness.PRESENT and self.value is not None:
            raise ValueError(
                f"state={self.state.value} must not carry a value "
                f"(got {self.value!r}); absence of a value is the point"
            )
        return self

    # --- constructors -------------------------------------------------------

    @classmethod
    def present(cls, value: float, **kw: Any) -> Quantity:
        return cls(value=value, state=Missingness.PRESENT, **kw)

    @classmethod
    def category(cls, raw: str | bool, **kw: Any) -> Quantity:
        """A present value whose anchor is categorical or boolean.

        Such inputs have no meaningful float: `protein_domain_context` is
        "active_site", not 0.9. The anchor turns the category into the number, so
        storing a float here would mean anchoring twice.
        """
        return cls(state=Missingness.PRESENT, raw=raw, **kw)

    @classmethod
    def absent(cls, state: Missingness, **kw: Any) -> Quantity:
        if state is Missingness.PRESENT:
            raise ValueError("use Quantity.present() for present values")
        return cls(state=state, **kw)

    @classmethod
    def searched_not_found(cls, **kw: Any) -> Quantity:
        return cls.absent(Missingness.SEARCHED_NOT_FOUND, **kw)

    @classmethod
    def not_searched(cls, **kw: Any) -> Quantity:
        return cls.absent(Missingness.NOT_SEARCHED, **kw)

    @classmethod
    def not_applicable(cls, **kw: Any) -> Quantity:
        return cls.absent(Missingness.NOT_APPLICABLE, **kw)

    @classmethod
    def failed(cls, **kw: Any) -> Quantity:
        return cls.absent(Missingness.FAILED_RETRIEVAL, **kw)

    # --- scoring semantics --------------------------------------------------

    @property
    def is_present(self) -> bool:
        return self.state is Missingness.PRESENT

    @property
    def contributes_zero(self) -> bool:
        """Searched and genuinely nothing there: weak negative information.

        This is the ONLY absent state that contributes a number. It is also the
        only absent state that may lower saturation, because it is the only one
        that constitutes evidence about the variant rather than about the run.
        """
        return self.state is Missingness.SEARCHED_NOT_FOUND

    @property
    def enters_weighted_sum(self) -> bool:
        return self.is_present or self.contributes_zero

    @property
    def scoring_value(self) -> float | None:
        if self.is_present:
            return self.value
        if self.contributes_zero:
            return 0.0
        return None

    @property
    def anchor_input(self) -> Any:
        """What gets handed to the anchor: the raw value when there is one.

        Numeric inputs carry `value`; categorical and boolean inputs carry only
        `raw`. Routing both through one property keeps the scorer from having to
        know which kind of anchor it is about to hit.
        """
        if not self.is_present:
            raise ValueError(
                f"anchor_input requested for a non-present Quantity "
                f"(state={self.state.value}); absence is handled by the caller"
            )
        return self.raw if self.raw is not None else self.value

    @property
    def blocks_panel(self) -> bool:
        return self.state is Missingness.FAILED_RETRIEVAL


class Provenance(Strict):
    """Spec section 5.4. Attached to every claim that can move a rank."""

    source: str
    source_id: str | None = None
    source_url: str | None = None
    source_release: str | None = None
    retrieved_at: datetime
    query: str | None = None
    lookup_key: str | None = None
    response_hash: str | None = None
    document_hash: str | None = None
    extractor_model: str | None = None
    extractor_prompt_version: str | None = None
    license: str | None = None


class SourceSnapshot(Strict):
    """Spec sections 30, 67. One entry per source that contributed to a run."""

    source: str
    release: str | None = None
    access: str = Field(
        default="api",
        description="api | local_table | local_tsv | local_vcf_slice | release_file",
    )
    accessed_at: datetime | None = None
    response_hash: str | None = None
    license: str | None = None
    redistribution_of_scores: bool | None = None


class DateStamp(Strict):
    as_of: date

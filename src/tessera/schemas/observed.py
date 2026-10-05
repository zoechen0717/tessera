"""Observed[T]: a value together with why it is or is not there (SCHEMA §2.1).

`null` never means zero, false, wild type or benign. The status says whether the
value was measured, searched for and absent, never searched, or could not be
retrieved -- and those cases rank differently.
"""

from __future__ import annotations

from enum import StrEnum
from typing import Generic, Self, TypeVar

from pydantic import Field, model_validator

from tessera.schemas.base import Strict

T = TypeVar("T")


class ObservationStatus(StrEnum):
    VALUE_PRESENT = "value_present"
    NOT_SEARCHED = "not_searched"
    SEARCHED_NOT_FOUND = "searched_not_found"
    NOT_AVAILABLE = "not_available"
    ACCESS_DENIED = "access_denied"
    NOT_APPLICABLE = "not_applicable"
    UNSUPPORTED = "unsupported"
    FAILED = "failed"


_NEEDS_REASON = {
    ObservationStatus.FAILED,
    ObservationStatus.UNSUPPORTED,
    ObservationStatus.NOT_APPLICABLE,
    ObservationStatus.ACCESS_DENIED,
}


class Observed(Strict, Generic[T]):
    status: ObservationStatus
    value: T | None = None
    source_snapshot_ids: list[str] = Field(default_factory=list)
    reason_code: str | None = None
    retrieval_event_id: str | None = None

    @model_validator(mode="after")
    def _status_consistency(self) -> Self:
        present = self.status is ObservationStatus.VALUE_PRESENT
        if present and self.value is None:
            raise ValueError("status=value_present requires a value")
        if not present and self.value is not None:
            raise ValueError(f"status={self.status.value} must not carry a value")
        if self.status in _NEEDS_REASON and not self.reason_code:
            raise ValueError(f"status={self.status.value} requires reason_code")
        if self.status is ObservationStatus.SEARCHED_NOT_FOUND and not self.retrieval_event_id:
            raise ValueError("searched_not_found requires a logged retrieval_event_id")
        if (
            present or self.status is ObservationStatus.SEARCHED_NOT_FOUND
        ) and not self.source_snapshot_ids:
            raise ValueError(
                f"status={self.status.value} asserts something about a source; "
                "source_snapshot_ids must be nonempty"
            )
        return self

    @property
    def is_present(self) -> bool:
        return self.status is ObservationStatus.VALUE_PRESENT

    @classmethod
    def not_searched(cls) -> "Observed[T]":
        return cls(status=ObservationStatus.NOT_SEARCHED)

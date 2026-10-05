"""Source snapshots and supporting locators (SCHEMA §5)."""

from __future__ import annotations

from typing import Annotated, Literal, Self

from pydantic import Field, model_validator

from tessera.schemas.base import SCHEMA_VERSION, Strict


class SourceSnapshot(Strict):
    schema_version: str = SCHEMA_VERSION
    snapshot_id: str = Field(description="snap_ + SHA-256 of raw content")
    source: str
    source_record_id: str
    source_version: str
    source_uri: str | None = None
    retrieved_at: str
    media_type: Literal["application/json", "text/plain"]
    raw_content_digest: str
    storage_policy: Literal["retained", "metadata_only", "restricted"] = "retained"
    access_status: Literal[
        "available",
        "abstract_only",
        "full_text_unavailable",
        "supplement_unavailable",
        "controlled_access",
        "failed",
    ] = "available"
    license_status: str = "unknown"
    is_synthetic: bool

    @model_validator(mode="after")
    def _addressed_by_content(self) -> Self:
        if self.snapshot_id != "snap_" + self.raw_content_digest:
            raise ValueError("snapshot_id must be snap_ + raw_content_digest")
        return self


class JsonFieldLocator(Strict):
    kind: Literal["json_field"] = "json_field"
    snapshot_id: str
    pointer: str = Field(description="RFC 6901 JSON Pointer")
    expected_value_digest: str


class TextSpanLocator(Strict):
    """Zero-based code-point offsets, start inclusive, end exclusive."""

    kind: Literal["text_span"] = "text_span"
    snapshot_id: str
    text_digest: str
    parser_version: str
    start: int = Field(ge=0)
    end: int = Field(ge=0)
    excerpt: str

    @model_validator(mode="after")
    def _ordered(self) -> Self:
        if self.end <= self.start:
            raise ValueError("text span end must be greater than start")
        return self


SupportingLocator = Annotated[JsonFieldLocator | TextSpanLocator, Field(discriminator="kind")]

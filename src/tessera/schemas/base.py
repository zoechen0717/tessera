"""Strict base model, canonical serialization and content digests.

Every persisted object is serialized through `canonical_json`, so the same
content always produces the same bytes and the same digest. Replay equality
(R-16) is defined on these digests.
"""

from __future__ import annotations

import hashlib
import json
from enum import StrEnum
from typing import Any

from pydantic import BaseModel, ConfigDict

SCHEMA_VERSION = "2.1.0-draft"


class Strict(BaseModel):
    """Unknown fields are an error, not a silent pass-through (SCHEMA §1)."""

    model_config = ConfigDict(extra="forbid", frozen=True, strict=False)


def canonical_json(obj: Any) -> str:
    """UTF-8 JSON with sorted keys and no whitespace."""
    if isinstance(obj, BaseModel):
        obj = obj.model_dump(mode="json")
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def sha256_hex(data: bytes | str) -> str:
    if isinstance(data, str):
        data = data.encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def content_digest(obj: Any) -> str:
    return sha256_hex(canonical_json(obj))


class ValidationState(StrEnum):
    """SCHEMA §2.4."""

    DRAFT = "draft"
    MECHANICALLY_VALID = "mechanically_valid"
    ACCEPTED_AUTO = "accepted_auto"
    ACCEPTED_HUMAN = "accepted_human"
    QUARANTINED = "quarantined"
    REJECTED = "rejected"
    SUPERSEDED = "superseded"

    @property
    def is_accepted(self) -> bool:
        return self in (ValidationState.ACCEPTED_AUTO, ValidationState.ACCEPTED_HUMAN)

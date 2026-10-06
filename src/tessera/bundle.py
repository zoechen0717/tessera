"""Content-addressed run bundle (DEC-13).

A bundle directory holds frozen objects as canonical JSONL. Everything ranking
needs is in the bundle, so replay reads only the bundle: no network, no model,
no live config.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, TypeVar

from pydantic import BaseModel

from tessera.schemas.base import canonical_json, sha256_hex

M = TypeVar("M", bound=BaseModel)

# Frozen inputs to ranking. Their combined digest is the evidence digest.
FROZEN_FILES = (
    "candidate_mentions.jsonl",
    "identity_decisions.jsonl",
    "candidate_variants.jsonl",
    "memberships.jsonl",
    "annotations.jsonl",
    "source_snapshots.jsonl",
    "claims.jsonl",
    "evidence_links.jsonl",
    "validation.jsonl",
    "gene_mechanism.json",
    "observation_search.json",
)
# Derived outputs whose bytes must be identical under replay (R-16).
DERIVED_FILES = ("features.jsonl", "ranked_variants.csv", "strata_report.json", "report.md")


def write_jsonl(path: Path, objs: Iterable[BaseModel | dict]) -> None:
    lines = [canonical_json(o) for o in objs]
    path.write_text("".join(line + "\n" for line in lines), encoding="utf-8")


def read_jsonl(path: Path, model: type[M]) -> list[M]:
    if not path.exists():
        raise FileNotFoundError(f"bundle artifact missing: {path.name}")
    return [model.model_validate_json(line) for line in path.read_text(encoding="utf-8").splitlines() if line]


def write_json(path: Path, obj) -> None:
    path.write_text(canonical_json(obj) + "\n", encoding="utf-8")


def read_json(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"bundle artifact missing: {path.name}")
    return json.loads(path.read_text(encoding="utf-8"))


def file_digest(path: Path) -> str:
    return sha256_hex(path.read_bytes())


def digest_files(directory: Path, names: Iterable[str]) -> dict[str, str]:
    return {n: file_digest(directory / n) for n in names}


def combined_digest(digests: dict[str, str]) -> str:
    return sha256_hex(canonical_json(sorted(digests.items())))

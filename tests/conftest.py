from __future__ import annotations

import json
from pathlib import Path

import pytest

from tessera.pipeline import run
from tessera.synthetic import build


@pytest.fixture
def fixture_input(tmp_path: Path) -> Path:
    return build(tmp_path / "input")


@pytest.fixture
def run_bundle(fixture_input: Path, tmp_path: Path) -> Path:
    out = tmp_path / "run"
    run(fixture_input, out)
    return out


def load_features(bundle: Path) -> dict[str, dict]:
    return {f["label"]: f for f in map(json.loads, (bundle / "features.jsonl").read_text().splitlines())}


def load_ranked(bundle: Path) -> list[dict]:
    import csv

    with (bundle / "ranked_variants.csv").open() as fh:
        return list(csv.DictReader(fh))


def mention_label(bundle: Path, mention_id: str) -> str | None:
    decisions = map(json.loads, (bundle / "identity_decisions.jsonl").read_text().splitlines())
    variants = {
        v["variant_id"]: f"{v['chrom']}:{v['pos_1based']}:{v['ref']}>{v['alt']}"
        for v in map(json.loads, (bundle / "candidate_variants.jsonl").read_text().splitlines())
    }
    for d in decisions:
        if d["mention_id"] == mention_id and d["variant_ids"]:
            return variants[d["variant_ids"][0]]
    return None

"""LLM layer: recording, replay, budgets, schema repair, secrets (no network)."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import pytest
from pydantic import BaseModel

from tessera.llm import (
    Budget, LLMClient, LLMError, ModelOutputStore, RawCompletion, _model_core, load_profile, parse_json,
)

CONFIG = Path(__file__).resolve().parents[1] / "config" / "llm.yaml"


class Out(BaseModel):
    allele: str
    direction: str


@dataclass
class FakeBackend:
    replies: list[str]
    model: str = "anthropic.claude-opus-4-8"
    calls: list = field(default_factory=list)

    def complete(self, system, user):
        self.calls.append(user)
        return RawCompletion(self.replies[len(self.calls) - 1], self.model, {"input_tokens": 10, "output_tokens": 5}, "stop")


def _client(tmp_path, replies, mode="live", budget=None, **kw):
    os.environ.setdefault("TESSERA_LLM_BASE_URL", "https://example.invalid/v1")
    prof = load_profile(CONFIG, "claude_aws_proxy")
    return LLMClient(prof, ModelOutputStore(tmp_path / "mo"), budget or Budget(10, 10_000, 10_000),
                     mode=mode, backend=FakeBackend(replies, **kw) if mode == "live" else None)


def test_fenced_json_is_parsed():
    assert parse_json('```json\n{"a": 1}\n```') == {"a": 1}


def test_record_then_replay_without_backend(tmp_path):
    c = _client(tmp_path, ['{"allele": "p.Arg990Ter", "direction": "reduced"}'])
    r = c.structured(task_id="t1", prompt_version="v1", system="s", user="u", schema=Out)
    assert r.value.allele == "p.Arg990Ter" and not r.outputs[0].model_mismatch
    replay = _client(tmp_path, [], mode="replay")
    assert replay.structured(task_id="t1", prompt_version="v1", system="s", user="u", schema=Out).value == r.value


def test_replay_fails_when_record_missing(tmp_path):
    with pytest.raises(LLMError, match="snapshot_missing"):
        _client(tmp_path, [], mode="replay").structured(task_id="t", prompt_version="v1", system="s", user="new", schema=Out)


def test_one_schema_repair_then_failure_is_reported(tmp_path):
    c = _client(tmp_path, ["not json", '{"allele": "x"}'])
    r = c.structured(task_id="t", prompt_version="v1", system="s", user="u", schema=Out)
    assert r.value is None and r.error.startswith("schema_failure") and len(r.outputs) == 2


def test_repair_counts_against_budget(tmp_path):
    c = _client(tmp_path, ["bad", "bad"], budget=Budget(1, 10_000, 10_000))
    with pytest.raises(LLMError, match="budget_exhausted"):
        c.structured(task_id="t", prompt_version="v1", system="s", user="u", schema=Out)


def test_returned_model_mismatch_is_flagged(tmp_path):
    c = _client(tmp_path, ['{"allele": "a", "direction": "b"}'], model="claude-haiku-4-5")
    assert c.structured(task_id="t", prompt_version="v1", system="s", user="u", schema=Out).outputs[0].model_mismatch


def test_api_key_never_recorded(tmp_path, monkeypatch):
    monkeypatch.setenv("TESSERA_LLM_API_KEY", "sk-SECRET-should-not-appear")
    c = _client(tmp_path, ['{"allele": "a", "direction": "b"}'])
    c.structured(task_id="t", prompt_version="v1", system="s", user="u", schema=Out)
    for p in (tmp_path / "mo").glob("*.json"):
        assert "sk-SECRET" not in p.read_text()


def test_gateway_model_names_compare_equal():
    assert _model_core("[AWS]claude-opus-4-8") == _model_core("anthropic.claude-opus-4-8")
    assert _model_core("[k-逆向按次-8k]claude-opus-5-5") == "claude-opus-5-5"

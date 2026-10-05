"""Replay: R-16, INV-13, INV-14."""

from __future__ import annotations

import json
import socket

import pytest

from tessera.pipeline import replay


def test_replay_is_byte_identical_without_network(run_bundle, tmp_path, monkeypatch):
    def no_network(*a, **k):
        raise AssertionError("replay attempted a network connection")

    monkeypatch.setattr(socket, "socket", no_network)
    monkeypatch.setattr(socket, "create_connection", no_network)
    result, same = replay(run_bundle, tmp_path / "replay")
    assert same
    for name, digest in result.artifact_digests.items():
        assert json.loads((run_bundle / "manifest.json").read_text())["derived_artifacts"][name] == digest


def test_replay_fails_on_modified_frozen_bundle(run_bundle, tmp_path):
    p = run_bundle / "annotations.jsonl"
    p.write_text(p.read_text().replace("0.91", "0.99", 1))
    with pytest.raises(ValueError, match="INV-13"):
        replay(run_bundle, tmp_path / "replay")


def test_replay_fails_on_missing_artifact(run_bundle, tmp_path):
    (run_bundle / "claims.jsonl").unlink()
    with pytest.raises(FileNotFoundError):
        replay(run_bundle, tmp_path / "replay")


def test_manifest_references_policy_and_evidence(run_bundle):
    m = json.loads((run_bundle / "manifest.json").read_text())
    assert m["policy_status"] == "provisional" and m["evidence_digest"] and m["policy_digest"]
    assert m["network_calls"] == 0 and m["llm_calls"] == 0 and m["is_synthetic"] is True
    for f in map(json.loads, (run_bundle / "features.jsonl").read_text().splitlines()):
        assert f["evidence_digest"] == m["evidence_digest"] and f["policy_digest"] == m["policy_digest"]

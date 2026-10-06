"""Cross-run request cache (no network: urlopen is replaced)."""

from __future__ import annotations

import io
import json

import pytest

from tessera.sources import http as H


class _Resp(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


@pytest.fixture
def calls(monkeypatch):
    seen = []

    def fake(req, timeout=None):
        seen.append(req.full_url)
        return _Resp(json.dumps({"n": len(seen)}).encode())

    monkeypatch.setattr(H.urllib.request, "urlopen", fake)
    return seen


def test_ttl_hit_and_live_search(tmp_path, calls):
    f = H.Fetcher(tmp_path / "s", cache_dir=tmp_path / "c")
    a = f.get("src", "r", "v1", "https://x.invalid/a", ttl_days=30)
    b = f.get("src", "r", "v1", "https://x.invalid/a", ttl_days=30)
    assert a.content == b.content and len(calls) == 1 and f.cache_hits == 1
    f.get("src", "q", "v1", "https://x.invalid/search")       # no ttl: always live
    f.get("src", "q", "v1", "https://x.invalid/search")
    assert len(calls) == 3


def test_new_source_version_misses_cache(tmp_path, calls):
    f = H.Fetcher(tmp_path / "s", cache_dir=tmp_path / "c")
    f.get("ensembl", "r", "116", "https://x.invalid/vep", ttl_days=float("inf"))
    f.get("ensembl", "r", "117", "https://x.invalid/vep", ttl_days=float("inf"))
    assert len(calls) == 2


def test_cache_hit_keeps_original_retrieval_time_and_index_merges(tmp_path, calls):
    f1 = H.Fetcher(tmp_path / "s", cache_dir=tmp_path / "c")
    f1.get("src", "r1", "v", "https://x.invalid/one", ttl_days=30)
    f1.write_index()
    first = json.loads((tmp_path / "s" / "snapshots.jsonl").read_text().splitlines()[0])
    f2 = H.Fetcher(tmp_path / "s", cache_dir=tmp_path / "c")
    f2.get("src", "r1", "v", "https://x.invalid/one", ttl_days=30)
    f2.get("src", "r2", "v", "https://x.invalid/two", ttl_days=30)
    f2.write_index()
    idx = [json.loads(l) for l in (tmp_path / "s" / "snapshots.jsonl").read_text().splitlines()]
    assert len(idx) == 2
    assert next(m for m in idx if m["source_record_id"] == "r1")["retrieved_at"] == first["retrieved_at"]

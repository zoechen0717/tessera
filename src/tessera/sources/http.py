"""HTTP client that records every response as a content-addressed snapshot.

Only adapters in tessera.sources use the network, and only through this class.
Freeze/derive/replay never import it. Retries cover transient failures only,
with bounded backoff (ARCH §11); every physical attempt is counted.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from tessera.schemas.base import sha256_hex

USER_AGENT = "tessera/0.2 (allele nomination research tool)"
_TRANSIENT = {429, 500, 502, 503, 504}


class SourceError(Exception):
    def __init__(self, source: str, code: str, detail: str):
        super().__init__(f"{source}: {code}: {detail}")
        self.source, self.code, self.detail = source, code, detail


@dataclass
class Response:
    content: bytes
    snapshot_id: str

    def json(self):
        return json.loads(self.content)

    def text(self) -> str:
        return self.content.decode("utf-8")


@dataclass
class Fetcher:
    snapshot_dir: Path
    timeout: float = 60.0
    max_retries: int = 3
    min_interval: dict[str, float] = field(default_factory=dict)
    cache_dir: Path | None = None
    attempts: int = 0
    logical_requests: int = 0
    cache_hits: int = 0
    _last: dict[str, float] = field(default_factory=dict)
    _meta: list[dict] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.snapshot_dir.mkdir(parents=True, exist_ok=True)
        if self.cache_dir:
            (self.cache_dir / "requests").mkdir(parents=True, exist_ok=True)
            (self.cache_dir / "blobs").mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------- requests

    def get(self, source: str, record_id: str, version: str, url: str, params: dict | None = None,
            media: str = "application/json", ttl_days: float | None = None) -> Response:
        """`ttl_days`: reuse a cached response younger than this (math.inf = immutable
        content such as an archived full text). None = always fetch (searches)."""
        full = url + ("?" + urllib.parse.urlencode(params) if params else "")
        req = urllib.request.Request(full, headers={"User-Agent": USER_AGENT, "Accept": media if media == "application/json" else "*/*"})
        return self._send(source, record_id, version, req, {"method": "GET", "url": full}, media, ttl_days)

    def post_json(self, source: str, record_id: str, version: str, url: str, payload: dict,
                  ttl_days: float | None = None) -> Response:
        body = json.dumps(payload, sort_keys=True).encode()
        req = urllib.request.Request(
            url, data=body, method="POST",
            headers={"User-Agent": USER_AGENT, "Content-Type": "application/json", "Accept": "application/json"},
        )
        request = {"method": "POST", "url": url, "body_sha256": sha256_hex(body)}
        return self._send(source, record_id, version, req, request, "application/json", ttl_days)

    # ----------------------------------------------------------------- cache

    def _cache_key(self, request_meta: dict) -> str:
        return sha256_hex(json.dumps(request_meta, sort_keys=True))  # includes source/version, see _send

    def _cache_get(self, request_meta: dict, ttl_days: float) -> tuple[bytes, dict] | None:
        idx = self.cache_dir / "requests" / f"{self._cache_key(request_meta)}.json"
        if not idx.exists():
            return None
        meta = json.loads(idx.read_text())
        age = (datetime.now(timezone.utc) - datetime.fromisoformat(meta["retrieved_at"])).total_seconds() / 86400
        blob = self.cache_dir / "blobs" / meta["digest"]
        if age > ttl_days or not blob.exists():
            return None
        content = blob.read_bytes()
        return (content, meta) if sha256_hex(content) == meta["digest"] else None

    def _cache_put(self, request_meta: dict, content: bytes, retrieved_at: str) -> None:
        digest = sha256_hex(content)
        (self.cache_dir / "blobs" / digest).write_bytes(content)
        (self.cache_dir / "requests" / f"{self._cache_key(request_meta)}.json").write_text(
            json.dumps({"digest": digest, "retrieved_at": retrieved_at, "url": request_meta["url"]}))

    def _send(self, source, record_id, version, req, request_meta, media, ttl_days=None) -> Response:
        self.logical_requests += 1
        # a new source release (e.g. Ensembl 117) must not be served from an older cache entry
        request_meta = {**request_meta, "_source": source, "_version": version}
        if self.cache_dir and ttl_days is not None:
            hit = self._cache_get(request_meta, ttl_days)
            if hit:
                self.cache_hits += 1
                return self._record(source, record_id, version, request_meta, media, hit[0], hit[1]["retrieved_at"])
        wait = self.min_interval.get(source, 0.0) - (time.monotonic() - self._last.get(source, 0.0))
        if wait > 0:
            time.sleep(wait)
        last_err = ""
        for attempt in range(self.max_retries + 1):
            self.attempts += 1
            try:
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    content = resp.read()
                self._last[source] = time.monotonic()
                now = datetime.now(timezone.utc).isoformat(timespec="seconds")
                if self.cache_dir:
                    self._cache_put(request_meta, content, now)
                return self._record(source, record_id, version, request_meta, media, content, now)
            except urllib.error.HTTPError as e:
                last_err = f"HTTP {e.code}"
                if e.code not in _TRANSIENT:
                    raise SourceError(source, "http_error", f"{last_err} for {request_meta['url']}") from e
                retry_after = e.headers.get("Retry-After") if e.headers else None
                delay = min(30.0, float(retry_after)) if retry_after and retry_after.isdigit() else min(30.0, 2.0**attempt)
            except (urllib.error.URLError, TimeoutError) as e:
                last_err = f"transport: {e}"
                delay = min(30.0, 2.0**attempt)
            if attempt < self.max_retries:
                time.sleep(delay)
        raise SourceError(source, "transient_network", f"{last_err} after {self.max_retries + 1} attempts")

    # ------------------------------------------------------------- snapshots

    def _record(self, source, record_id, version, request_meta, media, content: bytes, retrieved_at: str) -> Response:
        digest = sha256_hex(content)
        snapshot_id = "snap_" + digest
        ext = {"application/json": ".json", "application/zip": ".zip"}.get(media, ".txt")
        fname = f"{source}_{digest[:16]}{ext}"
        (self.snapshot_dir / fname).write_bytes(content)
        if not any(m["snapshot_id"] == snapshot_id for m in self._meta):
            self._meta.append({
                "file": fname,
                "snapshot_id": snapshot_id,
                "source": source,
                "source_record_id": record_id,
                "source_version": version,
                "source_uri": request_meta["url"],
                "retrieved_at": retrieved_at,  # original retrieval time, also for cache hits
                "media_type": media if media in ("application/json", "application/zip") else "text/plain",
                "raw_content_digest": digest,
                "license_status": "unknown",
                "is_synthetic": False,
            })
        return Response(content=content, snapshot_id=snapshot_id)

    def seed_from_snapshots(self, snapshot_dir: Path, record_prefixes: tuple[str, ...]) -> int:
        """Populate the cache from an earlier run's GET snapshots (same source/version)."""
        if not self.cache_dir or not (snapshot_dir / "snapshots.jsonl").exists():
            return 0
        n = 0
        for line in (snapshot_dir / "snapshots.jsonl").read_text().splitlines():
            m = json.loads(line)
            if not m["source_record_id"].startswith(record_prefixes) or not m.get("source_uri"):
                continue
            content = (snapshot_dir / m["file"]).read_bytes()
            if sha256_hex(content) != m["raw_content_digest"]:
                continue
            meta = {"method": "GET", "url": m["source_uri"], "_source": m["source"], "_version": m["source_version"]}
            self._cache_put(meta, content, m["retrieved_at"])
            n += 1
        return n

    def write_index(self) -> None:
        """Merge with an existing index so reusing a directory never drops earlier snapshots."""
        idx = self.snapshot_dir / "snapshots.jsonl"
        merged = {}
        if idx.exists():
            for line in idx.read_text().splitlines():
                m = json.loads(line)
                if (self.snapshot_dir / m["file"]).exists():
                    merged[m["snapshot_id"]] = m
        for m in self._meta:
            merged.setdefault(m["snapshot_id"], m)
        lines = sorted(merged.values(), key=lambda m: m["snapshot_id"])
        (self.snapshot_dir / "snapshots.jsonl").write_text(
            "".join(json.dumps(m, sort_keys=True) + "\n" for m in lines)
        )

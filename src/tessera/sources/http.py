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
    attempts: int = 0
    logical_requests: int = 0
    _last: dict[str, float] = field(default_factory=dict)
    _meta: list[dict] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.snapshot_dir.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------- requests

    def get(self, source: str, record_id: str, version: str, url: str, params: dict | None = None,
            media: str = "application/json") -> Response:
        full = url + ("?" + urllib.parse.urlencode(params) if params else "")
        req = urllib.request.Request(full, headers={"User-Agent": USER_AGENT, "Accept": media})
        return self._send(source, record_id, version, req, {"method": "GET", "url": full}, media)

    def post_json(self, source: str, record_id: str, version: str, url: str, payload: dict) -> Response:
        body = json.dumps(payload, sort_keys=True).encode()
        req = urllib.request.Request(
            url, data=body, method="POST",
            headers={"User-Agent": USER_AGENT, "Content-Type": "application/json", "Accept": "application/json"},
        )
        request = {"method": "POST", "url": url, "body_sha256": sha256_hex(body)}
        return self._send(source, record_id, version, req, request, "application/json")

    def _send(self, source, record_id, version, req, request_meta, media) -> Response:
        self.logical_requests += 1
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
                return self._record(source, record_id, version, request_meta, media, content)
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

    def _record(self, source, record_id, version, request_meta, media, content: bytes) -> Response:
        digest = sha256_hex(content)
        snapshot_id = "snap_" + digest
        fname = f"{source}_{digest[:16]}" + (".json" if media == "application/json" else ".txt")
        (self.snapshot_dir / fname).write_bytes(content)
        if not any(m["snapshot_id"] == snapshot_id for m in self._meta):
            self._meta.append({
                "file": fname,
                "snapshot_id": snapshot_id,
                "source": source,
                "source_record_id": record_id,
                "source_version": version,
                "source_uri": request_meta["url"],
                "retrieved_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "media_type": "application/json" if media == "application/json" else "text/plain",
                "raw_content_digest": digest,
                "license_status": "unknown",
                "is_synthetic": False,
            })
        return Response(content=content, snapshot_id=snapshot_id)

    def write_index(self) -> None:
        lines = sorted(self._meta, key=lambda m: m["snapshot_id"])
        (self.snapshot_dir / "snapshots.jsonl").write_text(
            "".join(json.dumps(m, sort_keys=True) + "\n" for m in lines)
        )

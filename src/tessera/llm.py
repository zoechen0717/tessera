"""Provider-neutral LLM client for bounded semantic workers (SPEC §6, ARCH §6, §10).

- Providers are user-configured profiles (config/llm.yaml); keys come from the
  environment or a git-ignored .env and are never written to outputs or logs.
- Every call is recorded as a content-addressed model output (request digest,
  requested and returned model, usage, raw text). `ReplayClient` serves those
  records without any network access, so extraction can be replayed exactly
  and re-extraction variability measured separately (R-16).
- Host budgets are reserved before dispatch, so retries and schema repairs
  cannot exceed them.
- Output is parsed as JSON (code fences stripped) and validated against a
  Pydantic model; one repair attempt is allowed. A model is never trusted to
  have followed the format.
"""

from __future__ import annotations

import json
import os
import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Protocol, TypeVar

import yaml
from pydantic import BaseModel, ValidationError

from tessera.schemas.base import canonical_json, sha256_hex

T = TypeVar("T", bound=BaseModel)
_TRANSIENT = {408, 429, 500, 502, 503, 504, 529}


class LLMError(Exception):
    def __init__(self, code: str, detail: str):
        super().__init__(f"{code}: {detail}")
        self.code, self.detail = code, detail


# ------------------------------------------------------------------ config


def load_dotenv(path: Path) -> None:
    """Populate os.environ from KEY=VALUE lines without overriding existing vars."""
    if not path.exists():
        return
    for line in path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#") and "=" in line:
            k, v = line.split("=", 1)
            os.environ.setdefault(k.strip(), v.strip())


@dataclass(frozen=True)
class Profile:
    name: str
    api_format: str
    model: str
    api_key_env: str
    base_url: str | None
    max_tokens: int
    timeout_s: float
    effort: str | None = None

    @property
    def public(self) -> dict:
        """What may be recorded: never the key, never the key's value."""
        return {"profile": self.name, "api_format": self.api_format, "model": self.model,
                "base_url": self.base_url, "max_tokens": self.max_tokens, "effort": self.effort}


def load_profile(config_path: Path, name: str | None = None, env_file: Path | None = None) -> Profile:
    if env_file:
        load_dotenv(env_file)
    cfg = yaml.safe_load(config_path.read_text())
    name = name or cfg["default_profile"]
    p = cfg["profiles"][name]
    if p["api_format"] not in ("openai_compatible", "anthropic"):
        raise LLMError("invalid_config", f"unknown api_format {p['api_format']}")
    base = p.get("base_url") or (os.environ.get(p["base_url_env"]) if p.get("base_url_env") else None)
    if p["api_format"] == "openai_compatible" and not base:
        raise LLMError("invalid_config", f"profile {name} needs base_url or base_url_env")
    return Profile(name=name, api_format=p["api_format"], model=p["model"], api_key_env=p["api_key_env"],
                   base_url=base.rstrip("/") if base else None, max_tokens=int(p.get("max_tokens", 8000)),
                   timeout_s=float(p.get("timeout_s", 180)), effort=p.get("effort"))


# ------------------------------------------------------------------ budget


@dataclass
class Budget:
    max_requests: int
    max_input_tokens: int
    max_output_tokens: int
    requests: int = 0
    input_tokens: int = 0
    output_tokens: int = 0

    def reserve(self) -> None:
        if self.requests >= self.max_requests:
            raise LLMError("budget_exhausted", f"max_requests={self.max_requests}")
        if self.input_tokens >= self.max_input_tokens or self.output_tokens >= self.max_output_tokens:
            raise LLMError("budget_exhausted", "token budget reached")
        self.requests += 1

    def charge(self, usage: dict) -> None:
        self.input_tokens += int(usage.get("input_tokens") or 0)
        self.output_tokens += int(usage.get("output_tokens") or 0)


# ------------------------------------------------------------------ records


class ModelOutput(BaseModel):
    """One recorded model call. Contains no credentials."""

    output_id: str
    request_digest: str
    task_id: str
    prompt_version: str
    schema_name: str
    provider: dict
    requested_model: str
    returned_model: str | None
    model_mismatch: bool
    usage: dict
    finish_reason: str | None
    text: str
    created_at: str
    attempt: int


@dataclass
class ModelOutputStore:
    directory: Path
    _index: dict[str, ModelOutput] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.directory.mkdir(parents=True, exist_ok=True)
        for p in sorted(self.directory.glob("*.json")):
            rec = ModelOutput.model_validate_json(p.read_text())
            self._index[f"{rec.request_digest}:{rec.attempt}"] = rec

    def put(self, rec: ModelOutput) -> None:
        (self.directory / f"{rec.output_id}.json").write_text(canonical_json(rec) + "\n")
        self._index[f"{rec.request_digest}:{rec.attempt}"] = rec

    def get(self, request_digest: str, attempt: int) -> ModelOutput | None:
        return self._index.get(f"{request_digest}:{attempt}")


# ------------------------------------------------------------------ backends


@dataclass
class RawCompletion:
    text: str
    returned_model: str | None
    usage: dict
    finish_reason: str | None


class Backend(Protocol):
    def complete(self, system: str, user: str) -> RawCompletion: ...


@dataclass
class OpenAICompatibleBackend:
    profile: Profile
    max_retries: int = 3

    def complete(self, system: str, user: str) -> RawCompletion:
        key = os.environ.get(self.profile.api_key_env)
        if not key:
            raise LLMError("missing_credentials", f"environment variable {self.profile.api_key_env} is not set")
        body = {"model": self.profile.model, "max_tokens": self.profile.max_tokens,
                "messages": [{"role": "system", "content": system}, {"role": "user", "content": user}]}
        req = urllib.request.Request(
            f"{self.profile.base_url}/chat/completions", data=json.dumps(body).encode(), method="POST",
            headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json",
                     "User-Agent": "tessera/0.2"})
        last = ""
        for attempt in range(self.max_retries + 1):
            try:
                with urllib.request.urlopen(req, timeout=self.profile.timeout_s) as resp:
                    d = json.load(resp)
                choice = d["choices"][0]
                u = d.get("usage") or {}
                return RawCompletion(
                    text=choice["message"].get("content") or "",
                    returned_model=d.get("model"),
                    usage={"input_tokens": u.get("prompt_tokens"), "output_tokens": u.get("completion_tokens")},
                    finish_reason=choice.get("finish_reason"),
                )
            except urllib.error.HTTPError as e:
                last = f"HTTP {e.code}"
                body = e.read()[:300].decode("utf-8", errors="replace")
                if e.code in (401, 402, 403):
                    # credentials, quota or balance: retrying cannot help, and every
                    # further request would fail the same way
                    raise LLMError("provider_quota_or_auth", f"{last}: {body}") from e
                if e.code not in _TRANSIENT:
                    raise LLMError("provider_error", f"{last}: {body}") from e
            except (urllib.error.URLError, TimeoutError) as e:
                last = f"transport: {e}"
            if attempt < self.max_retries:
                time.sleep(min(30.0, 2.0**attempt))
        raise LLMError("transient_network", f"{last} after {self.max_retries + 1} attempts")


@dataclass
class AnthropicBackend:
    """Anthropic Messages API through the official SDK (imported lazily)."""

    profile: Profile

    def complete(self, system: str, user: str) -> RawCompletion:
        try:
            import anthropic
        except ImportError as e:
            raise LLMError("missing_dependency", "pip install anthropic to use api_format: anthropic") from e
        client = anthropic.Anthropic(api_key=os.environ.get(self.profile.api_key_env),
                                     timeout=self.profile.timeout_s)
        kwargs = {}
        if self.profile.effort:
            kwargs["output_config"] = {"effort": self.profile.effort}
        try:
            with client.messages.stream(model=self.profile.model, max_tokens=self.profile.max_tokens,
                                        thinking={"type": "adaptive"}, system=system,
                                        messages=[{"role": "user", "content": user}], **kwargs) as stream:
                msg = stream.get_final_message()
        except anthropic.APIStatusError as e:
            raise LLMError("provider_error", f"HTTP {e.status_code}") from e
        except anthropic.APIConnectionError as e:
            raise LLMError("transient_network", str(e)) from e
        if msg.stop_reason == "refusal":
            raise LLMError("refusal", "model declined the request")
        text = "".join(b.text for b in msg.content if b.type == "text")
        return RawCompletion(text=text, returned_model=msg.model,
                             usage={"input_tokens": msg.usage.input_tokens, "output_tokens": msg.usage.output_tokens},
                             finish_reason=msg.stop_reason)


def backend_for(profile: Profile) -> Backend:
    return OpenAICompatibleBackend(profile) if profile.api_format == "openai_compatible" else AnthropicBackend(profile)


# ------------------------------------------------------------------ client


_FENCE = re.compile(r"^\s*```(?:json)?\s*(.*?)\s*```\s*$", re.S)


def parse_json(text: str) -> object:
    m = _FENCE.match(text)
    return json.loads(m.group(1) if m else text.strip())


@dataclass
class StructuredResult:
    value: BaseModel | None
    outputs: list[ModelOutput]
    error: str | None


@dataclass
class LLMClient:
    """mode='live' calls the backend and records; mode='replay' only reads records."""

    profile: Profile
    store: ModelOutputStore
    budget: Budget
    mode: str = "live"
    backend: Backend | None = None
    max_schema_repairs: int = 1

    def __post_init__(self) -> None:
        if self.mode not in ("live", "replay"):
            raise ValueError("mode must be live or replay")
        if self.mode == "live" and self.backend is None:
            self.backend = backend_for(self.profile)

    def structured(self, *, task_id: str, prompt_version: str, system: str, user: str, schema: type[T]) -> StructuredResult:
        outputs: list[ModelOutput] = []
        prompt = user
        last_err = None
        for attempt in range(self.max_schema_repairs + 1):
            digest = sha256_hex(canonical_json({"provider": self.profile.public, "prompt_version": prompt_version,
                                                "schema": schema.__name__, "system": system, "user": prompt}))
            rec = self.store.get(digest, attempt)
            if rec is None:
                if self.mode == "replay":
                    raise LLMError("snapshot_missing", f"no recorded output for task {task_id} attempt {attempt}")
                self.budget.reserve()
                raw = self.backend.complete(system, prompt)
                self.budget.charge(raw.usage)
                returned = raw.returned_model
                rec = ModelOutput(
                    output_id="mo_" + sha256_hex(f"{digest}:{attempt}:{raw.text}")[:32],
                    request_digest=digest, task_id=task_id, prompt_version=prompt_version,
                    schema_name=schema.__name__, provider=self.profile.public,
                    requested_model=self.profile.model, returned_model=returned,
                    model_mismatch=bool(returned) and _model_core(returned) != _model_core(self.profile.model),
                    usage=raw.usage, finish_reason=raw.finish_reason, text=raw.text,
                    created_at=datetime.now(timezone.utc).isoformat(timespec="seconds"), attempt=attempt)
                self.store.put(rec)
            outputs.append(rec)
            try:
                return StructuredResult(value=schema.model_validate(parse_json(rec.text)), outputs=outputs, error=None)
            except (json.JSONDecodeError, ValidationError) as e:
                last_err = f"{type(e).__name__}: {str(e)[:400]}"
                prompt = (f"{user}\n\nYour previous reply could not be used ({last_err}). "
                          "Reply again with only the JSON object, matching the schema exactly.")
        return StructuredResult(value=None, outputs=outputs, error=f"schema_failure: {last_err}")


def _model_core(name: str) -> str:
    """Compare model names across gateway decorations: '[AWS]claude-opus-4-8' ~ 'anthropic.claude-opus-4-8'."""
    name = re.sub(r"^\[[^\]]*\]", "", name).split("/")[-1]
    return re.sub(r"^(anthropic|us|eu|global)\.", "", name).strip()

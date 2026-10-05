"""Absolute anchored scales (spec 19.6).

Every raw input consumed by any score component is mapped to [0, 1] by an entry
in anchors.yaml. No component may use min-max, z-score, quantile, or rank
normalization across the candidate set.

Why this is a correctness requirement and not a style preference: under
candidate-set-dependent normalization, adding one variant to the candidate
universe changes every other variant's score. That makes the ranking a function
of candidate-set composition, which violates Guarantee 4 (identical snapshot ->
identical score) and makes the "known variants retained" regression check
unstable for reasons unrelated to any code change.

A missing anchor raises at config load. It never defaults.
"""

from __future__ import annotations

import bisect
import math
from pathlib import Path
from typing import Any

import yaml
from pydantic import Field, model_validator

from tessera.enums import AnchorType
from tessera.schemas.common import Strict


class MissingAnchorError(KeyError):
    """Raised when a raw input has no anchor entry. Never downgraded to a default."""


class InvalidAnchorError(ValueError):
    pass


class Anchor(Strict):
    """One raw input -> [0, 1] mapping."""

    name: str
    type: AnchorType
    note: str | None = None

    # categorical / boolean
    map: dict[str, float] | None = None
    unmapped: float | None = Field(
        default=None,
        description="Value for a category absent from `map`. None = raise.",
    )

    # piecewise
    breakpoints: list[float] | None = None
    values: list[float] | None = None
    log_scale: bool = False

    # linear_clamp / identity
    lo: float = 0.0
    hi: float = 1.0
    invert: bool = False

    @model_validator(mode="after")
    def _shape_matches_type(self) -> "Anchor":
        t = self.type
        if t is AnchorType.CATEGORICAL or t is AnchorType.BOOLEAN:
            if not self.map:
                raise InvalidAnchorError(f"anchor {self.name!r}: {t.value} requires `map`")
        elif t is AnchorType.PIECEWISE:
            if not self.breakpoints or not self.values:
                raise InvalidAnchorError(
                    f"anchor {self.name!r}: piecewise requires `breakpoints` and `values`"
                )
            if len(self.breakpoints) != len(self.values):
                raise InvalidAnchorError(
                    f"anchor {self.name!r}: breakpoints ({len(self.breakpoints)}) and "
                    f"values ({len(self.values)}) must be the same length"
                )
            if list(self.breakpoints) != sorted(self.breakpoints):
                raise InvalidAnchorError(f"anchor {self.name!r}: breakpoints must ascend")
        elif t in (AnchorType.LINEAR_CLAMP, AnchorType.IDENTITY):
            if self.hi == self.lo:
                raise InvalidAnchorError(f"anchor {self.name!r}: hi must differ from lo")
        for v in (self.values or []):
            if not 0.0 <= v <= 1.0:
                raise InvalidAnchorError(f"anchor {self.name!r}: value {v} outside [0,1]")
        for v in (self.map or {}).values():
            if not 0.0 <= v <= 1.0:
                raise InvalidAnchorError(f"anchor {self.name!r}: map value {v} outside [0,1]")
        return self

    # --- application --------------------------------------------------------

    def apply(self, raw: Any) -> float:
        if raw is None:
            raise InvalidAnchorError(
                f"anchor {self.name!r}: got None. Absence is a Missingness state, "
                f"not a raw value -- do not route it through an anchor."
            )
        out = self._apply_inner(raw)
        if self.invert:
            out = 1.0 - out
        return min(1.0, max(0.0, out))

    def _apply_inner(self, raw: Any) -> float:
        match self.type:
            case AnchorType.BOOLEAN | AnchorType.CATEGORICAL:
                key = self._normalize_key(raw)
                if key in self.map:  # type: ignore[operator]
                    return self.map[key]  # type: ignore[index]
                if self.unmapped is not None:
                    return self.unmapped
                raise MissingAnchorError(
                    f"anchor {self.name!r}: category {key!r} is not in the map and no "
                    f"`unmapped` fallback is declared. Add the category explicitly -- "
                    f"an undeclared category must not silently become a number."
                )
            case AnchorType.IDENTITY:
                return float(raw)
            case AnchorType.LINEAR_CLAMP:
                return (float(raw) - self.lo) / (self.hi - self.lo)
            case AnchorType.PIECEWISE:
                return self._piecewise(float(raw))
        raise InvalidAnchorError(f"anchor {self.name!r}: unhandled type {self.type}")

    @staticmethod
    def _normalize_key(raw: Any) -> str:
        if isinstance(raw, bool):
            return "true" if raw else "false"
        if isinstance(raw, float) and raw.is_integer():
            return str(int(raw))
        return str(raw).strip().lower().replace(" ", "_").replace("/", "_")

    def _piecewise(self, x: float) -> float:
        bps, vals = self.breakpoints, self.values
        assert bps is not None and vals is not None
        if x <= bps[0]:
            return vals[0]
        if x >= bps[-1]:
            return vals[-1]
        i = bisect.bisect_right(bps, x) - 1
        x0, x1 = bps[i], bps[i + 1]
        y0, y1 = vals[i], vals[i + 1]
        if self.log_scale:
            floor = 1e-12
            lx, lx0, lx1 = (math.log10(max(v, floor)) for v in (x, x0, x1))
            frac = 0.0 if lx1 == lx0 else (lx - lx0) / (lx1 - lx0)
        else:
            frac = 0.0 if x1 == x0 else (x - x0) / (x1 - x0)
        return y0 + frac * (y1 - y0)


class AnchorTable(Strict):
    """The loaded anchors.yaml. Lookup is strict by construction."""

    version: str
    anchors: dict[str, Anchor]

    @classmethod
    def load(cls, path: str | Path) -> "AnchorTable":
        raw = yaml.safe_load(Path(path).read_text())
        version = raw.pop("version", None)
        if version is None:
            raise InvalidAnchorError(f"{path}: anchors file must declare `version`")
        anchors = {
            name: Anchor(name=name, **(spec or {})) for name, spec in raw.get("anchors", {}).items()
        }
        if not anchors:
            raise InvalidAnchorError(f"{path}: no anchors defined")
        return cls(version=str(version), anchors=anchors)

    def __getitem__(self, name: str) -> Anchor:
        try:
            return self.anchors[name]
        except KeyError:
            raise MissingAnchorError(
                f"no anchor for raw input {name!r}. Every raw input mapped into a score "
                f"component must have an entry in anchors.yaml (spec 19.6). Declared: "
                f"{sorted(self.anchors)}"
            ) from None

    def apply(self, name: str, raw: Any) -> float:
        return self[name].apply(raw)

    def require(self, names: list[str]) -> None:
        """Assert every listed raw input is anchored. Called at config load so a
        gap is a startup failure rather than a wrong number much later."""
        missing = [n for n in names if n not in self.anchors]
        if missing:
            raise MissingAnchorError(
                f"scoring config references {len(missing)} unanchored raw input(s): "
                f"{missing}. Add them to anchors.yaml with an explicit rationale."
            )

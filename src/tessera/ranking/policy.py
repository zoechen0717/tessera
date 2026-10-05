"""Load and validate a nomination ranking policy.

Everything checkable is checked at load, so a malformed policy fails before any
row is ranked rather than producing a plausible-looking order.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from tessera.schemas.annotation import Mechanism
from tessera.schemas.base import content_digest
from tessera.schemas.features import Concordance

KEYS = {"K", "O", "P"}
# Predictors known to contain another predictor; never both as keys (DEC-12).
_CONTAINS = {"avi": {"alphamissense"}}


class PolicyError(ValueError):
    pass


@dataclass(frozen=True)
class SubclassRule:
    feature: str
    values: dict[str, Concordance]
    missing: Concordance


@dataclass(frozen=True)
class NominationPolicy:
    policy_id: str
    status: str
    key_order: tuple[str, ...]
    characterized_policy: str
    control_roles: frozenset[str]
    stratum_rules: tuple[tuple[str, frozenset[str]], ...]
    default_stratum: str
    predictors: dict[str, str | None]
    table_id: str
    stratum_level: dict[Mechanism, dict[str, Concordance]]
    subclass: dict[Mechanism, dict[str, SubclassRule]]
    digest: str
    raw: dict[str, Any]

    @property
    def strata(self) -> list[str]:
        out: list[str] = []
        for s, _ in self.stratum_rules:
            if s not in out:
                out.append(s)
        if self.default_stratum not in out:
            out.append(self.default_stratum)
        return out

    def stratum_for(self, terms: list[str]) -> str:
        present = set(terms)
        for stratum, rule_terms in self.stratum_rules:
            if present & rule_terms:
                return stratum
        return self.default_stratum

    @classmethod
    def load(cls, path: Path) -> "NominationPolicy":
        raw = yaml.safe_load(path.read_text())
        table_path = path.parent / raw["concordance_table"]
        table = yaml.safe_load(table_path.read_text())
        return cls.from_dicts(raw, table)

    @classmethod
    def from_dicts(cls, raw: dict[str, Any], table: dict[str, Any]) -> "NominationPolicy":
        key_order = tuple(raw["key_order"])
        if set(key_order) - KEYS or len(set(key_order)) != len(key_order):
            raise PolicyError(f"key_order must be a permutation of a subset of {sorted(KEYS)}")
        if raw["characterized_policy"] not in ("route_to_control", "keep_in_discovery"):
            raise PolicyError("characterized_policy must be route_to_control or keep_in_discovery")

        rules = tuple((r["stratum"], frozenset(r["terms"])) for r in raw["stratum_rules"])
        policy_strata = {s for s, _ in rules} | {raw["default_stratum"]}

        predictors = dict(raw["predictors"])
        if set(predictors) != policy_strata:
            raise PolicyError(
                f"predictors must name every stratum exactly: missing "
                f"{sorted(policy_strata - set(predictors))}, extra {sorted(set(predictors) - policy_strata)}"
            )
        used = {p for p in predictors.values() if p}
        for outer, inner in _CONTAINS.items():
            if outer in used and used & inner:
                raise PolicyError(f"{outer} contains {sorted(inner)}; they cannot both be keys")

        stratum_level: dict[Mechanism, dict[str, Concordance]] = {}
        subclass: dict[Mechanism, dict[str, SubclassRule]] = {}
        for mech in Mechanism:
            row = table["stratum_level"].get(mech.value)
            if row is None or set(row) != policy_strata:
                raise PolicyError(f"stratum_level for {mech.value} must cover strata {sorted(policy_strata)}")
            stratum_level[mech] = {s: Concordance(v) for s, v in row.items()}
            subclass[mech] = {
                s: SubclassRule(
                    feature=r["feature"],
                    values={str(k).lower(): Concordance(v) for k, v in r["values"].items()},
                    missing=Concordance(r["missing"]),
                )
                for s, r in (table["subclass"].get(mech.value) or {}).items()
            }
        if any(c is not Concordance.NOT_ASSESSABLE for c in stratum_level[Mechanism.UNKNOWN].values()) or subclass[
            Mechanism.UNKNOWN
        ]:
            raise PolicyError("unknown mechanism must be not_assessable everywhere")

        return cls(
            policy_id=raw["policy_id"],
            status=raw["status"],
            key_order=key_order,
            characterized_policy=raw["characterized_policy"],
            control_roles=frozenset(raw.get("control_roles", [])),
            stratum_rules=rules,
            default_stratum=raw["default_stratum"],
            predictors=predictors,
            table_id=table["table_id"],
            stratum_level=stratum_level,
            subclass=subclass,
            digest=content_digest({"policy": raw, "table": table}),
            raw={"policy": raw, "table": table},
        )

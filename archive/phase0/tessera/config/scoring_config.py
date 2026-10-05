"""Scoring configuration loader.

Everything checkable is checked here, at load, so that a misconfiguration is a
startup failure rather than a plausible-looking wrong number in a ranked table.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import yaml
from pydantic import Field, model_validator

from tessera.config.anchors import AnchorTable, MissingAnchorError
from tessera.enums import ConsequenceClass, Mechanism, RankingMode, ScoreComponent
from tessera.schemas.common import Strict

WEIGHT_SUM_TOLERANCE = 1e-9


class ScoringConfigError(ValueError):
    pass


class ModeWeights(Strict):
    """The four weighted components, plus the novelty multiplier.

    `novelty_alpha` is NOT a weight and is excluded from the sum-to-one check --
    that is the whole point of spec 19.5.
    """

    disease: float
    functional: float
    mechanism: float
    editability: float

    def as_dict(self) -> dict[ScoreComponent, float]:
        return {
            ScoreComponent.DISEASE: self.disease,
            ScoreComponent.FUNCTIONAL: self.functional,
            ScoreComponent.MECHANISM: self.mechanism,
            ScoreComponent.EDITABILITY: self.editability,
        }

    @model_validator(mode="after")
    def _sum_to_one(self) -> "ModeWeights":
        total = sum(self.as_dict().values())
        if not math.isclose(total, 1.0, abs_tol=1e-6):
            raise ScoringConfigError(
                f"component weights must sum to 1.0, got {total:.6f}. If you are "
                f"trying to add a fifth component, read spec 19.5 first: novelty "
                f"belongs in novelty_alpha, not in the weighted sum."
            )
        return self


class PanelComposition(Strict):
    discovery_fraction: float
    validation_fraction: float
    control_fraction: float

    @model_validator(mode="after")
    def _fractions_sum_to_one(self) -> "PanelComposition":
        total = self.discovery_fraction + self.validation_fraction + self.control_fraction
        if not math.isclose(total, 1.0, abs_tol=1e-6):
            raise ScoringConfigError(f"panel composition fractions must sum to 1.0, got {total}")
        return self


class ModeConfig(Strict):
    name: RankingMode
    weights: ModeWeights
    novelty_alpha: float = Field(ge=0.0, le=1.0)
    composition: PanelComposition | None = None


class DiseaseSubWeights(Strict):
    fine_mapping: float
    clinical_evidence: float
    recurrence: float
    rarity: float
    segregation: float

    @model_validator(mode="after")
    def _sum_to_one(self) -> "DiseaseSubWeights":
        total = (
            self.fine_mapping
            + self.clinical_evidence
            + self.recurrence
            + self.rarity
            + self.segregation
        )
        if not math.isclose(total, 1.0, abs_tol=1e-6):
            raise ScoringConfigError(f"disease sub-weights must sum to 1.0, got {total}")
        return self


class FunctionalConfig(Strict):
    aggregation: str = "max"
    replication_bonus: float = 0.0

    @model_validator(mode="after")
    def _max_only(self) -> "FunctionalConfig":
        if self.aggregation != "max":
            raise ScoringConfigError(
                f"functional aggregation must be 'max', got {self.aggregation!r}. "
                f"Summing evidence tiers reintroduces paper-counting after it was "
                f"excluded from recurrence (spec 19.2)."
            )
        return self


class EditabilitySubWeights(Strict):
    predicted_efficiency: float
    design_count: float
    bystander: float
    off_target: float
    indel_length: float

    @model_validator(mode="after")
    def _sum_to_one(self) -> "EditabilitySubWeights":
        total = sum(
            (
                self.predicted_efficiency,
                self.design_count,
                self.bystander,
                self.off_target,
                self.indel_length,
            )
        )
        if not math.isclose(total, 1.0, abs_tol=1e-6):
            raise ScoringConfigError(f"editability sub-weights must sum to 1.0, got {total}")
        return self


class MechanismBranches(Strict):
    """M per consequence class (spec 19.3). Each branch's weights sum to 1.0, and
    every raw input named must be anchored."""

    version: str
    branches: dict[ConsequenceClass, dict[str, float]]

    @model_validator(mode="after")
    def _each_branch_sums_to_one(self) -> "MechanismBranches":
        for cls_, weights in self.branches.items():
            total = sum(weights.values())
            if not math.isclose(total, 1.0, abs_tol=1e-6):
                raise ScoringConfigError(
                    f"mechanism branch {cls_.value!r} weights sum to {total}, not 1.0"
                )
        return self

    def raw_inputs(self) -> list[str]:
        seen: list[str] = []
        for weights in self.branches.values():
            for key in weights:
                if key not in seen:
                    seen.append(key)
        return seen

    def for_class(self, cls_: ConsequenceClass) -> dict[str, float]:
        if cls_ not in self.branches:
            raise ScoringConfigError(
                f"no M branch defined for consequence class {cls_.value!r}. Every class "
                f"needs an explicit branch -- falling back to another class's predictors "
                f"is the category error this split exists to prevent."
            )
        return self.branches[cls_]


class MechanismGating(Strict):
    version: str
    gating: dict[Mechanism, dict[str, float]]

    def multiplier(self, mechanism: Mechanism, variant_gate_key: str) -> tuple[float, str]:
        """Return (multiplier, key_used). Falls back to the mechanism's `default`
        only if declared; otherwise 1.0 with the reason recorded."""
        table = self.gating.get(mechanism, {})
        if variant_gate_key in table:
            return table[variant_gate_key], variant_gate_key
        if "default" in table:
            return table["default"], "default"
        return 1.0, "ungated"


class SaturationConfig(Strict):
    version: str
    levels: dict[str, float]
    missingness: dict[str, dict[str, Any]]

    def level(self, key: str) -> float:
        if key not in self.levels:
            raise ScoringConfigError(
                f"unknown saturation level {key!r}. Declared: {sorted(self.levels)}"
            )
        return self.levels[key]


class ScoringConfig(Strict):
    """Fully loaded, fully validated scoring configuration."""

    version: str
    config_dir: Path
    modes: dict[RankingMode, ModeConfig]
    disease: DiseaseSubWeights
    functional: FunctionalConfig
    editability: EditabilitySubWeights
    branches: MechanismBranches
    gating: MechanismGating
    saturation: SaturationConfig
    anchors: AnchorTable
    forbidden_scoring_fields: list[str]

    # ------------------------------------------------------------------ load

    @classmethod
    def load(cls, config_dir: str | Path, *, scoring_file: str = "scoring.yaml") -> "ScoringConfig":
        d = Path(config_dir)
        raw = yaml.safe_load((d / scoring_file).read_text())

        anchors = AnchorTable.load(d / raw.get("anchors", "anchors.yaml"))
        branches = MechanismBranches(
            **yaml.safe_load((d / raw["mechanism_branches"]).read_text())
        )
        gating = MechanismGating(**yaml.safe_load((d / raw["mechanism_gating"]).read_text()))
        saturation = SaturationConfig(**yaml.safe_load((d / raw["saturation"]).read_text()))

        modes: dict[RankingMode, ModeConfig] = {}
        mode_specs: dict[str, dict] = raw["ranking_modes"]
        for name, spec in mode_specs.items():
            spec = dict(spec)
            if "inherits" in spec:
                parent = mode_specs[spec.pop("inherits")]
                merged = {k: v for k, v in parent.items() if k != "composition"}
                merged.update(spec)
                spec = merged
            modes[RankingMode(name)] = ModeConfig(
                name=RankingMode(name),
                weights=ModeWeights(**spec["weights"]),
                novelty_alpha=spec["novelty_alpha"],
                composition=(
                    PanelComposition(**spec["composition"]) if "composition" in spec else None
                ),
            )

        comp = raw["components"]
        cfg = cls(
            version=str(raw["version"]),
            config_dir=d,
            modes=modes,
            disease=DiseaseSubWeights(**comp["disease"]),
            functional=FunctionalConfig(**comp["functional"]),
            editability=EditabilitySubWeights(**comp["editability"]),
            branches=branches,
            gating=gating,
            saturation=saturation,
            anchors=anchors,
            forbidden_scoring_fields=list(raw["validation"]["forbid_fields_in_scoring"]),
        )
        cfg._require_all_anchors()
        return cfg

    # ---------------------------------------------------------- consistency

    def _require_all_anchors(self) -> None:
        """Spec 19.6: every raw input referenced by any component must be anchored.

        A gap here is a startup failure, not a wrong number discovered later.
        """
        required = set(self.branches.raw_inputs())
        required |= {
            # disease
            "gwas_pip",
            "gwas_log10_p",
            "clinvar_classification",
            "clinvar_review_stars",
            "disease_match_type",
            "cohort_normalized_recurrence",
            "segregation_support",
            "gnomad_af",
            # functional
            "evidence_directness",
            "model_system_relevance",
            "replication_count",
            "mave_functional_class",
            # editability
            "pegrna_predicted_efficiency",
            "pegrna_design_count",
            "bystander_risk",
            "off_target_risk",
            "indel_length",
        }
        self.anchors.require(sorted(required))

    def mode(self, mode: RankingMode) -> ModeConfig:
        if mode not in self.modes:
            raise ScoringConfigError(
                f"ranking mode {mode.value!r} is not configured. Available: "
                f"{sorted(m.value for m in self.modes)}"
            )
        return self.modes[mode]

    def fingerprint(self) -> dict[str, str]:
        """Goes into the run manifest. Every version that can change a number."""
        return {
            "scoring_version": self.version,
            "anchors_version": self.anchors.version,
            "branches_version": self.branches.version,
            "gating_version": self.gating.version,
            "saturation_version": self.saturation.version,
        }

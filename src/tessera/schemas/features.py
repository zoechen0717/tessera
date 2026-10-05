"""Features and ranked rows for nomination_v0.1 (SPEC §8, SCHEMA §12).

K and O are ordinal categories; P is a raw predictor value compared only within
one stratum and one predictor. There is no global score field.
"""

from __future__ import annotations

from enum import StrEnum

from pydantic import Field

from tessera.schemas.base import SCHEMA_VERSION, Strict
from tessera.schemas.observed import Observed


class Concordance(StrEnum):
    CONCORDANT = "concordant"
    PARTIAL = "partial"
    NOT_ASSESSABLE = "not_assessable"
    DISCORDANT = "discordant"

    @property
    def order(self) -> int:
        """Higher is better. not_assessable sits above discordant: unknown is not
        evidence against, but is not rewarded either."""
        return {"concordant": 3, "partial": 2, "not_assessable": 1, "discordant": 0}[self.value]


class ObservationTier(StrEnum):
    CASE_CONTROLLED = "case_controlled"
    CASE_UNCONTROLLED = "case_uncontrolled"
    RELATED_CONDITION = "related_condition"
    NONE = "none"

    @property
    def order(self) -> int:
        return {"case_controlled": 3, "case_uncontrolled": 2, "related_condition": 1, "none": 0}[
            self.value
        ]


class Channel(StrEnum):
    DISCOVERY = "discovery"
    CONTROL = "control"


class Characterization(StrEnum):
    UNCHARACTERIZED = "uncharacterized"
    CHARACTERIZED_ALTERED = "characterized_altered"
    CHARACTERIZED_OTHER = "characterized_other"


class FeatureRow(Strict):
    schema_version: str = SCHEMA_VERSION
    variant_id: str
    label: str
    policy_id: str
    policy_digest: str
    policy_status: str
    evidence_digest: str
    primary_channel: Channel
    roles: list[str]
    stratum: str
    consequence_terms: list[str]
    characterization: Characterization
    k: Concordance
    k_basis: str
    o: ObservationTier
    o_search_status: str = Field(description="Whether observation sources were searched")
    p_predictor: str | None
    p: Observed[float] | None
    contributing_link_ids: dict[str, list[str]] = Field(default_factory=dict)
    flags: list[str] = Field(default_factory=list)


class RankedRow(Strict):
    schema_version: str = SCHEMA_VERSION
    variant_id: str
    label: str
    primary_channel: Channel
    stratum: str
    tier_rank: int = Field(ge=1, description="Dense rank; equal keys share a rank")
    tie_group_id: str
    tie_group_size: int = Field(ge=1)
    display_index: int = Field(ge=1, description="Coordinate order within tie group; no priority meaning")
    ordering_tuple: list[str]


class KeyDegeneracy(Strict):
    key: str
    distinct_values: int
    degenerate: bool


class StratumReport(Strict):
    channel: Channel
    stratum: str
    n: int
    n_tiers: int
    largest_tie_group: int
    class_level_concordance: Concordance
    keys: list[KeyDegeneracy]
    p_predictor: str | None
    p_missing: int

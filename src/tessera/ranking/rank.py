"""Stratified lexicographic ranking with explicit ties (SPEC §8.2–8.3, DEC-05).

No rank ever compares alleles from different (channel, stratum) groups
(INV-18). Rows equal on every key share a dense tier rank and a tie group; the
coordinate order inside a tie group is for display and carries no priority
(INV-19). Each key reports whether it was constant -- a degenerate key did no
work in that stratum and must not be described as having shaped the order.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import groupby

from tessera.ranking.policy import NominationPolicy
from tessera.schemas.annotation import GeneMechanismAssessment
from tessera.schemas.features import (
    Channel,
    FeatureRow,
    KeyDegeneracy,
    RankedRow,
    StratumReport,
)

_CHROM_SORT = lambda c: (0, int(c), "") if c.isdigit() else (1, 0, c)  # noqa: E731


def _key_value(row: FeatureRow, key: str) -> tuple:
    if key == "K":
        return (row.k.order,)
    if key == "O":
        return (row.o.order,)
    if key == "P":
        if row.p is not None and row.p.is_present:
            return (1, float(row.p.value))
        return (0, 0.0)
    raise KeyError(key)


def _key_label(row: FeatureRow, key: str) -> str:
    if key == "K":
        return f"K={row.k.value}"
    if key == "O":
        return f"O={row.o.value}"
    if row.p is None:
        return "P=no_predictor"
    if row.p.is_present:
        return f"P={row.p.value:g}"
    return f"P=missing:{row.p.status.value}"


def rank(
    rows: list[FeatureRow],
    coords: dict[str, tuple[str, int, str, str]],
    policy: NominationPolicy,
    mechanism: GeneMechanismAssessment,
) -> tuple[list[RankedRow], list[StratumReport]]:
    groups: dict[tuple[Channel, str], list[FeatureRow]] = defaultdict(list)
    for r in rows:
        groups[(r.primary_channel, r.stratum)].append(r)

    def coord_key(r: FeatureRow):
        chrom, pos, ref, alt = coords[r.variant_id]
        return (_CHROM_SORT(chrom), pos, ref, alt)

    def sort_key(r: FeatureRow):
        return tuple(_key_value(r, k) for k in policy.key_order)

    ranked: list[RankedRow] = []
    reports: list[StratumReport] = []
    stratum_order = {s: i for i, s in enumerate(policy.strata)}
    for (channel, stratum) in sorted(groups, key=lambda g: (g[0].value, stratum_order.get(g[1], 99))):
        members = groups[(channel, stratum)]
        ordered = sorted(members, key=lambda r: (tuple(-x for t in sort_key(r) for x in t), coord_key(r)))
        display = 0
        tiers = 0
        largest = 0
        for tier_rank, (_, tie) in enumerate(groupby(ordered, key=sort_key), start=1):
            tie = list(tie)
            tiers = tier_rank
            largest = max(largest, len(tie))
            for r in tie:
                display += 1
                ranked.append(
                    RankedRow(
                        variant_id=r.variant_id,
                        label=r.label,
                        primary_channel=channel,
                        stratum=stratum,
                        tier_rank=tier_rank,
                        tie_group_id=f"{channel.value}:{stratum}:t{tier_rank}",
                        tie_group_size=len(tie),
                        display_index=display,
                        ordering_tuple=[_key_label(r, k) for k in policy.key_order],
                    )
                )
        keys = []
        for k in policy.key_order:
            distinct = len({_key_value(r, k) for r in members})
            keys.append(KeyDegeneracy(key=k, distinct_values=distinct, degenerate=distinct <= 1))
        predictor = policy.predictors.get(stratum)
        reports.append(
            StratumReport(
                channel=channel,
                stratum=stratum,
                n=len(members),
                n_tiers=tiers,
                largest_tie_group=largest,
                class_level_concordance=policy.stratum_level[mechanism.mechanism][stratum],
                keys=keys,
                p_predictor=predictor,
                p_missing=sum(1 for r in members if predictor and not (r.p and r.p.is_present)),
            )
        )
    return ranked, reports

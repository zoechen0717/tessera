"""Deterministic CSV and Markdown rendering. No timestamps: replay compares bytes."""

from __future__ import annotations

import csv
import io
from collections import Counter

from tessera.ranking.policy import NominationPolicy
from tessera.schemas.annotation import GeneMechanismAssessment
from tessera.schemas.features import FeatureRow, RankedRow, StratumReport

CSV_COLUMNS = [
    "primary_channel", "stratum", "tier_rank", "tie_group_id", "tie_group_size", "display_index",
    "variant_id", "label", "K", "O", "P_predictor", "P", "ordering_tuple", "characterization",
    "roles", "consequence_terms", "k_basis", "o_search_status", "editing_status", "flags",
    "contributing_link_ids", "policy_id", "policy_status",
]


def render_csv(ranked: list[RankedRow], features: dict[str, FeatureRow]) -> str:
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(CSV_COLUMNS)
    for r in ranked:
        f = features[r.variant_id]
        p = "" if f.p is None else (f"{f.p.value:g}" if f.p.is_present else f.p.status.value)
        links = ";".join(f"{k}:{','.join(v)}" for k, v in sorted(f.contributing_link_ids.items()))
        w.writerow([
            r.primary_channel.value, r.stratum, r.tier_rank, r.tie_group_id, r.tie_group_size,
            r.display_index, r.variant_id, r.label, f.k.value, f.o.value, f.p_predictor or "",
            p, "|".join(r.ordering_tuple), f.characterization.value, ";".join(f.roles),
            ";".join(f.consequence_terms), f.k_basis, f.o_search_status, "not_assessed",
            ";".join(f.flags), links, f.policy_id, f.policy_status,
        ])
    return buf.getvalue()


def render_report(
    *,
    scope: dict,
    policy: NominationPolicy,
    mechanism: GeneMechanismAssessment,
    features: list[FeatureRow],
    ranked: list[RankedRow],
    strata: list[StratumReport],
    exclusions: list[dict],
    n_mentions: int,
    link_states: dict[str, str],
    obs_search: dict,
) -> str:
    gene = scope["gene"]
    L: list[str] = []
    add = L.append
    add(f"# tessera nomination report — {gene['symbol']} / {scope['disease'].get('input')}")
    add("")
    if scope["is_synthetic"]:
        add("> **SYNTHETIC FIXTURE RUN.** Reference, alleles, sources and evidence are synthetic. "
            "Nothing here is an observation about any real gene or variant.")
        add("")
    add(f"> Ranking profile `{policy.policy_id}` is **{policy.status}**: a proposed policy, not a "
        "validated predictor of effect size, pathogenicity or experimental success. "
        "This table is a nomination for review, not an optimized or ready-to-order panel.")
    add("")

    add("## Scope and coverage")
    add("")
    add(f"- Gene: `{gene['gene_id']}` ({gene['symbol']}), assembly `{scope['assembly']}`, "
        f"reference bundle `{scope['reference_bundle_id']}`")
    add(f"- Candidate route: manual input only (M0). Mentions: {n_mentions}; "
        f"resolved alleles: {len(features)}; excluded mentions: {len(exclusions)}")
    for status, n in sorted(Counter(e["status"] for e in exclusions).items()):
        add(f"  - excluded `{status}`: {n}")
    add(f"- Observation sources: `{obs_search['status']}`"
        + (f" ({obs_search.get('scope')})" if obs_search.get("scope") else ""))
    states = Counter(link_states.values())
    add(f"- Evidence links: {len(link_states)} — "
        + ", ".join(f"{k} {v}" for k, v in sorted(states.items())))
    add("- Not assessed in this run: editing feasibility (all `not_assessed`), GEO datasets, "
        "literature discovery (M1b).")
    if scope.get("requested_panel_size") is not None:
        add(f"- Requested panel size {scope['requested_panel_size']}: **unsupported** — panel step "
            "not implemented.")
    add("")

    add("## Gene mechanism (K)")
    add("")
    add(f"- Mechanism: **{mechanism.mechanism.value}** (status `{mechanism.status}`, decision "
        f"`{mechanism.decision}`, policy `{mechanism.policy_id}`)")
    for b in mechanism.basis:
        add(f"  - {b.source} `{b.record_id}` {b.field} = {b.value}")
    add("- K is a class-level inference from this mechanism and the concordance table "
        f"`{policy.table_id}`. It is never evidence about an individual allele.")
    add("")

    add("## Strata")
    add("")
    add("Ranks compare alleles only within one channel and consequence stratum. "
        "Equal keys share a tier; order inside a tier is genomic coordinate and means nothing.")
    add("")
    add("| channel | stratum | n | tiers | largest tie | class-level K | degenerate keys | P predictor | P missing |")
    add("|---|---|---|---|---|---|---|---|---|")
    for s in strata:
        deg = ", ".join(k.key for k in s.keys if k.degenerate) or "—"
        add(f"| {s.channel.value} | {s.stratum} | {s.n} | {s.n_tiers} | {s.largest_tie_group} | "
            f"{s.class_level_concordance.value} | {deg} | {s.p_predictor or '—'} | {s.p_missing} |")
    add("")
    warnings = []
    for s in strata:
        for k in s.keys:
            if k.key == "P" and s.p_predictor is None:
                continue
            if k.degenerate and s.n > 1:
                warnings.append(f"`{s.channel.value}/{s.stratum}`: key {k.key} is constant "
                                f"across {s.n} alleles and did not influence the order.")
        if s.n > 1 and s.n_tiers == 1:
            warnings.append(f"`{s.channel.value}/{s.stratum}`: all {s.n} alleles are tied; "
                            "this stratum has no ranking.")
    no_pred = [f"{s.channel.value}/{s.stratum}" for s in strata if s.p_predictor is None]
    if no_pred:
        warnings.append("No P predictor is declared for " + ", ".join(f"`{x}`" for x in no_pred)
                        + " (policy choice, SCI-10).")
    if warnings:
        add("**Degeneracy warnings**")
        add("")
        L.extend(f"- {w}" for w in warnings)
        add("")

    feat = {f.variant_id: f for f in features}
    add("## Ranked alleles")
    add("")
    current = None
    for r in ranked:
        key = (r.primary_channel.value, r.stratum)
        if key != current:
            if current is not None:
                add("")
            current = key
            add(f"### {key[0]} / {key[1]}")
            add("")
            add("| tier | allele | K | O | P | characterization | flags |")
            add("|---|---|---|---|---|---|---|")
        f = feat[r.variant_id]
        p = "—" if f.p is None else (f"{f.p.value:g}" if f.p.is_present else f"missing ({f.p.status.value})")
        tier = f"{r.tier_rank}" + (f" (tie ×{r.tie_group_size})" if r.tie_group_size > 1 else "")
        add(f"| {tier} | `{r.label}` | {f.k.value} | {f.o.value} | {p} | {f.characterization.value} | "
            f"{', '.join(f.flags) or '—'} |")
    add("")

    add("## Excluded mentions")
    add("")
    if exclusions:
        add("| mention | status | reason |")
        add("|---|---|---|")
        for e in exclusions:
            add(f"| {e['mention_id']} | {e['status']} | {e['reason']} |")
    else:
        add("None.")
    add("")
    add("## Limitations")
    add("")
    add("- O=`none` with observation search status `not_searched` means unknown, not absent.")
    add("- Missing predictions sort after all scored alleles in their stratum; they are not zero.")
    add("- Characterized alleles are routed to the control channel (DEC-02); a no-effect result is "
        "not evidence of benignity.")
    add("- Cross-stratum allocation belongs to the panel step, which does not exist yet.")
    return "\n".join(L) + "\n"

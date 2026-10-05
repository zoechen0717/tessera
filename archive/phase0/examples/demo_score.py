"""Phase 0 demo: score a small hand-built candidate set, with no network and no LLM.

Deliberately built to reproduce the SETD1A + schizophrenia situation described in
spec 47.0 -- gene-level disease evidence that is identical across candidates -- so
that the degenerate-component diagnostic fires and the report says which
components actually drove the ordering.

Run:  ./.venv/bin/python examples/demo_score.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "tests")]

from conftest import make_annotations, make_evidence, make_variant  # noqa: E402

from tessera.config.scoring_config import ScoringConfig  # noqa: E402
from tessera.enums import (  # noqa: E402
    ConsequenceClass,
    DiseaseMatchType,
    EvidenceLevel,
    Mechanism,
    MechanismConfidence,
    RankingMode,
)
from tessera.ranking.ranker import diagnose, rank_stratified, score_variant  # noqa: E402
from tessera.schemas.gene import GeneMechanism  # noqa: E402

CANDIDATES = [
    # label,            class,                       alphamissense, domain,           evidence
    ("p.Gln191His",     ConsequenceClass.MISSENSE,   0.93, "catalytic_domain",        "direct"),
    ("p.Arg1392Cys",    ConsequenceClass.MISSENSE,   0.71, "annotated_domain",        None),
    ("p.Ser512Leu",     ConsequenceClass.MISSENSE,   0.22, "disordered_region",       None),
    ("p.Tyr1499*",      ConsequenceClass.PTV,        None, "catalytic_domain",        None),
    ("c.4582+1G>A",     ConsequenceClass.SPLICE,     None, "annotated_domain",        None),
    ("p.Leu88=",        ConsequenceClass.SYNONYMOUS, None, "outside_annotated_feature", None),
]


def build(*, variant_level_evidence: bool) -> list[tuple]:
    out = []
    for i, (label, cls_, am, domain, ev_kind) in enumerate(CANDIDATES):
        vid = f"ga4gh:VA.demo{i}"
        variant = make_variant(vrs_id=vid, consequence_class=cls_, hgvs_p=label)
        ann = make_annotations(
            vrs_id=vid,
            alphamissense=am,
            domain=domain,
            spliceai=0.92 if cls_ is ConsequenceClass.SPLICE else 0.03,
            af=2e-6,
            # Gene-level only: every candidate carries the SAME disease evidence,
            # which is what makes D degenerate here.
            clinvar_class=None,
            stars=None,
            efficiency=0.70 - 0.06 * i,
        )
        # Every candidate inherits the SAME gene-level disease evidence, whose
        # condition is BROADER than the query -- the SETD1A + schizophrenia
        # situation. D is therefore constant and contributes nothing to the order.
        evidence = [
            make_evidence(
                vrs_id=None,
                evidence_id=f"ev_gene_{i}",
                level=EvidenceLevel.GENE_LEVEL,
                directness="gene_level_relevant_cell",
                model_system="related_human_neuronal",
                match=DiseaseMatchType.BROADER,
                proband=None,
                source_id="PMID:burden",
            )
        ]
        if ev_kind == "direct" and variant_level_evidence:
            evidence.append(
                make_evidence(
                    vrs_id=vid,
                    evidence_id=f"ev_{i}",
                    match=DiseaseMatchType.BROADER,
                    model_system="related_human_neuronal",
                )
            )
        out.append((label, variant, ann, evidence))
    return out


def run_scenario(cfg, mode, mechanism, title: str, *, variant_level_evidence: bool) -> None:
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")

    rows, scores = [], []
    for label, variant, ann, evidence in build(variant_level_evidence=variant_level_evidence):
        s = score_variant(cfg, mode, variant, ann, evidence, mechanism)
        scores.append(s)
        rows.append((label, s))

    hdr = (
        f"{'variant':<16}{'class':<12}{'D':>7}{'F':>7}{'M':>7}{'E':>7}"
        f"{'sat':>7}{'mult':>7}{'final':>8}"
    )
    print(hdr)
    print("-" * len(hdr))
    for label, s in rows:

        def fmt(c) -> str:
            return "    n/a" if c.value is None else f"{c.value:7.3f}"

        sat = "    n/a" if s.saturation.value is None else f"{s.saturation.value:7.2f}"
        print(
            f"{label:<16}{s.consequence_class.value:<12}"
            f"{fmt(s.disease)}{fmt(s.functional)}{fmt(s.mechanism)}{fmt(s.editability)}"
            f"{sat}{s.novelty_multiplier:7.3f}{s.final_score:8.4f}"
        )

    print("\nRanked within consequence class (cross-class pooling is not meaningful):")
    for cls_, group in rank_stratified(scores).items():
        labels = [next(l for l, s in rows if s.vrs_id == g.vrs_id) for g in group]
        print(f"  {cls_.value:<12} " + " > ".join(labels))

    diag = diagnose(cfg, mode, scores)
    print(f"\nComponent discrimination over {diag.n_candidates} candidates:")
    print(
        f"  {'component':<14}{'weight':>8}{'variance':>11}{'eff.contrib':>13}"
        f"{'coverage':>10}  verdict"
    )
    for row in diag.discrimination:
        verdict = "DEGENERATE" if row.degenerate else ("sparse" if row.sparse else "-")
        print(
            f"  {row.component.value:<14}{row.configured_weight:8.2f}{row.variance:11.5f}"
            f"{row.effective_contribution:13.4f}{row.coverage:10.2f}  {verdict}"
        )
    for w in diag.warnings:
        print(f"\n  WARNING: {w}")


def main() -> None:
    cfg = ScoringConfig.load(ROOT / "config")
    mechanism = GeneMechanism(
        gene_symbol="SETD1A",
        mechanism=Mechanism.HAPLOINSUFFICIENCY,
        confidence=MechanismConfidence.HIGH,
    )
    mode = RankingMode.DISCOVERY

    print(f"scoring fingerprint: {cfg.fingerprint()}")
    print(f"mode: {mode.value}  alpha={cfg.mode(mode).novelty_alpha}")
    print(f"gene mechanism: {mechanism.mechanism.value} ({mechanism.confidence.value})")

    run_scenario(
        cfg,
        mode,
        mechanism,
        "Scenario A -- one variant has direct functional evidence",
        variant_level_evidence=True,
    )
    print(
        "\n  Note: p.Gln191His scores LOWER on D than candidates with no recurrence\n"
        "  data at all. That is intended. A single proband in a 4,000-person cohort is\n"
        "  near chance expectation, so the recurrence term enters at a low value --\n"
        "  whereas for the others it is absent and its weight is renormalized away.\n"
        "  Weak evidence is not the same as no evidence, and the scorer must not make\n"
        "  a variant look better by knowing less about it."
    )

    run_scenario(
        cfg,
        mode,
        mechanism,
        "Scenario B -- gene-level evidence only (the SETD1A + schizophrenia case)",
        variant_level_evidence=False,
    )
    print(
        "\n  This is spec 47.0. Every candidate inherits the same gene-level burden\n"
        "  evidence for a BROADER condition, so D is identical across the set and its\n"
        "  configured weight of 0.32 does no work. The ordering is driven by M and E.\n"
        "  The run says so rather than presenting a disease-weighted score."
    )


if __name__ == "__main__":
    main()

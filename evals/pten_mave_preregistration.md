# Pre-registration: PTEN / autism — nomination_v0.1 vs MAVE readouts

Written 2026-10-06, **before** running tessera on PTEN or looking at any overlap with MAVE data.
Changes after this point are logged as amendments with a reason; the analysis below is not tuned
on the outcome.

## Question

Within the discovery missense stratum, does the `nomination_v0.1` order put functionally
abnormal alleles ahead of normal ones, and does it beat simpler baselines?

## Inputs (tessera run)

- Request: `requests/pten_autism.yaml` — gene PTEN, disease autism spectrum disorder
  (MONDO:0005258). Routes: ClinVar, LitVar2. **No literature deep read** and no screening of the
  MAVE papers, so no MAVE result can enter the ranking as functional evidence. No case/control
  browser (O is expected to be uninformative; reported as degenerate if so).
- Policy and concordance table unchanged from SETD1A (`nomination_v0.1`,
  `mechanism_concordance_v0.1`).

## Ground truth (never a ranking input)

| Role | MaveDB score set | Readout | Label |
|---|---|---|---|
| Primary | urn:mavedb:00000054-a-1 (PMID 29706350) | lipid phosphatase activity | **abnormal** = score ≤ −2.13 ("Truncation-like", IGVF Coding Variant Focus Group calibration); **normal** = score > −1.11; between = intermediate (excluded from binary metrics, kept for rank correlation) |
| Secondary | urn:mavedb:00000102-0-1 (PMID 34649609) | VAMP-seq abundance (combined) | continuous only (no calibration published in MaveDB) |
| Secondary | urn:mavedb:00001244-a-1 | morphological impact in hiPSC-derived neurons (VIS-seq) | continuous only |

Joining: MaveDB `hgvs_pro` single amino-acid substitutions ↔ tessera MANE `hgvs_p`. Only
alleles that tessera ranked (discovery, missense) and that have a MAVE score are evaluated.
Synonymous/nonsense MAVE entries are not used.

## Metrics

1. Spearman ρ between tessera tier rank and MAVE score (direction: better tier ↔ lower
   function), for each score set.
2. Primary binary task (abnormal vs normal, intermediate excluded): AUROC of the tessera order;
   precision@k and recall@k at k = 10, 25, 50; tie groups handled by expected value over random
   order within the tie.
3. 95% CIs by bootstrap over alleles (2,000 resamples, seed 42).

## Baselines (same allele set)

- Random order (expected values).
- AlphaMissense alone (the stratum's P key without K).
- Domain membership alone (K without P).
- ClinVar lookup: P/LP first, then VUS, then B/LB — a lookup baseline, not independent validation.

## Known limitations declared in advance

- AlphaMissense was benchmarked (not trained) on DMS data including PTEN; any AlphaMissense-driven
  result is partly a re-statement of that benchmark.
- The candidate set is ascertained (ClinVar/LitVar alleles), not all possible missense; ClinVar
  submissions may themselves cite the MAVE data.
- One gene: a case study, not a generalization claim (EVAL §8.1a).
- O will likely be degenerate (no case/control source wired for autism in this run).

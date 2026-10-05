# Variant Evidence System — REVIEW of 2.0.0-draft

Date: 2026-10-05 | Reviewer: Claude (with project scientist) | Reviewed: SPEC, EVIDENCE_SCHEMA, ARCHITECTURE, EVAL, IMPLEMENTATION_PLAN, REVIEW_GUIDE

## 1. Decision

**Ready with concrete fixes.** The engineering contracts (allele identity, claim/link separation, observation status, replay modes, evaluation denominators) are sound and are kept. The scientific objective is not: `evidence_priority_v0.1` answers "which alleles already have the most evidence", while the owner's decision is "which *new* alleles to nominate for the next screen". The fixes are recorded in [DECISIONS.md](DECISIONS.md) and applied in version 2.1.0-draft. ISSUE-001–004 may proceed under the revised plan.

## 2. Findings

| ID | Sev | File / section | Violated requirement | Concrete failure | Fix (decision) | Affected tests |
|---|---|---|---|---|---|---|
| S1 | P1 | SPEC §8.2 | §1 objective (owner: nomination) | F ranks any exact-allele experiment above none, regardless of outcome. Alleles already measured fill the top of the discovery channel; the screen re-tests known alleles | Replace default profile with `nomination_v0.1`; already-characterized alleles route to the control channel (DEC-01, DEC-02) | ranking policy tests |
| S2 | P1 | SPEC §8.2 | R-12 (meaningful ordering) | For a gene whose disease evidence is gene-level burden (e.g. ultra-rare PTV burden), nearly every allele is D=F=C=0; order is decided by a SHA-256 tie-break, which also reshuffles on any `identity_version` bump | Dense tier ranks, explicit tie groups, coordinate display order, degenerate-key report (DEC-05) | tie/degeneracy tests |
| S3 | P1 | SPEC §8, SCHEMA §2.3 | R-09 left without a lawful outlet | Gene-mechanism knowledge (LoF vs GoF vs dominant-negative) cannot affect rank at all, so a synonymous allele with one cohort observation outranks a PTV in a haploinsufficient gene; for a GoF gene the PTV would be wrongly favoured if anything else were added naively | Class-level mechanism-concordance key K, labelled inferred/class-level, never exact-allele evidence; mechanism resolved per gene with `unknown` (DEC-03) | INV-04, INV-17 |
| S4 | P1 | SPEC §8.2 vs EVAL §8.2 | Cross-document consistency | SPEC pools all consequence classes in one tuple; EVAL requires consequence strata | Rank within consequence stratum; cross-class allocation is the panel step (DEC-04) | INV-18 |
| S5 | P2 | SPEC §8.2 D=2 | R-07 | A single carrier in a case cohort without a control denominator is near base rate for ultra-rare alleles, yet outranks an expert-panel classification | Observation key O distinguishes case-with-control-denominator from uncontrolled observation (DEC-07) | O-tier tests |
| E1 | P1 | SPEC §3, PLAN | Prior plan: "deterministic real-data run is the milestone worth protecting" | M1 bundles structured sources with the literature/GEO agents; the first real-gene table arrives only after all LLM work, so S2-type defects surface late | Split M1 into M1a (structured, no LLM) and M1b (literature/GEO workers) (DEC-10) | — |
| E2 | P1 | PLAN §0, ARCH §13 | — | Plan says "no repository built"; tessera Phase 0 exists (weighted composite, 32 tests). Builder would either rewrite or run two codebases | Evolve the `tessera` package; archive the weighted composite (DEC-09) | — |
| E3 | P1 | PLAN ISSUE-002 | M0 "no network" | Reference-aware HGVS libraries need real transcript resources (UTA/SeqRepo-style); they cannot run against a synthetic reference | M0 accepts genomic (VCF-style) alleles only; HGVS c./p. mapping moves to M1a (DEC-11) | identity fixtures |
| E4 | P2 | SPEC §8.2 tie-break | R-16 stability | `variant_id` includes `identity_version`; a normalization version bump changes every tied row's order with no biological cause | Covered by DEC-05 | replay test |
| E5 | P2 | ARCH §4 | R-08 | AlphaMissense (local tables already available) and splice predictors are absent from the source table although PredictionRecord names them; splice prediction matters for LoF genes | Add AlphaMissense, SpliceAI-class predictor, ClinGen dosage, gnomAD constraint, G2P to source plan (DEC-12) | — |
| E6 | P2 | SCHEMA §13, PLAN | Scope | 16 invariants + append-only revisions + SQLite checkpoints in M0 for a single-researcher project delays any usable output | M0 store is a content-addressed file bundle; SQLite/checkpoints arrive with M1b where retries and resumable LLM stages exist (DEC-13) | — |

No P0 finding: nothing in the package creates false identity or evidence as written.

## 3. Scientific decisions for the owner

Resolved in this review (see DECISIONS.md): objective (DEC-01), handling of characterized alleles (DEC-02), mechanism axis (DEC-03), stratification (DEC-04).

Still open, with build defaults:

| ID | Question | Default in 2.1 |
|---|---|---|
| SCI-09 | Key order inside a stratum: mechanism → observation → prediction? | `(K, O, P)`; configurable, provisional |
| SCI-10 | Which predictor per consequence class? | missense: AlphaMissense; splice/synonymous: SpliceAI-class Δ; PTV: none in v0.1 (NMD status is a K sub-class feature, DEC-14). One predictor per class; AVI is an alternative, never added on top of AlphaMissense |
| SCI-11 | Should common alleles be excluded from discovery? | Exported and flagged only; no filter until owner sets a threshold |
| SCI-12 | Panel quotas per class and controls | Not implemented in M0; panel step is its own issue |
| SCI-13 | Is the SETD1A scPRIME batch usable as ground truth? | Case study only (historically selected panel, low cells/allele, no biological replicates); cross-gene MAVE/DMS benchmark is primary |

## 4. Source / API feasibility audit

| Source | Checked | Finding |
|---|---|---|
| AlphaGenome Atlas / AVI | 2026-10-05, DeepMind blog (2026-09-08) and secondary coverage | Exists. All possible SNVs plus ~100M observed indels; web portal + API; **non-commercial academic license**; AVI integrates AlphaGenome + AlphaMissense + conservation, with SHAP-style decomposition. Consequence: AVI and AlphaMissense are not independent (DEC-12). Endpoint paths and release definitions not yet verified against API docs — builder must do this in M1a |
| Other adapters | not checked | ClinVar, gnomAD, VEP, UniProt/InterPro, ClinGen, G2P, MaveDB, Europe PMC, GEO must be verified at M1a/M1b implementation time |

## 5. Schema/policy changes and migration

- FeatureRow D/F/C → K/O/P plus `characterized` status (SCHEMA §12). No data exists yet; no migration.
- New object `GeneMechanismAssessment` (SCHEMA §2.5).
- New invariants INV-17 (mechanism concordance never creates exact-allele evidence), INV-18 (no cross-stratum rank), INV-19 (ties reported, not hidden).
- tessera Phase 0 weighted-composite modules archived; `Quantity` superseded by `Observed[T]`.

## 6. First implementation tasks

ISSUE-001–004 as revised in IMPLEMENTATION_PLAN 2.1: contracts → genomic identity → file bundle store → offline nomination slice with exports and replay. Cuts: SQLite, HGVS mapping, LLM workers, panel optimisation.

## 7. Remaining uncertainty

- Key order (K, O, P) is a judgement call; it should be compared against alternatives on the DMS benchmark before any utility claim.
- Mechanism resolution from structured sources will be `unknown` for many newer risk genes; then K is constant and the run must say so.
- DMS readouts measure molecular/cellular function in specific assays, not neuronal transcriptomic effects; agreement there is necessary, not sufficient.

Source: [AlphaGenome Atlas announcement](https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/)

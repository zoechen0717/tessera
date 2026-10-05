# Variant Evidence and Screen Prioritization System — SPEC

Version: 2.1.0-draft | Date: 2026-10-05 | Owner: project scientist

Status: reviewable engineering specification; proposed scientific policies are not validated biological conclusions. This document defines WHAT must be built. Read ARCHITECTURE.md for execution, EVIDENCE_SCHEMA.md for contracts, and EVAL.md for verification.

## 1. Scientific objective

Given a disease, gene, and experimental context, nominate new alleles for a variant-level functional screen as an auditable, one-genomic-allele-per-row table (DEC-01). The tool must transfer across genes and diseases; nothing gene-specific is hard-coded. The immediate use case is SETD1A, schizophrenia, human iPSC-derived neurons, and scPRIME. This is a starting benchmark, not evidence of generalization.

Two parallel project tracks exist:

| Track | Input → output | Version 2 scope |
|---|---|---|
| A | Disease → risk genes | Future; share evidence infrastructure but do not implement gene ranking now |
| B | Disease + gene → screen candidate variants | Primary build target |

The system distinguishes three questions:

1. What evidence connects this exact allele to the disease?
2. What molecular or cellular effects have been measured or predicted?
3. Can this allele be installed and measured in the proposed experiment?

Disease relevance, expected transcriptomic perturbation, and editing feasibility are separate outputs. No claim that one score is a probability of pathogenicity or experimental success is permitted without a separately validated calibration model.

## 2. Contract language and ownership

MUST = required behavior; SHOULD = default with a documented reason for deviation; MAY = optional. Requirement IDs are stable across documents.

| Decision | Owner | Implementer may do |
|---|---|---|
| Target scientific objective and evidence interpretation | Scientist | Propose alternatives and show consequences |
| Disease expansion, transcript policy, ranking profile, control definition | Scientist | Implement declared policy and document ambiguity |
| Tool wrappers, database, retries, exports, orchestration | Engineer/coding agent | Choose routine implementation details within contracts |
| Gold labels and biological acceptance | Scientist/reviewers | Build tooling; cannot invent truth labels |
| Provider/library choice | Engineer | Verify current official APIs and document capabilities |

Numerical examples and the baseline ordering below are **proposed policies**, not scientific facts. Implementation can proceed with them tagged `provisional`. A research claim of ranking utility requires independent evaluation. No additional user approval is needed to build or test the provisional implementation.

## 3. Scope by milestone

### M0: offline vertical slice

Manual genomic (VCF-style) allele input → reference validation (DEC-11) → fixture annotations → typed evidence → validation ledger → deterministic ranking → CSV/JSONL/report. No network or model credentials required. Synthetic examples must be prominently labeled; never masquerade as SETD1A observations.

### M1a: structured real-gene run (no LLM)

Disease + gene and optional manual variants → ClinVar gene-level candidates → canonical registry with HGVS c./p. mapping on real GRCh38/MANE → ClinVar, gnomAD (frequency and constraint), consequence, protein-domain, AlphaMissense and splice-prediction annotations → gene mechanism assessment (ClinGen dosage, gnomAD constraint, Gene2Phenotype) → `nomination_v0.1` → exports. First real SETD1A table. No model credentials required.

### M1b: research workers

Bounded PubMed/Europe PMC and GEO research, claim extraction and validation, one discovery expansion round, SQLite store and stage checkpoints, AlphaGenome Atlas adapter (license: non-commercial academic), editing adapter interface and design imports.

The following are required capabilities across M1a/M1b, but successful live retrieval is conditional on source access:

- Gene/disease resolution and exact allele normalization.
- ClinVar adapter with condition and submission-level provenance.
- gnomAD population annotation with release and callable/coverage information when available.
- Consequence and transcript mapping using pinned resources.
- UniProt/InterPro protein context with isoform mapping.
- PubMed/Europe PMC literature and GEO accession research.
- AlphaGenome Atlas annotation adapter. If unavailable, use an explicit status, not fabricated values.
- Editing adapter interface and imports of externally computed designs. M1 may return `not_assessed` and must label the output as nomination rather than a ready-to-order panel.

Optional extensions: GWAS/fine-mapping, Open Targets, regulatory variants outside coding/splice scope, broad web fallback, structural interpretation, editing design engines, panel optimization, cross-gene comparison. Missing optional adapters must not silently change the candidate scope.

### Out of scope

Clinical diagnosis, automatic ACMG classification, unrestricted computer/code execution by research workers, downloading and reprocessing all GEO/SRA reads, deploying a web application, training a predictor, and automatic experimental execution.

## 4. Input contract

Required input: disease text/ID, gene text/ID, species, assembly, cell context, assay, and objective profile. Defaults must be explicit in the resolved request. Human GRCh38, coding/splice SNVs and small indels are the proposed initial scope; exact indel-size bounds belong in the versioned capability config.

```yaml
spec_version: 2.0.0-draft
gene: {input: SETD1A}
disease:
  input: schizophrenia
  expansion_policy: exact_and_approved_synonyms
screen_context:
  species: Homo sapiens
  assembly: GRCh38
  cell_type: iPSC-derived neuron
  assay: scPRIME
  editing_mode: prime_editing
  genotype_target: heterozygous
  host_background: null
candidate_scope:
  consequences: [coding, splice]
  include_controls: true
  manual_variants_file: null
  locus_intervals: []
objective:
  ranking_profile: evidence_priority_v0.1
  policy_status: provisional
  requested_panel_size: null
execution:
  mode: refresh
  allow_partial_optional_sources: true
```

Cell type is insufficient to define prime editing: host genotype, local haplotype, editor/PAM configuration, design settings, and target zygosity may materially alter feasibility. Missing values are preserved.

R-01: resolve aliases from versioned structured resources. The LLM can propose ontology matches but cannot commit an ambiguous identity. Related developmental disorders are not automatically equivalent to schizophrenia. Record relationships: `exact`, `approved_synonym`, `related`, `unrelated`, `unknown`.

R-02: preserve all original variant descriptions, publication context, assembly, transcript version, and mapping alternatives. Protein shorthand, rsIDs, and truncated frameshift notation are not unique genomic identifiers. Unresolved identity is a valid output, but cannot enter the canonical ranking.

## 5. Candidate universe and completeness

R-03: define the universe BEFORE reporting recall or ranking. Required candidate routes:

| Route | Behavior | Important limitation |
|---|---|---|
| Manual | Import user alleles, retain provenance and role | Inclusion is not evidence of disease relevance |
| ClinVar gene-level retrieval | Enumerate within declared scope; classify disease match afterward | Filtering only disease-matched ClinVar records can miss useful candidates |
| Gene–disease literature | Discover allele mentions and exact variant studies | Full text/supplements may be inaccessible |
| Controls | Explicit user candidates or approved policy | Never infer benignity from rarity, synonymous consequence, or negative search |
| Prediction-driven enumeration, optional | Enumerate only predeclared intervals and allele types | No genome-wide unconstrained generation |
| Disease genetics, optional | Import exact variants and gene links | Burden evidence is gene/mutation-class evidence; proximity is not causality |

For every route, report query, filters, pagination/truncation, scope, retrieval status, and counts. Record candidates discovered, identity resolved, excluded, annotated, and deeply researched. Retain every exclusion reason.

Literature can reveal candidates after the initial registry is built. M1 permits one bounded discovery expansion round: normalize new mentions, add resolvable alleles, annotate, then research within remaining budget. Mentions arriving after the configured closure are stored in `candidate_mentions.jsonl` with `deferred_next_run`; they are not silently lost or inserted as unannotated ranked rows.

No publication found means `searched_not_found` under a documented search scope. It does not establish absence of evidence or novelty.

## 6. Agent versus software boundary

R-04: all model actions are bounded by host-owned state and tools.

| Component | Model may decide | Host/software must enforce |
|---|---|---|
| Literature scout | Which permitted query/tool addresses a gap | Templates, aliases, domains, counters, time limits |
| Extractor | Interpret supplied text into claim proposals | Schema, source locator, identity links, validation gates |
| Dataset scout | Interpret metadata relevance; request permitted details | Accession existence, sample/file metadata, retrieval logging |
| Semantic reviewer | Flag mismatch, unsupported inference, conflict | Cannot overrule invalid reference alleles or invent citations |
| Mechanistic interpreter | Explain evidence and explicit hypotheses | Separate observations, predictions, and inference; no rank edits |
| Normalization/annotation/ranking/export | No authority to alter calculations | Deterministic, versioned code |

An agent is optional for choosing searches; an LLM call that extracts text is a bounded semantic worker even if it has no autonomous tool loop. Multiple roles do not require multiple running agents or a multi-agent framework.

## 7. Evidence acceptance

R-05: accepted ranking-relevant evidence must reference a real source snapshot and an inspectable field or text location. Bibliographic existence alone is insufficient. A PMID is preferred for PubMed records, but DOI, PMCID, GEO accession, and versioned database record IDs are legitimate identifiers. “No PMID = no evidence” is explicitly rejected.

R-06: identity/directness, disease match, experiment quality, and model relevance are independent dimensions. Exact variant in a biochemical assay is direct evidence; a gene knockout in neurons is gene-level evidence. Their quality depends on the question, not a universal assay hierarchy.

R-07: database classification and patient recurrence are distinct. SCV counts, independent submitting organizations, publication counts, and deduplicated individuals/cohorts must never share one `recurrence` field. ClinVar classifications for a related condition are not disease-specific proof.

R-08: predictions remain predictions. Store raw and transformed values, source model/release, supported allele types, target feature/tissue, calibration/reference set, and missingness. AlphaGenome molecular tracks, Atlas AVI, and protein-impact prediction channels must be distinguished. Do not assume every AVI component or scoring method is exposed by every endpoint; inspect the release definition. No arbitrary rescaling of unrelated raw scores to a shared probability.

R-09: same-residue, same-domain, mutation-class, gene, and pathway evidence cannot be promoted to exact-allele evidence. Gene-level burden results do not assign association significance to every allele. Inferred LoF is not experimentally demonstrated LoF; terminal-exon/isoform/NMD limitations must be retained when known.

R-10: retain contradictory outcomes, no-detected-effect observations, retractions/corrections, inaccessible sources, and possible cohort overlap. A no-effect assay is not universal benign evidence; record its sensitivity and tested endpoint.

R-11: deterministic checks cannot prove semantic entailment. Export automatic and human validation separately. A second LLM provides another check, not ground truth. All evidence influencing the top review set must be inspectable by the scientist.

## 8. Ranking policy: nomination baseline

R-12: ranking is a pure function of frozen accepted evidence, annotations, explicit policy, and scope. LLM prose/confidence does not directly set the rank. The extracted evidence still contains uncertainty; determinism does not make it true.

The default profile is `nomination_v0.1` (DEC-01): which **new** alleles are most worth installing in a variant-level screen. It is provisional and lexicographic. The 2.0 profile `evidence_priority_v0.1` survives only as an evaluation baseline for known-evidence recovery.

### 8.1 Eligibility and channels

- Normalized, reference-validated, in-scope candidates are eligible.
- Unresolved/out-of-scope records remain in exclusions or mentions, with no rank.
- Channels: `discovery` and `control`. Explicit user controls are `control`.
- An allele with an accepted exact-allele functional experiment is `characterized` (DEC-02). Default primary channel `control`, role `positive_control_candidate` when the outcome is `altered`, otherwise `characterized_other`. `characterized_policy: keep_in_discovery` overrides this per run and is reported.
- Conflicting evidence does not remove a candidate; it adds a visible review flag.
- One ranked row per allele; roles may be a list; the primary channel is documented.

### 8.2 Strata

Ranking happens within `(channel, consequence_stratum)` (DEC-04). Strata in M0: `ptv`, `splice`, `missense`, `inframe_indel`, `synonymous`, `other`. A cross-stratum order is never produced by the ranker; allocation across strata is the panel step (§8.5).

### 8.3 Ordered keys within a stratum

Sort descending on `(K, O, P)` (SCI-09; key order is configurable). Keys are categories except P; never average them.

| Key | Values (best first) | Definition |
|---|---|---|
| K mechanism concordance | `concordant` > `partial` > `not_assessable` > `discordant` | From the gene's `GeneMechanismAssessment` and the declared concordance table. Within a stratum the variant class is constant, so K varies only through declared **sub-class features**: PTV NMD status (`triggers`/`escapes`), missense/in-frame functional-domain membership. Missing feature → `not_assessable`. Class-level concordance (e.g. PTV vs missense under LoF) is reported per stratum and feeds panel quotas (§8.5), not within-stratum order. Class-level inference, never exact-allele evidence (INV-17). Unknown mechanism ⇒ `not_assessable` everywhere |
| O observation | `case_controlled` > `case_uncontrolled` > `related_condition` > `none` | Exact allele observed in an approved disease cohort with a control denominator available / without one / observed only in a related condition / no accepted observation. `none` records whether observation sources were searched |
| P predicted impact | raw value of the stratum's declared predictor, higher = more disruptive; missing last | One predictor per stratum (DEC-06). Missing, unsupported or failed predictions sort after all present values in a separate `P=missing` group, never as 0 |

Ties (DEC-05): rows equal on every key share a dense `tier_rank` and a `tie_group_id`; display order inside a tie group is genomic coordinate, which carries no priority meaning. The run reports, per stratum, tie-group sizes and whether each key was degenerate (constant over the stratum). A degenerate key must not be described as having influenced the order.

### 8.4 Exported but not ranked

Population frequency, ClinVar classification and review status, domain context, other predictions (including AVI when it is not the stratum's declared predictor), editing status and research coverage are exported per row and do not change the order. AVI and AlphaMissense never both act as keys (DEC-12).

### 8.5 Panel selection

Separate operation: per-stratum quotas informed by class-level mechanism concordance, control counts, editing constraints and total capacity. Not implemented in M0; a requested panel size is recorded and reported as unsupported until the panel policy exists.

## 9. Experimental feasibility

R-13: distinguish `design_found`, `no_design_under_config`, `unsupported`, `not_assessed`, `failed`, and `measured`. A design with a PAM is not an efficiency estimate. Off-target/bystander results require a declared tool/config and genomic background; do not synthesize numbers.

`no_design_under_config` describes the selected engine/settings, not impossibility under all editors. Measured efficiency includes assay, host cell, target genotype, replicate count, uncertainty and batch metadata. A final experimentally eligible ordering requires a separate declared selection policy; M1 does not fabricate it when host/design parameters are missing.

## 10. Output contract

R-14: required run artifacts:

| Artifact | Purpose |
|---|---|
| input.yaml, resolved_entities.json | Requested and resolved scope |
| candidate_variants.jsonl, candidate_mentions.jsonl, exclusions.jsonl | Universe and identity decisions |
| annotations.jsonl | Versioned structured source results |
| evidence.jsonl, validation.jsonl | Accepted and rejected claims with decisions |
| datasets.jsonl | Accession and sample/genotype relevance |
| features.jsonl, ranked_variants.csv | Features, K/O/P keys, stratum tier ranks and tie groups |
| provenance.json, manifest.json | Queries, hashes, versions, capabilities, budgets, run status |
| report.md | Evidence cards, conflicts, coverage, limitations, and review queue |

Optional Parquet mirrors are allowed; JSONL remains the M0 interchange contract. Export `primary_channel`, roles, stratum, tier rank, tie group, normalized allele, transcript effects, K/O/P, characterized status, ranking profile/status, editing status, research coverage, contributing evidence IDs, warnings and source snapshot IDs. Transcript/protein aliases may be multiple: no duplicate allele rows merely for transcript effects.

Every report must state: scope, route completeness, cutoff, unavailable sources, number deeply researched, provisional policy status, and whether editing was actually assessed. Explanations cite evidence IDs and locator links and distinguish source observations from interpretation.

## 11. Failure and reproducibility contract

R-15: source absence, retrieval failure, unsupported capability, and negative search are separate statuses. Required identity/reference failures stop ranking of affected alleles; optional source failures allow `partial` runs with explicit coverage. No silent fallback to invented annotations.

R-16: deterministic replay reuses saved normalized objects and accepted model outputs; it invokes neither network nor LLM and reproduces biological content, ordering, and stable content hashes exactly. Timestamps/run IDs may differ and are excluded from content equality. Re-extraction on saved documents calls the LLM again and is measured for stability; it has no promise of exact identity. Refresh may change both sources and model outputs.

R-17: log input/config/code/resource hashes, provider/model returned identifiers, prompt/schema versions, query parameters, document/snapshot hashes, budgets, validation decisions, and source capability checks. If exact model weights are not selectable, record that limitation; temperature/seed are not an exact reproducibility guarantee.

## 12. Scientific decisions remaining

| ID | Question | Build default / effect |
|---|---|---|
| SCI-01 | Evidence triage versus prospective effect discovery? | Resolved: nomination (DEC-01) |
| SCI-02 | Include related neurodevelopmental conditions? | Annotate as related, never upgrade disease match automatically |
| SCI-03 | Transcript/isoform exceptions? | MANE-oriented display when available; preserve every supported effect; identity is genomic |
| SCI-04 | How compare coding, splice and regulatory mechanisms? | Coding/splice MVP; no universal prediction score |
| SCI-05 | Which editing engine/background? | Import designs or return not_assessed |
| SCI-06 | Which controls and panel quotas? | Explicit user candidates only until a policy is declared |
| SCI-07 | Association thresholds and conflict resolution? | O=`case_controlled` requires a control denominator; retain unresolved conflicts |
| SCI-08 | What counts as experimental success? | Separate effect, editing/QC, and design utility in EVAL.md |
| SCI-09 | Key order within a stratum | `(K, O, P)`, configurable |
| SCI-10 | Predictor per consequence stratum | missense AlphaMissense; splice/synonymous splice Δ; PTV none in v0.1 (NMD status enters K) |
| SCI-11 | Exclude common alleles from discovery? | Flag only; no threshold yet |
| SCI-12 | Panel quotas and controls | Panel step not in M0 |
| SCI-13 | SETD1A scPRIME batch as ground truth? | Case study; MAVE/DMS cross-gene benchmark is primary |

The builder must report these choices and their implementation impact; infrastructure work can proceed without pretending they are resolved.

## 13. Acceptance gates

R-18: M0 passes identity adversarial tests, nomination policy tests (strata, ties, degeneracy, characterized routing), schema invariants, rejection ledger tests, duplicate-invariance tests, deterministic replay, and offline export checks. M1 additionally demonstrates real-source adapters and benchmark execution. Engineering gates are required; proposed numerical scientific targets are not evidence of achieved performance. EVAL.md owns metric definitions and denominators.

R-19: public-data MVP may process public variants and public aggregate evidence. Private scPRIME results used for local evaluation stay outside public retrieval/model prompts unless explicitly authorized. No credentials in outputs or logs.

R-20: coding agents must not silently modify scientific semantics. Policy/schema changes require a written proposal, migration implications, and the affected evaluation. Review instructions and implementation issues are in IMPLEMENTATION_PLAN.md and REVIEW_GUIDE.md.

# Variant Evidence System — REVIEW_GUIDE

Version: 2.1.0-draft | Date: 2026-10-05

Purpose: give another coding/scientific review agent a falsifiable task. This package is a draft specification, not an implemented or externally validated system. Review the actual files rather than conversational summaries.

## 1. Package reading order and authority

| File | Owns | Expected reviewer question |
|---|---|---|
| SPEC.md | Scope, scientific boundaries, ranking definition, requirement IDs | Are we solving the correct experimental decision? |
| EVIDENCE_SCHEMA.md | Objects, statuses, identity/claim semantics, invariants | Can the data distinguish the scientifically different cases? |
| ARCHITECTURE.md | Host workflow, permissions, state, retries, snapshots | Can code enforce those distinctions? |
| EVAL.md | Units, labels, metrics, leakage, release gates | Could the tests disprove utility or correctness? |
| IMPLEMENTATION_PLAN.md | Milestones, issue order and acceptance | Can a builder complete the first slice without inventing policy? |

SPEC owns scientific behavior; schema owns object types; architecture owns execution; eval owns denominators. If these disagree, report a contradiction rather than selecting the convenient file. Source adapters must verify current official documentation independently.

## 2. Required review output

Return a written REVIEW.md with:

1. Decision: ready for M0 / ready with concrete fixes / blocked by specified contradictions.
2. Findings table: ID, severity, file/section, violated requirement, concrete failure example, proposed fix and affected tests.
3. Scientific decisions needing owner input, distinguished from delegated engineering choices.
4. Source/API feasibility audit: verified URLs/releases/capabilities/date; do not infer an endpoint exists from a proposed interface.
5. Minimal schema/policy changes with migration implications.
6. Suggested first implementation tasks and scope cuts.
7. Remaining uncertainty and assumptions even if no blocker is found.

Do not automatically implement or rewrite the package during a review-only task. Produce proposed diffs where useful. The owner can subsequently assign implementation.

## 3. Severity rubric

| Level | Meaning | Example |
|---|---|---|
| P0 | Can create false identity/evidence, violate permissions, or invalidate reproducibility | Protein shorthand guessed as exact genomic allele |
| P1 | Can systematically distort rankings/evaluation or prevents required implementation | Candidate filtering hides missed gold alleles; inconsistent schema status |
| P2 | Maintainability/performance/clarity issue with a defined workaround | Inefficient repeated source queries |
| P3 | Optional preference or future feature | Framework or UI suggestion |

Each P0/P1 finding needs a concrete scenario, not just “this could be risky.” Absence of current measured results is expected; invented certainty is a finding.

## 4. Scientific challenges

- Resolved in 2.1 (DEC-01): the objective is nomination. Challenge instead the `(K, O, P)` key order, the mechanism concordance table, and whether `characterized` routing to controls is right.
- Are the proposed D tiers sensible? Direct cohort observation can outrank a related-condition clinical assertion, but a single observation is not association. Challenge this policy explicitly.
- Can exact-allele functional results be accepted when patient background confounds causality? What should the comparator/quality predicates require?
- Is “functional evidence exists” different from “the variant impairs function” everywhere in the schema, rank and report?
- Does gene-level haploinsufficiency/burden evidence become variant-level support accidentally?
- Are measured negative effects preserved without implying benignity? Are consequence-only predicted LoF and demonstrated LoF kept distinct?
- Are related disease, cell type and developmental context matches represented honestly?
- Should eventual discovery ranking use prediction/context strata, Pareto selection, explicit quotas or a validated learned model? The baseline should not preclude these.
- Could literature coverage or search budget masquerade as novelty? Does exploration allocation mitigate but not conceal that bias?
- Are patient counts, cohort independence, clinical submitters and paper duplicates disentangled?
- Does a prediction composite reuse another predictor and risk double counting? Are raw prediction scores ever compared across unsupported contexts/classes?
- Are controls defined by evidence/experimental role rather than assumptions about synonymous/common variants?
- Can existing scPRIME data evaluate allele outcomes independently of edit efficiency and panel ascertainment?
- What scientific claim would a successful held-out benchmark support, and what would remain unproven?

## 5. Engineering challenges

- Does canonical identity encode exact reference accession and normalization convention? Are HGVS and VCF conventions distinguished?
- Are ambiguous transcript/protein/rsID mappings represented one-to-many without false merging?
- Is variant-aware assembly mapping reference-validated on both sides?
- Can literature discovery add alleles safely within a closed, bounded expansion policy?
- Are stage checkpoints invalidated by changed dependencies and resume operations idempotent?
- Does each tool declare actual supported variants/builds/modalities and retrieval status?
- Are mechanical validation and semantic entailment separate? Can a second model overrule host reference integrity?
- Do accepted claims reference actual source snapshot locations? Are parser offsets reproducible?
- Is model output replay independent of live provider state? Does re-extraction expose variability?
- Can transient retries, schema repairs, or duplicate calls evade budgets?
- Are denied/failed/inaccessible sources distinct from successful empty results?
- Does a rejected claim remain in the audit ledger while being excluded from features?
- Does sample-level GEO evidence establish genotype rather than just gene mentions?
- Can editing status distinguish found design from measured efficiency and no-design-under-settings from general impossibility?
- Can every feature predicate be traced to an accepted link and policy version?
- Do reports distinguish automatic acceptance, human review, inference and prediction?

## 6. Evaluation challenges

- Are candidate universe, source cutoff, consequence scope and accessible corpus fixed before recall is calculated?
- Does ranking evaluation include discovery omissions rather than only user-supplied candidates?
- Are gold labels independently reviewed and aligned with the intended objective?
- Are unjudged variants kept out of false-negative counts?
- Are publication/cohort/query clusters kept within splits and confidence intervals?
- Does a ClinVar-based gold task explicitly identify ClinVar lookup as a baseline?
- Are agentic query selection and fixed-query extraction compared under matched budgets?
- Are pre-verification unsupported proposals reported alongside accepted-claim precision?
- Are deterministic replay, extraction variability and live source drift evaluated separately?
- Would any pass criterion allow a fabricated benchmark claim or label a one-gene case study as general validation?

## 7. Copyable reviewer prompt

> Review all files in this package as a scientific/engineering specification for a disease+gene variant screening evidence system. Do not implement code yet. Try to break allele identity, evidence directness, disease/cell-context matching, recurrence deduplication, source provenance, dataset genotype assignment, editing status, bounded agent execution, replay and evaluation validity. Challenge the provisional nomination_v0.1 keys, strata and mechanism concordance table rather than treating them as scientific truth. Return REVIEW.md with P0–P3 findings, file/section and requirement/invariant IDs, concrete counterexamples, proposed fixes and affected tests. Separate unresolved scientific choices from routine engineering choices. Independently verify any claimed API/source capability using primary documentation. Do not invent variants, citations, biological scores, labels or performance. State whether ISSUE-001–004 can proceed and list the smallest changes needed.

## 8. Response to review

After review, create a decision log recording accepted/rejected recommendations with rationale. Version the package when scientific behavior or contracts change. Do not tune scientific policies just to make synthetic tests pass; rewrite the expected policy tests only when the actual policy changes. Preserve the prior documents and benchmark versions so a reviewer can reconstruct what changed.

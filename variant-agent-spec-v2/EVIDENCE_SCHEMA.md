# Variant Evidence System — EVIDENCE_SCHEMA

Version: 2.1.0-draft | Date: 2026-10-05

This file owns object semantics, enums, identity and cross-object invariants. Implement with strict Pydantic models or equivalent, generate JSON Schema, reject unknown fields, and export UTF-8 JSON/JSONL. No implementation or validated biological records are included in this document.

## 1. Schema rules

- Every persisted object has `schema_version`, a stable object ID, and version/revision when mutable concepts change. Revisions are append-only.
- Times use ISO-8601 UTC with timezone; resource versions must be explicit or `unknown` with reason.
- IDs are strings; counts/positions are integers; probabilities have declared bounds; effect sizes include units/scales. Avoid implicit coercion of numeric strings.
- `null` means no observed value, with a separate status. It does not mean zero, false, wild type, absence of phenotype, or benign.
- Arrays default to empty only for real empty sets; an unqueried set has a status showing it was unqueried.
- Supporting excerpts are concise for audit; store source locations/content digests rather than reproducing whole copyrighted papers in reports. Full snapshots are retained only when permitted.
- All enum changes and evidence policy changes are versioned. Unknown source values preserve raw strings and map to `unknown`, not a convenient known category.

## 2. Shared types

### 2.1 Observation status

`value_present | not_searched | searched_not_found | not_available | access_denied | not_applicable | unsupported | failed`

`value_present` can contain a measured zero or a source-reported false. `searched_not_found` requires a logged successful query with scope and retrieval event. `not_available` requires a documented missing field/resource, not a network exception.

Generic `Observed[T]`:

| Field | Type / requirement |
|---|---|
| status | ObservationStatus; required |
| value | T or null; non-null iff status=value_present, except explicitly nullable source payload types |
| source_snapshot_ids | list[str]; nonempty when a source-derived value or successful negative query is asserted |
| reason_code | str or null; required for failed, unsupported, not_applicable, access_denied |
| retrieval_event_id | str or null; required for searched_not_found |

A biological statement is not represented by missingness. “No phenotype detected in assay X” is a value-bearing experiment outcome.

### 2.2 Identity match

`exact_genomic | transcript_supported_exact | protein_equivalent_only | same_residue | same_domain | mutation_class | gene_only | locus_only | pathway_only | ambiguous | no_match`

Only exact_genomic and transcript_supported_exact can establish an exact genomic allele link. A validated genomic allele derived from transcript HGVS becomes exact after genomic normalization. Protein equivalence alone cannot.

### 2.3 Claim scope and evidence kind

Scope: `variant | residue | domain | mutation_class | gene | locus | pathway`.

Kind: `patient_observation | association | segregation | clinical_classification | functional_experiment | computational_prediction | gene_mechanism | dataset_description | review_summary | hypothesis`.

These are independent. A prediction can concern an exact variant without being an experiment. A gene experiment can be high quality without becoming variant-specific. Hypotheses may appear only in interpretation outputs, never as accepted ranking evidence.

### 2.4 Validation state

`draft | mechanically_valid | accepted_auto | accepted_human | quarantined | rejected | superseded`

Effective ranking state: accepted_auto or accepted_human, except where SPEC explicitly requires human confirmation. Rejected/quarantined/superseded claims never contribute. Mechanical validity alone is insufficient.

Human_review: `unreviewed | accepted | rejected | needs_revision`; independent of automatic check results. Effective decision resolution records which reviewer/version is authoritative.

### 2.5 GeneMechanismAssessment (DEC-03)

| Field | Type / semantics |
|---|---|
| gene_id | HGNC/Ensembl ID |
| disease_concept | ontology ID + relationship to the query disease |
| mechanism | `loss_of_function | gain_of_function | dominant_negative | mixed | unknown` |
| status | ObservationStatus; `unknown` mechanism with `searched_not_found` differs from `not_searched` |
| basis | list of {source, record ID/version, field, value, locator}: e.g. ClinGen haploinsufficiency score, gnomAD LOEUF, Gene2Phenotype allelic requirement/mechanism, curated literature claim |
| decision | `rule` (declared mapping policy), `human`, or `imported`; policy ID/version |
| concordance_table_id | ID/digest of the mechanism × variant-class table used for K |

A mechanism is a gene–disease property. It sets K only through the concordance table and is never copied into an allele's evidence links. Constraint metrics alone (e.g. low LOEUF) indicate intolerance, not disease mechanism; a rule that maps them to `loss_of_function` must be declared and is reported as such.

## 3. RunRequest / ResolvedScope / RunManifest

RunRequest mirrors SPEC.md input. Required subobjects: gene, disease, screen_context, candidate_scope, objective, execution. Record the original request untouched, then materialize defaults into ResolvedScope.

ResolvedScope:

| Field | Type / semantics |
|---|---|
| gene | GeneEntity |
| disease | DiseaseEntity |
| reference_bundle_id | str; pinned assembly accession + FASTA checksum |
| transcript_bundle_id | str; pinned annotation/MANE resource |
| allowed_consequences | list[str] from declared ontology/policy |
| allowed_allele_types | list[str] capability-scoped |
| cell_context | species, cell_type text/ontology ID, developmental stage, tissue, culture/model details; unknowns explicit |
| editing_context | editor/version, PAM, host genotype/haplotype, target zygosity, design config ID; values may be unassessed |
| candidate_routes | list[route configs] with enabled/required flag and scope |
| disease_expansion_policy_id | str |
| ranking_policy_id/status | str / provisional or evaluated |

GeneEntity: stable HGNC/Ensembl or species-appropriate IDs, current symbol, approved aliases, organism, chromosome/strand, transcript IDs with versions, selected display transcript and selection reason. Genome membership and gene link are not identical for distal regulatory variants.

DiseaseEntity: original input, resolved ontology ID/release/name, approved synonyms, related concepts and matching policy. Each match has `relationship`, `mapping_basis`, and optional human decision. Ambiguous input produces candidates and `needs_resolution`, never a guessed match.

RunManifest requires run ID, spec/schema versions, code commit or source-tree digest, dependency lock digest, input/resolved-scope digest, source capability audit, prompt/model/schema/config digests, evidence/candidate snapshot digests, counters, stop reasons, stage statuses, terminal run status, privacy scope and artifact content hashes. Provider fields include requested and returned model identifiers, configuration actually supported, unavailable pinning information and request IDs where provided.

## 4. Variant and alias contracts

### 4.1 CanonicalVariant

| Field | Type / requirement |
|---|---|
| variant_id | str; hash of declared canonical identity material |
| identity_version | str; normalization/hash convention version |
| assembly | GRCh38 for M1 canonical output |
| reference_accession | versioned sequence accession |
| reference_bundle_id | pinned resource ID/checksum |
| chrom | declared contig convention; not the sole identity |
| pos_1based | int >=1, VCF convention |
| ref, alt | nonempty uppercase literal DNA alleles, unequal; M1 excludes symbolic ALT |
| allele_type | snv, deletion, insertion, delins, mnv; capabilities gate accepted types |
| normalization_status | reference_validated only in ranked registry |
| normalization_tool/version/config_digest | required |
| gene_links | list[GeneLink]; includes mapping method and evidence, not only symbol |
| transcript_effects | list[TranscriptEffect] |
| alias_ids | list[str] |

Variant ID convention: SHA-256 of canonical UTF-8 JSON with sorted keys, no whitespace, keys `identity_version`, `assembly`, `reference_accession`, `pos_1based`, `ref`, `alt`. Prefix `var_`. Reference_bundle_id is provenance; reference_accession anchors identity, and the exact reference bundle must validate it. Normalization convention changes require migration rather than silently reusing IDs.

HGVS is a representation/alias, not the primary key. TranscriptEffect fields: transcript ID/version, gene ID, HGVS c/p, consequence list, exon/intron, isoform/protein accession, mapping tool/config, selected_display flag, and status. Effects can disagree across isoforms; retain all and display selected policy.

### 4.2 VariantMention and IdentityDecision

VariantMention stores mention ID, source snapshot/locator, raw text, source-reported gene/assembly/transcript, parser output, role suggestion, and discovery route. Role suggestions from the model are not accepted controls without a declared policy.

IdentityDecision stores mention ID, status (`resolved | ambiguous | insufficient | ref_mismatch | unsupported | unmapped`), candidate variant IDs/representations, mapping method, resource versions, match justification, validation errors, decision maker and timestamp. `resolved` requires one supported concrete genomic allele, not one arbitrarily selected alternative.

Alias stores raw value, alias type, publication-specific context, transcript/assembly where applicable, mapping relation and decision ID. An rsID may map to multiple ALT alleles. Aliases can be one-to-many, so enforce uniqueness on canonical allele identity rather than alias string.

### 4.3 CandidateMembership

Fields: run ID, variant ID or unresolved mention ID, route ID, roles (`discovery | positive_control | negative_control | comparison | user_priority`), primary_channel (`discovery | control`), inclusion/exclusion reason, scoped eligibility, discovered stage, annotation/research status, and selection probability/stratum if sampled for research.

Each route logs counts and exclusions. A control label carries its own justification and source/human selection; pathogenic/benign terminology is never inferred solely from experiment role.

## 5. SourceSnapshot and supporting locator

SourceSnapshot:

| Field | Type / requirement |
|---|---|
| snapshot_id | str, content addressed |
| source | str, e.g. ClinVar, EuropePMC, GEO, AlphaGenomeAtlas, synthetic_fixture |
| source_record_id/version | str and explicit version or unknown reason |
| source_uri | stable public URI/accession; no credentials or signed token URLs |
| retrieved_at | UTC timestamp |
| request | canonical query/parameters, source query translation if available |
| adapter_version, capability_version | required |
| raw_content_digest | SHA-256 or explicit restricted-content status |
| text_content_digest, text_parser_version | required for normalized text locators |
| storage_policy | retained, metadata_only, restricted |
| access_status | available, abstract_only, full_text_unavailable, supplement_unavailable, controlled_access, failed |
| license_status | known terms/release or unknown; do not assume permission |
| is_synthetic | boolean; required, cannot mix synthetic and scientific evaluation silently |

SupportingLocator is a discriminated union:

| kind | Required fields |
|---|---|
| json_field | snapshot_id, JSON Pointer, expected value digest |
| text_span | snapshot_id, normalized text digest, parser version, start/end character offsets, concise excerpt |
| table_cell | snapshot_id, table/sheet ID, row/column coordinates, cell digest, parser version |
| figure_panel | snapshot_id, page/figure/panel locator, image digest, caption/text span when available; semantic/human review needed |

Offsets are zero-based Unicode code-point indices, start inclusive and end exclusive, in the pinned normalized text. Do not mix byte and character offsets. PDF page labels and rendered page indices are stored separately where needed. A URL without a specific source location does not satisfy the claim contract.

Bibliographic metadata: PMID/PMCID/DOI and title/authors/year/version when available. Non-PubMed sources remain legitimate. Link preprint/published/correction/retraction relationships; never treat publication versions as independent replication by default.

## 6. Claim, link and validation: normalized evidence model

One underlying gene-level finding can be relevant to many candidate alleles. Store it once as a Claim; use EvidenceLink to express applicability. This prevents repeated copies from becoming replication.

### 6.1 Claim

| Field | Type / semantics |
|---|---|
| claim_id, revision | stable claim ID + positive int |
| evidence_kind, scope | enums above |
| subject_gene_ids | list[str] |
| source_reported_variants | list[mention IDs]; may be empty for gene claims |
| disease_concepts | list[ontology/text IDs + relationships] |
| claim_text | concise factual statement attributable to the source |
| outcome | structured direction/result or null; no implicit priority polarity |
| context | organism, model, cell type, tissue, developmental stage, host genotype |
| functional_assay | FunctionalAssay or null |
| association | AssociationRecord or null |
| clinical_assertion | ClinicalAssertion or null |
| prediction | PredictionRecord or null |
| locators | nonempty SupportingLocator list |
| source_snapshot_ids | nonempty list[str] |
| cohort_ids, experiment_ids | list[str], overlap flags when unknown |
| publication_group_id | str or null |
| extraction | method human/deterministic/llm, model request/output IDs, prompt/schema versions, optional self-reported certainty |
| supersedes | claim revision or null |

Self-reported model certainty is diagnostic only: it cannot become an evidence strength, biological probability, or numeric ranking feature.

### 6.2 EvidenceLink

Fields: link ID, claim ID/revision, variant ID (nullable only for unassigned claims), identity_match, identity_decision IDs, applicability scope, disease_relationship, cell_context_relationship (`exact | approved | related | other | unknown`), justification locators, inference flag and effective validation decision ID.

An exact experiment link requires the source to establish the tested allele, not merely cite a paper that mentions it elsewhere. Tested construct, patient allele and source-reported protein change may all need separate mapping. Isogenic editing versus patient-derived carrier comparison is stored in FunctionalAssay, not inferred from the word “variant”.

### 6.3 ValidationDecision

Fields: decision ID, target claim/link revision, check policy/version, per-check status/reasons, mechanical/identity/citation/semantic results, reviewer type, reviewer model/request or human pseudonymous ID, human_review status, effective state, timestamp and predecessor decision.

Per-check enum: pass, fail, unknown, not_applicable. Semantic unknown→quarantined, never accepted by default. Export unsupported-claim proposals and reasons as well as accepted evidence, so hallucination evaluation does not inspect only the survivors.

## 7. FunctionalAssay

Required for functional_experiment claims:

| Field | Meaning |
|---|---|
| assay_name/type | Reporter, biochemistry, imaging, transcriptomics, electrophysiology, etc. |
| tested_allele_mentions | Identity-linked descriptions of what was actually tested |
| perturbation_type | exact_allele_edit, allele_construct, carrier_comparison, gene_knockout, knockdown, overexpression, other |
| genotype/zygosity | HET/WT/HOM or explicit unknown; phase where reported |
| comparator | WT/isogenic/rescue/vehicle/other, with description |
| endpoint | Measured quantity, units, direction and timepoint |
| outcome_type | altered, no_detected_effect, rescue, mixed, inconclusive |
| effect_size/CI/statistics | Observed fields with scale, n and test/multiplicity where reported |
| biological_replicates/technical_replicates | Separate observed counts; cells are not automatically biological replicates |
| experimental_context | Cell/model/tissue/developmental stage, donor and batch where public |
| confound_flags | Patient genetic background, edit efficiency, clone/batch, viability, differentiation, multiple edits, etc. |
| quality_review | Comparator adequacy, sensitivity, reporting completeness and assessor/version |

Gene knockout results cannot produce exact-allele functional evidence. A carrier-vs-control experiment can be exact-allele-associated observation but is not necessarily a causal installation experiment. Comparison and confounds must be exposed. Mechanism inferred from consequence, e.g. predicted LoF, belongs to prediction/annotation rather than a measured result.

## 8. AssociationRecord, recurrence and ClinicalAssertion

AssociationRecord stores study design, trait definition, cohort/sample counts, allele counts, allele frequency, effect/statistic/CI/test/multiplicity, segregation/de-novo details, source-reported thresholds, analysis unit (`allele | mutation_class | gene | locus`), independence/overlap and whether the exact allele was analyzed.

RecurrenceSummary keeps separate fields:

- unique_individuals, independent_families, independent_cohorts: Observed counts with deduplication basis and uncertainty.
- source_reported_observations: raw counts by study, not summed without independence checks.
- clinvar_scv_count and submitting_organizations: database metadata.
- publication_count and publication_groups: document metadata.

Do not infer absent denominators, independence or penetrance. Multiple records may refer to one patient. Burden-test p-values remain mutation-class/gene associations; each constituent allele does not inherit that p-value as allele-specific evidence.

ClinicalAssertion fields: VCV/RCV/SCV accessions and versions where available, condition labels/IDs, disease match, assertion type (`germline | oncogenicity | somatic_clinical_impact`), raw and mapped classification, textual review status, stars as auxiliary metadata, submitter IDs, assertion method, last evaluated, contributing status, aggregate conflict, source citations/locators and release.

Do not compress review status to stars alone or pool germline and somatic assertions. The adapter records both aggregate condition-specific records and submissions; provenance links classifications to their supporting sources without counting them again as independent experiments.

## 9. AnnotationRecord and prediction channels

AnnotationRecord: annotation ID, variant ID, source snapshot, annotation kind, status, raw payload/version, normalized payload, transformation/config digest, and transcript/feature context.

Population annotation fields: dataset/release, exome/genome or joint scope, assembly, AC, AN, AF, group-specific AC/AN/AF, filtering/callability/coverage, homozygote/hemizygote counts where available, and source field definitions. Absence from an uncallable position does not justify AF=0.

Protein/domain annotation: source protein accession/version/isoform, transcript mapping, residue/range, domain IDs and boundaries, residue mapping status, disorder/structure confidence where available. Low structural confidence is not pathogenicity, and domain membership alone is not a measured effect.

PredictionRecord:

| Field | Type / semantics |
|---|---|
| predictor/channel | AlphaGenome molecular modality, Atlas AVI, AlphaMissense, splice model, etc. |
| model_version/resource_release | required or explicit unavailable reason |
| supported_allele_class | declared source capability |
| feature/tissue/ontology IDs | raw and normalized context; not inferred target-cell match |
| input_interval/assembly/reference | full supported scoring context |
| score_name/raw_value/units | source-defined; may be multiple track outputs |
| score_direction | source-defined interpretation, unknown if unavailable |
| transformed_value | optional; requires transform ID and calibration/reference distribution |
| component_dependencies | known constituent predictors/features; unknown explicitly recorded |
| affected_gene_assignment | method, target ID, uncertainty; nearest gene not proof |
| status | observation status |

AVI is retained as a source-defined prediction with its release-specific formula. If an aggregate includes a protein predictor, adding that predictor again as an independent term risks double counting. No promise that AVI is available for every indel or that modality outputs and composite values share calibration.

## 10. DatasetRecord and sample-level links

DatasetRecord fields: source/accession, title, organism, study design, assay, sample count, disease/model/cell metadata, public availability of raw/processed files, controlled-access status, related publications, source snapshots, metadata retrieval status, dataset-level relevance and its basis.

SampleVariantLink fields: dataset/sample accession, variant ID/mention, genotype description, evidence locator, match status, control/comparator availability, sample cell/tissue context, provenance and validation. Genotype-confirmed carrier samples do not imply the allele was experimentally introduced.

Dataset relevance classes: `exact_allele_confirmed | gene_perturbation | disease_context | gene_expression_only | unrelated | unknown`. “SETD1A expression was measured” is gene_expression_only, not SETD1A perturbation. Raw availability, processed matrix availability and access permission use Observed booleans; unknown is not false.

A dataset contributes target-context exact-allele evidence only with validated sample linkage and target-context matching, not merely a GSE title or a gene mention. Dataset availability does not imply a functional result unless an actual analysis/source reports one.

## 11. EditingAssessment

Fields: variant ID, assessment status (`design_found | no_design_under_config | unsupported | not_assessed | failed | measured`), design engine/version, editor/PAM, reference/local haplotype/host genotype, target genotype, search parameters, returned design IDs, pegRNA/PBS/RTT design fields where generated, tool-defined scores with units/calibration, off-target/bystander method, measured efficiency/uncertainty/context where present, provenance and limitation codes.

`design_found` requires at least one engine-produced or imported design and complete design provenance. `measured` requires actual experimental source/context. A failed assessment cannot emit expected_efficiency=0. A found design cannot emit a calibrated success probability without a declared validated model.

## 12. FeatureRow and RankedRow

FeatureRow required fields: variant ID, frozen evidence digest, ranking policy ID/digest/status, primary_channel, roles, consequence stratum, characterized status, K (category), O (category), P (raw value + predictor ID or missing status), per-feature evidence IDs and satisfied/unmet predicates, source coverage/missingness summary, conflict flags, consequence group, prediction annotation references, editing assessment ID/status.

RankedRow additionally has stratum-specific dense `tier_rank` (1-based), ordering tuple, `tie_group_id` and size, coordinate display index (no priority meaning), comparison_scope ID, explanation template/version, review flags and evidence card ID. No `final_probability` or fake global score field in baseline M1.

Coverage reports which sources and deep-research tasks ran. “O=none with not_searched” and “O=none with searched_not_found” must be distinguishable. Ranked rows cannot cite an evidence ID absent from the frozen bundle or a decision that is not accepted.

## 13. Cross-object invariants: required tests

| ID | Invariant |
|---|---|
| INV-01 | Ranked allele has a unique normalized identity and matching reference sequence |
| INV-02 | Exact-allele accepted link includes an identity decision and source-specific tested/observed allele support |
| INV-03 | Each accepted factual claim has inspectable supporting locator(s) to the pinned source snapshot |
| INV-04 | Gene/domain/mutation-class links never produce exact-allele functional/association features |
| INV-05 | Clinical classification retains condition and assertion type; submitting organizations never become patient counts |
| INV-06 | Prediction stores score semantics; unsupported/missing values cannot become numerical zero |
| INV-07 | No-detected-effect observation retains tested assay and comparator; cannot become universal benignity |
| INV-08 | Dataset exact-allele status requires sample/genotype locators |
| INV-09 | Rejected/quarantined evidence cannot contribute feature predicates |
| INV-10 | Duplicate paper/version/claim insertion cannot increase baseline tiers or replication counts |
| INV-11 | O=`case_controlled` requires an accepted exact-allele observation in an approved disease cohort and a control denominator with status value_present |
| INV-12 | Synthetic fixtures never appear as real-source benchmark evidence without explicit synthetic labeling |
| INV-13 | Source/prompt/policy revisions cannot overwrite the bundle used by a completed run |
| INV-14 | Every exported rank references policy and evidence digests; every source failure remains observable |
| INV-15 | Multiple transcript consequences produce one genomic-allele ranked row |
| INV-16 | Every meaningful null/unknown or negative search has a declared status and scope |
| INV-17 | Mechanism concordance (K) never creates or upgrades an exact-allele evidence link |
| INV-18 | No rank compares alleles across (channel, consequence stratum) |
| INV-19 | Rows equal on all keys share tier rank and tie group; degenerate keys are reported per stratum |

## 14. Synthetic example, never scientific evidence

This abbreviated proposal demonstrates orthogonal dimensions. It is intentionally not a fully validated canonical variant. A fixture builder must supply the real synthetic reference and matching locators before accepting it.

```json
{
  "schema_version": "2.0.0-draft",
  "claim_id": "fixture_claim_001",
  "revision": 1,
  "evidence_kind": "functional_experiment",
  "scope": "gene",
  "subject_gene_ids": ["FIXTURE_GENE_1"],
  "source_reported_variants": [],
  "disease_concepts": [],
  "claim_text": "Gene knockout changed the reporter signal in fixture cells.",
  "outcome": {"outcome_type": "altered", "direction": "decreased"},
  "context": {"organism": "synthetic", "cell_type": "fixture_cells"},
  "functional_assay": {
    "assay_name": "synthetic_reporter",
    "tested_allele_mentions": [],
    "perturbation_type": "gene_knockout",
    "genotype": "unknown",
    "comparator": {"type": "wild_type", "description": "fixture WT"},
    "endpoint": {"name": "signal", "units": "arbitrary"},
    "outcome_type": "altered",
    "confound_flags": []
  },
  "association": null,
  "clinical_assertion": null,
  "prediction": null,
  "locators": [{
    "kind": "json_field",
    "snapshot_id": "fixture_snapshot_001",
    "pointer": "/results/0",
    "expected_value_digest": "<fixture-builder-computed-digest>"
  }],
  "source_snapshot_ids": ["fixture_snapshot_001"],
  "cohort_ids": [],
  "experiment_ids": ["fixture_experiment_001"],
  "publication_group_id": null,
  "extraction": {"method": "human"},
  "supersedes": null
}
```

Expected outcome: accepted gene-level claim may inform a gene mechanism summary, but must not mark any individual allele `characterized` or change its O/P keys. A synthetic exact-allele counterpart must name its tested allele and validate its identity link.

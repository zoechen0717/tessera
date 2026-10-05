# Variant Evidence System — EVAL

Version: 2.1.0-draft | Date: 2026-10-05

Evaluation separates engineering correctness, evidence extraction, candidate discovery, ranking utility and experimental outcomes. Passing deterministic tests is necessary but does not demonstrate biological utility. No results or gold labels are claimed here.

## 1. Questions and evaluation units

| Question | Unit | Primary measurement |
|---|---|---|
| Can identity be trusted? | Allele mention→genomic allele decision | Exact mapping precision, recall, ambiguity handling |
| Did discovery find known candidates? | Candidate allele in declared universe | Candidate recall by route/class |
| Does a cited source support a claim? | Atomic claim/link | Citation support precision and unsupported proposal rate |
| Was directness classified correctly? | Evidence link | Exact-allele precision and confusion matrix |
| Is GEO genotype/context relevance correct? | Dataset/sample link | Exact-allele dataset precision and class accuracy |
| Does ranking fit the stated objective? | Gene/disease query + candidate set | NDCG/utility, blinded panel preference, expert rubric |
| Is execution reproducible? | Run pair | Replay equality; re-extraction stability |
| Does an experiment show an effect? | Installed allele with biological replicates | Pre-specified assay outcome conditional on editing/QC |

These units must not be exchanged. A gene burden hit is not a gold positive label for all alleles. A biological effect is not automatically disease causality.

## 2. Benchmark stages

### B0: adversarial engineering fixtures

Use an explicit synthetic reference bundle and synthetic papers/JSON records. Include:

- Valid SNV, normalized indel representations, multiallelic rsID and conflicting transcript effects.
- Minus-strand gene with forward genomic alleles.
- Protein-only change compatible with multiple DNA edits; abbreviated frameshift without an allele.
- GRCh37 input with successful mapping, REF mismatch, failed and ambiguous mapping.
- Gene knockout paper that mentions a different exact variant only in introduction.
- Same-residue/different-allele result and same-domain result.
- ClinVar multiple SCVs describing one underlying study; related-condition classification and conflicting aggregate.
- Real-looking citation without a source, real source with unsupported claim, wrong text offsets and altered source snapshot.
- Patient-derived carrier comparison versus isogenic installation; functional no-effect outcome versus no paper found.
- Gene-expression-only GEO metadata versus genotype-confirmed variant samples.
- Atlas unsupported allele and unavailable prediction channel; gnomAD uncallable versus measured AC=0.
- Preprint/published duplicate, shared cohort, contradictory assays, retraction/correction.
- Prompt injection in a supplement and model tool request outside the allowlist.
- Rate limit, transport failure, invalid model JSON, exhausted budget, checkpoint interruption, replay artifact missing.

Expected properties and labels are manually specified from fixture contents. Synthetic IDs must remain visibly synthetic. M0 cannot be described as a real SETD1A demonstration.

### B1: SETD1A/scPRIME development benchmark

The owner supplies a versioned genomic allele list, source-backed labels, experimental metadata, and scope. Do not reconstruct exact genomic alleles from remembered protein shorthands or invent labels. Existing panel choices are useful recall targets, not necessarily an optimal-ranked gold set.

Potential inputs, each independently verified:

- Publicly documented disease alleles and source loci.
- Previously selected experimental alleles and explicit control role.
- Verified direct functional studies, including negative/contradictory findings.
- Search/access labels for real literature and GEO data.
- Editing efficiencies, genotype QC, cell counts and transcriptomic outcomes for retrospective analysis when authorized.

The reader must know whether labels measure disease relevance, perturbation effect, editing success, or historical inclusion. They are different tasks.

### B2: held-out generalization

After development, evaluate additional gene/disease queries spanning consequence classes and mechanisms. Choose dataset size using confidence interval/power needs. A one-gene demonstration supports a case study, not general performance claims.

Split by gene/query, publication family, cohort and temporal cutoff where appropriate. Avoid random row splits that put related claims from the same experiment into train/development and test. Freeze test labels and policy before evaluating alternatives.

## 3. Gold dataset contracts

| File | Required fields |
|---|---|
| gold_variants.csv | allele identity, alias/source locator, eligibility, disease relationship, class, control role, label purpose |
| gold_mentions.jsonl | mention locator, valid allele alternatives, expected mapping/ambiguity status |
| gold_evidence.jsonl | atomic claim, source snapshot/locator, kind, scope, identity match, disease/cell relation, outcome, confounds |
| gold_datasets.jsonl | accession/sample, genotype locator, relevance class, access status |
| gold_ranking_labels.csv | query ID, allele ID, utility rubric dimensions, grade or unjudged, assessor/version |
| split_manifest.json | development/validation/test assignments, cutoff, grouping keys, digest |
| adjudication.jsonl | reviewer disagreements, reason, final decision and revision |

Human reviewers inspect original sources, not just model summaries. For high-impact labels, use two independent reviews and adjudication where feasible. A coding agent may draft annotations, but a draft model output cannot be the gold answer used to declare that same system accurate.

Partial judgment is expected: unjudged is not irrelevant. Report number/coverage of judged candidates. Review accessible and inaccessible sources separately; a paywall affects retrieval coverage, not biological truth.

## 4. Identity metrics

Let gold-resolvable mentions be those with a unique reference-supported allele under the benchmark policy.

- Exact mapping precision = correctly resolved mentions / all mentions the system resolved.
- Exact mapping recall = correctly resolved gold-resolvable mentions / all gold-resolvable mentions.
- Ambiguity handling accuracy = ambiguous/insufficient mentions correctly left unresolved / all gold ambiguous/insufficient mentions.
- False merge rate = distinct gold alleles incorrectly unified / gold allele pairs or merge decisions under the declared denominator.
- False split rate = duplicate representations not unified / gold equivalence groups.

Report both allele-level and mention-level results. Variant-level evidence assigned to a wrong allele is a critical failure even if gene/disease match is correct.

## 5. Discovery metrics

CandidateRecall = gold eligible alleles included in the resolved candidate set / gold eligible alleles within declared scope and cutoff.

Report by candidate route, consequence class, clinical/literature status and source accessibility. Also report recall before and after filtering, unresolved mention count, and top-k recall separately. Candidate discovery failure cannot be hidden by evaluating ranking only among alleles already supplied by the user.

Compare disease+gene discovery with manual-list annotation as separate tasks. Gold unknown outside queried sources prevents universal recall claims; declare benchmark ascertainment and list which candidate routes the gold set came from.

## 6. Claim support and directness metrics

Citation existence precision = proposed source IDs that resolve to the intended source / proposed source IDs. This does NOT measure support.

Citation support precision = judged accepted atomic claims whose cited location supports the claim at its stated identity/scope/context / all judged accepted claims. Support requires more than a matching quote. Report exact denominator and unjudged counts.

Unsupported proposal rate = unsupported judged draft claims / all judged draft claims, before validation. Unsupported accepted rate uses accepted claims only. Record both so a rejecting verifier does not hide extraction problems. Measure rejection recall for unsupported proposals, acceptance precision, and loss of supported proposals.

Directness confusion matrix: gold versus predicted identity/applicability scope, with particular attention to gene→exact allele, protein-equivalent→exact genomic and same-residue→tested allele errors. Exact-allele functional precision = correct exact-allele functional links / all judged accepted exact-allele functional links.

Functional extraction field accuracy covers tested allele, assay/comparator, outcome, model, cell context, biological n, and confounds. Do not score omitted unknowns as wrong when the source does not report the field. Include an explicit unknown/error analysis.

## 7. Dataset metrics

Measure relevance-class confusion matrix, accession validity and exact-allele sample-link precision. Test genotype availability claims against inspectable metadata/sample files. Report inaccessible data separately. A dataset with disease cells is not an exact-allele experiment, and gene expression measurements are not gene perturbations.

## 8. Ranking evaluation

### 8.1 Label definition precedes the metric

For nomination_v0.1, the deterministic ordering can be tested against the key rules. Scientific utility requires independent labels: measured per-allele effects (§8.3), not the K/O/P keys themselves. Do not label alleles with the policy's own keys and claim that sorting them validates utility.

Evaluate alternative objectives separately:

1. Recovery of established disease-related alleles.
2. Selection of understudied alleles with later measured effects.
3. Scientist-rated usefulness of a screen proposal.
4. Editing-compatible panel yield and diversity.

The proposed baseline may be good for (1) and poor for (2); the evaluation must allow that result.

### 8.1a Primary nomination benchmark (DEC-08)

Units: gene (held out) × consequence stratum. Labels: per-allele functional readouts from MAVE/DMS and saturation genome editing (MaveDB and primary publications), with the assay's own effect threshold or a predeclared one, assay identity and direction recorded. Question: within a held-out gene and stratum, does the top-k by `nomination_v0.1` contain more effect-positive alleles than (a) random, (b) the stratum predictor alone, (c) `evidence_priority_v0.1`, (d) frequency-only?

- Restrict to alleles the run could nominate (in scope, normalizable); report coverage.
- Group splits by gene; never tune K/O/P order or the concordance table on test genes.
- Leakage: record whether the predictor saw DMS data in training or calibration; report such genes separately. AlphaMissense and AVI training/calibration sets must be checked at build time.
- DMS assays measure specific molecular/cellular functions; record assay-to-mechanism relevance. A DMS hit is evidence of perturbation, not disease causality.

SETD1A scPRIME (B1) is a retrospective case study: historically selected panel, few cells per allele, no biological replicates; report effect estimates with uncertainty and do not pool it with the DMS benchmark.

### 8.2 Defined metrics

- Recall@k: relevant gold eligible alleles in top k / all relevant gold eligible alleles in scope.
- Precision@k: relevant alleles in top k / k, ONLY when all top-k entries are judged. With partial labels, report judged precision and judged coverage, or use an explicit sampling estimator; never mark unjudged entries negative.
- NDCG@k: use predeclared utility grades with gain `2^grade - 1` and discount `log2(rank+1)`. IDCG is defined over the same fixed candidate/gold universe. With zero total relevance report undefined, not a perfect score.
- Expert review: blinded, randomized output order and fixed source access. Capture rubric dimensions, disagreement and review time, not just preference percentages.
- Panel yield: number installed/passing QC divided by attempted variants, and number interpretable effect-positive alleles divided by successfully installed/QC-passing alleles. Include uncertainty and failed attempts.

Baseline tiers are categorical; correlation with effect size should not imply calibrated distances. Separate discovery/control channels and consequence strata. Report tie counts and k-boundary sensitivity rather than letting a hash tie-breaker suggest scientific superiority.

## 9. Reproducibility and uncertainty

### Exact replay

On a frozen accepted-evidence bundle, code/config/resources fixed, require identical candidate identity set, accepted evidence references, feature predicates, ordering tuples and ranked content digest. Ignore timestamps/run IDs and unstable serialization metadata. No network or model calls may occur. Detect missing artifacts as failure.

### Model re-extraction

Repeat extraction on frozen documents, same configured model/prompts at least five times for an initial variability study. Exact output is not assumed. Compare accepted-claim support, directness classifications, normalized claim signatures, candidate/evidence set overlap and top-k overlap. Measure extraction and validation variability separately.

Evidence signatures exclude arbitrary model IDs/timestamps but include source snapshot/locator, canonical subject, assay/outcome and scope. Stable signatures need a declared canonicalization policy; fuzzy semantic similarity alone can conceal meaning changes.

Top-k Jaccard = intersection / union. Report k relative to available candidate count and primary channel. Rank correlations require a common declared set and a missing-item convention; do not quietly discard lost candidates. Kendall tau-b handles ties in scientific tiers, while exported deterministic ordinal ranks are separately reported.

### Refresh and attribution

Compare saved runs for source drift, new/removed candidates, extraction/policy changes and source access changes. To attribute effects to a prompt, hold source bundle, candidate universe and policy fixed. To test a source update, hold extraction/policy where feasible. Multi-factor changes are reported as such.

## 10. Baselines and ablations

| Condition | What it tests |
|---|---|
| Structured-source only | Added value of unstructured evidence discovery |
| Fixed-query retrieval + same extractor | Added value of agent query choice, rather than simply LLM extraction |
| Bounded adaptive retrieval + extractor | Full research workflow |
| Without semantic verifier | Whether verifier improves support precision and at what recall/cost loss |
| Without disease relationship checks | Effect of related-condition evidence leakage |
| Without deduplication | Bias from repeated publications/cohorts |
| Baseline tiers versus declared alternative ranker | Utility of the ranking policy, separately from retrieval |
| Prediction omitted versus included in alternative policy | Incremental value of predictions under class-appropriate calibration |

Hold candidate budgets, accessible corpus, information cutoff and evaluation labels comparable. Report cost, tool/LLM requests, tokens, wall time, rate-limit incidents, human review time, abstention and coverage. More retrieval budget is not proof of better agent reasoning.

## 11. Leakage and circularity controls

- If ClinVar labels define a gold pathogenicity task, using those same labels as ranking inputs is a lookup baseline, not independent prediction validation.
- Existing chosen panel membership reflects historical selection and feasibility; it is not exhaustive truth about useful variants.
- Gene-level information shared by all alleles is not allele-specific discrimination.
- Do not tune weights/thresholds on held-out scPRIME outcomes and then report that same set as prospective validation.
- Publication/cohort duplicates must not straddle split boundaries.
- Retrieval date differs from publication/data availability date. A claimed historical/prospective benchmark must exclude post-cutoff sources and labels.
- Models may have pretraining exposure; record this limitation. A time cutoff on tool sources alone does not prove no model knowledge leakage. Require source-grounded claims nonetheless.
- Do not send withheld experimental labels to scouts or the extractor.

## 12. scPRIME-specific retrospective evaluation

Separate editing efficiency, target genotype certainty, assay sensitivity and biological phenotype. A negative transcriptomic result can reflect insufficient edited cells, low expression/detection, batch, differentiation state or limited gene panel coverage.

Predeclare effect definitions, quality filters, cell/biological replicate handling, multiple-testing policy and variant-level summary. Do not select genes or regulators based on the same outcomes later used to test prediction. Mechanistic evidence and disease association are evaluated independently from perturbation magnitude. Controls can reveal technical bias; never assume a synonymous edit is biologically inert.

Use existing scPRIME data to ask whether evidence features predict an interpretable perturbation after editing/QC, with confidence intervals and sample-size limitations. This cannot establish disease causality by itself. Future active learning requires a separate prospective design/evaluation.

## 13. Gates and proposed numerical targets

### Engineering release gates: required

- All INV-01 through INV-16 adversarial fixtures pass.
- No invalid identity, unsupported citation, forbidden tool action, or accepted rejected-claim feature in fixtures.
- Exact offline replay equality and checkpoint resume equivalence.
- Schema/export referential integrity and deterministic tier policy tests.
- Missing source, truncation, unresearched candidate and unjudged label are visible.

### Scientific development targets: proposals, not achieved results

Candidate recall >=0.90 within the declared gold universe; citation support precision >=0.95 on judged accepted claims; exact-allele functional precision >=0.98; zero adjudicated wrong-allele accepted links in the reviewed top set. These numbers need owner review and confidence intervals; they are not guarantees and must not be used to declare a tiny benchmark conclusive.

Report point estimates and uncertainty at an appropriate independent unit. Bootstrap by query/study/cohort for correlated claims; use binomial intervals only where independence is defensible. Zero observed errors does not mean zero population risk. Ranking utility has no invented universal pass threshold: agree on the objective and held-out comparator before making a scientific claim.

### Publication/generalization gate

Require B2 held-out evaluation, frozen protocol, source/label audit, uncertainty, matched baselines and error analysis. Disclose unavailable sources, model pinning limitations and human review workload. Engineering completion alone is not this gate.

## 14. Evaluation artifact contract

`eval_manifest.json`, metric tables with numerators/denominators and confidence intervals, identity/citation/directness mismatch JSONL, coverage/truncation tables, run-to-run comparison, cost/runtime table, blinded review/adjudication record and `evaluation_report.md`.

Each result cites run/content digests, gold/split versions, label purpose, objective profile, cutoff and source access. Report `not_evaluated` when gold labels are absent. The builder must never fill expected-performance placeholders with invented benchmark outcomes.

# tessera work log

Running log of what was built, what was found, and what is still open. Newest entries at the
bottom. Design decisions are recorded in [variant-agent-spec-v2/DECISIONS.md](variant-agent-spec-v2/DECISIONS.md)
(DEC-xx); this log records the work and the evidence behind it.

Conventions: dates are absolute; costs are LLM spend on the configured gateway; "tests" is the
offline pytest suite (`./.venv/bin/python -m pytest tests -q`).

---

## Baseline — 2025-09-25: Phase 0 (v1 spec)

- Weighted-composite scorer (D/F/M/E, anchors, novelty multiplier), 32 invariant tests,
  6-candidate synthetic demo. Built against the single-file v1 architecture
  (`../variant_prioritization_agent_architecture.md`).
- Known issue at the time: in the demo, D was constant across candidates (gene-level evidence),
  so its weight did no work.

## 2026-10-05 — Spec v2 review → 2.1.0-draft

- Reviewed the v2 package (SPEC, EVIDENCE_SCHEMA, ARCHITECTURE, EVAL, IMPLEMENTATION_PLAN).
  Wrote `REVIEW.md` (3 scientific + 6 engineering findings, no P0) and `DECISIONS.md`.
- Owner decisions: the product is **allele nomination** (new alleles for the next screen),
  transferable across genes (DEC-01); evolve the `tessera` package (DEC-09).
- Ranking replaced: `evidence_priority_v0.1` → `nomination_v0.1` — within (channel,
  consequence stratum), lexicographic (K mechanism concordance, O patient observation,
  P one predictor per stratum); characterized alleles → control channel; ties reported, no hash
  tie-break (DEC-01…07).
- Correction found while implementing: with stratified ranking, class-level concordance is
  constant inside a stratum → K uses sub-class features (PTV NMD status, missense domain
  membership) (DEC-14).
- M1 split into M1a (structured, no LLM) and M1b (literature/GEO) (DEC-10).

## 2026-10-05 — M0 offline vertical slice

- Phase 0 archived to `archive/phase0/`. New package: strict schemas (`Observed[T]`, identity,
  claims/links/decisions), reference bundle + left-align normalization, host mechanical
  validation, stratified ranking, content-addressed bundle, replay.
- Synthetic EVAL-B0 fixture (22 mentions incl. REF mismatch, protein-only, rsID-only, wrong
  assembly, multiallelic, non-left-aligned indel, gene claim posing as exact, intro-only mention,
  wrong offsets, human rejection, carrier comparison).
- 41 tests; replay byte-identical with sockets disabled.
- Bug caught by a test: snapshot raw bytes were not re-hashed on import.

## 2026-10-05 — M1a structured real-gene run (SETD1A / schizophrenia)

- Sources verified live: Ensembl REST 116 (MANE, RefSeq accession, sequence, VEP with
  AlphaMissense/SpliceAI/LOFTEE), gnomAD r4, SCHEMA browser REST (per-allele case/control allele
  counts), ClinVar E-utilities, UniProt, ClinGen dosage (DEC-15).
- Mechanism rule: ClinGen HI = 3 (disease MONDO:0005090) → loss_of_function (DEC-16).
- First real table: 1,889 ranked alleles. VEP AlphaMissense = local AlphaMissense table on all
  1,097 missense alleles.
- Findings: every observed SETD1A PTV is NMD-triggering → K degenerate in the PTV stratum
  (only 2 tiers); SCHEMA `in_analysis` meaning unverified, does not gate O (DEC-18).
- Robustness: optional-source failure (Ensembl xrefs HTTP 500) → partial run, not abort.

## 2026-10-05 — M1b layer 1: deterministic literature route

- LitVar2 gene variants resolved with Ensembl variant_recoder; unique → candidate, several →
  ambiguous, none → unmapped, frameshift protein shorthand → unsupported (DEC-21).
- Recoder pitfalls handled: batch results matched by `input` (not order); `c.4582-2delAG`
  silently mapped to a 1-bp deletion → stated-deletion length check.
- LitVar mis-assigned other genes' variants (e.g. SET/TAF1 p.V92A); rejected by the reference
  check.

## 2026-10-05 — LLM layer

- Provider-neutral client (`config/llm.yaml`, `openai_compatible` / `anthropic`), keys only via
  env / git-ignored `.env`, every call recorded with requested **and returned** model, replay
  without network, host budgets, one schema repair (DEC-20).
- Gateway pricing (checked via `/api/pricing`): several routes bill **per call** —
  `[AWS]claude-opus-4-8` $0.80/call, `[1000k按次计费]gemini-3.1-pro-preview` $0.05/call.
  Gemini returns `gemini-3.1-pro-low` (flagged as model mismatch).
- Incident: first screening run used the $0.80/call route; 17 calls succeeded, then the account
  ran out of balance (HTTP 403 pre-charge failure). Fixes: quota/auth errors stop the run
  immediately; default profile switched to Gemini.

## 2026-10-05 — M1b layers 2–3: screening and deep reading (DEC-22, DEC-23, DEC-24)

- Pool: Europe PMC 517 ∪ LitVar2 368 ∪ PubTator3 720 hits → **901 unique papers**.
- Only **79** name SETD1A in title/abstract → LLM screening only for those (2 Gemini calls);
  the rest get a free full-text scan (660 open-access papers).
- Rule fixes after first screening: common-SNP association / somatic cancer / engineered-mutant
  papers no longer P1; LitVar links must be confirmed in full text.
- Screening result: P1 24, P2 29, P3 138, P4 691, P5 19; 42 deep-readable; 11 P1 papers need
  manual access (e.g. Singh 2016, Takata 2014). Review sheet:
  `runs/setd1a/literature/review_screening.csv`.
- Cost controls: cross-run HTTP cache with per-source TTL; per-paper LLM result caches;
  `screen --refresh-pool` for incremental updates. Identical rerun cost: $0.
- Deep read (42 papers, 5 Gemini calls): supplements scanned deterministically; passages
  selected and packed several papers per call. 32 claims passed host checks, 27 rejected
  (20 non-verbatim quotes, 6 other-gene alleles, 1 allele not in passage).
- Safety rule added: a notation is attributed to the gene only if the nearest preceding gene
  symbol is the target — a PDF supplement listed SCN3A p.Val1532Leu beside SETD1A, and SETD1A
  residue 1532 is also Val, so the reference check alone would have accepted it.
- Literature notation normalizer (HTML entities, typographic characters, "/+", multiple alleles,
  "c.N-2delAG" shorthand expanded only within one exon/intron side).
- Verified against the paper (my initial doubt was wrong): PMID 40962831 and 41422157 used
  CRISPR isogenic knock-ins → c.4582-2_4582-1del and c.4596dup are `characterized_altered`
  → control channel.
- Panel coverage (protein-level, frameshift ±2 residues): 15/18 sites reached; 12 same allele,
  3 different allele at the same site (L93fs, E757fs, D424fs); not found: Y42fs, G705fs,
  E988E (designed control).
- Total LLM spend on Gemini for SETD1A literature: ~$0.45.

## 2026-10-06 — Provenance in the ranked table

- `ranked_variants.csv` gains `hgvs_c`, `hgvs_p`, `discovery_routes`, `source_refs` (ClinVar
  VCV, SCHEMA counts, PMIDs from LitVar/text/supplements) and `literature_evidence`; report rows
  show sources. Mentions and identity decisions are now part of the frozen evidence digest.
- SETD1A: 1,975 ranked alleles; 270 with ≥1 PMID; 9 with accepted literature evidence.
- Limitation noted: literature patient reports do not feed O (no counts extracted).

## 2026-10-06 — Status and open items

- Tests: 90 passing. Replay identical for the current SETD1A bundle.
- Open scientific decisions: SCI-09 (O before P?), PTV stratum has almost no ordering signal,
  low `case_controlled` threshold, literature patient reports vs O, human review of
  LLM-extracted claims, burden p threshold in `mechanism_rule_v0.1`.
- Open engineering: manual full-text input for inaccessible papers; GEO datasets; SQLite /
  checkpoints; panel step; Track A (disease → genes) not started.
- Next (agreed 2026-10-06): run Track B on a gene that has MAVE/DMS data so the ranking can be
  evaluated against measured effects.

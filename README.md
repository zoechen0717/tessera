# tessera

Allele nomination for variant-level functional screens (scPRIME and similar).
Given a disease, a risk gene and a screen context, tessera nominates **new**
alleles to install, as an auditable one-genomic-allele-per-row table. It is
meant to transfer across genes and diseases; SETD1A / schizophrenia is the
first development case, not a hard-coded assumption.

A *tessera* is one tile of a mosaic: each nominated allele is a tile, the
screen panel is the mosaic.

Specification: [`variant-agent-spec-v2/`](variant-agent-spec-v2/) — version
2.1.0-draft. Start with [REVIEW.md](variant-agent-spec-v2/REVIEW.md) and
[DECISIONS.md](variant-agent-spec-v2/DECISIONS.md); SPEC §8 defines the ranking.

## Status: M0 offline vertical slice (ISSUE-001–004)

| Milestone | Scope | State |
|---|---|---|
| **M0** | Contracts, genomic identity, frozen bundle, `nomination_v0.1`, exports, replay — synthetic data only | **done** |
| **M1a** | Real GRCh38/MANE, ClinVar, SCHEMA, gnomAD, VEP (consequence, AlphaMissense, SpliceAI, LOFTEE), UniProt, ClinGen; first real SETD1A table; no LLM | **done** |
| M1b | Literature / GEO workers, claim extraction, SQLite + checkpoints, AlphaGenome Atlas | — |
| Eval | MAVE/DMS cross-gene benchmark; SETD1A scPRIME case study | — |
| Panel | Per-stratum quotas, controls, editing constraints | — |

There is no LLM anywhere. Network access exists only in `tessera fetch`
(`src/tessera/sources/`); freeze, ranking and replay are offline.

## Ranking in one paragraph

Alleles are ranked only within `(channel, consequence stratum)` — PTV, splice,
missense, in-frame, synonymous, other — by `(K, O, P)`:
**K** mechanism concordance (the gene's mechanism, e.g. loss-of-function vs
dominant-negative, applied through sub-class features like NMD status and
domain membership; class-level inference, never allele evidence),
**O** patient observation tier, **P** the stratum's single declared predictor.
Alleles already characterized by an exact-allele experiment go to the control
channel. Ties are reported as ties, and keys that were constant in a stratum
are flagged as having done no work. All of this is **provisional policy**.

## Quick start

```bash
./.venv/bin/python -m pytest tests -q            # 51 tests
# real gene (network, ~1 min): fetch → run → replay
PYTHONPATH=src ./.venv/bin/python -m tessera.cli fetch --request requests/setd1a_schizophrenia.yaml --output runs/setd1a/input
PYTHONPATH=src ./.venv/bin/python -m tessera.cli run --input runs/setd1a/input/input.yaml --output runs/setd1a/run
PYTHONPATH=src ./.venv/bin/python -m tessera.cli replay --bundle runs/setd1a/run --output runs/setd1a/replay
./.venv/bin/python examples/synthetic_demo.py    # build fixture → run → replay
PYTHONPATH=src ./.venv/bin/python -m tessera.cli run --input <input.yaml> --output <dir>
PYTHONPATH=src ./.venv/bin/python -m tessera.cli replay --bundle <dir> --output <dir2>
```

Use `.venv` (native arm64 Python 3.13); the system `/usr/local/bin/python3.12`
is an x86_64 build under Rosetta.

## Layout

```
config/
  nomination_v0.1.yaml              key order, strata, one predictor per stratum
  mechanism_concordance_v0.1.yaml   mechanism × class / sub-class → concordance
src/tessera/
  schemas/      Observed[T], identity, sources/locators, claims/links/decisions,
                annotations + gene mechanism, features/ranked rows
  identity/     pinned reference bundle, left-align normalization, registry
  evidence/     host mechanical validation, effective-state resolution
  ranking/      policy loader, per-allele features, stratified ranking
  sources/      M1a adapters (Ensembl/VEP, gnomAD, SCHEMA, ClinVar, UniProt, ClinGen)
  m1a.py        fetch: sources → snapshots → input directory
  pipeline.py   freeze → derive; run; replay
  bundle.py     canonical JSONL bundle and digests
  reporting.py  ranked_variants.csv, report.md (deterministic)
  synthetic.py  synthetic EVAL-B0 fixture builder
  cli.py
requests/       per-gene M1a requests (setd1a_schizophrenia.yaml)
tests/          identity, evidence, ranking, replay, source parsers (real excerpts)
archive/phase0/ superseded weighted-composite scorer (not on the import path)
```

## What M0 enforces (tests)

| Spec ID | Check |
|---|---|
| INV-01, INV-15 | Equivalent indel representations merge; multiallelic alleles stay distinct; one row per allele |
| R-02 | REF mismatch, protein shorthand, bare rsID, wrong assembly are excluded with reasons, never dropped |
| INV-02 | A paper that mentions an allele outside the experiment cannot characterize it |
| INV-03 | Wrong text offsets or an altered snapshot reject the claim / fail the import |
| INV-04 | A gene-scoped claim posing as exact-allele evidence is rejected |
| INV-09 | Rejected and human-rejected evidence contributes nothing |
| INV-10 | Duplicated evidence does not change tiers |
| INV-12 | Synthetic and real sources cannot be mixed |
| INV-13, R-16 | Replay is byte-identical with networking disabled; a modified bundle fails |
| INV-17 | Mechanism never appears as allele evidence |
| INV-18, INV-19 | No cross-stratum rank; ties share a tier; degenerate keys reported |
| DEC-02 | Characterized alleles route to control (configurable) |
| DEC-03/14 | Switching LoF → dominant-negative reorders the PTV stratum; unknown mechanism makes K degenerate and the report says so |
| DEC-12 | Policy load fails if AVI and AlphaMissense are both keys |

# Implementation notes

Decisions made while building Phase 0 that the specification did not settle, plus
one thing the code found that the spec got wrong.

## 1. `Quantity` instead of `float | None`

Every numeric annotation is a `Quantity` carrying a `Missingness` state. Verbose,
and the reason is spec §42.1: a bare `float | None` cannot distinguish "not a
coding variant" from "table lookup failed" from "we never looked", and those three
must score differently. Making the states unrepresentable-as-numbers is what stops
`not_searched` from silently becoming 0.0.

Categorical inputs (`protein_domain_context`, `nmd_predicted`) carry only `raw`,
not `value`, because there is no meaningful float for "active_site" — the anchor
turns the category into a number, and storing a float would anchor twice. Hence
`Quantity.category()` and the `anchor_input` property.

## 2. Anchors receive raw values, not pre-scaled ones

Building the test fixtures immediately surfaced this: `_weighted` was passing
`q.value` to anchors, which breaks for every categorical anchor. The rule is now
that anchors always receive `q.anchor_input`. Worth recording because it is the
kind of mismatch that would otherwise produce plausible-looking wrong numbers
rather than an error.

## 3. `sparse` is separated from `degenerate` in the diagnostics

The first version of `diagnose()` flagged any component with zero variance as
degenerate and warned that it "contributed nothing to the ordering". The demo
showed that was an overclaim: a component present on only one candidate has zero
variance yet still moves that candidate's score relative to the rest. So:

- `degenerate` — present on ≥2 candidates and constant. Genuinely contributes
  nothing to the order.
- `sparse` — present on <2 candidates. Variance cannot be assessed; the asymmetry
  is a coverage gap, not a finding.

The diagnostic exists to prevent overclaiming, so it should not overclaim.

## 4. The LLM critic may downgrade but never promote

`Evidence.effective_level` takes the LLM critic's proposed level only if it is
*less* direct than the extracted level. Spec §14.2 says the critic is advisory;
this makes the asymmetry mechanical. The model may express doubt about evidence,
but the category an item belongs to is established deterministically.

## 5. Recurrence returns `None` rather than a number when it cannot be normalized

`_cohort_normalized_recurrence` returns `None` when no cohort size is available,
so the term is excluded and the remaining disease sub-weights renormalize. The
alternative — falling back to a raw proband count — is the gene-length artifact
spec §19.1.1 prohibits. The null model inside it is a placeholder (`cohort × 1e-4`);
Phase 2 replaces it with a per-gene, per-class mutation-rate model. The interface
is deliberately the *ratio*, not the count, so swapping the null model does not
require re-anchoring.

## 6. `_best_disease_match` defaults to `related`, not `exact`

Absent any evidence about how well the source condition matches the query
disease, the multiplier defaults to the `related` anchor (0.30). Assuming an exact
match in the absence of evidence is precisely the error spec §12.1.1 warns about,
and for SETD1A specifically it is the difference between a defensible ranking and
a confidently wrong one.

## 7. Environment

`/usr/local/bin/python3.12` on this machine is an x86_64 build under Rosetta;
polars detects the missing CPU features and warns that it will crash. The venv
uses `~/miniconda3/bin/python3.13`, which is native arm64.

## Not yet built, and deliberately so

- No network adapters. Phase 0 is offline by construction.
- No storage layer yet (`src/tessera/storage/` is an empty package). DuckDB +
  Parquet per spec §55; it is not needed until there are real annotations to persist.
- No panel selection yet (Phase 7), though `panel_eligible` and the blocking-flag
  machinery it depends on are in place and tested.
- No CLI yet. `pyproject.toml` declares the `tessera` entry point for Phase 3.

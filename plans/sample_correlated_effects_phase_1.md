# Correlated effects — Phase 1 results

**Spec:** `plans/sample_correlated_effects.md` (v3.1), §5.1 and §8 Phase 1.
**Branch:** `feat/genome-effects-v49`.
**Follows:** `plans/sample_correlated_effects_phase_0.md`.
**Status:** complete. Full suite green.
**Date:** 2026-09-20.

Phase 1 commits the schema. `ind_phenotype` carries the two residual-realization
columns the feature needs plus the two columns that were previously bolted on by
on-demand `ALTER TABLE`; `phenotype_meta` carries `condition_change_action` and
`define_phenotype()` sets it. **Nothing writes `residual_value` or
`residual_condition_level` yet** — both are `NULL` on every row until Phase 5,
and the new tests assert exactly that. There is no behaviour change for any
existing call.

---

## What shipped

| File | Change |
|---|---|
| `R/define_trait.R` | base `ind_phenotype` DDL gains `liability_value DOUBLE`, `cat_name VARCHAR`, `residual_value DOUBLE`, `residual_condition_level VARCHAR` |
| `R/open_pop.R` | `phenotype_meta` DDL gains `condition_change_action VARCHAR DEFAULT 'error'` |
| `R/add_phenotype.R` | the two `dbListFields()` + `ALTER TABLE ADD COLUMN` blocks for `liability_value` / `cat_name` are deleted; the `store_liability` and `cat_names` paths assign into `records` and the base columns receive them |
| `R/define_phenotype.R` | new argument `condition_change_action = c("error", "independent")`, `match.arg()`-validated, written to the `phenotype_meta` row; roxygen documents when it applies (D2) and that every phenotype in one residual block must agree (D6) |
| `R/sql_utils.R` | `TABLE_RESERVED_COLS$ind_phenotype` gains `residual_value`, `residual_condition_level`; `$phenotype_meta` gains `condition_change_action` |
| `R/schema.R` | `.sm_col()` descriptions for all five new columns, so `describe_table()` never prints `(no description)` for them |
| `tests/testthat/test-phenotype_schema.R` | **new** — 21 expectations, see below |
| `tests/testthat/test-schema-registries.R` | `liability_value` / `cat_name` leave `DEFERRED_COLS` — they exist from `open_pop()` now, so the registry test checks them like any other column |
| `CLAUDE.md` | `ind_phenotype` and `phenotype_meta` schema tables updated |
| `man/define_phenotype.Rd` | regenerated |

## Why the shape is what it is

**Base columns, not `ALTER TABLE`.** §4(a) of the plan: pre-1.0.0 there is no
migration path, so a column either exists from `CREATE TABLE` or it does not
exist. The old on-demand `ALTER` for `liability_value` / `cat_name` meant the
same table had a different column set depending on which phenotype types had
been recorded so far — a fact `mutate_table()`'s reserved-column check had to
anticipate (`DEFERRED_COLS` in the registry test) and `archive_replicate()` had
to tolerate. Now the column set is fixed at `open_pop()`.

**`residual_value` is liability-scale.** It is the `resid` term of
`liability <- pheno_mean + covariate_contrib + tbv + resid`, before any
count/categorical conversion — the only scale on which conditioning across
phenotypes is valid. `liability_value` (the whole liability) and
`residual_value` (its residual component) are different quantities; both are
kept.

**`residual_condition_level` stores the selected stratum, not the raw column
value** (Codex B6b). `NULL` means "the unconditional `R` was used", including
for an animal whose level matched no conditional stratum. That is what lets D2
be enforced: a later call can see exactly which `R` a realization was drawn
under.

**`condition_change_action` lives on `phenotype_meta`, per phenotype**, exactly
like `missing_component_action` — same table, same defaulting, same "configuration
belongs in a table so `restore_pop()` is complete" rule. D6 (agreement across a
block) is a Phase 2 / Phase 5 check; Phase 1 only stores the value.

**Appending with fewer columns than the table has.** The three
`dbWriteTable(..., append = TRUE)` calls in `add_phenotype()` write `records`
without the residual columns (and, for non-categorical phenotypes, without
`liability_value` / `cat_name`). DuckDB's R driver appends by column name, so
the omitted columns take `NULL`. This is the same mechanism that already let a
call without `...` follow a call with `...`; it is now also asserted directly
(the continuous-phenotype test checks all four are `NA`). Phase 5 replaces these
calls with register + `INSERT` anyway.

## Tests (`test-phenotype_schema.R`, 21 expectations)

1. A freshly opened population's `ind_phenotype` has exactly the nine base
   columns, in DDL order; `phenotype_meta` has `condition_change_action`; every
   one of them is in `TABLE_RESERVED_COLS` and has a non-empty
   `describe_table()` description.
2. A categorical phenotype with `store_liability = TRUE` and `cat_names`:
   `dbListFields()` is **identical** before and after `add_phenotype()` (no
   `ALTER TABLE`); `liability_value` is non-`NA` on every row; `cat_name` is one
   of the declared labels; liability and category agree (every "Dead" liability
   exceeds every "Alive" one); `residual_value` and `residual_condition_level`
   are `NA` on every row.
3. A continuous phenotype leaves `liability_value` and `cat_name` `NA`;
   `mutate_table(residual_value = 1)` is refused as reserved.
4. `define_phenotype()` stores `'error'` by default and `'independent'` on
   request, rejects any other value, and `overwrite = TRUE` replaces the stored
   value.

## Verification

- `test-phenotype_schema.R`: 21 passed.
- Full `tests/testthat` suite via `pkgload::load_all()` + `testthat::test_dir()`:
  **2932 passed, 0 failed, 0 errors, 1 skipped** (2911 before Phase 1 + the 21 new expectations; the one skip is pre-existing).

## Plan bookkeeping

- `plans/sample_correlated_effects.md`: header status → Phases 0–1 shipped;
  §5.1 retitled "shipped" and rewritten in the past tense with the
  `define_phenotype()` argument and the test file named; §8 Phase 1 marked
  shipped; `file:line` references into `add_phenotype.R`, `define_trait.R`,
  `open_pop.R`, `sql_utils.R` refreshed.
- `NEWS.md`: Phase 1 entry under **0.71.0 (in development)**.
- `CLAUDE.md`: both schema tables updated (this is Phase 8 work in the plan,
  pulled forward because CLAUDE.md must describe the schema that exists).

## Next

**Phase 2 — centralize covariance definitions.** `validate_phenotype_cov_block()`
implementing D1 (pair-row block discovery, `N == U`, per-stratum completeness,
one condition column per block, PSD) and the D3 realization lock
(`residual_value IS NOT NULL` — the column now exists); the three
`phenotype_var_comp` writers wrapped in transactions and calling it; D5 for the
two diagonal writers; distribution (blocks ≥ 2 only) and source checks; D6
agreement at definition time.

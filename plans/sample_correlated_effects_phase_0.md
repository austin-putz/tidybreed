# Correlated effects — Phase 0 results

**Spec:** `plans/sample_correlated_effects.md` (v3.1), §8 Phase 0.
**Branch:** `feat/genome-effects-v49`.
**Status:** complete. No behaviour change. Full suite green.
**Date:** 2026-09-20.

Phase 0 removes the legacy code left in the phenotype layer after the v0.63.x /
v0.64.0 migration cleanup: dead table-existence guards, duplicate `CREATE TABLE`
DDL, a dead fallback branch, a stale string, and an unused parameter. It writes
no new code and touches no schema. Its purpose is to leave Phases 1–7 editing
one copy of everything.

---

## What shipped

| File | Change |
|---|---|
| `R/add_phenotype.R` | the "backward-compat fallback" second `get_phenotype_var()` lookup in the independent residual branch is deleted (both lookups read the same unconditional diagonal of `phenotype_var_comp`); the `phenotype_components` existence guard around step 4 is deleted; the no-match warning names `phenotype_var_comp` instead of the long-gone `phenotype_residual_cov`; `get_residual_cov()` is called without the unused subset argument |
| `R/phenotype_helpers.R` | `get_residual_cov()` loses its `subset_df` parameter and its `phenotype_var_comp` existence guard |
| `R/define_effect_cov_matrix.R` | `ensure_trait_var_comp()` and `ensure_phenotype_var_comp()` are **deleted** with their five call sites; the `dbListTables()` guards in `get_trait_var()`, `get_phenotype_var()`, `load_trait_cov()`, `load_phenotype_cov()` are deleted |
| `R/define_residual_cov.R` | the `ensure_phenotype_var_comp()` call is deleted |
| `R/define_trait.R` | `ensure_trait_tables()`'s DDL list loses its `phenotype_meta`, `phenotype_components`, `phenotype_var_comp` entries — byte-identical copies of the DDL `open_pop()` runs a few lines before calling it |
| `R/effect_helpers.R` | `delete_existing_effect()` deletes stale `phenotype_random_effects` rows unconditionally |
| `man/` | `ensure_phenotype_var_comp.Rd`, `ensure_trait_var_comp.Rd` deleted; `get_residual_cov.Rd` regenerated |

9 files, +17 / −199.

## Why each removal is safe

Every guard removed asked "does table X exist?" for a table that `open_pop()`
creates unconditionally — `trait_var_comp`, `phenotype_meta`,
`phenotype_components`, `phenotype_var_comp` directly
([open_pop.R:279-344](../R/open_pop.R#L279-L344)), and `trait_meta`,
`phenotype_effects`, `phenotype_random_effects`, `ind_phenotype`, `ind_tbv`,
`ind_tgv`, `ind_ebv`, `ind_index`, `ind_true_index` through the
`ensure_trait_tables()` call at [open_pop.R:179](../R/open_pop.R#L179).
`restore_pop()` refuses a file without `ind_meta` and otherwise opens whatever
`open_pop()` wrote. There is no path to a `tidybreed_pop` whose connection lacks
these tables, so every guard had exactly one branch that could run.

The two `ensure_*_var_comp()` functions and the three duplicate entries in
`ensure_trait_tables()` were the drift risk the plan named: two copies of a
`CREATE TABLE` that only one of could ever execute. Phase 1 changes the
`ind_phenotype` and `phenotype_meta` DDL; after Phase 0 each has one home.

The `add_phenotype()` fallback was equivalent code: `residual_var_unconditional[t]`
is the diagonal of the unconditional residual rows for `t`
([phenotype_helpers.R:301-314](../R/phenotype_helpers.R#L301-L314)); the fallback
`get_phenotype_var(pop, "residual", t)` selects the same cell with the same
`condition_column IS NULL OR = ''` predicate. If the first is `NA`, so is the
second.

`pop$tables` is unaffected: `open_pop()` lists the four directly-created tables
in `tables_created` before `ensure_trait_tables()` appends its own names through
`unique()`, so removing later duplicates changes neither membership nor order
(which `schema()` reads).

## Scope note

The plan listed six stragglers. Five more of the same two kinds were found
while removing them and are included above: three further `dbListTables()`
guards in `define_effect_cov_matrix.R`, one in `add_phenotype()`, one in
`effect_helpers.R`, and the three duplicate DDL entries in
`ensure_trait_tables()`. Nothing outside the phenotype layer was touched.

**Left alone, on purpose.** `ensure_trait_tables()` is still called from
`define_trait()`, `define_phenotype()` and `define_residual_cov()` as well as
from `open_pop()`. Those three calls are no-ops for the same reason the guards
were, but the function is also the lazy DDL block Phase 1 edits, and its
`existing` check is one mechanism rather than a scattered set of guards.
Collapsing it to a single unconditional call from `open_pop()` is a reasonable
follow-up after this plan lands; it is not needed for it.

## Verification

- `test-define_phenotype.R`, `test-phenotype_composite.R`, `test-add_phenotype.R`,
  `test-schema-print.R`, `test-open_pop.R`, `test-restore_pop.R`,
  `test-define_trait.R`: all pass (253 expectations, 0 failures).
- Full `tests/testthat` suite via `pkgload::load_all()` + `testthat::test_dir()`: **2911 passed, 0 failed, 0 errors, 1 skipped** (the one skip is pre-existing and unrelated).
- `grep -rn "ensure_phenotype_var_comp\|ensure_trait_var_comp\|phenotype_residual_cov" R/ man/ tests/` → no hits; `subset_df` survives only as a local variable name inside `add_phenotype()`'s TBV section, unrelated to the removed parameter.

## Plan bookkeeping

- `plans/sample_correlated_effects.md`: header status → "implementation in
  progress, Phase 0 shipped"; §8 Phase 0 rewritten as a shipped table; every
  `file:line` reference into the five edited R files refreshed to the
  post-Phase-0 positions; the `phenotype_meta` DDL reference moved from
  `define_trait.R` to `open_pop.R`, which is now its only home.
- `NEWS.md`: entry under **0.71.0 (in development)**. `DESCRIPTION` is bumped
  in Phase 8 with the rest of the feature, matching how the genome-effects
  phases were versioned.

## Next

**Phase 1 — schema.** Base `ind_phenotype` gains `liability_value`, `cat_name`,
`residual_value`, `residual_condition_level`; the two on-demand `ALTER TABLE`
blocks in `add_phenotype()` go; `phenotype_meta` gains
`condition_change_action` and `define_phenotype()` the argument that sets it;
`TABLE_RESERVED_COLS` and the `_schema_meta` column descriptions updated.

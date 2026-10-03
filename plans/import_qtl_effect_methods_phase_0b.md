# Import QTL-effect methods — Step 0b results

**Spec:** `plans/import_qtl_effect_methods.md`, §10 Step 0b (B-1, B-2); found
while reviewing `plans/import_qtl_effect_methods_codex_review.md` (findings 2–3).
**Version:** 0.71.2.
**Status:** complete. Full suite: 66 files, 963 tests, 3505 passed, 0 failed, 0 errors, 1 skipped (run before the review-pass fixes). After those fixes, the eight phenotype and determinism test files were re-run: 278 expectations, 0 failed.
**Date:** 2026-10-02.

Two bugs in the code as it stood, unrelated to the import itself. They were fixed
ahead of step 1 so that each one ships with its own regression test. Steps 3 and 5
later retarget the same code from `ind_tbv` to `ind_tgv_total`, and they keep these
tests.

## What shipped

**B-1: group-mate sums accumulate exactly.**
- `.group_mate_tbv()` (`R/contributor_tbv.R`) summed mates' TBVs with a plain
  `SUM(t.tbv_value)`.
- It now uses `CAST(SUM(CAST(t.tbv_value AS GEV_ACC_TYPE)) AS DOUBLE)`, the
  evaluator's exact accumulator. Errors go through `.gev_accumulator_error()`, so an
  out-of-range sum gives the evaluator's message.
- `COUNT` and the R-side `total / n_mates` mean are unchanged.
- This covers every group path: `components` with `contributor_type = "group"`,
  and `formula_tbv`'s `group_sum()` / `group_mean()`. Both route through this one
  helper.

**B-2: no silent `Va = 0` under `prevalence`.**
- `.ap_liability_records()` looked up `get_trait_var(pop, "gen_add", t)` with the
  phenotype name and turned `NA` into 0.
- `define_phenotype()` now refuses `prevalence` together with `components` or
  `formula_tbv`. The message comes from `.prevalence_composite_msg()` and points to
  `thresholds =`.
- `add_phenotype()` gets a new PLAN-stage check, `.ap_check_prevalence()`
  (`R/add_phenotype_stages.R`). It runs before TBV materialisation, any draw, or any
  write.
  - It refuses a composite phenotype with a prevalence. This is a second line of
    defence, for a database edited by other means.
  - It refuses a simple trait with no stored `gen_add` diagonal. The message names
    `define_effect_cov_matrix()`, `define_trait(target_add_var = )` and
    `thresholds =`.
- The `NA → 0` fallback is gone. An `NA` there is now an internal error.

## What a user sees differently

- A composite categorical phenotype with `prevalence` errors in `define_phenotype()`.
  Before, it was accepted and its threshold ignored all genetic variance.
- A simple categorical trait with `prevalence` whose effects were written without a
  stored target (manual `effects`, or `scale_to_target = FALSE` with no target)
  errors in `add_phenotype()`. Nothing is written and no RNG is consumed.
- SGE phenotypes on large pens are now bit-identical across DuckDB thread counts.

## Tests

- **New: `tests/testthat/test-group-contributor-determinism.R`.**
  - `.group_mate_tbv()` sum and mean are `expect_identical()` at 1 and 8 threads.
    They match an R-side sum over the other pen members to 1e-12.
  - A seeded `formula_tbv` phenotype with `group_sum()` and `group_mean()` is
    `expect_identical()` between a 1-thread and an 8-thread run.
- **`test-add_phenotype.R`.**
  - `prevalence` with `formula_tbv` is refused, and nothing is written.
  - A hand-inserted `phenotype_components` row with `prevalence` is refused in PLAN,
    with no RNG use and no `ind_tbv` write.
  - `prevalence` with no stored target errors in `add_phenotype()`. Afterwards
    `.Random.seed`, `ind_phenotype` and `ind_tbv` are unchanged. A `user_values`
    call on the same trait still works.
  - The same trait with `thresholds` records phenotypes.
- **`test-phenotype_composite.R` §13.** The poultry SGE mortality test was the B-2
  case itself: a composite with `prevalence = 0.1`. It now asserts the refusal and
  runs with `thresholds = 1.3`.

**Checked against the old code** by reverting the `R/contributor_tbv.R` change:
- The SGE phenotype test **fails** on the old code.
- The direct helper comparison happened not to diverge on this fixture. It stays
  for its sum/mean semantics, and the file header says which test guards what.

## Review pass (2026-10-02)

A review of the finished diff found one regression and one coverage gap. Both are fixed.

- **Regression: `user_values` calls.** `.ap_check_prevalence()` ran for every call.
  But `add_phenotype(user_values = )` never evaluates the model and never places a
  threshold, so a prevalence trait with no target would have errored for nothing. The
  check now runs only when `user_values` is `NULL`. The no-target test adds a
  `user_values` call that succeeds.
- **Coverage: the composite backstop.** `define_phenotype(overwrite = TRUE)` deletes
  old `phenotype_components` rows. The only way to reach `.ap_check_prevalence()`'s
  composite branch is therefore a hand-written row. A new test inserts one and checks
  the refusal, that `.Random.seed` is unchanged, and that nothing is written to
  `ind_tbv`.
- **Checked, no change needed.**
  - Alignment of `has_components` / `has_formula_tbv` with `pheno_meta` rows: both
    are in `phenos` order before the topological sort.
  - `.gev_accumulator_error()` re-signals any non-accumulator error unchanged.
  - `define_trait_simple()` forwards `prevalence` but never `components`, so it is
    unaffected.
  - `test-phenotype_schema.R`'s prevalence test has a stored target and still passes.

## Deviations from the plan

- **Fixture size.** B-1's test needed pens of 200; the plan said "at least five
  mates". A probe found that the old `SUM()` never diverged at pens of 10, but
  diverged on 10 of 10 runs with 2 pens of 200 on a one-trait fixture. Whether
  DuckDB parallelises the aggregate depends on the join size and the query plan. A
  small fixture would have passed on the broken code.
- **Where the composite refusal lives.** It is in `define_phenotype()` as well as
  in `add_phenotype()`, so the error comes at definition time. The main plan's §6A
  and §10 record this.

## Plan bookkeeping

`plans/import_qtl_effect_methods.md`:
- the §10 table marks 0b done;
- the Step 0b section gains an "As built" paragraph;
- §6A's prevalence bullet notes that step 3 extends `.ap_check_prevalence()` to the
  active-block rule rather than adding a new check.

Also updated:
- `.claude/skills/tidybreed-api/SKILL.md` (`define_phenotype()` key arguments:
  `prevalence`, and the exact group sum);
- `man/define_phenotype.Rd`, plus two internal Rd files, via `devtools::document()`;
- `DESCRIPTION` and `NEWS.md`.

## Next

Step 1, the rename-only release (0.72.0).

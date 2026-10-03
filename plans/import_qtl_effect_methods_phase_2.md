# Import QTL-effect methods — Step 2 results

**Spec:** `plans/import_qtl_effect_methods.md`, §10 Step 2. Also §6B rule 5, §6C, §7.1–§7.6,
and gates A0–A22 in §11 (A12 was withdrawn by Q21).
**Version:** 0.73.0.
**Status:** complete.
- Full suite: 68 files, 1003 tests, 3692 expectations passed, 0 failed, 0 errors, 0 skipped (run with `NOT_CRAN=true`, so the dual-anchor tests ran); 127 warnings, 116 of them the new pool-expectation diagnostic (Q22).
- Before the step there were 66 files and 967 tests.

**Date:** 2026-10-03.

Part A brings in the exact multi-trait calibration and `anchor =`. It also brings in the
§6C rules for generation targets. Under Q21, `effects` and `scale_to_target` stay until
step 3; step 2 adds only the A19 refusal.

## What shipped

### The method: `R/qtl_congruence.R` (new, internal)

Ported from `simulate_qtl_effects/R/qtl_effects.R`, the current file (never the frozen
`qtl_effects_paper.R`).
- **Ported unchanged:** the validators, `.qtl_psd_eigen()`, `.qtl_relative_spectrum()`
  and `.qtl_spectrum_text()`.
- **`.qtl_congruence(B0, target, anchor)`:** the source's `.congruence` closure, lifted
  out as a function. Its algebra is unchanged.
- **A second rank error, checked first:** `rank(G) > rank(M)`, "the anchor cannot carry
  the target". After it comes the source's `rank(B0'MB0) < rank(G)`, "the drawn
  architecture cannot". The source had only the second (A5).
- **Anchor objects, so a generator never needs a genotype matrix for the genic case:**
  - `.qtl_anchor_diag(w)`;
  - `.qtl_anchor_design(Xc, n - 1)`: `O(nmk)`, never forms an m × m matrix, and gets its
    rank from the SVD of `Xc`.
- **`.qtl_sim_effects()`:** an internal mirror of the source's `sim_qtl_effects()`, with
  genic, realised, reference, `marginal_observed` and dual anchors. It has no samplers
  and no `marginal =` (Q5). It exists so the source suite ports nearly verbatim.

### `define_additive_effects()` (rewritten)

- **One flow for k = 1 and k ≥ 2:**
  1. Validate.
  2. Resolve the target (`.dae_resolve_target()`).
  3. Resolve loci, base and anchor data.
  4. Seed, then draw (`.draw_additive_architecture()`).
  5. Calibrate (`.qtl_congruence()`).
  6. Commit terms and target in one transaction.
  7. Run diagnostics.
  8. Print the messages.
- **New arguments:** `trait_var_comp_tbl`, `anchor = c("genic", "realised")` and
  `warn_bounds = c(0.8, 1.25)`. `G` may now be a single number for one trait.
- **Exact:** `B' M B = G` to ~1e-15 for `method = "shared"`, correlations included. In
  the smoke test, a target r_g of 0.4 was delivered as exactly 0.4.
- **`"union"`** scales each trait on its own (a 1 × 1 congruence). For a non-zero
  covariance it warns "approximate" and gives the delivered covariance and correlation.
- **`anchor = "realised"`:** `Cov(X)` of the base individuals.
  - `.dae_collect_dosages()` collects them as integer SQL sums, `id_ind` then `locus_id`
    ordered.
  - It refuses partial genotypes (where `extract_genotypes()` would read a missing copy
    as 0).
  - It refuses a base above `QTL_REALISED_MAX_CELLS = 2e7` cells.
  - It refuses pool and copy bases and scoped calls (A6).
- **Atomic:** `.ge_commit()` gained a `before_commit` hook. The passed `G` is written by
  `.tvc_write_block()` inside the term transaction (A10).
- **RNG:** `seed` is applied after every check. No refusal touches `.Random.seed` (A22).
- **Diagnostics (§7.4):**
  - pool expectation `2 Cov(H)` for a founder base;
  - observed `Cov(X)` for an individuals base;
  - the genic limit under `"realised"`.
  - A relative spectrum outside `warn_bounds` is a warning. Nothing is stored.
- **Messages:** the closing message says "exact" / "approximate" / "not calibrated" and
  gives the delivered covariance. A line-scoped call adds the §7.6 line-mean message
  (A13).

### Targets: `trait_var_comp` (§6C)

- **DDL:** `line_name VARCHAR` in the base `CREATE TABLE` (`R/open_pop.R`). Registered in
  `R/schema.R` and `R/sql_utils.R`.
- **One writer, `.tvc_write_block()`** (in `R/define_effect_cov_matrix.R`):
  - checks finite and PSD;
  - refuses on a whole-table, NULL-safe check of `effect_name` × any trait × `line_name`,
    even for an identical matrix. The error gives a working `remove_rows()` call (A15)
    for the connected block, `.tvc_block_traits()`;
  - inserts `%.17g` literals through `dbExecute()`, never `dbWriteTable()`;
  - opens no transaction of its own.
- **`define_effect_cov_matrix()`:**
  - gains `line_name` (genetic effects only);
  - accepts a number for one name;
  - checks dimnames and never relabels them (`.check_cov_dimnames()`);
  - wraps the genetic write in a transaction.
- **Readers:** `get_trait_var()` and `load_trait_cov()` gain `line_name = NULL`. A line
  uses its own rows, else falls back to the population-wide rows (per `effect_name`);
  lines are never mixed. The BLUPF90 parameter file and the prevalence threshold use the
  default, which is the population-wide rows.
- **Vocabulary:** `GENETIC_EFFECT_NAMES`, `GENETIC_EFFECT_NAMES_FUTURE` and
  `DERIVED_EFFECT_NAMES`, checked by `.check_effect_name_input()`.
  - `define_effect_cov_matrix()` refuses the future and derived names.
  - `define_effect_random()`, `define_effect_fixed_class()`, `define_effect_fixed_cov()`
    and `write_phenotype_cov_block()` refuse all three sets (A20).

### One entry path

- `define_trait()` loses `target_add_var` and `target_add_mean`, and `trait_meta` loses
  its column. `write_trait_var_diag()` is deleted.
- `define_trait_simple()` is deleted: R file, Rd, `NAMESPACE`, `_pkgdown.yml`, README,
  `package_summary` (md, R and html), the vignette table and the skill section.
- `restore_pop()` refuses a file whose `trait_var_comp` has no `line_name` (this always
  fires on an older file), or whose `trait_meta` still has `target_add_mean`.
- `rescale_effects_to_target()` is deleted (see Deviations).

## Files

- **New:**
  - `R/qtl_congruence.R`;
  - `tests/testthat/test-qtl-congruence.R`, `test-define_additive_effects-anchor.R`,
    `helper-qtl-effects.R`;
  - internal Rd pages for the new helpers and constants.
- **Deleted:** `R/define_trait_simple.R`, `man/define_trait_simple.Rd`,
  `man/write_trait_var_diag.Rd`.
- **Changed, R:** `define_additive_effects.R`, `define_effect_cov_matrix.R`,
  `define_trait.R`, `define_genome_effect_terms.R` (the `.ge_commit()` hook),
  `open_pop.R`, `restore_pop.R`, `schema.R`, `sql_utils.R`, `define_effect_random.R`,
  `define_effect_fixed_class.R`, `define_effect_fixed_cov.R`, `phenotype_cov_block.R`,
  `define_phenotype.R` (examples), `add_phenotype_stages.R` (error string),
  `chr_meta_helpers.R` (doc).
- **Changed, docs:** `CLAUDE.md`, both skills, `README.md`, `package_summary.md`,
  `dev/package_summary/*`, `dev/benchmarks/*` (3 files), the introduction vignette, the
  swine vignette script, `_pkgdown.yml`, `NEWS.md` and `DESCRIPTION`.
  - `CLAUDE.md`: the Two-Layer paragraph, and a new naming rule 6 (the reserved
    `effect_name` vocabulary).
  - Introduction vignette: the ADG + BF block is now written once, through
    `remove_rows()` + `G =`.

## Tests

- **Churn:**
  - 135 `define_trait(..., target_add_var = v)` call sites across 33 test files, plus
    `helper-parity.R`, were rewritten mechanically to the test-only helper
    `with_additive_target(pop, trait, var)` (`helper-pop.R`), which is `define_trait()`
    plus `define_effect_cov_matrix()`. Neither draws from the RNG, so the parity goldens
    are unaffected.
  - 29 calls at `G =` sites, in 5 files, were changed to plain `define_trait()`,
    because `G` now writes the target.
  - Rewritten tests:
    - `test-define_trait.R`: no target arguments or column; replacing a target goes
      through `remove_rows()`.
    - `test-define_additive_effects.R`: the two-`G` union/shared test.
    - `test-genome-effects-writer.R`: gate 44.
    - `test-mutate_table_defaults.R`: its `define_trait_simple()` call.
  - Deleted: the `define_trait_simple()` test in `test-phenotype_composite.R`.
- **New gate tests:** A0, A2, A3, A5, A6, A7, A8, A9, A10 (twice: a failure inside the
  transaction via a mocked `.tvc_write_block()`, and an infeasible rank), A11, A13–A22,
  and paper-12 through `add_offspring()`. The paper-12 TBV variance goes from 1.05 at
  generation 0 to 3.05 at generation 6, against a target of 4.
- **New internal tests:** paper-1 to 9 (paper-10 is inside paper-9), paper-11, A1
  (three seeds × G ∈ {1e-6, 4, 1e6}), A2 against a dense `C^-1/2 G^1/2` oracle (three
  seeds × m ∈ {2, 200} × full and rank-1), A4, A5, and the design anchor.
  - Dual tests are `skip_on_cran()`. They ran with `NOT_CRAN=true`.
- **Mutation spot-checks:**
  - Putting back the per-trait rescale for `"shared"` fails A2, A3, A5, A7, A10 and A18.
  - Moving `set.seed()` before validation fails A22.
  - Both were reverted and re-run green.
- **Vignette:** the introduction vignette, purled and sourced under `load_all()`, runs
  end to end. Its one warning is the pre-existing `farm` column notice.
- `devtools::document()` runs clean. `pkgdown::check_pkgdown()` finds no problems.

## Deviations from the plan

Recorded in the main plan's Step 2 "As built" paragraph.

1. **k = 1 uses the congruence too, and `rescale_effects_to_target()` is gone.**
   - One calibration path. A1 shows it equals the scalar rescale (~1e-16 in practice).
   - Behaviour change: a calibrated QTL set with no segregating locus is now an error,
     where it used to warn and write the draw unscaled.
   - The parity goldens compare at 1e-8 and still pass.
2. **`G` together with `trait_var_comp_tbl` is refused outright** ("not both").
3. **`anchor = "realised"` refuses manual `effects` and `scale_to_target = FALSE`.**
4. **`define_effect_cov_matrix()` accepts a number for one name.**
5. **Reserved names are refused for fixed effects** and in `write_phenotype_cov_block()`,
   not only for random effects.
6. **Diagnostics** skip an `ind_haplotype` base, and skip with a message above the cell
   limit or when the base genotypes cannot be collected. They run for the common scope
   only.
7. **Paper tests:**
   - 10 and 11 test the dual anchor, which is internal, so they stay on the internals
     (§1A table updated).
   - Only 12 goes through `add_offspring()`.
8. **`seed`** is validated as an integer scalar.

## Found while building

- **The default `warn_bounds` fires on most small founder pools.** It produced 116 of
  the suite's 127 warnings, on pools of 20–100 haplotypes. Sampling LD in `2 Cov(H)` is
  larger than the ±25% band. This is recorded as **Q22** (open, recommendation (b): a
  message for the pool comparison, a warning only for observed individuals). It is a
  user decision.
- `extract_genotypes()` reads a missing allele copy as dosage 0 and does not sort its
  rows. The realised anchor therefore collects dosages itself rather than through it.
- `trait_var_comp` has no uniqueness constraint. The block key is enforced in
  `.tvc_write_block()`; a `NULL` `line_name` would defeat a SQL `UNIQUE` anyway.

## Deferred

- **Q22:** the diagnostics default.
- **Step 3:** removing `effects` / `scale_to_target` (Q21), and the Q21 refusal in
  `define_effect_cov_matrix()`, scoped by line. Step 3 also rewords the "already stored"
  error so it says remove, then regenerate with `G =`.
- **Part B:** `extract_genetic_variance()`. Until then, the diagnostics are the only
  delivered-versus-target report.

## Plan bookkeeping

`plans/import_qtl_effect_methods.md`:
- §10 table: step 2 marked done.
- Step 2 gains an "As built" paragraph.
- §1A test table: as-built row added.
- §7.4: a pointer to Q22.
- New Q22.
- Review log entry 16.

## Next

Q22, then step 3 (consolidation, P2, Q18, and the Q21 removal), at 0.74.0.

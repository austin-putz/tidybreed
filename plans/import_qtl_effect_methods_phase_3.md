# Import QTL-effect methods — Step 3 results

**Spec:** `plans/import_qtl_effect_methods.md`, §10 Step 3; §6 (P1), §6A (P2), §6B
(value half), Q18, Q21, and gates T1–T9, PH1–PH8 (§11).
**Plan:** `plans/import_qtl_effect_methods_phase_3_plan.md` (written before the step).
**Shipped in three parts**, each its own commit and review (decided 2026-10-04):

| Part | Version | Content | Status |
|---|---|---|---|
| 3a | 0.74.0 | Consolidation (P1) + phenotypes read the total (P2) + value names + active-block prevalence rule | **done** 2026-10-04 |
| 3b | 0.74.1 | Q21: `effects` / `scale_to_target` removed, owner rule, `define_effect_cov_matrix()` refusal | **done** 2026-10-04 |
| 3b review | 0.74.2 | Fallback-line scope in the cov refusal; `G =` refuses/warns for stranded scopes | **done** 2026-10-04 |
| 3c | 0.74.3 | Q18: `formula_tbv` → `formula_tgv`, DSL `component =` / `table =`, Stage 1 ids out of SQL | **done** 2026-10-04 |
| 3c follow-ups | 0.74.4 | `mutate_derived()` ids out of SQL, constant warning dropped, case hints | **done** 2026-10-05 |
| Codex review | 0.74.5 | Nine findings verified and fixed; `remove_generated_effects()`; additive true-index warning | **done** 2026-10-05 |

---

## 3a — Consolidation and P2 (0.74.0)

**Status:** complete.
- Full suite (`NOT_CRAN=true`): 70 files, 1021 tests, 3817 expectations; 0 failed, 0 errors, 0 skipped. 11 warnings, all pre-existing (unchanged from 0.73.2).
- Before the step: 68 files, 1003 tests at 0.73.0 (0.73.2 added the R1–R9 gates).
- **Date:** 2026-10-04.

`ind_tgv` is now the only table of true genetic values. Every phenotype path
reads the **total** genetic value, so dominance, A×A and indicator terms reach
phenotypes. Under an additive-only model nothing changes: the total is the
breeding value, bit for bit.

### What shipped

**`ind_tbv` and `add_tbv()` are removed.**
- The breeding value is `ind_tgv` with `component_name = "additive"`.
- `R/add_tbv.R`, its Rd pages, the `NAMESPACE` export and the `_pkgdown.yml` entry
  are deleted, and so is the `CREATE TABLE ind_tbv`.
- `ind_tbv` is gone from every registry: reserved columns, primary keys, row
  keys, `SYSTEM_TABLES`, the `id_ind` table list, `schema.R` descriptions and
  table order, the `archive_replicate()` defaults and `print(pop)`.
- The evaluator helpers that existed only for it are deleted:
  `.gev_reserved_additive()`, the `tbv` branch of `.gev_require_terms()`, and
  `.gev_warn_tbv_stale()` with its warning.

**`add_tgv()` inherits from `add_tbv()`.** New signature:
`add_tgv(tbl, trait_name, index_names, weight_type, component_name = "additive",
overwrite_index, ...)`.
- `weight_type` is the old `type`, renamed for naming rule 1.
- `component_name` chooses which value a true index weights: `"additive"`
  (default), another component, or `"total"`. It is validated against
  `TGV_COMPONENT_NAMES` + `"total"`.
- `...` gives scalar custom columns on `ind_tgv`.
- Unknown indices are refused before anything is written.
- True-index rows are written in one transaction. Ids are registered as views,
  never pasted into SQL; the old `add_tbv()` pasted them.
- **Re-evaluation upserts.** Components the model no longer produces are deleted
  and the rest are updated in place on `(id_ind, trait_name, component_name)`, so
  custom columns survive. The old delete-and-reinsert would have wiped them on
  every `add_phenotype()` call.

**`ind_true_index` gains `component_name`.** Its row key is now
`(id_ind, index_name, weight_type, component_name)`, so an index on the
breeding values and one on the total can coexist.

**Value components renamed (§6B):**
- `order1_additive` → `additive`;
- `order1_dominance` → `dominance`;
- `order1_other` → `indicator`;
- `interaction` is unchanged.

The closed set is `TGV_COMPONENT_NAMES`, and `.gev_component()` asserts it.

**One reader of genetic values.**
- `.tgv_read(conn, ids, traits, components)` (`R/add_tgv.R`) serves phenotypes
  and true indices. It joins the selected individuals first and aggregates
  after, so a call never aggregates the whole table.
- `.tgv_by_id()` and `.group_mate_tgv()` (`R/contributor_tbv.R`) replace
  `.tbv_by_id()` and `.group_mate_tbv()`.
- `components` is `"total"` or a listed set. A listed component the individual
  has no row for contributes 0. An individual with no row at all is missing, and
  `missing_component_action` applies.

**Phenotypes read the total (P2).**
- Every path reads it: simple phenotypes, every composite contributor (self,
  dam, sire, group) and the `formula_tbv` DSL.
- Stage 1 calls `add_tgv()`.
- `phenotype_components.component_names`:
  - it defaults to `"total"` (in the DDL and in `define_phenotype()`) and is no
    longer reserved;
  - it accepts a comma-separated list of components;
  - `define_phenotype()` refuses unknown names, `"total"` mixed with others,
    and duplicates.

**`ind_tgv_total` is an ordered sum:** `list_sum(list(tgv_value ORDER BY
component_name))`. It is bit-identical at any thread count, and a one-component
total equals its row bit for bit. Listed components are summed the same way.
Group-mate sums keep `GEV_ACC_TYPE`.

**Simple-phenotype precheck.** The trait needs at least one term, of any kind
and any owner. A model written only with `define_genome_effect_terms()`
records phenotypes.

**Prevalence: the active-block rule** (`.ap_prevalence_genetic_var()`).
- The threshold uses the sum of the stored population-wide `additive`,
  `dominance` and `additive_by_additive` diagonals, each counted only if the
  model has terms of that kind.
- It errors, naming `thresholds =`, in two cases:
  - the model has an indicator surface or another kind of interaction;
  - a kind of term has no stored target.
- `.ap_check_prevalence()` (PLAN, before any write or draw) and
  `.ap_liability_records()` call the same helper.

**`restore_pop()` refuses pre-0.74 files.** It refuses:
- an `ind_tbv` table;
- an `ind_true_index` without `component_name`;
- stored `order1_*` names in `ind_tgv` or `phenotype_components`.

**`add_index()`** auto-detects `tgv_value` on `ind_tgv`. Its existing duplicate
error stops an unfiltered multi-component table from being summed.

### Files

- **Deleted:**
  - `R/add_tbv.R`;
  - `man/add_tbv.Rd`, `upsert_ind_tbv.Rd`, `dot-tbv_by_id.Rd`, `dot-group_mate_tbv.Rd`.
- **New tests:**
  - `test-tgv-consolidation.R` (T3, T4, T7–T9, custom-column survival, the
    reader-vs-view check, the unknown-index refusal);
  - `test-phenotype-total-genetic-value.R` (PH1, PH3–PH7).
- **Renamed tests (`git mv`):**
  - `test-add_tbv.R` → `test-add_tgv_breeding_value.R`;
  - `test-add_tbv_index.R` → `test-add_tgv_index.R`.
- **Changed, R:**
  - core: `add_tgv.R`, `genome_effects_eval.R`, `genome_effects_helpers.R`,
    `contributor_tbv.R`;
  - phenotypes: `add_phenotype.R`, `add_phenotype_stages.R`, `formula_helpers.R`,
    `define_phenotype.R`;
  - schema, registries and restore: `define_trait.R`, `open_pop.R`, `schema.R`,
    `sql_utils.R`, `restore_pop.R`;
  - other: `add_index.R`, `archive_replicate.R`, `tidybreed_pop.R`;
  - roxygen only: `remove_rows.R`, `define_index.R`, `add_dosage.R`,
    `define_additive_effects.R`.
- **Changed, docs:**
  - `CLAUDE.md`: the D7 sentence, "One evaluator, one table of genetic values",
    design principle 4, naming-rule examples, the Two-Layer bullet, Roadmap and
    coverage notes;
  - both skills;
  - `README.md`, the introduction vignette and the swine vignette script;
  - `dev/benchmarks/*` (3 files) and `dev/package_summary/render_package_summary.R`;
  - `NEWS.md` and `DESCRIPTION`;
  - regenerated `man/`.
- **Parity golden:** `tests/testthat/parity_golden/tbv.rds`. Only the column was
  renamed (`tbv_value` → `tgv_value`), and its values were kept. The new code
  matches it to 6.7e-16; the parity tolerance is 1e-8.

### Tests

**Mechanical churn:**
- 108 `add_tbv(` calls became `add_tgv(`.
- Reads of `ind_tbv` / `tbv_value` moved to `ind_tgv` filtered to `additive`,
  through the test helper `tgv_additive()` (`helper-pop.R`) or inline SQL.
- `type =` became `weight_type =`.

**Rewritten:**
- the old stale-warning tests were replaced by T5: a custom additive term is in
  `additive`;
- `add_index()` on `ind_tbv` became T6, which errors unfiltered and works
  filtered to `additive`;
- the D7 contract test now targets `add_tgv()`;
- the group determinism test reads both `"total"` and `"additive"`;
- the remove_rows, schema and open_pop table lists.

**New gates:**

| Gate | What it checks |
|---|---|
| T3 | one `additive` row, identical to the total |
| T4 | every component of a mixed model against an independent R oracle |
| T7 | additive and total true indices coexist, and the values match |
| T8 | archive, `remove_rows()`, and three `restore_pop()` refusals |
| T9 | schema |
| PH1 | simple phenotype is `mean + total + residual`, identical on an additive-only model |
| PH3 | components: the total, listed components, a missing component contributing 0, refusals, and atomicity under `overwrite` |
| PH4 | `mean` is an intercept |
| PH5 | total and phenotype identical at 1 and 8 threads, 400 animals with four components |
| PH6 | custom-owner-only trait records phenotypes |
| PH7 (3a part) | A, A+D and A+D+A×A sums; unused stored targets ignored; refusals leave tables and RNG untouched; an end-to-end threshold |

### Deviations from the plan

1. **`ind_tgv_total` uses an ordered floating sum, not `GEV_ACC_TYPE`.** The
   `DOUBLE → DECIMAL(38, 18) → DOUBLE` round trip moved about 10% of
   one-component totals by one ulp, which failed T3. The main plan's §6A was
   rewritten to match.
2. **`add_tgv()` re-evaluation upserts** instead of deleting and re-inserting, as
   described under "What shipped".
3. **`.tbv_by_id()` / `.group_mate_tbv()` were renamed in 3a**, because their
   bodies changed anyway. The other Q18 internal renames stay in 3c.
4. **The unknown-index check moved before the `ind_tgv` write.** `add_tbv()`
   wrote TBVs before failing on a missing index.

### Found while building

- **Bug in `define_phenotype()` (fixed).** It checked `components` *after*
  inserting the `phenotype_meta` row. Under `overwrite = TRUE` it checked
  after deleting the old definition too. A refused call therefore left a
  half-defined phenotype, or deleted the existing one. All component checks now
  run before any write, and PH3 pins this.
- **`restore_pop()` check fixed during review.** It matched `order1_` with SQL
  `LIKE`, where `_` is a wildcard. It now uses `starts_with()` / `contains()`.
- **Pre-existing, not fixed:** Stage 1 builds its contributor subset with
  `filter(id_ind %in% ids)`. That puts individual ids in SQL text, against
  CLAUDE.md's hard rule. Candidate for 3b.
- **Pre-existing, unchanged:** with `overwrite_index = FALSE`, individuals that
  already have a true-index row are skipped even when their genetic values were
  just recomputed. `add_tbv()` did the same.

### Verification

- Targeted test files pass after the review fixes.
- The introduction vignette, purled and sourced under `load_all()`, runs end to
  end. Its one warning is the old `farm` notice.
- The swine script parses, but cannot run here because it needs BLUPF90.
- `devtools::document()` and `pkgdown::check_pkgdown()` are clean.
- Grep for removed names (PH8 subset): the only hits are the `restore_pop()`
  refusal and the tests that pin it.
- Benchmark (`benchmark_phenotype_scale.R`, n = 250 / 1000, 0.73.2 vs 0.74.0 on the same machine): the PLAN stage is ~10–15% slower (e.g. correlated n = 1000: 0.176 → 0.200 s; named_effect n = 1000 call 2: 0.191 → 0.241 s), a fixed ~15–30 ms per call — ms/ind still falls with n, so no per-individual cost was added. It is the dearer `add_tgv()` write (upsert with NOT EXISTS) and the ordered-sum read. RESOLVE and COMMIT are unchanged.

### Plan bookkeeping

- `plans/import_qtl_effect_methods.md`:
  - §10 table: step 3 is split, and 3a is marked done.
  - Step 3 records the three 2026-10-04 decisions and gains the "As built, 3a"
    paragraph.
  - §6.3 task 3 now says `component_name`.
  - §6A's exact-sum paragraph is rewritten.
- `plans/import_qtl_effect_methods_phase_3_plan.md`: status line, and the
  3a.4 note on the view.

---

## 3b — Q21: generated means calibrated (0.74.1)

**Status:** complete.
- Full suite (`NOT_CRAN=true`): 70 files, 1024 tests, 3827 expectations; 0 failed, 0 errors, 0 skipped. 11 warnings, all pre-existing (unchanged from 3a).
- After 3a: 70 files, 1021 tests, 3817 expectations.
- **Date:** 2026-10-04.

Before 3b, `"generated"` meant "written by a generator", not "calibrated". Manual
`effects =` and `scale_to_target = FALSE` wrote generator-owned terms that matched
no target, and the prevalence threshold trusted the target anyway. Now the
generator always calibrates, the threshold accepts only generated terms, and a
target cannot be written under terms calibrated to a different one.

### What shipped

**`define_additive_effects()` loses `effects` and `scale_to_target`.**
- Every call samples and calibrates. Removed: the manual branch, the
  `need = "none" / "sigma"` modes of `.dae_resolve_target()`, the A19 refusal,
  the unscaled union warning, the "not calibrated" message and the matching
  branch of `.dae_check_realised()`.
- Known coefficients go through `define_genome_effect_terms()` under a user
  owner. The roxygen gains a *Generated means calibrated* section, and
  `add_tgv()` gains a writer example (line-specific known effects with
  `base_tbl` filling the centres).
- The sex-linked QTL error (`assert_qtl_autosomal()`) no longer says
  "pass `scale_to_target = FALSE`". It points at the writer.
- Seeded output of calibrated calls is unchanged: their code path is the same.

**The owner rule** (`.ap_prevalence_genetic_var()`).
- A trait with any term not owned by `"generated"` is refused for
  `prevalence =`, naming the owners and `thresholds =`.
- It runs before the kind check, in PLAN (no write, no draw) and in the
  liability stage, through the same helper.
- The `define_phenotype(prevalence = )` roxygen replaces the 0.72.3 caveat
  with the rule. `mean =` is documented as an intercept, with the recipe for a
  realised base mean (§6A "Mean").

**`define_effect_cov_matrix()` refuses a target under generated terms**
(`.tvc_refuse_under_generated()`).
- It refuses a genetic block when any trait of the block has `"generated"`
  terms of that kind at the block's scope: population-wide terms (common or
  parent-only) for `line_name = NULL`, and line-C terms for `"C"`.
- It is still refused after the old block is removed.
- The error gives the route: `remove_rows()` (when a block is stored), then
  `define_additive_effects(G = )`.
- A line's target before that line's effects is accepted (A17).
- Only the exported function checks. A generator's `G =` writes through
  `.tvc_write_block()` with its terms, as before.
- The term classification is shared with the threshold:
  - `.gev_target_kind()` maps a term to the `trait_var_comp` block that
    describes it;
  - `.gev_term_line()` gives its line scope.

**Reworded errors.**
- `.tvc_write_block()`'s "already stored" error gives the full sequence
  (remove, then re-run the generator).
- So does `define_additive_effects(G = )`'s version. The text before the
  `remove_rows()` call still ends "remove it first:", which A15 parses.

### Files

- **R:**
  - `define_additive_effects.R` (the removal; the roxygen);
  - `add_phenotype_stages.R` (the owner rule);
  - `define_effect_cov_matrix.R` (the refusal, its helpers, the reworded error);
  - `genome_effects_eval.R` (`.gev_target_kind()`, `.gev_term_line()`);
  - `chr_meta_helpers.R` (the error text);
  - `define_phenotype.R` and `add_tgv.R` (roxygen).
- **Tests:**
  - new helpers `with_additive_terms()` and `plant_generated_additive()`
    (`helper-pop.R`);
  - `additive_flat` (`helper-genome-effects-db.R`) now covers every owner's
    order-one additive terms and has an `effect_owner` column;
  - migrated: `test-add_tgv_breeding_value.R`, `test-genome-effects-eval.R`,
    `test-extract_allele_freq.R`, `test-add_tgv_index.R`,
    `test-define_additive_effects.R`, `test-define_additive_effects-anchor.R`,
    `test-genome-effects-writer.R`, `test-add_phenotype.R` and
    `test-parity.R` (a comment);
  - new gates in `test-phenotype-total-genetic-value.R`.
- **Docs:**
  - `CLAUDE.md` gains a "Generated means calibrated" hard rule;
  - both skills;
  - `README.md`, where the stale `base = "current_pop"` became `base_tbl`;
  - the swine script, `dev/benchmarks/benchmark_tgv_scale.R` and
    `benchmark_phenotype_scale.R`;
  - `NEWS.md`, `DESCRIPTION` and regenerated `man/`.

### Tests

**Call-site census.** At the start of 3b, 99 non-comment lines in 8 test files
matched `effects =` / `scale_to_target =`. That is about the plan's 93 calls:
some calls span two matching lines, and a few test names mention the argument.
They moved by category:

| Category | Files | Migration |
|---|---|---|
| (a) known values for a breeding-value oracle (40) | `add_tgv_breeding_value` 20, `genome-effects-eval` 17, `extract_allele_freq` 2, `add_tgv_index` 1 (Y loci) | `with_additive_terms()`, same arguments; no assertion changed |
| (d) filler while testing the generator (~38) | `define_additive_effects`, `genome-effects-writer` | the generator, calibrated and seeded; tests assert on scopes, centres and loci, not values |
| (b) union per-trait sets (3) | `define_additive_effects-anchor` A7, R5, R5b; one union fixture in `define_additive_effects` | `plant_generated_additive()` |
| (c) tests of the removed arguments | A19, A22's `G + effects` line, "accepts manual effects" | deleted; one "unused argument" test for both arguments |

**Notes on the migration:**
- The breeding-value oracle (`independent_tbv()`) and `additive_flat` read
  every owner now, matching `add_tgv()`.
- `make_two_line_pop()` fixes each line at `p = 0` or `1`, where no variance can
  be calibrated. Its line-scoped centring checks use the helper, which runs the
  generator's own `.dae_resolve_base()`. Its pooled checks (`p = 0.5`) still
  call the generator, including the Wahlund warnings.
- Gate 42 needs one segregating locus per one-locus call, so it now uses a
  30-locus fixture and picks loci with `0 < p < 1` in pool A.
- A15: the `define_effect_cov_matrix()` call after generation now hits the new
  refusal, not "already stored".
- PH7's kind refusals (indicator surface; dominance with no target) plant their
  terms under `"generated"`, so the owner rule does not mask them.
- `test-add_phenotype.R`'s "no stored target" case is now a real user path:
  generate with `G`, then `remove_rows()` the target.

**New gates (PH7, 3b part):**

| Gate | What it checks |
|---|---|
| Codex finding 1 | target 1, ten `ad_terms()` effects of 10, `prevalence = 0.1` → error; RNG, `ind_phenotype`, `phenotype_random_effects` and `ind_tgv` untouched; `thresholds =` works |
| owner rule | generated terms matching the target plus one tiny user term → error naming the owner |
| cov refusal | refused with a stored block (message names `remove_rows()` and `define_additive_effects(G = )`); refused after removal; another kind not blocked; the stated route succeeds and the threshold follows the new target |
| line scope | a line-A target before line-A effects is accepted; after line-A generation it is refused, even after removal; line B is still free |
| unused arguments | `effects =` and `scale_to_target =` are unused-argument errors, and nothing is written |

### Found while building

- **The store-then-select route was a dead end under existing terms.** When a
  `dominance` or `additive_by_additive` target is stored,
  `define_additive_effects(G = )` refused, and told the user to store `G` with
  `define_effect_cov_matrix()` and then select it with `trait_var_comp_tbl`.
  Once the trait has generated additive terms, that store is now refused, so the
  advice led nowhere. The error now checks the state. If generated additive
  terms exist, it says: remove the non-additive block, re-run with `G`, then
  store the block again. The cov-refusal gate pins this. Only the message
  changed; the 2026-10-04 "no new argument" decision stands.
- **Pre-existing, still open:** Stage 1 builds its contributor subset with
  `filter(id_ind %in% ids)`, which puts ids in SQL text. It was not touched in
  3b (it is unrelated to Q21). Scheduled for 3c (user decision, 2026-10-04;
  plan 3c.2b).

### Verification

- Every migrated file passes on its own after the change.
- **Mutation check:** with the owner check removed (in memory), the Codex
  finding-1 gate fails (5 failures). With the check restored, it passes.
- The introduction vignette, purled and sourced under `load_all()`, runs to the
  end; its one warning is the old `farm` notice. The swine script parses.
- `devtools::document()` and `pkgdown::check_pkgdown()` are clean.
- Both migrated benchmarks run: `benchmark_phenotype_scale.R --check` and
  `benchmark_tgv_scale.R` at n = 200.
- Grep: `scale_to_target` and `effects =` on the generator survive only in the
  unused-argument test, past NEWS, closed plans and `tools/quarto/legacy/`.

### Plan bookkeeping

- `plans/import_qtl_effect_methods.md`:
  - §7.1's signature is marked as built;
  - A19 is marked moot;
  - the §10 table marks 3b done;
  - the Q21 step-list bullet records the closed route;
  - an "As built, 3b" paragraph is added.
- `plans/import_qtl_effect_methods_phase_3_plan.md`: the status line.

### 3b review (0.74.2)

A full review of 3b after its commit (`1e00e39`) found two ways to leave
generated terms describing a target other than the stored one. Both were fixed,
and the generator-side behaviour was decided by the user.

1. **`define_effect_cov_matrix()` missed fallback line terms (bug).** Take
   line-A terms generated with no line-A target. They were calibrated to the
   population-wide target through the generator's `line → NULL` fallback. The
   refusal only looked at population-wide terms, so after `remove_rows()` a
   population-wide target of 50 was accepted for terms that deliver about 1,
   and the prevalence threshold trusted it. Reproduced, then fixed:
   - `.tvc_generated_terms()` counts, for a `line_name = NULL` block, the terms
     of every line with no block of its own.
   - This is sound because a line block cannot be added under existing line
     terms. A line block present now was therefore present when its terms were
     generated.
2. **`define_additive_effects(G = )` could strand other scopes (design,
   decided 2026-10-04).** A call replaces only its own scope. Take common
   terms with line-A fallback terms. Removing the target and re-running the
   common scope with `G = 50` succeeded, while the line-A terms still
   delivered 1. `.dae_target_dependents()` now acts before the seed:
   - It **refuses** when a fallback line depends on the target. The route is to
     give that line its own `G` first (`line_name = , G =`), then re-run.
   - It **warns** for terms of another `parent_origin` at the same target
     scope, naming them. Targets are per line, not per parent, so a refusal
     would deadlock: each scope would block the other, and generated terms
     cannot be removed any other way. The user chose this split over refusing
     everything with a new reset path.

**Gates** (`test-phenotype-total-genetic-value.R`):
- "line terms that fell back … block a new one";
- "G = refuses when a fallback line depends on the target, warns for other
  parents". It checks that nothing is written or drawn on the refusal, that the
  route succeeds, and that the warning's re-run is silent.

**Mutation check:** with the fallback-line clause removed (in memory), 7
expectations fail across both gates.

**Reviewed and found sound:**
- the generator diff, a pure removal (`exact` can never be NA);
- the owner rule's placement, before the kind check, in both PLAN and the
  liability stage;
- `.gev_target_kind()` / `.gev_term_line()`; `'any'` and `'unknown'` origins
  count as population-wide;
- the test helpers and the `additive_flat` change;
- the A15 parse anchor;
- docs, NEWS and the skills.

**Suite after the fixes** (`NOT_CRAN=true`): 70 files, 1026 tests, 3840
expectations; 0 failed, 0 errors, 0 skipped. 11 warnings, all pre-existing.
`pkgdown::check_pkgdown()` is clean.

**Version shift:** 3c becomes 0.74.3.

## 3c — Q18: `formula_tgv` and the DSL arguments (0.74.3)

Plan: `import_qtl_effect_methods_phase_3_plan.md` §3c, with 3c.2b added after
3a. **Breaking:** the argument and the `phenotype_meta` column are renamed, and
the DSL's positional `table` is removed. Seeded output is unchanged.

### What shipped

1. **Rename, no alias.** These are now `formula_tgv`:
   - `define_phenotype(formula_tbv = )`;
   - the `phenotype_meta.formula_tbv` column (DDL, `TABLE_RESERVED_COLS`,
     `schema()` text);
   - every internal name: `.validate_/.eval_formula_tgv()`,
     `.walk_formula_tgv_ast()`, `.build_tgv_env()`, `.ap_materialize_tgvs()`,
     `.assemble_composite_tgv()`, `tgv_kind`, `.FORMULA_TGV_DSL_FUNS`, and the
     plan entry's `tbv` field (now `tgv`, the total);
   - `R/contributor_tbv.R` became `R/contributor_tgv.R` (`git mv`), along with
     its Rd topic.

   `restore_pop()` refuses a file with the old column. "Composite TBV" in
   roxygen, messages and `schema()` text became "composite genetic value".
2. **The DSL** (`R/formula_helpers.R`, `.FORMULA_TGV_DSL_ARGS`):
   - `self` / `dam` / `sire` take one positional trait;
   - `group_sum` / `group_mean` take `trait, col`;
   - every call takes a named-only `component =`, which is one of
     `TGV_COMPONENT_NAMES` or `"total"` (the default);
   - the group calls take a named-only `table =`, default `"ind_meta"`.

   Trait, column and table may be symbols or strings. `col` and `table` must
   pass `validate_sql_identifier()`. Each reference reads its own component
   through `.tgv_by_id()` / `.group_mate_tgv()`.
3. **One pass.** `.walk_formula_tgv_ast()` validates each reference and
   replaces it with a `.tgv_<n>` placeholder in the same traversal, returning
   the substituted expression. `.substitute_tbv_ast()` is deleted.
4. **Define-time validation** (`.validate_formula_tgv(conn, formula)`). All of
   the following are checked before any write:
   - the parse, with exactly one expression;
   - the grammar;
   - the traits, with `agrep()` suggestions;
   - for each group reference, that the table exists, has `id_ind` and has the
     column.

   The old "validated at `add_phenotype()` time" message is gone.
5. **3c.2b, Stage 1 ids out of SQL.** `add_tgv()` is split into these parts:
   - the exported function, which resolves the individuals and writes the true
     index;
   - `.tgv_compute()`, which evaluates and writes `ind_tgv`;
   - `.tgv_compute_ids()`, used by Stage 1 for its contributor sets (dams,
     sires, group-mates). It registers the ids as a view, joins `ind_meta`, and
     drops `NA`, repeated and unknown ids, as `resolve_subset_ids()` does.

   Simple phenotypes still pass the user's own table to `add_tgv()`.

### Found while building

- **Arbitrary code ran from a stored formula (security bug).** The old walker
  recursed into any call, and the evaluator ran the expression in an
  environment whose parent is `baseenv()`. `formula_tbv = "system('...')"`
  therefore passed `define_phenotype()`, was stored, and executed at every
  `add_phenotype()`. A formula may now use only:
  - the DSL calls;
  - `+ - * / ^` and parentheses;
  - numbers;
  - the math whitelist.

  A string constant outside a call, a function name used as a trait, and a
  second expression (`"T; T"`) are refused too. Named arguments of the math
  functions (`round(x, digits = 2)`) still work.
- **The planned "duplicate-ref" bug was not live.** The old
  `.substitute_tbv_ast()` ignored `table` when matching. But the walk and the
  substitution visited references in the same depth-first order, so each node
  always took its own reference. Checked against 0.74.2 with formulas whose
  group terms differ only in `table`, in both orders. NEWS and the plan say so;
  the one-pass design now makes it hold by construction.
- **No quote can reach an id through the API.** `add_founders()` validates
  `line_name`, so ids are identifiers. The 3c.2b gate therefore passes ids with
  quotes straight to `.tgv_compute_ids()` (they are dropped as unknown) and
  compares 700 ids, plus `NA` and repeats, with an `ind_meta` filter, using
  `expect_identical()`.

### Deviations from the plan

- `.substitute_tgv_ast()` was dropped, not fixed (item 3), so the planned
  mutation "drop `table` from ref matching" has no code to mutate. The mutation
  run instead ignored a named `table =`. Three expectations failed: the
  table-reading gate, and the refusals of `table = "nope"` and of a table
  without `id_ind`.
- **Group column validation moved to define time (design consequence).** The
  plan asked for it, and it has a cost: a group column must now exist before
  `define_phenotype()`, so the config-first order no longer works for it.
  `formula =` (derived) keeps its warn-only check. `components =` groups are
  still checked only at `add_phenotype()`, as before; that is out of 3c's scope.

### Tests

`test-formula_tgv_dsl.R` (new; fixture: generated additive plus user dominance,
30 offspring, two pen groupings):

| Gate | What it checks |
|---|---|
| PH2 components | `T + dam(T)` reads totals; `dam(T, component = "additive")` reads the dam's additive value, which differs from her total; `sire(T, component = 'dominance')` |
| PH2 table | `group_sum(T, pen) + 2 * group_sum(T, pen, table = "pens")` matches a hand sum over each grouping (and the groupings differ); `group_mean(..., table = pens, component = "additive")` |
| PH2 refusals | 19 bad formulas: bogus or non-string component, extra positional, unknown or repeated named argument, nested call, positional `table`, too few arguments, non-identifier `col` / `table`, missing table, table without `id_ind`, missing column, `system()`, a string constant, a function name as trait, two expressions, unknown trait. `phenotype_meta` stays empty |
| 3c.2b | `.tgv_compute_ids()` on 700 ids (reversed, repeated, `NA`, two quoted unknowns) gives the same `ind_tgv` rows as `filter(id_ind %in% ids) |> add_tgv()` |

Also:
- `test-formula_phenotype.R`: fC5 now expects the define-time error and an
  empty `phenotype_meta`;
- `test-tgv-consolidation.R` T8: the `formula_tbv` refusal;
- the 36 `formula_tbv` uses in tests are renamed.

### Verification

- **Suite** (`NOT_CRAN=true`): 71 files, 1032 tests, 3891 expectations, after the review
  fixes; 0 failed, 0 errors, 0 skipped. There are 11 warnings, all
  pre-existing (the same count as 0.74.2).
- **Mutation check:** the walker was edited to ignore a named `table =`. Three
  expectations in `test-formula_tgv_dsl.R` failed. The file was restored and
  checked afterwards.
- The touched files pass on their own, and so does T8 with the new
  `formula_tbv` refusal.
- `devtools::document()` and `pkgdown::check_pkgdown()` are clean.
- The introduction vignette, purled and sourced under `load_all()`, runs to the
  end; its one warning is the old `farm` notice. The swine script parses (its
  one `formula_tgv` is `"WWD + dam(WWM)"`).
- **PH8 grep** over `R/`, `tests/`, `man/`, `vignettes/`, `dev/`, `NAMESPACE`,
  `CLAUDE.md`, `README.md`, `_pkgdown.yml`, `package_summary.md` and the
  skills, for `ind_tbv`, `add_tbv`, `tbv_value`, `id_tbv`, `formula_tbv`,
  `order1_`, `.gev_reserved_additive` and `.gev_warn_tbv_stale`. Hits remain
  only in the code and tests that refuse or check for these names:
  - `restore_pop()`'s pre-0.74.0 and pre-0.74.3 refusals;
  - T8's old-shape fixtures;
  - `test-open_pop.R`'s absence check;
  - PH3's `order1_additive` refusal.

### 3c review (before commit)

A full review of 3c found two bugs and two stale docs, all fixed in 0.74.3
before the commit.

1. **3c.2b was incomplete (bug).** `.ap_plan()` still read the phenotyped
   subset with `filter(id_ind %in% !!subset_ids)`, which renders the ids. It now
   uses `.ap_read_by_id()`.
   - The first 3c.2b gate called `.tgv_compute_ids()` directly, so it could not
     see this.
   - A new gate traces DuckDB's `dbSendQuery` method, which carries every
     statement, dbplyr's included. It runs `add_phenotype()` on a filtered
     subset, with simple, `components` (sire) and `formula_tgv` (dam,
     component, group) phenotypes, and asserts that no quoted id appears in any
     of the statements.
   - Mutation check: with the old `.ap_plan()` filter restored, the gate fails;
     with the fix, it passes.
2. **Derived `formula =` could run arbitrary R (security; predates 3c).**
   Confirmed end to end: `formula = "ADG + nchar(system(...))"` was accepted
   and the command ran during `add_phenotype()`.
   - `.check_derived_formula()` now gives the same closed grammar as the DSL
     (phenotype names, numbers, operators, the math whitelist).
   - It runs in `.validate_derived_formula()` and again in
     `.eval_derived_formula()`, which now evaluates with `enclos = baseenv()`.
   - The gate checks the refusals at define time, and that a malicious formula
     written straight into `phenotype_meta` is refused at `add_phenotype()`. It
     sets an environment variable as its payload, and the variable stays unset;
     no record is written.
3. **Stale docs.**
   - The `add_phenotype()` roxygen now says DSL references can read one
     component.
   - CLAUDE.md's "One evaluator" rule names `component =`, and a new hard rule
     says stored formulas have a closed grammar.

**Left for 0.74.4** (agreed with the user):
- `mutate_derived()` renders join ids into SQL;
- the scalar-constant warning fires on ordinary weights such as
  `0.5 * dam(WWM)`;
- the group-table existence check is case-sensitive, while DuckDB is not.

### Plan bookkeeping

- `plans/import_qtl_effect_methods.md`: Step 3 is marked done, the §10 table
  marks step 3 done, and an "As built, 3c" paragraph is added.
- `plans/import_qtl_effect_methods_phase_3_plan.md`: the status line, plus an
  as-built note on the duplicate-ref item.
- Skills: `formula_tgv` in the `phenotype_meta` table (schema), the DSL
  grammar under `define_phenotype()` (api), and `.tgv_compute_ids()` in the
  contributor note.

### 3c follow-ups (0.74.4)

These are the three items the 3c review left for after the commit (`67eb28f`).

1. **`mutate_derived()` ids out of SQL.** The `join_table` read and the
   destination-key read used `filter(.data[[join_by]] %in% !!ids)`.
   - Both now go through `.md_rows_by_key()`, which registers the key values as
     a view, joins on the quoted `join_by`, and drops `NA` keys, as SQL `IN`
     did.
   - `dbQuoteIdentifier()` is used rather than `validate_sql_identifier()`. The
     column already has to exist, and refusing SQL keywords would break a
     legitimate `join_by` such as `date`.
   - `record_sql()` / `leaked_ids()` moved to `helper-sql.R`, shared by both
     id gates.
   - Mutation check: the 0.74.3 `mutate_derived.R` fails the new gate.
2. **The scalar-constant warning is dropped.** Numbers are part of both
   grammars. The warning fired on every maternal weight, and was the source of
   five of the suite's eleven standing warnings.
3. **Case hint, not case folding.** The alternatives were:
   - match names case-insensitively, as DuckDB does, which would have meant
     changing every exact `%in%` check on `pop$tables` and `dbListFields()`;
   - or keep exact matching everywhere and explain the slip.

   The hint was chosen. `.case_hint()` adds "Did you mean 'pens'? Names are
   case-sensitive." in `.validate_formula_tgv()` (define time) and in
   `.read_one_per_id()` (the `add_phenotype()` lookup for `components` groups
   and covariates), so both places treat a case slip the same way.

**Suite** (`NOT_CRAN=true`): 71 files, 1033 tests, 3892 expectations; 0
failed, 0 errors, 0 skipped. There are 6 warnings, down from 11; the five
removed were the scalar-constant warning. The remaining six are pre-existing:
- `add_founders` named pool;
- `add_phenotype` on `ind_ebv`;
- `genome_map` BIGINT;
- parity;
- the composite covariate test.

`pkgdown::check_pkgdown()` is clean.

## Codex review of step 3 (0.74.5)

**Review:** `plans/import_qtl_effect_methods_phase_3_codex_review.md` (of 0.74.4,
`76aa0df`). **Date:** 2026-10-05.

Every one of the nine findings was reproduced through the public API before
it was fixed, with the review's own fixture (seed 3401, 12 loci, 30 + 30
founders). The measured numbers matched the review exactly: genic 100 vs
threshold variance 1 (F1); 6 surviving U terms delivering 1 under a target of
0 (F2); 2 vs 1 (F3); affected fraction 0.3442, liability variance 10.127 (F4);
60/60 non-finite records (F6); the `setNames()` length error (F7); 0 rows left
(F8); 33 ties, 0.833 vs 0.283 (F9); an all-zero true index (F5). **All nine
are agreed with; none was rejected.**

The common thread: Q21's "every term is `generated`" proves that each
*variant* was calibrated to *a* target. It does not prove that the stored
target the threshold reads is the one used (F1, F2), that the variants
together have that variance (F3), or that the genetic variance is the whole
liability's (F4). The threshold now checks each of those, and refuses rather
than approximate when a check cannot be made.

### Fixes

| # | Fix | Where |
|---|---|---|
| 1 | `trait_var_comp_tbl` must select the block the call's scope reads (`line → NULL`, `.tvc_resolve_line()`); it chooses effect and traits, never scope. Refused before the seed | `.dae_resolve_target()` |
| 2 | Every trait's old variant at the scope is deleted, including a zero-target `union` trait that gets no new terms. The message says so. A trait left with no terms takes the ordinary "no genome effects" error downstream | `define_additive_effects()` step 7 |
| 3 | The threshold refuses when one kind at one line scope has generated variants for more than one parent scope (reciprocal, or common + parent-only fallback). New `.gev_term_parent()` | `.ap_prevalence_genetic_var()` |
| 4 | `V` includes every named random effect's stored variance (`normal`, `uniform`). Refused: a `gamma` effect (mean `sqrt(v)`), a residual with conditional strata (the unconditional row is only a fallback), and `V = 0`. Checked in PLAN | `.ap_prevalence_env_var()`, `.ap_liability_records()` |
| 5 | *(user decision: warn)* An `"additive"` true index warns when an index trait has `indicator` terms or hand-written interactions. The index is still written; `"total"` is not offered as a breeding-value substitute | `.tgv_warn_structural_additive()` |
| 6 | A non-finite `formula_tgv` result (`Inf`, `NaN`, not an `NA` from a missing contributor) is an error in PLAN naming the individuals; no draw, no record | `.eval_formula_tgv()` |
| 7 | A constant `formula_tgv` is broadcast | `.eval_formula_tgv()` |
| 8 | `mean` (one finite number) and `thresholds` are validated before anything is written, and the replacement (delete, insert, residual, components) is one transaction. The residual goes through the transaction-free `.pvc_write_block()` | `define_phenotype()` |
| 9 | Strict exceedance: a liability on a cutpoint stays in the lower category (`findInterval(left.open = TRUE)`), matching "the fraction above" and PH7's oracle. `thresholds` must be finite and strictly ascending (`c(1, NA)` used to declare three categories and classify into two) | `liability_to_categorical()`, `define_phenotype()` |

**The parent-fallback warning advertised routes that are refused** (the
review's last qualification). *User decision:* a new exported
`remove_generated_effects(pop, trait_name, line_name = NULL, parent_origin =
NULL)` deletes the generated terms at exactly one scope, in one validated
transaction, and errors (deleting nothing, all-or-nothing across traits) when a
trait has none there. Targets and `ind_tgv` are untouched. The warning now
names it, and says prevalence is refused while both variants stand.

*Follow-up decision (2026-10-05, option A):* it removes **every kind** of
generated term at the scope, not only additive ones, so step 5's jointly
calibrated common-scope model is removed whole in one call, and never one
component at a time. Alternatives considered and rejected: additive-only with
a refusal for non-additive terms, and a per-kind argument (it would leave a
joint calibration half-removed). This was already the code's behaviour; the
roxygen, the error text and a gate (planted generated dominance and A×A at the
common scope removed with the common additive terms, a line-A variant
surviving) now say so.

**Scientific qualifications**, now in the `define_phenotype(prevalence = )`
roxygen and the helper's: the threshold is a Gaussian approximation at an
HWE/LE reference (a few large QTL miss the prevalence even there); summing
component targets assumes orthogonal components; fixed effects are not in `V`
(the prevalence is for records whose fixed effects are 0); exact fractions in a
known population come from `thresholds`. The review's other qualifications
(deterministic summation is not exact arithmetic; cached true indices are
caches) were already stated and needed no change.

### Seeded output

Unchanged for every model the suite exercises: ties occur only for discrete
liabilities with no residual, and the new refusals and variance terms apply
only to models the threshold used to get wrong. A categorical phenotype with
`prevalence` **and** a named random effect now gets a different (correct)
threshold.

### Tests

| File | Gates |
|---|---|
| `test-prevalence-threshold.R` (new) | F3 reciprocal (PLAN refusal: RNG and `ind_phenotype` untouched) and common + fallback, then the `remove_generated_effects()` route; F4 analytic cutoff `qnorm(0.9) * sqrt(10)` checked record by record, uniform counted, gamma and conditional strata refused; F9 ties at unit and end-to-end level; `V = 0` refused |
| `test-define_additive_effects.R` | F1: common call with line rows, line call with another line's rows, line call with population rows despite its own block (RNG and terms untouched), same scope accepted and calibrated, line → population fallback accepted. F2: the review's state, asserted on the stored terms and an independent `Σ n p q a²`, not the call's masked `delivered` |
| `test-formula_tgv_dsl.R` | F6: `T / 0`, overflow, `log` of a negative, and a categorical phenotype; F7: `"2"` and `"T + 2"` over 30 offspring |
| `test-define_phenotype.R` | F8: three invalid `mean`s and an injected failure after the delete (mocked `next_int_id()`), each leaving `phenotype_meta`, `phenotype_components` and `phenotype_var_comp` identical; F9 threshold validation |
| `test-add_tgv_index.R` | F5: the dosage surface warns and still writes 0; `"total"` and a generated-only trait do not warn |
| `test-remove_generated_effects.R` (new) | one scope removed and nothing else; empty scope and bad input refused with nothing deleted; user-owner terms untouched; every kind at the scope removed together |
| `test-tgv-consolidation.R` | T7's mixed fixture (hand-written indicator + interaction) now expects the F5 warning |

### Verification

- Every finding re-run through the review's probes after the fix: each now
  refuses, or gives the right number (F4: affected fraction 0.1013 for a
  requested 0.1 in 10,000).
- **Mutation checks:** with the env-variance term zeroed, the parent-scope
  check disabled, the inclusive classifier restored or the F5 warning
  removed, the matching gates fail (2/2, 2/2, 1/1, 1/1 tests).
- **Full suite** (`NOT_CRAN=true`): 0 failed, 0 errors, 0 skipped. The six
  standing warnings are unchanged; two new ones from T7 were the F5 warning
  on its mixed fixture, now asserted (and the warning no longer fires on a
  call that skips every individual). The two touched files were re-run after
  that change.
- `devtools::document()` and `pkgdown::check_pkgdown()` are clean.

### Plan bookkeeping

- `plans/import_qtl_effect_methods.md`: §10's table and the Step 3 heading;
  an "As built, step-3 review" paragraph; Q21 refined; requirements for steps 4
  and 5.
- Skills (`tidybreed-api`, `tidybreed-schema`), CLAUDE.md's "Generated means
  calibrated" rule, `NEWS.md`, `DESCRIPTION`, `_pkgdown.yml`, regenerated
  `man/`.

### Next

Step 3 is complete. The next step in `plans/import_qtl_effect_methods.md` §10
is step 4 (Part B, 0.75.0).

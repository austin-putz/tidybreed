# Import QTL-effect methods — Step 3 results

**Spec:** `plans/import_qtl_effect_methods.md`, §10 Step 3; §6 (P1), §6A (P2), §6B
(value half), Q18, Q21, and gates T1–T9, PH1–PH8 (§11).
**Plan:** `plans/import_qtl_effect_methods_phase_3_plan.md` (written before the step).
**Shipped in three parts**, each its own commit and review (decided 2026-10-04):

| Part | Version | Content | Status |
|---|---|---|---|
| 3a | 0.74.0 | Consolidation (P1) + phenotypes read the total (P2) + value names + active-block prevalence rule | **done** 2026-10-04 |
| 3b | 0.74.1 | Q21: `effects` / `scale_to_target` removed, owner rule, `define_effect_cov_matrix()` refusal | planned |
| 3c | 0.74.2 | Q18: `formula_tbv` → `formula_tgv`, DSL `component =` / `table =` | planned |

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

### Next

3b (0.74.1): remove `effects` / `scale_to_target` and migrate the 93 call sites;
add the owner rule; add the `define_effect_cov_matrix()` refusal; reword the
"already stored" error.

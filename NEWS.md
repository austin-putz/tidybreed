# tidybreed 0.36.1 (2026-06-06)

## Breaking changes

* `add_schema_description()` renamed to `define_schema_description()` to match
  the package naming convention (`define_*` = model configuration / metadata).
* Signature changed from `(pop, table_name, description, column_name, notes)` to
  `(tbl, column_name = NULL, description, notes = NULL)` where `tbl` is a
  `tidybreed_table` from `get_table()`, consistent with `define_chip()` and all
  other action functions. Column name is now the first positional argument so
  `define_schema_description("sex", "Sex of individual")` works without named args.

# tidybreed 0.36.0 (2026-06-06)

## New features

* **In-database schema documentation** — descriptions for all tables and
  columns are now stored in a `_schema_meta` system table inside the DuckDB
  file. Descriptions persist across `restore_pop()` restores and travel with
  the `.duckdb` file.
* `schema(pop)` — returns an aligned tibble of all user-visible tables with
  row counts, column counts, and descriptions.
* `describe_table(pop, "table_name")` — prints column names, types, and
  descriptions for any table; returns an invisible tibble.
* `add_schema_description(pop, table_name, description, column_name, notes)` —
  upserts a description for any table or column, including user-defined ones
  added via `mutate_table()` or `...` in `add_founders()`.
* `migrate_schema_meta(pop)` — adds `_schema_meta` to databases created before
  v0.36.0 without recreating the database.
* `print.tidybreed_pop()` now shows a hint to use `schema()` and
  `describe_table()`, and hides the system `_schema_meta` table from the table
  list.
* `summary.tidybreed_pop()` now pulls table descriptions from `_schema_meta`
  rather than a hard-coded internal vector and shows them below each table
  header. The `_schema_meta` table is excluded from the summary output.
* `define_chip()` now auto-registers a column description in `_schema_meta`
  for the chip membership flag it creates.
* `restore_pop()` now warns when opening a database that predates `_schema_meta`
  and suggests running `migrate_schema_meta()`.

## Internal

* `SYSTEM_TABLES` in `sql_utils.R` extended to include `_schema_meta`,
  `founder_haplotypes`, `phenotype_meta`, `phenotype_components`,
  `phenotype_residual_cov`, and `trait_random_effects`.

# tidybreed 0.35.1 (2026-06-06)

## Bug fixes

* `add_phenotype()` now validates `...` arguments against common misspellings
  of `phenotype_name` (e.g., `trait`, `phenotype`, `trait_name`, `pheno_name`).
  Previously, `phenotype_name = "ADG"` mistyped as `trait = "ADG"` would silently
  be treated as an extra column, leaving the selector NULL and adding **all**
  phenotypes instead of just ADG. Now raises a clear error directing users to the
  correct argument name.

# tidybreed 0.35.0 (2026-06-05)

## Breaking changes

* `mutate_table()` no longer accepts plain R vectors (length > 1) for setting
  per-row values. Positional vector assignment was unsafe: row order is not
  guaranteed and can silently drift when rows are added or reordered, causing
  wrong values to land on wrong rows with no error.

  Replace with an explicit tibble carrying the primary key column:
  ```r
  # Before (unsafe):
  mutate_table(gen = c(1, 2, 3))

  # After (safe):
  mutate_table(gen = tibble::tibble(
    id_ind = c("A_1", "A_2", "A_3"),
    gen = c(1L, 2L, 3L)
  ))
  ```

  Scalar broadcasting is unchanged and unaffected.

# tidybreed 0.34.2 (2026-06-05)

## Breaking changes

* `add_phenotype()` argument renamed from `trait_name` to `phenotype_name`.
  The old name contradicted the two-layer design (traits = genetic components,
  phenotypes = observed records) and caused silent bugs: passing
  `phenotype_name = "ADG"` fell into `...` as an extra column, leaving the
  selector `NULL` and generating records for *all* phenotypes. Update all
  named calls: `trait_name = "ADG"` → `phenotype_name = "ADG"`. Positional
  calls (`add_phenotype("ADG")`) are unaffected.

# tidybreed 0.34.1 (2026-06-03)

## Breaking changes

* `define_phenotype()` argument `trait_type` renamed to `type`. The old name
  conflicts with the genetic layer (`define_trait()` world); the observation
  layer should use the shorter, unambiguous name. Update all calls:
  `trait_type = "continuous"` → `type = "continuous"`.

* `phenotype_meta` DB column renamed from `trait_type` to `type` to match the
  argument name (per naming-consistency rule: argument names must match the
  column they populate). Existing databases must be recreated.

## Improvements

* `define_phenotype()` man page now fully documents `formula_tbv`, `formula`,
  and `type = "derived_formula"`, with examples for continuous, count,
  categorical, maternal composite (both `components` data frame and `formula_tbv`
  shorthand), SGE, and derived-formula patterns.

* `define_trait_simple()` argument `trait_type` likewise renamed to `type`.

* `sql_utils.R`: removed stale observation-layer columns (`trait_type`,
  `repeatable`, etc.) from `trait_meta` reserved cols (they moved to
  `phenotype_meta` in v0.31.0); added `phenotype_meta` reserved cols entry.

# tidybreed 0.34.0 (2026-05-24)

## New features

* `define_phenotype()` gains two new arguments for formula-based phenotype
  specification, reducing boilerplate for common composite and derived patterns:

  * **`formula_tbv`** — a DSL string that composes the true breeding value from
    named traits using a small set of functions:
    * Bare symbol or `self(trait)` — the focal individual's own TBV
    * `dam(trait)` — the dam's TBV (NA for founders; handled by `missing_component_action`)
    * `sire(trait)` — the sire's TBV
    * `group_sum(trait, col, table = "ind_meta")` — sum of group-mates' TBVs (self excluded)
    * `group_mean(trait, col, table = "ind_meta")` — mean of group-mates' TBVs
    * Standard R arithmetic operators and a whitelist of math functions (`sqrt`,
      `log`, `exp`, `abs`, `round`, etc.)

    Example — maternal weaning weight:
    ```r
    define_phenotype(pop, "WW", trait_type = "continuous",
                     mean = 230, residual_var = 180,
                     formula_tbv = "WWD + 0.5 * dam(WWM)")
    ```

    Example — SGE average daily gain:
    ```r
    define_phenotype(pop, "ADG_obs", trait_type = "continuous",
                     mean = 850, residual_var = 300,
                     formula_tbv = "ADG_direct + group_sum(ADG_SGE, pen_id)")
    ```

    The existing `components` data frame path is retained unchanged for advanced
    cases (covariate weights, Legendre polynomial weights).

  * **`formula`** — an arithmetic expression over previously simulated
    `ind_phenotype` records, used with `trait_type = "derived_formula"`. No
    residual variance is drawn; the value is computed directly. Topological sort
    in `add_phenotype()` ensures dependencies are evaluated first.

    Example — feed conversion ratio:
    ```r
    define_phenotype(pop, "FCR", trait_type = "derived_formula",
                     formula = "ADFI / ADG")
    ```

    Supports chains (`FCR_pct = "FCR * 100"`), scalar coefficients
    (`"ADFI - 0.036*ADG - 0.0072*MBW"`), and whitelisted math functions
    (`"ADG / MBW^0.75"`). Division by zero and `sqrt(negative)` produce NA
    with a warning.

* `add_phenotype()` now topologically sorts phenotypes when any
  `derived_formula` phenotype is present, guaranteeing dependency order even
  when `add_phenotype()` is called with no phenotype name.

* `phenotype_meta` gains two new `VARCHAR` columns: `formula_tbv` and
  `formula`. Existing databases are migrated automatically on first use via
  `ensure_trait_tables()`.

* `define_phenotype()` validates `formula_tbv` trait names against `trait_meta`
  at definition time, with close-match suggestions (using `agrep`) for
  unknown names. Group column existence is validated at `add_phenotype()` time
  with a clear error message and a hint to use the `table=` argument.

# tidybreed 0.33.0 (2026-05-23)

## Breaking changes

* `add_residual_cov()` renamed to `define_residual_cov()` — residual covariance
  is model configuration, not simulation output.
* `define_effect_int()` renamed to `define_effect_intercept()` — abbreviation
  violated the no-abbreviation naming rule.
* `compute_derived()` renamed to `mutate_derived()` — aligns with the `mutate_`
  prefix for column-update operations.
* `delete_rows()` renamed to `drop_rows()` — `delete_` is not a defined prefix
  in the naming convention.
* `add_trait_covariate()` removed (was deprecated since v0.6.0).

# tidybreed 0.32.0 (2026-05-20)

## New features

* `contributor_type = "group"` is now fully implemented in composite phenotypes.
  The SGE (Social Genetic Effects / Bijma) model is supported: each focal
  individual's social component is the sum (or mean, via `aggregation`) of
  pen-/cage-/litter-mates' social TBVs. Self is always excluded from the group
  aggregation. Singletons (no group-mates) receive a social contribution of 0
  and are not excluded. Relevant use-cases include pig ADG in nucleus pens,
  poultry cage mortality, litter effects, and any competition/cooperation model.

* `define_phenotype()` gains a `missing_component_action` argument
  (`"skip"` by default, or `"error"`). This is stored in `phenotype_meta` and
  applied uniformly across **all** missing composite components (group
  assignment, dam/sire TBV, etc.) on every `add_phenotype()` call.
  `"skip"` returns no phenotype record for the affected individual and emits a
  warning with the count and up to 5 example IDs. `"error"` stops with an
  informative message.

## Schema changes

* `phenotype_meta` gains a `missing_component_action VARCHAR DEFAULT 'skip'`
  column. Existing databases are migrated automatically by `ensure_trait_tables()`.

# tidybreed 0.31.0 (2026-05-19)

## Breaking changes

* `define_trait()` no longer accepts phenotype-level arguments (`trait_type`,
  `residual_var`, `expressed_sex`, `mean`, `min_value`, `max_value`,
  `prevalence`, `thresholds`, `cat_values`, `cat_names`, `store_liability`,
  `index_weight`, `economic_value`, `recorded_on`). Pass these to the new
  `define_phenotype()` function instead.

* `ind_phenotype.trait_name` is renamed to `phenotype_name`. SQL queries and
  direct table writes against `ind_phenotype` must use the new column name.
  Same rename applies to `trait_effects.phenotype_name` and
  `trait_random_effects.phenotype_name`.

* `add_tbv()` no longer applies `expressed_sex` filtering (that column was
  removed from `trait_meta`). Sex filtering is now handled by `add_phenotype()`
  via `phenotype_meta.expressed_sex`.

## New functions

* `define_phenotype()` — registers an observed phenotype in the new
  `phenotype_meta` table. Accepts all phenotype-level metadata (mean,
  `expressed_sex`, `trait_type`, `repeatable`, categorical thresholds,
  prevalence, etc.) and an optional `residual_var` that writes to the new
  `phenotype_residual_cov` table.

* `add_residual_cov()` — writes conditional or unconditional residual
  (co)variance entries to `phenotype_residual_cov`, enabling heterogeneous
  residuals by sex, herd, or any `ind_meta` group.

## New database tables

* `phenotype_meta` — one row per observed phenotype (decoupled from
  `trait_meta`).
* `phenotype_components` — assembly rules for composite phenotypes
  (maternal, social, reaction-norm). Supports `contributor_type` values
  `"self"`, `"dam"`, and `"sire"`.
* `phenotype_residual_cov` — residual (co)variance store; replaces the
  `"residual"` rows formerly in `trait_effect_cov` for phenotype residuals.

## Changes

* `define_effect_cov_matrix(pop, "residual", R)` now routes to
  `phenotype_residual_cov` via `add_residual_cov()`.
* `define_trait_simple()` chains through `define_phenotype()` automatically.
* `add_phenotype()` reads all phenotype metadata from `phenotype_meta` and
  assembles composite TBVs from `phenotype_components` when present.
* `add_index()` detects `phenotype_name` column in `ind_phenotype`
  automatically (no user action required).

# tidybreed 0.30.1 (2026-05-19)

## New

* Added `inst/CITATION` with two entries: an `@Manual` entry for the R package
  and an `@Misc` entry for the Quarto quick-start manual. Users can retrieve
  them with `citation("tidybreed")` or `toBibtex(citation("tidybreed"))`.
  A formal DOI will be assigned at v1.0.0.

# tidybreed 0.30.0 (2026-05-19)

## Breaking changes

* `set_qtl_effects_multi()` is removed from the public API. Pass a character
  vector to `define_additive_effects()` instead.

  ```r
  # old
  tbl |> set_qtl_effects_multi(trait_names = c("ADG", "BW"), G = G)

  # new
  tbl |> define_additive_effects(c("ADG", "BW"), G = G)
  ```

  A deprecated wrapper remains and will forward calls with a warning, but will
  be fully removed in a future release.

## Enhancements

* `define_additive_effects()` now accepts `trait_name` as a character vector
  (length >= 2) for correlated multi-trait effect sampling. New parameters:
  * `G` — additive-genetic (co)variance matrix; `NULL` reads from
    `trait_effect_cov`.
  * `method` — `"shared"` (default, all traits share the filtered QTL set) or
    `"union"` (per-trait QTL sets from existing `genome_effects` rows).
  Single-trait behaviour is unchanged.

# tidybreed 0.29.0 (2026-05-19)

## Breaking changes — rename `add_*` configuration functions to `define_*`

Eight functions renamed from `add_*` to `define_*` to match their
metadata/configuration semantics (pre-v1 API cleanup):

* `add_trait()` → `define_trait()`
* `add_additive_effects()` → `define_additive_effects()`
* `add_trait_simple()` → `define_trait_simple()`
* `add_effect_fixed_class()` → `define_effect_fixed_class()`
* `add_effect_fixed_cov()` → `define_effect_fixed_cov()`
* `add_effect_random()` → `define_effect_random()`
* `add_effect_cov_matrix()` → `define_effect_cov_matrix()`
* `add_effect_int()` → `define_effect_int()`

The rule: `add_*` now exclusively means inserting simulation output rows
(individuals, phenotypes, TBVs, EBVs, index values). `define_*` means writing
model configuration that shapes how the simulation runs (trait specs, effect
definitions, variance matrices).

# tidybreed 0.28.0 (2026-05-18)

## Enhancements

* `add_index()` is now fully generalized: any `tidybreed_table` with `id_ind`,
  `trait_name`, and a numeric value column can be passed as the first argument —
  not just `ind_ebv`. This enables phenotypic indexes (`ind_phenotype`), TBV-based
  indexes (`ind_tbv`), and user-defined tables.
  - New `value_col` argument (default `NULL`): auto-detected from the table name
    (`ind_ebv` → `"ebv_value"`, `ind_phenotype` → `"pheno_value"`,
    `ind_tbv` → `"tbv_value"`). Pass explicitly for unknown or user-defined tables.
  - A uniqueness check errors if any `(id_ind, trait_name)` combination has more
    than one row after filtering, preventing silent mis-computation from repeat
    phenotype records or multiple EBV evaluations.
  - A warning is issued when no filter is applied, reminding users to narrow to a
    single record per `(id_ind, trait_name)` before calling `add_index()`.

# tidybreed 0.27.1 (2026-05-18)

## Administrative

* License changed from MIT to **GPL-3**. Any modifications or derivative works
  distributed publicly must also be released under GPL-3, protecting the
  project's open-source contributions.

# tidybreed 0.27.0 (2026-05-16)

## New features

* New `define_table()` function. Creates a user-defined custom table inside
  the tidybreed DuckDB database and registers it in `pop$tables` so that
  `get_table()` and `mutate_table()` work on it immediately. Column names and
  types are declared with named typed-NA `...` arguments (same convention as
  `mutate_table()` schema pre-declaration). Optional `primary_key` argument
  names a column as the `PRIMARY KEY`. Use `DBI::dbAppendTable()` to insert
  rows after creation.

# tidybreed 0.26.0 (2026-05-16)

## New features

* `add_tbv()` gains three new arguments for computing true selection index
  values from TBVs:
  - `index_names`: character vector of named indices (from `define_index()`) for
    which a true index value should be computed. When `NULL` (default), only
    `ind_tbv` is written — existing behaviour is unchanged.
  - `type`: which weight column from `index_meta` to use — `"index"` (default,
    `index_weight`), `"economic"` (`economic_weight`), or `"both"` (writes two
    rows per individual distinguished by the `weight_type` column).
  - `overwrite_index`: when `FALSE` (default), individuals already present in
    `ind_true_index` for a given `(index_name, weight_type)` are skipped —
    avoids redundant recomputation across generations. Set to `TRUE` to
    recompute when index weights change.
* New table `ind_true_index` (`id_true_index`, `id_ind`, `index_name`,
  `weight_type`, `true_index_value`). Created automatically by
  `ensure_trait_tables()` on the next call that triggers it (e.g.,
  `add_trait()`, `initialize_genome()`). Existing databases are migrated
  transparently.

# tidybreed 0.25.0 (2026-05-16)

## New features

* `extract_genotypes()` gains an `effects_tbl` argument. Pass a filtered
  `tidybreed_table` from `get_table(pop, "genome_effects")` to extract
  genotypes at QTL loci (selected by trait, effect type, effect size, line,
  or any other filter). `chip_name` and `effects_tbl` may be used together;
  locus sets are unioned and deduplicated. `chip_name` now has a `NULL`
  default (backward compatible — positional calls like
  `extract_genotypes(tbl, "50k")` are unchanged).

# tidybreed 0.24.1 (2026-05-16)

## Naming consistency fixes

* `index_meta.economic_value` renamed to `index_meta.economic_weight` to match
  the `index_weight` column and the `{prefix}_weight` convention. Existing
  databases are migrated automatically.
* `define_index()` parameter `economic_values` renamed to `economic_wts` to
  match the `index_wts` parameter naming convention.
* `CLAUDE.md` updated: `index_meta` and `ind_index` table schemas documented;
  `define_index()` / `add_index()` added to the Implemented Functions section.

# tidybreed 0.24.0 (2026-05-16)

## New features

* `index_meta` gains an `economic_value DOUBLE` column. Existing databases are
  migrated automatically on the next call to `ensure_trait_tables()`.
* `add_trait()` now writes a global economic-value entry to `index_meta`
  (`index_name = NULL`) alongside its `trait_meta` row. When `overwrite = FALSE`
  (the default), an existing entry is preserved. When `overwrite = TRUE`, the
  `economic_value` is updated in place.
* `define_index()` gains two new parameters:
  - `economic_values` — optional numeric vector (same length as `trait_names`)
    written to the `index_meta.economic_value` column. Some values may be `0`.
  - `overwrite = FALSE` — re-calling for an existing `(index_name, trait_name)`
    pair is now a no-op by default. Set `overwrite = TRUE` to update weights
    and economic values in place (the previous unconditional-upsert behaviour).

# tidybreed 0.23.0 (2026-05-15)

## Breaking changes — genome_effects table replaces genome_meta effect columns

* New `genome_effects` table stores all QTL effect data in long format. The
  dynamically-added wide columns `add_{trait}`, `is_QTL_{trait}`, and
  `base_allele_freq_{trait}` no longer exist in `genome_meta`. Existing
  databases created with prior versions are incompatible.

* `define_qtl()` is removed. Its role is absorbed by `add_additive_effects()`,
  which now accepts a filtered `tidybreed_table` from `get_table("genome_meta")`
  as its first argument. QTL selection and effect assignment happen in one step.

  **Migration:**
  ```r
  # Old (v0.22.x):
  pop |> get_table("genome_meta") |> filter(...) |> define_qtl("T") |> add_additive_effects("T")

  # New (v0.23.0):
  pop |> get_table("genome_meta") |> filter(...) |> add_additive_effects("T")
  ```

* `add_additive_effects()` first argument is now always a `tidybreed_table`
  from `get_table("genome_meta")` (never a bare `tidybreed_pop`). The
  `base = "current_pop"` individual filter is now supplied via the new
  `base_tbl` argument instead of piping from an `ind_meta` table.

* `set_qtl_effects_multi()` first argument is now a `tidybreed_table` from
  `get_table("genome_meta")` (defining the shared/union loci set) instead of
  a `tidybreed_pop`. Add `line_name` and `base_tbl` arguments.

## Improvements

* `add_phenotype()` no longer duplicates TBV computation. It now calls
  `add_tbv()` internally and reads TBVs from `ind_tbv` for the phenotype model.
* `genome_effects` table supports a `line_name` column for future line-specific
  QTL effects (`NULL` = population-wide, the default).

# tidybreed 0.22.1 (2026-05-15)

## Breaking changes — rename `set_qtl_effects()` to `add_additive_effects()`

`set_qtl_effects()` has been renamed to `add_additive_effects()` to align with
the `add_` naming convention and to prepare for future `add_dominance_effects()`
and other effect-type functions. All callers (tests, vignettes, documentation)
have been updated.

**Migration:** replace every `set_qtl_effects(` with `add_additive_effects(`.
`set_qtl_effects_multi()` is unchanged.

# tidybreed 0.22.0 (2026-05-15)

## Breaking changes — eliminate implicit locus ordering

All positional locus selection methods have been removed. Loci are now always
selected via the `get_table("genome_meta") |> filter(...) |> fn()` pattern,
which is order-safe and consistent with individual selection throughout the
package.

### `define_chip()` — signature change

Old interface (removed):
```r
define_chip(pop, chip_name, n, method, locus_tf, locus_ids, locus_names)
```

New interface — pipe a filtered `genome_meta` table:
```r
pop |> get_table("genome_meta") |> filter(...) |> define_chip(chip_name)
```

The helper functions `select_by_n()`, `select_by_locus_ids()`, and
`select_by_locus_names()` in `chip_helpers.R` have been removed entirely.

### `define_qtl()` — signature change

Old interface (removed):
```r
define_qtl(pop, trait_name, n, method, locus_tf, locus_ids, locus_names)
```

New interface — pipe a filtered `genome_meta` table; `trait_name` is now optional:
```r
pop |> get_table("genome_meta") |> filter(...) |> define_qtl("ADG")
pop |> get_table("genome_meta") |> filter(...) |> define_qtl(c("ADG", "BW"))
pop |> get_table("genome_meta") |> filter(...) |> define_qtl()  # all traits
```

Shared QTL pleiotropy pattern:
```r
pop |> get_table("genome_meta") |> filter(is_QTL_ADG == TRUE) |> define_qtl("BW")
```

`define_qtl()` still returns `tidybreed_pop`, so it can pipe into
`set_qtl_effects()`.

### `add_trait_simple()` — `qtl_method` argument removed

The `qtl_method` parameter no longer exists. QTL are always placed randomly
inside `add_trait_simple()`. For non-random placement, use `define_qtl()`
directly after filtering `genome_meta`.

### `add_phenotype()` and `add_tbv()` — `trait_name` now optional

When omitted, all traits in `trait_meta` are used (in `id_trait` order).

```r
pop |> get_table("ind_meta") |> add_phenotype()  # all traits
pop |> get_table("ind_meta") |> add_tbv()         # all traits
```

# tidybreed 0.21.0 (2026-05-15)

## Breaking changes — naming consistency

All column and argument names have been standardised. Scripts that reference
any of the old names must be updated.

### Database column renames

| Table | Old column | New column |
|---|---|---|
| `ind_meta` | `line` | `line_name` |
| `ind_phenotype` | `value` | `pheno_value` |
| `ind_tbv` | `tbv` | `tbv_value` |
| `ind_ebv` | `ebv` | `ebv_value` |
| `trait_effect_cov` | `trait_1` | `trait_name_1` |
| `trait_effect_cov` | `trait_2` | `trait_name_2` |
| `trait_effect_cov` | `cov` | `cov_value` |
| `index_meta` | `index_wt` | `index_weight` |

### Public function parameter renames

* `add_phenotype(tbl, trait, ...)` — `trait` → `trait_name`
* `add_tbv(tbl, trait, ...)` — `trait` → `trait_name`
* `add_ebv(tbl, trait_name, ...)` — already correct; now consistent throughout
* Internal helpers (`build_data_file`, `parse_blupf90_solutions`,
  `update_covars_from_blupf90`, `compute_covariate_contribution`,
  `next_pheno_numbers`, `load_effect_cov`) — `trait` → `trait_name`

### CLAUDE.md — new naming consistency rules

Five rules now documented in `CLAUDE.md`:
1. Function argument names must match the database column they populate exactly.
2. All primary numeric value columns follow the `{prefix}_value` pattern.
3. All name/label columns end in `_name`.
4. All ID columns start with `id_`.
5. No abbreviations when the full word is unambiguous.

# tidybreed 0.20.1 (2026-05-14)

## Documentation & cleanup

* `NAMESPACE` now exports `compute_derived()` (was missing from v0.20.0 commit).
* `man/compute_derived.Rd` added (roxygen docs for `compute_derived()`).
* `man/add_trait.Rd` and `man/add_phenotype.Rd` updated to reflect v0.19.0
  binary → categorical unification (`cat_values`, `cat_names`, `store_liability`
  params; removed `"binary"` from `trait_type` choices).
* Added swine basic replication/generation YAML vignette script.

# tidybreed 0.20.0 (2026-05-14)

## New features

* `compute_derived()` — new action function that joins a filtered primary
  table with an optional secondary table, applies a user-supplied R function
  to produce a derived vector, and writes the result to one or more destination
  columns in any combination of tables. Reduces ~25-line dplyr workflows for
  computed date fields (e.g. `puberty_date = birth_date + AP_value`) to a
  single declarative pipe step. Accepts `write_to = c(table = "col")` syntax
  to write the same computed value under different column names in different
  tables.

# tidybreed 0.19.0 (2026-05-14)

## Breaking change

* `trait_type = "binary"` has been **removed**. Use `trait_type = "categorical"`
  with `prevalence = <p>` (for a 2-category trait) or `thresholds = c(<value>)`
  and add `cat_values = c(0, 1)` to recover 0/1 encoding. Binary was a strict
  subset of categorical — this eliminates the duplication.

## New features

* `add_trait()` gains three new parameters for categorical traits:
  * `cat_values` — numeric vector (length = n_categories) giving the phenotype
    value stored in `ind_phenotype` for each category (e.g. `c(0, 1)` for
    standard binary, or `c(1, 2, 3, 4)` for a 4-level score). Defaults to
    `1, 2, ..., K` when `NULL`.
  * `cat_names` — character vector of human-readable labels per category
    (e.g. `c("Alive", "Dead")`), written to the reserved `cat_name` column
    in `ind_phenotype`. Must not contain commas.
  * `store_liability` — logical. When `TRUE`, the raw liability value on the
    continuous scale is written to the reserved `liability_value` column in
    `ind_phenotype` alongside the observed category value.

* `add_trait()` now accepts `prevalence` for `categorical` traits (previously
  only accepted for `binary`). A single prevalence value defines a 2-category
  threshold; the actual threshold is computed at phenotype time from
  `target_add_mean + qnorm(1 - prevalence) * sqrt(VA + VR)`.

* `add_trait()` validates that exactly one of `thresholds` or `prevalence` is
  supplied for categorical traits; supplying both is an error.

* Two reserved column names are added to `ind_phenotype`:
  `liability_value` (DOUBLE) and `cat_name` (VARCHAR). Both are added
  dynamically via `ALTER TABLE` on first use, so existing databases are not
  affected until a categorical trait with these options is phenotyped.

* `trait_meta` gains three new columns: `cat_values VARCHAR`,
  `cat_names VARCHAR`, `store_liability BOOLEAN`. Existing databases are
  migrated automatically on the next call to `ensure_trait_tables()`.

# tidybreed 0.18.3 (2026-05-14)

## Bug fix / Enhancement

* `add_phenotype()` now enforces the `repeatable` field from `trait_meta`.
  If `repeatable = FALSE` and individuals in the candidate set already have a
  phenotype record for that trait, they are silently skipped and a warning is
  issued reporting the number rejected and the number that will receive a new
  record. Individuals with no prior record are phenotyped as usual.

# tidybreed 0.18.2 (2026-05-14)

## Bug fix

* `add_phenotype()` binary trait now uses a **fixed theoretical threshold**
  (`target_add_mean + qnorm(1 - prevalence) * sqrt(VA + VR)`) rather than the
  empirical quantile of the current batch. The old approach mechanically forced
  exactly the target prevalence; the new approach correctly shows sampling
  variation around the target, matching the standard threshold model from
  quantitative genetics.
* Updated the existing testthat binary-trait test to use explicit absolute
  tolerance (`abs(rate - 0.1) <= 0.04`) instead of relative tolerance.
* Added `tests/test_binary_mortality.R` — a stand-alone demo script that sets
  up 10% mortality with VA = 1, VR = 9 (h² = 0.10) and verifies the observed
  rate is within 1.5 pp of the target.

# tidybreed 0.18.1 (2026-05-14)

## Documentation

* Removed all references to removed functions (`mutate_ind_meta()`,
  `mutate_genome_meta()`, `set_residual_cov()`) from roxygen examples,
  `@seealso` sections, vignettes, and test scripts; replaced with current API.
* Replaced `[add_trait_covariate()]` in `@seealso` of `add_trait()` and
  `add_phenotype()` with the current `add_effect_*()` functions.
* Tightened `@param` type documentation across all exported functions to
  explicitly state scalar vs. vector, valid choices, and R type suffixes.
* Updated vignette subtitle to v0.18.1 and corrected table name `ind_phenotype`.

# tidybreed 0.18.0 (2026-05-14)

## New features

* `close_pop()` gains a `results_dir` argument. When provided, the `.duckdb`
  file is moved to that directory after the connection is closed. The directory
  is created automatically if it does not exist. An error is raised if a file
  with the same name already exists at the destination. In-memory databases
  are silently ignored.

# tidybreed 0.17.0 (2026-05-14)

## Breaking changes

* `pop$metadata` has been removed from the `tidybreed_pop` S3 object. Any
  code that reads `pop$metadata$n_loci`, `pop$metadata$n_chr`,
  `pop$metadata$chr_len_Mb`, `pop$metadata$chr_names`, or
  `pop$metadata$n_individuals` will need to query the database directly
  (e.g. `SELECT COUNT(*) FROM genome_meta`).

## Improvements

* **Database is the single source of truth.** Genome shape (`n_loci`,
  chromosome lengths, etc.) is no longer cached on the pop object at
  initialization. All functions that need this information query `genome_meta`
  live, so they automatically reflect any changes made to the table.

* **Dynamic genomes are now possible.** Because no genome summary is cached,
  users can add rows to `genome_meta` and columns to `genome_haplotype` /
  `genome_genotype` after initialization (e.g. for novel mutations, gene
  editing, or adding newly discovered QTL). `add_offspring()` and `add_tbv()`
  will automatically pick up the updated locus set without any metadata
  refresh step.

* `build_chr_info()` no longer takes a `chr_len_Mb` argument. Chromosome
  length for recombination is now derived as `MAX(pos_Mb)` per chromosome
  from the in-memory `genome_meta_df` slice — genomically equivalent since
  crossovers beyond the last locus position cannot affect observed alleles.

* `restore_pop()` simplified: no metadata reconstruction block. Function now
  just opens the connection, verifies `genome_meta` is non-empty, infers
  `pop_name`, and returns the pop object.

# tidybreed 0.16.0 (2026-05-14)

## New features

* `restore_pop(db_path)` — reconnects to an existing `.duckdb` file and
  reconstructs a fully operational `tidybreed_pop` object. Allows simulations
  to be resumed after `close_pop()`, a crash, or in a fresh R session without
  any stored metadata tables. All runtime metadata (`n_loci`, `n_chr`,
  `chr_len_Mb`, `chr_names`) is derived on the fly from the current state of
  `genome_meta`. An optional `pop_name` argument overrides the name inferred
  from the filename.

# tidybreed 0.15.2 (2026-05-13)

## New features

* `add_ebv()`: new optional `phenotype` argument accepts a `tidybreed_table`
  from `get_table("ind_phenotype") |> filter(...)`. When supplied, only
  phenotype records matching the filter are included in the BLUPF90 data file;
  the pedigree (and which animals receive EBVs) is unaffected. This solves the
  "future phenotype" problem where records are sampled early (e.g., number
  weaned, age at puberty) but should be excluded from evaluations until
  actually observed. Ignored with a warning in `parent_avg` mode.

# tidybreed 0.15.1 (2026-05-13)

## Bug fixes

* `print.tidybreed_table()`: fixed crash when `select()` is called before
  printing. Previously, the print method derived column names from
  `DBI::dbListFields()` (the full physical table schema) and tried to select
  them from `x$tbl` via `dplyr::all_of()`. After a `select()` call narrows
  the lazy tbl to a subset of columns, those original column names no longer
  exist in `x$tbl`, causing `dplyr::all_of()` to error. The fix uses
  `colnames(x$tbl)` instead, which always reflects the current visible columns.
  The header field count now also matches the narrowed column set.

# tidybreed 0.15.0 (2026-05-13)

## Breaking changes

* `id_ind` format changed from `{line_name}-{n}` (e.g. `Libra-1020`) to
  `{line_name}_{n}` (e.g. `Libra_1020`). Underscore separator avoids SQL
  query issues and DML-string construction with glue/paste. **Existing `.ddb`
  files are incompatible — recreate databases with this version.**

* `line_name` no longer accepts hyphens. Only letters, digits, and underscores
  are allowed (regex `^[a-zA-Z][a-zA-Z0-9_]*$`).

* `ind_phenotype`: `id_record VARCHAR` replaced by `id_phenotype INTEGER`.
  The new column is a global auto-incrementing integer primary key rather than
  a per-trait formatted string like `"ADG-529"`.

* New integer primary-key columns added to all major tables:

  | Table           | New column           |
  |-----------------|----------------------|
  | `ind_tbv`       | `id_tbv INTEGER`     |
  | `ind_ebv`       | `id_ebv INTEGER`     |
  | `ind_index`     | `id_index INTEGER`   |
  | `trait_meta`    | `id_trait INTEGER`   |
  | `trait_effect_cov` | `id_trait_effect_cov INTEGER` |
  | `index_meta`    | `id_index_name INTEGER` |

  Each column holds a globally unique auto-incrementing integer, making it
  easy to reference specific rows without composing multi-column keys. Former
  composite primary keys are now enforced at the application level (unique
  constraint retained where relevant).

* New internal helper `next_int_id(conn, table, id_col)` in `R/sql_utils.R`
  returns `MAX(id_col) + 1` for any table, used by all insert functions.

# tidybreed 0.14.0 (2026-05-13)

## New features

* New `delete_rows(tbl, tables = NULL, dry_run = FALSE, verbose = TRUE)` —
  delete rows from one or more population tables using the
  `get_table() |> filter() |> delete_rows()` pipeline.
  - `tables = NULL` (default): delete from the current table only, using the
    full composite row key (e.g. filters on `ind_tbv` by `trait_name` delete
    only that trait's rows, not all rows for those animals).
  - `tables = "all"`: delete matching individuals from every `ind_*` table
    plus `genome_haplotype` and `genome_genotype` that exist in the population
    (includes `ind_meta` for a complete "hard delete").
  - `tables = c(...)`: delete from an explicit subset of tables; each must
    have an `id_ind` column.
  - `dry_run = TRUE`: reports what would be deleted without modifying the DB.
  - Uses a temp-table JOIN internally (consistent with `mutate_table()`);
    safe for large ID sets.
* New internal constants in `R/sql_utils.R`:
  - `TABLE_ROW_KEYS`: composite row key columns per table, used for exact
    single-table deletion.
  - `IND_TABLE_ID_IND_COLS`: character vector of all tables with `id_ind`.

# tidybreed 0.13.2 (2026-05-12)

## Performance

* `add_offspring()` / `make_gamete()`: ~3x speedup for gamete simulation.
  - Replaced the per-locus `for` / `while` loop in `make_gamete()` with
    `findInterval()`, which counts crossovers before each locus position in
    vectorized C code. Haplotype assignment becomes a single arithmetic
    expression: `(current_hap - 1 + n_toggles %% 2) %% 2 + 1`.
  - Added `build_chr_info()` helper that pre-computes per-chromosome locus
    indices, positions, and lengths once per `add_offspring()` call.
    Previously, `genome_meta_df$chr == chr_id` (an O(n_loci) scan) was
    repeated for every chromosome × every gamete. The result is passed into
    each `make_gamete()` call, eliminating redundant work.

# tidybreed 0.13.1 (2026-05-11)

## Bug fixes

* `format_sql_value()`: added class-based guards for `Date` and `POSIXct`/`POSIXlt`
  values **before** the `db_type` string-matching branches. If a `Date` object
  loses its class attribute (e.g. via R's `for`-loop iterator stripping class on
  some R builds), `infer_duckdb_type()` could fall through to `"DOUBLE"`, causing
  `as.character(date)` = `"2026-02-26"` to be embedded **unquoted** in the UPDATE
  SQL. DuckDB then interpreted the bare expression as integer arithmetic
  (`2026 − 02 − 26 = 1998`) and threw
  `Conversion Error: Unimplemented type for cast (INTEGER -> DATE)`. The fix
  ensures Date and TIMESTAMP values always produce `DATE '...'` /
  `TIMESTAMP '...'` literals regardless of what `db_type` was inferred.

# tidybreed 0.13.0 (2026-05-10)

## Bug fixes

* `infer_duckdb_type()`: moved `inherits(value, "Date")` and
  `inherits(value, "POSIXct"/"POSIXlt")` checks **before** `is.numeric()`.
  `Date` objects are stored internally as doubles, so `is.numeric()` could
  return `TRUE` for them in some configurations, causing the type to be
  inferred as `DOUBLE` instead of `DATE`. The symptom was
  `SET mate_date = 2026-02-26` (unquoted) in the UPDATE SQL, which DuckDB
  interpreted as integer arithmetic and threw
  `Conversion Error: Unimplemented type for cast (INTEGER -> DATE)`.

* `format_sql_value()`: date and timestamp literals now use the explicit SQL
  type-keyword syntax (`DATE '2026-02-26'`, `TIMESTAMP '...'`) instead of bare
  quoted strings. This removes any remaining ambiguity about literal types
  regardless of column context.

# tidybreed 0.12.0 (2026-05-03)

## New features

* `summary.tidybreed_pop()` — detailed per-table summary with frequency tables
  for low-cardinality columns, 5-number summaries for numeric columns, date
  ranges for date/timestamp columns, and safe handling of wide genome tables
  (locus columns are not SUMMARIZE'd). Dispatches via `summary(pop)`.
* `print.tidybreed_summary()` — tidyverse-style console display with
  box-drawing separators and aligned column output. Accepts optional `tables`
  and `max_values` parameters to control scope and display thresholds.

# tidybreed 0.11.3 (2026-05-02)

## Bug fixes

* `add_index()`: replaced `stats::reshape()` + `sub()` rename with a direct
  matrix-based EBV pivot. `stats::reshape()` behaves differently on tibbles
  (returned by `dplyr::collect()` in the filter branch) vs plain data frames
  (returned by `DBI::dbGetQuery()` in the no-filter branch), causing column
  name mangling that made `wide[[tr]]` return `NULL` and threw
  `vapply … values must be length 1, but FUN(X[[1]]) result is length 0`.
  The fix builds an `n_ind × n_traits` matrix using direct row lookup and
  `match()`, then computes index values via `%*%`. No behaviour change for the
  no-filter path; the filter path now works correctly for any number of traits.

# tidybreed 0.11.2 (2026-05-02)

## Bug fixes

* `add_index()`: renamed `replace_index` parameter to `overwrite_index` to
  prevent R's partial argument matching from silently consuming user-defined
  extra column values named `rep`. Previously, passing `rep = 1L` as a custom
  column would be matched to `replace_index` (since `"rep"` is a unique prefix
  of `"replace_index"`), causing a `stopifnot(is.logical(...))` error. Now
  consistent with `overwrite_trait` in `add_ebv()`.

# tidybreed 0.11.1 (2026-05-02)

## Bug fixes

* `add_ebv()`: renamed `replace_trait` parameter to `overwrite_trait` to
  prevent R's partial argument matching from silently consuming user-defined
  extra column values named `rep`. Previously, passing `rep = 1L` as a custom
  column would be matched to `replace_trait` (since `"rep"` is a unique prefix
  of `"replace_trait"`), causing the value to never reach `...` and leaving the
  `rep` column at its declared default (`NA`).

# tidybreed 0.11.0 (2026-05-01)

## New features

* New `define_index(pop, index_name, trait_names, index_wts, ...)` registers a
  named selection index in the new `index_meta` table. Not all population traits
  need to be indexed — specify only the traits with non-zero economic weights.
  Supports extra user-defined columns via `...` (same pattern as other
  `add_*()` / `define_*()` functions). Re-calling with the same
  `(index_name, trait_name)` updates the weights in place (upsert).

* New `add_index(tbl, index_name, replace_index = FALSE, delete_all = FALSE, ...)`
  computes a selection index by multiplying EBVs by the weights defined in
  `define_index()` and appends results to the new `ind_index` table. Requires
  piping through `get_table("ind_ebv")` (the only `add_*()` function that takes
  `ind_ebv` as its starting table, not `ind_meta`). Each run increments an
  `index_number` column per individual, mirroring `ind_ebv`'s `eval_number`
  pattern. `replace_index = TRUE` clears prior runs for the named index;
  `delete_all = TRUE` clears the entire `ind_index` table. Issues a warning when
  no filter is applied (auto-selects the latest `eval_number` per individual per
  trait). Errors immediately if any individual is missing an EBV for any required
  index trait.

* Two new database tables created by `initialize_genome()`:
  - `index_meta (index_name, trait_name, index_wt)` — index definitions; supports
    user-defined extra columns via `define_index(...)`
  - `ind_index (id_ind, index_name, index_number, index_value)` — computed index
    values; supports user-defined extra columns via `add_index(...)`

# tidybreed 0.10.0 (2026-05-01)

## Breaking changes

* `ind_tbv` no longer has a `date_calc` column. The column was removed from the
  schema, the `add_tbv()` function signature (`date_calc` parameter dropped),
  and the internal `add_phenotype()` TBV write path. Users who need a date can
  add a custom column via `mutate_table()`.

* `ind_ebv` no longer has a `date_calc` column. The `add_ebv()` `date_calc`
  parameter has been removed. The `parse_blupf90_solutions()` internal function
  no longer accepts `date_calc`; its last argument is now `eval_nums` (named
  integer vector).

* `ind_ebv` primary key changed from `(id_ind, trait_name, model)` to
  `(id_ind, trait_name, model, eval_number)`. Existing databases are
  automatically migrated: `date_calc` is dropped and `eval_number INTEGER` is
  added (defaulting to 1 for any pre-existing rows).

## New features

* `ind_ebv` gains an `eval_number INTEGER` column. Each call to `add_ebv()` for
  a given trait increments the counter by 1 (global per `trait_name`, across
  all models). The latest evaluation can be found with
  `ORDER BY eval_number DESC LIMIT 1`.

* `add_ebv()` gains two new logical parameters:
  - `replace_trait = TRUE` — deletes all `ind_ebv` rows for the traits being
    added before inserting; new rows receive `eval_number = 1`.
  - `delete_all = TRUE` — deletes **all** rows in `ind_ebv` before inserting;
    new rows receive `eval_number = 1`. Takes precedence over `replace_trait`.

* `add_ebv(..., parent_avg = TRUE)` now queries parent EBVs using
  `MAX(eval_number)` per parent per `(trait_name, model)`. A warning is printed
  when any parent has more than one evaluation for a given trait, noting which
  `eval_number` is used.

# tidybreed 0.9.6 (2026-05-01)

## Bug fixes

* `upsert_ind_tbv()` / `upsert_ind_ebv()`: replaced the DELETE + `dbWriteTable`
  pattern with a DuckDB-native `INSERT … ON CONFLICT DO UPDATE SET` UPSERT.
  Previously, calling `add_tbv()` (or `add_phenotype()`) after a prior
  `add_tbv()` that had written user-defined extra columns (e.g. `rep`) would
  silently reset those columns to `NULL` in existing rows for other traits or
  on re-runs without the extra-column argument. The new UPSERT only updates
  columns that are explicitly present in the incoming data frame, leaving all
  other columns in existing rows untouched.

# tidybreed 0.9.5 (2026-05-01)

## Bug fixes

* `add_tbv()`: removed the `existing_ids` skip guard that prevented re-computing
  TBVs for individuals already in `ind_tbv`. The guard was redundant with the
  DELETE + INSERT logic inside `upsert_ind_tbv()` and blocked custom column
  updates (e.g. `rep`) when looping over replicates.

# tidybreed 0.9.4 (2026-04-29)

## New features

* `mutate_table()` gains `.set_default` parameter to create SQL DEFAULT
  constraints on new columns. When `.set_default = TRUE`, future INSERT
  operations from `add_founders()`, `add_phenotype()`, and other `add_*()`
  functions automatically use the default value when the column is not
  explicitly specified. This enables easy schema pre-declaration and
  consistent metadata across individuals.

# tidybreed 0.9.2 (2026-04-29)

## Documentation

* Added `vignettes/tidybreed-introduction.Rmd` — a placeholder introduction
  vignette directing users to the QMD-based quick-start in `vignettes/qmd/`
  until the package reaches v1.0.0.

# tidybreed 0.9.1 (2026-04-28)

## Internal

* `genome_haplotype` and `genome_genotype` locus columns now use `UTINYINT`
  instead of `INTEGER`, reducing per-locus memory from 4 bytes to 1 byte (4×
  reduction). `UTINYINT` range (0–255) covers biallelic haplotypes (0/1),
  diploid genotypes (0/1/2), and polyploid up to 8n or 16n.

# tidybreed 0.9.0 (2026-04-27)

## New Features

* **Eager table initialization**: `initialize_genome()` now creates all core
  database tables upfront (`ind_meta`, `ind_phenotype`, `ind_tbv`, `ind_ebv`,
  `trait_meta`, `trait_effects`, `trait_effect_cov`, `trait_random_effects`).
  Users can call `get_table()` and `mutate_table()` on any table immediately
  after `initialize_genome()`, before any data has been added.

* **Custom field forwarding in `add_*` functions**: `add_founders()`,
  `add_phenotype()`, `add_tbv()`, and `add_ebv()` now accept `...` for custom
  column values that are written atomically with the new rows. Column types are
  inferred from the R type (`0L` → INTEGER, `0` → DOUBLE, `"text"` → VARCHAR,
  `TRUE` → BOOLEAN). Scalars are broadcast; vectors must match the inserted
  row count. Reserved column names are blocked.

  ```r
  # Founders with custom columns in a single call
  pop <- pop |>
    add_founders(n_males = 10, n_females = 100, line_name = "A",
                 gen = 0L, farm = "Iowa")
  ```

* **`mutate_table()` on empty tables**: calling `mutate_table()` on an empty
  table now creates the column schema (via `ALTER TABLE ADD COLUMN`) instead
  of warning and returning early. This enables pre-declaring typed column
  schemas using typed NAs before any rows exist:

  ```r
  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(gen = NA_integer_, farm = NA_character_)
  ```

* **`prepare_extra_cols()` internal helper** (in `R/sql_utils.R`): shared
  validation + type-inference + ALTER TABLE logic used by all `add_*`
  functions with `...` support.

* **`TABLE_RESERVED_COLS`** and **`TABLE_PRIMARY_KEYS`** in `R/sql_utils.R`
  expanded to cover `ind_tbv`, `ind_ebv`, `trait_meta`, `trait_effects`, and
  `trait_effect_cov`.

---

# tidybreed 0.8.2 (2026-04-24)

## Bug Fixes

* **`add_ebv()` / `parse_blupf90_solutions()`**: fixed two bugs when reading
  back blupf90+ solutions with `OPTION origID` active:
  - The parser was looking for a file named `solutions`, but blupf90+ writes
    `solutions.orig` when `OPTION origID` is present. The file path is now
    correct.
  - `original_id` values that look numeric (e.g. `10`, `20`) were parsed as
    integers by `read.table`, causing the `%in% all_ped_ids` comparison against
    a character vector to silently return zero matches. The column is now
    coerced to character immediately after parsing.
* Added `tests/testthat/test-parse_blupf90_solutions.R` covering: correct
  parsing of the aligned/header format, numeric-looking IDs, missing-file
  error, and no-match warning.

---

# tidybreed 0.8.1 (2026-04-24)

## Bug Fixes

* **`add_ebv()` / `write_renum_par()`**: fixed three errors in the generated
  `renum.par` parameter file for BLUPF90:
  - The intercept (mu) effect was declared as `cross` (class variable) instead
    of `cov` (covariate). The mu column is always 1 for every row, so it must
    be `cov`.
  - The animal random effect was missing its own `EFFECT` block before `RANDOM
    animal`. Without it, `renumf90` was associating the last fixed-effect
    `EFFECT` block (e.g. sex at column 3) with the animal random effect,
    assigning the wrong column number instead of column 2 (`id_ind`).
  - As a consequence of the missing animal `EFFECT` block, the fixed-effect
    column for sex was not written correctly. With the animal `EFFECT` block
    restored at column 2, the sex `EFFECT` block at column 3 now lands
    correctly.

---

# tidybreed 0.8.0 (2026-04-24)

## New Features

* **`add_ebv()`** — new function to populate the `ind_ebv` table with estimated
  breeding values. Two modes:
  - `software = "blupf90"`: runs `renumf90` + `blupf90+` from PATH. Builds
    pedigree, data, and (optionally) genotype files automatically from the
    database. Parameter file (`renum.par`) is auto-generated from `trait_effects`
    and `trait_effect_cov`. All input/output files are written to a timestamped
    subfolder of `run_dir` for full reproducibility. Supports BLUP and ssGBLUP
    (when `chip_name` is supplied). After the run, solutions are parsed and EBVs
    are stored in `ind_ebv`.
  - `parent_avg = TRUE`: computes EBVs as the simple average of parent EBVs
    already in `ind_ebv` for the given `model`. Returns `NA` (with a warning)
    for animals whose parents lack EBVs.
* `add_ebv()` accepts a `tidybreed_table` (from `get_table()` + optional
  `filter()`) as its first argument, following the same calling convention as
  `add_phenotype()` and `add_tbv()`.

---

# tidybreed 0.7.2 (2026-04-24)

## Breaking Changes

* **`ind_phenotype` schema simplified**: removed `env`, `rep`, and
  `date_measured` columns. These are replaced by a single `pheno_number
  INTEGER` column that auto-increments per individual × trait (1 = first
  record, 2 = second, etc.). Users can add their own columns via
  `get_table("ind_phenotype") |> mutate_table(...)`.

* **`add_phenotype()` parameters removed**: `env`, `rep`, and `date_measured`
  arguments have been dropped from the function signature. Remove them from
  any existing calls.

* **`ind_phenotype` now registered for `mutate_table()`**: primary key
  (`id_record`) and reserved columns (`id_record`, `id_ind`, `trait_name`,
  `value`, `pheno_number`) are registered, enabling filtered column additions
  via `get_table("ind_phenotype") |> filter(...) |> mutate_table(my_col = ...)`.

---

# tidybreed 0.7.1 (2026-04-24)

## Bug Fixes

* `add_tbv()` now checks for existing TBV records before computing. Individuals
  that already have a TBV for a requested trait are skipped with an informative
  message rather than silently overwritten.

# tidybreed 0.7.0 (2026-04-24)

## Breaking Changes

* **`add_effect_cov_matrix(pop, effect_name, cov_matrix, trait_names, tol)`** —
  new unified function for storing variance/covariance matrices. Use
  `effect_name = "gen_add"` for additive genetic (co)variances and
  `effect_name = "residual"` for residual (co)variances. Any named random effect
  (e.g. `"litter"`, `"dam"`) is also supported. Replaces `set_residual_cov()` and
  `set_random_effect_cov()`.

* **`set_residual_cov()` removed** — use
  `add_effect_cov_matrix(pop, "residual", R)` instead.

* **`set_random_effect_cov()` removed** — use
  `add_effect_cov_matrix(pop, effect_name, R)` instead.

* **`trait_meta` columns removed**: `target_add_var` and `residual_var` are no
  longer stored in `trait_meta`. Supply them via `add_trait(target_add_var = ...)` /
  `add_trait(residual_var = ...)` (values are written to `trait_effect_cov`) or
  call `add_effect_cov_matrix()` directly.

* **`trait_effects` column removed**: `variance` column removed from
  `trait_effects`. Random effect variances are now stored exclusively in
  `trait_effect_cov`.

* **`set_qtl_effects_multi()` — `G` parameter is now optional** (default `NULL`).
  If omitted, the additive genetic covariance matrix is read from
  `trait_effect_cov` (stored via `add_effect_cov_matrix("gen_add", ...)`).

* **`add_effect_random()` — `variance` parameter is now optional** (default
  `NULL`). If a diagonal entry for `(effect_name, trait_name)` already exists in
  `trait_effect_cov`, it is used automatically. Only required when no stored
  value exists.

## New Tables

* **`trait_effect_cov`** — unified variance/covariance table replacing
  `trait_residual_cov` and `trait_random_effect_cov`. Schema:
  `(effect_name VARCHAR, trait_1 VARCHAR, trait_2 VARCHAR, cov DOUBLE)`.
  Both `(i,j)` and `(j,i)` pairs stored for symmetric lookup.

## Bug Fixes / Internal Changes

* All writes to `trait_effect_cov` use `DBI::dbExecute()` with raw SQL instead
  of `DBI::dbWriteTable()` to avoid consuming R's RNG state (DuckDB's
  `dbWriteTable` internally touches the RNG, which would shift simulation
  reproducibility).

---

# tidybreed 0.6.0 (2026-04-23)

## New Features

* **`add_effect_int(pop, trait_name, mean)`** — sets the intercept
  (`target_add_mean`) for a trait; convenience alternative to setting it at
  `add_trait()` time.

* **`add_effect_fixed_class(pop, trait_name, effect_name, source_column, levels, source_table, overwrite)`**
  — discrete fixed effect. Maps levels of a grouping column to numeric shifts.
  Errors (not silent 0) at phenotyping time if an individual's level is not in
  `levels`. Replaces the `effect_class = "fixed"` path of `add_trait_covariate()`.

* **`add_effect_fixed_cov(pop, trait_name, effect_name, source_column, slope, center, source_table, overwrite)`**
  — continuous covariate regression term. Contribution = `slope * (x - center)`.
  When `center = NULL` the mean of `source_column` is computed from the current
  `source_table` and stored for reproducibility.

* **`add_effect_random(pop, trait_name, effect_name, source_column, variance, distribution, source_table, overwrite)`**
  — random group effect. Drawn values are now persisted in a new
  `trait_random_effects` table so the same group receives the same shift on
  repeated calls to `add_phenotype()`, without requiring a fixed `seed`.

* **`set_random_effect_cov(pop, effect_name, traits, R, tol, overwrite)`** —
  stores a covariance matrix enabling joint MVN draws of a named random effect
  across multiple traits. Analogous to `set_residual_cov()`.

* All `add_effect_*()` functions accept a `source_table` parameter (default
  `"ind_meta"`). Any database table with an `id_ind` column can serve as the
  source for effect levels, enabling future repeated-measures and multi-table
  workflows.

## Schema Changes

* `trait_effects` gains three nullable columns: `source_table VARCHAR`,
  `slope DOUBLE`, `center DOUBLE`. Existing rows (written by
  `add_trait_covariate()`) are backward-compatible — `source_table` defaults to
  `"ind_meta"` when `NULL`.

* Two new tables: `trait_random_effects` (stores drawn group-level random effect
  values) and `trait_random_effect_cov` (covariance structure for correlated
  random effects across traits).

## Deprecations

* `add_trait_covariate()` now emits a deprecation warning. It continues to work
  for backward compatibility but users should migrate to `add_effect_fixed_class()`
  or `add_effect_random()`.

---

# tidybreed 0.5.0 (2026-04-23)

## Breaking Changes

* `add_trait()`: parameter `mean` renamed to `target_add_mean`. The column in
  `trait_meta` is also renamed. Existing serialized databases with the old
  column name will not load correctly — re-initialize when upgrading.

* `set_qtl_effects()` and `set_qtl_effects_multi()`: first argument renamed from
  `pop` to `x` and now accepts either a `tidybreed_pop` or a `tidybreed_table`
  (for `base = "current_pop"` to specify which individuals define the base).

## New Features

* `set_qtl_effects()` / `set_qtl_effects_multi()`: new `base` parameter
  (`"founder_haplotypes"` default, or `"current_pop"`). Controls which allele
  frequencies are used for effect scaling and TBV centering.

* Effect scaling now uses the **Falconer formula**
  (`V_A = Σ 2·p·(1−p)·α²`) instead of the realized `var(G %*% alpha)`. This
  gives a theoretically grounded `target_add_var` guarantee rather than an
  empirical one.

* `set_qtl_effects()` writes a `base_allele_freq_{trait}` column to
  `genome_meta` recording which allele frequencies were used. This column is
  read by `add_tbv()` and `add_phenotype()` to center TBVs:
  `TBV_i = (G_i − 2·p_base) · α`, ensuring `E[TBV] ≈ 0` for the base
  population.

* `base = "founder_haplotypes"` computes allele frequencies from the actual
  rows in the `founder_haplotypes` table (the pool used to sample founders),
  not from the stored theoretical `founder_allele_freq` column.

---

# tidybreed 0.4.3 (2026-04-23)

## Internal / Cleanup

* Action functions (`add_phenotype()`, `add_tbv()`, `add_genotypes()`,
  `extract_genotypes()`) fully migrated to the `tidybreed_table`-first calling
  convention announced in 0.4.1. Removed the internal `resolve_pending_filter()`
  helper; each function now collects and extracts `id_ind` directly from the
  `tidybreed_table`.

* `define_chip()` and `define_qtl()` now call `mutate_table()` directly instead
  of relying on the removed `mutate_genome_meta()`.

* `add_offspring()` now uses `validate_sql_identifier()` for extra column
  validation, consistent with the rest of the package.

* Source files for deleted functions (`R/mutate_genome_meta.R`,
  `R/mutate_ind_meta.R`, `R/add_ebv.R`) and their associated man pages and
  tests are now physically removed from the repository.

* Expanded roxygen documentation: new `man/` pages for many exported and
  internal functions, and updated `get_table()` / `infer_duckdb_type()` docs.

* Test suite updated throughout to use the new `get_table() |> filter() |>
  action()` pattern.

# tidybreed 0.4.2 (2026-04-23)

## New Features

* `get_table()` chains now support the full `slice_*` family: `slice_max()`,
  `slice_min()`, `slice_head()`, `slice_tail()`, and `slice_sample()`.
  These can be used directly on a `tidybreed_table` before passing to action
  functions such as `add_phenotype()`.

# tidybreed 0.4.1 (2026-04-23)

## Breaking Changes

* `filter()` can no longer be called directly on a `tidybreed_pop` object.
  You must now call `get_table()` first to identify which table to filter.
  This applies to `add_phenotype()`, `add_tbv()`, `add_genotypes()`, and
  `extract_genotypes()`.
  ```r
  # Old (no longer works):
  pop |> dplyr::filter(sex == "F") |> add_phenotype("ADG")

  # New (required):
  pop |> get_table("ind_meta") |> dplyr::filter(sex == "F") |> add_phenotype("ADG")

  # Any table with id_ind now works:
  pop |> get_table("ind_phenotype") |> dplyr::filter(value > 500) |> add_phenotype("ADG2")
  ```

* `add_phenotype()`, `add_tbv()`, `add_genotypes()`, and `extract_genotypes()`
  now require a `tidybreed_table` as the first argument (from `get_table()`),
  not a `tidybreed_pop`.

* `add_ebv()` has been removed. A proper evaluation runner will replace it in
  a future version.

# tidybreed 0.4.0 (2026-04-23)

## Breaking Changes

* `mutate_ind_meta()` and `mutate_genome_meta()` have been **removed**. Use the
  new generic `mutate_table()` instead:
  ```r
  # Old:
  pop <- mutate_ind_meta(pop, gen = 1L)
  # New:
  pop <- pop |> get_table("ind_meta") |> mutate_table(gen = 1L)
  ```

## New Functions

* `mutate_table(tbl_obj, ...)` — generic column add/update for any table.
  Chain after `get_table()` (and optionally `filter()`) to add new columns or
  update existing ones. Returns `pop` invisibly. Supports scalar and vector
  values, type inference, reserved-column blocking, and informative messages.
  ```r
  # All rows:
  pop <- pop |> get_table("ind_meta") |> mutate_table(gen = 1L)

  # Filtered rows only (females get NULL for a new column):
  pop <- pop |>
    get_table("ind_meta") |>
    filter(sex == "M") |>
    mutate_table(gen = 2L)
  ```

## Changes

* `get_table()` now returns a `tidybreed_table` S3 object instead of a raw
  `tbl_duckdb_connection`. Backward-compatible via `collect.tidybreed_table`,
  `filter.tidybreed_table`, `select.tidybreed_table`, `arrange.tidybreed_table`,
  `pull.tidybreed_table`, and `count.tidybreed_table` S3 methods — all existing
  `get_table(...) |> dplyr::filter(...) |> dplyr::collect()` patterns continue
  to work unchanged.
* `infer_duckdb_type()` moved to `R/sql_utils.R` and is now a shared internal
  utility available to all mutation helpers.
* `define_chip()`, `define_qtl()`, `set_qtl_effects()`, and
  `set_qtl_effects_multi()` updated to call `mutate_table()` internally.

---

# tidybreed 0.3.0 (2026-04-21)

## New Functions

* `add_genotypes(pop, chip_name)` — marks a filtered subset of animals as
  genotyped on a named SNP chip by writing a `has_<chip_name>` BOOLEAN column
  to `ind_meta`. Follows the same `filter()` -> action pipe pattern as
  `add_phenotype()`. Operation is additive: animals already marked TRUE remain
  TRUE across multiple calls. Chip must exist in `genome_meta` (via
  `define_chip()`) before calling.
* `extract_genotypes(pop, chip_name)` — returns a tibble of genotypes
  (0/1/2 encoding) for animals marked as genotyped, restricted to chip loci.
  The returned set is the intersection of animals with `has_<chip_name> == TRUE`,
  any pending `filter()` predicates, and loci with `is_<chip_name> == TRUE`.
  Intended for use immediately before GBLUP/GWAS evaluation.

---

# tidybreed 0.2.2 (2026-04-20)

## Performance

- `add_founders()`: eliminated nested R loop (O(n_founders × n_loci) iterations)
  that caused multi-minute runtimes for large genomes. Haplotype and genotype
  frames are now built via vectorized matrix indexing and addition in C, reducing
  frame construction from minutes to ~0.18 s regardless of genome size. Also
  switched genome table writes to `duckdb_register` + `INSERT SELECT` for a
  further ~2× speedup on the DB write step. Typical runtimes: 0.38 s (1k loci)
  → 2.7 s (5k loci) → 8 s (10k loci) for 215 founders.

# tidybreed 0.2.1 (2026-04-20)

## API

- Renamed `name` argument to `trait_name` in `add_trait()` and `add_trait_simple()` for consistency
- Renamed `trait` argument to `trait_name` in `define_qtl()` and `set_qtl_effects()`
- Renamed `traits` argument to `trait_names` in `set_qtl_effects_multi()`
- Renamed `name` argument to `chip_name` in `define_chip()`

# tidybreed 0.2.0 (2026-04-20)

## New: Trait Architecture & Phenotype Simulation

Trait definition, QTL selection, additive-effect sampling (single and
correlated multi-trait), fixed and random covariates, phenotype generation
for continuous / count / binary / categorical traits, imprinting support,
and storage of true and estimated breeding values.

### New exported functions

* `add_trait()` — insert a row in `trait_meta` with target variance
  components, trait type (continuous / count / binary / categorical),
  expression rules (sex-limited, parent-of-origin), and index/economic
  weights. Creates the six trait-layer tables on first call.
* `define_qtl()` — mirror of `define_chip()` for QTL loci; writes an
  `is_QTL_{trait}` BOOLEAN column in `genome_meta`. Reuses the same six
  selection methods (by count + `random` / `even` / `chromosome_even`, by
  logical vector, by locus ids, by locus names).
* `set_qtl_effects()` — write the `add_{trait}` additive-effect column.
  Supports manual effects or sampled effects (`normal` / `gamma`) with
  automatic rescaling to hit `trait_meta$target_add_var`.
* `set_qtl_effects_multi()` — draw correlated additive effects across
  multiple traits from `MVN(0, G)` via `MASS::mvrnorm`; supports `"shared"`
  (pleiotropy) and `"union"` strategies.
* `set_residual_cov()` — store a residual covariance matrix `R` across
  traits in a new `trait_residual_cov` table. Consumed by `add_phenotype()`
  when multiple traits share the same filtered subset.
* `add_trait_covariate()` — append fixed or random covariate rows to a
  `trait_effects` table. Fixed-effect levels are serialised to a JSON-style
  VARCHAR; random effects store distribution + variance.
* `add_phenotype()` — the workhorse. Generates phenotypes for a subset of
  individuals for one or more traits and writes rows to `ind_phenotype`. Also
  computes and stores the underlying TBV in `ind_tbv`. Joint MVN residual
  draws when multiple traits share the subset and `R` is stored.
* `add_tbv()` — compute and store TBV without generating phenotype records.
* `add_ebv()` — ingest externally computed estimated breeding values into
  `ind_ebv`, tagged with a user-supplied model label.
* `add_trait_simple()` — one-shot wrapper chaining `add_trait()` +
  `define_qtl()` + `set_qtl_effects()`.

### filter() on tidybreed_pop

New S3 method `filter.tidybreed_pop()` stashes dplyr predicates on the
population object. The next `add_phenotype()` / `add_tbv()` call applies and
clears them. Multiple `filter()` calls stack with AND semantics.

```r
pop |>
  dplyr::filter(sex == "F", gen == 1L) |>
  add_phenotype("ADG")
```

### New tables

`trait_meta`, `trait_effects`, `trait_residual_cov`, `ind_phenotype`,
`ind_tbv`, `ind_ebv`. Also new columns written to `genome_meta`:
`is_QTL_{trait}` (BOOLEAN) and `add_{trait}` (DOUBLE) per trait.

### Other

* Added `MASS` to Imports for multivariate normal sampling. MASS is a
  Recommended R package that ships with base R distributions — no extra
  install burden.

---

# tidybreed 0.1.0 (2026-04-17)

## New Functions

* `add_offspring()` — core mating function. User supplies a `matings` tibble
  (one row per offspring) with required columns `id_parent_1`, `id_parent_2`,
  `sex`, and `line`. Gametes are produced via chromosomal crossover simulation
  (Haldane map: crossovers per chromosome ~ Poisson(chr_len_Mb / 100)).
  New `ind_meta`, `genome_haplotype`, and `genome_genotype` rows are written
  for all offspring. Animal-breeder aliases `id_sire` / `id_dam` are accepted.
  Any extra columns in `matings` (e.g. `gen = 2L`) are validated and written
  to `ind_meta`, with automatic `ALTER TABLE` if the column is new.
* `make_gamete()` — internal recombination helper (`R/recombination_helpers.R`).
  Not exported; used by `add_offspring()`.

---

# tidybreed 0.0.3 (2026-04-17)

## Breaking Changes

* `ind_meta`, `genome_haplotype`, `genome_genotype`: column `ind_id` renamed to
  `id_ind`; `parent_1` → `id_parent_1`; `parent_2` → `id_parent_2`. The `id_`
  prefix groups all identifier columns together for clarity.

**Migration for existing databases:**
```r
conn <- DBI::dbConnect(duckdb::duckdb(), "your_database.duckdb")
DBI::dbExecute(conn, "ALTER TABLE ind_meta RENAME COLUMN ind_id TO id_ind")
DBI::dbExecute(conn, "ALTER TABLE ind_meta RENAME COLUMN parent_1 TO id_parent_1")
DBI::dbExecute(conn, "ALTER TABLE ind_meta RENAME COLUMN parent_2 TO id_parent_2")
DBI::dbExecute(conn, "ALTER TABLE genome_haplotype RENAME COLUMN ind_id TO id_ind")
DBI::dbExecute(conn, "ALTER TABLE genome_genotype RENAME COLUMN ind_id TO id_ind")
DBI::dbDisconnect(conn)
```

## Bug Fixes

* `infer_duckdb_type()`: bare `NA` (untyped logical) now warns and returns
  `VARCHAR` instead of silently returning `BOOLEAN`
* `mutate_ind_meta()`: type check now runs before length check, so passing an
  unsupported type (e.g. a `list`) produces "Unsupported type" rather than a
  misleading length-mismatch error
* Fixed pre-existing test bugs: `expect_error(expr, NA)` inverted assertion in
  `test-add_founders.R`; null `db_conn` in `test-mutate_genome_meta.R` now uses
  a real in-memory connection

# tidybreed 0.0.2 (2026-04-15)

## Bug Fixes

* `infer_duckdb_type()`: fixed incorrect VARCHAR inference when a typed vector
  starts with `NA` (e.g. `c(NA_real_, rnorm(n))`). Type is now determined from
  the R class of the full vector before inspecting individual elements, so
  `NA_real_` → `DOUBLE`, `NA_integer_` → `INTEGER`, `NA` → `BOOLEAN`. A bare
  all-`NA` vector (class `logical`) now warns with a message pointing users to
  typed NA constants.

# tidybreed 0.0.1 (2026-04-14)

## Breaking Changes

### Terminology: "population" → "line"

* `add_founders()` parameter renamed: `pop_name` → `line_name`
* `ind_meta` table column renamed: `population` → `line`

**Rationale:** "population" refers to one complete genome build from
`initialize_genome()`. "Line" refers to distinct genetic lines within that
population (e.g., line A selected for growth, line B for egg quality).

**Migration for existing databases:**
```r
conn <- DBI::dbConnect(duckdb::duckdb(), "your_database.duckdb")
DBI::dbExecute(conn, "ALTER TABLE ind_meta RENAME COLUMN population TO line")
DBI::dbDisconnect(conn)
```

## New Functions

* `mutate_genome_meta()` — add or update user columns in `genome_meta`;
  reserved columns (`locus_id`, `locus_name`, `chr`, `chr_name`, `pos_Mb`)
  are blocked
* `define_chip()` — convenience wrapper that marks loci as members of a named
  SNP chip by writing a `BOOLEAN` column `is_{chip_name}` to `genome_meta`;
  supports `"random"`, `"even"`, and `"chr_even"` selection methods

## Other Changes

* Removed stale design and planning documents (DESIGN.md, QUICKSTART.md,
  IMPLEMENTATION_STATUS.md, IMPLEMENTATION_SUMMARY.md, TODO_finalize.md,
  MUTATE_IND_META_IMPLEMENTATION.md)
* Moved manual test scripts to `tests/` directory
* Added CLAUDE.md with developer and AI context for the project

---

# tidybreed 0.0.0

* Initial package setup and framework

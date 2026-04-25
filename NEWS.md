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

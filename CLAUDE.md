# tidybreed — Developer & AI Context

## What This Package Does

`tidybreed` simulates animal and plant breeding programs. It provides a
pipe-friendly API where all genomic and individual data is stored in a
file-based DuckDB database (not in R memory). Users build a `tidybreed_pop`
object step by step, query tables with `dplyr`, and eventually run selection
and mating cycles.

## Design Principles

1. **Database-first** — all data lives in DuckDB, not R objects
2. **File-based by default** — enables simulations larger than RAM; supports
   resumable runs and sharing
3. **Lazy evaluation** — use `dplyr::tbl()` / `get_table()` and filter before
   `collect()`-ing into R
4. **Pipe-friendly** — most exported functions accept a `tidybreed_pop` and
   return a `tidybreed_pop`; action functions (`add_phenotype`, `add_tbv`,
   `add_genotypes`, `extract_genotypes`) accept a `tidybreed_table` from
   `get_table()` and return `tidybreed_pop`
5. **Type-safe** — all table columns have explicit DuckDB types; user-added
   columns are inferred via `infer_duckdb_type()`

## Function Naming Convention

| Prefix        | Meaning                                        |
|---------------|------------------------------------------------|
| `initialize_` | Creates the DuckDB database and core tables    |
| `add_`        | Inserts new rows (individuals, traits, etc.)   |
| `mutate_`     | Adds or updates columns in an existing table   |
| `define_`     | Convenience wrapper for common `mutate_` tasks |

## Database Schema

### `genome_meta`

Locus-level metadata. One row per locus.

| Column     | Type    | Notes                          |
|------------|---------|--------------------------------|
| locus_id   | INTEGER | Primary key (1, 2, …, n_loci)  |
| locus_name | VARCHAR | e.g. "Locus_1", "rs12345"      |
| chr        | INTEGER | Chromosome number              |
| chr_name   | VARCHAR | Chromosome name string         |
| pos_Mb     | DOUBLE  | Position in megabases          |
| *user cols*| any     | Added via `mutate_genome_meta()` or `define_chip()` |

**Reserved** (cannot be modified): `locus_id`, `locus_name`, `chr`, `chr_name`, `pos_Mb`

Example user columns: `is_50K BOOLEAN`, `is_HD BOOLEAN`, `effect_ADG DOUBLE`

### `genome_haplotype`

Phased haplotypes. Two rows per individual (parent_origin 1 = paternal, 2 = maternal).

| Column        | Type    |
|---------------|---------|
| id_ind        | VARCHAR |
| parent_origin | INTEGER |
| locus_1 … locus_n | INTEGER (0 or 1) |

### `genome_genotype`

Genotypes in 0/1/2 encoding. One row per individual.

| Column        | Type    |
|---------------|---------|
| id_ind        | VARCHAR |
| locus_1 … locus_n | INTEGER (0, 1, or 2) |

Both haplotypes and genotypes are stored (not computed on demand) because disk
is cheap and recomputation during mating is expensive.

### `ind_meta`

Individual-level metadata. Created by `add_founders()`.

| Column      | Type    | Notes                          |
|-------------|---------|--------------------------------|
| id_ind      | VARCHAR | Primary key, format `{line}-{n}` |
| id_parent_1 | VARCHAR | NA for founders                |
| id_parent_2 | VARCHAR | NA for founders                |
| line        | VARCHAR | Genetic line name              |
| sex         | VARCHAR | "M" or "F"                     |
| *user cols* | any     | Added via `mutate_ind_meta()`  |

**Reserved**: `id_ind`, `id_parent_1`, `id_parent_2`, `line`, `sex`

### `trait_meta`

One row per trait. Populated by `add_trait()`.

| Column             | Type    | Notes                                                         |
|--------------------|---------|---------------------------------------------------------------|
| trait_name         | VARCHAR | Primary key                                                   |
| description        | VARCHAR | Free text                                                     |
| units              | VARCHAR | e.g. "kg"                                                     |
| trait_type         | VARCHAR | `"continuous"` / `"count"` / `"binary"` / `"categorical"`     |
| repeatable         | BOOLEAN | Repeated measures allowed?                                    |
| recorded_on        | VARCHAR | `"self"` / `"dam"` / `"sire"` / `"offspring_mean"`            |
| expressed_sex      | VARCHAR | `"both"` / `"M"` / `"F"` for sex-limited traits               |
| expressed_parent   | VARCHAR | `"both"` / `"parent_1"` / `"parent_2"` for imprinting         |
| mean               | DOUBLE  | Overall mean / intercept on liability scale                   |
| min_value          | DOUBLE  | For count traits; clip (NA = no limit)                        |
| max_value          | DOUBLE  | For count traits; clip                                        |
| prevalence         | DOUBLE  | For binary traits                                             |
| thresholds         | VARCHAR | For categorical traits; comma-separated cutpoints             |
| index_weight       | DOUBLE  | Weight in selection index                                     |
| economic_value     | DOUBLE  | Economic value per unit                                       |

### `trait_effects`

Non-additive-genetic, non-residual terms in the phenotype model (fixed and
random effects). One row per (trait × effect).

| Column        | Type    | Notes                                                    |
|---------------|---------|----------------------------------------------------------|
| trait_name    | VARCHAR |                                                          |
| effect_name   | VARCHAR | e.g. "sex", "gen", "litter"                              |
| effect_class  | VARCHAR | `"fixed_class"`, `"fixed_cov"`, or `"random"`            |
| source_column | VARCHAR | Column in source table used as grouping variable         |
| source_table  | VARCHAR | Table containing `source_column` (default `"ind_meta"`)  |
| distribution  | VARCHAR | For random effects: `"normal"`, `"gamma"`, `"uniform"`   |
| levels_json   | VARCHAR | For fixed_class effects: JSON `{"M":30,"F":0}`           |
| slope         | DOUBLE  | For fixed_cov effects: regression coefficient            |
| center        | DOUBLE  | For fixed_cov effects: centering value                   |
| value         | DOUBLE  | Rarely used scalar                                       |

### `trait_effect_cov`

Unified variance/covariance table for all random effects (additive genetic,
residual, and named random effects). One row per (effect_name, trait_1, trait_2).
Both `(i,j)` and `(j,i)` pairs stored. Populated by `add_effect_cov_matrix()`,
`add_trait()`, and `add_effect_random()`.

Reserved `effect_name` values: `"gen_add"` (additive genetic G matrix),
`"residual"` (residual R matrix). Any other name maps to a user-defined random
effect matching `effect_name` in `trait_effects`.

| Column      | Type    | Notes                                              |
|-------------|---------|----------------------------------------------------|
| effect_name | VARCHAR | `"gen_add"`, `"residual"`, or random effect name   |
| trait_1     | VARCHAR | Primary key with effect_name + trait_2             |
| trait_2     | VARCHAR |                                                    |
| cov         | DOUBLE  | Variance (diagonal) or covariance (off-diagonal)   |

### `ind_phenotype`

Phenotype records in long format. Populated by `add_phenotype()`.

| Column      | Type    | Notes                                             |
|-------------|---------|---------------------------------------------------|
| id_record   | VARCHAR | Auto `{trait}-{n}`                                |
| id_ind      | VARCHAR |                                                   |
| trait_name  | VARCHAR |                                                   |
| value       | DOUBLE  | Phenotype value                                   |
| pheno_number| INTEGER | 1 = first record for this individual × trait, etc.|
| *user cols* | any     | Added via `mutate_table()`                        |

### `ind_tbv`

True breeding values (simulation ground truth). Populated by
`add_phenotype()` and `add_tbv()`. Composite key `(id_ind, trait_name)`.

| Column     | Type    |
|------------|---------|
| id_ind     | VARCHAR |
| trait_name | VARCHAR |
| tbv        | DOUBLE  |
| date_calc  | DATE    |

### `ind_ebv`

Estimated breeding values from external BLUP / GBLUP runs. Composite key
`(id_ind, trait_name, model)`. Will be populated by a future `add_ebv()`
replacement (current version removed).

| Column     | Type    | Notes                          |
|------------|---------|--------------------------------|
| id_ind     | VARCHAR |                                |
| trait_name | VARCHAR |                                |
| model      | VARCHAR | User label, e.g. "ssGBLUP_v1"  |
| ebv        | DOUBLE  |                                |
| acc        | DOUBLE  | Optional accuracy              |
| se         | DOUBLE  | Optional standard error        |
| date_calc  | DATE    |                                |

## Implemented Functions (Phase 1 Complete)

### `initialize_genome()`

`R/initialize_genome.R`

Creates the DuckDB file (or in-memory DB) and populates three tables:
`genome_meta`, `genome_haplotype` (empty), `genome_genotype` (empty).
Returns a `tidybreed_pop` S3 object.

Key params: `pop_name`, `n_loci`, `n_chr`, `chr_len_Mb`, `db_path`

### `add_founders()`

`R/add_founders.R`

Samples haplotypes for each founder individual using per-locus allele
frequencies. Populates `ind_meta` (core 5 cols), `genome_haplotype` (2 rows
each), and `genome_genotype` (1 row each). ID format: `{line_name}-{n}`.

Key params: `n_males`, `n_females`, `line_name`, `allele_freq`

### `mutate_table()`

`R/mutate_table.R`

Generic replacement for the removed `mutate_ind_meta()` / `mutate_genome_meta()`
functions. Adds or updates columns in **any** database table. Chain after
`get_table()` (and optionally `filter()`) then call `mutate_table(col = value)`.
Scalar values are broadcast to all (or filtered) rows; vectors must match the
effective row count. Type is inferred via `infer_duckdb_type()`. Reserved
columns are blocked via `TABLE_RESERVED_COLS` in `R/sql_utils.R`. Returns `pop`
invisibly.

```r
# All rows
pop <- pop |> get_table("ind_meta") |> mutate_table(gen = 1L)

# Filtered rows (unmatched rows get NULL for new columns)
pop <- pop |>
  get_table("ind_meta") |>
  filter(sex == "M") |>
  mutate_table(gen = 2L)
```

### `define_chip()`

`R/define_chip.R`

Convenience wrapper. Marks loci as members of a named chip by writing a
`BOOLEAN` column `is_{chip_name}` to `genome_meta`. Four selection methods:

- `n` + `method = "random"` — randomly sample `n` loci
- `n` + `method = "even"` — evenly spaced across all loci by position
- `n` + `method = "chromosome_even"` — proportional to chromosome length
- `locus_tf` — logical vector (same length as `genome_meta` rows); TRUE = on chip. Useful for passing the complement of an existing chip (`!existing_tf`).
- `locus_ids` — integer vector of specific locus IDs
- `locus_names` — character vector of specific locus names

### `get_table()` / `close_pop()` / `print.tidybreed_pop()`

`R/tidybreed_pop.R`

`get_table(pop, "table_name")` returns a `tidybreed_table` S3 object that
carries the pop reference, table name, lazy dplyr tbl, and pending filter.
Supports `filter()`, `collect()`, `select()`, `arrange()`, `pull()`, `count()`,
and `mutate_table()`. `close_pop()` safely closes the DuckDB connection.

**Subset selection for action functions** (`add_phenotype`, `add_tbv`,
`add_genotypes`, `extract_genotypes`) requires `get_table()` as the first step.
`filter()` is called on the `tidybreed_table`, not on the pop directly.
The unique `id_ind` values from the collected filtered table are used as the
candidate set. Any table that has an `id_ind` column can be used (e.g.
`ind_meta`, `ind_phenotype`, `genome_haplotype`, `genome_genotype`).

```r
# All individuals
pop |> get_table("ind_meta") |> add_phenotype("ADG")

# Filtered by ind_meta
pop |>
  get_table("ind_meta") |>
  dplyr::filter(sex == "F", gen == 1L) |>
  add_phenotype("ADG")

# Pre-select top performers from a prior phenotype
pop |>
  get_table("ind_phenotype") |>
  dplyr::filter(value > 500) |>
  add_phenotype("ADG2")
```

### `add_trait()` / `define_qtl()` / `set_qtl_effects()` / `set_qtl_effects_multi()`

`R/add_trait.R`, `R/define_qtl.R`, `R/set_qtl_effects.R`

- `add_trait()` — one row in `trait_meta`. Creates trait tables on first call.
  `target_add_var` and `residual_var` params write to `trait_effect_cov` (not
  `trait_meta`).
- `define_qtl()` — mirror of `define_chip()` for QTL. Writes
  `is_QTL_{trait}` BOOLEAN. Reuses all `chip_helpers.R` selection methods.
- `set_qtl_effects()` — writes `add_{trait}` DOUBLE. Manual or sampled
  (normal/gamma) with optional rescale. Reads `target_add_var` from
  `trait_effect_cov` (`effect_name = "gen_add"`).
- `set_qtl_effects_multi()` — correlated effects from `MVN(0, G)` across
  multiple traits via `MASS::mvrnorm`. `G = NULL` reads from `trait_effect_cov`.
  `method = "shared"` (pleiotropy) or `"union"`.

### `add_effect_cov_matrix()` / `add_effect_random()` / `add_effect_fixed_class()` / `add_effect_fixed_cov()` / `add_effect_int()`

`R/add_effect_cov_matrix.R`, `R/add_effect_random.R`, `R/add_effect_fixed_class.R`, `R/add_effect_fixed_cov.R`, `R/add_effect_int.R`

- `add_effect_cov_matrix(pop, effect_name, cov_matrix)` — **single entry
  point for all variance/covariance data**. Stores symmetric matrix in
  `trait_effect_cov`. Use `effect_name = "gen_add"` or `"residual"` for
  reserved effects. Can be called before `add_trait()` or `add_effect_random()`.
- `add_effect_random()` — `variance` optional if already in `trait_effect_cov`.
- `add_effect_fixed_class()` — discrete level → shift mapping.
- `add_effect_fixed_cov()` — linear regression term (`slope * (x - center)`).
- `add_effect_int()` — sets intercept (`target_add_mean`) for a trait.

### `add_trait_covariate()` *(deprecated)*

`R/add_trait_covariate.R`

- Deprecated since v0.6.0. Use `add_effect_fixed_class()`, `add_effect_fixed_cov()`,
  or `add_effect_random()` instead.

### `add_phenotype()` / `add_tbv()`

`R/add_phenotype.R`, `R/add_tbv.R`

Both functions accept a `tidybreed_table` (from `get_table()` + optional
`filter()`) as their first argument and return `tidybreed_pop`.

- `add_phenotype()` — the workhorse. Extracts unique `id_ind` from the
  filtered table, intersects with `expressed_sex`, computes TBV
  (genotype × effects, or haplotype dose for imprinted traits), adds
  fixed/random covariate contributions, samples residuals (joint `MVN(0, R)`
  when multiple traits share the subset and `R` is stored; otherwise
  independent). Converts liability to phenotype per `trait_type`. Writes
  `ind_phenotype` rows and updates `ind_tbv`.
- `add_tbv()` — TBV-only; no phenotype records.

### `add_trait_simple()`

`R/add_trait_simple.R`

Convenience wrapper that chains `add_trait()` → `define_qtl()` →
`set_qtl_effects()` for a single uncorrelated trait.

## Roadmap

### Longer-Term

- `select_parents()` — selection index or truncation selection
- Export: PLINK `.bed/.bim/.fam`, VCF
- Visualization helpers
- `mutate_ind_phenotype()` — user columns on `ind_phenotype`
- Dominance and epistasis effects (currently only additive)

## Versioning Policy

Use **three-part semantic versioning**: `MAJOR.MINOR.PATCH` (e.g. `0.0.1`).
Do **not** use the four-part devtools convention (`0.0.0.9000`).

| Part  | Bump when…                                           |
|-------|------------------------------------------------------|
| PATCH | Bug fixes, doc updates, minor internal changes       |
| MINOR | New exported functions or non-breaking feature additions |
| MAJOR | Breaking API changes                                 |

**Before every commit + push, update:**
1. `DESCRIPTION` — `Version:` field
2. `NEWS.md` — add an entry under the new version heading

## Key File Locations

| Path                       | Contents                          |
|----------------------------|-----------------------------------|
| `R/`                       | All exported functions            |
| `tests/testthat/`          | Formal testthat test suite        |
| `tests/test_*.R`           | Manual/dev test scripts           |
| `man/`                     | Roxygen-generated documentation   |

## Design Rationale

**Why DuckDB?** Columnar, embedded (no server), SQL via dbplyr, handles
datasets larger than RAM, excellent R integration.

**Why 0/1/2 encoding?** Standard in genomics (PLINK, VCF). Easy to interpret
(count of alternate allele). Efficient integer storage.

**Why store both haplotypes and genotypes?** Haplotypes are required for
recombination and phased exports. Genotypes are required for GWAS and genomic
prediction. Computing genotypes on the fly during every query would be wasteful.

## Development Environment

### Running R Commands

The R executable path is platform-specific. When running R or Rscript via the
Bash tool, use the appropriate path:

**Windows (Hendrix Genetics AVD):**
- Check if working directory contains "Hendrix"
- R: `"/c/Program Files/R/R-4.5.1/bin/x64/R.exe"`
- Rscript: `"/c/Program Files/R/R-4.5.1/bin/x64/Rscript.exe"`

**Mac/Linux:**
- Use standard shell commands: `R` or `r` and `Rscript`

Example test command:
```bash
"/c/Program Files/R/R-4.5.1/bin/x64/Rscript.exe" -e "print('Hello from R')"
```





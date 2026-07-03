# tidybreed — Developer & AI Context

## What This Package Does

`tidybreed` simulates breeding programs. It provides a
pipe-friendly API where all genomic and individual data is stored in a
file-based DuckDB database (not in R memory). Users build a `tidybreed_pop`
object step by step, query tables with `dplyr`, and eventually run selection
and mating cycles.

## Design Principles

1. **Database-first** — all data lives in DuckDB, not R objects
2. **File-based by default** — enables simulations larger than RAM; supports
   replicates, resumable runs, and sharing easily
3. **Lazy evaluation** — use `dplyr::tbl()` / `get_table()` and filter before
   `collect()`-ing into R
4. **Pipe-friendly** — most exported functions accept a `tidybreed_pop` and
   return a `tidybreed_pop`; action functions (`add_phenotype`, `add_tbv`,
   `add_genotypes`, `extract_genotypes`, `define_chip`, `define_additive_effects`)
   accept a `tidybreed_table` from `get_table()` and return `tidybreed_pop`
5. **Type-safe** — all table columns have explicit DuckDB types; user-added
   columns are inferred via `infer_duckdb_type()`
6. **Disdain and intolerance for storing metadata** - storing data such as
   n_loci for loci count is silly when the function that needs that data
   can run a basic SQL query to pull out info and makes the program truly modular
7. **No implicit ordering** — never accept positional vectors (logical TRUE/FALSE
   or integer indices) to select rows from a database table. Row order is not
   guaranteed and changes silently when rows are added or reordered. All
   selection flows through named identifiers (`id_ind`, `locus_name`, `locus_id`
   as PK) or explicit SQL filter predicates via `get_table() |> filter()`.

If metadata is required, it needs to be considered carefully and likely 
should be stored in a table and never in the R object itself, for restarts
and for `restore_pop()` to be complete without having to recreate it or
given as options again. 

## Naming Consistency Rules

1. **Function argument names must match the database column they populate exactly.**
   Never use a different name for the same thing (e.g. argument `line_name` → column
   `line_name`, not `line`).

2. **All primary numeric value columns follow the `{prefix}_value` pattern.**
   Examples: `pheno_value` in `ind_phenotype`, `tbv_value` in `ind_tbv`,
   `ebv_value` in `ind_ebv`, `cov_value` in `trait_var_comp`,
   `index_value` in `ind_index`.

3. **All name/label columns end in `_name`.**
   Examples: `trait_name`, `locus_name`, `chr_name`, `line_name`, `effect_name`,
   `index_name`. Never abbreviate to just `trait`, `line`, `chr`, etc. when
   used as a column or function parameter that maps to one of these columns.

4. **All ID foreign-key columns start with `id_`.**
   Examples: `id_ind`, `id_trait`, `id_ebv`, `id_tbv`. No `phenotype_id`-style names.

5. **No abbreviations in column names when the full word is unambiguous.**
   `index_weight` not `index_wt`; `trait_name_1`/`trait_name_2` not `trait_1`/`trait_2`.

## Two-Layer Phenotype Design (v0.31.0+)

The model is split into two distinct layers with a strict boundary between them:

**Genetic component layer** — managed by `define_trait()`:
- One row in `trait_meta` per underlying genetic quantity (e.g. `ADG_direct`, `ADG_social`, `WWD`, `WWM`)
- Has QTL effects in `genome_effects`, TBVs in `ind_tbv`, additive variance in `trait_var_comp`
- Arguments: `target_add_var`, `target_add_mean`, `expressed_parent`, `description`, `units`
- No phenotype-level information at all — no mean, no residual, no type, no expressed_sex

**Observation layer** — managed by `define_phenotype()`:
- One row in `phenotype_meta` per observed phenotype individuals receive records for (e.g. `ADG`, `WW`, `mortality`)
- For simple traits, `phenotype_name` equals the `trait_name` of its single genetic component
- For composite traits (maternal, SGE), `phenotype_name` is new and one or more `trait_meta` rows feed into it via `phenotype_components`
- Arguments: `type`, `mean`, `expressed_sex`, `repeatable`, `min_value`, `max_value`, `prevalence`, `thresholds`, `cat_values`, `cat_names`, `store_liability`, `residual_var`, `components`, `formula_tbv`, `formula`, `missing_component_action`

**The rule**: if an argument describes the genetics (variance, QTL structure, parent-of-origin), it belongs in `define_trait()`. If it describes what observers record (mean, distribution, sex expression, residual noise, how to assemble from components), it belongs in `define_phenotype()`. Never put observation-layer arguments on `define_trait()` or genetic-layer arguments on `define_phenotype()`.

## Function Naming Convention

| Prefix        | Meaning                                        |
|---------------|------------------------------------------------|
| `initialize_` | Creates the DuckDB database and core tables    |
| `add_`        | Inserts simulation output rows                 |
| `define_`     | Writes model configuration / metadata          |
| `mutate_`     | Adds or updates columns in an existing table   |

**`add_*` vs `define_*` rule**: if the function writes rows that represent
simulation *output* (data produced by running the model), use `add_`. If the
function writes rows that configure *how* the model runs (parameters, weights,
effect definitions), use `define_`.

Examples: `add_founders()`, `add_phenotype()`, `add_tbv()`, `add_ebv()`,
`add_index()` — all write simulation output.  
`define_trait()`, `define_additive_effects()`, `define_effect_cov_matrix()`,
`define_chip()`, `define_index()` — all write model configuration.

## Database Schema

### `genome_meta`

Locus-level metadata. One row per locus.

| Column     | Type    | Notes                          |
|------------|---------|--------------------------------|
| locus_id   | INTEGER | Primary key (1, 2, …, n_loci)  |
| locus_name | VARCHAR | e.g. "Locus_1", "rs12345"; validated unique |
| chr        | INTEGER | Chromosome number              |
| chr_name   | VARCHAR | Chromosome name string         |
| pos_Mb     | DOUBLE  | Position in megabases          |
| introduced_gen | INTEGER | NULL = founding locus; generation a novel mutation was added (Stage 5) |
| *user cols*| any     | Added via `mutate_table()` or `define_chip()` |

**Reserved** (cannot be modified): `locus_id`, `locus_name`, `chr`, `chr_name`, `pos_Mb`, `introduced_gen`

Example user columns: `is_50K BOOLEAN`, `is_HD BOOLEAN`

**Note**: QTL effects are **not** stored as columns in `genome_meta`. They live in
the `genome_effects` table (see below). There are no `add_{trait}`, `is_QTL_{trait}`,
or `base_allele_freq_{trait}` columns.

### `genome_effects`

QTL effect data. One row per (locus × trait × effect type × line). Populated by
`define_additive_effects()` (single or multi-trait).

| Column             | Type    | Notes                                                      |
|--------------------|---------|------------------------------------------------------------|
| id_genome_effect   | INTEGER | Primary key (auto-incrementing)                            |
| locus_name         | VARCHAR | FK to `genome_meta.locus_name`                             |
| line_name          | VARCHAR | NULL = population-wide; set for line-specific effects      |
| trait_name         | VARCHAR | FK to `trait_meta.trait_name`                              |
| genome_effect_type | VARCHAR | `"additive"` now; `"dominance"` and others later           |
| genome_value       | DOUBLE  | Effect size                                                |
| base_allele_freq   | DOUBLE  | Base allele frequency used for TBV centering (Falconer)    |

**Reserved**: all columns (the table is managed exclusively by `define_additive_effects()`).

QTL membership is **implicit**: a locus is a QTL for a trait if it has a row in
`genome_effects` for that `(trait_name, genome_effect_type)`. No separate boolean
flag is stored.

### `ind_haplotype`

Phased haplotypes in **long** format. One row per (individual × haplotype ×
locus). Populated by `add_founders()` and `add_offspring()`.

| Column        | Type     | Notes                                                       |
|---------------|----------|-------------------------------------------------------------|
| id_ind        | VARCHAR  | FK to `ind_meta.id_ind`; part of composite PK               |
| parent_origin | UTINYINT | 1 = from parent_1 (sire), 2 = from parent_2 (dam); PK part  |
| strand        | UTINYINT | Copy within a parent's contribution; always 1 for diploids; PK part |
| line_origin   | VARCHAR  | Founding line this allele traces to; used by `add_tbv()` for line-specific crossbreeding TBV |
| locus_id      | INTEGER  | FK to `genome_meta.locus_id`; physical sort/PK key          |
| locus_name    | VARCHAR  | FK to `genome_meta.locus_name`; denormalized for direct `genome_effects` joins |
| allele        | UTINYINT | 0 or 1 (phased)                                             |

**Primary key**: `(id_ind, parent_origin, strand, locus_id)`.

### `ind_genotype`

Genotype dosages in **long** format, 0/1/2 encoding. One row per (individual ×
locus). **On-demand cache** — starts empty and is populated only by
`add_dosage()` (never by `add_founders()`/`add_offspring()`). May be empty or
partial.

| Column       | Type     | Notes                                        |
|--------------|----------|----------------------------------------------|
| id_ind       | VARCHAR  | FK to `ind_meta.id_ind`; part of composite PK |
| locus_id     | INTEGER  | FK to `genome_meta.locus_id`; part of composite PK |
| locus_name   | VARCHAR  | FK to `genome_meta.locus_name`               |
| dosage_value | UTINYINT | Sum of alleles across strands (0/1/2 diploid) |

**Primary key**: `(id_ind, locus_id)`. Populated via `INSERT OR REPLACE`
(idempotent).

Haplotypes are the source of truth; dosage is derived on demand (SUM of alleles)
rather than auto-stored, because dosage is cheap to recompute and only needed for
specific downstream analyses (MAS, GBLUP export, allele frequencies).

### `chr_meta`

Per-chromosome inheritance rules. One row per chromosome. Created by
`define_genome()` with default diploid-autosome rows. `define_chr()` and
non-default rules (sex chromosomes, organelles, polyploidy) are Stage 4.

| Column      | Type     | Notes                                                    |
|-------------|----------|----------------------------------------------------------|
| chr_name    | VARCHAR  | Primary key; matches `genome_meta.chr_name`              |
| copies_M    | UTINYINT | Copies carried by males (default 2)                      |
| copies_F    | UTINYINT | Copies carried by females (default 2)                    |
| hemi_parent | VARCHAR  | For single-copy chromosomes: `"parent_1"`/`"parent_2"`; NULL for diploids |
| recombines  | BOOLEAN  | TRUE if recombination occurs (default TRUE)              |

### `founder_haplotypes`

Founder haplotype pool in **long** format. Created by
`define_founder_haplotypes()`, sampled by `add_founders()`.

| Column       | Type     | Notes                                                    |
|--------------|----------|----------------------------------------------------------|
| line_name    | VARCHAR  | NULL = shared pool; set for line-specific pools; logical key part |
| haplotype_id | INTEGER  | Sequential within the pool (unique per `line_name`); logical key part |
| locus_name   | VARCHAR  | FK to `genome_meta.locus_name`; logical key part         |
| allele       | INTEGER  | 0 or 1                                                   |

Logical key `(line_name, haplotype_id, locus_name)` enforced in R (nullable
`line_name`), matching the `index_meta` convention.

### `ind_meta`

Individual-level metadata. Created empty by `initialize_genome()`; rows
populated by `add_founders()` and `add_offspring()`.

| Column      | Type    | Notes                               |
|-------------|---------|-------------------------------------|
| id_ind      | VARCHAR | Primary key, format `{line_name}_{n}` (e.g. `Libra_1020`) |
| id_parent_1 | VARCHAR | NA for founders                     |
| id_parent_2 | VARCHAR | NA for founders                     |
| line_name   | VARCHAR | Genetic line name                   |
| sex         | VARCHAR | "M" or "F"                          |
| *user cols* | any     | Added via `mutate_table()` or `...` in `add_founders()` |

**Reserved**: `id_ind`, `id_parent_1`, `id_parent_2`, `line_name`, `sex`

### `trait_meta`

One row per **genetic component trait**. Populated by `define_trait()`.
Contains only genetic-layer information — no phenotype-level metadata.
Observation-layer metadata lives in `phenotype_meta`.

| Column          | Type    | Notes                                                              |
|-----------------|---------|--------------------------------------------------------------------|
| id_trait        | INTEGER | Primary key (auto-incrementing)                                    |
| trait_name      | VARCHAR | Unique identifier; equals `phenotype_name` for simple traits       |
| description     | VARCHAR | Free text                                                          |
| units           | VARCHAR | e.g. `"kg"`, `"g/day"`                                             |
| expressed_parent| VARCHAR | `"both"` (default), `"parent_1"` (paternal), `"parent_2"` (maternal) — imprinting |
| target_add_mean | DOUBLE  | TBV centering mean for the base population; default `0`            |

**What does NOT belong here** (all moved to `phenotype_meta` in v0.31.0):
`type`, `expressed_sex`, `repeatable`, `mean`, `min_value`, `max_value`,
`prevalence`, `thresholds`, `cat_values`, `cat_names`, `residual_var`,
`index_weight`, `economic_value`.

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

### `trait_var_comp`

Genetic-layer variance components. One row per (effect_name, trait_name_1, trait_name_2).
Both `(i,j)` and `(j,i)` pairs stored. Populated by `define_effect_cov_matrix()` and
`define_trait()`. Stores **only** genetic effects — no phenotype-level variances.

Valid `effect_name` values: `"gen_add"` (additive genetic G matrix);
future: `"dominance"`, `"epistasis"`. Named random effects (HYS, litter, pen)
go to `phenotype_var_comp`, not here.

| Column           | Type    | Notes                                              |
|------------------|---------|----------------------------------------------------|
| id_trait_var_comp| INTEGER | Primary key (auto-incrementing)                    |
| effect_name      | VARCHAR | `"gen_add"`; future: `"dominance"`, `"epistasis"`  |
| trait_name_1     | VARCHAR |                                                    |
| trait_name_2     | VARCHAR |                                                    |
| cov_value        | DOUBLE  | Variance (diagonal) or covariance (off-diagonal)   |

### `phenotype_meta`

Observed phenotype definitions. One row per phenotype name. Populated by
`define_phenotype()`. Analogous to `trait_meta` but for the observation layer —
simple traits have the same name in both tables; composite phenotypes (e.g. WW,
SGE ADG) appear only here.

| Column                   | Type    | Notes                                                         |
|--------------------------|---------|---------------------------------------------------------------|
| id_phenotype_meta        | INTEGER | Primary key (auto-incrementing)                               |
| phenotype_name           | VARCHAR | Unique. Equals `trait_name` for simple traits.                |
| type                     | VARCHAR | `"continuous"`, `"count"`, `"categorical"`, `"derived_formula"` |
| mean                     | DOUBLE  | Phenotypic population mean / liability intercept              |
| expressed_sex            | VARCHAR | `"both"`, `"M"`, or `"F"`                                     |
| repeatable               | BOOLEAN | Repeated records allowed?                                     |
| min_value / max_value    | DOUBLE  | Clipping bounds for count traits                              |
| prevalence               | DOUBLE  | For 2-category categorical traits                             |
| thresholds               | VARCHAR | Comma-separated liability cutpoints for K-category traits     |
| cat_values               | VARCHAR | Comma-separated numeric phenotype values per category         |
| cat_names                | VARCHAR | Comma-separated labels per category                           |
| store_liability          | BOOLEAN | Write raw liability to `ind_phenotype.liability_value`        |
| missing_component_action | VARCHAR | `"skip"` (default) or `"error"` — what to do when any component of a composite phenotype cannot be resolved for an individual |

**Reserved**: all columns (managed by `define_phenotype()`).

### `phenotype_components`

Component definitions for composite phenotypes. One row per (phenotype ×
component). Populated by `define_phenotype(..., components = ...)`. Simple
(non-composite) phenotypes have no rows here.

| Column             | Type    | Notes                                                              |
|--------------------|---------|--------------------------------------------------------------------|
| id_phenotype_comp  | INTEGER | Primary key (auto-incrementing)                                    |
| phenotype_name     | VARCHAR | FK to `phenotype_meta.phenotype_name`                              |
| source_trait_name  | VARCHAR | FK to `trait_meta.trait_name` — the genetic component trait        |
| contributor_type   | VARCHAR | `"self"`, `"dam"`, `"sire"`, or `"group"`                          |
| group_column       | VARCHAR | Column in `group_table` that holds group membership (required for `"group"`) |
| group_table        | VARCHAR | Table containing `group_column`; default `"ind_meta"`              |
| aggregation        | VARCHAR | `"sum"` (default) or `"mean"` — how group-mates' TBVs are combined |
| weight             | DOUBLE  | Scalar multiplier; default `1.0`                                   |
| weight_type        | VARCHAR | `"fixed"` (default), `"covariate"`, `"legendre"`, `"raw_poly"`     |
| covariate_name     | VARCHAR | Covariate column when `weight_type = "covariate"`                  |
| covariate_table    | VARCHAR | Table containing `covariate_name`                                  |
| poly_order         | INTEGER | Polynomial basis order                                             |
| poly_scale_min/max | DOUBLE  | Legendre scaling bounds                                            |
| genome_effect_types| VARCHAR | Default `"additive"`                                               |
| missing_action     | VARCHAR | Per-component fallback (currently unused; use `phenotype_meta.missing_component_action`) |
| contributor_filter | VARCHAR | Reserved for spatial/neighborhood lookup — not yet implemented     |

**Note on SGE (Social Genetic Effects / Bijma model)**: for `contributor_type = "group"`,
`add_phenotype()` aggregates group-mates' TBVs (excluding self). A singleton (no
group-mates) receives a social contribution of 0 and is not excluded. An individual
with no group assignment receives `NA` and is handled by `missing_component_action`.

### `phenotype_var_comp`

Phenotype-layer variance components. One row per (effect_name, phenotype₁, phenotype₂,
condition). Both (i,j) and (j,i) pairs stored for off-diagonal entries. Populated by
`define_phenotype(..., residual_var = ...)`, `define_residual_cov()`, and
`define_effect_random()`.

`effect_name = 'residual'` is reserved for residual noise. All named random effects
(HYS, litter, pen, etc.) use their own string (e.g. `'hys'`, `'litter'`).
The `condition_column` / `condition_level` columns are used only for `'residual'`
to model heterogeneous residual variance by sex, group, etc.

| Column               | Type    | Notes                                                              |
|----------------------|---------|--------------------------------------------------------------------|
| id_phenotype_var_comp| INTEGER | Primary key (auto-incrementing)                                    |
| effect_name          | VARCHAR | `'residual'` or any named random effect (e.g. `'hys'`, `'litter'`)|
| phenotype_name_1     | VARCHAR |                                                                    |
| phenotype_name_2     | VARCHAR |                                                                    |
| cov_value            | DOUBLE  |                                                                    |
| condition_column     | VARCHAR | NULL = unconditional; used only for `effect_name = 'residual'`     |
| condition_table      | VARCHAR | Default `"ind_meta"`                                               |
| condition_level      | VARCHAR | Value of `condition_column` for this row                           |
| weight_type          | VARCHAR | Default `"fixed"`                                                  |
| poly_order           | INTEGER | Polynomial order for `"legendre"` weight type                      |

### `ind_phenotype`

Phenotype records in long format. Populated by `add_phenotype()`.

| Column         | Type    | Notes                                             |
|----------------|---------|---------------------------------------------------|
| id_phenotype   | INTEGER | Primary key (global auto-incrementing integer)    |
| id_ind         | VARCHAR |                                                   |
| phenotype_name | VARCHAR | FK to `phenotype_meta.phenotype_name`             |
| pheno_value    | DOUBLE  | Phenotype value                                   |
| pheno_number   | INTEGER | 1 = first record for this individual × trait, etc.|
| *user cols*    | any     | Added via `mutate_table()` or scalar `...` in `add_phenotype()` |

### `ind_tbv`

True breeding values (simulation ground truth). Populated by
`add_phenotype()` and `add_tbv()`. Logical key `(id_ind, trait_name)` unique.

| Column     | Type    | Notes                                  |
|------------|---------|----------------------------------------|
| id_tbv     | INTEGER | Primary key (auto-incrementing)        |
| id_ind     | VARCHAR |                                        |
| trait_name | VARCHAR |                                        |
| tbv_value  | DOUBLE  |                                        |

### `ind_ebv`

Estimated breeding values from external BLUP / GBLUP runs. Logical key
`(id_ind, trait_name, model, eval_number)` unique. Populated by `add_ebv()`.

| Column      | Type    | Notes                                                   |
|-------------|---------|--------------------------------------------------------------|
| id_ebv      | INTEGER | Primary key (auto-incrementing)                              |
| id_ind      | VARCHAR |                                                              |
| trait_name  | VARCHAR |                                                              |
| model       | VARCHAR | User label, e.g. "ssGBLUP_v1"                               |
| ebv_value   | DOUBLE  |                                                              |
| acc         | DOUBLE  | Optional accuracy                                            |
| se          | DOUBLE  | Optional standard error                                      |
| eval_number | INTEGER | Auto-incrementing counter per trait (global across models); 1 = first evaluation |

### `index_meta`

Selection index definitions. One row per (index × trait). Populated by
`define_index()`. A special row with `index_name = NULL` is written by
`define_trait()` to record the global economic weight for each trait.

| Column          | Type    | Notes                                                         |
|-----------------|---------|---------------------------------------------------------------|
| id_index_name   | INTEGER | Primary key (auto-incrementing)                               |
| index_name      | VARCHAR | NULL = global/default entry written by `define_trait()`; named index written by `define_index()` |
| trait_name      | VARCHAR | FK to `trait_meta.trait_name`                                 |
| index_weight    | DOUBLE  | Selection index weight (NULL for global rows)                 |
| economic_weight | DOUBLE  | Economic value per unit of the trait                          |
| *user cols*     | any     | Added via `...` in `define_index()`                           |

**Reserved**: `id_index_name`, `index_name`, `trait_name`, `index_weight`, `economic_weight`

**Unique constraint**: `(index_name, trait_name)` (NULL `index_name` uniqueness enforced in R code).

### `ind_index`

Computed selection index values. One row per (individual × index × run).
Populated by `add_index()`.

| Column       | Type    | Notes                                          |
|--------------|---------|------------------------------------------------|
| id_index     | INTEGER | Primary key (auto-incrementing)                |
| id_ind       | VARCHAR |                                                |
| index_name   | VARCHAR | FK to `index_meta.index_name`                  |
| index_number | INTEGER | Auto-incrementing run counter per individual   |
| index_value  | DOUBLE  | Computed index value                           |
| *user cols*  | any     | Added via `...` in `add_index()`               |

### `ind_true_index`

True selection index values computed from TBVs (simulation ground truth).
Populated by `add_tbv()` when `index_names` is supplied.

| Column           | Type    | Notes                                                        |
|------------------|---------|--------------------------------------------------------------|
| id_true_index    | INTEGER | Primary key (auto-incrementing)                              |
| id_ind           | VARCHAR |                                                              |
| index_name       | VARCHAR | FK to `index_meta.index_name`                                |
| weight_type      | VARCHAR | `"index"` (uses `index_weight`) or `"economic"` (uses `economic_weight`) |
| true_index_value | DOUBLE  | Weighted sum: `sum(weight_i * tbv_i)` across index traits    |

**Reserved**: all columns (managed exclusively by `add_tbv()` when `index_names` is supplied).

Logical row key `(id_ind, index_name, weight_type)` — one true index value per
individual × index × weight type. No SQL `UNIQUE` constraint; uniqueness enforced
in R via DELETE + INSERT when `overwrite_index = TRUE`.

## Implemented Functions (Phase 1 Complete)

### `initialize_genome()`

`R/initialize_genome.R`

Creates the DuckDB file (or in-memory DB) and populates **all** core tables
eagerly, including genome tables and all individual/trait tables:

- Genome: `genome_meta`, `ind_haplotype` (empty), `ind_genotype` (empty), `chr_meta` (default autosome rows)
- Optionally: `founder_haplotypes` (if `n_haplotypes` provided)
- Individual/trait (all empty): `ind_meta`, `ind_phenotype`, `ind_tbv`,
  `ind_ebv`, `trait_meta`, `trait_effects`, `trait_var_comp`,
  `trait_random_effects`

All tables are registered in `pop$tables` immediately. Users can call
`get_table()` and `mutate_table()` on any table right after init, before any
data is inserted.

Key params: `pop_name`, `n_loci`, `n_chr`, `chr_len_Mb`, `db_path`

### `add_founders()`

`R/add_founders.R`

Samples haplotypes for each founder individual from the `founder_haplotypes`
pool. Appends rows to `ind_meta` (core 5 cols) and `ind_haplotype`
(long: 2 haplotypes × n_loci rows each, with `line_origin` = the founder's line
and `strand = 1`). Does **not** write `ind_genotype` (on-demand via
`add_dosage()`). ID format: `{line_name}_{n}` (e.g. `Libra_1`).

Accepts `...` for custom `ind_meta` columns written atomically with the new
rows (see **Custom field forwarding** below).

Key params: `n_males`, `n_females`, `line_name`, then `...` for custom fields.

### Custom field forwarding in `add_*` functions

`add_founders()`, `add_phenotype()`, `add_tbv()`, and `add_ebv()` all accept
`...` for optional custom columns written to the target table at the same time
the new rows are inserted. This avoids a redundant second `mutate_table()` step.

```r
# Single call — gen and farm written with the founders
pop <- pop |>
  add_founders(n_males = 10, n_females = 100, line_name = "A",
               gen = 0L, farm = "Iowa")

# add_phenotype: scalar only (broadcast to all phenotype records)
pop <- pop |>
  get_table("ind_meta") |>
  add_phenotype("ADG", test_env = "barn_A")
```

**Argument disambiguation**: R's standard matching routes explicit formal params
(`n_males`, `n_females`, `line_name`, etc.) to their positions; anything else
falls into `...` and is treated as a custom column. Reserved column names
(`id_ind`, `sex`, `line`, etc.) are blocked with an error.

**Type safety** — column types are inferred from the R value via
`infer_duckdb_type()`. Use R's type suffixes to get the right DuckDB type:

| R value | DuckDB type |
|---------|-------------|
| `0L`, `NA_integer_` | `INTEGER` |
| `0`, `0.0`, `NA_real_` | `DOUBLE` |
| `TRUE`/`FALSE`, `NA` (bare) | `BOOLEAN` |
| `"text"`, `NA_character_` | `VARCHAR` |
| `as.Date(...)` | `DATE` |
| `Sys.time()`, `as.POSIXct(...)` | `TIMESTAMP` |

Common pitfall: `gen = 0` gives DOUBLE, not INTEGER. Use `gen = 0L`.

**Pre-declaring a column schema** before data exists (typed-NA workflow):

```r
# 1. After initialize_genome(), declare column types on the empty table
pop <- pop |>
  get_table("ind_meta") |>
  mutate_table(gen = NA_integer_, farm = NA_character_)

# 2. add_founders fills values in; types already match
pop <- pop |>
  add_founders(n_males = 10, n_females = 100, line_name = "A",
               gen = 0L, farm = "Iowa")
```

The shared internal helper `prepare_extra_cols()` in `R/sql_utils.R` handles
validation, type inference, ALTER TABLE, and scalar broadcast for all `add_*`
functions.

**Scalar vs. vector**:
- `add_founders()` / `add_offspring()`: scalars broadcast; vectors must have
  length `n_males + n_females` (or `n_offspring`).
- `add_phenotype()` / `add_tbv()` / `add_ebv()`: scalar only (row count per
  trait varies). Use `mutate_table()` afterwards for per-record vectors.

### `mutate_table()`

`R/mutate_table.R`

Generic function. Adds or updates columns in **any** database table. Chain after
`get_table()` (and optionally `filter()`) then call `mutate_table(col = value)`.
Scalar values are broadcast to all (or filtered) rows; vectors must match the
effective row count. Type is inferred via `infer_duckdb_type()`. Reserved
columns are blocked via `TABLE_RESERVED_COLS` in `R/sql_utils.R`. Returns `pop`
invisibly.

When called on an **empty table**, `mutate_table()` still creates the column
schema via `ALTER TABLE ADD COLUMN` (no rows are updated). This is the mechanism
for pre-declaring typed column schemas before data arrives.

```r
# All rows
pop <- pop |> get_table("ind_meta") |> mutate_table(gen = 1L)

# Filtered rows (unmatched rows get NULL for new columns)
pop <- pop |>
  get_table("ind_meta") |>
  filter(sex == "M") |>
  mutate_table(gen = 2L)

# Pre-declare schema on empty table
pop <- pop |>
  get_table("ind_ebv") |>
  mutate_table(model_version = NA_character_)
```

### `define_chip()`

`R/define_chip.R`

Marks the loci in a filtered `genome_meta` table as members of a named chip,
writing a `BOOLEAN` column `is_{chip_name}` to `genome_meta`. All other loci
receive `FALSE`. Accepts a `tidybreed_table` from `get_table("genome_meta")`
(optionally filtered) as its first argument; returns `tidybreed_pop`.

```r
# Filter first, then define chip
pop |> get_table("genome_meta") |> filter(chr %in% 1:5) |> define_chip("chr1to5")

# Random chip: sample locus names, then filter
sel <- pop |> get_table("genome_meta") |> collect() |>
  slice_sample(n = 500) |> pull(locus_name)
pop |> get_table("genome_meta") |> filter(locus_name %in% sel) |> define_chip("50K")

# Complement of an existing chip
pop |> get_table("genome_meta") |> filter(is_50K == FALSE) |> define_chip("non50K")
```

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
`ind_meta`, `ind_phenotype`, `ind_haplotype`, `ind_genotype`).

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

### `define_trait()` / `define_additive_effects()`

`R/define_trait.R`, `R/define_additive_effects.R`

- `define_trait()` — **genetic layer only**. Writes one row to `trait_meta` and
  a global `(index_name = NULL, trait_name, economic_weight = 0)` row to
  `index_meta`. Accepted arguments: `trait_name`, `target_add_var` (writes to
  `trait_var_comp`), `target_add_mean`, `expressed_parent`, `description`,
  `units`, `overwrite`. **Never** pass observation-layer arguments here
  (`type`, `mean`, `expressed_sex`, `residual_var`, etc.) — those belong
  in `define_phenotype()`. `overwrite = FALSE` (default) errors if the trait
  already exists; `overwrite = TRUE` replaces both the `trait_meta` row and its
  `index_meta` entry.
- `define_additive_effects()` — accepts a `tidybreed_table` from
  `get_table("genome_meta")` (optionally filtered) as its **first argument**.
  `trait_name` accepts a scalar **or vector** of trait names:
  - **Single trait** — manual (`effects` vector) or sampled (`distribution =
    "normal"/"gamma"`) with optional Falconer rescale. Reads `target_add_var`
    from `trait_var_comp`.
  - **Multiple traits** — draws correlated effects from `MVN(0, G)` via
    `MASS::mvrnorm`. `G = NULL` reads from `trait_var_comp`. `method =
    "shared"` (all traits use the filtered loci) or `"union"` (per-trait QTL
    sets from existing `genome_effects` rows, restricted to the filtered loci).

  Re-calling for the same trait replaces existing rows in `genome_effects`.

  ```r
  # Single trait
  pop |> get_table("genome_meta") |> filter(chr %in% 1:5) |> define_additive_effects("ADG")

  # Multiple correlated traits (shared QTL set)
  G <- matrix(c(0.25, 0.10, 0.10, 0.30), 2, 2,
              dimnames = list(c("ADG", "BW"), c("ADG", "BW")))
  pop |> get_table("genome_meta") |> filter(chr %in% 1:5) |>
    define_additive_effects(c("ADG", "BW"), G = G)

  # Use generation-0 animals to define base allele frequencies
  gen0 <- get_table(pop, "ind_meta") |> filter(gen == 0L)
  pop |> get_table("genome_meta") |> filter(...) |>
    define_additive_effects("ADG", base = "current_pop", base_tbl = gen0)
  ```

### `define_effect_cov_matrix()` / `define_effect_random()` / `define_effect_fixed_class()` / `define_effect_fixed_cov()` / `define_effect_int()`

`R/define_effect_cov_matrix.R`, `R/define_effect_random.R`, `R/define_effect_fixed_class.R`, `R/define_effect_fixed_cov.R`, `R/define_effect_int.R`

- `define_effect_cov_matrix(pop, effect_name, cov_matrix)` — **single entry
  point for all variance/covariance data**. Routes by `effect_name`:
  genetic effects (`"gen_add"`, `"dominance"`, `"epistasis"`) → `trait_var_comp`;
  `"residual"` → `define_residual_cov()` → `phenotype_var_comp`;
  any other name → `phenotype_var_comp` with that `effect_name`.
  Can be called before `define_trait()` or `define_effect_random()`.
- `define_effect_random()` — `variance` optional if already in `phenotype_var_comp`.
- `define_effect_fixed_class()` — discrete level → shift mapping.
- `define_effect_fixed_cov()` — linear regression term (`slope * (x - center)`).
- `define_effect_int()` — sets intercept (`target_add_mean`) for a trait.

### `define_phenotype()` / `define_residual_cov()`

`R/define_phenotype.R`, `R/define_residual_cov.R`

- `define_phenotype(pop, phenotype_name, type, mean, expressed_sex, repeatable, ...)` —
  registers an observed phenotype in `phenotype_meta`. For simple traits
  `phenotype_name` matches the `trait_name` already in `trait_meta`. For
  composite phenotypes (e.g. weaning weight, SGE ADG) the name is new and no
  prior `define_trait()` call is needed for it.

  Key arguments:
  - `residual_var` — scalar; writes one unconditional diagonal entry to
    `phenotype_var_comp` (with `effect_name = 'residual'`). For correlated or
    heterogeneous residuals use `define_residual_cov()` afterwards.
  - `components` — data frame with columns `source_trait_name` and
    `contributor_type` (`"self"`, `"dam"`, `"sire"`, `"group"`). Optional
    columns: `weight`, `weight_type`, `aggregation`, `group_column`,
    `group_table`, `covariate_name`, etc. Writes to `phenotype_components`.
    `NULL` (default) = simple single-self trait.
  - `missing_component_action` — `"skip"` (default) or `"error"`. Stored in
    `phenotype_meta` and applied uniformly by `add_phenotype()` for **any**
    missing composite piece (missing group assignment, missing dam/sire TBV,
    etc.). `"skip"` excludes the individual and warns with a count + up to 5
    example IDs. `"error"` stops immediately.

- `define_residual_cov(pop, phenotype_names, cov_matrix, condition_column = NULL, ...)` —
  writes conditional or unconditional residual (co)variance entries to
  `phenotype_var_comp` (always with `effect_name = 'residual'`). Supply a named
  matrix for multi-phenotype correlated residuals, or call once per sex/group
  level with `condition_column = "sex"` and `condition_level = "M"` / `"F"` for
  heterogeneous residuals.

### `add_trait_covariate()` *(deprecated)*

`R/add_trait_covariate.R`

- Deprecated since v0.6.0. Use `define_effect_fixed_class()`, `define_effect_fixed_cov()`,
  or `define_effect_random()` instead.

### `add_phenotype()` / `add_tbv()`

`R/add_phenotype.R`, `R/add_tbv.R`

Both functions accept a `tidybreed_table` (from `get_table()` + optional
`filter()`) as their first argument and return `tidybreed_pop`.

- `add_phenotype()` — the workhorse. `phenotype_name` (formerly `trait_name`)
  defaults to all phenotypes in `phenotype_meta` when omitted. Internally calls
  `add_tbv()` first for all required source traits (including all composite
  components and group-member IDs), then assembles the composite TBV via
  `.assemble_composite_tbv()`. For group contributors (SGE model), all
  group-member TBVs are pre-fetched so group aggregation is a single in-memory
  pass. Adds fixed/random covariate contributions, samples residuals (joint
  `MVN(0, R)` when multiple phenotypes share the subset and `R` is stored;
  otherwise independent). Converts liability to phenotype per `type`.
  Writes `ind_phenotype` rows and updates `ind_tbv`.
- `add_tbv()` — TBV-only; no phenotype records. Computes centered TBV by joining
  `ind_haplotype` to `genome_effects` (`genome_effect_type = "additive"`) on
  `(locus_name, line_origin)`: a line-specific effect row (`genome_effects.line_name
  = ind_haplotype.line_origin`) is preferred, falling back per-locus to the
  population-wide row (`line_name IS NULL`) only when no line-specific row exists
  for that locus/line. This is what makes crossbreeding TBV correct (e.g. a Duroc ×
  Landrace F1 centered against each parent line's own QTL effects and base allele
  frequency). Imprinted traits (`expressed_parent = "parent_1"`/`"parent_2"`)
  restrict the join to that parent's `parent_origin` before the same line-matching
  logic applies. `trait_name` also defaults to all traits in `trait_meta` when
  omitted. Optional arguments for true index computation:
  - `index_names` — character vector of named indices; when supplied, multiplies
    per-trait TBVs by the index weights and writes results to `ind_true_index`.
    `NULL` (default) skips index computation.
  - `type` — `"index"` (default, uses `index_weight`), `"economic"` (uses
    `economic_weight`), or `"both"` (writes two rows per individual distinguished
    by `weight_type`).
  - `overwrite_index = FALSE` — when `FALSE`, skips individuals that already have
    a value in `ind_true_index` for the given `(index_name, weight_type)`. Set
    `TRUE` to recompute (e.g. after updating index weights).

### `define_trait_simple()`

`R/define_trait_simple.R`

Convenience wrapper that chains `define_trait()` and `define_additive_effects()` for
a single uncorrelated trait. QTL are always placed randomly (n = `n_qtl`).
For non-random QTL placement or correlated multi-trait effects, use the
functions individually with `get_table("genome_meta") |> filter(...)  |>
define_additive_effects()`.

### `define_index()` / `add_index()`

`R/define_index.R`, `R/add_index.R`

- `define_index(pop, index_name, trait_names, index_wts, economic_wts = NULL, overwrite = FALSE, ...)` —
  registers a named selection index in `index_meta`. `overwrite = FALSE` (default)
  is a no-op when `(index_name, trait_name)` already exists; `overwrite = TRUE`
  updates weights and economic weights in place. `economic_wts` is an optional
  numeric vector (same length as `trait_names`; some values may be 0) written to
  `index_meta.economic_weight`. Extra `...` columns are broadcast or per-trait.
- `add_index(tbl, index_name, value_col = NULL, overwrite_index = FALSE, delete_all = FALSE, ...)` —
  accepts a `tidybreed_table` from `get_table()` (optionally filtered). Any table
  with `id_ind`, `trait_name`, and a numeric value column is accepted: `ind_ebv`,
  `ind_phenotype`, `ind_tbv`, or a user-defined table. `value_col` is auto-detected
  from the table name (`ind_ebv` → `"ebv_value"`, `ind_phenotype` → `"pheno_value"`,
  `ind_tbv` → `"tbv_value"`); supply it explicitly for unknown tables.
  Multiplies each individual's values by the index weights in `index_meta` and
  appends to `ind_index`. Every individual must have exactly one value per index
  trait — an error is thrown if duplicates are found (filter to a single model /
  `eval_number` / `pheno_number` first). Issues a warning when no filter is applied.
  `overwrite_index = TRUE` clears prior runs for the named index; `delete_all = TRUE`
  clears all of `ind_index`.

## Roadmap

### Longer-Term

- `select_parents()` — selection index or truncation selection
- Export: PLINK `.bed/.bim/.fam`, VCF
- Visualization helpers
- Dominance and epistasis effects (currently only additive)

## License

`tidybreed` is licensed under the **MIT License**. Use `License: MIT + file LICENSE`
in `DESCRIPTION`. The `LICENSE` file must contain the standard MIT license text with
the year and copyright holder.

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





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
4. **Pipe-friendly** — every exported function accepts a `tidybreed_pop` and
   returns a `tidybreed_pop`
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
| ind_id        | VARCHAR |
| parent_origin | INTEGER |
| locus_1 … locus_n | INTEGER (0 or 1) |

### `genome_genotype`

Genotypes in 0/1/2 encoding. One row per individual.

| Column        | Type    |
|---------------|---------|
| ind_id        | VARCHAR |
| locus_1 … locus_n | INTEGER (0, 1, or 2) |

Both haplotypes and genotypes are stored (not computed on demand) because disk
is cheap and recomputation during mating is expensive.

### `ind_meta`

Individual-level metadata. Created by `add_founders()`.

| Column   | Type    | Notes                          |
|----------|---------|--------------------------------|
| ind_id   | VARCHAR | Primary key, format `{line}-{n}` |
| parent_1 | VARCHAR | "0" for founders               |
| parent_2 | VARCHAR | "0" for founders               |
| line     | VARCHAR | Genetic line name              |
| sex      | VARCHAR | "M" or "F"                     |
| *user cols* | any  | Added via `mutate_ind_meta()`  |

**Reserved**: `ind_id`, `parent_1`, `parent_2`, `line`, `sex`

### `ind_phenotype`

Phenotypic records in long format. Populated by `add_phenotype()` (not yet implemented).

We need to think through this table more before implementing. 

| Column        | Type    | Notes          |
|---------------|---------|----------------|
| ind_id        | VARCHAR |                |
| trait         | VARCHAR | Trait name     |
| value         | DOUBLE  | Phenotypic value |
| env           | VARCHAR | Optional       |
| rep           | INTEGER | Optional       |
| date_measured | DATE    | Optional       |

### `trait_meta`

Trait genetic architecture. Populated by `add_trait()` (not yet implemented).

We need to think through this function and table before implementing. 

| Column      | Type    | Notes                                    |
|-------------|---------|------------------------------------------|
| trait_name  | VARCHAR | e.g. "ADG"                               |
| effect_name | VARCHAR | e.g. "QTL_1", "residual"                 |
| effect_type | VARCHAR | "QTL", "polygenic", "environmental"      |
| level       | VARCHAR | "additive", "dominance", "epistasis"     |
| locus_id    | INTEGER | Associated locus (NULL for residual)     |
| true_effect | DOUBLE  | True effect size                         |
| distribution| VARCHAR | "normal", "uniform", …                   |
| mean        | DOUBLE  | Distribution mean in founder pop         |
| variance    | DOUBLE  | Distribution variance in founder pop     |

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

### `mutate_ind_meta()`

`R/mutate_ind_meta.R`

Adds or updates user columns in `ind_meta`. Scalar values are broadcast to all
rows; vectors assign per-individual. Type is inferred via
`infer_duckdb_type()`. Reserved columns are blocked.

### `mutate_genome_meta()`

`R/mutate_genome_meta.R`

Same pattern as `mutate_ind_meta()` but targets `genome_meta`. Used internally
by `define_chip()` to write chip membership columns.

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

`get_table(pop, "table_name")` returns a lazy `dplyr` tibble (`tbl()`).
`close_pop()` safely closes the DuckDB connection.

## Roadmap

### Next: Trait Architecture

**`add_trait()`** — creates/populates `trait_meta`:

- User specifies trait name, number of QTLs, additive variance target,
  residual variance
- Function samples QTL loci from `genome_meta`, draws effect sizes scaled to
  hit `target_add_var`, writes rows to `trait_meta`
- Should also write a `BOOLEAN` column `is_QTL_{trait}` to `genome_meta` for
  easy filtering

### Next: Phenotype Simulation

**`add_phenotype(trait)`** — populates `ind_phenotype`:

- Reads `trait_meta` to get QTL loci and effects
- Reads `genome_genotype` to compute true breeding value (TBV) per individual
- Adds residual drawn from N(0, res_var)
- Writes `ind_id`, `trait`, `value` rows to `ind_phenotype`
- Should accept a filtered subset of individuals (via pipe) so users can
  phenotype only e.g. females in a given generation

### `filter()` Integration

Users should be able to pipe `dplyr::filter()` on `ind_meta` before calling
`add_phenotype()` or genotyping functions to restrict which animals are
processed:

```r
pop |>
  filter(sex == "F", gen == 1L) |>
  add_phenotype(trait = "ADG")
```

### Longer-Term

- `mate()` — recombination engine, produces offspring haplotypes/genotypes
- `select_parents()` — selection index or truncation selection
- Export: PLINK `.bed/.bim/.fam`, VCF
- Visualization helpers

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





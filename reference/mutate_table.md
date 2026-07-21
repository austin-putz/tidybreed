# Add or modify columns in any population database table

Generic replacement for the table-specific `mutate_ind_meta()` /
`mutate_genome_meta()` family. Chain after
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
(and optionally
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html))
to add new columns or update existing ones. `value` can be a scalar
broadcast to all (or filtered) rows, or a tibble keyed on the table's
primary key for per-row values — see `@examples` below for both, plus
the schema-pre-declaration pattern for empty tables.

## Usage

``` r
mutate_table(tbl_obj, ..., .set_default = FALSE)
```

## Arguments

- tbl_obj:

  A `tidybreed_table` object returned by
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md).

- ...:

  Named arguments of the form `column_name = value`. `value` can be a
  scalar (applied to all affected rows) or a tibble/data frame
  containing the primary key column and a column named after the field.
  Plain vectors (length \> 1) are not allowed; use a tibble instead.

- .set_default:

  Logical; if `TRUE`, creates new columns with a SQL DEFAULT constraint
  set to the provided value. The DEFAULT applies to future INSERT
  operations (e.g., from
  [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md),
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md))
  when the column is not explicitly specified in `...`. Only valid for
  scalar values; an error is raised if used with tibbles. Has no effect
  on columns that already exist (DuckDB does not support modifying
  existing column defaults). Default: `FALSE`.

## Value

The parent `tidybreed_pop` object (invisibly).

## Details

**Type inference** (`logical` → BOOLEAN, `integer` → INTEGER, `numeric`
→ DOUBLE, `Date` → DATE, `POSIXct` → TIMESTAMP, `character` → VARCHAR)
is performed automatically via
[`infer_duckdb_type()`](https://austin-putz.github.io/tidybreed/reference/infer_duckdb_type.md).

**Per-row values (tibble path)**: to set different values for different
rows, pass a tibble/data frame containing the primary key column plus a
column named after the field being updated:

    mutate_table(gen = tibble::tibble(id_ind = c("A_1", "A_2"),
                                       gen = c(1L, 2L)))

The tibble is entirely self-contained; any upstream
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
is ignored for that field (but continues to apply to scalar fields in
the same call). All IDs in the tibble must exist in the table; an error
is raised for unknown IDs.

**Filtering**: when a
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
is applied upstream, scalar values broadcast to all matching rows. New
columns created in this context will be `NULL` for all non-matching
rows.

**Default constraints**: When `.set_default = TRUE`, new columns are
created with a SQL DEFAULT constraint. This means future INSERT
operations from
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md),
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
or other `add_*()` functions will automatically use the default value
when the column is not explicitly provided in `...`.

For populated tables, existing rows are still updated via the standard
UPDATE mechanism; the DEFAULT only affects subsequent INSERT operations.
For empty tables (e.g., right after
[`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md)
/
[`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md),
before any founders are added), the DEFAULT is set at the schema level
with no data operations.

Note: DEFAULT constraints can only be added when creating new columns.
If a column already exists, `.set_default` is ignored (but the UPDATE
proceeds normally).

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- open_pop(pop_name = "sim", db_name = ":memory:") |>
  define_genome(n_loci = 100, n_chr = 5, chr_len_Mb = 50) |>
  define_founder_haplotypes(n_haplotypes = 100, line_name = "A")
pop <- pop |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 10, n_females = 100, line_name = "A")

# Broadcast: assign gen = 1 to every individual
pop <- pop |>
  get_table("ind_meta") |>
  mutate_table(gen = 1L)

# Filtered subset: bump gen only for males
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(sex == "M") |>
  mutate_table(gen = 2L)

# Per-row values via a tibble keyed on the primary key column
pop <- pop |>
  get_table("ind_meta") |>
  mutate_table(
    gen = tibble::tibble(
      id_ind = c("A_1", "A_2"),
      gen    = c(5L, 6L)
    )
  )

# Pre-declare a typed column schema before any data exists
pop2 <- open_pop(pop_name = "sim2", db_name = ":memory:") |>
  define_genome(n_loci = 100, n_chr = 5, chr_len_Mb = 50)
pop2 <- pop2 |>
  get_table("ind_meta") |>
  mutate_table(farm = NA_character_)

# Pre-declare WITH a SQL DEFAULT so future add_founders() rows auto-fill
pop3 <- open_pop(pop_name = "sim3", db_name = ":memory:") |>
  define_genome(n_loci = 100, n_chr = 5, chr_len_Mb = 50) |>
  define_founder_haplotypes(n_haplotypes = 100, line_name = "A")
pop3 <- pop3 |>
  get_table("ind_meta") |>
  mutate_table(gen = 0L, active = TRUE, .set_default = TRUE)

# Future founders automatically get gen = 0L and active = TRUE
pop3 <- pop3 |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 10, n_females = 100, line_name = "A")
pop3 |> get_table("ind_meta") |> dplyr::collect()  # all have defaults

# Override the default by supplying an explicit value
pop3 <- pop3 |>
  get_table("founder_haplotypes") |>
  dplyr::filter(line_name == "A") |>
  add_founders(n_males = 5, n_females = 50, line_name = "B", gen = 1L)

# Line B: gen = 1L (explicit), active = TRUE (default)
pop3 |>
  get_table("ind_meta") |>
  dplyr::filter(line_name == "B") |>
  dplyr::collect()
} # }
```

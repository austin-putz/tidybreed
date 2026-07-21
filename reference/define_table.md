# Define a custom table in the population database

Creates a new user-defined table inside the tidybreed DuckDB database
and registers it in `pop$tables` so that
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
and
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
work on it immediately.

Column names and types are specified as named `...` arguments using the
same typed-NA convention as
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
schema pre-declaration:

    pop |> define_table(
      "sim_timing",
      run_id       = NA_integer_,
      duration_sec = NA_real_,
      label        = NA_character_
    )

Types are inferred by
[`infer_duckdb_type()`](https://austin-putz.github.io/tidybreed/reference/infer_duckdb_type.md).
Use R's typed NA values (`NA_integer_`, `NA_real_`, `NA_character_`,
`as.Date(NA)`, etc.) to get the intended DuckDB column type. There is no
typed logical NA in R, so for a `BOOLEAN` column pass a concrete
placeholder (`FALSE` or `TRUE`) instead. Bare `NA` defaults to `VARCHAR`
with a warning.

After creation, add rows via
`DBI::dbAppendTable(pop$db_conn, table_name, my_tibble)` and query via
`get_table(pop, table_name) |> dplyr::collect()`.

## Usage

``` r
define_table(pop, table_name, ..., primary_key = NULL, overwrite = FALSE)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- table_name:

  Character scalar. Name for the new table. Must be a valid SQL
  identifier and must not conflict with a system-managed tidybreed
  table.

- ...:

  Named column definitions. Each name becomes a column; the value's R
  type determines the DuckDB column type. Use typed NAs to declare
  columns without supplying real data (e.g. `run_id = NA_integer_`).

- primary_key:

  Character scalar (optional). Name of one of the `...` columns to
  declare as the `PRIMARY KEY`. When supplied,
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  filtered-row updates will work on this table (DuckDB enforces
  uniqueness). Defaults to `NULL` (no primary key).

- overwrite:

  Logical. If `FALSE` (default), an error is raised when `table_name`
  already exists. If `TRUE`, the existing table is dropped and recreated
  (all data in it will be lost).

## Value

The `tidybreed_pop` (invisibly). Assign the result back.

## See also

[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md),
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- open_pop(pop_name = "demo", db_name = ":memory:") |>
  define_genome(n_loci = 200, n_chr = 2, chr_len_Mb = 100)

# Create a custom table for tracking simulation run metadata
pop <- pop |> define_table(
  "sim_timing",
  run_id       = NA_integer_,
  duration_sec = NA_real_,
  label        = NA_character_,
  primary_key  = "run_id"
)

# Insert rows directly via DBI
DBI::dbAppendTable(
  pop$db_conn, "sim_timing",
  tibble::tibble(run_id = 1L, duration_sec = 42.3, label = "baseline")
)

# Query via the tidy interface
get_table(pop, "sim_timing") |> dplyr::collect()

# Add more columns later
pop <- pop |> get_table("sim_timing") |> mutate_table(notes = NA_character_)

close_pop(pop)
} # }
```

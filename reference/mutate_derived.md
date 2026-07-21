# Compute and write a derived column from a cross-table join

A pipe-friendly action function that:

1.  Collects the current (optionally filtered) primary table

2.  Left-joins an optional secondary table by `join_by`

3.  Applies a user-supplied `compute` function to produce a derived
    vector

4.  Writes the result to one or more destination columns across any
    tables

This reduces ~25-line dplyr join +
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
workflows for computed fields (e.g.
`puberty_date = birth_date + AP_value`) to a single pipe step.

## Usage

``` r
mutate_derived(
  tbl_obj,
  compute,
  join_table = NULL,
  join_by = "id_ind",
  write_to
)
```

## Arguments

- tbl_obj:

  A `tidybreed_table` returned by
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  with an optional
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
  applied.

- compute:

  A function `function(df) -> vector`. `df` is the joined tibble
  (primary table left-joined to `join_table`). The returned vector must
  be a scalar (broadcast to all rows) or have the same length as
  `nrow(df)`. Primary table column names take precedence on name
  collisions; any colliding column from `join_table` is suffixed with
  `".<join_table>"`.

- join_table:

  Character scalar or `NULL`. Name of a secondary table to left-join
  with the primary table. Must exist in `pop$tables`. `NULL` (default)
  skips the join; `compute` receives only the primary table rows.

- join_by:

  Character scalar. Column name to join on; must exist in both tables
  when `join_table` is not `NULL`. Default `"id_ind"`.

- write_to:

  Named character vector. Names are destination table names; values are
  destination column names. Example:
  `c(ind_meta = "puberty_date", ind_phenotype = "pheno_date")`. Each
  destination table must be registered in the package primary-key map
  (`TABLE_PRIMARY_KEYS`). Reserved columns are blocked.

## Value

The parent `tidybreed_pop` object (invisibly).

## Details

### Writing back to the primary table

Uses the table's own primary key for exact-row matching. A filtered
`ind_phenotype` (e.g. `filter(phenotype_name == "AP")`) only updates
those specific rows — other phenotype records for the same individual
are untouched.

### Writing to a secondary table

Requires at most one computed value per `join_by` ID. If the primary
table has multiple rows for the same individual that produce *different*
computed values, `mutate_derived()` errors with a message asking you to
pre-filter. Identical values for the same `join_by` ID are silently
deduplicated (first occurrence wins).

### Column name collisions during the join

If the primary and secondary tables share a column name other than
`join_by`, the secondary table's copy gains a `".<join_table>"` suffix
inside `compute`. Primary table columns are always unambiguous (no
suffix).

### join_table row expansion

If `join_table` itself contains multiple rows per `join_by` value, the
left-join silently expands the primary table rows. Pre-filter
`join_table` externally before passing it if this is a concern.

### New column type inference

When a destination column does not yet exist, its DuckDB type is
inferred from `result_vec` via
[`infer_duckdb_type()`](https://austin-putz.github.io/tidybreed/reference/infer_duckdb_type.md).
R's `Date` arithmetic preserves the `Date` class, so
`birth_date + integer` correctly produces a `DATE` column. Use typed
`NA` vectors (e.g. `NA_real_`) if `result_vec` is all-NA.

## Examples

``` r
if (FALSE) { # \dontrun{
# Age-at-puberty: puberty_date = birth_date (ind_meta) + AP value (days)
pop <- pop |>
  get_table("ind_phenotype") |>
  dplyr::filter(phenotype_name == "AP", pheno_number == 1L) |>
  mutate_derived(
    compute    = \(df) df$birth_date + as.integer(df$pheno_value),
    join_table = "ind_meta",
    join_by    = "id_ind",
    write_to   = c(ind_meta = "puberty_date", ind_phenotype = "pheno_date")
  )

# Compute entirely from the primary table (no join needed)
pop <- pop |>
  get_table("ind_phenotype") |>
  dplyr::filter(phenotype_name == "BW") |>
  mutate_derived(
    compute  = \(df) df$pheno_value * 2.205,   # kg -> lb
    write_to = c(ind_phenotype = "pheno_value_lb")
  )
} # }
```

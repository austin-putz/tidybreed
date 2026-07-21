# Validate a tidybreed_table, col_name, and overwrite state before any RNG.

Returns a named list with all needed state, or NULL when the filter
matched zero rows (callers must propagate the early return).

## Usage

``` r
.resolve_group_target(tbl_obj, col_name, output_type, overwrite)
```

## Arguments

- tbl_obj:

  A `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md).

- col_name:

  Character scalar. Target column name.

- output_type:

  `"INTEGER"` or `"VARCHAR"`.

- overwrite:

  Logical. If `FALSE`, error when filtered rows already have non-NULL
  values in `col_name`.

## Value

A named list with elements `pop`, `conn`, `table_name`, `pk_col`, `pks`
(sorted primary key values of the filtered rows), `n_total`, and
`col_exists`; or `NULL` when the filter matched zero rows (callers must
propagate the early return by returning `invisible(tbl_obj$pop)`).

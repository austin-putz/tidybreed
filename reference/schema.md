# View table-level descriptions for a tidybreed population

Returns a tibble of all user-visible tables with row counts, column
counts, and descriptions from `_schema_meta`. Print the result for an
aligned overview; use `describe_table(pop, "name")` to drill into a
specific table.

## Usage

``` r
schema(pop)
```

## Arguments

- pop:

  A `tidybreed_pop` object, or a `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (its `pop` reference is used).

## Value

A tibble of class `tidybreed_schema` with columns `table_name`,
`n_rows`, `n_cols`, and `description`. Printed via
[`print.tidybreed_schema()`](https://austin-putz.github.io/tidybreed/reference/print.tidybreed_schema.md).

## See also

[`describe_table()`](https://austin-putz.github.io/tidybreed/reference/describe_table.md),
[`define_schema_description()`](https://austin-putz.github.io/tidybreed/reference/define_schema_description.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- open_pop(pop_name = "MySim", db_name = ":memory:") |>
  define_genome(n_loci = 500, n_chr = 5, chr_len_Mb = 100)
schema(pop)
} # }
```

# View column-level descriptions for a tidybreed table

Returns a tibble of all columns in `table_name` with their DuckDB types
and descriptions from `_schema_meta`. Any wide table with `locus_<n>`
columns only lists its non-locus metadata columns to avoid printing
thousands of locus columns.

## Usage

``` r
describe_table(pop, table_name)
```

## Arguments

- pop:

  A `tidybreed_pop` object, or a `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (its `pop` reference is used).

- table_name:

  Character. Name of the table to describe.

## Value

A tibble of class `tidybreed_table_desc` with columns `column_name`,
`column_type`, `description`, and `notes`. Printed via
[`print.tidybreed_table_desc()`](https://austin-putz.github.io/tidybreed/reference/print.tidybreed_table_desc.md).

## See also

[`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md),
[`define_schema_description()`](https://austin-putz.github.io/tidybreed/reference/define_schema_description.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- open_pop(pop_name = "MySim", db_name = ":memory:") |>
  define_genome(n_loci = 500, n_chr = 5, chr_len_Mb = 100)
describe_table(pop, "ind_meta")
describe_table(pop, "genome_effects")

# Also accepts a tidybreed_table, so it chains directly after get_table()
pop |> get_table("ind_meta") |> describe_table("ind_meta")
} # }
```

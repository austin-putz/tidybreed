# Define or update a description for a table or column

Upserts a description into `_schema_meta`. Use this to document custom
columns added via
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
or `...` in
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
etc. Pipe a `tidybreed_table` from
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
as the first argument.

## Usage

``` r
define_schema_description(tbl, column_name = NULL, description, notes = NULL)
```

## Arguments

- tbl:

  A `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md).
  The table name is inferred from `tbl$table_name`.

- column_name:

  Character or `NULL`. When `NULL` (default), the description applies to
  the table itself. When a column name is given, the description applies
  to that specific column.

- description:

  Character. Human-readable description.

- notes:

  Character or `NULL`. Optional supplementary context.

## Value

The `tidybreed_table`, invisibly, enabling back-to-back chained calls
without repeating
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md).
The underlying database is modified in place via the DBI connection — do
not assign the result back to `pop` (use as a side-effect or extract
with `tbl$pop`).
[`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md)
and
[`describe_table()`](https://austin-putz.github.io/tidybreed/reference/describe_table.md)
accept either `tidybreed_pop` or `tidybreed_table`.

## See also

[`describe_table()`](https://austin-putz.github.io/tidybreed/reference/describe_table.md),
[`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Chain multiple column descriptions — no need to repeat get_table()
pop |>
  get_table("ind_meta") |>
  define_schema_description("id_ind",    "Unique individual identifier") |>
  define_schema_description("sex",       "Sex of individual (M or F)")   |>
  define_schema_description("line_name", "Genetic line name")
describe_table(pop, "ind_meta")  # pop still valid; DBI conn is a reference

# Table-level description (column_name = NULL)
pop |>
  get_table("ind_meta") |>
  define_schema_description(description = "Individual metadata table")
} # }
```

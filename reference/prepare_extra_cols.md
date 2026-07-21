# Validate, type-infer, and expand extra user columns for pre-insert attachment

Called by
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md),
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md),
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md),
and
[`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
to process the `...` (extra field) arguments before they are attached to
an insertion data frame.

## Usage

``` r
prepare_extra_cols(extra_cols, n_rows, table_name, conn)
```

## Arguments

- extra_cols:

  Named list of values (captured from `list(...)`).

- n_rows:

  Integer. Number of rows in the insertion data frame.

- table_name:

  Character. Target table name.

- conn:

  DuckDB connection.

## Value

Named list of expanded value vectors, one per field.

## Details

Scalars are broadcast to `n_rows`. Vectors must already have length
`n_rows`. Any column that does not yet exist in the target table is
added via `ALTER TABLE ADD COLUMN` before the caller writes the data
frame.

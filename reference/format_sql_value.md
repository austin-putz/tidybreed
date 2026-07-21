# Format a scalar R value as a SQL literal

Format a scalar R value as a SQL literal

## Usage

``` r
format_sql_value(value, db_type)
```

## Arguments

- value:

  A scalar R value (length-1 vector)

- db_type:

  DuckDB type string from
  [`infer_duckdb_type()`](https://austin-putz.github.io/tidybreed/reference/infer_duckdb_type.md)

## Value

A character string safe to embed in a SQL statement

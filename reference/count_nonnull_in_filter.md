# Count non-NULL values in a filtered set of rows

Used by
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
to report how many existing values will be replaced vs. filled in before
executing an UPDATE.

## Usage

``` r
count_nonnull_in_filter(conn, table_name, field_name, pk_col, filter_ids)
```

## Arguments

- conn:

  DuckDB connection

- table_name:

  Target table name

- field_name:

  Column to inspect

- pk_col:

  Primary key column name

- filter_ids:

  Vector of PK values identifying the affected rows

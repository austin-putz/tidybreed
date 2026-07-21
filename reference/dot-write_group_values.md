# Write PK → value pairs to a table column via a temp-table UPDATE JOIN.

Thin wrapper around
[`mutate_table_vector()`](https://austin-putz.github.io/tidybreed/reference/mutate_table_vector.md)
(defined in mutate_table.R).

## Usage

``` r
.write_group_values(conn, table_name, col_name, pk_col, pks, values)
```

## Arguments

- conn:

  DuckDB connection.

- table_name:

  Target table name.

- col_name:

  Column to write.

- pk_col:

  Primary key column name.

- pks:

  Primary key values, in the same order as `values`.

- values:

  Vector of new values, same length and order as `pks`.

## Value

`NULL`, invisibly.

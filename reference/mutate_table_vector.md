# Execute a vector UPDATE on a table column via a temp-table JOIN

Execute a vector UPDATE on a table column via a temp-table JOIN

## Usage

``` r
mutate_table_vector(conn, table_name, field_name, pks_ordered, values, pk_col)
```

## Arguments

- conn:

  DuckDB connection

- table_name:

  Target table name

- field_name:

  Column to update

- pks_ordered:

  PK values in canonical row order (from `get_pks_in_order()`)

- values:

  Vector of new values, same length and order as `pks_ordered`

- pk_col:

  Primary key column name

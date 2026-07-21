# Execute a scalar UPDATE on a table column

Execute a scalar UPDATE on a table column

## Usage

``` r
mutate_table_scalar(
  conn,
  table_name,
  field_name,
  value,
  db_type,
  filter_ids,
  pk_col
)
```

## Arguments

- conn:

  DuckDB connection

- table_name:

  Target table name

- field_name:

  Column to update

- value:

  Scalar R value

- db_type:

  DuckDB type string

- filter_ids:

  NULL = update all rows; vector of PK values = update only those rows

- pk_col:

  Primary key column name (required when filter_ids is not NULL)

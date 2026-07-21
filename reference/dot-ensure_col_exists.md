# Create a column in a table if it does not already exist.

Create a column in a table if it does not already exist.

## Usage

``` r
.ensure_col_exists(conn, table_name, col_name, output_type)
```

## Arguments

- conn:

  DuckDB connection.

- table_name:

  Target table name.

- col_name:

  Column to create if absent.

- output_type:

  DuckDB type string (e.g. `"INTEGER"`, `"VARCHAR"`) used for the
  `ADD COLUMN` statement when the column does not yet exist.

## Value

`NULL`, invisibly.

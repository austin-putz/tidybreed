# Register table/column descriptions in \_schema_meta (upsert)

Register table/column descriptions in \_schema_meta (upsert)

## Usage

``` r
register_schema_meta(conn, entries)
```

## Arguments

- conn:

  A DBI connection to the DuckDB database.

- entries:

  A data.frame with columns: object_type, table_name, column_name (NA
  for table-level), description, notes (NA ok).

# Per-table on-disk size, in bytes

DuckDB has no direct per-table byte count.
`duckdb_tables().estimated_size` looks right but is estimated
*cardinality*, not bytes. The workable route is `PRAGMA storage_info()`,
counting the distinct blocks a table's segments occupy and multiplying
by the database block size.

## Usage

``` r
.schema_table_bytes(conn, tables)
```

## Arguments

- conn:

  A DBI connection.

- tables:

  Character vector of table names.

## Value

Numeric vector of bytes, aligned to `tables`; `NA` where the pragma
could not be read.

## Details

The caller must `CHECKPOINT` first: before one, recently written data
carries `block_id = -1` and every table reports zero. That is a write,
which is why `sizes` is opt-in on
[`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md)
rather than always on.

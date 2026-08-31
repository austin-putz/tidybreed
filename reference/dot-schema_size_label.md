# Whole-database size label for the `schema()` header

`PRAGMA database_size` returns one row of already-formatted strings, so
no byte arithmetic is needed. Two cases:

## Usage

``` r
.schema_size_label(conn, db_path)
```

## Arguments

- conn:

  A DBI connection.

- db_path:

  Character. `pop$db_path` — the literal `":memory:"` for an in-memory
  population.

## Value

A single string such as `"4.0 MiB on disk"`, or `NA_character_` if the
pragma is unavailable.

## Details

- **File-backed**: report `database_size`, plus `wal_size` when the WAL
  is non-empty. The WAL is only folded into the file by a `CHECKPOINT`,
  so a user looking at the `.duckdb` on disk sees both; reporting only
  `database_size` under-reports.
  [`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md)
  must not checkpoint — that is a write.

- **In-memory** (`db_name = ":memory:"`): `database_size` is `"0 bytes"`
  and `block_size` is `0` because there is no file. Report
  `memory_usage` and label it as memory, not disk.

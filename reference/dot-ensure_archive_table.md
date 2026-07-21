# Ensure an archive table exists and its schema matches the source.

Creates the archive table from the source schema on first call. On
subsequent calls checks that the schemas are identical (strict policy —
errors on any drift, including user-added columns that appear only in
one DB). The `replicate` column added by archive_replicate() is excluded
from the comparison.

## Usage

``` r
.ensure_archive_table(conn, tbl_name, add_replicate)
```

## Value

Character vector of source column names (in schema order), invisibly.

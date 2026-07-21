# Get the next ID for a BIGINT primary-key sequence (overflow-safe).

Like
[`next_int_id()`](https://austin-putz.github.io/tidybreed/reference/next_int_id.md)
but returns a **double** rather than an `integer`, so it is safe for
`BIGINT` PK columns whose values can exceed `.Machine$integer.max`
(`2^31 - 1`). [`as.integer()`](https://rdrr.io/r/base/integer.html)
would return `NA` past that ceiling; a double is exact to `2^53`. Used
for `ind_crossover.id_crossover`, which accumulates far more rows than
any `INTEGER`-keyed table. DuckDB already returns `BIGINT` `MAX()` as an
R double, so no precision is lost on the read.

## Usage

``` r
next_row_id(conn, table, id_col)
```

## Arguments

- conn:

  A DBI connection.

- table:

  Character. Table name.

- id_col:

  Character. BIGINT PK column name.

## Value

Numeric (double) scalar.

## Details

For multi-row inserts, pre-generate a range:
`seq.int(next_row_id(conn, tbl, col), length.out = n)`

# One value of `column` from `table` per id, in `ids` order

Requires exactly one row per id; zero or several rows is an error naming
the table, the column and up to five example ids. The column is returned
in its native type (`NA` for `NULL`).

## Usage

``` r
.read_one_per_id(conn, table, column, ids, what)
```

## Arguments

- what:

  Prefix for error messages, e.g. `"Residual condition lookup"`.

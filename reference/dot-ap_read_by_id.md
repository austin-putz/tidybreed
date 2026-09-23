# Rows of `table` for a set of individuals, by registered-view join

The ids never enter the SQL text. Returns whatever rows match — zero,
one or several per id — so the caller applies its own contract.

## Usage

``` r
.ap_read_by_id(conn, table, ids, columns, where = NULL)
```

## Arguments

- columns:

  Character vector of columns to return besides `id_ind`.

- where:

  Optional extra predicate on `t` (already quoted SQL).

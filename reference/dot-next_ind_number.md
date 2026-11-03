# First individual-id number to use per line

Individual ids are exactly `line_name || '_' || n`, where `n` is a
positive decimal integer with no leading zero and
`1 <= n <= .Machine$integer.max`. For each line in `n_new`, reads the
`ind_meta` rows whose `line_name` column is that line (never an id
prefix, so line `A` does not see line `A_B`), takes everything after the
line's own `line_name_` prefix as `n`, and returns `max(n) + 1` (`1` for
a line with no rows).

## Usage

``` r
.next_ind_number(conn, n_new)
```

## Arguments

- conn:

  A DBI connection.

- n_new:

  Whole numbers named by line: how many new ids each line needs. Names
  must be unique.

## Value

Integer vector named like `n_new`: the first `n` for each line.

## Details

Line names reach DuckDB through a registered view, never SQL text. The
helper only registers a view and reads, so it does not advance R's RNG.

Errors, naming the line, when a stored id of the line is malformed (only
a hand-edited table can hold one) or when the whole allocation
`next .. next + n_new - 1` would leave the supported range. "Current
maximum" is over the rows `ind_meta` holds now: numbers freed by
deletion can be reused, and
[`archive_replicate()`](https://austin-putz.github.io/tidybreed/reference/archive_replicate.md)
restarts numbering.

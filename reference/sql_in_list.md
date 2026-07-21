# Quote and escape a character vector for a SQL `IN (...)` list

Escapes embedded single quotes (quote-doubling, matching
[`format_sql_value()`](https://austin-putz.github.io/tidybreed/reference/format_sql_value.md)'s
VARCHAR branch) and returns the comma-joined, quoted fragment only —
callers still write `paste0("... IN (", x, ")")` around the result.

## Usage

``` r
sql_in_list(values, what = "value")
```

## Arguments

- values:

  Character vector to escape and format for a SQL `IN (...)` list. Must
  be non-empty and free of `NA` (a bare `NA` would otherwise silently
  become the literal string `'NA'` in SQL).

- what:

  Human-readable role label used in error messages (e.g. `"locus name"`,
  `"individual ID"`).

## Value

A single character string: comma-separated, single-quoted, escaped
values.

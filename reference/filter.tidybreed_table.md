# Filter method for tidybreed_table

Stashes filter predicates on the `tidybreed_table` object (for
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
to use in the SQL WHERE clause) AND applies them eagerly to the
underlying lazy dplyr tbl (so
[`collect()`](https://dplyr.tidyverse.org/reference/compute.html)
continues to work correctly for read-only queries).

## Usage

``` r
# S3 method for class 'tidybreed_table'
filter(.data, ..., .preserve = FALSE)
```

## Arguments

- .data:

  A `tidybreed_table` object returned by
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md).

- ...:

  Unquoted predicate expressions.

- .preserve:

  Not used; present for S3 signature compatibility.

## Value

The `tidybreed_table` object with the filter stashed and applied.

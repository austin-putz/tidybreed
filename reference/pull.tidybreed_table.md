# Pull method for tidybreed_table

Pulls a single column as a vector, delegating to the underlying lazy
tbl.

## Usage

``` r
# S3 method for class 'tidybreed_table'
pull(.data, var = -1, name = NULL, ...)
```

## Arguments

- .data:

  A `tidybreed_table` object.

- var:

  Column to pull (tidy-select).

- name:

  Optional column to use as names.

- ...:

  Passed to
  [`dplyr::pull()`](https://dplyr.tidyverse.org/reference/pull.html).

## Value

A vector.

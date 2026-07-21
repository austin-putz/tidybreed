# Collect method for tidybreed_table

Delegates to the underlying lazy dplyr tbl, preserving backward
compatibility for all existing `get_table(...) |> collect()` patterns.

## Usage

``` r
# S3 method for class 'tidybreed_table'
collect(x, ...)
```

## Arguments

- x:

  A `tidybreed_table` object.

- ...:

  Passed to
  [`dplyr::collect()`](https://dplyr.tidyverse.org/reference/compute.html).

## Value

A tibble with the (filtered) table contents.

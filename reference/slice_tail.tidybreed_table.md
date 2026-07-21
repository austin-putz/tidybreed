# slice_tail method for tidybreed_table

Applies
[`dplyr::slice_tail()`](https://dplyr.tidyverse.org/reference/slice.html)
to the underlying lazy tbl and returns a modified `tidybreed_table` so
that further dplyr chains continue to work.

## Usage

``` r
# S3 method for class 'tidybreed_table'
slice_tail(.data, ...)
```

## Arguments

- .data:

  A `tidybreed_table` object.

- ...:

  Passed to
  [`dplyr::slice_tail()`](https://dplyr.tidyverse.org/reference/slice.html)
  (e.g. `n = 10`).

## Value

A `tidybreed_table` with the slice applied.

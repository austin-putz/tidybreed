# slice_min method for tidybreed_table

Applies
[`dplyr::slice_min()`](https://dplyr.tidyverse.org/reference/slice.html)
to the underlying lazy tbl and returns a modified `tidybreed_table` so
that further dplyr chains continue to work.

## Usage

``` r
# S3 method for class 'tidybreed_table'
slice_min(.data, order_by, ...)
```

## Arguments

- .data:

  A `tidybreed_table` object.

- order_by:

  Column or expression to order by (tidy-eval).

- ...:

  Passed to
  [`dplyr::slice_min()`](https://dplyr.tidyverse.org/reference/slice.html).

## Value

A `tidybreed_table` with the slice applied.

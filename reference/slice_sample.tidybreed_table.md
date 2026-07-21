# slice_sample method for tidybreed_table

Applies
[`dplyr::slice_sample()`](https://dplyr.tidyverse.org/reference/slice.html)
to the underlying lazy tbl and returns a modified `tidybreed_table` so
that further dplyr chains continue to work.

## Usage

``` r
# S3 method for class 'tidybreed_table'
slice_sample(.data, ...)
```

## Arguments

- .data:

  A `tidybreed_table` object.

- ...:

  Passed to
  [`dplyr::slice_sample()`](https://dplyr.tidyverse.org/reference/slice.html)
  (e.g. `n = 10`, `prop = 0.1`).

## Value

A `tidybreed_table` with the slice applied.

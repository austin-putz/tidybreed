# Arrange method for tidybreed_table

Applies row ordering to the underlying lazy tbl and returns a modified
`tidybreed_table` so that further dplyr chains continue to work.

## Usage

``` r
# S3 method for class 'tidybreed_table'
arrange(.data, ...)
```

## Arguments

- .data:

  A `tidybreed_table` object.

- ...:

  Ordering expressions.

## Value

A `tidybreed_table` with the ordering applied.

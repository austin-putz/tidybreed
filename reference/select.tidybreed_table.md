# Select method for tidybreed_table

Applies column selection to the underlying lazy tbl and returns a
modified `tidybreed_table` so that further dplyr chains continue to
work.

## Usage

``` r
# S3 method for class 'tidybreed_table'
select(.data, ...)
```

## Arguments

- .data:

  A `tidybreed_table` object.

- ...:

  Column selection (tidy-select).

## Value

A `tidybreed_table` with the selection applied.

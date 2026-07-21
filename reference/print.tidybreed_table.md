# Print method for tidybreed_table

Prints a summary header (rows × fields) plus a data preview of up to `n`
rows and `c` columns. Any active filter is shown below the header.

## Usage

``` r
# S3 method for class 'tidybreed_table'
print(x, n = 10, c = 10, ...)
```

## Arguments

- x:

  A `tidybreed_table` object.

- n:

  Maximum number of rows to preview (default 10).

- c:

  Maximum number of columns to preview (default 10).

- ...:

  Ignored.

## Value

`x` invisibly.

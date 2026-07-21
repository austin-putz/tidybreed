# Close tidybreed population database connection

Close tidybreed population database connection

## Usage

``` r
close_pop(pop, results_dir = NULL)
```

## Arguments

- pop:

  A tidybreed_pop object

- results_dir:

  Optional path to a directory. When provided, the `.duckdb` file is
  moved there after the connection is closed. The directory is created
  if it does not yet exist. Errors if the destination file already
  exists. Ignored silently for in-memory databases.

## Value

`NULL`, invisibly.

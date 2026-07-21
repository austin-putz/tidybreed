# Count method for tidybreed_table

Counts rows (optionally by group), delegating to the underlying lazy
tbl.

## Usage

``` r
# S3 method for class 'tidybreed_table'
count(x, ..., wt = NULL, sort = FALSE, name = "n")
```

## Arguments

- x:

  A `tidybreed_table` object.

- ...:

  Grouping columns (tidy-select).

- wt:

  Optional weighting column.

- sort:

  Logical; sort by descending count?

- name:

  Name for the count column.

## Value

A lazy tibble.

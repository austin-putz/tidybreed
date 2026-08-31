# Abbreviate a row count for display

`ind_haplotype` alone otherwise dictates the width of the `Rows` column
for every table in the printout.

## Usage

``` r
.schema_format_rows(n)
```

## Arguments

- n:

  Integer vector of row counts.

## Value

Character vector: `"?"` for `NA`, `"1.10M"` style above a million,
comma-grouped below it.

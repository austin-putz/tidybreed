# Largest-remainder integer rounding.

Rounds `weights * n` to integers summing exactly to `n`. Tie-breaking is
stable left-to-right (earlier groups win).

## Usage

``` r
.largest_remainder(weights, n)
```

## Arguments

- weights:

  Numeric vector of proportional weights (must sum to 1).

- n:

  Integer. Total count to distribute.

## Value

Integer vector of the same length as `weights`.

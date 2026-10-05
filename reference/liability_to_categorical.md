# Convert a liability vector to ordered integer categories

A record is in category `k + 1` when its liability is strictly **above**
cutpoint `k`; a liability exactly on a cutpoint stays in the lower
category. That is what `prevalence` ("the fraction above the threshold")
means, and it matters for discrete genetic values with no residual.

## Usage

``` r
liability_to_categorical(liability, thresholds)
```

## Arguments

- liability:

  Numeric vector.

- thresholds:

  Numeric vector of cutpoints (ascending).

## Value

Integer vector of category indices (1-based).

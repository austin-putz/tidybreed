# Clip a count vector to `[min_value, max_value]`

Clip a count vector to `[min_value, max_value]`

## Usage

``` r
clip_count(x, min_value = NA_real_, max_value = NA_real_)
```

## Arguments

- x:

  Numeric vector.

- min_value, max_value:

  Numeric scalars or `NA`.

## Value

Integer vector (rounded) clipped to bounds.

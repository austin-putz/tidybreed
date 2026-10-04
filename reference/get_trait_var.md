# Get the variance (diagonal) for one trait from trait_var_comp

Get the variance (diagonal) for one trait from trait_var_comp

## Usage

``` r
get_trait_var(pop, effect_name, trait_name, line_name = NULL)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- effect_name:

  Character.

- trait_name:

  Character.

- line_name:

  `NULL` (population-wide rows) or a line, which falls back to the
  population-wide rows when it has none of its own.

## Value

Numeric scalar, or `NA_real_` if not found.

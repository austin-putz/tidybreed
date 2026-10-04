# Load a full covariance matrix from trait_var_comp

Load a full covariance matrix from trait_var_comp

## Usage

``` r
load_trait_cov(pop, effect_name, trait_names, line_name = NULL)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- effect_name:

  Character.

- trait_names:

  Character vector of trait names.

- line_name:

  `NULL` (population-wide rows) or a line, which falls back to the
  population-wide rows when it has none of its own. Lines are never
  mixed.

## Value

Named numeric matrix, or `NULL` if any entry is missing.

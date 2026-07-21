# Load a full covariance matrix from trait_var_comp

Load a full covariance matrix from trait_var_comp

## Usage

``` r
load_trait_cov(pop, effect_name, trait_names)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- effect_name:

  Character.

- trait_names:

  Character vector of trait names.

## Value

Named numeric matrix, or `NULL` if any entry is missing.

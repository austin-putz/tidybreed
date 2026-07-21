# Load a full covariance matrix from phenotype_var_comp for a named random effect

Load a full covariance matrix from phenotype_var_comp for a named random
effect

## Usage

``` r
load_phenotype_cov(pop, effect_name, phenotype_names)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- effect_name:

  Character. The random effect name (not "residual").

- phenotype_names:

  Character vector of phenotype names.

## Value

Named numeric matrix, or `NULL` if any entry is missing.

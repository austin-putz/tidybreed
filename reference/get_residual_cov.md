# Look up residual covariance information from phenotype_var_comp

Look up residual covariance information from phenotype_var_comp

## Usage

``` r
get_residual_cov(pop, phenotype_names, subset_df = NULL)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- phenotype_names:

  Character vector of phenotype names.

- subset_df:

  Currently unused inside this function (accepted for a future
  per-subset lookup); callers such as
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  separately look up each individual's `condition_column` value and
  match it against `R_by_level`'s names to pick the right (co)variance
  matrix.

## Value

A list:

- `R_unconditional`: numeric matrix (the unconditional residual R) or
  NULL.

- `condition_column`: character scalar or NULL.

- `condition_table`: character scalar or NULL.

- `R_by_level`: named list of matrices keyed by condition_level, or
  NULL.

- `residual_var_unconditional`: named numeric vector of diagonal
  variances.

# Compute the covariate contribution for each individual for a single phenotype

Compute the covariate contribution for each individual for a single
phenotype

## Usage

``` r
compute_covariate_contribution(pop, phenotype_name, subset_df)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- phenotype_name:

  Character. Phenotype name to look up in `phenotype_effects`.

- subset_df:

  Data frame: the per-phenotype subset of `ind_meta` (already
  sex-filtered). Must contain `id_ind` and any `ind_meta` columns
  referenced by effects.

## Value

A list:

- `contribution`: numeric vector of length `nrow(subset_df)`. `NA` for
  individuals excluded by `null_class_action = "skip"`.

- `n_skipped`: integer count of NAs introduced.

# Stage-1 covariate terms of one phenotype: fixed contributions and the random-effect level every record touches

No RNG, no writes. Fixed-class and fixed-covariate effects are evaluated
to a per-record contribution, with `null_class_action = "skip"` marking
the record `NA`. Random effects are **not** drawn here: the function
records which level of the grouping column each record falls in so that
Stage 2 can resolve the draws once the record list is final (see
[add_phenotype_stages](https://austin-putz.github.io/tidybreed/reference/add_phenotype_stages.md)).

## Usage

``` r
.ap_covariate_terms(pop, phenotype_name, subset_df)
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

- `fixed`: numeric vector of length `nrow(subset_df)`; `NA` for
  individuals excluded by `null_class_action = "skip"`.

- `random`: a list, one element per random effect in `effect_name`
  order, each `list(effect_name, distribution, level)` where `level` is
  a character vector of length `nrow(subset_df)` (`NA` = no level).

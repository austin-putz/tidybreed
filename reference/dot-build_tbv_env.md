# Pre-fetch all TBV vectors needed by a formula_tbv expression.

Calls .fetch_contributor_tbvs() or .fetch_group_tbvs() for each
trait_ref, returning a named list ready to be used as an eval()
environment.

## Usage

``` r
.build_tbv_env(pop, trait_refs, subset_df, phenotype_name)
```

## Arguments

- pop:

  A tidybreed_pop object.

- trait_refs:

  List from .walk_formula_tbv_ast()\$trait_refs.

- subset_df:

  Data frame: sex-filtered ind_meta rows for focal individuals.

- phenotype_name:

  Character. Used in error messages.

## Value

Named list: placeholder → numeric vector (length = nrow(subset_df)).

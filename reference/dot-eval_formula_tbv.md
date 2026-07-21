# Evaluate a formula_tbv string for a set of individuals.

Orchestrates: parse → AST walk → TBV pre-fetch → AST substitution →
eval().

## Usage

``` r
.eval_formula_tbv(pop, formula_tbv, subset_df, phenotype_name)
```

## Arguments

- pop:

  A tidybreed_pop object.

- formula_tbv:

  Character. DSL formula string from phenotype_meta.

- subset_df:

  Data frame: sex-filtered ind_meta rows.

- phenotype_name:

  Character. Used in error messages.

## Value

Named numeric vector (names = id_ind). NA marks excluded individuals
(missing dam/sire TBV, or NA group membership).

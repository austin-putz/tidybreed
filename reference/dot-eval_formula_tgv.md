# Evaluate a formula_tgv string for a set of individuals.

Orchestrates: parse → AST walk (references replaced by placeholders) →
genetic-value pre-fetch → eval().

## Usage

``` r
.eval_formula_tgv(pop, formula_tgv, subset_df, phenotype_name)
```

## Arguments

- pop:

  A tidybreed_pop object.

- formula_tgv:

  Character. DSL formula string from phenotype_meta.

- subset_df:

  Data frame: sex-filtered ind_meta rows.

- phenotype_name:

  Character. Used in error messages.

## Value

Named numeric vector (names = id_ind). NA marks excluded individuals (a
missing dam/sire genetic value, or NA group membership). A constant
expression is broadcast to every individual. An `Inf`, `-Inf` or `NaN`
from the arithmetic itself (division by zero, overflow, a domain error)
is an error: it is a fault of the model, not a missing component.

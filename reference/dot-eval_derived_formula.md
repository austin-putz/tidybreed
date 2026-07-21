# Evaluate a derived formula over existing ind_phenotype records.

Evaluate a derived formula over existing ind_phenotype records.

## Usage

``` r
.eval_derived_formula(pop, formula, ids, phenotype_name)
```

## Arguments

- pop:

  A tidybreed_pop object.

- formula:

  Character. Formula string from phenotype_meta.

- ids:

  Character vector. Individual IDs to compute for.

- phenotype_name:

  Character. Name of the derived phenotype (for messages).

## Value

Numeric vector, same length as ids. NA propagates naturally; Inf/NaN
converted to NA with a warning.

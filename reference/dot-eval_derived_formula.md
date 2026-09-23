# Evaluate a derived formula over ind_phenotype records.

Evaluate a derived formula over ind_phenotype records.

## Usage

``` r
.eval_derived_formula(pop, formula, ids, phenotype_name, pending = NULL)
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

- pending:

  Optional data frame of records planned earlier in the same
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  call (`id_ind`, `phenotype_name`, `pheno_value`, `pheno_number`) that
  are not on disk yet. They are treated exactly as if they had been
  written, so a derived phenotype can consume a feeder phenotype from
  the same call.

## Value

Numeric vector, same length as ids. NA propagates naturally; Inf/NaN
converted to NA with a warning.

# Plan the records of one phenotype

Applies the covariate skip and the genetic-value exclusions, reads the
genetic value, and assigns `pheno_number`. Emits the same warnings and
messages the exclusions always have. No RNG, no writes.

## Usage

``` r
.ap_plan_phenotype(
  pop,
  t,
  m,
  subset_df,
  path,
  tgv_kind,
  formula_tgv,
  formula,
  comp_rows,
  user_values,
  n_phenos
)
```

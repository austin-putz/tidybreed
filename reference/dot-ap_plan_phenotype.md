# Plan the records of one phenotype

Applies the covariate skip and the TBV exclusions, reads the TBV, and
assigns `pheno_number`. Emits the same warnings and messages the
exclusions always have. No RNG, no writes.

## Usage

``` r
.ap_plan_phenotype(
  pop,
  t,
  m,
  subset_df,
  path,
  tbv_kind,
  formula_tbv,
  formula,
  comp_rows,
  user_values,
  n_phenos
)
```

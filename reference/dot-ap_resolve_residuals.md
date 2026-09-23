# The residual adapter: every model-path residual of the call

Implements `plans/sample_correlated_effects.md` §5.3–§5.5 for
`effect_name = 'residual'`. Every model-path phenotype with at least one
planned record must belong to a residual covariance block (else "No
residual variance found") unless its residuals are all fixed by
`user_residual`, in which case nothing is drawn or conditioned for it.
Blocks are the plan's, processed in
[`find_covariance_blocks()`](https://austin-putz.github.io/tidybreed/reference/find_covariance_blocks.md)
order; within a block the entity is `(id_ind, pheno_number)`, the
coordinates are the block's phenotypes, and each entity's residual is
drawn from the stratum its condition value selects, conditional on what
it has already realized — stored residuals of *any* block member at the
same `pheno_number` (subject to D2) plus the `user_residual` values
fixed in this call.

## Usage

``` r
.ap_resolve_residuals(plan, fixed = list())
```

## Arguments

- plan:

  The Stage-1 plan.

- fixed:

  The validated `user_residual` list from
  [`.ap_fixed_residuals()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_fixed_residuals.md).

## Value

A list named by model-path phenotype (those with planned records). Each
element has `value` (numeric per planned record), `level` (the
`residual_condition_level` per record: the selected stratum, `NA` for
the unconditional `R`) and `var_unconditional` (the phenotype's
unconditional residual variance, `NA` if no unconditional stratum is
stored).

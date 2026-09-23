# Resolve the residuals of one covariance block

Resolve the residuals of one covariance block

## Usage

``` r
.ap_residual_block(plan, b, targets, fixed)
```

## Arguments

- plan:

  The Stage-1 plan.

- b:

  One block from
  [`find_covariance_blocks()`](https://austin-putz.github.io/tidybreed/reference/find_covariance_blocks.md).

- targets:

  The model-path phenotypes of the call with planned records.

- fixed:

  The validated `user_residual` list.

## Value

A list named by the block's in-call phenotypes, each
`list(value, level)` per planned record.

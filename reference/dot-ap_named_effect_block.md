# Resolve the draws of one named-effect covariance block

Resolve the draws of one named-effect covariance block

## Usage

``` r
.ap_named_effect_block(plan, b, targets)
```

## Arguments

- plan:

  The Stage-1 plan.

- b:

  One block from
  [`find_covariance_blocks()`](https://austin-putz.github.io/tidybreed/reference/find_covariance_blocks.md)
  for a named effect.

- targets:

  The model-path phenotypes of the call with planned records and a
  random term for `b$effect_name`.

## Value

A list with `contribution` (named by the block's in-call phenotypes, one
value per planned record) and `pending` (new `phenotype_random_effects`
rows: every level drawn here, by coordinate then level).

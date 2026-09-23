# Stage 2: every random draw and the in-memory record assembly

Order of RNG consumption, fixed regardless of database row order: the
named-effect adapter first
([`.ap_resolve_named_effects()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_resolve_named_effects.md):
effects in byte-sorted `effect_name` order, blocks in
[`find_covariance_blocks()`](https://austin-putz.github.io/tidybreed/reference/find_covariance_blocks.md)
order, one resolver call per sample-set group in sorted group order),
then the residual adapter
([`.ap_resolve_residuals()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_resolve_residuals.md):
one block at a time, in loader order; within a block one resolver call
per `(stratum, sample set)` group in sorted group order). Record
assembly (derived formulas, liability, type conversion) follows in plan
order and draws nothing. Nothing is written.

## Usage

``` r
.ap_resolve(plan, user_residual = NULL)
```

## Arguments

- plan:

  The Stage-1 plan.

- user_residual:

  The `user_residual` argument, or `NULL`.

## Value

A list with `records` (one tibble per phenotype, in plan order, possibly
empty) and `random_effects` (a data frame of new
`phenotype_random_effects` rows, possibly empty).

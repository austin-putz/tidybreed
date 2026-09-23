# The named-effect adapter: every random-effect draw of the call

Implements `plans/sample_correlated_effects.md` §5.6 and §5.8 for every
`effect_name` other than `'residual'`. The entity is the *level* — a
pen, a herd, an `id_ind` for a permanent-environment effect — and a
level's draw is realized once and reused by every record that ever
touches it, in this call or any later one. Effects are processed in
byte-sorted order, each through its blocks (from the plan, in
[`find_covariance_blocks()`](https://austin-putz.github.io/tidybreed/reference/find_covariance_blocks.md)
order): a level draws its planned coordinates conditional on the
coordinates already stored in `phenotype_random_effects` for the block's
other phenotypes. Every model-path phenotype with a random term for the
effect must be in a block (else "No variance stored"); the §5.6 checks
are re-run here as the backstop; a 1 x 1 block whose effect is `gamma`
or `uniform` keeps its marginal sampler.

## Usage

``` r
.ap_resolve_named_effects(plan)
```

## Arguments

- plan:

  The Stage-1 plan.

## Value

A list with `contribution` (named by model-path phenotype with planned
records: the summed random-effect value per record, `0` for a record
whose level is `NULL`) and `pending` (the new `phenotype_random_effects`
rows for Stage 3).

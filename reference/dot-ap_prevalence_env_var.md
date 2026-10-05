# The named random-effect variance a prevalence threshold uses

The liability carries every named random effect of the phenotype, so the
threshold's variance includes their stored variances (each effect's
diagonal in `phenotype_var_comp`). `normal` and `uniform` effects are
centred with that variance and enter it; a `gamma` effect (shape 1) has
mean `sqrt(variance)`, which no threshold from `mean` accounts for, and
is refused.

## Usage

``` r
.ap_prevalence_env_var(pop, t)
```

## Value

The summed random-effect variance (a number, `0` with none).

## Details

Also refuses a residual block with conditional strata: the unconditional
stratum is then only the fallback for records whose level has no stratum
of its own, not the population's marginal residual variance, which would
need the levels' frequencies.

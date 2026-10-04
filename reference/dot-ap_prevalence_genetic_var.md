# The genetic variance a prevalence threshold uses (the active-block rule)

The liability carries the trait's total genetic value, so the threshold
uses the variance of the **active** model
(plans/import_qtl_effect_methods.md §6A): the population-wide
(`line_name IS NULL`) stored diagonal of `additive`, `dominance` and
`additive_by_additive`, each counted only if the trait's model has terms
of that kind (any scope, any owner). A stored target for a kind the
model lacks — one a generator was told to leave out — does not enter.
Errors, naming `define_phenotype(thresholds = )`, when the model has
terms outside the three blocks (an `indicator` surface, any interaction
other than additive-by-additive) or terms of a kind with no stored
target.

## Usage

``` r
.ap_prevalence_genetic_var(pop, t)
```

## Value

The summed variance (a number).

## Details

The result is the *target* at the reference population, an approximation
for a selected or line-scoped population, as the roxygen of
[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
says.

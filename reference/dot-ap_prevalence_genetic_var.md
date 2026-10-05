# The genetic variance a prevalence threshold uses (the active-block rule)

The liability carries the trait's total genetic value, so the threshold
uses the variance of the **active** model
(plans/import_qtl_effect_methods.md §6A): the population-wide
(`line_name IS NULL`) stored diagonal of `additive`, `dominance` and
`additive_by_additive`, each counted only if the trait's model has terms
of that kind (any scope). A stored target for a kind the model lacks —
one a generator was told to leave out — does not enter.

## Usage

``` r
.ap_prevalence_genetic_var(pop, t)
```

## Value

The summed variance (a number).

## Details

**The owner rule (Q21).** Every term must be owned by `"generated"`:
only a generator writes that owner, and a generator always calibrates
its terms to the stored target, so the target is known to describe them.
A term from
[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md)
(any user owner) carries values nobody checked against the target, and
the threshold would silently miss the requested prevalence.

Errors, naming `define_phenotype(thresholds = )`, when any term is not
`"generated"`, when the model has terms outside the three blocks (an
`indicator` surface, any interaction other than additive-by-additive),
or when it has terms of a kind with no stored target. Also errors when
one kind at one line scope has generated variants for more than one
parent-of-origin scope: each variant is calibrated to the target alone,
so "generated" proves a variant's calibration, not the whole trait's
variance.

The sum of diagonals assumes orthogonal components (statistical coding
at one HWE/LE base). The result is the *target* at the reference
population, an approximation for a selected or line-scoped population,
as the roxygen of
[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
says.

# Error if any of the given loci sit on a non-"full" chromosome

[`rescale_effects_to_target()`](https://austin-putz.github.io/tidybreed/reference/rescale_effects_to_target.md)'s
Falconer variance formula (`2 * p * (1-p) * a^2`) assumes every QTL is
diploid/autosomal — it would silently overstate the variance
contribution of any QTL on a hemizygous (`copy_mode != "full"`) locus,
since a single-copy Bernoulli(p) draw has genic variance `p*(1-p)*a^2`,
not `2*p*(1-p)*a^2`. Rather than generalizing the formula (an
approximation that would need to average over the sexes actually present
in the population),
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
fails loudly instead when `scale_to_target = TRUE` and any selected
locus is sex-linked/organelle.

## Usage

``` r
assert_qtl_autosomal(conn, locus_names)
```

## Arguments

- conn:

  A DBI connection.

- locus_names:

  Character vector of locus names to check.

## Value

Invisible `NULL` on success; errors otherwise.

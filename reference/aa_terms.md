# Build `terms` rows for additive-by-additive pairs

Expands one coefficient per locus pair into the two-member
`additive x additive` term
[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md)
takes: two rows per pair, both `contrast_name = "additive"`, carrying
the same coefficient.

- `coding = "functional"` — the pair contributes
  `e * (g_1 - 1) * (g_2 - 1)` (centres 0.5), where `g` is the dosage of
  allele 1. Values are absolute, so the term has a non-zero mean.

- `coding = "cockerham"` — the pair contributes
  `e * (g_1 - 2 p_1) * (g_2 - 2 p_2)` (centres `p_1`, `p_2`), which has
  mean 0 under HWE and linkage equilibrium.

Within a pair the two loci are put in C-locale order (each `p` moves
with its locus), so `(L10, L44)` and `(L44, L10)` give the same rows.
`e` is never rescaled. A pair with `e = 0` is dropped, and every `e`
being zero is an error, as
[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md)
does for `a` and `d`. Every input is validated before zero pairs are
dropped.

## Usage

``` r
aa_terms(
  locus_name_1,
  locus_name_2,
  e,
  p_1,
  p_2,
  coding = c("functional", "cockerham"),
  effect_name = NULL,
  report = TRUE
)
```

## Arguments

- locus_name_1, locus_name_2:

  Character vectors of the same length, one pair per position. A pair's
  two loci must differ.

- e:

  Numeric coefficient per pair, recycled.

- p_1, p_2:

  Allele-1 frequency of each pair's first and second locus, recycled.
  Required under both codings: the centres under Cockerham, the basis of
  the reported mean under functional.

- coding:

  `"functional"` (default) or `"cockerham"`.

- effect_name:

  Optional label carried onto every term.

- report:

  Logical; print the implied genetic mean (functional coding). Default
  `TRUE`.

## Value

A `terms` data frame with the fixed builder column set: two rows per
pair with a non-zero `e`.

## Every builder returns the same columns

[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md),
`aa_terms()` and
[`genotype_terms()`](https://austin-putz.github.io/tidybreed/reference/genotype_terms.md)
return exactly `term_id`, `locus_name`, `contrast_name`, `center_value`,
`copy_count_value`, `dosage_value`, `genome_value`, `effect_name`, in
that order and with those types (character, character, character,
double, integer, integer, double, character), with a typed `NA` where a
column does not apply. So their outputs
[`rbind()`](https://rdrr.io/r/base/cbind.html) in any combination. A
`term_id` encodes the builder and the term's loci unambiguously,
whatever characters the locus names contain, so bound outputs never
share a `term_id` unless they describe the same term on the same loci.
(The writer may still refuse overlapping definitions on the same loci;
that is its family rule, not an id collision.)

## Functional `a` is not Cockerham `alpha` once pairs exist

`ad_terms(coding = "cockerham")` takes the average effect `alpha`, not
the functional `a`. Already without pairs `alpha = a + (q - p) d`; with
pairs it also needs `sum_l e_jl (2 p_l - 1)`, which
[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md)
never sees. Feeding functional `a` to the Cockerham branch therefore
changes the **genotypic values**, not just how the variance is split
between components. For a hand-written model with pairs, use functional
coding for both builders.

## The reported mean is not written anywhere

Under functional coding a
[`message()`](https://rdrr.io/r/base/message.html) reports each pair's
share of the implied genetic mean, `e (2 p_1 - 1)(2 p_2 - 1)` (its
expectation under HWE and linkage equilibrium), and their running total.
It is written to no table, for the reason given in
[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md).

## See also

[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md),
[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md),
[`genotype_terms()`](https://austin-putz.github.io/tidybreed/reference/genotype_terms.md).

## Examples

``` r
if (FALSE) { # \dontrun{
pop |>
  define_genome_effect_terms(
    trait_name = "ADG",
    terms = rbind(
      ad_terms(locus_name = c("Locus_10", "Locus_44"),
               a = c(0.30, -0.12), d = c(0.10, 0.05), p = c(0.35, 0.60)),
      aa_terms(locus_name_1 = "Locus_10", locus_name_2 = "Locus_44",
               e = 0.08, p_1 = 0.35, p_2 = 0.60)))
} # }
```

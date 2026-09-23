# Build `terms` rows for a functional or Cockerham (a, d) pair

Expands one locus's additive and dominance coefficients into the two
member rows
[`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md)
takes. There is no third table and no separate storage mode — functional
coding is `additive` with `center_value = 0.5` plus an `indicator` on
the heterozygous state, Cockerham coding is `additive` plus `dominance`,
both centred at `p`.

- `coding = "functional"` — the genotypic value is
  `a * (g - 1) + d * 1[g = 1]`, where `g` is the dosage of allele 1.
  Values are **absolute**, so the genetic component has a non-zero mean
  by construction.

- `coding = "cockerham"` — `a` is the average effect `alpha` and `d` the
  dominance deviation `delta` on the HWE-orthogonal contrast (`-2p^2`,
  `2pq`, `-2q^2`). Values are deviations, with mean 0.

## Usage

``` r
ad_terms(
  locus_name,
  a,
  d = 0,
  p,
  coding = c("functional", "cockerham"),
  effect_name = NULL,
  report = TRUE
)
```

## Arguments

- locus_name:

  Character vector of locus names.

- a:

  Numeric additive coefficient(s), recycled to `locus_name`.

- d:

  Numeric dominance coefficient(s), recycled. Default `0`.

- p:

  Numeric allele-1 frequency per locus, recycled. Required for
  `"cockerham"`; used for the reported mean under `"functional"`.

- coding:

  `"functional"` (default) or `"cockerham"`.

- effect_name:

  Optional label carried onto both members of each locus.

- report:

  Logical; print the implied genetic mean. Default `TRUE`.

## Value

A `terms` data frame: two rows per locus with a non-zero coefficient,
`term_id` `"<locus>_a"` / `"<locus>_d"`.

## The reported mean is not written anywhere

For functional coding the implied mean of the **genetic component** is
`mu = a(p - q) + 2pq*d`. It is reported and stored nowhere — in
particular not in `phenotype_meta.mean`, which would double-count it
once non-additive genetic values reach the phenotype layer, because the
raw genetic values already have expectation `mu`. The running total
across loci is reported only because every term here is a single-locus
main effect; no such total exists for an epistatic term, whose
expectation depends on joint genotype frequencies and LD.

## See also

[`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md),
[`genotype_terms()`](https://austin-putz.github.io/tidybreed/reference/genotype_terms.md).

## Examples

``` r
if (FALSE) { # \dontrun{
tt <- ad_terms("Locus_10", a = 0.4, d = 0.2, p = 0.3)
pop <- pop |> define_genome_effects("ADG", tt, effect_owner = "functional")
} # }
```

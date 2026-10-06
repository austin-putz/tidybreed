# Build `terms` rows from a table of genotype values

Turns a genotype-by-value table into `indicator` terms: one term per row
of `genotypes`, one member per locus column. This is indicator
completeness made concrete — an arbitrary function of a finite set of
genotype states *is* a sum of indicator terms, so a hand-entered surface
needs no second representation, only rows.

A cell you leave out is a term you did not write, contributing zero. Ask
for the opposite with
`define_genome_effect_terms(require_complete = TRUE)`, which then
demands every reachable `(copy_count, dosage)` state.

## Usage

``` r
genotype_terms(
  genotypes,
  value,
  copy_count = NULL,
  drop_zero = TRUE,
  effect_name = NULL
)
```

## Arguments

- genotypes:

  A data frame whose columns are named by `locus_name` and hold
  **dosages** of allele 1, one row per genotype combination. A
  `copy_count` list may accompany it for variable-copy loci.

- value:

  Numeric vector of genotypic values, one per row of `genotypes`.

- copy_count:

  Optional named list or vector giving `copy_count_value` per locus
  column. Omitted, it is inferred by the writer at diploid-autosomal
  loci; a variable-copy locus needs it, or a `copy_count` column-shaped
  data frame matching `genotypes`.

- drop_zero:

  Logical; drop rows whose `value` is exactly 0, since they contribute
  nothing. Default `TRUE`.

- effect_name:

  Optional label carried onto every term.

## Value

A `terms` data frame with the fixed builder column set (see
[`aa_terms()`](https://austin-putz.github.io/tidybreed/reference/aa_terms.md))
and `nrow(genotypes) * ncol(genotypes)` rows (before `drop_zero`). Each
row of `genotypes` is one term; its `term_id` encodes the surface's loci
and the row number, so surfaces over different loci never share a
`term_id` when bound together.

## See also

[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md),
[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md),
[`aa_terms()`](https://austin-putz.github.io/tidybreed/reference/aa_terms.md).

## Examples

``` r
if (FALSE) { # \dontrun{
cells <- expand.grid(Locus_10 = 0:2, Locus_44 = 0:2)
vals  <- c(0, 0, 0,  0, 1.4, 2.1,  0, 2.1, 3.6)
pop <- pop |> define_genome_effect_terms(
  "ADG", genotype_terms(cells, vals), effect_owner = "epistasis_AxA")
} # }
```

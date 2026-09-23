# Extract per-locus allele frequencies from a filtered table

Computes the frequency of allele 1 at every locus in `genome_meta` from
the allele copies a filtered `tidybreed_table` selects. The table's
identity says *what kind of thing* is being selected; the user's
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
says *which ones*:

|  |  |  |
|----|----|----|
| Piped table | Meaning | Copies counted |
| `founder_haplotypes` | the founder pool | the filtered pool rows |
| `ind_haplotype` | these allele copies | the filtered haplotype rows |
| any other table with `id_ind` | these individuals | every `ind_haplotype` row of the **distinct** selected `id_ind` |

The frequency is computed in one SQL statement with the filter rendered
as a subquery; nothing but the per-locus result is collected into R. For
an `id_ind` table the frequency depends only on which individuals are
selected, never on how many rows each has — `ind_phenotype` with five
records per animal gives the same answer as `ind_meta` for the same
animals.

This is the single place tidybreed turns a population selection into
`p`.
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
and
[`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md)
both call it for their `base_tbl`, so a base selection means the same
population in either. Users call it to obtain `p` for
[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md).

## Usage

``` r
extract_allele_freq(tbl)
```

## Arguments

- tbl:

  A `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md),
  optionally filtered. One of the three shapes above. Columns removed
  with
  [`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
  are checked before any SQL runs.

## Value

A tibble with exactly one row per `genome_meta` locus in ascending
`locus_id` order: `locus_id` (integer), `locus_name` (character),
`allele_freq` (double). `allele_freq` is `NA` at a locus the selection
has no copies for; it is never `0` in that case, and an observed
frequency of exactly `0` or `1` is a real value and is kept. Errors if
no locus has a copy at all. Never warns and never writes.

## Two notions of line

`ind_meta.line_name` is a pedigree label — an F1 is whatever it was
called; `ind_haplotype.line_origin` is the founding line each allele
copy traces to, carried through every
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
call. They coincide for purebreds and diverge for crosses.
`get_table(pop, "ind_meta") |> filter(line_name == "F1")` counts every
copy those animals carry, whichever line it came from;
`get_table(pop, "ind_haplotype") |> filter(line_origin == "Duroc")`
counts Duroc copies wherever they sit, at any cross depth.

## See also

[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md),
[`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md),
[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Founder pool of one line
p <- pop |> get_table("founder_haplotypes") |>
  dplyr::filter(line_name == "Duroc") |> extract_allele_freq()

# Generation-0 animals
p <- pop |> get_table("ind_meta") |> dplyr::filter(gen == 0L) |>
  extract_allele_freq()

# Duroc copies wherever they sit, including inside crossbreds
p <- pop |> get_table("ind_haplotype") |>
  dplyr::filter(line_origin == "Duroc") |> extract_allele_freq()

# Feed p to ad_terms()
loci <- c("Locus_10", "Locus_44")
tt <- ad_terms(loci, a = c(0.4, 0.1), d = c(0.2, 0.05),
               p = p$allele_freq[match(loci, p$locus_name)],
               coding = "cockerham")
} # }
```

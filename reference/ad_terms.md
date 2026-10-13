# Build `terms` rows for a functional or Cockerham (a, d) pair

Turns each locus's additive and dominance coefficients into the `terms`
rows
[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md)
writes: one row per non-zero coefficient. Functional coding gives an
`additive` row centred at `0.5` and an `indicator` row on the
heterozygote; Cockerham coding gives `additive` and `dominance` rows
centred at `p`.

- `coding = "functional"` (default): the genotypic value is
  `a * (g - 1) + d * 1[g = 1]`, where `g` is the dosage of allele 1 (0,
  1 or 2). So `a` is half the difference between the two homozygotes and
  `d` the heterozygote's deviation from their midpoint. Values are
  **absolute**, so the genetic values have a non-zero mean.

- `coding = "cockerham"`: `a` is the average effect `alpha` and `d` the
  dominance effect on the HWE dominance contrast, which is `-2p^2`,
  `2pq`, `-2q^2` at dosages 0, 1, 2. Values are deviations with mean 0
  under HWE at `p`.

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

  Character vector of distinct locus names (from `genome_meta`).

- a:

  Numeric additive coefficient(s), one per locus or one recycled.

- d:

  Numeric dominance coefficient(s), one per locus or one recycled.
  Default `0` (no dominance row).

- p:

  Allele-1 frequency per locus, in `[0, 1]`, one per locus or one
  recycled. Always required: the centre under Cockerham coding, and the
  basis of the reported mean under functional coding (where it does not
  change the stored rows). Use the frequency of the population the model
  describes, e.g. from
  [`extract_allele_freq()`](https://austin-putz.github.io/tidybreed/reference/extract_allele_freq.md).

- coding:

  `"functional"` (default) or `"cockerham"`; see *Description*.

- effect_name:

  Optional label carried onto every row.

- report:

  Logical; print the implied genetic mean as a message. Default `TRUE`.

## Value

A `terms` data frame with the fixed builder column set (see
[`aa_terms()`](https://austin-putz.github.io/tidybreed/reference/aa_terms.md)):
one row per non-zero coefficient. `term_id` encodes the builder, the
locus and `a` / `d`, so it never collides with another builder's ids
when outputs are [`rbind()`](https://rdrr.io/r/base/cbind.html)-ed.

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

[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md),
[`aa_terms()`](https://austin-putz.github.io/tidybreed/reference/aa_terms.md),
[`genotype_terms()`](https://austin-putz.github.io/tidybreed/reference/genotype_terms.md).

## Examples

``` r
# The rows only; nothing is written
ad_terms(c("Locus_10", "Locus_20"), a = c(0.4, -0.2), d = c(0.1, 0),
         p = c(0.3, 0.5))
#> Implied genetic mean (functional coding), reported only — written to no table: Locus_10 mu = -0.118; Locus_20 mu =  0.000. Running total over these single-locus main effects: -0.118.
#>           term_id locus_name contrast_name center_value copy_count_value
#> 1 ad:8:Locus_10#a   Locus_10      additive          0.5               NA
#> 2 ad:8:Locus_20#a   Locus_20      additive          0.5               NA
#> 3 ad:8:Locus_10#d   Locus_10     indicator           NA                2
#>   dosage_value genome_value effect_name
#> 1           NA          0.4        <NA>
#> 2           NA         -0.2        <NA>
#> 3            1          0.1        <NA>

# The same model's statistical (Cockerham) coefficients at p = 0.3
# (Locus_10 only): alpha = a + (q - p) d
ad_terms("Locus_10", a = 0.4 + (0.7 - 0.3) * 0.1, d = 0.1, p = 0.3,
         coding = "cockerham")
#> Implied genetic mean (cockerham coding), reported only — written to no table: Locus_10 mu = 0. Running total over these single-locus main effects: 0.
#>           term_id locus_name contrast_name center_value copy_count_value
#> 1 ad:8:Locus_10#a   Locus_10      additive          0.3               NA
#> 2 ad:8:Locus_10#d   Locus_10     dominance          0.3               NA
#>   dosage_value genome_value effect_name
#> 1           NA         0.44        <NA>
#> 2           NA         0.10        <NA>

# Write them for a trait
set.seed(1)
pop <- open_pop(pop_name = "ad_example", db_name = ":memory:") |>
  define_genome(n_loci = 60, n_chr = 2, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 100)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpYogPXV/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> Opened in-memory population 'ad_example'
#> Defined genome: 60 loci across 2 chromosomes | chr lengths (Mb): all equal to 100 Mb
#>   Tables written: genome_meta, genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_inheritance, chr_recombination, genome_effects, genome_effect_members, genome_effect_member_origins
#>   Views created: genome_effect_terms, genome_effect_loci
#> Generating 100 founder haplotypes (method = "uniform")...
#>   Created founder_haplotypes (100 haplotypes x 60 loci)
pop <- pop |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 10, n_females = 10, line_name = "A")
#> Added 20 founders (10 males, 10 females) to line 'A'
pop <- define_trait(pop, "ADG")
#> Added trait 'ADG'.
pop <- pop |>
  define_genome_effect_terms("ADG",
    ad_terms("Locus_10", a = 0.4, d = 0.1, p = 0.3),
    effect_owner = "qtl_map")
#> Implied genetic mean (functional coding), reported only — written to no table: Locus_10 mu = -0.118. Running total over these single-locus main effects: -0.118.
#> Wrote 2 genome-effect terms (2 member rows) for trait 'ADG' under owner 'qtl_map' [append].
close_pop(pop)
```

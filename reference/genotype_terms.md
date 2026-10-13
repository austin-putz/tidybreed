# Build `terms` rows from a table of genotype values

Turns a table of genotypic values into `indicator` terms: one term per
row of `genotypes`, adding `value` to every individual with exactly that
genotype at those loci. Any function of the genotypes at a few loci can
be written this way, e.g. a two-locus epistasis table read from a paper.

A genotype you leave out contributes zero. To be told when one is
missing, write with
`define_genome_effect_terms(require_complete = TRUE)`, which then
requires every reachable genotype state.

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

  A data frame with one column per locus, named by `locus_name`, holding
  the **dosage** of allele 1 (0, 1 or 2 at a diploid locus); one row per
  genotype combination.
  [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html) builds a
  full table.

- value:

  Numeric vector of genotypic values, one per row of `genotypes`.

- copy_count:

  Optional `copy_count_value` (copies present) per locus column: a named
  list or vector, or a data frame shaped like `genotypes`. Leave it
  `NULL` at ordinary diploid autosomal loci, where it is 2; a locus
  whose copy number varies (e.g. sex chromosomes) needs it.

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
# A two-locus table: value only when both loci carry allele 1
cells <- expand.grid(Locus_10 = 0:2, Locus_44 = 0:2)
vals  <- c(0, 0, 0,  0, 1.4, 2.1,  0, 2.1, 3.6)
genotype_terms(cells, vals)
#>                        term_id locus_name contrast_name center_value
#> 1 geno:8:Locus_10|8:Locus_44#5   Locus_10     indicator           NA
#> 2 geno:8:Locus_10|8:Locus_44#6   Locus_10     indicator           NA
#> 3 geno:8:Locus_10|8:Locus_44#8   Locus_10     indicator           NA
#> 4 geno:8:Locus_10|8:Locus_44#9   Locus_10     indicator           NA
#> 5 geno:8:Locus_10|8:Locus_44#5   Locus_44     indicator           NA
#> 6 geno:8:Locus_10|8:Locus_44#6   Locus_44     indicator           NA
#> 7 geno:8:Locus_10|8:Locus_44#8   Locus_44     indicator           NA
#> 8 geno:8:Locus_10|8:Locus_44#9   Locus_44     indicator           NA
#>   copy_count_value dosage_value genome_value effect_name
#> 1               NA            1          1.4        <NA>
#> 2               NA            2          2.1        <NA>
#> 3               NA            1          2.1        <NA>
#> 4               NA            2          3.6        <NA>
#> 5               NA            1          1.4        <NA>
#> 6               NA            1          2.1        <NA>
#> 7               NA            2          2.1        <NA>
#> 8               NA            2          3.6        <NA>

set.seed(1)
pop <- open_pop(pop_name = "genotype_example", db_name = ":memory:") |>
  define_genome(n_loci = 60, n_chr = 2, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 100)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpOJ05Ue/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> Opened in-memory population 'genotype_example'
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
  define_genome_effect_terms("ADG", genotype_terms(cells, vals),
    effect_owner = "epistasis_table")
#> Wrote 4 genome-effect terms (8 member rows) for trait 'ADG' under owner 'epistasis_table' [append].
close_pop(pop)
```

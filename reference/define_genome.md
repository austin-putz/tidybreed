# Define the genome of a population

Lays out the loci of a population opened with
[`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md):
how many, on which chromosomes, at which physical and genetic positions.
It is the second step of every simulation, and runs once per population:

    pop <- open_pop(pop_name = "pigs", db_name = ":memory:") |>
      define_genome(n_loci = 50000, n_chr = 18, chr_len_Mb = 100)

Next come
[`define_chromosome()`](https://austin-putz.github.io/tidybreed/reference/define_chromosome.md)
for any sex chromosomes or organelles,
[`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)
for the founder allele pool, and
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md).

## Usage

``` r
define_genome(
  pop,
  n_loci,
  n_chr,
  chr_len_Mb,
  cM_per_Mb = 1,
  locus_names = NULL,
  chr_names = NULL,
  recombines_M = TRUE,
  recombines_F = TRUE
)
```

## Arguments

- pop:

  A `tidybreed_pop` from
  [`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md)
  with no genome yet; `define_genome()` is refused on a population that
  already has one.

- n_loci:

  Whole number: the total number of loci (SNPs) across the genome, at
  least `n_chr`.

- n_chr:

  Whole number: the number of chromosomes, including any sex chromosomes
  or organelles you will set up with
  [`define_chromosome()`](https://austin-putz.github.io/tidybreed/reference/define_chromosome.md).

- chr_len_Mb:

  Chromosome length(s) in megabases: one number for every chromosome, or
  a vector of `n_chr` values in chromosome order. Sets the physical
  positions, and with `cM_per_Mb` the genetic map. It does not change
  how many loci a chromosome gets.

- cM_per_Mb:

  Recombination rate in centiMorgans per megabase: one number, or a
  vector of `n_chr` values. Must be positive. Default `1`, about one
  crossover per 100 Mb per meiosis.

- locus_names:

  Optional character vector of `n_loci` unique, non-empty names, in
  genome order (chromosome 1 first). `NULL` (default) gives `"Locus_1"`,
  `"Locus_2"`, ...

- chr_names:

  Optional vector of `n_chr` unique, non-empty chromosome names, in
  order, e.g. `c(1:18, "X", "Y")` (stored as text). `NULL` (default)
  gives `"1"`, `"2"`, ... Name sex chromosomes and organelles here so
  [`define_chromosome()`](https://austin-putz.github.io/tidybreed/reference/define_chromosome.md)
  can refer to them.

- recombines_M, recombines_F:

  Logical: whether a male (`M`) or female (`F`) parent recombines its
  chromosomes when making gametes, genome-wide. Both default `TRUE`. Set
  one to `FALSE` for a species where one sex has no crossing over at
  all, e.g. `recombines_M = FALSE` for *Drosophila* males or
  `recombines_F = FALSE` for silkworm females. For single chromosomes (a
  Y or W), use
  [`define_chromosome()`](https://austin-putz.github.io/tidybreed/reference/define_chromosome.md)
  instead.

## Value

The `tidybreed_pop`, with the genome tables added.

## Details

**Where the loci go.** Loci are divided as evenly as possible across the
chromosomes, whatever their lengths (50 loci on 3 chromosomes gives 17,
16, 17), and spaced evenly within each chromosome, never at either end.
Each gets an id (`locus_id`, 1 to `n_loci` in genome order), a name, a
chromosome and a physical position `pos_bp`, stored in `genome_meta`. To
make some of them QTL or markers later, filter
`get_table(pop, "genome_meta")` and pipe it into
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md),
[`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md)
or
[`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md).

**Genetic map and recombination.** Each locus's genetic position is
`pos_cM = pos_bp / 1e6 * cM_per_Mb`, stored in `genome_map`. In
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
a chromosome's number of crossovers per meiosis is Poisson with mean
equal to its genetic length in Morgans (the last locus's
`pos_cM / 100`), placed uniformly in cM (the Haldane model). So the
default 1 cM/Mb gives a 100 Mb chromosome about one crossover per
meiosis.

**What is created.** `define_genome()` writes `genome_meta`,
`genome_map`, the seeded chromosome rules `chr_inheritance` and
`chr_recombination`, and the empty tables later functions fill:
`ind_haplotype` (alleles, written by
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
and
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)),
`ind_genotype` (dosages, filled on demand by
[`add_dosage()`](https://austin-putz.github.io/tidybreed/reference/add_dosage.md)),
`ind_crossover` and the three genome-effect tables with their views
`genome_effect_terms` and `genome_effect_loci`. Every chromosome starts
as a diploid autosome: one copy from each parent, recombining in both
sexes unless `recombines_M` / `recombines_F` say otherwise.

The call is **atomic**: everything is written in one transaction, so an
error (an invalid argument, loci too dense to place) leaves the database
as it was and you can simply call again with corrected arguments.

## See also

[`define_chromosome()`](https://austin-putz.github.io/tidybreed/reference/define_chromosome.md),
[`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md),
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md);
`get_table(pop, "genome_meta")` and `get_table(pop, "genome_map")` to
look at the result.

## Examples

``` r
# 1,000 loci on 10 chromosomes of 100 Mb each
pop <- open_pop(pop_name = "A", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp7yuhM6/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> Opened in-memory population 'A'
#> Defined genome: 1000 loci across 10 chromosomes | chr lengths (Mb): all equal to 100 Mb
#>   Tables written: genome_meta, genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_inheritance, chr_recombination, genome_effects, genome_effect_members, genome_effect_member_origins
#>   Views created: genome_effect_terms, genome_effect_loci
get_table(pop, "genome_meta")
#> <tidybreed_table: genome_meta>  [1000 rows × 5 fields]
#> # A tibble: 10 × 5
#>    locus_id locus_name   chr chr_name  pos_bp
#>       <int> <chr>      <int> <chr>      <dbl>
#>  1        1 Locus_1        1 1         990099
#>  2        2 Locus_2        1 1        1980198
#>  3        3 Locus_3        1 1        2970297
#>  4        4 Locus_4        1 1        3960396
#>  5        5 Locus_5        1 1        4950495
#>  6        6 Locus_6        1 1        5940594
#>  7        7 Locus_7        1 1        6930693
#>  8        8 Locus_8        1 1        7920792
#>  9        9 Locus_9        1 1        8910891
#> 10       10 Locus_10       1 1        9900990
get_table(pop, "genome_map")
#> <tidybreed_table: genome_map>  [1000 rows × 7 fields]
#> # A tibble: 10 × 7
#>    id_genome_map locus_id locus_name sex   line_name map_name pos_cM
#>            <int>    <int> <chr>      <chr> <chr>     <chr>     <dbl>
#>  1             1        1 Locus_1    NA    NA        default   0.990
#>  2             2        2 Locus_2    NA    NA        default   1.98 
#>  3             3        3 Locus_3    NA    NA        default   2.97 
#>  4             4        4 Locus_4    NA    NA        default   3.96 
#>  5             5        5 Locus_5    NA    NA        default   4.95 
#>  6             6        6 Locus_6    NA    NA        default   5.94 
#>  7             7        7 Locus_7    NA    NA        default   6.93 
#>  8             8        8 Locus_8    NA    NA        default   7.92 
#>  9             9        9 Locus_9    NA    NA        default   8.91 
#> 10            10       10 Locus_10   NA    NA        default   9.90 
close_pop(pop)

# Chromosomes of different lengths, a higher recombination rate and
# your own names (the first three cattle autosomes plus X)
pop <- open_pop(pop_name = "B", db_name = ":memory:") |>
  define_genome(n_loci = 400, n_chr = 4,
                chr_len_Mb = c(158, 137, 121, 139),
                cM_per_Mb  = 1.2,
                chr_names  = c("1", "2", "3", "X"))
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp7yuhM6/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> Opened in-memory population 'B'
#> Defined genome: 400 loci across 4 chromosomes | chr lengths (Mb): 1:158, 2:137, 3:121, X:139
#>   Tables written: genome_meta, genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_inheritance, chr_recombination, genome_effects, genome_effect_members, genome_effect_member_origins
#>   Views created: genome_effect_terms, genome_effect_loci
close_pop(pop)

# The usual next steps: a founder pool, then founders
set.seed(1)
pop <- open_pop(pop_name = "C", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 100,
                            min_allele_freq = 0.05,
                            max_allele_freq = 0.95)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp7yuhM6/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> Opened in-memory population 'C'
#> Defined genome: 1000 loci across 10 chromosomes | chr lengths (Mb): all equal to 100 Mb
#>   Tables written: genome_meta, genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_inheritance, chr_recombination, genome_effects, genome_effect_members, genome_effect_member_origins
#>   Views created: genome_effect_terms, genome_effect_loci
#> Generating 100 founder haplotypes (method = "uniform")...
#>   Created founder_haplotypes (100 haplotypes x 1000 loci)
pop <- pop |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 10, n_females = 10, line_name = "A")
#> Added 20 founders (10 males, 10 females) to line 'A'
close_pop(pop)

# No crossing over in males on any chromosome (Drosophila); a single
# non-recombining chromosome goes through define_chromosome() instead
pop <- open_pop(pop_name = "fly", db_name = ":memory:") |>
  define_genome(n_loci = 400, n_chr = 4, chr_len_Mb = c(23, 46, 50, 1.4),
                recombines_M = FALSE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp7yuhM6/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> Opened in-memory population 'fly'
#> Defined genome: 400 loci across 4 chromosomes | chr lengths (Mb): 1:23, 2:46, 3:50, 4:1.4
#>   Tables written: genome_meta, genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_inheritance, chr_recombination, genome_effects, genome_effect_members, genome_effect_member_origins
#>   Views created: genome_effect_terms, genome_effect_loci
get_table(pop, "chr_recombination")
#> <tidybreed_table: chr_recombination>  [8 rows × 4 fields]
#> # A tibble: 8 × 4
#>   chr_name parent_sex line_name recombines
#>   <chr>    <chr>      <chr>     <lgl>     
#> 1 1        M          NA        FALSE     
#> 2 1        F          NA        TRUE      
#> 3 2        M          NA        FALSE     
#> 4 2        F          NA        TRUE      
#> 5 3        M          NA        FALSE     
#> 6 3        F          NA        TRUE      
#> 7 4        M          NA        FALSE     
#> 8 4        F          NA        TRUE      
close_pop(pop)
```

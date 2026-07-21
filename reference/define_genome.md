# Define the genome structure of a breeding population

Adds genome tables (`genome_meta`, `genome_map`, `ind_haplotype`,
`ind_genotype`, `chr_meta`) to a population opened with
[`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md).
Pipe-friendly — accepts a `tidybreed_pop` and returns a `tidybreed_pop`.
Physical position (`pos_bp`, base pairs) lives in `genome_meta`; the
genetic map (`pos_cM`, centiMorgans) lives in the long `genome_map`
table. Haplotypes are stored in long format; `ind_genotype` is an
on-demand dosage cache populated by
[`add_dosage()`](https://austin-putz.github.io/tidybreed/reference/add_dosage.md).

This is the second step in setting up a simulation:

    pop <- open_pop(...) |>
      define_genome(n_loci = 50000, n_chr = 18, chr_len_Mb = 100)

To create the founder haplotype pool needed by
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md),
call
[`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)
after this function.

## Usage

``` r
define_genome(
  pop,
  n_loci,
  n_chr,
  chr_len_Mb,
  cM_per_Mb = 1,
  locus_names = NULL,
  chr_names = NULL
)
```

## Arguments

- pop:

  A `tidybreed_pop` object from
  [`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md).

- n_loci:

  Integer scalar. Total number of loci to simulate.

- n_chr:

  Integer scalar. Number of chromosomes.

- chr_len_Mb:

  Numeric scalar or numeric vector of length `n_chr`. Chromosome
  length(s) in megabases. A single scalar applies the same length to all
  chromosomes; a vector specifies each chromosome separately. Determines
  the physical coordinate `genome_meta.pos_bp`.

- cM_per_Mb:

  Numeric scalar or numeric vector of length `n_chr`. Genetic map rate
  in centiMorgans per megabase, used to derive the genetic map position
  `genome_map.pos_cM` from `pos_bp` (`pos_cM = pos_bp/1e6 * cM_per_Mb`).
  Default `1.0`. Must be finite and strictly positive.

- locus_names:

  Character vector of length `n_loci` or `NULL`. Custom locus names.
  When `NULL` (default), names are auto-generated as `"Locus_1"`,
  `"Locus_2"`, etc.

- chr_names:

  Character vector of length `n_chr` or `NULL`. Custom chromosome names.
  When `NULL` (default), chromosomes are numbered `1, 2, ..., n_chr`.

## Value

The input `tidybreed_pop` with genome tables added.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple setup
pop <- open_pop(pop_name = "A", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100)

# Different chromosome lengths (cattle)
pop <- open_pop(pop_name = "Cattle", db_name = ":memory:") |>
  define_genome(
    n_loci = 50000,
    n_chr  = 30,
    chr_len_Mb = c(158, 137, 121, 120, 121, 119, 112, 113, 105, 104,
                   107,  91,  84,  84,  85,  81,  75,  66,  64,  72,
                    71,  61,  52,  62,  42,  51,  45,  46,  51,  42)
  )

# With founder haplotypes
pop <- open_pop(pop_name = "B", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 100,
                            min_allele_freq = 0.05,
                            max_allele_freq = 0.95)
} # }
```

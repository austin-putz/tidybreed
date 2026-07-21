# Produce a gamete from a parent's two haplotypes via chromosomal crossovers

Simulates chromosomal recombination using the Haldane map function,
drawing from this gamete's own dqrng stream (seeded here from
`base_seed` + `sid`). Crossover count per chromosome ~
Poisson(chr_len_cM / 100), i.e. Poisson(genetic length in Morgans);
positions uniform in genetic (cM) distance. All chromosomes of the
gamete are drawn sequentially from the one seeded stream in
ascending-chromosome order.

## Usage

``` r
make_gamete(hap_matrix, chr_info, base_seed, sid, store_crossovers = FALSE)
```

## Arguments

- hap_matrix:

  2 x n_loci integer matrix. Row 1 = haplotype from parent_origin 1, row
  2 = parent_origin 2.

- chr_info:

  List from
  [`build_chr_info()`](https://austin-putz.github.io/tidybreed/reference/build_chr_info.md)
  — pre-computed locus indices, positions, and lengths per chromosome.

- base_seed:

  Integer base seed for the
  [`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
  call.

- sid:

  Integer stream id (see
  [`.gamete_stream_id()`](https://austin-putz.github.io/tidybreed/reference/dot-gamete_stream_id.md),
  kind = 1 autosome).

- store_crossovers:

  Logical; emit the crossover buffer when `TRUE`.

## Value

List: `allele` and `homolog` (length-`n_loci` integer vectors, the
gamete alleles and the donating homolog per locus), plus `xover_chr_idx`
(the `chr_info` position of each crossover event) and `xover_pos_cM`
(its cM position) — both `length(0)` unless `store_crossovers`.

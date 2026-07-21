# Pre-compute per-chromosome locus metadata for gamete simulation

Called once per
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
invocation. The returned list is passed to every
[`make_gamete()`](https://austin-putz.github.io/tidybreed/reference/make_gamete.md)
call, avoiding O(n_loci × n_chr) masking work inside the per-offspring
loop.

## Usage

``` r
build_chr_info(locus_map)
```

## Arguments

- locus_map:

  Data frame with columns `chr` and `pos_cM` — the resolved genetic map
  (from `resolve_genome_map()`), one row per locus, ordered by
  `locus_id`.

## Value

Named list, one element per chromosome. Each element is a list with
`locus_idx` (integer indices into the locus dimension), `pos_cM`
(genetic positions in centiMorgans), and `chr_len` (derived as
`MAX(pos_cM)` for the chromosome; crossovers beyond the last locus are
invisible to inheritance).

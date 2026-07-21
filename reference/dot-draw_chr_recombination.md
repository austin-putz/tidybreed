# Draw one chromosome's recombination pattern from the current dqrng stream

The per-chromosome core shared by autosomes
([`make_gamete()`](https://austin-putz.github.io/tidybreed/reference/make_gamete.md))
and recombining special chromosomes (the
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
special path). Assumes the gamete's stream is ALREADY seeded — draws
sequentially in the fixed canonical order (count, then start homolog,
then sorted positions) so R↔C++ parity holds. Purely a homolog-index
computation: it never touches allele values, so the caller gathers
alleles/`line_origin` from its own matrices via `hap_idx_vec`.

## Usage

``` r
.draw_chr_recombination(chr_pos, chr_len, store_crossovers = FALSE)
```

## Arguments

- chr_pos:

  Numeric locus genetic positions (cM) for this chromosome.

- chr_len:

  Numeric chromosome genetic length (cM) = `max(chr_pos)`.

- store_crossovers:

  Logical; when `TRUE` also return the drawn crossover positions (cM)
  for `ind_crossover`.

## Value

List with `hap_idx_vec` (length `length(chr_pos)`, values in `{1,2}` —
which homolog donates each locus) and `xover_pos_cM` (sorted crossover
positions in cM, or `numeric(0)`).

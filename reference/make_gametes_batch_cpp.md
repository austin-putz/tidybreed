# Autosome gamete kernel (C++), one resolved-map group

Bit-for-bit port of
[`make_gametes_batch_r()`](https://austin-putz.github.io/tidybreed/reference/make_gametes_batch_r.md).
Same arguments, same list return; see the R reference for the contract.
Not called directly — dispatched to by the internal selector
[`make_gametes_batch()`](https://austin-putz.github.io/tidybreed/reference/make_gametes_batch.md).

## Usage

``` r
make_gametes_batch_cpp(
  parent_allele,
  parent_lo_code,
  gamete_parent_idx,
  gamete_o,
  gamete_origin,
  chr_start,
  chr_end,
  chr_pos_cM,
  chr_len_cM,
  base_seed,
  store_crossovers
)
```

## Arguments

- parent_allele:

  Integer matrix (2\*n_parents) x n_autosome_loci.

- parent_lo_code:

  Integer matrix, same shape: line_origin codes.

- gamete_parent_idx:

  Integer length G: 1-based packed parent row.

- gamete_o:

  Integer length G: global offspring index per gamete.

- gamete_origin:

  Integer length G: parent_origin (1/2) per gamete.

- chr_start, chr_end:

  Integer per-chromosome 1-based inclusive ranges.

- chr_pos_cM:

  Numeric per-locus positions (cM), length n_autosome_loci.

- chr_len_cM:

  Numeric per-chromosome length (cM).

- base_seed:

  Integer base seed.

- store_crossovers:

  Logical.

## Value

Named list matching make_gametes_batch_r().

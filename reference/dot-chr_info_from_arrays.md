# Reconstruct the per-chromosome `chr_info` list from flat kernel arrays

The C++-ready seam
([`make_gametes_batch_r()`](https://austin-putz.github.io/tidybreed/reference/make_gametes_batch_r.md))
receives the autosome map as flat arrays
(`chr_start`/`chr_end`/`chr_pos_cM`/`chr_len_cM`) rather than the nested
`chr_info` list
[`make_gamete()`](https://austin-putz.github.io/tidybreed/reference/make_gamete.md)
consumes. This rebuilds that list so the unchanged per-gamete core can
be reused verbatim. `chr_start`/`chr_end` are 1-based inclusive
locus-index ranges into `1..n_autosome_loci`.

## Usage

``` r
.chr_info_from_arrays(chr_start, chr_end, chr_pos_cM, chr_len_cM)
```

## Arguments

- chr_start, chr_end:

  Integer per-chromosome locus-index ranges.

- chr_pos_cM:

  Numeric per-locus genetic positions (cM), length `n_autosome_loci`.

- chr_len_cM:

  Numeric per-chromosome genetic length (cM).

## Value

A `chr_info`-shaped list (one element per chromosome).

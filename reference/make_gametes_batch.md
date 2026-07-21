# Dispatch the autosome gamete kernel to the C++ or R implementation

The single call site
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
uses. Prefers the compiled kernel
[`make_gametes_batch_cpp()`](https://austin-putz.github.io/tidybreed/reference/make_gametes_batch_cpp.md)
and falls back to the pure-R reference / parity oracle
[`make_gametes_batch_r()`](https://austin-putz.github.io/tidybreed/reference/make_gametes_batch_r.md)
when the compiled object is unavailable or explicitly disabled. Both
take byte-identical inputs and return byte-identical output for a given
`base_seed` (the R↔C++ parity contract), so the choice never changes
results.

## Usage

``` r
make_gametes_batch(
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
  store_crossovers = FALSE
)
```

## Arguments

- parent_allele:

  Integer matrix `(2 * n_parents) x n_autosome_loci`, homolog-major:
  rows `2p-1` and `2p` are parent `p`'s homolog 1 and 2.

- parent_lo_code:

  Integer matrix, same shape: `line_origin` dictionary codes (the caller
  decodes back to `VARCHAR`).

- gamete_parent_idx:

  Integer length `G`: 1-based packed parent row that produces each
  gamete.

- gamete_o:

  Integer length `G`: global offspring index `o` per gamete (original
  `matings` row); drives the per-gamete stream.

- gamete_origin:

  Integer length `G`: `parent_origin` (1 = parent_1/sire, 2 =
  parent_2/dam) per gamete.

- chr_start, chr_end:

  Integer per-chromosome locus-index ranges (1-based inclusive,
  contiguous, into `1..n_autosome_loci`) for THIS group's map.

- chr_pos_cM:

  Numeric per-locus genetic positions (cM), length `n_autosome_loci`,
  for this group's map.

- chr_len_cM:

  Numeric per-chromosome genetic length (cM) for this group.

- base_seed:

  Integer base seed for the
  [`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
  call.

- store_crossovers:

  Logical; emit the ragged crossover buffer when `TRUE`.

## Value

The list returned by the selected kernel.

## Details

Selection: `getOption("tidybreed.kernel")`, else the `TIDYBREED_KERNEL`
environment variable, else `"auto"`. `"r"` forces the R reference (for
debugging / no-compiler environments); anything else uses C++ when
available.

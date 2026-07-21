# Generate autosome gametes for one resolved-map group (the swappable seam)

The single seam the C++ kernel (Stage 3) drops into: a dependency-free
numeric function over **packed integer arrays** — no `tidybreed_pop`, no
DBI, no `data.frame`, no character `line_origin`, no string-keyed lists
— so the R reference and the C++ kernel take byte-identical inputs and
their outputs can be compared directly for parity. This is the R
reference / parity oracle;
[`make_gametes_batch()`](https://austin-putz.github.io/tidybreed/reference/make_gametes_batch.md)
dispatches to it (or to the compiled kernel).

## Usage

``` r
make_gametes_batch_r(
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

List of long parallel vectors over the `G` gametes in the given order
(gamete -\> locus ascending): `parent_origin`, `locus_idx`
(1..n_autosome_loci), `allele`, `line_origin_code`; plus the ragged
crossover buffer `xover_gamete_o` (global `o`), `xover_parent_origin`,
`xover_chr` (index into the autosome chromosome order), `xover_pos_cM`
(length 0 unless `store_crossovers`).

## Details

**Gamete-flat / one-map-per-call (group-by-map seam).** Each gamete may
use a different genetic map (sex-specific, achiasmy, line-specific), so
the caller groups gametes by resolved-map key and calls this once per
group with that group's `chr_*` arrays. The unit is a **gamete**,
identified by `(gamete_o, gamete_origin)`; each draws from its own dqrng
stream keyed on its **global** offspring index `gamete_o` (autosome kind
= 1), so output is independent of batch size, iteration order, and how
gametes are partitioned into map-groups or calls.

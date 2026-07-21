# Compute per-locus allele frequencies from the base population

Compute per-locus allele frequencies from the base population

## Usage

``` r
compute_base_allele_freq(pop, base, base_ids = NULL)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- base:

  Character. `"founder_haplotypes"` or `"current_pop"`.

- base_ids:

  Optional character vector of `id_ind` for `"current_pop"`.

## Value

Numeric vector of allele frequencies, length `n_loci`, in `locus_id`
order.

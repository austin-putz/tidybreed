# Assert that a set of individuals all have ploidy = 2

Ploidy-aware dosage and genotype extraction for real polyploidy (ploidy
\> 2) is not yet supported. Sex-linked and organelle `chr_inheritance`
rules (reduced `from_parent_1`/`from_parent_2` copy counts) are fine
here — `SUM(allele)` over however many `ind_haplotype` rows exist per
locus already produces correct dosage regardless of row count, so this
guard targets actual organism ploidy, not chromosome copy-count shape.

## Usage

``` r
assert_ploidy_2(pop, ids = NULL)
```

## Arguments

- pop:

  A `tidybreed_pop`.

- ids:

  Optional character vector of `id_ind` to check; `NULL` (default)
  checks every individual in `ind_meta`.

## Value

Invisible `NULL` on success; errors otherwise.

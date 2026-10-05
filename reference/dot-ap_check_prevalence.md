# Refuse a prevalence threshold that has no genetic variance to use

Runs
[`.ap_prevalence_genetic_var()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_prevalence_genetic_var.md)
for every categorical phenotype placed by `prevalence`, in Stage 1,
before any `ind_tgv` write or draw. A composite phenotype is refused
outright: its liability combines several traits and contributors, which
no stored diagonal describes.

## Usage

``` r
.ap_check_prevalence(pop, pheno_meta, composite)
```

## Arguments

- pheno_meta:

  The call's `phenotype_meta` rows.

- composite:

  Logical, per row: has `phenotype_components` or `formula_tgv`.

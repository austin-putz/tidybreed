# Refuse a prevalence threshold that has no genetic variance to use

A categorical phenotype defined by `prevalence` gets its liability
threshold from the stored additive target of its trait. A composite
phenotype has no such target, and neither does a simple trait whose
effects were never calibrated to one; both used to fall back silently to
zero genetic variance.

## Usage

``` r
.ap_check_prevalence(pop, pheno_meta, composite)
```

## Arguments

- pheno_meta:

  The call's `phenotype_meta` rows.

- composite:

  Logical, per row: has `phenotype_components` or `formula_tbv`.

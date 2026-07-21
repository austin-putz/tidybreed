# Topologically sort phenotypes for safe evaluation order.

Derived formula phenotypes depend on other phenotypes being present in
ind_phenotype first. formula_tbv and components phenotypes have no
inter-phenotype dependencies (they depend on trait_meta, not
phenotype_meta).

## Usage

``` r
.topo_sort_phenotypes(pheno_meta)
```

## Arguments

- pheno_meta:

  Data frame from phenotype_meta (all columns).

## Value

Character vector: phenotype_name in safe evaluation order.

## Details

Uses Kahn's BFS algorithm. Detects cycles and stops with an informative
error.

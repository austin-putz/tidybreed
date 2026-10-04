# Genetic value of a trait per id (`NA` for a `NA` id or no `ind_tgv` row)

Genetic value of a trait per id (`NA` for a `NA` id or no `ind_tgv` row)

## Usage

``` r
.tgv_by_id(conn, trait_name, ids, components = "total")
```

## Arguments

- components:

  `"total"` (the `ind_tgv_total` view, the default for every phenotype
  path) or a character vector of `ind_tgv.component_name` values, summed
  exactly; a listed component the individual has no row for
  contributes 0. See `.tgv_read()`.

# Assemble the composite genetic value of one phenotype from `phenotype_components`

Sums `weight * contributor genetic value` over the phenotype's component
rows, one contributor lookup per row (see
[`?contributor_tgv`](https://austin-putz.github.io/tidybreed/reference/contributor_tgv.md)).
Each row reads its `component_names` (`"total"` by default) from
`ind_tgv`, which must already hold the source traits for every
contributor
([`.ap_materialize_tgvs()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_materialize_tgvs.md)).
A missing piece — a `NULL` dam or sire, a contributor with no `ind_tgv`
row, a `NULL` group value, a `NULL` covariate — makes the individual's
composite `NA`.

## Usage

``` r
.assemble_composite_tgv(
  pop,
  phenotype_name,
  comp_rows,
  subset_df,
  missing_component_action
)
```

## Arguments

- comp_rows:

  The phenotype's `phenotype_components` rows.

- subset_df:

  The planned `ind_meta` rows (needs `id_ind`, `id_parent_1`,
  `id_parent_2`).

- missing_component_action:

  `"skip"` warns and returns `NA` for the excluded individuals;
  `"error"` stops. Both name the count and up to five ids.

## Value

Numeric vector named by `id_ind`; `NA` marks an excluded individual.

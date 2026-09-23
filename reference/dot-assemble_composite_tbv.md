# Assemble the composite TBV of one phenotype from `phenotype_components`

Sums `weight * contributor TBV` over the phenotype's component rows, one
contributor lookup per row (see
[`?contributor_tbv`](https://austin-putz.github.io/tidybreed/reference/contributor_tbv.md)).
`ind_tbv` must already hold the source traits for every contributor
([`.ap_materialize_tbvs()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_materialize_tbvs.md)).
A missing piece — a `NULL` dam or sire, a contributor with no TBV, a
`NULL` group value, a `NULL` covariate — makes the individual's
composite `NA`.

## Usage

``` r
.assemble_composite_tbv(
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

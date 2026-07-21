# Assemble composite TBV from phenotype_components for a single phenotype

Reads contributor TBVs from `ind_tbv` (which must already be populated
by
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
for all relevant source traits and contributor IDs) and multiplies by
component weights, summing across all components.

## Usage

``` r
.assemble_composite_tbv(
  pop,
  phenotype_name,
  comp_rows,
  subset_df,
  missing_component_action = "skip"
)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- phenotype_name:

  Character. The composite phenotype name.

- comp_rows:

  Data frame. Rows from `phenotype_components` for this phenotype.

- subset_df:

  Data frame. Sex-filtered (and skip-masked) `ind_meta` rows.

- missing_component_action:

  Character. `"skip"` (default) to warn and exclude individuals with any
  missing component (dam/sire/group TBV unavailable, no group
  assignment, etc.); `"error"` to stop immediately.

## Value

A list with `composite_tbv`: named numeric vector (NA = excluded).

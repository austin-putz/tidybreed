# Fetch and aggregate group-mate TBVs for a group_sum / group_mean term.

Validates that grp_table and grp_col exist (errors with clear message).
Singleton pens (no group-mates) receive 0. Individuals with NA group get
NA.

## Usage

``` r
.fetch_group_tbvs(
  pop,
  trait_safe,
  focal_ids,
  grp_col,
  grp_table,
  agg_fn,
  phenotype_name
)
```

## Arguments

- pop:

  A tidybreed_pop object.

- trait_safe:

  Character. Trait name, already SQL-escaped (single quotes doubled).

- focal_ids:

  Character vector of focal individual IDs.

- grp_col:

  Character. Column in `grp_table` holding group membership (e.g. pen or
  litter ID).

- grp_table:

  Character. Table containing `grp_col` and `id_ind` (default
  `"ind_meta"` at the call site).

- agg_fn:

  Function. `sum` or `mean`, applied to group-mates' TBVs (excluding the
  focal individual itself).

- phenotype_name:

  Character. Name of the phenotype being computed; used only in error
  messages.

## Value

Numeric vector, same length as `focal_ids`: the aggregated group-mate
TBV (0 for singletons with no group-mates), or `NA` for individuals with
no group assignment.

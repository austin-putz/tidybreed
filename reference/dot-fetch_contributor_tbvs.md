# Fetch TBV values for a set of contributor IDs, aligned to focal_ids. Returns NA where contributor is absent or has no TBV record.

Fetch TBV values for a set of contributor IDs, aligned to focal_ids.
Returns NA where contributor is absent or has no TBV record.

## Usage

``` r
.fetch_contributor_tbvs(pop, trait_safe, contr_ids, focal_ids)
```

## Arguments

- pop:

  A tidybreed_pop object.

- trait_safe:

  Character. Trait name, already SQL-escaped (single quotes doubled).

- contr_ids:

  Character vector of contributor IDs (e.g. `id_parent_1` or
  `id_parent_2`), same length and order as `focal_ids`.
  `NA`/`"NA"`/empty values are treated as missing.

- focal_ids:

  Character vector of focal individual IDs.

## Value

Numeric vector, same length as `focal_ids`, with the contributor's
`tbv_value` or `NA` where the contributor is missing or has no TBV.

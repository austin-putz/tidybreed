# Aggregated group-mate TBV per focal individual

Aggregated group-mate TBV per focal individual

## Usage

``` r
.group_mate_tbv(
  conn,
  trait_name,
  focal_ids,
  group_column,
  group_table,
  aggregation,
  what
)
```

## Arguments

- aggregation:

  `"sum"` or `"mean"` over the mates that have a TBV.

## Value

Numeric per focal: the aggregate, `0` with no such mates, `NA` with no
group value.

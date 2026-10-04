# Aggregated group-mate genetic value per focal individual

Aggregated group-mate genetic value per focal individual

## Usage

``` r
.group_mate_tgv(
  conn,
  trait_name,
  focal_ids,
  group_column,
  group_table,
  aggregation,
  what,
  components = "total"
)
```

## Arguments

- aggregation:

  `"sum"` or `"mean"` over the mates that have a genetic value for the
  trait.

- components:

  As in
  [`.tgv_by_id()`](https://austin-putz.github.io/tidybreed/reference/dot-tgv_by_id.md).

## Value

Numeric per focal: the aggregate, `0` with no such mates, `NA` with no
group value.

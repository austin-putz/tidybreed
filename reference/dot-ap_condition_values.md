# The residual condition value of every planned record

[`.read_one_per_id()`](https://austin-putz.github.io/tidybreed/reference/dot-read_one_per_id.md)
on the condition table (exactly one row per planned id), returned as
character (`NA` for `NULL`) to match how
`phenotype_var_comp.condition_level` is stored.

## Usage

``` r
.ap_condition_values(conn, condition_table, condition_column, ids)
```

# Reserved columns per table (cannot be modified by the user)

A table absent from this list is unprotected:
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
will happily overwrite any column on it. Tables written and owned
exclusively by a `define_*`/`add_*` function (`founder_haplotypes`,
`phenotype_components`, `phenotype_random_effects`, `genome_effects`,
...) therefore list **every** column. Only entity-shaped tables that
users legitimately annotate — `ind_meta`, `genome_meta`,
`ind_phenotype`, `index_meta` — leave room for user columns.

## Usage

``` r
TABLE_RESERVED_COLS
```

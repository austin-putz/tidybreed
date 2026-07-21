# S3 class wrapping a table reference for generic mutation

Returned by
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md).
Carries the parent population object, the table name, a lazy dplyr
tibble for reading, and a pending filter list. Supports
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html),
[`collect()`](https://dplyr.tidyverse.org/reference/compute.html), and
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md).

## Usage

``` r
new_tidybreed_table(pop, table_name)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- table_name:

  Character scalar. Name of the table to reference.

## Value

A `tidybreed_table` S3 object.

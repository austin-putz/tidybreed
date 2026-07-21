# Get a table reference from a tidybreed population

Returns a `tidybreed_table` object that can be piped into
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
and
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md),
or collected with
[`dplyr::collect()`](https://dplyr.tidyverse.org/reference/compute.html).
All existing
`get_table(pop, x) |> dplyr::filter(...) |> dplyr::collect()` patterns
continue to work unchanged.

## Usage

``` r
get_table(pop, table_name)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- table_name:

  Name of the table to retrieve.

## Value

A `tidybreed_table` object.

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- open_pop(pop_name = "A", db_name = ":memory:") |>
  define_genome(n_loci = 100, n_chr = 2, chr_len_Mb = 100)

# Read-only query
get_table(pop, "genome_meta") |> dplyr::collect()

# Mutate all rows
pop <- pop |> get_table("ind_meta") |> mutate_table(gen = 1L)

# Mutate a filtered subset
pop <- pop |>
  get_table("ind_meta") |>
  filter(sex == "M") |>
  mutate_table(gen = 2L)
} # }
```

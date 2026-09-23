# Mark animals as genotyped on a SNP chip

Records which individuals have been genotyped on a named SNP chip by
writing or updating a BOOLEAN column `has_<chip_name>` in `ind_meta`.
Pipe a `tidybreed_table` (from
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
and optionally
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html))
as the first argument to restrict which animals are marked.

The operation is **additive**: animals already marked `TRUE` remain
`TRUE`. Only new animals are flipped. This mirrors real life — once an
animal is genotyped it stays genotyped.

## Usage

``` r
add_genotypes(tbl, chip_name, col_name = paste0("has_", chip_name))
```

## Arguments

- tbl:

  A `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md),
  optionally piped through
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html).
  Any table with an `id_ind` column is accepted; the individuals acted
  on are the distinct `id_ind` values present in the (filtered) table.
  An unfiltered `ind_meta` selects every individual; an unfiltered
  `ind_ebv`, `ind_index`, `ind_genotype`, ... selects only the
  individuals that have rows there. A table without `id_ind` is an
  error.

- chip_name:

  Character. Name of an existing SNP chip (must have an `is_<chip_name>`
  column in `genome_meta`, created by
  [`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md)).

- col_name:

  Character. Name of the BOOLEAN column to write in `ind_meta`. Default:
  `paste0("has_", chip_name)`.

## Value

The modified `tidybreed_pop` (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
# Genotype all females in generation 1 on the 50k chip
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(sex == "F", gen == 1L) |>
  add_genotypes("50k")

# Also genotype all generation 2 animals (additive — gen 1 females stay TRUE)
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(gen == 2L) |>
  add_genotypes("50k")
} # }
```

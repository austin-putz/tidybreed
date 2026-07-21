# Define a SNP chip

Marks the loci in a filtered `genome_meta` table as members of a SNP
chip, creating a `BOOLEAN` column (default: `is_{chip_name}`) in
`genome_meta`.

Pipe `get_table("genome_meta")` — optionally piped through
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html) —
into this function. The unique `locus_id` values in the collected result
determine chip membership. All other loci receive `FALSE`.

This replaces the old positional selection methods (`n`/`method`,
`locus_tf`, `locus_ids`, `locus_names`). Selection is now done entirely
by the caller using standard dplyr verbs on the `genome_meta` table,
which is order-safe and compatible with dynamic genomes.

## Usage

``` r
define_chip(tbl, chip_name, col_name = paste0("is_", chip_name))
```

## Arguments

- tbl:

  A `tidybreed_table` from `get_table("genome_meta")` (optionally piped
  through
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)).
  Must contain a `locus_id` column.

- chip_name:

  Character scalar. Name of the SNP chip; used in messages and to derive
  the default column name.

- col_name:

  Character scalar. Column name created in `genome_meta`. Default:
  `paste0("is_", chip_name)`. Must be a valid SQL identifier.

## Value

The modified `tidybreed_pop` object (invisibly).

## See also

[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md),
[`add_genotypes()`](https://austin-putz.github.io/tidybreed/reference/add_genotypes.md),
[`extract_genotypes()`](https://austin-putz.github.io/tidybreed/reference/extract_genotypes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# All loci on chromosomes 1-5
pop <- pop |>
  get_table("genome_meta") |>
  dplyr::filter(chr %in% 1:5) |>
  define_chip("chr1to5")

# Random chip — sample locus names, then filter
selected <- pop |>
  get_table("genome_meta") |>
  dplyr::collect() |>
  dplyr::slice_sample(n = 500) |>
  dplyr::pull(locus_name)

pop <- pop |>
  get_table("genome_meta") |>
  dplyr::filter(locus_name %in% selected) |>
  define_chip("50K")

# Complement of an existing chip
pop <- pop |>
  get_table("genome_meta") |>
  dplyr::filter(is_50K == FALSE) |>
  define_chip("non50K")
} # }
```

# Extract genotype data for individuals, by chip and/or QTL loci

Returns a tibble of genotypes (0/1/2 encoding) for a set of individuals,
restricted to loci selected by a chip definition, a filtered
`genome_effects` table, or both. Pipe a `tidybreed_table` (from
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
and optionally
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html))
as the first argument to restrict which individuals are included.

At least one of `chip_name` or `effects_tbl` must be provided. When both
are given the locus sets are **unioned** (deduplicated, ordered by
`locus_id`).

**Chip path** — the returned individual set is the intersection of:

- Animals with `has_<chip_name> == TRUE` in `ind_meta`

- Animals matching any pending
  [`filter()`](https://rdrr.io/r/stats/filter.html) predicates on `tbl`

- Loci with `is_<chip_name> == TRUE` in `genome_meta`

**QTL path (`effects_tbl`)** — the returned individual set is:

- Animals matching any pending
  [`filter()`](https://rdrr.io/r/stats/filter.html) predicates on `tbl`
  (or all individuals in `ind_haplotype` when no filter is applied)

- Loci whose `locus_name` appears in the collected `effects_tbl`

## Usage

``` r
extract_genotypes(
  tbl,
  chip_name = NULL,
  effects_tbl = NULL,
  loci_tbl = NULL,
  col_name = NULL
)
```

## Arguments

- tbl:

  A `tidybreed_table` object from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (optionally piped through
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)).
  The table must contain an `id_ind` column when a filter is applied.

- chip_name:

  Character or `NULL`. Name of a chip previously defined via
  [`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md)
  and applied to animals via
  [`add_genotypes()`](https://austin-putz.github.io/tidybreed/reference/add_genotypes.md).
  When `NULL` the chip path is skipped.

- effects_tbl:

  A `tidybreed_table` from `get_table(pop, "genome_effects")`
  (optionally filtered), or `NULL`. The collected table must contain a
  `locus_name` column. Use
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
  to restrict by `trait_name`, `genome_effect_type`, `genome_value`,
  `line_name`, etc. When `NULL` the QTL path is skipped.

- loci_tbl:

  A `tidybreed_table` from `get_table(pop, "genome_meta")` (optionally
  filtered), or `NULL`. A general locus filter independent of chips and
  QTL sets — e.g. `filter(!chr_name %in% c("X", "Y", "MT"))` to restrict
  to autosomes. The collected table's `locus_name` values become the
  locus set. Unioned with `chip_name`/`effects_tbl` when combined.

- col_name:

  Character. Name of the BOOLEAN column in `ind_meta` that records chip
  genotyping status. Default: `paste0("has_", chip_name)`. Ignored when
  `chip_name` is `NULL`.

## Value

A tibble with one row per individual and one column per selected locus
(`id_ind` + `locus_N` columns). Locus columns use 0/1/2 integer encoding
and are ordered by `locus_id`.

## Examples

``` r
if (FALSE) { # \dontrun{
# All genotyped animals on the 50k chip (unchanged usage)
geno <- pop |> get_table("ind_meta") |> extract_genotypes("50k")

# Only females genotyped on the HD chip
geno <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(sex == "F") |>
  extract_genotypes("HD")

# QTL loci for a trait (all individuals)
geno <- pop |>
  get_table("ind_meta") |>
  extract_genotypes(
    effects_tbl = get_table(pop, "genome_effects") |>
      dplyr::filter(trait_name == "ADG")
  )

# Large-effect QTL only, females only
geno <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(sex == "F") |>
  extract_genotypes(
    effects_tbl = get_table(pop, "genome_effects") |>
      dplyr::filter(trait_name == "ADG", abs(genome_value) > 0.15)
  )

# Chip loci + QTL loci unioned
geno <- pop |>
  get_table("ind_meta") |>
  extract_genotypes(
    chip_name   = "50k",
    effects_tbl = get_table(pop, "genome_effects") |>
      dplyr::filter(trait_name == "ADG")
  )
} # }
```

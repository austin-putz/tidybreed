# Compute and store true breeding values without writing phenotypes

Computes the true breeding value (TBV) for each individual in the
current subset and each requested trait, and writes them to `ind_tbv`.
This is the exact function
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
calls internally (once for every source trait it needs) before
assembling phenotype records — there is no separate "TBV math"
duplicated elsewhere.

TBV is the Falconer-centered sum, across every `ind_haplotype` row (one
per allele copy, not genotype dosage) for the individual, of:


      TBV_i = sum over haplotype rows of (allele - base_allele_freq) * genome_value

`genome_value` and `base_allele_freq` are read from `genome_effects`
(`genome_effect_type = "additive"`). For each haplotype row, a
**line-specific** effect (`genome_effects.line_name` matching that row's
`line_origin`) is preferred; the **population-wide** effect
(`genome_effects.line_name IS NULL`) is used only when no line-specific
row exists for that locus/line. This per-locus fallback is what makes
crossbreeding TBV correct — e.g. a Duroc x Landrace F1 is centered
against each parent line's own QTL effects and base allele frequency
(see the "Crossbreeding TBV" example below). For **imprinted** traits
(`trait_meta.expressed_parent` = `"parent_1"` or `"parent_2"`), only
haplotype rows from that parent's `parent_origin` are summed before the
same line-matching logic applies.

Optionally computes true selection index values by multiplying per-trait
TBVs by weights from named indices defined with
[`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md),
and writes them to `ind_true_index`.

Pipe a `tidybreed_table` (from
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
and optionally
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html))
as the first argument to select individuals. The `expressed_sex` rule
from `trait_meta` is applied on top.

Useful for tracking genetic trend across generations without collecting
phenotypes.

## Usage

``` r
add_tbv(
  tbl,
  trait_name = NULL,
  index_names = NULL,
  type = c("index", "economic", "both"),
  overwrite_index = FALSE,
  ...
)
```

## Arguments

- tbl:

  A `tidybreed_table` object from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (optionally piped through
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)).
  The table must contain an `id_ind` column.

- trait_name:

  Character vector of trait name(s). When `NULL` (default), all traits
  currently in `trait_meta` are used (in `id_trait` order).

- index_names:

  Character vector of named index(es) from `index_meta` for which true
  index values should be computed from TBVs and written to
  `ind_true_index`. When `NULL` (default), no true index computation is
  performed. All index traits must be included in `trait_name` (or all
  traits when `trait_name = NULL`).

- type:

  Which weight column from `index_meta` to use: `"index"` uses
  `index_weight`, `"economic"` uses `economic_weight`, `"both"` computes
  and stores both (distinguished by the `weight_type` column in
  `ind_true_index`). Defaults to `"index"`.

- overwrite_index:

  Logical. When `FALSE` (default), individuals that already have a true
  index value in `ind_true_index` for the given
  `(index_name, weight_type)` combination are skipped — avoids redundant
  recomputation across generations. When `TRUE`, existing rows are
  deleted and recomputed (use when index weights have changed).

- ...:

  Optional extra columns written to `ind_tbv` (scalars only; broadcast
  to all records).

## Value

The modified `tidybreed_pop` (invisibly).

## See also

[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
[`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md),
[`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Crossbreeding TBV: line-specific additive effects for two pure lines, then
# a Duroc x Landrace F1 centered against each parent line's own effects and
# base allele frequency (see the Description above for the matching rule)
pop <- pop |>
  get_table("genome_meta") |>
  define_additive_effects("ADG", effects = duroc_effects, line_name = "Duroc")
pop <- pop |>
  get_table("genome_meta") |>
  define_additive_effects("ADG", effects = landrace_effects, line_name = "Landrace")
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(line_name == "F1") |>
  add_tbv("ADG")

# TBVs only, for a generation subset
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(gen == 2L) |>
  add_tbv(c("ADG", "BW"))

# TBVs + true index values (both index and economic weights) written to
# ind_true_index, distinguished by weight_type
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(gen == 2L) |>
  add_tbv(c("ADG", "BW"), index_names = "terminal", type = "both")
} # }
```

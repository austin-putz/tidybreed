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


      TBV_i = sum over allele copies of (allele - center_value) * genome_value

`genome_value` and `center_value` come from the **order-one `additive`
terms** under the reserved effect owner `generated_additive_tbv`, the
terms
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
writes. `center_value` is that variant's base allele frequency.

**This is one filtered call into the same evaluator
[`add_tgv()`](https://austin-putz.github.io/tidybreed/reference/add_tgv.md)
uses**, not a second implementation: `add_tbv()` is
[`add_tgv()`](https://austin-putz.github.io/tidybreed/reference/add_tgv.md)
restricted to the reserved owner and to single-member additive terms.

The filter is deliberate and not merely conservative. Under functional
\\(a, d)\\ input the stored coefficient is \\a\\, while the
breeding-value coefficient in a diploid HWE base is \\\alpha = a + d(q -
p)\\; under epistasis, average effects depend on other loci and on LD.
So arbitrary terms written through
[`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md)
contribute to `ind_tgv` but **never silently redefine the breeding
value**, additive members appearing inside interactions are ignored, and
`ind_tbv` keeps its exact meaning. Deriving average effects from a
general non-additive model is a separate calculation.

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

## When the stored coefficients stop being average effects

Ignoring those terms is right, but it stops giving *the model's*
breeding value as soon as one of them contributes to the additive
component or shifts the coefficients this function reads. `add_tbv()`
warns once per trait in exactly that case: a non-reserved order-one
`additive` term (it is part of A and is skipped), an `indicator` surface
(raw functional coding — at a locus that also carries a generated
additive term the stored `a` is no longer the average effect, \\\alpha =
a + d(q - p)\\), or an interaction (whose additive projection depends on
other loci and on LD, so there is no local correction).

An order-one `dominance` term centred where the additive term is centred
is the **exception and stays silent**: Cockerham coding is
HWE-orthogonal, so it contributes nothing to A and leaves the additive
coefficient alone — `tbv_value` is still exact. Warning there would cry
wolf on the common case.

Each allele copy takes the **most specific** variant whose origin
predicate matches its `(line_origin, parent_origin)` label, falling back
per copy to the common variant. This per-copy fallback is what makes
crossbreeding TBV correct — e.g. a Duroc x Landrace F1 is centered
against each parent line's own effects and base allele frequency (see
the "Crossbreeding TBV" example below). **Imprinting** is a property of
the effect, not of the trait: a term scoped to one `parent_origin` (see
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md))
reads only that parent's allele copies, per locus and per line.

Optionally computes true selection index values by multiplying per-trait
TBVs by weights from named indices defined with
[`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md),
and writes them to `ind_true_index`.

Pipe a `tidybreed_table` (from
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
and optionally
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html))
as the first argument to select individuals. Every individual in that
subset receives a TBV for every requested trait — unlike
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
no sex-expression rule is applied here (`expressed_sex` is an
observation-layer property of `phenotype_meta`, not of a genetic
component trait).

Useful for tracking genetic trend across generations without collecting
phenotypes.

## See also

[`add_tgv()`](https://austin-putz.github.io/tidybreed/reference/add_tgv.md)
for every component of the genetic value,
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

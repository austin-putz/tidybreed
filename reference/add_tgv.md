# Compute and store true genetic values

Evaluates the stored effect model for each individual in the current
subset and writes the result to `ind_tgv`, one row per (individual x
trait x **component**). The total is the derived view `ind_tgv_total`,
never a stored row — a stored total would make every `SUM(tgv_value)`
double-count.

Unlike
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md),
this evaluates **every** term of the trait: additive, dominance,
hand-entered genotype surfaces and multi-locus interactions, under every
effect owner. `component_name` records how each term was *declared*:

|  |  |
|----|----|
| `component_name` | Terms it collects |
| `"order1_additive"` | single-member terms whose contrast is `additive` |
| `"order1_dominance"` | single-member `dominance` terms |
| `"order1_other"` | single-member `indicator` terms (a hand-entered surface) |
| `"interaction"` | any term with two or more members |

These are **model-structure components, not variance components.** A
functional A x A term contributes to \\V_A\\, \\V_D\\ *and* \\V_I\\ in
the statistical sense; the names carry the declared order precisely so
they cannot be read as an orthogonal decomposition.

`ind_tgv` stores the **raw sum of the stored terms — no mean is added.**
A pure Cockerham model yields centered deviations; a raw `indicator`
surface yields absolute genotypic values with a non-zero mean by
construction, and is not silently re-centered. See
[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md),
which reports the implied genetic mean and writes it nowhere.

Writes are idempotent: an individual's rows for a trait are replaced, so
re-evaluating after changing the effect model never leaves stale
components behind.

## Usage

``` r
add_tgv(tbl, trait_name = NULL)
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
  in `trait_meta` are used (in `id_trait` order).

## Value

The modified `tidybreed_pop` (invisibly).

## Resource guard

Evaluation enumerates **label-vectors** — one distinct
`(line_origin, parent_origin)` label per member — not allele copies per
individual, so its cost is a function of the model rather than of
population size. The count is estimated per fallback family before any
work runs; it warns above
`getOption("tidybreed.label_vector_warn", 1e4)` and stops above
`getOption("tidybreed.label_vector_max", 1e6)`, so an accidental
high-order scoped term fails loudly instead of appearing to hang.

## See also

[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
for the breeding value,
[`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md)
and
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
for writing the terms this evaluates.

## Examples

``` r
if (FALSE) { # \dontrun{
# Every component of the model, for one generation
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(gen == 2L) |>
  add_tgv("ADG")

# The components, and the derived total
pop |> get_table("ind_tgv") |> dplyr::collect()
pop |> get_table("ind_tgv_total") |> dplyr::collect()
} # }
```

# Compute and store true genetic values

Evaluates the stored effect model for each individual in the current
subset and writes the result to `ind_tgv`, one row per (individual x
trait x **component**). The total is the derived view `ind_tgv_total`,
never a stored row — a stored total would make every `SUM(tgv_value)`
double-count. This is the one table of true genetic values:
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
calls this function for every trait a phenotype reads, and reads the
total (or the components `phenotype_components.component_names` lists)
from it.

Every term of the trait is evaluated: additive, dominance, hand-entered
genotype surfaces and multi-locus interactions, under every effect
owner. `component_name` records how each term was declared:

|  |  |
|----|----|
| `component_name` | Terms it collects |
| `"additive"` | one-locus terms whose contrast is `additive` |
| `"dominance"` | one-locus `dominance` terms |
| `"indicator"` | one-locus `indicator` terms (a hand-entered genotype surface) |
| `"interaction"` | any term over two or more loci |

**The breeding value is `component_name = "additive"`** for effects
written in statistical coding at a single base allele frequency per
locus, which is what the generators
([`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md))
write. For hand-written functional terms
([`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md)
with `coding = "functional"`, or an `indicator` surface) `additive` is
the functional additive effect, not the breeding value: under functional
coding the average effect is \\\alpha = a + d(q - p)\\, and under
epistasis it depends on other loci and on LD. A true index on
`"additive"` therefore **warns** when an index trait has `indicator`
terms or hand-written interaction terms: a 0/1/2 dosage surface entered
with
[`genotype_terms()`](https://austin-putz.github.io/tidybreed/reference/genotype_terms.md)
is genetically additive but has no `additive` rows, and its additive
index is 0. `"total"` is a different objective (selection on genetic
value), not a breeding-value projection.

`ind_tgv` stores the **raw sum of the stored terms — no mean is added.**
A pure Cockerham model yields centered deviations; a raw `indicator`
surface yields absolute genotypic values with a non-zero mean by
construction, and is not silently re-centered. See
[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md),
which reports the implied genetic mean and writes it nowhere.

Each allele copy takes the **most specific** variant whose origin
predicate matches its `(line_origin, parent_origin)` label, falling back
per copy to the common variant. This per-copy fallback is what makes
crossbred genetic values correct — e.g. a Duroc x Landrace F1 is
centered against each parent line's own effects and base allele
frequency. **Imprinting** is a property of the effect, not of the trait:
a term scoped to one `parent_origin` reads only that parent's allele
copies.

Re-evaluating replaces an individual's rows for a trait: a component the
model no longer has is deleted, so no stale component survives, and a
component that is still there is updated in place, keeping any custom
columns added to it.

Optionally computes true selection index values by multiplying per-trait
genetic values (`component_name`, default the breeding value) by weights
from named indices defined with
[`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md),
and writes them to `ind_true_index`.

Every individual in the subset receives a value for every requested
trait — unlike
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
no sex-expression rule is applied here (`expressed_sex` is a property of
`phenotype_meta`, not of a genetic trait).

## Usage

``` r
add_tgv(
  tbl,
  trait_name = NULL,
  index_names = NULL,
  weight_type = c("index", "economic", "both"),
  component_name = "additive",
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
  in `trait_meta` are used (in `id_trait` order).

- index_names:

  Character vector of named index(es) from `index_meta` for which true
  index values are computed and written to `ind_true_index`. When `NULL`
  (default), no true index is computed. Every index trait must have an
  `ind_tgv` row for every individual of the subset (computed by this
  call or an earlier one).

- weight_type:

  Which weight column from `index_meta` to use: `"index"` uses
  `index_weight`, `"economic"` uses `economic_weight`, `"both"` computes
  and stores both (distinguished by `ind_true_index.weight_type`).
  Defaults to `"index"`.

- component_name:

  Which genetic value the true index weights: one of `"additive"`
  (default — selection indices are on breeding values), `"dominance"`,
  `"indicator"`, `"interaction"`, or `"total"` (the `ind_tgv_total`
  view). A component the trait's model has no terms for contributes 0.
  Stored in `ind_true_index.component_name`.

- overwrite_index:

  Logical. When `FALSE` (default), individuals that already have a true
  index value for the given `(index_name, weight_type, component_name)`
  are skipped. When `TRUE`, existing rows are deleted and recomputed
  (use when index weights have changed).

- ...:

  Optional extra columns written to `ind_tgv` (scalars only; broadcast
  to every row written).

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

[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md)
and
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
for writing the terms this evaluates,
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
[`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md),
[`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Every component of the model, for one generation
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(gen == 2L) |>
  add_tgv("ADG")

# The breeding values, and the derived total
pop |> get_table("ind_tgv") |>
  dplyr::filter(component_name == "additive") |> dplyr::collect()
pop |> get_table("ind_tgv_total") |> dplyr::collect()

# Known coefficients (a QTL map, GWAS estimates), line-specific, through
# the writer rather than the generator; add_tgv() evaluates every owner.
# center_value is filled from the Duroc founders' allele frequencies.
pop <- define_genome_effect_terms(pop, "ADG",
  data.frame(term_id = 1:2, locus_name = c("Locus_10", "Locus_44"),
             contrast_name = "additive", genome_value = c(0.4, -0.2)),
  origin = list(line_match_type = "exact", line_name = "Duroc"),
  base_tbl = get_table(pop, "founder_haplotypes") |>
    dplyr::filter(line_name == "Duroc"),
  effect_owner = "qtl_map")

# Genetic values + true index values (index and economic weights) on the
# breeding values, written to ind_true_index
pop <- pop |>
  get_table("ind_meta") |>
  dplyr::filter(gen == 2L) |>
  add_tgv(c("ADG", "BW"), index_names = "terminal", weight_type = "both")
} # }
```

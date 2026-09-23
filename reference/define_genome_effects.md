# Define genome effects as terms over one or more loci

Writes genome effects in the term / member / origin form: one
coefficient (`genome_effects`) over one or more loci
(`genome_effect_members`), each locus optionally scoped to allele copies
of a given line and/or parent of origin
(`genome_effect_member_origins`).

`terms` is a **long data frame, one row per (term x locus)** — the same
shape as every table
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
returns. Scope is supplied separately in `origin` so the common case
stays flat.

`define_genome_effects()` writes any effect you supply;
`define_*_effects()` functions such as
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
sample effects of one shape and write them through the same path.

## Usage

``` r
define_genome_effects(
  pop,
  trait_name,
  terms,
  effect_owner = "custom",
  mode = c("append", "replace_scope", "replace_owner", "replace_trait"),
  origin = NULL,
  base_tbl = NULL,
  require_complete = FALSE,
  allow_reserved_owner = FALSE
)
```

## Arguments

- pop:

  A `tidybreed_pop`.

- trait_name:

  Character scalar. Must exist in `trait_meta`.

- terms:

  Long data frame; see **The `terms` data frame**.

- effect_owner:

  Character scalar naming the writer that owns these rows, **for
  replacement only**. Owners always sum and are never selected between,
  because `effect_owner` is part of the fallback-family signature.
  Default `"custom"`.

- mode:

  One of `"append"`, `"replace_scope"`, `"replace_owner"`,
  `"replace_trait"`.

- origin:

  Scope; see **Scope (`origin`)**.

- base_tbl:

  Optional `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (optionally filtered) selecting the allele copies whose frequencies
  fill `center_value` on any `additive` or `dominance` member that has
  none: `founder_haplotypes`, `ind_haplotype`, or any table with an
  `id_ind` column, from the same `pop` — see
  [`extract_allele_freq()`](https://austin-putz.github.io/tidybreed/reference/extract_allele_freq.md).
  An explicit `center_value` is never overwritten; `indicator` members
  are never touched; the base is queried only if some centre is actually
  missing. `NULL` (default) fills nothing, and a missing centre is an
  error. The fill is Cockerham `p` only — functional coding writes `0.5`
  explicitly
  ([`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md)
  does). One `base_tbl` gives one `p` per locus per call, so line-scoped
  surfaces are written one line at a time with `mode = "append"`.

- require_complete:

  Logical. When `TRUE`, an indicator surface must name every reachable
  `(copy_count, dosage)` state on every member — including
  `copy_count_value = 0` where a chromosome can be absent. Default
  `FALSE` (sparse: a cell you do not write contributes zero).

- allow_reserved_owner:

  Logical. Permit writing under a package-reserved `effect_owner`.
  Default `FALSE`; the package's own generator
  ([`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md))
  writes under its reserved owner through the same engine.

## Value

The `tidybreed_pop`, invisibly.

## The `terms` data frame

|  |  |  |
|----|----|----|
| Column | Required | Meaning |
| `term_id` | yes, unless one term | Groups rows into one term. **User-facing only** — never stored; the writer replaces it with an `id_genome_effect`. Any atomic type |
| `genome_value` | yes | The term's coefficient. Constant within a `term_id` |
| `effect_name` | no | Per-term label. Constant within a `term_id`; no mathematical meaning |
| `locus_name` | yes | Resolved to `locus_id`; each locus at most once per term |
| `contrast_name` | yes | `"additive"`, `"dominance"` or `"indicator"` |
| `center_value` | non-indicator | `p` for Cockerham coding, `0.5` for functional |
| `copy_count_value`, `dosage_value` | indicator | The local genotype state. `copy_count_value` is inferred at ordinary diploid-autosomal loci |

A single-term call may omit `term_id` entirely.

## Scope (`origin`)

- `NULL` (default) — the common scope: no origin rows, matches every
  allele copy.

- A **named scalar list**, e.g.
  `list(line_name = "Duroc", parent_origin = 1)` — one scope applied to
  every member of every term. Accepted names are `line_match_type`,
  `line_name`, `parent_origin` and `copy_count`.

- A **data frame** with columns `term_id`, `locus_name`,
  `line_match_type`, `line_name`, `parent_origin`, `copy_count` —
  per-member scopes, needed for the exact multisets a `dominance` or
  `indicator` member takes. Keyed by `locus_name`, so you never touch
  canonical slot order.

An `additive` member takes at most one origin row; a genotype member
takes an exact multiset whose `copy_count`s sum to the state's copy
count (2 for `dominance`, `copy_count_value` for `indicator`). A scalar
list therefore scopes `additive` members but is usually not enough for a
genotype member — the validator says so, naming the locus.

## Replacement modes

- `"append"` — insert; a term duplicating an existing family + scope
  identity is rejected.

- `"replace_scope"` — delete only the variants in
  `(trait_name, effect_owner)` whose origin predicate **equals** the
  supplied scope, then insert. This is what
  [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
  uses, so successive common / line-A / line-B calls each replace only
  their own variant. Requires `origin` to be `NULL` or a scalar list (a
  per-member `origin` data frame has no single scope to key on).

- `"replace_owner"` — replace everything under
  `(trait_name, effect_owner)`.

- `"replace_trait"` — clear every owner's terms for the trait. Never
  implied.

## See also

[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md)
and
[`genotype_terms()`](https://austin-putz.github.io/tidybreed/reference/genotype_terms.md)
build `terms`;
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
for generated additive QTL effects.

## Examples

``` r
if (FALSE) { # \dontrun{
# One dominance term, Cockerham coding at p = 0.3
pop <- pop |> define_genome_effects(
  trait_name = "ADG",
  terms = data.frame(locus_name    = "Locus_10",
                     contrast_name = "dominance",
                     center_value  = 0.3,
                     genome_value  = 0.8)
)

# A hand-entered 3x3 A x A surface: nine cells, nine terms, two members each
cells <- expand.grid(g1 = 0:2, g2 = 0:2)
cells$value <- c(0, 0, 0,  0, 1.4, 2.1,  0, 2.1, 3.6)
surface <- rbind(
  data.frame(term_id = seq_len(nrow(cells)), locus_name = "Locus_10",
             contrast_name = "indicator", dosage_value = cells$g1,
             genome_value  = cells$value),
  data.frame(term_id = seq_len(nrow(cells)), locus_name = "Locus_44",
             contrast_name = "indicator", dosage_value = cells$g2,
             genome_value  = cells$value)
)
pop <- pop |> define_genome_effects("ADG", surface[surface$genome_value != 0, ],
                                    effect_owner = "epistasis_AxA")

# Reciprocal dominance: the F1 value depends on which parent gave which line
pop <- pop |> define_genome_effects(
  "ADG",
  terms = data.frame(term_id = 1L, locus_name = "Locus_10",
                     contrast_name = "dominance",
                     center_value = 0.3, genome_value = 1.2),
  origin = data.frame(term_id = 1L, locus_name = "Locus_10",
                      line_match_type = "exact",
                      line_name     = c("Duroc", "Landrace"),
                      parent_origin = c(1L, 2L),
                      copy_count    = c(1L, 1L)),
  effect_owner = "reciprocal"
)

# Let the writer fill Cockerham p from a base population: leave
# center_value out and pass base_tbl (see extract_allele_freq()).
pop <- pop |> define_genome_effects(
  "ADG",
  data.frame(locus_name = "Locus_10", contrast_name = "dominance",
             genome_value = 0.8),
  base_tbl = get_table(pop, "founder_haplotypes") |>
    dplyr::filter(line_name == "A")
)
} # }
```

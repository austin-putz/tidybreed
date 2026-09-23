# Define additive QTL effects for one or more traits

Selects QTL from a filtered `genome_meta` table and writes one order-one
`additive` term per locus through the same engine as
[`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md),
under the reserved effect owner `"generated_additive_tbv"`.
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
reads order-one `additive` variants from that owner and nothing else, so
effects written here and effects a user writes with
[`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md)
can never be confused for one another.

**Single trait** (`trait_name` length 1) — two modes:

- **Manual**: pass `effects`, a numeric vector of length `n_qtl` (number
  of filtered loci) in ascending `locus_id` order.

- **Sampled**: draw effects from `distribution` (`"normal"` or
  `"gamma"`). If `scale_to_target = TRUE`, effects are rescaled using
  the Falconer formula so the expected additive variance in the base
  population equals the `target_add_var` stored for this trait.

**Multiple traits** (`trait_name` length \>= 2) — effects are drawn
jointly from a multivariate normal distribution keyed by the
additive-genetic covariance matrix `G`. Two locus-selection methods:

- `method = "shared"` — the loci in `tbl` become the QTL set of every
  trait, and each locus receives one joint draw.

- `method = "union"` — the loci in `tbl` form the candidate pool;
  per-trait membership is read from the terms already stored at this
  scope, and a locus draws jointly only for the traits it is a QTL for.

[`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md)
writes any effect you supply; `define_*_effects()` functions such as
this one sample effects of one shape and write them through the same
path.

## Usage

``` r
define_additive_effects(
  tbl,
  trait_name,
  effects = NULL,
  distribution = c("normal", "gamma"),
  G = NULL,
  method = c("shared", "union"),
  base_tbl = NULL,
  line_name = NULL,
  parent_origin = NULL,
  scale_to_target = TRUE,
  seed = NULL
)
```

## Arguments

- tbl:

  A `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)`("genome_meta")`
  (with an optional
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)).
  The filtered rows determine which loci are QTL.

- trait_name:

  Character scalar **or** vector. Name(s) of existing traits in
  `trait_meta`. When length \>= 2, effects are drawn jointly from
  `MVN(0, G)` and `G` / `method` become active.

- effects:

  Optional numeric vector of length `n_qtl` (manual mode, single trait
  only), in ascending `locus_id` order. Error if
  `length(trait_name) > 1`.

- distribution:

  Character. `"normal"` (default) or `"gamma"`, used when `effects` is
  `NULL` and `length(trait_name) == 1`. Ignored for multi-trait.

- G:

  Optional numeric matrix of additive-genetic (co)variances (multi-trait
  only). Must be square and symmetric with side length
  `length(trait_name)`. When supplied, stored to `trait_var_comp` under
  `"gen_add"`. When `NULL`, read from `trait_var_comp`.

- method:

  Character. `"shared"` (default) or `"union"`. Multi-trait only.
  `"shared"` — all listed traits use the filtered loci as their shared
  QTL set. `"union"` — per-trait QTL sets are read from existing
  `genome_effects` rows, restricted to the filtered loci.

- base_tbl:

  Optional `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (optionally filtered) selecting the allele copies that define base
  allele frequencies: `founder_haplotypes`, `ind_haplotype`, or any
  table with an `id_ind` column. Must come from the same `pop` as `tbl`.
  `NULL` (default) resolves to the founder pool of the line the effect
  applies to — see *Which population centers the effects* above.

- line_name:

  Optional character. When set, effects are scoped to allele copies of
  this genetic line: a copy whose `line_origin` matches takes these
  values, and falls back per copy to the common variant where no
  line-specific one exists. Also selects the default `base_tbl`. `NULL`
  (default) means the common scope, matching every copy.

- parent_origin:

  Optional `1` (sire / parent_1) or `2` (dam / parent_2) — imprinting,
  restricting the term to copies inherited from that parent. `NULL`
  (default) means both parents' copies. **Per trait**: a scalar is
  recycled, a vector must match `trait_name` positionally, or name its
  entries by trait. A call mixing origins across traits while supplying
  `G` is rejected — under random mating the paternal and maternal copies
  at a locus are independent, so the requested genetic covariance
  between a paternal-only and a maternal-only trait is zero and cannot
  be realized. For imprinting that varies locus by locus, write the
  terms with
  [`define_genome_effects()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effects.md).

- scale_to_target:

  Logical. If `TRUE`, rescale effects so the expected additive variance
  equals the stored `target_add_var`:
  `V_A = sum_j n_eligible,j * p_j q_j a_j^2`, where `n_eligible` is 2
  for an unparented term and 1 for a parent-qualified one.

- seed:

  Optional integer for reproducibility.

## Value

The modified `tidybreed_pop` (invisibly).

## Which population centers the effects

Base allele frequencies center the true breeding value (the Falconer
`allele - p` term) and set the `2pq` denominator used by
`scale_to_target`. They come from `base_tbl`, a filtered
`tidybreed_table` whose identity says *what kind of thing* is selected
and whose
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
says *which* (see
[`extract_allele_freq()`](https://austin-putz.github.io/tidybreed/reference/extract_allele_freq.md)
for the three accepted shapes: the founder pool, allele copies in
`ind_haplotype`, or individuals from any table with `id_ind`).

When `base_tbl = NULL` the base is **the population the effect applies
to**, resolved with the same `line -> NULL` precedence as
`resolve_genome_map()`: a line-specific effect (`line_name = "A"`)
centers on line A's own founder pool, or on the shared
(`line_name = NULL`) pool when no named pool exists; a population-wide
effect (`line_name = NULL`) centers on the whole founder table. Only
that last case warns when the founder table holds more than one pool:
pooling divergent lines overstates within-line heterozygosity — the
Wahlund effect. Two lines fixed for opposite alleles each have zero
within-line variance, but pool to `p = 0.5` and an apparent `2pq = 0.5`;
the inflated denominator then makes `scale_to_target` **under**-scale
the effects, and realized within-line additive variance falls short of
`target_add_var`. An explicit `base_tbl` is an intentional selection and
never warns — pass `base_tbl = get_table(pop, "founder_haplotypes")` to
pool on purpose, which is how the common fallback variant of a
crossbreeding model is defined.

A selected QTL locus with no allele copies in the base is an error,
never silently centered at `p = 0`.

The centering constant is stored per member as
`genome_effect_members.center_value` and travels with its
`genome_value`, so evaluation applies each allele copy's own line's
centering — a crossbred animal's line-A alleles are centered on line A
and its line-B alleles on line B.

## Scope, and what a re-run replaces

`line_name` and `parent_origin` compose into the single origin row an
`additive` member is allowed:

|             |                 |                                               |
|-------------|-----------------|-----------------------------------------------|
| `line_name` | `parent_origin` | Stored scope                                  |
| `NULL`      | `NULL`          | no origin rows (the common scope)             |
| `"A"`       | `NULL`          | `('exact', 'A', parent NULL, copy_count = 1)` |
| `NULL`      | `1` / `2`       | `('any', NULL, parent, copy_count = 1)`       |
| `"A"`       | `1` / `2`       | `('exact', 'A', parent, copy_count = 1)`      |

Re-running replaces **only the variant at the same scope**
(`mode = "replace_scope"`), so successive common / line-A / line-B calls
each keep the others: the per-copy fallback that makes crossbred
breeding values correct depends on all of them standing. It also means
changing `parent_origin` on a re-run **adds** a variant rather than
replacing one — the two are in a containment relation and both apply, to
different copies. That is legal and rarely intended, so the function
warns on exactly that case.

## See also

[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md),
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Single trait — all loci on chr 1-5 become QTL; scale to target variance
pop <- pop |>
  define_trait("ADG", target_add_var = 0.25) |>
  get_table("genome_meta") |>
  dplyr::filter(chr %in% 1:5) |>
  define_additive_effects("ADG", distribution = "normal")

# Multiple correlated traits — shared QTL set, joint MVN draw
G <- matrix(c(0.25, 0.10, 0.10, 0.30), 2, 2,
            dimnames = list(c("ADG", "BW"), c("ADG", "BW")))
pop <- pop |>
  define_effect_cov_matrix("gen_add", G) |>
  get_table("genome_meta") |>
  dplyr::filter(chr %in% 1:5) |>
  define_additive_effects(c("ADG", "BW"), G = G)

# Generation-0 individuals define the base allele frequencies
pop <- pop |>
  get_table("genome_meta") |>
  dplyr::filter(chr %in% 1:5) |>
  define_additive_effects("ADG",
    base_tbl = get_table(pop, "ind_meta") |> dplyr::filter(gen == 0L))

# Crossbreeding: three variants. The common fallback names its base to say
# "yes, pool"; each line's variant centers on its own founder pool by default.
gm <- pop |> get_table("genome_meta") |> dplyr::filter(chr %in% 1:5)
pop <- gm |> define_additive_effects("ADG",
               base_tbl = get_table(pop, "founder_haplotypes"))
pop <- gm |> define_additive_effects("ADG", line_name = "Duroc")
pop <- gm |> define_additive_effects("ADG", line_name = "Landrace")

# Duroc allele copies wherever they sit, including inside crossbreds
pop <- gm |> define_additive_effects("ADG", line_name = "Duroc",
  base_tbl = get_table(pop, "ind_haplotype") |>
    dplyr::filter(line_origin == "Duroc"))
} # }
```

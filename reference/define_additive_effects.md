# Define additive QTL effects for one or more traits

Selects QTL from a filtered `genome_meta` table and writes additive
effects to the `genome_effects` table.

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

- `method = "shared"` — the loci in `tbl` become the shared QTL set for
  all traits. Loci that are QTL for only a subset of traits in
  `genome_effects` also receive independent draws (with the diagonal
  variance of `G` for that trait).

- `method = "union"` — the loci in `tbl` form the candidate pool;
  per-trait membership is determined from existing rows in
  `genome_effects`.

The `base` argument controls which allele frequencies are used:

- `"founder_haplotypes"` (default) — computes allele frequencies
  directly from the `founder_haplotypes` table (requires
  [`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)
  was called). Restrict to one founder pool with `base_line_name`. (This
  does **not** read `genome_meta.founder_allele_freq`, which is
  informational only.)

- `"current_pop"` — computes allele frequencies from the current
  `ind_haplotype` table. Pass a filtered `tidybreed_table` via
  `base_tbl` to restrict which individuals define the base population.

## Usage

``` r
define_additive_effects(
  tbl,
  trait_name,
  effects = NULL,
  distribution = c("normal", "gamma"),
  G = NULL,
  method = c("shared", "union"),
  base = c("founder_haplotypes", "current_pop"),
  base_tbl = NULL,
  base_line_name = NULL,
  line_name = NULL,
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

- base:

  Character. `"founder_haplotypes"` (default) or `"current_pop"`.

- base_tbl:

  Optional `tidybreed_table` (from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  on any table with an `id_ind` column) used when `base = "current_pop"`
  to restrict which individuals define the base allele frequencies. When
  `NULL`, all individuals in `ind_haplotype` are used. Ignored (with a
  warning) when `base = "founder_haplotypes"` — use `base_line_name`
  there.

- base_line_name:

  Optional character, `base = "founder_haplotypes"` only. Which founder
  pool defines the base allele frequencies. **Defaults to `line_name`**,
  so line-specific effects are centered on their own line; pass `NULL`
  explicitly to pool every line instead. Errors if no
  `founder_haplotypes` rows carry that line. See *Which population
  centers the effects* above.

- line_name:

  Optional character. When set, effects are tagged to this genetic line:
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
  then prefers these rows for alleles whose `line_origin` matches,
  falling back per-locus to the population-wide rows. Also becomes the
  default for `base_line_name`. `NULL` (default) means population-wide
  effects.

- scale_to_target:

  Logical. If `TRUE`, rescale effects using the Falconer formula so the
  expected additive variance equals the stored `target_add_var`.

- seed:

  Optional integer for reproducibility.

## Value

The modified `tidybreed_pop` (invisibly).

## Which population centers the effects

Base allele frequencies center the true breeding value (the Falconer
`allele - p` term) and set the `2pq` denominator used by
`scale_to_target`. By default they come from **the population the effect
applies to**: `base_line_name` inherits `line_name`, so a line-specific
effect is centered on that line's own founder pool and a population-wide
effect (`line_name = NULL`) on the whole founder base.

This matters because pooling divergent lines overstates within-line
heterozygosity — the Wahlund effect. Two lines fixed for opposite
alleles each have zero within-line variance, but pool to `p = 0.5` and
an apparent `2pq = 0.5`; the inflated denominator then makes
`scale_to_target` **under**-scale the effects, and realized within-line
additive variance falls short of `target_add_var`. Pass
`base_line_name = NULL` explicitly to force pooling anyway.

The centering constant is stored per row in
`genome_effects.base_allele_freq` and travels with its `genome_value`,
so
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
applies each allele's own line's centering — a crossbred animal's line-A
alleles are centered on line A and its line-B alleles on line B.

Calling this function again for the same
`(trait_name, genome_effect_type, line_name)` replaces the existing rows
in `genome_effects`.

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

# current_pop: use generation-0 individuals to define base allele frequencies
gen0_tbl <- get_table(pop, "ind_meta") |> dplyr::filter(gen == 0L)
pop <- pop |>
  get_table("genome_meta") |>
  dplyr::filter(chr %in% 1:5) |>
  define_additive_effects("ADG", base = "current_pop", base_tbl = gen0_tbl)

# Crossbreeding: each line's effects centered on its own founder pool.
# base_line_name inherits line_name, so nothing extra is needed.
pop <- pop |>
  get_table("genome_meta") |> dplyr::filter(chr %in% 1:5) |>
  define_additive_effects("ADG", line_name = "Duroc")
pop <- pop |>
  get_table("genome_meta") |> dplyr::filter(chr %in% 1:5) |>
  define_additive_effects("ADG", line_name = "Landrace")

# Line-specific effects, but deliberately centered on the pooled base
pop <- pop |>
  get_table("genome_meta") |> dplyr::filter(chr %in% 1:5) |>
  define_additive_effects("ADG", line_name = "Duroc", base_line_name = NULL)
} # }
```

# Define additive, dominance and epistatic QTL effects calibrated to targets

Selects QTL from a filtered `genome_meta` table, samples additive,
dominance and additive-by-additive (A x A) effects, and **calibrates**
them so that the stored targets `G_A`, `G_D` and `G_AA` are delivered
exactly under a named anchor. Writes them under the reserved effect
owner `"generated"`, replacing the traits' whole generated model. The
non-additive sibling of
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md):
the same pipe subject, target rules, base resolution and owner.

Use
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
for line-scoped or parent-of-origin additive models (crossbreeding
variants); use this function when dominance or epistasis may be part of
the model. With an additive-only target the two write identical rows for
the same seed.

## Usage

``` r
define_genome_effects(
  tbl,
  trait_name,
  G_A = NULL,
  G_D = NULL,
  G_AA = NULL,
  trait_var_comp_tbl = NULL,
  pairs = NULL,
  n_pairs = NULL,
  anchor = c("genic", "realised"),
  dominance_degree_mean = 0.19,
  dominance_degree_sd = 0.097,
  inbreeding_depression = NULL,
  base_tbl = NULL,
  warn_bounds = c(0.8, 1.25)
)
```

## Arguments

- tbl:

  A `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)`("genome_meta")`,
  optionally filtered. Its loci are the QTL of every trait.

- trait_name:

  Character vector of existing traits.

- G_A, G_D, G_AA:

  Optional `k x k` targets (a number for one trait), named in
  `trait_name` order or unnamed. Written with the effects; never over a
  stored block. `NULL` reads the block from `trait_var_comp_tbl`.

- trait_var_comp_tbl:

  Optional filtered `get_table(pop, "trait_var_comp")`: the stored rows
  to use for blocks not passed. `NULL` reads the stored population-wide
  rows.

- pairs:

  Optional data frame `locus_name_1`, `locus_name_2`.

- n_pairs:

  Optional number of random pairs, `1..floor(m / 2)`.

- anchor:

  `"genic"` (default) or `"realised"`.

- dominance_degree_mean, dominance_degree_sd:

  Mean and standard deviation of the dominance degrees that shape the
  dominance architecture. Defaults `0.19` and `0.097`.

- inbreeding_depression:

  Optional numeric vector named by trait.

- base_tbl:

  Optional `tidybreed_table` selecting the base population, as in
  [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md).
  Required under `"realised"`.

- warn_bounds:

  `c(lower, upper)` for the comparison with another population, or
  `NULL` to turn it off.

## Value

The `tidybreed_pop`, invisibly.

## Details

**The model.** In functional coding the genotypic value is
`g = sum_j a_j (x_j - 1) + sum_j d_j 1[x_j = 1] + sum_(k,l) e_kl (x_k - 1)(x_l - 1)`,
`x` the dosage of allele 1. Its exact statistical (NOIA) re-expression
has additive effects `alpha_j = a_j + b_j d_j + sum_l e_jl c_l`
(`c_l = 2 p_l - 1`), dominance effects `d` and A x A effects `e`. The
three targets are covariances of those three components.

**When the result is exact.** Every present block (`G_A`, and `G_D` /
`G_AA` when present) is delivered under the named anchor to a relative
tolerance of `1e-8` on the correlation scale. That is checked, never
assumed: a calibration that misses it is an error before anything is
written.

**How.** Highest order first: the A x A architecture is calibrated to
`G_AA`; the dominance architecture, `d = h |a|` with dominance degrees
`h ~ N(dominance_degree_mean, dominance_degree_sd)`, is calibrated to
`G_D`; then, with the coupling `b d + sum e c` the first two induce now
fixed, the additive architecture is solved so the statistical additive
effects deliver `G_A`. Each stage is a `k x k` right factor of its drawn
architecture, so the same seed gives the same architecture.

**The additive floor.** Dominance and epistasis already induce additive
variance (`b d` and `e c` are part of `alpha`). Part of it the additive
architecture can cancel, part it cannot; `G_A` below what it cannot
cancel is refused, and the error names that minimum per trait. The floor
belongs to the **sampled architecture**, not to `G_D` / `G_AA` alone:
another seed gives another floor, and with every base frequency at `0.5`
under `"genic"` it is zero however large `G_D` is.

**The anchor** is the reference covariance the calibration is exact for:

- `"genic"` (default): the random-mating (HWE + linkage-equilibrium)
  limit at the base allele frequencies. Use it for multi-generation
  studies.

- `"realised"`: the covariances of the individuals `base_tbl` selects,
  linkage disequilibrium included. The effects are still **stored** with
  HWE-referenced contrasts centred at the base frequencies, so the
  stored `additive` / `dominance` split is not the cohort's; the total
  is exact, and
  [`extract_genetic_variance()`](https://austin-putz.github.io/tidybreed/reference/extract_genetic_variance.md)`(anchor = "realised")`
  on the same individuals gives the targets back. `base_tbl` must select
  individuals with complete genotypes, and the in-memory design arrays
  have a size limit.

Under linkage disequilibrium the `additive` component of a realised
cohort is a contrast component, not a joint regression: calibration is
exact for the anchor and nothing more.

After calibration the delivered blocks are compared with what another
population sees: the genic limit under `"realised"`; the base
individuals' own covariances under `"genic"` when `base_tbl` selects
individuals; the founder pool's expectation of the additive block only
(dominance and A x A would need the pool's multi-locus LD). A departure
outside `warn_bounds` warns (the pool comparison is a message). Nothing
is stored.

## Targets

Each block (`additive`, `dominance`, `additive_by_additive`) comes from
exactly one source:

- a passed matrix (`G_A`, `G_D`, `G_AA`; a number for one trait),
  written to `trait_var_comp` in the same transaction as the effects. It
  is refused when that block is already stored for any of the traits,
  even as an identical matrix, and whether or not `trait_var_comp_tbl`
  shows it;

- otherwise the rows of `trait_var_comp_tbl` (a filtered
  `get_table(pop, "trait_var_comp")`), or, when it is `NULL`, the stored
  population-wide rows.

A block that comes from neither is **absent** from the model; a zero
matrix is an explicit exact zero, and writes no terms. An `additive`
block is required. A model with no non-zero dominance or A x A block
takes the additive-only route: the same draw, calibration and rows as
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md).
To add a block later, store or pass it and re-run: the additive
architecture draw is unchanged by the new block (same seed).

Non-zero dominance or A x A targets need a positive-definite `G_A`: this
release solves the additive stage for a full-rank factor. A correlation
of +/-1 or a zero additive variance is accepted for additive-only
models. Singular targets are reported with a message, since a typed `1`
meant as `0.99` looks the same.

## Pairs

A x A acts on pairs of QTL, used only with an `additive_by_additive`
block:

- `pairs`, a data frame `locus_name_1`, `locus_name_2`: your design. A
  locus may appear in several pairs (a hub gene); unknown loci, loci
  outside the filter, self-pairs and repeated pairs (in either order)
  are errors.

- `pairs = NULL`: a random matching, each QTL in at most one pair, as
  AlphaSimR does. `n_pairs = NULL` pairs every QTL once (`floor(m / 2)`
  pairs), and a message says so.

Within a pair the loci are ordered by name (C locale), and pairs are
written sorted by locus id, so the stored order does not depend on the
draw.

## Crossbreeding and scope

Effects are common to every line (no `line_name` or `parent_origin`;
line or parent-specific non-additive effects are a later release).
Common effects still give heterosis, through allele-frequency
differences between lines. The default base pools the founder table, so
the targets then hold at the **pooled** frequencies (with a Wahlund
warning); pass `base_tbl` to calibrate on one reference line.

## What a re-run replaces

Every `"generated"` term of the traits, at every scope: line-scoped
additive variants from
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
included (the message counts them). Terms of other owners are untouched.
Afterwards
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
refuses these traits while their generated model has dominance or A x A
terms; re-run this function with an additive-only target, or use
[`remove_generated_effects()`](https://austin-putz.github.io/tidybreed/reference/remove_generated_effects.md).

## Generated means calibrated

There is no option to write fixed or unscaled effects. Exact
coefficients go through
[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md)
with
[`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md)
/
[`aa_terms()`](https://austin-putz.github.io/tidybreed/reference/aa_terms.md)
under a user owner, for example:

    pop |>
      define_genome_effect_terms(
        trait_name = "ADG",
        terms = rbind(
          ad_terms(locus_name = c("Locus_10", "Locus_44"),
                   a = c(0.30, -0.12), d = c(0.10, 0.05), p = c(0.35, 0.60),
                   coding = "functional"),
          aa_terms(locus_name_1 = "Locus_10", locus_name_2 = "Locus_44",
                   e = 0.08, p_1 = 0.35, p_2 = 0.60, coding = "functional")))

## Inbreeding depression

`inbreeding_depression` (the drop in the mean per unit of inbreeding
`F`, positive = depression, `sum 2pq d`) sets each named trait's mean
dominance degree. It is exact for one trait. With several traits each
requested trait's mean is solved as if it stood alone, and the joint
calibration then mixes the columns, so the delivered depression is
approximate, with no closeness guarantee; a message gives requested and
delivered values. It is not stored.

## See also

[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md),
[`extract_genetic_variance()`](https://austin-putz.github.io/tidybreed/reference/extract_genetic_variance.md),
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md),
[`remove_generated_effects()`](https://austin-putz.github.io/tidybreed/reference/remove_generated_effects.md),
[`aa_terms()`](https://austin-putz.github.io/tidybreed/reference/aa_terms.md).

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- pop |> define_trait("ADG")
set.seed(1)
pop <- pop |>
  get_table("genome_meta") |>
  dplyr::filter(chr_name %in% c("1", "2")) |>
  define_genome_effects("ADG", G_A = 0.3, G_D = 0.1, G_AA = 0.05)

# Additive-only first, dominance added later: the additive draw is the same
set.seed(2)
pop <- pop |> get_table("genome_meta") |> define_genome_effects("BF", G_A = 1)
pop <- pop |> define_effect_cov_matrix("dominance", 0.2, trait_name = "BF")
set.seed(2)
pop <- pop |> get_table("genome_meta") |> define_genome_effects("BF")
} # }
```

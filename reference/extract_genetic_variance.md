# Measure the genetic covariance blocks of a cohort

Reports the genetic (co)variances a selected group of individuals
actually carries, in the shape of `trait_var_comp`: one row per
`(effect_name, trait_name_1, trait_name_2)`, every trait pair in both
orientations. So "did the model deliver its target?" is one
[`dplyr::inner_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)
against the stored targets. Read-only: nothing is written, and `ind_tgv`
is not touched.

`tbl` selects the individuals, exactly as for
[`add_tgv()`](https://austin-putz.github.io/tidybreed/reference/add_tgv.md):
any filtered table with `id_ind` (`ind_meta` by generation or date,
`ind_phenotype`, `ind_ebv`, ...).

## Usage

``` r
extract_genetic_variance(
  tbl,
  trait_name = NULL,
  base_tbl = NULL,
  anchor = c("realised", "genic")
)
```

## Arguments

- tbl:

  A `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md),
  optionally filtered, selecting the individuals.

- trait_name:

  Character vector of traits. `NULL` (default) is every trait with
  stored genome-effect terms, in `trait_meta` order. A named trait
  without terms is an error.

- base_tbl:

  Genic anchor only: a filtered `tidybreed_table` whose allele
  frequencies define the expectation (`founder_haplotypes`,
  `ind_haplotype` copies, or individuals). `NULL` (default) uses the
  selected individuals.

- anchor:

  `"realised"` (default) or `"genic"`.

## Value

A tibble with columns `effect_name`, `trait_name_1`, `trait_name_2`,
`cov_value`, `n_ind`, `decomposition`, `anchor`, sorted by effect
(`additive`, `dominance`, `additive_by_additive`, `unpartitioned`,
`between_components`, `total`) and trait order. A
[`message()`](https://rdrr.io/r/base/message.html) names the population
the estimates describe; the `anchor` column travels with the rows, but
which cohort or base produced them is only in that message.

## Two anchors, two definitions

- `anchor = "realised"` (default) measures the selected individuals:
  sample covariances (divisor `n - 1`) of each block's values. Under
  linkage disequilibrium (LD) or departures from Hardy-Weinberg
  equilibrium (HWE) the blocks correlate; those cross-block covariances
  are reported in `between_components`, so the block rows plus
  `between_components` sum exactly to `total`.

- `anchor = "genic"` is the HWE + linkage-equilibrium expectation at
  `base_tbl`'s allele frequencies. The blocks are orthogonal by
  construction, there is no `between_components` row, and `total` is the
  sum of the blocks. Available for fully decomposable models only (case
  1 below).

`base_tbl` is like the generators' `base_tbl` (see
[`extract_allele_freq()`](https://austin-putz.github.io/tidybreed/reference/extract_allele_freq.md))
with two differences: its `NULL` default is **the selected individuals**
(their whole-genotype frequencies, whatever table `tbl` is), not the
founder pool, so it follows the cohort as frequencies drift and does not
reproduce a generation target; and it is genic only (a non-`NULL`
`base_tbl` with `"realised"` is an error). An explicit `base_tbl` keeps
[`extract_allele_freq()`](https://austin-putz.github.io/tidybreed/reference/extract_allele_freq.md)'s
semantics: an `ind_haplotype` filter selects allele copies. To compare
with a generation target, pass the generation base.

## What the realised blocks are, under LD

The realised blocks are the sample covariances of the NOIA **contrast
components** at the cohort's frequencies (the source project's
`nonadd_decompose()`): each locus' heterozygosity is regressed on its
own dosage only, and the additive block is `(X - 2p) alpha` with
`alpha_j = a_j + b_j d_j + sum_l e_jl (2 p_l - 1)`. Under LE these are
the orthogonal least-squares components. Under LD they are not:
dominance at another locus, or a pair product, can still regress on the
dosages, so the `additive` row is **not** the variance of the cohort's
joint least-squares additive projection (the breeding values a
regression of the genetic value on all dosages would give), and
`between_components = 0` does not certify that it is. For example, one
pair `e (g_1 - 1)(g_2 - 1)` on a cohort with both loci at `p = 0.5` and
HWE margins but in LD reports `additive = 0`, while `lm(g ~ g_1 + g_2)`
explains part of its variance. `"full"` means every stored term has a
supported shape, not that the cohort's additive breeding-value variance
has been recovered.

## What the function can decompose

Each trait's stored model falls in one of three cases, reported in the
`decomposition` column.

1.  `"full"`: every term is common-scope (no line or parent-of-origin
    scope) on diploid-autosomal loci, and is an order-one `additive`,
    order-one `dominance`, order-one diploid `indicator` (any of the
    three genotype states), or two-member `additive x additive` term, of
    any owner and at any centre. These are converted to functional
    effects and re-projected on the NOIA (natural and orthogonal
    interactions) model at the cohort's frequencies, so both codings of
    one model give the same report.

2.  `"additive_only"`: every term is an order-one `additive` term and
    some have a line or parent-of-origin scope (crossbreeding,
    imprinting). The `additive` row is the **evaluated additive
    variance**, the covariance of the evaluated additive component, not
    a re-projection.

3.  `"partial"`: anything else. The covered terms are decomposed as in
    case 1; everything else goes to `unpartitioned`, the variance of
    those terms' value alone (not `residual`, which in this package is
    environmental noise). In this case the `additive` row is the
    additive contrast component of the covered terms, not the breeding
    value of the whole model.

Scope variants of one term compete (the most specific matching scope
wins), so a term family with any scoped variant is uncovered as a whole.
An off-diagonal row of two traits in different cases carries the less
complete label.

## Which rows appear

Rows follow the statistical decomposition, not the stored term kinds.
`additive` is reported for every trait with a decomposed term (a
one-locus genotype surface equal to the dosage is all additive;
dominance and interaction terms induce additive effects too).
`dominance` appears when the converted model has a non-zero dominance
coefficient, `additive_by_additive` when it has a non-zero pair. A
coefficient whose stored contributions cancel to within their
floating-point rounding (a genotype surface that is exactly linear,
written in any row order) counts as zero. A block whose variance is 0 in
this cohort is reported as 0; a missing row means the model has no such
block, so a stored target for it falls out of the `inner_join()` and
into the
[`dplyr::anti_join()`](https://dplyr.tidyverse.org/reference/filter-joins.html).
An off-diagonal block row appears only when both traits have the block.
`total` is always reported.

A locus fixed in the cohort contributes nothing to the dominance or
interaction blocks, but an interaction with a fixed partner is still a
real additive effect (`e (g_1 - 1)(g_2 - 1) = e (g_1 - 1)` when
`g_2 = 2`), and is kept in `additive`.

## Comparing with a target

The effect names are `trait_var_comp`'s (the interaction block is
`additive_by_additive`, not the `ind_tgv` component name `interaction`).
A mismatch with the stored target is information, never an error. The
comparison is like-for-like only when the trait's terms are all
`"generated"` **and** the anchor, the reference population, the measured
model and the calibrated scope match the generation call: a generated
line variant competes with common fallback terms, so a mixed cohort
measures the active combination; a restored or later-edited founder pool
is not the generation-time base. Stored diagonals are generation targets
under their calibration contract (genic or realised), not universally a
genic total.

`define_phenotype(prevalence = )` sums the stored diagonals. On a cohort
with LD, departures from HWE or drifted frequencies the realised `total`
differs from that sum by `between_components` and the drift in each
block, which this function shows. Even a matching variance does not
guarantee the prevalence: the threshold also assumes a near-normal
liability, which a skewed finite-locus model need not give.

## Precision

Each block is computed to an **absolute** precision of about machine
epsilon times the size of the model's coefficients, not relative to the
block itself. Re-projecting the stored terms adds the dominance and pair
coefficients into the additive one (`alpha = a + b d + sum e c`), so an
additive block many orders of magnitude below the dominance or A x A
block keeps few correct digits: a model stored with an additive variance
of exactly `1e-24` beside a dominance variance of 1 is reported as about
`1.0001e-24`. At ordinary ratios the error is at rounding level.

## See also

[`add_tgv()`](https://austin-putz.github.io/tidybreed/reference/add_tgv.md),
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md),
[`extract_allele_freq()`](https://austin-putz.github.io/tidybreed/reference/extract_allele_freq.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Realised blocks of generation 5
meas <- get_table(pop, "ind_meta") |>
  dplyr::filter(gen == 5L) |>
  extract_genetic_variance()

# Target vs measured, population-wide targets (a line: line_name == "A")
targets <- get_table(pop, "trait_var_comp") |>
  dplyr::filter(is.na(line_name)) |>
  dplyr::collect() |>
  dplyr::select(effect_name, trait_name_1, trait_name_2, target = cov_value)
dplyr::inner_join(targets, meas,
                  by = c("effect_name", "trait_name_1", "trait_name_2"))
dplyr::anti_join(targets, meas,      # targets with no measured block
                 by = c("effect_name", "trait_name_1", "trait_name_2"))

# Genic: a common-scope model generated with anchor = "genic" on the founder
# pool, measured at the founder pool's frequencies
get_table(pop, "ind_meta") |>
  extract_genetic_variance(anchor = "genic",
                           base_tbl = get_table(pop, "founder_haplotypes"))
} # }
```

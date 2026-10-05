# Define an observed phenotype

Registers an **observed phenotype**: the quantity written to
`ind_phenotype` when
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
is called. For simple (non-composite) traits, `phenotype_name` equals
`trait_name` and
[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md)
must already have been called for that name. For composite traits (e.g.
weaning weight assembled from direct and maternal components), no
[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md)
call is needed for the composite name itself.

Writes one row to `phenotype_meta`. If `residual_var` is supplied, also
writes an unconditional diagonal entry to `phenotype_var_comp`
(`effect_name = "residual"`). If `components` is supplied, writes one
row per component to `phenotype_components`.

## Usage

``` r
define_phenotype(
  pop,
  phenotype_name,
  type = c("continuous", "count", "categorical", "derived_formula"),
  mean = 0,
  expressed_sex = c("both", "M", "F"),
  repeatable = FALSE,
  min_value = NULL,
  max_value = NULL,
  prevalence = NULL,
  thresholds = NULL,
  cat_values = NULL,
  cat_names = NULL,
  store_liability = FALSE,
  residual_var = NULL,
  components = NULL,
  formula_tbv = NULL,
  formula = NULL,
  missing_component_action = c("skip", "error"),
  condition_change_action = c("error", "independent"),
  overwrite = FALSE
)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- phenotype_name:

  Character. Name of the observed phenotype; equals `trait_name` for
  simple (non-composite) traits.

- type:

  Character. One of `"continuous"`, `"count"`, `"categorical"`, or
  `"derived_formula"`. `"derived_formula"` phenotypes are computed at
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  time by evaluating the `formula` expression over already- recorded
  phenotype values for the same individuals; they have no TBV, no
  residual variance, and no QTL of their own.

- mean:

  Numeric. The intercept. Default `0`. A record is `mean` + the genetic
  value as stored + random effects + residual; nothing is added to make
  the realised mean hit `mean`. Generated effects have mean 0 in
  expectation at their base allele frequencies (Hardy-Weinberg and
  linkage equilibrium), so `mean` is the base phenotypic mean in
  expectation; a finite, selected or non-equilibrium base differs by its
  sample genetic mean, and hand-written functional terms by their
  implied mean (see
  [`ad_terms()`](https://austin-putz.github.io/tidybreed/reference/ad_terms.md)).
  For a particular realised base mean, set `mean` to the target minus
  the base individuals' mean of `ind_tgv_total`, measured after
  [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md).

- expressed_sex:

  Character. Who receives a phenotype record: `"both"` (default), `"M"`,
  or `"F"`.

- repeatable:

  Logical. Whether an individual can have multiple records for this
  phenotype (e.g. repeated test-day or litter-size records). Default
  `FALSE`.

- min_value, max_value:

  Numeric. Clipping bounds for count traits. `NULL` means no limit.

- prevalence:

  Numeric between 0 and 1. For categorical traits with one threshold
  (two categories), the fraction expected above the threshold. Mutually
  exclusive with `thresholds`. The liability carries the trait's total
  genetic value, so the threshold is placed from `mean`, the
  unconditional residual variance and the trait's stored genetic targets
  (`trait_var_comp`, population-wide rows): the sum of the `additive`,
  `dominance` and `additive_by_additive` diagonals, each counted only if
  the trait's model has terms of that kind.
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  errors when a kind of term the model has has no stored target, or when
  the model has terms outside those three kinds (an `indicator` surface,
  other interactions). The threshold uses the *target* at the reference
  population, so the realised prevalence of a selected or line-specific
  population differs. Every term of the trait must be owned by
  `"generated"`: a generator
  ([`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md))
  always calibrates its terms to the stored target, so the target
  describes them. Terms written with
  [`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md)
  carry values nothing checked against a target, so
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  errors for such a trait; give `thresholds`. Not valid for composite
  phenotypes (`components` or `formula_tbv`): their genetic liability
  combines several traits and contributors, which no stored variance
  describes. Give `thresholds` instead.

- thresholds:

  Numeric vector of length K−1 for K ordered categories. Liability
  cutpoints in ascending order. Mutually exclusive with `prevalence`.

- cat_values:

  Numeric vector of length K. Phenotype value stored in `ind_phenotype`
  for each category. Defaults to `1, 2, ..., K`.

- cat_names:

  Character vector of length K. Human-readable label per category, e.g.
  `c("Alive", "Dead")`. Must not contain commas.

- store_liability:

  Logical. When `TRUE`, the underlying liability value is written to the
  reserved `liability_value` column in `ind_phenotype`. Only meaningful
  for categorical traits.

- residual_var:

  Numeric or `NULL`. Scalar residual variance. When supplied, writes a 1
  × 1 unconditional residual block for this phenotype to
  `phenotype_var_comp` (`effect_name = "residual"`). It is an error when
  the phenotype already belongs to a multi-phenotype residual block —
  that block is redeclared as a whole with
  [`define_residual_cov()`](https://austin-putz.github.io/tidybreed/reference/define_residual_cov.md)
  — or when the phenotype's existing residual has realized draws in
  `ind_phenotype`. With `overwrite = TRUE` and no `residual_var`,
  `phenotype_var_comp` is left untouched. For heterogeneous or
  correlated residuals use
  [`define_residual_cov()`](https://austin-putz.github.io/tidybreed/reference/define_residual_cov.md).

- components:

  A data frame or `tibble` with one row per genetic component. Columns:

  - `source_trait_name` (required): component trait name in
    `trait_meta`.

  - `contributor_type` (required): `"self"`, `"dam"`, `"sire"`, or
    `"group"`.

  - `weight` (optional, default `1.0`): scalar multiplier.

  - `weight_type` (optional, default `"fixed"`): `"fixed"` or
    `"covariate"` (`weight * covariate`). Nothing else is implemented,
    and anything else is rejected here rather than at
    [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
    time.

  - `covariate_name` (optional): covariate key.

  - `covariate_table` (optional, default `"ind_meta"`): table containing
    the covariate column; it must have exactly one row per individual.

  - `poly_order` (optional): polynomial basis order.

  - `poly_scale_min`, `poly_scale_max` (optional): Legendre scaling
    bounds.

  - `component_names` (optional, default `"total"`): which genetic value
    of the contributor this row reads from `ind_tgv`. `"total"` is the
    total genetic value (the `ind_tgv_total` view: additive, dominance
    and every other component). A comma-separated list of
    `ind_tgv.component_name` values (`"additive"`, `"dominance"`,
    `"indicator"`, `"interaction"`, e.g. `"additive"` or
    `"additive,dominance"`) reads their sum; a listed component the
    trait's model has no terms for contributes 0. Anything else,
    `"total"` mixed with other names, or a duplicate is an error.

  - `group_column` (optional): column defining group membership.

  - `group_table` (optional, default `"ind_meta"`): table containing
    `group_column`.

  - `aggregation` (optional, default `"sum"`): `"sum"` or `"mean"` for
    group contributors.

  `NULL` (default) → simple single-self trait; `phenotype_components`
  not written. Mutually exclusive with `formula_tbv`.

- formula_tbv:

  Character. DSL shorthand for assembling a composite genetic value from
  component traits already in `trait_meta`. Every reference reads the
  contributor's **total** genetic value (`ind_tgv_total`). A bare trait
  symbol (e.g. `"WWD"`) means the individual's own (`"self"`) value;
  contributor roles can also be given explicitly as function calls:
  `self(trait)`, `dam(trait)`, `sire(trait)`, `group_sum(trait, col)`,
  and `group_mean(trait, col)` (`col` = grouping column in `ind_meta`,
  e.g. pen or litter). These are combined with the arithmetic operators
  `+`, `-`, `*`, `/`, and parentheses (e.g. `"WWD + dam(WWM)"`,
  `"ADG_direct + group_sum(ADG_social, pen_id)"`). Mutually exclusive
  with `components`. Not valid with `type = "derived_formula"`.

- formula:

  Character. Arithmetic expression evaluated over already- recorded
  phenotype values to produce a derived phenotype (e.g. `"ADFI / ADG"`
  for feed conversion ratio). Phenotype names reference `ind_phenotype`
  records (`pheno_number = 1`) produced earlier in the same
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  call — component phenotypes must already be recorded, and
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  topologically sorts multiple `derived_formula` phenotypes so
  dependencies are computed first. Operators `+`, `-`, `*`, `/`, `^` and
  parentheses are supported, plus the math functions `sqrt`, `log`,
  `log2`, `log10`, `exp`, `abs`, `round`, `ceiling`, `floor`, `sign`,
  `trunc`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`. Non-finite
  results (`Inf`/`NaN`) are converted to `NA` with a warning. Required
  when `type = "derived_formula"`. Not valid otherwise.

- missing_component_action:

  Character. What to do when an individual is missing one or more
  required composite components (e.g. no group assignment for a
  `"group"` contributor, or a missing dam/sire genetic value). `"skip"`
  (default) excludes the individual from `ind_phenotype` and emits a
  warning with a count. `"error"` stops with an informative message
  listing affected individuals. Stored in `phenotype_meta` so the
  behaviour is consistent across all
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  calls for this phenotype. Note: this is unrelated to
  `null_class_action` (set via
  [`define_effect_fixed_class()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_class.md)),
  which handles `NULL` levels for fixed-class covariate effects, and
  does not affect random-effect draws (new levels always get a fresh
  draw).

- condition_change_action:

  Character. Applies only when this phenotype is in a residual
  covariance block with a `condition_column` (see
  [`define_residual_cov()`](https://austin-putz.github.io/tidybreed/reference/define_residual_cov.md))
  and a correlated phenotype's residual was stored under a **different**
  condition level than the one the current record resolves to — e.g. an
  animal moved farms between the two records. `"error"` (default) stops,
  because no covariance is defined between the two strata.
  `"independent"` drops the incompatible stored residual from the
  conditioning set (stored residuals from the same stratum still
  condition the draw) and warns with a count. Stored in
  `phenotype_meta`; every phenotype in one residual block must carry the
  same value (D6), so this argument only sets it while the phenotype is
  still a block of one. Once the block has two or more members the value
  is block-scoped — change it with
  [`define_condition_change_action()`](https://austin-putz.github.io/tidybreed/reference/define_condition_change_action.md),
  which writes every member in one transaction and leaves the rest of
  their `phenotype_meta` rows alone. An immutable condition column such
  as `sex` never triggers either action.

- overwrite:

  Logical. If `TRUE` and a phenotype with the same name already exists,
  replace its rows in `phenotype_meta` and `phenotype_components`.
  Default `FALSE` errors on duplicate.

## Value

The modified `tidybreed_pop` (invisibly).

## See also

[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md),
[`define_residual_cov()`](https://austin-putz.github.io/tidybreed/reference/define_residual_cov.md),
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# ── Simple continuous trait ──────────────────────────────────────────────
pop <- pop |>
  define_trait("ADG") |>
  get_table("genome_meta") |>
  define_additive_effects("ADG", G = 100) |>
  define_phenotype("ADG",
    type         = "continuous",
    mean         = 850,
    residual_var = 120)

# ── Count trait with clipping bounds ────────────────────────────────────
pop <- pop |>
  define_trait("NW") |>
  get_table("genome_meta") |>
  define_additive_effects("NW", G = 2) |>
  define_phenotype("NW",
    type         = "count",
    mean         = 10,
    min_value    = 1,
    max_value    = 30,
    residual_var = 8)

# ── Categorical trait (binary via prevalence) ────────────────────────────
pop <- pop |>
  define_trait("mort") |>
  get_table("genome_meta") |>
  define_additive_effects("mort", G = 0.05) |>
  define_phenotype("mort",
    type       = "categorical",
    prevalence = 0.05,
    cat_names  = c("Alive", "Dead"))

# ── Maternal composite via components data frame: WW = WWD (self) + WWM (dam) ──
pop <- pop |>
  define_phenotype("WW",
    type         = "continuous",
    mean         = 230,
    residual_var = 180,
    components   = tibble::tribble(
      ~source_trait_name, ~contributor_type,
      "WWD",              "self",
      "WWM",              "dam"
    ))

# ── Maternal composite via formula_tbv shorthand (equivalent to above) ──
pop <- pop |>
  define_phenotype("WW2",
    type         = "continuous",
    mean         = 230,
    residual_var = 180,
    formula_tbv  = "WWD + dam(WWM)")

# ── SGE (social genetic effects): ADG = direct (self) + social group sum ──
pop <- pop |>
  define_phenotype("ADG_sge",
    type         = "continuous",
    mean         = 850,
    residual_var = 100,
    formula_tbv  = "ADG_direct + group_sum(ADG_social, pen_id)")

# ── Derived formula: FCR computed from already-recorded ADFI and ADG ────
# (Define ADFI and ADG first, then derive FCR — no TBV or residual needed)
pop <- pop |>
  define_phenotype("FCR",
    type    = "derived_formula",
    formula = "ADFI / ADG")
} # }
```

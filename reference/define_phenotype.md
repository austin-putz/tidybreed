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
  formula_tgv = NULL,
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
  phenotype values for the same individuals; they have no genetic value,
  no residual variance, and no QTL of their own.

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
  (two categories), the fraction expected strictly above the threshold.
  Mutually exclusive with `thresholds`. The threshold is
  `mean + qnorm(1 - prevalence) * sqrt(V)`, where `V` is the variance of
  the whole liability:

  - the trait's stored genetic targets (`trait_var_comp`,
    population-wide rows): the sum of the `additive`, `dominance` and
    `additive_by_additive` diagonals, each counted only if the trait's
    model has terms of that kind;

  - the stored variance of every named random effect
    ([`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md);
    `normal` and `uniform` are centred with that variance);

  - the unconditional residual variance.

  This is a **Gaussian approximation** at a reference population in
  Hardy-Weinberg and linkage equilibrium. It is exact only when the
  liability is normal: a few large QTL make the genetic value discrete,
  and the realised prevalence then differs even in an infinite reference
  population. Summing the component targets also assumes they are
  orthogonal (no covariance between additive, dominance and A x A
  values), which holds for statistical coding at one base, not under LD.
  Fixed effects shift the liability and are not included: the prevalence
  is for records whose fixed effects are 0. A selected or line-specific
  population also differs. For an exact fraction in a known population,
  compute cutpoints from it and pass `thresholds`.

  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  refuses the threshold, before any draw or write, when no stored target
  can describe the liability:

  - a term not owned by `"generated"`. A generator
    ([`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md))
    always calibrates its terms to the stored target; terms written with
    [`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md)
    carry values nothing checked against one;

  - a kind of term with no stored target, or terms outside the three
    kinds (an `indicator` surface, other interactions);

  - generated variants for two parent-of-origin scopes at one line
    (paternal-only plus maternal-only, or common plus a parent-only
    fallback). Each was calibrated to the target alone; together they
    have a different variance;

  - a `gamma` random effect (its mean is `sqrt(variance)`);

  - a residual with conditional strata (the marginal variance would
    depend on the levels' frequencies);

  - a total variance of 0, where no cutpoint gives a fraction.

  Not valid for composite phenotypes (`components` or `formula_tgv`):
  their genetic liability combines several traits and contributors,
  which no stored variance describes. Give `thresholds` instead.

- thresholds:

  Numeric vector of length K−1 for K ordered categories: finite
  liability cutpoints in strictly ascending order. A record is in
  category `k + 1` when its liability is strictly above cutpoint `k`; a
  liability exactly on a cutpoint stays in the lower category. Mutually
  exclusive with `prevalence`.

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
  not written. Mutually exclusive with `formula_tgv`.

- formula_tgv:

  Character. DSL shorthand for assembling a composite genetic value from
  component traits already in `trait_meta`. A bare trait symbol (e.g.
  `"WWD"`) is the individual's own (`"self"`) value; contributor roles
  are given as calls:

  - `self(trait)`, `dam(trait)`, `sire(trait)`: one positional trait;

  - `group_sum(trait, col)`, `group_mean(trait, col)`: the sum or mean
    over the individual's group-mates (the *other* individuals with the
    same value of `col`; a group of one gives `0`).

  Every call takes an optional named `component =`: which genetic value
  of the contributor to read from `ind_tgv` — `"total"` (the default,
  the `ind_tgv_total` view: additive, dominance and every other
  component), or one component, `"additive"` (the breeding value for
  generated effects), `"dominance"`, `"indicator"` or `"interaction"`; a
  component the trait's model has no terms for reads 0. The group calls
  also take a named `table =` (default `"ind_meta"`): the table holding
  `col`, with exactly one row per individual. Both must be named, e.g.
  `"WWD + dam(WWM, component = \"additive\")"` or
  `"ADG_direct + group_sum(ADG_social, pen, table = \"pens\")"`.

  References combine with `+`, `-`, `*`, `/`, `^`, parentheses, numbers
  and the math functions listed under `formula`. Anything else — another
  function, an unknown or positional extra argument, a component outside
  that list, a `col` or `table` that is not a plain identifier, or a
  table or column that does not exist yet — is an error here, before
  anything is written. Mutually exclusive with `components`. Not valid
  with `type = "derived_formula"`. A constant expression gives every
  individual that value. At
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
  a result that is `Inf`, `-Inf` or `NaN` (division by zero, overflow, a
  function outside its domain) is an error naming the individuals, and
  nothing is written; a missing contributor is
  `missing_component_action`'s business, as before.

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

# ── Maternal composite via formula_tgv shorthand (equivalent to above) ──
pop <- pop |>
  define_phenotype("WW2",
    type         = "continuous",
    mean         = 230,
    residual_var = 180,
    formula_tgv  = "WWD + dam(WWM)")

# ── SGE (social genetic effects): ADG = direct (self) + social group sum ──
pop <- pop |>
  define_phenotype("ADG_sge",
    type         = "continuous",
    mean         = 850,
    residual_var = 100,
    formula_tgv  = "ADG_direct + group_sum(ADG_social, pen_id)")

# ── Derived formula: FCR computed from already-recorded ADFI and ADG ────
# (Define ADFI and ADG first, then derive FCR — no genetic value or residual needed)
pop <- pop |>
  define_phenotype("FCR",
    type    = "derived_formula",
    formula = "ADFI / ADG")
} # }
```

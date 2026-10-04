# Define a variance-covariance matrix for any named effect

Single entry point for storing all variance and covariance data in
tidybreed. Routes to `trait_var_comp` for genetic effects and to
`phenotype_var_comp` for phenotype-level effects.

Common `effect_name` values:

- `"additive"` — additive genetic (co)variances (G matrix). Written to
  `trait_var_comp`. Used by
  [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
  when rescaling to target variance and as the sampling distribution for
  multi-trait draws.

- `"dominance"`, `"additive_by_additive"` — reserved genetic effects; no
  generator calibrates them yet. Written to `trait_var_comp`. Row/column
  names are trait names.

- `"residual"` — residual (co)variances (R matrix). Routed to
  `phenotype_var_comp` with `effect_name = "residual"`. Row/column names
  are phenotype names. Equivalent to calling
  [`define_residual_cov()`](https://austin-putz.github.io/tidybreed/reference/define_residual_cov.md)
  with `condition_column = NULL`. Use this for a multi-phenotype
  correlated residual matrix; for a single scalar residual use
  `residual_var` in
  [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
  instead.

- Any named random effect (`"hys"`, `"litter"`, `"pen"`, …) — written to
  `phenotype_var_comp`. Must match the `effect_name` used in
  [`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md).
  Row/column names are phenotype names. Each level of the effect (each
  pen) then carries one draw per phenotype with this covariance,
  realized sequentially: whichever phenotype
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  generates first for a level draws marginally, and the others are later
  drawn conditional on what the level has stored — however many calls
  apart (see
  [`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md)).

`define_effect_cov_matrix()` can be called **before**
[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md)
or
[`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md)
— no prior setup is required.

All n² pairs are stored. For `phenotype_var_comp` effects the names form
a *covariance block* that is declared in one call, as a complete matrix:
a call that names a fragment or a strict subset of an existing block is
an error, the matrix must be positive semi-definite, and a block cannot
be redefined once draws exist under it (the error gives the
[`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md)
call that clears them). In a block of two or more phenotypes every
[`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md)
row for the effect must use `distribution = "normal"` and read the same
`(source_column, source_table)`. See
[`define_residual_cov()`](https://austin-putz.github.io/tidybreed/reference/define_residual_cov.md)
for the full rules. A rejected call changes nothing.

**Genetic blocks are written once.** A genetic block (`"additive"`,
`"dominance"`, `"additive_by_additive"`) is validated as positive
semi-definite, stored at full double precision, and never overwritten:
if any row already exists for that `effect_name`, any of the named
traits and the same `line_name`, the call is an error, even when the
matrix is identical. The error gives the
[`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md)
call that clears the stored block. `trait_var_comp` is the single source
of generation targets; the effect generators read it and never overwrite
it either.

`"additive_by_dominance"` and `"dominance_by_dominance"` are reserved
for future generators and refused. `"total"`, `"unpartitioned"` and
`"between_components"` are output names of the variance extractor and
refused as input.

## Usage

``` r
define_effect_cov_matrix(
  pop,
  effect_name,
  cov_matrix,
  trait_name = NULL,
  line_name = NULL,
  tol = 1e-09
)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- effect_name:

  Character. Label for the variance component, e.g. `"additive"`,
  `"residual"`, `"hys"`.

- cov_matrix:

  A numeric square matrix, or a single number when one trait/phenotype
  is named. Must be symmetric within `tol`. A named matrix must carry
  the same names as `trait_name`, in the same order (it is never
  relabelled); an unnamed one is taken in `trait_name` order.

- trait_name:

  Character vector of trait/phenotype names (length ==
  `nrow(cov_matrix)`). Optional when the matrix has names.

- line_name:

  Character or `NULL` (default). Genetic effects only: the line whose
  generation target this is. `NULL` is the population-wide target, which
  a line without its own block falls back to.

- tol:

  Numeric. Tolerance for symmetry check (default `1e-9`).

## Value

The modified `tidybreed_pop` (invisibly).

## See also

[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md),
[`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md),
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md),
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Additive genetic covariance matrix → trait_var_comp
G <- matrix(c(100, -20, -20, 50), 2, 2,
            dimnames = list(c("ADG", "BF"), c("ADG", "BF")))
pop <- pop |>
  define_effect_cov_matrix("additive", G)

# One trait: a number is a 1 x 1 matrix
pop <- pop |>
  define_effect_cov_matrix("additive", 0.25, trait_name = "WW")

# Residual → phenotype_var_comp (effect_name = "residual")
R <- matrix(c(30, 5, 5, 10), 2, 2,
            dimnames = list(c("ADG", "BF"), c("ADG", "BF")))
pop <- pop |>
  define_effect_cov_matrix("residual", R)

# Multi-phenotype HYS covariance → phenotype_var_comp (effect_name = "hys")
R_hys <- matrix(c(0.2, 0.05, 0.05, 0.3), 2, 2,
                dimnames = list(c("ADG", "BF"), c("ADG", "BF")))
pop <- pop |>
  define_effect_cov_matrix("hys", R_hys)
} # }
```

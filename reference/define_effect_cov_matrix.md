# Define a variance-covariance matrix for any named effect

Single entry point for storing all variance and covariance data in
tidybreed. Routes to `trait_var_comp` for genetic effects and to
`phenotype_var_comp` for phenotype-level effects.

Common `effect_name` values:

- `"gen_add"` — additive genetic (co)variances (G matrix). Written to
  `trait_var_comp`. Used by
  [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
  when rescaling to target variance and as the sampling distribution for
  multi-trait draws.

- `"dominance"`, `"epistasis"` — future genetic effects. Written to
  `trait_var_comp`. Row/column names are trait names.

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
  Row/column names are phenotype names.

`define_effect_cov_matrix()` can be called **before**
[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md)
or
[`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md)
— no prior setup is required.

All n² pairs are stored. Previous entries for this `effect_name` × names
combination are replaced.

## Usage

``` r
define_effect_cov_matrix(
  pop,
  effect_name,
  cov_matrix,
  trait_names = NULL,
  tol = 1e-09
)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- effect_name:

  Character. Label for the variance component, e.g. `"gen_add"`,
  `"residual"`, `"hys"`.

- cov_matrix:

  A numeric square matrix. Must be symmetric within `tol`. Row and
  column names are used as trait/phenotype names when `trait_names` is
  not supplied.

- trait_names:

  Optional character vector of trait/phenotype names (length ==
  `nrow(cov_matrix)`). Overrides the matrix's `rownames` / `colnames`.

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
  define_effect_cov_matrix("gen_add", G)

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

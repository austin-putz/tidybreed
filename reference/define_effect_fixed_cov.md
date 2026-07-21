# Define a continuous covariate (regression) effect in a phenotype model

Inserts a row into `trait_effects` representing a polynomial regression
term. At phenotyping time the contribution per individual is
`slope * (source_column_value - center)^poly_order`. When `center` is
`NULL` the mean of `source_column` is computed and stored automatically.

## Usage

``` r
define_effect_fixed_cov(
  pop,
  phenotype_name,
  effect_name,
  source_column,
  slope,
  center = NULL,
  poly_order = 1L,
  source_table = "ind_meta",
  overwrite = FALSE
)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- phenotype_name:

  Character. Name of an existing phenotype in `phenotype_meta`.

- effect_name:

  Character. Unique label for this effect within the phenotype.

- source_column:

  Character. Column in `source_table` holding the continuous predictor
  values.

- slope:

  Numeric scalar. Regression coefficient.

- center:

  Numeric scalar or `NULL`. Value subtracted before raising to
  `poly_order`. If `NULL`, auto-computed as the column mean.

- poly_order:

  Integer ≥ 1. Polynomial order. Default `1` gives the standard linear
  `slope * (x - center)` term. Use `2` for quadratic, etc.

- source_table:

  Character. Database table containing `source_column`. Default
  `"ind_meta"`.

- overwrite:

  Logical. Replace an existing effect with the same name.

## Value

The modified `tidybreed_pop` (invisibly).

## See also

[`define_effect_intercept()`](https://austin-putz.github.io/tidybreed/reference/define_effect_intercept.md),
[`define_effect_fixed_class()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_class.md),
[`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md),
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Linear regression on age (center auto-computed)
pop <- pop |>
  define_effect_fixed_cov("ADG", "age",
    source_column = "age_days",
    slope = 2.5)

# Quadratic regression on DIM (dairy)
pop <- pop |>
  define_effect_fixed_cov("milk_td", "dim_sq",
    source_column = "dim",
    slope  = -0.02,
    center = 150,
    poly_order = 2L)
} # }
```

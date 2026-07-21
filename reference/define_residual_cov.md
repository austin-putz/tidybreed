# Define residual covariance entries for observed phenotypes

Writes rows to `phenotype_var_comp` (with `effect_name = "residual"`)
representing the residual (co)variance matrix for one or more
phenotypes, optionally conditioned on a group variable (e.g. farm, sex).
Both `(i,j)` and `(j,i)` pairs are stored.

Three typical call patterns:

1.  **Called by
    [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
    internally** when `residual_var` is supplied (scalar diagonal,
    unconditional).

2.  **Called by
    [`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md)**
    when `effect_name = "residual"` to store a full multi-trait
    unconditional R matrix.

3.  **Called directly** to add group-specific (heterogeneous) residual
    rows after
    [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
    has written the unconditional default.

## Usage

``` r
define_residual_cov(
  pop,
  phenotype_names,
  cov_matrix,
  condition_column = NULL,
  condition_table = "ind_meta",
  condition_level = NULL
)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- phenotype_names:

  Character vector of phenotype names (must match `rownames(cov_matrix)`
  and `colnames(cov_matrix)`). For a single phenotype, a scalar is
  accepted.

- cov_matrix:

  Numeric matrix. Must be symmetric with `dimnames(cov_matrix)` matching
  `phenotype_names`. For a single phenotype the matrix is `1×1`.

- condition_column:

  Character or `NULL`. Column in `condition_table` used to look up group
  membership at phenotype time. `NULL` (default) stores an unconditional
  default row (the fallback when no group-specific row exists).

- condition_table:

  Character. Table containing `condition_column`. Default `"ind_meta"`.

- condition_level:

  Character or `NULL`. Specific level of `condition_column` this row
  applies to. `NULL` (default) = unconditional fallback row.

## Value

The modified `tidybreed_pop` (invisibly).

## See also

[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md),
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Heterogeneous residual by sex for turkey BW
pop <- pop |>
  define_phenotype("BW_turkey",
    type         = "continuous",
    mean         = 8000,
    residual_var = 600) |>   # unconditional default
  define_residual_cov("BW_turkey",
    cov_matrix       = matrix(900, 1, 1, dimnames = list("BW_turkey","BW_turkey")),
    condition_column = "sex",
    condition_level  = "M") |>
  define_residual_cov("BW_turkey",
    cov_matrix       = matrix(400, 1, 1, dimnames = list("BW_turkey","BW_turkey")),
    condition_column = "sex",
    condition_level  = "F")
} # }
```

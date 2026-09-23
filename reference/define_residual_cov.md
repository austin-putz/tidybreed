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
    when `effect_name = "residual"` to store a full multi-phenotype
    unconditional R matrix.

3.  **Called directly** to declare group-specific (heterogeneous)
    residual strata, one call per `condition_level`.

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

  Numeric matrix with `dimnames(cov_matrix)` matching `phenotype_names`.
  For a single phenotype the matrix is `1×1`.

- condition_column:

  Character or `NULL`. Column in `condition_table` used to look up group
  membership at phenotype time. `NULL` (default) declares the
  unconditional stratum. Must be supplied together with
  `condition_level`.

- condition_table:

  Character. Table containing `condition_column`. Default `"ind_meta"`.
  Ignored (stored as `NULL`) for the unconditional stratum.

- condition_level:

  Character or `NULL`. Level of `condition_column` this stratum applies
  to. `NULL` (default) = unconditional stratum.

## Value

The modified `tidybreed_pop` (invisibly).

## A block is declared in one call, as a complete matrix

The phenotypes that share any stored residual pair row form a
*covariance block*; a pair written as an explicit `0` still joins its
two phenotypes. A call must name a whole block: declaring `{A, B}` and
then `{B, C}` is an error (the block would be `{A, B, C}` with
`Cov(A, C)` undeclared), and so is redeclaring `{A, B}` or `A` alone
once `{A, B, C}` exists. Redeclare the complete block instead, writing
`0` for uncorrelated pairs. The matrix must be symmetric, finite and
positive semi-definite; a rejected call changes nothing.

**Strata.** Conditional calls (`condition_column` + `condition_level`)
add strata to the block. Every stratum names the same phenotypes and a
block has one `condition_column`; to grow a block that already has
several strata, clear its rows from `phenotype_var_comp` with
[`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md)
and redeclare each stratum.

**How the block is sampled.**
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
draws each record's residual from this block conditional on the
residuals the same individual has already realized for the block's other
phenotypes at the same `pheno_number` — in the same call or any earlier
one — so the declared covariance holds whether the phenotypes are
recorded together or a hundred simulated days apart with culling in
between. With strata, each record draws from the stratum its
`condition_column` value selects, falling back to the unconditional
stratum when the value is `NULL` or matches none (an error if there is
no unconditional stratum).

**Realized draws lock the block.** Once any `ind_phenotype` row of a
member has a non-`NULL` `residual_value`, the block cannot be redefined;
the error gives the
[`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md)
call that clears those rows. Every phenotype in the block that is
already defined must carry the same `condition_change_action` (see
[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)).

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

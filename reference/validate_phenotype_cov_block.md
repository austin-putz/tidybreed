# Validate a `phenotype_var_comp` block declaration

Runs D1 (completeness, matrix), the residual strata rules, D6, the §5.6
named-effect checks and the D3 realization lock against the *current*
table state — call it before deleting anything. Returns the block
members (`== phenotype_names` when it returns at all).

## Usage

``` r
validate_phenotype_cov_block(
  conn,
  effect_name,
  phenotype_names,
  cov_matrix,
  condition_column = NULL,
  condition_table = NULL,
  condition_level = NULL,
  caller = "define_residual_cov()",
  tol = 1e-08
)
```

## Arguments

- conn:

  A DBI connection.

- effect_name:

  `"residual"` or a named random effect.

- phenotype_names:

  Character vector; the block being declared.

- cov_matrix:

  Numeric matrix with `dimnames` equal to `phenotype_names`.

- condition_column, condition_table, condition_level:

  The stratum (`NULL`s for the unconditional stratum). Residual only.

- caller:

  Prefix for error messages.

- tol:

  Relative symmetry / PSD tolerance.

## Value

The sorted block members, invisibly.

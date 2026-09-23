# Replace one stratum of a block: validate, delete, insert

No transaction management — the caller owns the transaction. Writes
through `duckdb_register()` + `INSERT`, which does not touch R's RNG.

## Usage

``` r
.pvc_write_block(
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

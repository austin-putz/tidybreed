# Write one stratum of a covariance block in its own transaction

The single write path for `phenotype_var_comp`, used by
[`define_residual_cov()`](https://austin-putz.github.io/tidybreed/reference/define_residual_cov.md),
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md)
and (inside its own transaction, via
[`.pvc_write_block()`](https://austin-putz.github.io/tidybreed/reference/dot-pvc_write_block.md))
[`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md).

## Usage

``` r
write_phenotype_cov_block(
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

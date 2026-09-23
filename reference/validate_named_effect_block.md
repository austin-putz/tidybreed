# §5.6: the `phenotype_effects` rows of a named-effect block are compatible

Applies to blocks of two or more phenotypes. Every row for
`(effect_name, member)` must be a `random` effect with
`distribution = "normal"` and the same `(source_column, source_table)`.

## Usage

``` r
validate_named_effect_block(
  conn,
  effect_name,
  phenotype_names,
  pending = NULL,
  caller = "define_effect_cov_matrix()"
)
```

## Arguments

- conn:

  A DBI connection.

- effect_name:

  Character scalar (not `"residual"`).

- phenotype_names:

  The block members.

- pending:

  Optional one-row data frame with columns `phenotype_name`,
  `effect_class`, `source_column`, `source_table`, `distribution` for a
  `phenotype_effects` row about to be written; it replaces any stored
  row of the same phenotype in the comparison.

- caller:

  Prefix for the error message.

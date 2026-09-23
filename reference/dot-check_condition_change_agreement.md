# D6: `condition_change_action` agrees across a residual block

D6: `condition_change_action` agrees across a residual block

## Usage

``` r
.check_condition_change_agreement(
  conn,
  phenotype_names,
  pending = NULL,
  caller = "define_residual_cov()"
)
```

## Arguments

- conn:

  A DBI connection.

- phenotype_names:

  The block members.

- pending:

  Optional `list(phenotype_name =, condition_change_action =)` for a
  `phenotype_meta` row about to be written; it replaces any stored row
  of the same name in the comparison.

- caller:

  Prefix for the error message.

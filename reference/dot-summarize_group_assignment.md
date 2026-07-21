# Emit group-level diagnostics after assignment.

No-op when `quiet = TRUE`.

## Usage

``` r
.summarize_group_assignment(values, n_total, table_name, col_name, quiet)
```

## Arguments

- values:

  Vector of assigned group values (`NA` for unassigned rows), in the
  order they were written.

- n_total:

  Total number of filtered rows considered for assignment.

- table_name:

  Target table name, used only in the message text.

- col_name:

  Target column name, used only in the message text.

- quiet:

  Logical. If `TRUE`, suppress all messages (no-op).

## Value

`NULL`, invisibly.

# Assign each table to a display group and a pipeline sort position

Assign each table to a display group and a pipeline sort position

## Usage

``` r
.schema_group_of(table_names)
```

## Arguments

- table_names:

  Character vector of table names.

## Value

A data.frame with `table_group` and `sort_key` (integer), aligned to
`table_names`.

# Walk the AST of a formula_tbv expression, returning structured trait refs.

Each occurrence of a bare symbol or DSL function call gets a unique
placeholder name (e.g. ".tbv_self_WWD_1", ".tbv_dam_WWM_2") so that the
same trait can appear multiple times with distinct pre-fetched vectors.

## Usage

``` r
.walk_formula_tbv_ast(expr)
```

## Arguments

- expr:

  Parsed R expression (from `parse()[[1]]`).

## Value

A list: \$trait_refs: list of lists, each with: - trait: character trait
name - type: "self", "dam", "sire", "group_sum", or "group_mean" - col:
group column name (NA for non-group types) - table: group table name (NA
for non-group types) - placeholder: unique R symbol name for the
pre-fetched vector \$has_scalar_constant: logical

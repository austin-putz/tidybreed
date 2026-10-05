# Walk a formula_tgv expression: validate it, collect its references, and replace each with a placeholder

One depth-first pass. Each bare trait symbol or DSL call becomes one
reference and is replaced, in the returned expression, by that
reference's unique placeholder symbol (`.tgv_1`, `.tgv_2`, ...), so the
same trait can appear several times with different contributors,
components or tables.

## Usage

``` r
.walk_formula_tgv_ast(expr)
```

## Arguments

- expr:

  Parsed R expression (from `.parse_formula_tgv()`).

## Value

A list: \$trait_refs: list of lists, each with: - trait: character trait
name - type: "self", "dam", "sire", "group_sum", or "group_mean" - col:
group column name (NA for non-group types) - table: group table name (NA
for non-group types) - component: "total" or one
ind_tgv.component_name - placeholder: unique R symbol name for the
pre-fetched vector - call: the reference as written, for messages
\$expr: `expr` with every reference replaced by its placeholder

## Details

The DSL: a bare symbol is `self(trait)`; `self(trait)`, `dam(trait)` and
`sire(trait)` take one positional trait; `group_sum(trait, col)` and
`group_mean(trait, col)` take a trait and a group column. All five take
an optional named `component =` (one of
[TGV_COMPONENT_NAMES](https://austin-putz.github.io/tidybreed/reference/TGV_COMPONENT_NAMES.md)
or `"total"`, the default); the group functions also take a named
`table =` (default `"ind_meta"`). Trait, column and table are symbols or
strings; `col` and `table` must be SQL identifiers. Any other argument,
any call outside the DSL, the arithmetic operators and the math
whitelist, and any constant that is not a number is an error naming the
offending call.

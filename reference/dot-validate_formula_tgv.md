# Validate a formula_tgv string at define_phenotype() time.

Parses the DSL formula and walks it with
[`.walk_formula_tgv_ast()`](https://austin-putz.github.io/tidybreed/reference/dot-walk_formula_tgv_ast.md),
which refuses any call outside the DSL, the arithmetic operators and the
math whitelist, and any DSL call with an argument it does not take. Then
checks every referenced trait against `trait_meta` (with close-match
suggestions via [`agrep()`](https://rdrr.io/r/base/agrep.html)), and
every group `table` / `col` against the database: the table must exist
and hold `id_ind` and the column.

## Usage

``` r
.validate_formula_tgv(conn, formula_tgv)
```

## Arguments

- conn:

  DBI connection.

- formula_tgv:

  Character. The DSL formula string.

## Value

Invisible NULL on success. Stops on error; warns for scalar constants.

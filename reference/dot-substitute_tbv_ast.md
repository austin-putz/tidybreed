# Replace DSL calls and bare symbols in the formula_tbv AST with placeholders.

Traverses in the same depth-first order as .walk_formula_tbv_ast(),
matching each node to its pre-assigned placeholder from trait_refs
(using sequential .used flags so repeated occurrences are handled
correctly).

## Usage

``` r
.substitute_tbv_ast(expr, trait_refs)
```

## Arguments

- expr:

  Parsed R expression.

- trait_refs:

  List from .walk_formula_tbv_ast()\$trait_refs.

## Value

Modified expression with placeholder symbols.

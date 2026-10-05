# Check a derived formula's grammar and return the phenotypes it reads

The expression is [`eval()`](https://rdrr.io/r/base/eval.html)ed at
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
time, so only phenotype names, numbers, the arithmetic operators and the
math whitelist are accepted. Any other call or constant is an error
naming it.

## Usage

``` r
.check_derived_formula(expr, formula)
```

## Arguments

- expr:

  Parsed R expression.

- formula:

  The formula string, for messages.

## Value

Character vector of the phenotype names referenced (unique).

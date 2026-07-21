# Validate a formula_tbv string at define_phenotype() time.

Parses the DSL formula, walks the AST to collect trait references, and
checks that every trait name is present in trait_meta. Provides
close-match suggestions via agrep() for typos.

## Usage

``` r
.validate_formula_tbv(formula_tbv, known_traits)
```

## Arguments

- formula_tbv:

  Character. The DSL formula string.

- known_traits:

  Character vector of trait names from trait_meta.

## Value

Invisible NULL on success. Stops on error; warns for scalar constants.

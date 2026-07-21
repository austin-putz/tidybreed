# Validate a derived formula string at define_phenotype() time.

Symbols that are not in phenotype_meta generate a warning (not an
error), allowing config-first workflows where components are defined
before their dependents. A hard error at add_phenotype() time fires if
still missing.

## Usage

``` r
.validate_derived_formula(formula, known_phenos)
```

## Arguments

- formula:

  Character. The arithmetic formula string.

- known_phenos:

  Character vector of phenotype names from phenotype_meta.

## Value

Invisible NULL on success.

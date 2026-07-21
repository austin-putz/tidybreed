# Rescale QTL effects to hit a target additive variance using Falconer formula

Rescale QTL effects to hit a target additive variance using Falconer
formula

## Usage

``` r
rescale_effects_to_target(qtl_tf, qtl_effects, target_add_var, p_base)
```

## Arguments

- qtl_tf:

  Logical mask of QTL loci, length `n_loci`.

- qtl_effects:

  Numeric effects at QTL loci, length `sum(qtl_tf)`.

- target_add_var:

  Target additive variance.

- p_base:

  Numeric vector of base allele frequencies, length `n_loci`.

## Value

Rescaled `qtl_effects` vector.

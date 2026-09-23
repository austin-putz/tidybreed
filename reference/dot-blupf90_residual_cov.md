# The unconditional residual (co)variance matrix of a set of traits

Assembled from the residual covariance blocks in `phenotype_var_comp`
([`find_covariance_blocks()`](https://austin-putz.github.io/tidybreed/reference/find_covariance_blocks.md)):
traits in different blocks are independent, so the matrix is
block-diagonal with explicit zeros between blocks. A trait in no block,
or in a block with only conditional strata, is an error — BLUPF90 takes
one residual matrix.

## Usage

``` r
.blupf90_residual_cov(pop, trait)
```

## Value

Numeric `trait` x `trait` matrix with dimnames.

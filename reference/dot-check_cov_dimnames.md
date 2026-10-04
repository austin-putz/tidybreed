# Check a covariance matrix's names against `trait_name`, never relabelling

A single number is a 1 x 1 matrix when one name is given. A named matrix
(row names, and column names when present) must equal `trait_name` in
order; an unnamed one is taken in `trait_name` order. With no
`trait_name`, the row names are the names. Before 0.73.0 the names were
assigned over whatever the matrix carried, so a matrix named
`c("BF", "ADG")` passed with `trait_name = c("ADG", "BF")` was stored
with its rows the wrong way round.

## Usage

``` r
.check_cov_dimnames(x, trait_name, arg = "cov_matrix")
```

## Value

The matrix with `dimnames = list(names, names)`.

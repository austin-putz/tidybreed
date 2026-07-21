# Sample residuals for one or more phenotypes

Sample residuals for one or more phenotypes

## Usage

``` r
sample_residuals(n_ind, residual_var, R = NULL)
```

## Arguments

- n_ind:

  Integer. Number of individuals.

- residual_var:

  Named numeric vector of per-phenotype residual variances.

- R:

  Optional covariance matrix with row/col names matching
  `names(residual_var)`.

## Value

Numeric matrix of shape `n_ind × length(residual_var)`, columns named by
phenotype.

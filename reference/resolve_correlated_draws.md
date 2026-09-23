# Draw coordinates of a Gaussian block conditional on observed coordinates

For `n` entities sharing one covariance matrix `covariance` over a block
of coordinates, draws the `sample_coordinates` of every entity from the
multivariate normal conditional on that entity's observed coordinates.
An entity's observed set is whatever is non-`NA` in its row of
`observed`; entities with different observed sets are grouped by pattern
and the conditional mean coefficients and covariance are computed once
per pattern. Coordinates that are neither observed nor sampled are
latent and are never realized.

The function is pure apart from RNG: it never reads or writes a database
and knows nothing about individuals, strata or tables.

## Usage

``` r
resolve_correlated_draws(
  covariance,
  sample_coordinates,
  entity_keys,
  observed = NULL,
  tolerance = NULL
)
```

## Arguments

- covariance:

  Numeric matrix with dimnames: one stratum of a block, as returned by
  [`find_covariance_blocks()`](https://austin-putz.github.io/tidybreed/reference/find_covariance_blocks.md).

- sample_coordinates:

  Character vector of the coordinates to draw, a subset of the block's
  names with no duplicates. May be empty.

- entity_keys:

  Opaque entity identifiers: a vector or list (one entity per element)
  or a data frame (one entity per row). Only the count is used, and it
  fixes `n`. Rows of the result follow its order.

- observed:

  `NULL`, or a numeric matrix / data frame with `n` rows whose column
  names are block coordinates not in `sample_coordinates`. `NA` means
  "not observed for this entity"; every non-`NA` value must be finite.
  Stored and caller-fixed values both go here.

- tolerance:

  Relative tolerance; see the numerical contract. `NULL` uses the
  default.

## Value

A numeric `n` x `length(sample_coordinates)` matrix with the sample
coordinates as column names, rows in `entity_keys` order.

## Numerical contract

- `tolerance` is *relative*. The absolute tolerance is
  `tolerance * lambda_max(covariance)`, so the same relative
  perturbation is treated the same way at variance `1e-6` and `1e6`. The
  default is `nrow(covariance) * sqrt(.Machine$double.eps)`.

- `covariance` must be square, finite, symmetric within tolerance, with
  identical unique row and column names, and positive semi-definite
  within tolerance. Perfect correlation and zero-variance coordinates
  are valid.

- The observed block `R_oo` is inverted by Cholesky when it is positive
  definite and by an eigen pseudoinverse when it is only semi-definite.
  In the singular case every observed vector must lie in the support of
  its Gaussian: the component of `e_o` outside the range of `R_oo` must
  be within `tolerance * max(||e_o||, sqrt(lambda_max))` of zero. A
  vector outside the support (a non-zero value for a zero-variance
  coordinate, or two perfectly correlated coordinates that disagree) has
  probability zero and defines no conditional distribution, so it is an
  error rather than a pseudoinverse guess.

- The conditional covariance is symmetrized; eigenvalues in
  `[-tolerance * lambda_max, 0)` are projected to zero and anything more
  negative is an error. It is factored by Cholesky when positive
  definite (platform-deterministic) and by `V sqrt(D)` otherwise. A
  coordinate with zero conditional variance is returned at its
  conditional mean exactly.

- **RNG.** Every check runs before the first random number is drawn, so
  a rejected call leaves `.Random.seed` untouched. A successful call
  consumes exactly `n * length(sample_coordinates)` standard normals
  from [`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html), in
  entity order and `sample_coordinates` order within an entity, whatever
  the observed patterns and even when some conditional variances are
  zero. Zero entities or zero sample coordinates consume nothing.

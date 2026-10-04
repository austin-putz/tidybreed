# Write one genetic covariance block to `trait_var_comp`

The single write path for generation targets, shared by
[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md)
and the generators' `G =`. It

- validates the block as finite, symmetric and positive semidefinite, on
  the correlation scale so the check does not depend on the traits'
  units;

- refuses when **any** row exists for `effect_name`, any of the block's
  traits and the same `line_name` (NULL-safe), checked on the whole
  table, even for an identical matrix;

- inserts all n^2 rows with `%.17g` literals, which round-trip a double
  exactly, through `dbExecute()` (never `dbWriteTable()`, which advances
  the RNG).

## Usage

``` r
.tvc_write_block(conn, effect_name, G, line_name = NULL)
```

## Arguments

- G:

  Named square matrix (`dimnames` = trait names).

## Details

It opens **no** transaction: the caller owns one, so a generator commits
the target together with its terms.

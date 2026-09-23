# Load the covariance block(s) touching a set of phenotypes

A *block* is a connected component of the graph whose edges are stored
pair rows in `phenotype_var_comp` for `effect_name` (a `cov_value` of
`0` is an edge). The requested phenotypes may fall into several
components; each is returned with every stratum stored for it. A
phenotype with no row at all belongs to no block and is simply absent
from the result — callers decide what that means (for the residual
effect it is "no residual variance").

## Usage

``` r
find_covariance_blocks(conn, effect_name, phenotype_names)
```

## Arguments

- conn:

  A DBI connection.

- effect_name:

  `"residual"` or a named random effect.

- phenotype_names:

  Character vector of the phenotypes whose blocks are wanted.

## Value

A list of blocks, ordered by their first member; possibly empty. Each
block is a list with:

- `effect_name`:

  as supplied.

- `phenotypes`:

  sorted character vector of members.

- `condition_table`, `condition_column`:

  the block's condition column, or `NULL` when it has only the
  unconditional stratum.

- `unconditional`:

  the unconditional matrix, or `NULL` if none is stored.

- `conditional`:

  a list of matrices named by `condition_level` (empty when there are no
  conditional strata).

Every matrix is `phenotypes` x `phenotypes` with dimnames, exactly
symmetric, in the sorted member order.

## Details

The writers in `R/phenotype_cov_block.R` guarantee that every stratum of
a block is a complete, symmetric matrix over the same phenotypes and
that a block has one `(condition_table, condition_column)`. Rows can
still be removed by hand with
[`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md),
so the loader re-checks those invariants and errors, naming the block,
when they no longer hold.

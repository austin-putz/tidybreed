# Members of the block(s) touching a set of phenotypes

The transitive closure over pair-row existence for `effect_name`, across
every stratum. Always contains `phenotype_names` itself, so a phenotype
with no stored rows is its own (empty) block.

## Usage

``` r
.pvc_block_members(conn, effect_name, phenotype_names)
```

## Arguments

- conn:

  A DBI connection.

- effect_name:

  Character scalar.

- phenotype_names:

  Character vector.

## Value

Character vector of block members, sorted.

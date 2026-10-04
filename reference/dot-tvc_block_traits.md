# The traits of the stored block(s) touching `traits`

A block is found from the rows themselves: the traits linked by
off-diagonal rows within one `effect_name` x `line_name`. No block id is
stored. Returns the connected closure of `traits`, restricted to traits
that have at least one stored row.

## Usage

``` r
.tvc_block_traits(conn, effect_name, line_name, traits)
```

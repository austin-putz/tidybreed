# The `line_name` whose rows a reader should use

`NULL` reads the population-wide rows. A named line reads its own rows
when it has any for `effect_name` and these traits, and otherwise falls
back to the population-wide rows. The fallback is decided per
`effect_name`.

## Usage

``` r
.tvc_resolve_line(conn, effect_name, traits, line_name)
```

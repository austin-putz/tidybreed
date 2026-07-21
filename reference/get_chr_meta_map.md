# Load `chr_meta` as a named list keyed by `chr_name`

Convenience reader used by
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
and
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
to avoid re-querying `chr_meta` per chromosome inside a loop.

## Usage

``` r
get_chr_meta_map(conn)
```

## Arguments

- conn:

  A DBI connection.

## Value

Named list (by `chr_name`), each element a 1-row data.frame with
`copy_mode_M`, `copy_mode_F`, `hemi_parent`, `recombines_M`,
`recombines_F`.

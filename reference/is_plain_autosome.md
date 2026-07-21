# Is a chromosome "plain" (both sexes full copy, recombines in both sexes)?

A chromosome is handled by the fast, unchanged autosome path in
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)/[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
iff both sexes get a full complement and recombination occurs in both
sexes — i.e. `chr_meta` has its default row for it. Any other
combination (sex-linked copy mode, or achiasmy in either sex) routes
through the special-chromosome path.

## Usage

``` r
is_plain_autosome(chr_row)
```

## Arguments

- chr_row:

  One row of `chr_meta` (as returned by
  [`get_chr_meta_map()`](https://austin-putz.github.io/tidybreed/reference/get_chr_meta_map.md)).

## Value

Logical scalar.

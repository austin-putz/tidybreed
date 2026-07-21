# Resolve a chromosome's copy count for a given sex and ploidy

Translates `chr_meta`'s relative `copy_mode`
(`"full"`/`"half"`/`"none"`) into an absolute copy count for one
individual, given their own ploidy. Shared by
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
and
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
so both agree on expected row counts per chromosome per sex —
duplicating this mapping inline in both files would risk them silently
disagreeing.

## Usage

``` r
resolve_chr_copy_count(copy_mode, ploidy = 2L)
```

## Arguments

- copy_mode:

  Character scalar: `"full"`, `"half"`, or `"none"`.

- ploidy:

  Integer scalar. The individual's own ploidy (always `2L` in this
  version of tidybreed).

## Value

Integer copy count (`0L`, `1L`, or `2L` in this version).

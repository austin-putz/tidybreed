# The traits with generated terms of one kind calibrated to one target scope

A `line_name = "C"` target covers line-C terms. A `line_name = NULL`
target covers the population-wide terms (common or parent-only) **and**
the line-scoped terms of every line that has no stored block of its own
for that kind and trait: the generator resolves a line's target with the
`line -> NULL` fallback
([`.tvc_resolve_line()`](https://austin-putz.github.io/tidybreed/reference/dot-tvc_resolve_line.md)),
so those terms were calibrated to the population-wide target. A line
block cannot be added under existing line terms (the refusal below), so
a line block present now was present when its terms were generated.

## Usage

``` r
.tvc_generated_traits(conn, effect_name, traits, line_name)
```

## Value

Sorted character vector of trait names (empty when none).

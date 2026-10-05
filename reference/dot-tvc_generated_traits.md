# The traits with generated terms of one kind at one target scope

A `line_name = NULL` target covers the population-wide terms (common or
parent-only); a `line_name = "C"` target covers line-C terms.

## Usage

``` r
.tvc_generated_traits(conn, effect_name, traits, line_name)
```

## Value

Sorted character vector of trait names (empty when none).

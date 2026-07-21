# Recursively extract all name nodes from a parsed R expression. Skips the call head (function name) and numeric/logical literals. Used by derived-formula validation and topo sort.

Recursively extract all name nodes from a parsed R expression. Skips the
call head (function name) and numeric/logical literals. Used by
derived-formula validation and topo sort.

## Usage

``` r
.extract_all_symbols(e)
```

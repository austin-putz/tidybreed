# Parse a memory size (bytes number or a string like "512MB", "2GB")

Parse a memory size (bytes number or a string like "512MB", "2GB")

## Usage

``` r
parse_mem_size(x)
```

## Arguments

- x:

  Numeric bytes, or a string with an optional k/m/g/t suffix
  (case-insensitive; a trailing "B" is optional). Powers of 1024.

## Value

Numeric bytes.

# Composite row keys per table (used for exact single-table DELETE)

For tables with a single PK this is length-1. For composite-PK tables
all key columns are listed so that
[`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md)
can target exact rows rather than over-deleting by `id_ind` alone.

## Usage

``` r
TABLE_ROW_KEYS
```

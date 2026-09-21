# Tables where row deletion is not a meaningful operation

A table missing from
[TABLE_ROW_KEYS](https://austin-putz.github.io/tidybreed/reference/TABLE_ROW_KEYS.md)
is ambiguous: it may be a deliberate refusal or it may be an oversight,
and from the outside those look identical. Listing a table here makes
the refusal a decision with a reason attached, which
[`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md)
reports instead of the generic "not registered" error.

## Usage

``` r
TABLE_NO_ROW_DELETE
```

## Details

Names are table names; values are the reason shown to the user.

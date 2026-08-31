# Display groups for [`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md), in pipeline order

Option A from `plans/update_schema_print.md`: a hard-coded vector rather
than a name-prefix rule or a `_schema_meta` column. A table added later
and not registered here degrades *visibly* — it appears under **User
tables** — rather than being silently misfiled by a lexical rule that
cannot tell `ind_haplotype` (raw genome data) from `ind_tbv` (simulation
output).

## Usage

``` r
.schema_table_order()
```

## Value

A named list of character vectors, in display order.

## Details

In-group order is workflow order, not alphabetical. `"User tables"` is
an intentionally empty slot: unrecognized tables land there, sorted by
name, and `"System"` always prints last.

Must name the same tables as the `.\*_descriptions()` helpers above.

# All pre-built `_schema_meta` descriptions, in display-group order

Registered once by
[`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md).
Description rows are plain metadata and do not require their table to
exist yet, so registering the whole set up front is cheaper than
re-registering a subset from every `define_*` entry point — and it
removes the window in which a later `define_*` call would clobber a
user's
[`define_schema_description()`](https://austin-putz.github.io/tidybreed/reference/define_schema_description.md)
override.

## Usage

``` r
.all_schema_descriptions()
```

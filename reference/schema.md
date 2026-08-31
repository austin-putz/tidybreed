# View table-level descriptions for a tidybreed population

Returns a tibble of all user-visible tables with row counts, column
counts, and descriptions from `_schema_meta`. Print the result for an
aligned overview; use `describe_table(pop, "name")` to drill into a
specific table.

Tables are ordered by pipeline stage — the order in which a population
is built — and printed under group headings. The order does not depend
on how the population was opened:
[`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md)
populates `pop$tables` in creation order and
[`restore_pop()`](https://austin-putz.github.io/tidybreed/reference/restore_pop.md)
in DuckDB catalog order, so before this the same `.duckdb` file printed
differently depending on which one you used.

## Usage

``` r
schema(
  pop,
  order = c("pipeline", "name", "rows", "size"),
  show_empty = FALSE,
  include_system = FALSE,
  sizes = FALSE
)
```

## Arguments

- pop:

  A `tidybreed_pop` object, or a `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (its `pop` reference is used).

- order:

  One of `"pipeline"` (default, grouped by build stage with headings),
  `"name"` (flat alphabetical), `"rows"` (flat, largest row count first)
  or `"size"` (flat, largest on-disk size first). Group headings print
  only for `"pipeline"`; `table_group` stays in the returned tibble for
  all four. `"size"` requires `sizes = TRUE` and errors otherwise.

- show_empty:

  Logical. `FALSE` (default) collapses each group's zero-row tables into
  a single `+ n empty: ...` line — on a freshly built population most
  tables are empty and would otherwise bury the ones with data. `TRUE`
  prints every table on its own row.

- include_system:

  Logical. `FALSE` (default) hides the `_schema_meta` system table;
  `TRUE` lists it under a **System** group.

- sizes:

  Logical. `FALSE` (default). `TRUE` adds a per-table `size_bytes`
  column and prints a `Size` column with a mandatory caveat footnote.
  Opt-in because collecting it requires a `CHECKPOINT`, which writes to
  the database — see **Size reporting**.

## Value

A tibble of class `tidybreed_schema` with columns `table_name`,
`table_group`, `n_rows`, `n_cols`, and `description`. Printed via
[`print.tidybreed_schema()`](https://austin-putz.github.io/tidybreed/reference/print.tidybreed_schema.md).
`table_group` is a factor whose level order is the display order, so the
grouping is available as data for
[`filter()`](https://rdrr.io/r/stats/filter.html) /
[`split()`](https://rdrr.io/r/base/split.html) and not only as printed
text.

## Size reporting

`sizes = TRUE` adds a per-table `Size` column. Three things about it are
load-bearing:

1.  **It writes.** DuckDB has no direct per-table byte count; the
    workable route is `PRAGMA storage_info()`, and before a `CHECKPOINT`
    recently written data carries `block_id = -1` so every table reports
    zero. `schema()` therefore issues a `CHECKPOINT` when
    `sizes = TRUE`, and only then. This is the sole reason the argument
    is opt-in rather than always on.

2.  **It is quantized.** Size is on-disk storage attributed in whole 256
    KiB blocks, so small tables are not distinguishable: most read as
    0.25 MiB, and a table small enough to live in the catalog reads as
    0.00 MiB while still holding rows.

3.  **It does not add up.** Per-table sizes sum to less than the
    database total in the header, because catalog and header blocks are
    not attributed to any table. Blocks are not shared between tables,
    so the attribution is sound as far as it goes.

Caveats 2 and 3 are printed as a footnote under the table whenever the
column appears; the footnote has no suppression argument. For an
in-memory population there are no blocks at all, so `sizes = TRUE` warns
and omits the column.

## Database size

The printed header reports the size of the whole database, read from
`PRAGMA database_size`. For a file-backed population this is the
`.duckdb` file size, plus the write-ahead log when it is non-empty (the
WAL is only folded into the file by a `CHECKPOINT`, which `schema()`
deliberately does not issue — it is a write). For an in-memory
population there is no file, so the header reports memory usage instead
and says so.

## See also

[`describe_table()`](https://austin-putz.github.io/tidybreed/reference/describe_table.md),
[`define_schema_description()`](https://austin-putz.github.io/tidybreed/reference/define_schema_description.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- open_pop(pop_name = "MySim", db_name = ":memory:") |>
  define_genome(n_loci = 500, n_chr = 5, chr_len_Mb = 100)
schema(pop)

# Expand the collapsed empty tables, and show the system table
schema(pop, show_empty = TRUE, include_system = TRUE)

# The grouping is data, not just print formatting
subset(schema(pop), table_group == "Genome")

# Flat orderings answer "what is actually big here?"
schema(pop, order = "rows")

# On-disk sizes (issues a CHECKPOINT; file-backed populations only)
schema(pop, order = "size", sizes = TRUE)
} # }
```

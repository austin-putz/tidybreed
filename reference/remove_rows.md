# Remove rows from one or more population tables

Removes rows from the database based on a filtered `tidybreed_table`.
Chain after
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
and
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
to target specific rows. Returns the population invisibly — consistent
with all other action functions.

## Usage

``` r
remove_rows(
  tbl,
  tables = NULL,
  confirm_all = FALSE,
  dry_run = FALSE,
  verbose = TRUE
)
```

## Arguments

- tbl:

  A `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md),
  with a
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
  applied.

- tables:

  `NULL` (default) to delete from the current table only; `"all"` to
  delete from every `ind_*` table (including `ind_haplotype` and
  `ind_genotype`) that exists in the population; or a character vector
  of specific table names (each must have an `id_ind` column). When not
  `NULL`, the source table must also have an `id_ind` column.

- confirm_all:

  Logical. Set to `TRUE` to allow deletion of all rows when no filter
  has been applied. This is a deliberate safeguard — you must explicitly
  opt in to wiping an entire table. Default `FALSE`.

- dry_run:

  Logical. If `TRUE`, report what would be deleted without modifying the
  database. Default `FALSE`.

- verbose:

  Logical. If `TRUE`, print a message after each deletion reporting the
  row count removed. Default `TRUE`.

## Value

The `tidybreed_pop` object (invisibly).

## Details

A filter is **required** by default. Calling `remove_rows()` without a
prior
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
will stop with an error explaining how to opt in to full-table deletion
via `confirm_all = TRUE`.

**Single-table mode** (`tables = NULL`): uses the composite row key from
the internal `TABLE_ROW_KEYS` registry to delete exactly the rows
matched by the filter — no more, no less. For example, filtering
`ind_tbv` by `trait_name == "ADG"` deletes only the ADG rows, not all
TBV rows for those animals.

**Cross-table mode** (`tables != NULL`): extracts unique `id_ind` values
from the filtered table and issues a `DELETE ... WHERE id_ind IN (...)`
for each target table via a temp-table JOIN. `tables = "all"` targets
every `ind_*` table (including `ind_haplotype`, `ind_genotype`, and
`ind_true_index`) that currently exists in the population (including
`ind_meta`).

A temporary DuckDB table is used for both modes to avoid SQL injection
risk and to handle large ID sets efficiently.

## Examples

``` r
if (FALSE) { # \dontrun{
# Delete specific phenotype records for culled animals
pop |>
  get_table("ind_phenotype") |>
  dplyr::filter(id_ind %in% culled_ids, phenotype_name == "litter_size") |>
  remove_rows()

# Remove all data for culled animals across every individual table
pop |>
  get_table("ind_meta") |>
  dplyr::filter(id_ind %in% culled_ids) |>
  remove_rows(tables = "all")

# Delete from an explicit subset of tables only
pop |>
  get_table("ind_meta") |>
  dplyr::filter(id_ind %in% culled_ids) |>
  remove_rows(tables = c("ind_phenotype", "ind_tbv"))

# Preview before deleting
pop |>
  get_table("ind_meta") |>
  dplyr::filter(id_ind %in% culled_ids) |>
  remove_rows(tables = "all", dry_run = TRUE)

# Wipe an entire table (requires confirm_all = TRUE)
pop |>
  get_table("ind_phenotype") |>
  remove_rows(confirm_all = TRUE)
} # }
```

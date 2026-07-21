# Summarize a tidybreed population

Produces a detailed summary of all tables in a `tidybreed_pop` object.
For each table, shows row count, column count, and per-column
statistics: frequency tables for low-cardinality categorical columns,
5-number summaries for numeric columns, and date ranges for
date/timestamp columns. Any wide table with `locus_<n>` columns only
summarizes its non-locus metadata columns to avoid querying
50,000-column tables.

## Usage

``` r
# S3 method for class 'tidybreed_pop'
summary(object, tables = NULL, max_values = 12L, ...)
```

## Arguments

- object:

  A `tidybreed_pop` object.

- tables:

  Character vector of table names to include. `NULL` (default) includes
  all tables registered in `object$tables`.

- max_values:

  Integer. Columns with at most this many distinct values receive a
  frequency table display; default `12L`.

- ...:

  Additional arguments (not used; required by the S3 generic).

## Value

A `tidybreed_summary` object (a named list).

## See also

[`print.tidybreed_pop()`](https://austin-putz.github.io/tidybreed/reference/print.tidybreed_pop.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pop <- open_pop(pop_name = "MySim", db_name = ":memory:") |>
  define_genome(n_loci = 1000, n_chr = 5, chr_len_Mb = 100) |>
  define_founder_haplotypes(n_haplotypes = 50)
pop <- pop |>
  get_table("founder_haplotypes") |>
  add_founders(n_males = 50, n_females = 50, line_name = "A", gen = 0L)
summary(pop)
summary(pop, tables = c("ind_meta", "genome_meta"))
summary(pop, max_values = 5L)
} # }
```

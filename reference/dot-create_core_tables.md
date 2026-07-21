# Create all core (non-genome) database tables

Internal helper called by
[`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md).
Creates the `_schema_meta` system table and all non-genome metadata
tables. Genome tables (`genome_meta`, `ind_haplotype`, `ind_genotype`,
`chr_meta`) are created by
[`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md).

## Usage

``` r
.create_core_tables(db_conn)
```

## Arguments

- db_conn:

  An active DuckDB connection.

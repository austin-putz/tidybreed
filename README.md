# tidybreed

<!-- badges: start -->
<!-- badges: end -->

A pipe-friendly R package for breeding program simulation backed by DuckDB. Design large-scale genomic simulations without running out of memory — all data lives on disk in a DuckDB database and is queried lazily via `dplyr`.

## Installation

Install [pak](https://pak.r-lib.org/) from the Posit Package Manager, then install tidybreed from GitHub:

This is the new `repos` that is recommended by Posit. 

```r
# install pak package
install.packages("pak", repos = "https://packagemanager.posit.co/cran/latest")

# laod pak
library(pak)

# install this package
pak::pak("austin-putz/tidybreed")

# load tidybreed
library(tidybreed)
```

## Core Concept

Every function in tidybreed accepts and returns a `tidybreed_pop` object — a thin wrapper around a DuckDB connection. Chain operations with `|>` (or `%>%`):

```r
pop <- initialize_genome(...) |>
  add_founders(...)           |>
  mutate_ind_meta(...)        |>
  define_chip(...)
```

All tables are queryable at any point with `get_table()` and standard `dplyr` verbs.

Users should be able to add to tables at any time or mutate them. This allows
immense customizable data for the user and then query based on custom fields for
instance and go between filtering the individual data and the genome data. 

## Examples

### 1. Initialize a genome

```r
pop <- initialize_genome(
  pop_name   = "PigLine_A",
  n_loci     = 50000,
  n_chr      = 18,
  chr_len_Mb = 100
)

print(pop)
# <tidybreed population>
# Population: PigLine_A
# Loci: 50,000  Chromosomes: 18
# Tables: genome_meta, genome_haplotype, genome_genotype
```

> Remember: If you ever forget the table names, just print the `pop` object. 

### 2. Add founder animals

```r
pop <- pop |>
  add_founders(
    n_males   = 20,
    n_females = 200,
    line_name = "A"
  )

# Inspect individual data
pop %>% get_table("ind_meta)

# or return a tibble() with 
pop %>% get_table("ind_meta") |>
  dplyr::collect()   # collect will pull into R memory and create a tibble
```

Users can then summarize or manipulate the tibble() however they desire. 

This means they can pull genotypes, haplotypes and say QTL effects
and calculate it all on their own which means possibilities are endless. 

### 3. Add custom individual metadata

The real secret sauce of this package is to use user
defined data and fields into the database so you no longer
have to be responsible for this yourself and saving to disk saves
all the information for retrieval later. 

```r
pop <- pop |>
  mutate_ind_meta(
    gen        = 0L,       # add 0 as integer
    farm       = "Iowa",   # add farm name to all individuals
    date_birth = as.Date("2024-01-01")  # add your own date, now possible to run "real life" simulations without hacks
  )
```

Types are inferred automatically. Scalar values are broadcast to all current individuals; pass a vector to assign per-individual values.

### 4. Define a SNP chip

```r
# Evenly spaced 50K chip
pop <- pop |>
  define_chip(chip_name = "50K", n = 50000, method = "even")

# Add a denser HD chip proportional to chromosome length
pop <- pop |>
  define_chip(chip_name = "HD", n = 80000, method = "chr_even")
```

By default, this adds 'is_<chip_name>' to the 'genome_meta' table in the database. 

This provides the user the ability to filter the genotypes by QTL, SNP chip, or whatever they desire. 

```r
# Check chip assignment on genome_meta
get_table(pop, "genome_meta") |>
  dplyr::filter(is_50K == TRUE) |>
  dplyr::select(locus_name, chr, pos_Mb) |>
  dplyr::collect()
```

This filters the 'genome_meta' table to only the 50k chip loci and then selects af few columns
and then finally returns a `tibble()` (in memory) in R. 

### 5. Query with dplyr

All tables are lazy — filter and select before pulling data into R:

```r
# Males from generation 0
get_table(pop, "ind_meta") |>
  dplyr::filter(sex == "M", gen == 0L) |>
  dplyr::collect()

# Chromosome 1 loci on the 50K chip
get_table(pop, "genome_meta") |>
  dplyr::filter(chr == 1, is_50K == TRUE) |>
  dplyr::arrange(pos_Mb) |>
  dplyr::collect()
```

### 6. Coming soon: traits and phenotypes

```r
# Define trait architecture (QTL + residual variance)
pop <- pop |>
  add_trait(
    name           = "ADG",
    n_qtl          = 20,
    target_add_var = 1.0,
    res_var        = 0.5
  )

# Simulate phenotypes for selected animals
pop <- pop |>
  add_phenotype(trait = "ADG")
```

## Function Overview

| Function              | Purpose                                              |
|-----------------------|------------------------------------------------------|
| `initialize_genome()` | Create DuckDB database with genome tables            |
| `add_founders()`      | Add founder individuals with haplotypes/genotypes    |
| `mutate_ind_meta()`   | Add or update columns in `ind_meta`                  |
| `mutate_genome_meta()`| Add or update columns in `genome_meta`               |
| `define_chip()`       | Mark loci as on a named SNP chip                     |
| `get_table()`         | Return a lazy `dplyr` tibble from any table          |
| `close_pop()`         | Close the DuckDB connection                          |

## License

MIT

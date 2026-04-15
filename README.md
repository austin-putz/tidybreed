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
  pop_name   = "Pigs",         # name population (name individual genetic lines below)
  n_loci     = 100,            # total n loci to simulate
  n_chr      = 18,             # total n chromosomes to simulate
  chr_len_Mb = 100,            # length of each chromosome
  n_haplotypes = 200,          # number of haplotypes generated to sample from in founders
  overwrite    = TRUE          # overwrite the current duckdb (database) if it exists
)

print(pop)
# <tidybreed population>
# Population: PigLine_A
# Loci: 50,000  Chromosomes: 18
# Tables: genome_meta, genome_haplotype, genome_genotype
```

> Remember: If you ever forget the table names, just print the `pop` object. 

### 2. Add founder animals

We start by adding founder individuals. 

```r
pop <- pop |>
  add_founders(
    n_males   = 2,     # 2 male founders (sampled from haplotypes generated above)
    n_females = 8,     # 8 female founders
    line_name = "A"    # name of genetic line if you want to crossbred later
  )
```

This function adds 10 total individuals as founders and calls them line "A". 

Let's peak at that new table created. 

```r
# Inspect individual data
pop |> 
  get_table("ind_meta")
```

You can see we have 5 default columns that are the only ones required to simulate. 

```r
# or return a tibble() with the collect() function
pop |> 
  get_table("ind_meta") |>
  dplyr::collect()   # collect will pull into R memory and create a tibble
```

Users can then summarize or manipulate the `tibble()` however they desire. 

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
    date_birth = as.Date("2026-01-01")  # add your own date, now possible to run "real life" simulations without hacks
  )
```

Types are inferred automatically. Scalar values are broadcast to all current individuals; pass a vector to assign per-individual values.

### 4. Define a SNP chip

```r
# use mutate_genome_meta() to create your own custom chip 
pop <- pop %>%
  mutate_genome_meta(
    is_CustomSNPChip = sample(c(rep(TRUE,50), rep(FALSE,50)), size=100, replace=FALSE)    # generates exactly 50 SNP for this panel
  )

# Evenly spaced low-density chip
pop <- pop |>
  define_chip(name = "LowDensity", n = 30, method = "even")

# Denser chip proportional to chromosome length
pop <- pop |>
  define_chip(name = "HighDensity", n = 90, method = "random")

# Logical vector — define a chip as the complement of an existing one
ld_tf <- get_table(pop, "genome_meta") |> dplyr::pull(is_LowDensity)
pop <- pop |>
  define_chip(name = "NonLD", locus_tf = !ld_tf)
```

By default, this adds `is_<chip_name>` to the `genome_meta` table in the database.

The `locus_tf` argument accepts any logical vector of the same length as the number of loci,
making it easy to compose chips from existing membership columns or custom logic.

This provides the user the ability to filter the genotypes by QTL, SNP chip, or whatever they desire. 

```r
# Check chip assignment on genome_meta
get_table(pop, "genome_meta") |>
  dplyr::filter(is_LowDensity == TRUE) |>
  dplyr::select(locus_name, chr, pos_Mb) |>
  dplyr::collect()
```

This filters the 'genome_meta' table to only the 50k chip loci and then selects af few columns
and then finally returns a `tibble()` (in memory) in R. 

### 5. Query with dplyr

All tables are lazy — filter and select before pulling data into R:

```r
# Males from generation 0
pop %>%
  get_table("ind_meta") %>%
  dplyr::filter(sex == "M", gen == 0L) %>%
  dplyr::collect()

# Chromosome 1 loci on the 50K chip
get_table(pop, "genome_meta") |>
  dplyr::filter(chr == 1, is_HighDensity == TRUE) |>
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

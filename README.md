# `tidybreed` R Package

<!-- badges: start -->
<!-- badges: end -->

A pipe-friendly (`%>%` or `|>`) R package for breeding program simulation backed by [DuckDB](https://duckdb.org). 
Design large-scale genomic simulations without running out of memory — all data lives on disk in a DuckDB database 
and is queried lazily via `dplyr`. 

## Design

| Main | :arrow_right: | Reason or Description |
| ---- | ------------- | --------------------- |
| [R](https://www.r-project.org) | :arrow_right: | R is a standard for most scientists, allows easy design of custom (flexible) breeding programs
| **Modern Database**  | :arrow_right: | [DuckDB](https://duckdb.org) for efficient storage on disk, unlimited manipulation for users custom programs, nothing can store data like a database
| Pipe everything      | :arrow_right: | Animals, plants, insects, etc all need 'selected' for phenotyping, genotyping, evaluations, mating, selection, etc. Use `tidyverse` language to add what you want for each individual you need
| Flexible Base        | :arrow_right: | Users should be able to add everything custom for their breeding program using a `mutate_table()` function, then I can provide "helper functions" later on top (should allow maximum features before the author has time to add a helper function)

## Installation

Install [pak](https://pak.r-lib.org/) from the Posit Package Manager, then install `tidybreed` from [GitHub](https://github.com/austin-putz/tidybreed/):

This is the new `repos` that is recommended by Posit. 

```r
# install pak package
install.packages("pak", repos = "https://packagemanager.posit.co/cran/latest")

# laod pak
library(pak)

# install the tidybreed package
pak::pak("austin-putz/tidybreed")

# load tidybreed
library(tidybreed)
```

Pay attention to the `version` of the package as it may change rapidly pre-version `1.0.0`, you were warned... 

> DISCLAIMER: Pre `v1.0.0` packages are considered in 'beta' testing and subject to make many breaking changes

## Core Concept

Every function in `tidybreed` accepts and returns a `tidybreed_pop` object — a thin wrapper around a 
DuckDB connection. Chain operations with `|>` (or `%>%`):

```r
pop <- pop(...) |>        # pop contains the DuckDB connection
  get_table(...) |>       # identify which DB table you want to work with
  filter(...) |>          # filter out the individuals you want to work on or add to
  add_phenotype(...)      # use `add_*()` function to add rows to a table (often different table)
```

All tables are queryable at any point with `get_table()` and standard `dplyr` verbs, however
you need to use `collect()` before pulling into R memory and returning a `tibble()`. 

Users should be able to add to tables at any time or mutate them. This allows
**immense customizable data** for the user and then query based on custom fields for
instance and go between filtering the individual data and the genome data. No
more endless hack after hack to identify the individuals you want. 

Add phenotypes on 'selected' individuals

```r
pop %>%                               # pop points to the DuckDB connection
  get_table("ind_meta") %>%           # select your table to 'work on/with'
  filter(sex == "M", gen == 1L) %>%   # select what animals you want phenotyped for this trait
  add_phenotype("ADG")                # add phenotypes based on "trait_name" here "Average Daily Gain"
```

This allows users to easily add **phenotypes**, **genotypes**, **TBV**, **EBV**, **Index**, etc 
to any custom list of individuals. Making it extremely flexible and powerful. 

Furthermore, users can double query to select animals for mating using this pattern:

**Step 1:** Pull list of IDs for 'candidates'

```r
male_candidates <- pop %>%             # name list of IDs pulled below
  get_table("ind_meta") %>%            # pull from 'ind_meta' table (1 row / ind)
  filter(sex == "M", gen == 1L) %>%    # filter by sex + generation number
  pull(id_ind)                         # pulls a character vector of IDs for these individuals
```

**Step 2:** Choose what you want to select based on and filter IDs based on this list of candidates

```r
males_selected <- pop %>%          # now select males based on any table you want
  get_table("ind_ebv") %>%         # select any table to select animals based on (here EBV)
  filter(
    id_ind %in% male_candidates,   # subset EBV table by IDs pulled above in step 1
	trait_name == "ADG"            # filter by trait you are selecting for
  ) %>%
  slice_max(ebv, n=5) %>%          # select top 5 for ADG EBV
  pull(id_ind)                     # pull list of IDs into vector
```

Do not run this code as we didn't add column `gen` to `ind_meta` table or run `add_ebv()` to 
calculate EBVs yet...

Yet, this gives you a deep understand of the eventual power of `tidybreed` 

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

The real **secret sauce** of this package is to use user
defined data and fields into the database so you no longer
have to be responsible for this yourself and saving to disk saves
all the information for retrieval later. 

```r
pop <- pop |>
  mutate_table(
    gen        = 0L,       # add generation 0 as integer
    farm       = "Iowa",   # add farm name to all (character in R, VARCHAR in SQL)
    date_birth = as.Date("2026-01-01")  # add your own dates, now possible to run "real life" simulations without hacks
  )
```

Types are inferred automatically. Scalar values are broadcast to all current individuals; pass a vector to assign per-individual values.

### 4. Define a SNP chip

Use `mutate_table()` again to define a SNP chip:

```r
# use mutate_genome_meta() to create your own custom chip 
pop <- pop %>%
  mutate_table(
    is_50k = sample(c(rep(TRUE,50), rep(FALSE,50)), size=100, replace=FALSE)    # generates exactly 50 SNP for this panel
  )
```

Or use a 'helper-function' called `define_chip()`

```r
# Evenly spaced low-density chip
pop <- pop |>
  define_chip(chip_name = "LowDensity", n = 30, method = "even")
```

```r
# Denser chip proportional to chromosome length
pop <- pop |>
  define_chip(chip_name = "HighDensity", n = 90, method = "random")
```

We can extract the loci used in our 'LowDensity' SNP Chip to reuse. 

```r
# Logical vector — define a chip as the complement of an existing one
ld_tf <- pop |> 
  get_table("genome_meta") |> 
  dplyr::pull(is_LowDensity)

# define new chip called 'NonLD'
pop <- pop |>
  define_chip(chip_name = "NonLD", locus_tf = !ld_tf)
```

This will then identify all loci NOT on the 'LowDensity' chip as the 'NonLD' SNP chip. 

By default, this adds `is_<chip_name>` to the `genome_meta` table in the database.

The `locus_tf` argument accepts any logical vector of the same length as the number of loci,
making it easy to compose chips from existing membership columns or custom logic.

This provides the user the ability to filter the genotypes by QTL, SNP chip, or whatever they desire. 
No longer the need to have the software package do this for the user. 

```r
# Check chip assignment on genome_meta
pop %>%
  get_table("genome_meta") |>                  # table - 1 row per loci
  dplyr::filter(is_LowDensity == TRUE) |>      # filter by LowDensity SNP Chip loci
  dplyr::select(locus_name, chr, pos_Mb) |>    # pull only these 3 columns/fields from table
  dplyr::collect()                             # pull into R memory for further whatever...
```

This filters the `genome_meta` table to only the 50k chip loci and then selects af few columns
and then finally returns a `tibble()` (in memory) in R. 

### 5. Query with dplyr

All tables are lazy — filter and select before pulling data into R:

```r
# Males from generation 0
pop %>%
  get_table("ind_meta") %>%                   # select your table (ind_meta is 1 row per ind)
  dplyr::filter(sex == "M", gen == 0L) %>%    # select which animals
  dplyr::collect()                            # pull all data from 'ind_meta' for males in gen 0

# Chromosome 1 loci on the 50K chip
pop %>%
  get_table("genome_meta") |>                           # has 1 row per loci
  dplyr::filter(chr == 1, is_HighDensity == TRUE) |>    # only take chr 1 loci and on the HD chip
  dplyr::arrange(pos_Mb) |>                             # sort by position
  dplyr::collect()                                      # pull into R memory
```

### 6. Add Trait 

#### Add General Trait Info

Add a trait called **ADG** (Average Daily Gain):

We can first add basic trait information, can be used later by other functions for
**phenotyping**, adding QTL loci and effects, etc. 

```r
pop <- pop %>%
  add_trait(
    trait_name = "ADG",
	description = "Average Daily Gain",
	units = "g/d",
	trait_type = "continous",
	repeatable = FALSE,
	target_add_mean = 0,
	target_add_var = 1000,
	residual_var = 2000,
	index_weight = 0.05,
	economic_value = 0.02,
	overwrite = TRUE
  )
```

#### Add to Genome Info

Now we need **QTL loci and effects**:

```r
# Define trait architecture (QTL + residual variance)
pop <- pop |>
  # set QTL loci in genome for ADG
  define_qtl(
    trait_name = "ADG",
	n = 10,                  # 10 QTL loci
	method = "random"        # randomly select from genome_meta table
  ) |>
  # set QTL effects in genome for ADG
  set_qtl_effects(
    trait_name = "ADG", 
	distribution = "normal",   # effects sampled from normal distribution
	scale_to_target = TRUE,    # scale to target additive var (set above)
	base = "current_pop"       # use current population (founders added above, TBV = 0)
  )
```

Or simply use `mutate_table()` and do it all yourself! (encouraged)

Name new fields `is_QTL_<trait name>` and `add_<triat name>` in `genome_meta`

#### Add TBV

You can calculate and store the True Breeding Values (TBV) with
`add_tbv()`

First, we can add **custom user defined columns** in the `ind_tbv` table. 

```r
pop %>%
  get_table("ind_tbv") %>%
  mutate_table(
    type = NA_character_
  )
```

Make sure to understand `NA_character_` initilizes a `character` (R) or `VARCHAR` (SQL) variable
with nothing inside the field yet (no rows). 

Now we can calculate the TBV for all animals + add custom field value for these individuals:

```r
# calculate and store TBV for ADG trait
pop <- pop %>%
  get_table("ind_meta") %>% # here we specify the 'ind_meta' table so all animals will have their TBV calculated
  add_tbv(
    "ADG", 
    type = "additive"       # USER DEFINED Above
  )
```

Now view the table yourself, pull it with `collect()` and analyze with `dplyr` if 
desired:

```r
# now look at updated 'ind_tbv' table
pop %>% get_table("ind_tbv") %>% print()
```

#### Add Effects for Phenotyping 

**Add Overall Mean**

Use `add_effect_int()` for overall mean of phenotype: 

```r
pop %>%
  add_effect_int(
    "ADG",
	mean = 1000
  )
```

**Add fixed effect** with `add_effect_fixed_class()`:

```r
# add sex effect for ADG
pop %>%
  add_effect_fixed_class(
    trait_name = "ADG",
    effect_name = "sex",
    source_column = "sex", 
    levels = c(M=200, F=0),
    source_table = "ind_meta",
    overwrite = TRUE
  )
```

#### Add Phenotypes

Just like with `ind_tbv` we can first add a **custom field** to
`ind_phenotype` table:

```r
pop %>%
  get_table("ind_phenotype") %>%
  mutate_table(
    environment = NA_character_    # new column in 'ind_phenotype' table
  )
```

Now we can add phenotypes to all individuals:

```r
# phenotype all individuals for ADG
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
  add_phenotype(                # add rows to 'ind_phenotye' table
    trait = "ADG",              # trait name
    environment = "hot"         # USER DEFINED FIELD Above
  )  
```

**OR we can select which animals to phenotype!!**

```r
# phenotype ONLY males
pop %>%
  get_table("ind_meta") %>%     # will phenotype all individuals in this table with no filter
  filter(sex == "M") %>%        # select only males to phenotype for trait
  add_phenotype(                # add rows to 'ind_phenotye' table
    trait = "ADG",              # trait name
    environment = "hot"         # USER DEFINED FIELD Above
  )  
```

## 7. Run Evaluations (EBVs)

The package wouldn't be complete without some way to run evaluations to calculate
**EBVs** (Estimated Breeding Values):

Again, add a **custom field** to the new `ind_ebv` table:

```r
pop %>%
  get_table("ind_ebv") %>%
  mutate_table(
    option_set = NA_character_    # USER DEFINED - for parameters used maybe?
  )
```

Run [BLUPF90](https://nce.ads.uga.edu/wiki/doku.php?id=start) to calculate EBVs: 

```r
pop %>%
  get_table("ind_meta") %>%      # for all animals, no filters
  add_ebv(
    trait    = "ADG",
    software = "blupf90",        # Software package to use
    model    = "blup",           # BLUP model (no genomics)
    option_set = "2-gen-back"    # USER DEFINED FIELD
  )
```

Print `ind_ebv` table: 

```r
pop %>%
  get_table("ind_ebv") %>%
  print(n=50)
```

## Function Overview

| Function              | Purpose                                              |
|-----------------------|------------------------------------------------------|
| `initialize_genome()` | Create DuckDB database with genome tables            |
| `get_table()`         | Return a lazy `dplyr` tibble from any table          |
| `mutate_table()`      | Add or update columns in any table                   |
| `add_founders()`      | Add founder individuals with haplotypes/genotypes    |
| `define_chip()`       | Mark loci as on a named SNP chip                     |
| `define_qtl()`        | set QTL loci by trait                                | 
| `set_qtl_effects()`   | set QTL effects by trait                             | 
| `add_phenotype()`     | add phenotype by summing QTL + fixed/random + resid  | 
| `add_tbv()`           | sum QTL effects for True Breeding Value              | 
| `add_ebv()`           | run evaluation software or parent average            |
| `add_offspring()`     | add progeny/offspring based on tibble mating design  | 
| `close_pop()`         | Close the DuckDB connection                          |

## License

MIT (for now)




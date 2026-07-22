# Introduction to tidybreed

## What `tidybreed` is

`tidybreed` simulates animal and plant breeding programs. You build a
population step by step — a genome, founder animals, traits, phenotypes,
matings, selection — and the package tracks every individual, haplotype,
and record for you.

What makes it different from other simulators is where the data lives.
**Everything is stored in a [DuckDB](https://duckdb.org/) database, not
in R memory.** A simulation with a million animals does not need a
million animals’ worth of RAM. The database is a plain file on disk, so
runs are resumable, shareable, and queryable with ordinary SQL long
after the R session that made them has ended.

The API is pipe-friendly and built on `dplyr` verbs, so if you know the
tidyverse you already know most of the syntax.

This vignette walks through a complete single-generation simulation.
Every code block below actually runs when this page is built, so what
you see is real output.

## Installation

`tidybreed` is not on CRAN yet. Install from GitHub — `pak` handles the
compiled code and dependencies cleanly:

``` r

install.packages("pak")
pak::pak("austin-putz/tidybreed")
```

You will need a C++ compiler (Xcode command line tools on macOS, Rtools
on Windows). Then:

``` r

library(tidybreed)
library(dplyr)

set.seed(1)   # so this page renders the same numbers every time
```

## How the package works

Three ideas explain almost all of the API. They are worth reading before
the code.

### 1. R is the control layer, DuckDB is the storage layer

Your R session holds a small object with a database connection in it. It
does **not** hold your genotypes. When you ask for a table you get a
*lazy reference* — nothing is pulled into R until you call
[`collect()`](https://dplyr.tidyverse.org/reference/compute.html). This
is why you should filter first and collect last.

### 2. There are two object types

| Object | Created by | What it is |
|----|----|----|
| `tidybreed_pop` | [`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md) | The population: a database connection plus a table registry |
| `tidybreed_table` | [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md) | A lazy reference to one table, plus any pending filter |

Most functions take a `tidybreed_pop` and return a `tidybreed_pop`, so
they chain:

``` r

pop <- pop |> define_genome(...) |> define_founder_haplotypes(...)
```

But **action functions** — the ones that operate on a *subset of rows* —
take a `tidybreed_table` as their first argument instead. That is how
you tell them which rows to act on:

``` r

pop <- pop |>
  get_table("ind_meta") |>       # which table
  filter(sex == "F") |>          # which rows
  add_phenotype("ADG")           # what to do
```

[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md),
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md),
[`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md),
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md),
[`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md),
[`add_genotypes()`](https://austin-putz.github.io/tidybreed/reference/add_genotypes.md),
[`extract_genotypes()`](https://austin-putz.github.io/tidybreed/reference/extract_genotypes.md),
and
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
all work this way. Most of them return a `tidybreed_pop`, so the chain
continues; the `extract_*` family is the exception — it hands you data
back as a tibble, since pulling results into R is the whole point of it.

There is no positional row selection anywhere in `tidybreed` — no
`TRUE`/`FALSE` masks, no integer indices. Database row order is not
guaranteed, so selection always goes through names (`id_ind`,
`locus_name`) or a SQL filter.

### 3. Genetics and observations are separate layers

This one trips people up, so it is worth being explicit. A trait’s
**genetics** and the **phenotype someone writes on a clipboard** are
defined by two different functions:

| You want to specify | Function | Examples |
|----|----|----|
| Genetic architecture | [`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md) | `target_add_var`, `target_add_mean`, `expressed_parent`, `units` |
| What gets observed | [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md) | `mean`, `type`, `residual_var`, `expressed_sex`, `repeatable` |

For a simple trait the two share a name and you call both. The payoff
comes later: one observed phenotype (weaning weight) can be built from
several genetic components (direct + maternal), and the split makes that
natural instead of a special case.

## Open a population

[`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md)
creates the database. Here we use `":memory:"` so this vignette writes
nothing to disk; in a real run you would omit `db_name` and get a
`sim.duckdb` file inside an organised output folder.

``` r

pop <- open_pop(pop_name = "demo", db_name = ":memory:")
#> duckdb is keeping downloaded extensions in a temporary directory:
#> ℹ /tmp/Rtmp6WPL0L/duckdb/extensions
#> This is removed when the R session ends, so extensions are re-downloaded each session.
#> ℹ To keep them, point `options(duckdb.extension_directory =)` or the `DUCKDB_EXTENSION_DIRECTORY` environment variable at a permanent path.
#> Opened in-memory population 'demo'

pop
#> ── tidybreed population: demo ──────────────────────────────────────────────────
#>   Database   in-memory  [connected]
#> 
#>   schema(pop) · describe_table(pop, "name")
#> ────────────────────────────────────────────────────────────────────────────────
```

Printing a population gives you a live summary. Right now there is
almost nothing in it — sections appear as you add data.

## Define the genome

``` r

pop <- pop |>
  define_genome(
    n_loci     = 500,   # total loci (SNP and QTL both come from this pool)
    n_chr      = 5,     # chromosomes
    chr_len_Mb = 100    # length of each chromosome in megabases
  )
#> Defined genome: 500 loci across 5 chromosomes | chr lengths (Mb): all equal to 100 Mb
#>   Tables written: genome_meta, genome_map, ind_haplotype, ind_genotype, ind_crossover, chr_meta

pop |> get_table("genome_meta")
#> <tidybreed_table: genome_meta>  [500 rows × 5 fields]
#> # A tibble: 10 × 5
#>    locus_id locus_name   chr chr_name  pos_bp
#>       <int> <chr>      <int> <chr>      <dbl>
#>  1        1 Locus_1        1 1         990099
#>  2        2 Locus_2        1 1        1980198
#>  3        3 Locus_3        1 1        2970297
#>  4        4 Locus_4        1 1        3960396
#>  5        5 Locus_5        1 1        4950495
#>  6        6 Locus_6        1 1        5940594
#>  7        7 Locus_7        1 1        6930693
#>  8        8 Locus_8        1 1        7920792
#>  9        9 Locus_9        1 1        8910891
#> 10       10 Locus_10       1 1        9900990
```

`genome_meta` holds one row per locus with its **physical** position
(`pos_bp`). The **genetic** map lives in its own table, because a genome
can have several maps — male vs female, or line-specific:

``` r

pop |> get_table("genome_map")
#> <tidybreed_table: genome_map>  [500 rows × 7 fields]
#> # A tibble: 10 × 7
#>    id_genome_map locus_id locus_name sex   line_name map_name pos_cM
#>            <int>    <int> <chr>      <chr> <chr>     <chr>     <dbl>
#>  1             1        1 Locus_1    NA    NA        default   0.990
#>  2             2        2 Locus_2    NA    NA        default   1.98 
#>  3             3        3 Locus_3    NA    NA        default   2.97 
#>  4             4        4 Locus_4    NA    NA        default   3.96 
#>  5             5        5 Locus_5    NA    NA        default   4.95 
#>  6             6        6 Locus_6    NA    NA        default   5.94 
#>  7             7        7 Locus_7    NA    NA        default   6.93 
#>  8             8        8 Locus_8    NA    NA        default   7.92 
#>  9             9        9 Locus_9    NA    NA        default   8.91 
#> 10            10       10 Locus_10   NA    NA        default   9.90
```

Per-chromosome inheritance rules live in `chr_meta`. Every chromosome
starts as a normal diploid autosome;
[`define_chr()`](https://austin-putz.github.io/tidybreed/reference/define_chr.md)
is how you would declare an X, Y, or mitochondrial chromosome.

``` r

pop |> get_table("chr_meta") |> collect()
#> # A tibble: 5 × 6
#>   chr_name copy_mode_M copy_mode_F hemi_parent recombines_M recombines_F
#>   <chr>    <chr>       <chr>       <chr>       <lgl>        <lgl>       
#> 1 1        full        full        NA          TRUE         TRUE        
#> 2 2        full        full        NA          TRUE         TRUE        
#> 3 3        full        full        NA          TRUE         TRUE        
#> 4 4        full        full        NA          TRUE         TRUE        
#> 5 5        full        full        NA          TRUE         TRUE
```

## Explore the database

Two helpers exist so you never have to guess what is in the database.
`schema(pop)` lists every table:

``` r

schema(pop)
#> ── Schema: demo ────────────────────────────────────────────────────────────────
#>   Use describe_table(pop, "name") for column-level details.
#> 
#>   Table                 Rows  Cols  Description 
#>   ──────────────────────────────────────────────────────────────────────────────
#>   ind_meta                 0     6  Individual-level metadata. One row per in...
#>   trait_var_comp           0     5  Genetic variance component storage. One r...
#>   genome_effects           0     7  QTL effect data in long format. One row p...
#>   phenotype_meta           0    16  Observed phenotype definitions. One row p...
#>   phenotype_components     0    17  Component definitions for composite pheno...
#>   phenotype_var_comp       0    10  Phenotype-level variance component storag...
#>   trait_meta               0     6  Genetic component trait definitions. One ...
#>   trait_effects            0    12  Fixed and random effect configurations fo...
#>   trait_random_effects     0     5  Sampled random effect levels. One row per...
#>   ind_phenotype            0     5  Phenotype records in long format. One row...
#>   ind_tbv                  0     4  True breeding values (simulation ground t...
#>   ind_ebv                  0     8  Estimated breeding values from external B...
#>   index_meta               0     5  Selection index definitions. One row per ...
#>   ind_index                0     5  Computed selection index values. One row ...
#>   ind_true_index           0     5  True selection index values computed from...
#>   genome_meta            500     5  Locus-level metadata. One row per locus. ...
#>   genome_map             500     7  Genetic map in long format. One row per (...
#>   ind_haplotype            0     7  Phased haplotypes in long format. One row...
#>   ind_genotype             0     4  Genotype dosage cache in long format. One...
#>   ind_crossover            0     6  Crossover events in long format, one row ...
#>   chr_meta                 5     6  Per-chromosome inheritance rules. One row...
```

[`describe_table()`](https://austin-putz.github.io/tidybreed/reference/describe_table.md)
explains one table’s columns:

``` r

describe_table(pop, "ind_meta")
#> ── ind_meta ─────────────────────────────────────────────────── 0 rows · 6 cols 
#>   Individual-level metadata. One row per individual. Core columns are managed by the system; user-defined columns can be added via mutate_table() or the ... arguments of add_founders() and add_offspring().
#> 
#>   Column       Type      Description 
#>   ──────────────────────────────────────────────────────────────────────────────
#>   id_ind       VARCHAR   Primary key; format '{line_name}_{n}' (e.g. 'A_1',...
#>   id_parent_1  VARCHAR   Paternal parent id_ind; NA for founders
#>   id_parent_2  VARCHAR   Maternal parent id_ind; NA for founders
#>   line_name    VARCHAR   Genetic line name (e.g. 'A', 'Holstein')
#>   sex          VARCHAR   Sex of the individual: 'M' for male, 'F' for female
#>   ploidy       UTINYINT  Genome ploidy; declared at add_founders() time (mu...
```

Remember that
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
is lazy. The filter below runs *inside* DuckDB and only the answer comes
back to R:

``` r

pop |>
  get_table("genome_meta") |>
  filter(chr == 1L) |>
  count() |>
  collect()
#> # A tibble: 1 × 1
#>       n
#>   <dbl>
#> 1   100
```

## Create founder haplotypes

Founders are sampled from a pool of haplotypes. Tagging the pool with a
`line_name` is good practice — it is what lets you run crossbreeding
later.

``` r

pop <- pop |>
  define_founder_haplotypes(
    line_name       = "A",
    n_haplotypes    = 100,
    method          = "uniform",
    min_allele_freq = 0.05,
    max_allele_freq = 0.95
  )
#> Generating 100 founder haplotypes (method = "uniform")...
#>   Created founder_haplotypes for line 'A' (100 haplotypes x 500 loci)
```

Several methods control the allele-frequency spectrum and linkage
disequilibrium:

| `method` | What it gives you |
|----|----|
| `"uniform"` | Frequencies drawn uniformly between `min_allele_freq` and `max_allele_freq` |
| `"fixed"` | Every locus at the same `allele_freq` |
| `"beta"` | Beta-distributed frequencies (`beta_shape1`, `beta_shape2`) — U-shaped is realistic |
| `"balding_nichols"` | FST-based divergence between lines (`fst`, `mean_freq`) |
| `"mosaic"` | Haplotypes built from templates, producing simple LD (`n_templates`, `switch_rate`) |
| `"gaussian_copula"` | Correlated frequencies with distance-based decay (`decay_rate`) |

## Add founder individuals

This is the first action function, and it shows why the
`tidybreed_table` pattern exists: you pipe in the haplotype pool you
want to sample *from*.

Any extra named argument — like `gen` below — becomes a new column on
`ind_meta`, written in the same transaction as the animals themselves.

``` r

pop <- pop |>
  get_table("founder_haplotypes") |>
  filter(line_name == "A") |>
  add_founders(
    n_males   = 250,
    n_females = 250,
    line_name = "A",
    gen       = 0L      # custom column
  )
#> Added new column 'gen' (INTEGER) to `ind_meta`
#> Added 500 founders (250 males, 250 females) to line 'A'

pop |> get_table("ind_meta") |> collect() |> head()
#> # A tibble: 6 × 7
#>   id_ind id_parent_1 id_parent_2 line_name sex   ploidy   gen
#>   <chr>  <chr>       <chr>       <chr>     <chr>  <int> <int>
#> 1 A_1    NA          NA          A         M          2     0
#> 2 A_2    NA          NA          A         M          2     0
#> 3 A_3    NA          NA          A         M          2     0
#> 4 A_4    NA          NA          A         M          2     0
#> 5 A_5    NA          NA          A         M          2     0
#> 6 A_6    NA          NA          A         M          2     0
```

> **Watch the type suffix.** `gen = 0L` creates an `INTEGER` column;
> `gen = 0` would create a `DOUBLE`. The same applies to `NA_integer_`,
> `NA_character_`, and friends.

Haplotypes are stored in long format — one row per individual × parent ×
locus:

``` r

pop |> get_table("ind_haplotype")
#> <tidybreed_table: ind_haplotype>  [5e+05 rows × 7 fields]
#> # A tibble: 10 × 7
#>    id_ind parent_origin strand line_origin locus_id locus_name allele
#>    <chr>          <int>  <int> <chr>          <int> <chr>       <int>
#>  1 A_1                1      1 A                  1 Locus_1         0
#>  2 A_2                1      1 A                  1 Locus_1         0
#>  3 A_3                1      1 A                  1 Locus_1         0
#>  4 A_4                1      1 A                  1 Locus_1         1
#>  5 A_5                1      1 A                  1 Locus_1         0
#>  6 A_6                1      1 A                  1 Locus_1         0
#>  7 A_7                1      1 A                  1 Locus_1         0
#>  8 A_8                1      1 A                  1 Locus_1         0
#>  9 A_9                1      1 A                  1 Locus_1         0
#> 10 A_10               1      1 A                  1 Locus_1         1
```

## Add your own columns

[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
adds or updates columns on *any* table. This is the mechanism for
keeping simulation state in the database instead of in parallel R
objects.

``` r

# Broadcast a value to every row
pop <- pop |> get_table("ind_meta") |> mutate_table(farm = "Iowa")
#> Added new column 'farm' (VARCHAR) to `ind_meta`; 500 rows set

# Update only some rows
pop <- pop |>
  get_table("ind_meta") |>
  filter(sex == "M") |>
  mutate_table(farm = "AI_Stud")
#> Warning: 'farm': replaced 250 existing values in `ind_meta` [250 of 500 rows]

pop |> get_table("ind_meta") |> count(sex, farm) |> collect()
#> # A tibble: 2 × 3
#>   sex   farm        n
#>   <chr> <chr>   <dbl>
#> 1 F     Iowa      250
#> 2 M     AI_Stud   250
```

The warning above is deliberate:
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
tells you whenever it overwrites values that were already there, so you
notice an accidental clobber instead of discovering it three generations
later.

Calling
[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
on an *empty* table still creates the column, which lets you declare a
typed schema before any data arrives:

``` r

pop |> get_table("ind_ebv") |> mutate_table(model_version = NA_character_)
```

## Define a trait

The genetic layer. `target_add_var` is the additive genetic variance the
QTL effects will be scaled to hit.

``` r

pop <- pop |>
  define_trait("ADG", target_add_var = 0.25, units = "kg/day")
#> Added trait 'ADG'.

pop |> get_table("trait_meta") |> collect()
#> # A tibble: 1 × 6
#>   id_trait trait_name description units  expressed_parent target_add_mean
#>      <int> <chr>      <chr>       <chr>  <chr>                      <dbl>
#> 1        1 ADG        NA          kg/day both                           0
```

Note what is *not* here: no mean, no residual variance, no trait type.
Those are observation-layer properties and come later.

## Assign QTL effects

Which loci are QTL is decided by the filter you pipe in. Here
chromosomes 4 and 5 carry the QTL, leaving 1–3 free for a SNP chip.

``` r

pop <- pop |>
  get_table("genome_meta") |>
  filter(chr %in% c(4L, 5L)) |>
  define_additive_effects("ADG")
#> Set additive effects for 200 QTL on trait 'ADG' (base: founder_haplotypes).

pop |> get_table("genome_effects") |> collect() |> head()
#> # A tibble: 6 × 7
#>   id_genome_effect locus_name line_name trait_name genome_effect_type
#>              <int> <chr>      <chr>     <chr>      <chr>             
#> 1                1 Locus_301  NA        ADG        additive          
#> 2                2 Locus_302  NA        ADG        additive          
#> 3                3 Locus_303  NA        ADG        additive          
#> 4                4 Locus_304  NA        ADG        additive          
#> 5                5 Locus_305  NA        ADG        additive          
#> 6                6 Locus_306  NA        ADG        additive          
#> # ℹ 2 more variables: genome_value <dbl>, base_allele_freq <dbl>
```

QTL membership is **implicit**: a locus is a QTL for a trait if it has a
row in `genome_effects`. There is no `is_QTL` flag to keep in sync.
`base_allele_freq` is stored alongside each effect because it is what
centres breeding values.

## Define the phenotype

The observation layer. With `target_add_var = 0.25` and
`residual_var = 0.75`, this trait has a heritability of 0.25.

``` r

pop <- pop |>
  define_phenotype(
    "ADG",
    type         = "continuous",
    mean         = 1.0,
    residual_var = 0.75
  )
#> Added phenotype 'ADG' (type: continuous).

pop |> get_table("phenotype_meta") |> collect() |> select(1:6)
#> # A tibble: 1 × 6
#>   id_phenotype_meta phenotype_name type        mean expressed_sex repeatable
#>               <int> <chr>          <chr>      <dbl> <chr>         <lgl>     
#> 1                 1 ADG            continuous     1 both          FALSE
```

`type` also accepts `"count"`, `"categorical"` (with `prevalence` or
`thresholds`), and `"derived_formula"`.

## Add a fixed effect

Non-genetic terms attach to the *phenotype*, not the trait. Here males
grow 0.30 kg/day faster than females, with females as the reference
level.

``` r

pop <- pop |>
  define_effect_fixed_class(
    "ADG",
    effect_name   = "sex",
    source_column = "sex",
    levels        = c(M = 0.30, F = 0)
  )
#> Added fixed-class effect 'sex' to phenotype 'ADG' (2 levels).
```

[`define_effect_fixed_cov()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_cov.md)
adds a regression on a covariate, and
[`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md)
adds a named random effect such as pen or litter.

## Compute breeding values and phenotypes

[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
computes true breeding values — the simulation’s ground truth, which you
would never know in a real population.

``` r

pop <- pop |> get_table("ind_meta") |> add_tbv("ADG")
#> Computed TBV for 500 individuals on trait 'ADG'.

pop |> get_table("ind_tbv") |> collect() |> head()
#> # A tibble: 6 × 4
#>   id_tbv id_ind trait_name tbv_value
#>    <int> <chr>  <chr>          <dbl>
#> 1      5 A_5    ADG          -0.137 
#> 2     13 A_13   ADG           0.119 
#> 3     16 A_16   ADG           0.115 
#> 4     22 A_22   ADG           0.321 
#> 5     25 A_25   ADG           0.127 
#> 6     31 A_31   ADG          -0.0346
```

[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
builds the observed record: intercept, plus fixed and random effects,
plus the breeding value, plus a sampled residual. It calls
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
internally, so you can go straight to it.

``` r

pop <- pop |> get_table("ind_meta") |> add_phenotype("ADG")
#> Computed TBV for 500 individuals on trait 'ADG'.
#> Wrote 500 phenotype records for 'ADG'.

pop |> get_table("ind_phenotype") |> collect() |> head()
#> # A tibble: 6 × 5
#>   id_phenotype id_ind phenotype_name pheno_value pheno_number
#>          <int> <chr>  <chr>                <dbl>        <int>
#> 1            1 A_1    ADG                  0.466            1
#> 2            2 A_2    ADG                  2.24             1
#> 3            3 A_3    ADG                  1.16             1
#> 4            4 A_4    ADG                  2.96             1
#> 5            5 A_5    ADG                  1.41             1
#> 6            6 A_6    ADG                  0.550            1
```

The sex effect we defined should show up as a gap of roughly 0.30
between the group means:

``` r

pop |>
  get_table("ind_phenotype") |>
  collect() |>
  left_join(
    pop |> get_table("ind_meta") |> collect() |> select(id_ind, sex),
    by = "id_ind"
  ) |>
  group_by(sex) |>
  summarise(n = n(), mean_ADG = mean(pheno_value))
#> # A tibble: 2 × 3
#>   sex       n mean_ADG
#>   <chr> <int>    <dbl>
#> 1 F       250     1.00
#> 2 M       250     1.36
```

## A second, genetically correlated trait

Real programs select on several traits at once. Define the second trait,
then supply a genetic covariance matrix.

``` r

pop <- pop |> define_trait("BF", target_add_var = 0.30, units = "mm")
#> Added trait 'BF'.

G <- matrix(
  c(0.25, 0.08,
    0.08, 0.30),
  nrow = 2,
  dimnames = list(c("ADG", "BF"), c("ADG", "BF"))
)

pop <- pop |> define_effect_cov_matrix("gen_add", G)
#> Stored 'gen_add' covariance matrix for: ADG, BF.

pop |> get_table("trait_var_comp") |> collect()
#> # A tibble: 4 × 5
#>   id_trait_var_comp effect_name trait_name_1 trait_name_2 cov_value
#>               <int> <chr>       <chr>        <chr>            <dbl>
#> 1                 1 gen_add     ADG          ADG               0.25
#> 2                 2 gen_add     ADG          BF                0.08
#> 3                 3 gen_add     BF           ADG               0.08
#> 4                 4 gen_add     BF           BF                0.3
```

Covariance matrices are stored in the database, not passed around in R.
Now one call draws correlated effects for both traits from `MVN(0, G)`:

``` r

pop <- pop |>
  get_table("genome_meta") |>
  filter(chr %in% c(4L, 5L)) |>
  define_additive_effects(c("ADG", "BF"))
#> Set correlated additive effects for traits: ADG, BF (method: shared)

pop <- pop |>
  define_phenotype("BF", type = "continuous", mean = 12, residual_var = 0.70)
#> Added phenotype 'BF' (type: continuous).
```

> This **replaces** the ADG effects assigned earlier with a fresh
> correlated draw. Re-running
> [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
> for a trait always overwrites its rows in `genome_effects`, so
> recompute anything derived from them.

``` r

pop <- pop |> get_table("ind_meta") |> add_tbv()   # no trait_name = all traits
#> Computed TBV for 500 individuals on trait 'ADG'.
#> Computed TBV for 500 individuals on trait 'BF'.

pop |>
  get_table("ind_tbv") |>
  collect() |>
  group_by(trait_name) |>
  summarise(n = n(), mean_tbv = mean(tbv_value), var_tbv = var(tbv_value))
#> # A tibble: 2 × 4
#>   trait_name     n mean_tbv var_tbv
#>   <chr>      <int>    <dbl>   <dbl>
#> 1 ADG          500   0.0262   0.188
#> 2 BF           500  -0.0204   0.287
```

The realised variances are in the neighbourhood of the 0.25 and 0.30 we
asked for, but they are not exact — and that is worth understanding
rather than glossing over. `target_add_var` scales the QTL effects so
the additive variance comes out right **in the founder haplotype pool**.
The animals you actually created are a finite sample from that pool, so
their realised variance wobbles around the target. Larger founder
populations wobble less.

## Make some offspring

Matings are described by an ordinary tibble: one row per offspring, with
`id_parent_1`, `id_parent_2`, `sex`, and `line_name` required. Any
custom `ind_meta` column can ride along.

``` r

sires <- pop |> get_table("ind_meta") |> filter(sex == "M") |> collect() |> pull(id_ind)
dams  <- pop |> get_table("ind_meta") |> filter(sex == "F") |> collect() |> pull(id_ind)

matings <- tibble::tibble(
  id_parent_1 = rep(sires[1:5], each = 4),
  id_parent_2 = dams[1:20],
  sex         = rep(c("M", "F"), 10),
  line_name   = "A",
  gen         = 1L
)

head(matings)
#> # A tibble: 6 × 5
#>   id_parent_1 id_parent_2 sex   line_name   gen
#>   <chr>       <chr>       <chr> <chr>     <int>
#> 1 A_1         A_251       M     A             1
#> 2 A_1         A_252       F     A             1
#> 3 A_1         A_253       M     A             1
#> 4 A_1         A_254       F     A             1
#> 5 A_2         A_255       M     A             1
#> 6 A_2         A_256       F     A             1
```

[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
simulates meiosis — crossovers are drawn from the genetic map in
`genome_map`, respecting the per-chromosome rules in `chr_meta`.

``` r

pop <- add_offspring(pop, matings)
#> Added 20 offspring (base_seed = 409472004)

pop |> get_table("ind_meta") |> filter(gen == 1L) |> collect() |> head()
#> # A tibble: 6 × 8
#>   id_ind id_parent_1 id_parent_2 line_name sex   ploidy   gen farm 
#>   <chr>  <chr>       <chr>       <chr>     <chr>  <int> <int> <chr>
#> 1 A_501  A_1         A_251       A         M          2     1 NA   
#> 2 A_502  A_1         A_252       A         F          2     1 NA   
#> 3 A_503  A_1         A_253       A         M          2     1 NA   
#> 4 A_504  A_1         A_254       A         F          2     1 NA   
#> 5 A_505  A_2         A_255       A         M          2     1 NA   
#> 6 A_506  A_2         A_256       A         F          2     1 NA
```

Because you build the mating plan yourself, any design is possible —
nested, factorial, reciprocal crosses, assortative mating — without the
package needing a special argument for each.

``` r

pop
#> ── tidybreed population: demo ──────────────────────────────────────────────────
#>   Database   in-memory  [connected]
#> 
#>   Genome     500 loci · 5 chr · 495 Mb · 490 cM
#>              founder pool: 100 haplotypes
#>   Model      2 traits · 2 phenotypes · 200 QTL
#> 
#>   Individuals  520
#>     by sex     260 F · 260 M
#> 
#>   Records    phenotypes 500 · TBV 1,000
#> 
#>   schema(pop) · describe_table(pop, "name")
#> ────────────────────────────────────────────────────────────────────────────────
```

## SNP chips and genotypes

A chip is a named set of loci. Mark them with a filter:

``` r

pop <- pop |>
  get_table("genome_meta") |>
  filter(chr %in% c(1L, 2L, 3L)) |>
  define_chip("50K")
#> Defined chip '50K' with 300 SNPs in column 'is_50K'

pop |> get_table("genome_meta") |> count(is_50K) |> collect()
#> # A tibble: 2 × 2
#>   is_50K     n
#>   <lgl>  <dbl>
#> 1 FALSE    200
#> 2 TRUE     300
```

[`add_genotypes()`](https://austin-putz.github.io/tidybreed/reference/add_genotypes.md)
records *which animals were genotyped* — here, only the offspring:

``` r

pop <- pop |>
  get_table("ind_meta") |>
  filter(gen == 1L) |>
  add_genotypes("50K")
#> Chip '50K': 20 animal(s) now marked as genotyped.

pop |> get_table("ind_meta") |> count(has_50K) |> collect()
#> # A tibble: 2 × 2
#>   has_50K     n
#>   <lgl>   <dbl>
#> 1 FALSE     500
#> 2 TRUE       20
```

[`extract_genotypes()`](https://austin-putz.github.io/tidybreed/reference/extract_genotypes.md)
returns the familiar wide 0/1/2 dosage matrix, ready for GBLUP or
export:

``` r

geno <- pop |> get_table("ind_meta") |> extract_genotypes("50K")

dim(geno)
#> [1]  20 301
geno[1:5, 1:6]
#> # A tibble: 5 × 6
#>   id_ind locus_1 locus_2 locus_3 locus_4 locus_5
#>   <chr>    <int>   <int>   <int>   <int>   <int>
#> 1 A_501        0       2       1       1       0
#> 2 A_502        1       2       1       2       0
#> 3 A_503        1       2       2       1       0
#> 4 A_504        1       2       1       2       0
#> 5 A_505        0       1       2       1       0
```

Haplotypes are the source of truth and dosages are derived on demand, so
nothing is stored twice. If you need dosages repeatedly,
[`add_dosage()`](https://austin-putz.github.io/tidybreed/reference/add_dosage.md)
caches them in `ind_genotype`.

## Selection index

Register the weights:

``` r

pop <- pop |>
  define_index("terminal",
               trait_names = c("ADG", "BF"),
               index_wts   = c(1.2, -0.8))
#> Defined index 'terminal': ADG (1.2), BF (-0.8)

pop |> get_table("index_meta") |> collect()
#> # A tibble: 4 × 5
#>   id_index_name index_name trait_name index_weight economic_weight
#>           <int> <chr>      <chr>             <dbl>           <dbl>
#> 1             1 NA         ADG                NA                 0
#> 2             2 NA         BF                 NA                 0
#> 3             3 terminal   ADG                 1.2              NA
#> 4             4 terminal   BF                 -0.8              NA
```

The rows with a missing `index_name` are the global economic-weight
entries that
[`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md)
writes for every trait; the `terminal` rows are our index.

Because we know the true breeding values, we can compute the **true**
index — useful for measuring how well selection actually worked:

``` r

pop <- pop |> get_table("ind_meta") |> add_tbv(index_names = "terminal")
#> Computed TBV for 520 individuals on trait 'ADG'.
#> Computed TBV for 520 individuals on trait 'BF'.
#> Computed true index 'terminal' (index) for 520 individuals.

pop |> get_table("ind_true_index") |> collect() |> head()
#> # A tibble: 6 × 5
#>   id_true_index id_ind index_name weight_type true_index_value
#>           <int> <chr>  <chr>      <chr>                  <dbl>
#> 1             1 A_1    terminal   index                 0.483 
#> 2             2 A_10   terminal   index                 0.726 
#> 3             3 A_100  terminal   index                -0.362 
#> 4             4 A_101  terminal   index                -0.375 
#> 5             5 A_102  terminal   index                 0.0360
#> 6             6 A_103  terminal   index                 0.0964
```

[`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md)
does the same arithmetic on any table of values. In a real program you
would run it on `ind_ebv` after a BLUP evaluation; here we use the TBVs.
**Filter first** so there is exactly one value per individual per trait:

``` r

pop <- pop |>
  get_table("ind_tbv") |>
  filter(trait_name %in% c("ADG", "BF")) |>
  add_index("terminal")
#> Computing index 'terminal': ADG (wt=1.2), BF (wt=-0.8)
#> Added index 'terminal' (run #1) for 520 individuals

pop |>
  get_table("ind_index") |>
  collect() |>
  arrange(desc(index_value)) |>
  head(5)
#> # A tibble: 5 × 5
#>   id_index id_ind index_name index_number index_value
#>      <int> <chr>  <chr>             <int>       <dbl>
#> 1      460 A_512  terminal              1        1.56
#> 2      463 A_515  terminal              1        1.46
#> 3      404 A_462  terminal              1        1.36
#> 4      471 A_54   terminal              1        1.32
#> 5       70 A_161  terminal              1        1.20
```

Those top animals are your selection candidates. Pull their IDs, build a
new mating plan, call
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md),
and you have closed the loop on a breeding cycle.

## Close the population

``` r

close_pop(pop)
```

For an on-disk run,
[`close_pop()`](https://austin-putz.github.io/tidybreed/reference/close_pop.md)
flushes everything to the `.duckdb` file and you can pick the simulation
back up later:

``` r

pop <- restore_pop("tidybreed_output/sim.duckdb")
```

Nothing needs to be re-specified — the genome, traits, variance
components, and every animal are already in the file. This is the payoff
of keeping configuration in the database rather than in R.

## Where to go next

This vignette covered the main path. The package has more, grouped by
prefix:

| Prefix | Meaning |
|----|----|
| `open_` / `restore_` | Create or reopen a population |
| `define_` | Write model configuration — how the simulation should behave |
| `add_` | Write simulation output — data the model produced |
| `mutate_` | Add or update columns on an existing table |
| `get_` / `extract_` | Read data out, lazily or in analysis format |
| `remove_` / `archive_` | Delete rows, or stamp a finished replicate |

Functions not shown above, with a pointer to their help page:

| Function | What it does |
|----|----|
| [`?add_ebv`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md) | Run BLUPF90 (or parent average) and store estimated breeding values |
| [`?define_chr`](https://austin-putz.github.io/tidybreed/reference/define_chr.md) | Sex chromosomes, organelles, and achiasmatic meiosis |
| [`?define_phenotype`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md) | Composite traits: maternal effects and social genetic effects, via `components` |
| [`?define_residual_cov`](https://austin-putz.github.io/tidybreed/reference/define_residual_cov.md) | Correlated or heterogeneous residual (co)variances |
| [`?define_effect_random`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md) | Named random effects such as pen, litter, or herd-year-season |
| [`?define_effect_fixed_cov`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_cov.md) | Fixed regression on a covariate |
| [`?define_trait_simple`](https://austin-putz.github.io/tidybreed/reference/define_trait_simple.md) | Shortcut wrapper: [`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md) + [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md) |
| [`?archive_replicate`](https://austin-putz.github.io/tidybreed/reference/archive_replicate.md) | Collect many replicates into one archive database |
| [`?remove_rows`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md) | Delete rows safely across related tables |
| [`?mutate_group_seq`](https://austin-putz.github.io/tidybreed/reference/mutate_group_seq.md) | Litter and group utilities (`mutate_group_*`) |
| [`?define_table`](https://austin-putz.github.io/tidybreed/reference/define_table.md) | Create your own tables inside the population |

The project README covers installation, global options, the output
directory layout, and scenario YAML files in more detail.

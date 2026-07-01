# `tidybreed` R Package

<!-- badges: start -->
<!-- badges: end -->

A pipe-friendly (`%>%` or `|>`) R package for breeding program simulation backed by [DuckDB](https://duckdb.org). 
Design large-scale genomic simulations without running out of memory — all data lives on disk in a DuckDB database 
and is queried lazily via `dplyr`. 

## License

**MIT License**

`tidybreed` is released under the [MIT License](LICENSE). You are free to use,
modify, and distribute it, including in commercial and proprietary projects, provided
the copyright notice is retained.

> **NO WARRANTY** — This software is provided **as-is**, without warranty of any kind, express or
> implied. The authors accept **no liability** for any damages or losses arising from its use.

## Design

| Main | → | Reason or Description |
| ---- | - | --------------------- |
| [R](https://www.r-project.org) | → | Standard for most scientists; flexible design of custom breeding programs |
| **DuckDB** | → | Columnar, embedded, no server needed; handles datasets larger than RAM |
| Pipe everything | → | Filter individuals, add phenotypes/genotypes/EBVs/index with `tidyverse` verbs |
| Fully customizable | → | Add any column to any table with `mutate_table()`; query with standard `dplyr` |

## Installation

Install [pak](https://pak.r-lib.org/) then install `tidybreed` from [GitHub](https://github.com/austin-putz/tidybreed/):

```r
install.packages("pak", repos = "https://packagemanager.posit.co/cran/latest")
library(pak)
pak::pak("austin-putz/tidybreed")
library(tidybreed)
```

> **WARNING:** Pre-`v1.0.0` packages are considered beta and subject to breaking changes. **Pin your version** to avoid surprises.

```r
packageVersion("tidybreed")   # check installed version, e.g. '0.46.1'
pak::pak("austin-putz/tidybreed@v0.46.1")   # pin to a specific release
```

Browse all releases on the [GitHub Releases page](https://github.com/austin-putz/tidybreed/releases).

---

## Global Options

Set package-wide options at the top of your script (or in a startup / config
script) so simulations can be parameterized without changing function calls
throughout the codebase. Every option has a built-in default and is entirely
optional — set only the ones you need. Options prefixed `tidybreed.` are read
by `open_pop()` (folder layout) and `archive_replicate()` (multi-replicate
archiving).

```r
options(
  tidybreed.pop_name        = "my_project",                    # label stored on the pop object
  tidybreed.base_dir        = "~/path/to/project/",             # root folder (layer 1)
  tidybreed.output          = "tidybreed_output",                # output subfolder name (layer 2)
  tidybreed.scenario        = "baseline",                        # scenario subfolder name (layer 3)
  tidybreed.tools           = c("blupf90", "plink"),             # external tool subfolders (layer 4)
  tidybreed.db_name         = "sim.duckdb",                      # working DuckDB file name
  tidybreed.replicate       = 1L,                                # current replicate number
  tidybreed.archive_path    = "~/path/to/project/archive/",      # directory for the archive DuckDB file
  tidybreed.db_name_archive = "all_reps.duckdb",                 # archive DuckDB file name
  tidybreed.quiet           = FALSE                              # suppress the startup banner
)
```

| Option | Default | Used by | Description |
|--------|---------|---------|-------------|
| `tidybreed.pop_name` | `"sim"` | `open_pop()` | Optional label stored on the `tidybreed_pop` object (`pop$pop_name`); shown in `print(pop)`. Purely descriptive — does not affect file paths. |
| `tidybreed.base_dir` | `getwd()` | `open_pop()` | Root folder (layer 1) under which the output folder, scenario folder, and database file are created. Set this once per project so every script writes to the same place. |
| `tidybreed.output` | `"tidybreed_output"` | `open_pop()` | Name of the output subfolder (layer 2) created inside `base_dir`, e.g. `<base_dir>/tidybreed_output/`. Groups all simulation runs together, separate from other project files. |
| `tidybreed.scenario` | `NULL` | `open_pop()` | Name of the scenario subfolder (layer 3), e.g. `<base_dir>/<output>/<scenario>/`. When `NULL` (default), a `YYYYMMDD_HHMMSS` timestamp folder is auto-generated so every run is isolated and nothing gets overwritten by accident. Set explicitly (e.g. `"baseline"`) to reuse the same folder across runs, such as in an HPC array job. |
| `tidybreed.tools` | `NULL` | `open_pop()` | Character vector of external-tool subfolder names to create at layer 4 inside the scenario folder, e.g. `c("blupf90", "plink")` creates `.../blupf90/` and `.../plink/` for those tools' input/output files. `NULL` skips tool folder creation. |
| `tidybreed.db_name` | `"sim.duckdb"` | `open_pop()` | File name of the **working** DuckDB database, placed inside the scenario folder. Use `":memory:"` for an in-memory database (skips all folder creation entirely — useful for tests and quick prototyping). |
| `tidybreed.replicate` | `1L` | `archive_replicate()` | Integer replicate number stamped on every row written to the archive by `archive_replicate()`. Auto-increments by 1 after each successful archive call, so a loop over replicates does not need to manage a counter manually. |
| `tidybreed.archive_path` | `NULL` | `archive_replicate()` | Directory where the **archive** DuckDB file is written (as opposed to the working database set by `base_dir`/`output`/`scenario`/`db_name`). Combined with `tidybreed.db_name_archive` to form the full archive file path. When `NULL`, the archive lands next to the working database instead (see `tidybreed.db_name_archive`). |
| `tidybreed.db_name_archive` | `NULL` | `archive_replicate()` | File name of the **archive** DuckDB database that `archive_replicate()` appends each replicate's results to (multi-replicate analysis file, distinct from the per-run working database). When `NULL`, `archive_replicate()` skips writing an archive entirely and only resets the working database. |
| `tidybreed.quiet` | `FALSE` | package startup (`.onAttach`) | When `TRUE`, suppresses the `tidybreed` startup banner printed on `library(tidybreed)`. Useful for scripts, RMarkdown, and non-interactive/batch jobs. |

> **Archive path resolution** (used by `archive_replicate()`) — the first non-`NULL` result wins:
> 1. An explicit `archive_path` argument passed directly to `archive_replicate()`.
> 2. `file.path(tidybreed.archive_path, tidybreed.db_name_archive)` when both options are set.
> 3. `file.path(dirname(pop$db_path), tidybreed.db_name_archive)` when only `tidybreed.db_name_archive` is set — the archive lands next to the working database.
> 4. If `tidybreed.db_name_archive` is also `NULL`, no archive is written; only the reset phases of `archive_replicate()` run.

---

## Core Concept

Every function accepts and returns a `tidybreed_pop` object — a thin wrapper around a DuckDB connection. Chain operations with `|>` (or `%>%`):

```r
pop <- pop |>
  get_table("ind_meta") |>       # identify which DB table to work with
  filter(sex == "F", gen == 1L) |>  # filter to the individuals you want
  add_phenotype("ADG")           # add records to ind_phenotype
```

All tables are queryable at any time. Use `collect()` to pull results into R as a `tibble`.

---

## Simulation Workflow

### 1. Open a population

`open_pop()` creates or re-opens a DuckDB-backed population object. Use `clean = TRUE` to start fresh.

```r
pop <- open_pop(
  clean = TRUE   # wipe and recreate the database
)

print(pop)
```

### 2. Define the genome

```r
pop <- pop |>
  define_genome(
    n_loci     = 10000,  # total number of loci (SNP + QTL)
    n_chr      = 18,     # number of chromosomes
    chr_len_Mb = 50      # length of each chromosome in megabases
  )

pop |> get_table("genome_meta")   # 1 row per locus
```

### 3. Define founder haplotypes

Generate haplotype pools for each genetic line. Multiple methods are available to control allele frequency distributions.

```r
# Uniform allele frequencies between min and max
pop <- pop |>
  define_founder_haplotypes(
    line_name       = "A",
    n_haplotypes    = 1000,
    method          = "uniform",
    min_allele_freq = 0.01,
    max_allele_freq = 0.99
  )

# Fixed frequency (all loci p = 0.5)
pop <- pop |>
  define_founder_haplotypes(line_name = "B", n_haplotypes = 1000,
    method = "fixed", allele_freq = 0.5)

# Beta distribution
pop <- pop |>
  define_founder_haplotypes(line_name = "C", n_haplotypes = 1000,
    method = "beta", beta_shape1 = 0.5, beta_shape2 = 0.5)

# Balding-Nichols (FST-based)
pop <- pop |>
  define_founder_haplotypes(line_name = "D", n_haplotypes = 1000,
    method = "balding_nichols", fst = 0.1, mean_freq = 0.5)

# Mosaic (introduces simple LD)
pop <- pop |>
  define_founder_haplotypes(line_name = "E", n_haplotypes = 1000,
    method = "mosaic", n_templates = 32, switch_rate = 1.0)

# Gaussian copula
pop <- pop |>
  define_founder_haplotypes(line_name = "F", n_haplotypes = 1000,
    method = "gaussian_copula", decay_rate = 0.25)
```

### 4. Add founder individuals

Sample haplotypes for founders. Pass any custom `ind_meta` column as `...` arguments — they are written atomically with the new rows.

```r
pop <- pop |>
  get_table("founder_haplotypes") |>
  filter(line_name == "A") |>
  add_founders(
    n_males    = 400,
    n_females  = 1600,
    line_name  = "A",
    birth_date = sampled_birth_dates,   # Date vector (one per founder)
    alive      = TRUE,
    active     = FALSE
  )
```

### 5. Add custom columns with `mutate_table()`

The real power of `tidybreed` is freely adding columns to any table so your simulation state is stored in the database — no parallel R objects to maintain.

```r
# Add user-defined fields to ind_meta (declare schema before data arrives)
pop |>
  get_table("ind_meta") |>
  mutate_table(
    status        = NA_character_,   # VARCHAR: production status
    birth_date    = as.Date(NA),     # DATE
    puberty_date  = as.Date(NA),
    mate_date     = as.Date(NA),
    farrow_date   = as.Date(NA),
    wean_date     = as.Date(NA),
    cull_date     = as.Date(NA),
    alive         = TRUE,            # BOOLEAN with default value
    active        = FALSE,
    .set_default  = TRUE             # write this value as the column default
  )
```

Types are inferred from R values:

| R value | DuckDB type |
|---------|-------------|
| `0L`, `NA_integer_` | `INTEGER` |
| `0.0`, `NA_real_` | `DOUBLE` |
| `TRUE`/`FALSE` | `BOOLEAN` |
| `"text"`, `NA_character_` | `VARCHAR` |
| `as.Date(...)` | `DATE` |
| `Sys.time()` | `TIMESTAMP` |

Update filtered rows after data exists:

```r
# Mark animals that reached off-test today
pop |>
  get_table("ind_meta") |>
  filter(sex == "F", off_test_date == cur_date) |>
  mutate_table(status = "after-test-gilt")
```

Add per-individual descriptions to columns for documentation:

```r
pop |>
  get_table("ind_meta") |>
  define_schema_description("status",     "Production status (e.g. gestation, lactation)") |>
  define_schema_description("birth_date", "Date of birth") |>
  define_schema_description("alive",      "Is the animal alive?")

pop |> describe_table("ind_meta")   # view all column descriptions
schema(pop)                          # view all table schemas
```

### 6. Define a SNP chip

Filter `genome_meta` to the loci you want, then call `define_chip()`. Adds an `is_<chip_name>` boolean column to `genome_meta`.

```r
# Random 9,000-locus chip
pop |>
  get_table("genome_meta") |>
  slice_sample(n = 9000) |>
  define_chip(chip_name = "9k")

# QTL are defined as loci NOT on the chip
pop |>
  get_table("genome_meta") |>
  filter(is_9k != TRUE) |>          # non-chip loci become QTL candidates
  define_additive_effects(...)
```

### 7. Define variance components

Use a single entry point for all variance/covariance matrices. `effect_name = "gen_add"` routes to `trait_var_comp`; `"residual"` and named random effects (e.g. `"pen"`) route to `phenotype_var_comp`.

```r
# Additive genetic (co)variance — 3 traits (ADG, WWD, WWM)
vars.mat.add <- matrix(c(
  0.0045, 0.00, 0.00,
  0.00,   0.04, 0.00,
  0.00,   0.00, 0.13),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("ADG", "WWD", "WWM"),
                  c("ADG", "WWD", "WWM")))

pop <- pop |>
  define_effect_cov_matrix(effect_name = "gen_add", cov_matrix = vars.mat.add)

# Residual (co)variance — 2 phenotypes (ADG, WW)
vars.mat.res <- matrix(c(
  0.0067, 0.00,
  0.00,   0.45),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("ADG", "WW"),
                  c("ADG", "WW")))

pop <- pop |>
  define_effect_cov_matrix(effect_name = "residual", cov_matrix = vars.mat.res)

# Named random effect — pen effect on ADG only (1x1)
vars.mat.pen <- matrix(0.0005, nrow = 1, byrow = TRUE,
  dimnames = list("ADG", "ADG"))

pop <- pop |>
  define_effect_cov_matrix(effect_name = "pen", cov_matrix = vars.mat.pen)
```

### 8. Define traits and phenotypes (two-layer design)

**Genetic layer** (`define_trait`) — one row per underlying genetic quantity:

```r
pop <- pop |>
  define_trait(
    trait_name      = "ADG",
    description     = "Average Daily Gain",
    units           = "kg/d",
    target_add_mean = 0,      # TBV mean in base population
    overwrite       = TRUE
  )
```

**Observation layer** (`define_phenotype`) — what animals actually receive records for:

```r
pop <- pop |>
  define_phenotype(
    phenotype_name = "ADG",
    type           = "continuous",   # "continuous", "count", "categorical", "derived_formula"
    mean           = 0.92,           # phenotypic population mean
    expressed_sex  = "both",         # "both", "M", or "F"
    min_value      = 0,
    overwrite      = TRUE
  )
```

For composite phenotypes (e.g. weaning weight = direct + maternal):

```r
pop <- pop |>
  define_phenotype(
    phenotype_name           = "WW",
    type                     = "continuous",
    formula_tbv              = "WWD + dam(WWM)",   # DSL: self + dam component
    mean                     = 6.0,
    expressed_sex            = "both",
    missing_component_action = "skip",   # skip founders with no dam record
    overwrite                = TRUE
  )
```

Derived phenotypes computed from other phenotypes:

```r
pop <- pop |>
  define_phenotype(
    phenotype_name = "FCR",
    type           = "derived_formula",
    formula        = "ADFI / ADG",
    expressed_sex  = "both",
    overwrite      = TRUE
  )
```

Add QTL effects for a trait (or multiple correlated traits at once):

```r
# Single trait
pop |>
  get_table("genome_meta") |>
  filter(is_9k != TRUE) |>           # QTL are non-chip loci
  define_additive_effects(
    trait_name      = "ADG",
    distribution    = "normal",
    scale_to_target = TRUE,          # scale to target additive variance
    base            = "current_pop"  # standardize to current animals
  )

# Multiple correlated traits in one call (draws from MVN with G matrix)
pop |>
  get_table("genome_meta") |>
  filter(is_9k != TRUE) |>
  define_additive_effects(
    trait_name      = c("WWD", "WWM"),
    distribution    = "normal",
    scale_to_target = TRUE,
    base            = "current_pop"
  )
```

Add fixed effects:

```r
pop |>
  define_effect_fixed_class(
    "ADG",
    effect_name    = "sex",
    source_column  = "sex",
    levels         = c(M = 0.08, F = 0),   # male advantage in ADG
    source_table   = "ind_meta",
    overwrite      = TRUE
  )
```

### 9. Calculate True Breeding Values

```r
# All animals in ind_meta
pop <- pop |>
  get_table("ind_meta") |>
  add_tbv(trait_name = "ADG")

# Check means by trait
pop |>
  get_table("ind_tbv") |>
  collect() |>
  group_by(trait_name) |>
  summarise(MeanTBV = mean(tbv_value))
```

Define a selection index first, then compute true index values from TBVs (ground truth for monitoring genetic trend):

```r
# Define the index weights (must exist before add_tbv uses index_names)
pop |>
  define_index(
    index_name   = "maternal",
    trait_names  = c("ADG", "WWD", "WWM"),
    index_wts    = c(ADG = 0.5, WWD = 0.3, WWM = 0.2),
    economic_wts = c(ADG = 0.5, WWD = 0.3, WWM = 0.2)
  )

# Compute true index values from TBVs and write to ind_true_index
pop |>
  get_table("ind_meta") |>
  add_tbv(index_names = "maternal")

pop |> get_table("ind_true_index") |> collect()
```

### 10. Add Genotypes

```r
# Add 9k genotypes to all animals
pop |>
  get_table("ind_meta") |>
  add_genotypes(chip_name = "9k")

# Extract genotypes for downstream analysis
pop |>
  get_table("ind_meta") |>
  extract_genotypes(chip_name = "9k")
```

### 11. Add Phenotypes

```r
# Phenotype animals that reached off-test today
pop |>
  get_table("ind_meta") |>
  filter(off_test_date == cur_date) |>
  add_phenotype(
    phenotype_name = c("ADG", "BF"),
    pheno_date     = cur_date
  )

# Phenotype all females for age at puberty (sex-limited trait)
pop |>
  get_table("ind_meta") |>
  filter(sex == "F") |>
  add_phenotype(phenotype_name = "AP")
```

Check counts:

```r
pop |> get_table("ind_phenotype") |> filter(phenotype_name == "ADG") |> collect() |> count()
```

### 12. Run Evaluations (EBVs)

Run [BLUPF90](https://nce.ads.uga.edu/wiki/doku.php) to estimate breeding values:

```r
pop <- pop |>
  get_table("ind_meta") |>
  add_ebv(
    "ADG",
    software  = "blupf90",
    model     = "blup",
    phenotype = pop |>
      get_table("ind_phenotype") |>
      filter(pheno_date < cur_date | is.na(pheno_date)),  # exclude future records
    eval_date = cur_date
  )

pop |> get_table("ind_ebv") |> filter(trait_name == "ADG")
```

### 13. Define and Calculate a Selection Index

```r
# Define a multi-trait selection index
pop |>
  define_index(
    index_name   = "maternal",
    trait_names  = c("AP", "NW", "ADG", "ADFI", "BF", "WWD", "WWM"),
    index_wts    = b_adjusted,    # economic index weights (named vector)
    economic_wts = b_adjusted
  )

# Calculate index values from latest EBVs (run after add_ebv)
pop |>
  get_table("ind_ebv") |>
  filter(eval_date == cur_date) |>
  add_index("maternal", index_date = cur_date)

pop |> get_table("ind_index")
```

### 14. Select Parents

Pull candidate IDs and then select from any table:

```r
# Step 1: identify candidates
male_candidates <- pop |>
  get_table("ind_meta") |>
  filter(sex == "M", status %in% c("after-test-boar", "breeding-boar")) |>
  pull(id_ind)

# Step 2: select top animals by index
selected_males <- pop |>
  get_table("ind_index") |>
  filter(id_ind %in% male_candidates, index_date == latest_index_date) |>
  slice_max(index_value, n = 10) |>
  pull(id_ind)

# Update status
pop |>
  get_table("ind_meta") |>
  filter(id_ind %in% selected_males) |>
  mutate_table(status = "breeding-boar")
```

### 15. Add Offspring

Build a mating plan as a `tibble` (one row per offspring), then call `add_offspring()`:

```r
data.matings <- tibble(
  id_parent_1   = rep(sampled_sires, times = nw_per_dam),  # sire IDs
  id_parent_2   = rep(selected_dams, times = nw_per_dam),  # dam IDs
  line_name     = "A",
  sex           = "F",                  # force sex with sexed semen
  conc_date     = cur_date,
  birth_date    = cur_date + 116L,
  on_test_date  = cur_date + 116L + 70L,
  off_test_date = cur_date + 116L + 160L
)

pop |> add_offspring(data.matings)
```

### 16. Utility: mutate_group helpers

Add sequence numbers, named labels, or concatenated strings within groups:

```r
# Sequential number within each litter (group by dam ID)
pop |>
  get_table("ind_meta") |>
  filter(birth_date == cur_date) |>
  mutate_group_seq(
    group_col  = "id_parent_2",   # group by dam
    new_col    = "piglet_number"
  )

# Named labels within a group
pop |>
  get_table("ind_meta") |>
  mutate_group_named(
    group_col  = "id_parent_2",
    name_col   = "id_ind",
    new_col    = "litter_members"
  )
```

### 17. Archive and Restore

Save a replicate's DuckDB to an archive database for multi-replicate analysis:

```r
archive_replicate(pop, rep = 1L)
```

Restore a population from an existing DuckDB file (e.g. to resume a run):

```r
pop <- restore_pop(db_path = "~/path/to/project/tidybreed_output/sim.duckdb")
```

---

## Database Tables

| Table | Rows | Description |
|-------|------|-------------|
| `genome_meta` | 1 per locus | Locus positions; user chip columns added via `define_chip()` |
| `founder_haplotypes` | 1 per (haplotype × locus) | Haplotype pool sampled by `add_founders()` |
| `genome_haplotype` | 2 per individual | Phased haplotypes (paternal / maternal) |
| `genome_genotype` | 1 per individual | 0/1/2-encoded genotypes (SNP chip loci) |
| `genome_effects` | 1 per (locus × trait × effect type) | Additive QTL effect sizes |
| `ind_meta` | 1 per individual | Pedigree, sex, line; user date/status columns added via `mutate_table()` |
| `ind_phenotype` | 1 per (individual × phenotype record) | Long-format phenotype records |
| `ind_tbv` | 1 per (individual × trait) | True breeding values (simulation ground truth) |
| `ind_true_index` | 1 per (individual × index × weight type) | True index values from TBVs |
| `ind_ebv` | 1 per (individual × trait × evaluation) | Estimated breeding values from BLUP/GBLUP |
| `ind_index` | 1 per (individual × index × run) | Computed selection index values |
| `trait_meta` | 1 per genetic trait | Genetic-layer configuration (variance target, units) |
| `trait_var_comp` | 1 per (effect × trait pair) | Additive genetic (co)variance matrix entries |
| `trait_effects` | 1 per (trait × effect) | Fixed and random model effects |
| `phenotype_meta` | 1 per observed phenotype | Observation-layer config (mean, type, sex expression) |
| `phenotype_components` | 1 per (phenotype × component) | Composite phenotype wiring (self/dam/sire/group) |
| `phenotype_var_comp` | 1 per (effect × phenotype pair) | Residual and random effect (co)variance entries |
| `index_meta` | 1 per (index × trait) | Selection index weight definitions |

---

## Function Overview

### Population & Genome

| Function | Purpose |
|----------|---------|
| `open_pop()` | Open (or create) a DuckDB-backed population object |
| `define_genome()` | Define loci, chromosomes, and positions in `genome_meta` |
| `define_founder_haplotypes()` | Generate haplotype pools per line (`uniform`, `fixed`, `beta`, `balding_nichols`, `mosaic`, `gaussian_copula`) |
| `restore_pop()` | Reopen an existing population from a DuckDB file |
| `close_pop()` | Safely close the DuckDB connection |
| `print.tidybreed_pop()` | Print population summary |

### Individuals

| Function | Purpose |
|----------|---------|
| `add_founders()` | Sample haplotypes and create founder rows in `ind_meta` |
| `add_offspring()` | Add progeny given a mating-plan `tibble` (1 row per offspring) |

### Tables & Queries

| Function | Purpose |
|----------|---------|
| `get_table()` | Return a lazy `dplyr` tbl from any database table |
| `mutate_table()` | Add or update columns in any table (scalar or vector) |
| `remove_rows()` | Delete rows from a table by filter (with safety confirmation) |
| `schema()` | Print all table schemas |
| `describe_table()` | Print column descriptions for a table |
| `define_schema_description()` | Register a description string for a user-defined column |

### Genome & Chips

| Function | Purpose |
|----------|---------|
| `define_chip()` | Mark filtered loci as members of a named SNP chip |
| `add_genotypes()` | Write 0/1/2 genotypes to `genome_genotype` for a chip |
| `extract_genotypes()` | Pull genotypes into R for a chip or set of QTL effects |

### Traits & Model Configuration

| Function | Purpose |
|----------|---------|
| `define_trait()` | Register a genetic-layer trait in `trait_meta` |
| `define_trait_simple()` | Convenience wrapper: `define_trait()` + `define_additive_effects()` |
| `define_phenotype()` | Register an observed phenotype in `phenotype_meta` |
| `define_additive_effects()` | Assign QTL effects to filtered loci (single or correlated multi-trait) |
| `define_effect_cov_matrix()` | Load a (co)variance matrix into `trait_var_comp` or `phenotype_var_comp` |
| `define_effect_fixed_class()` | Add a discrete fixed-effect level-to-shift mapping |
| `define_effect_fixed_cov()` | Add a linear regression fixed covariate |
| `define_effect_random()` | Add a named random effect |
| `define_effect_intercept()` | Set the overall phenotypic mean (intercept) |
| `define_residual_cov()` | Write residual (co)variance entries to `phenotype_var_comp` |

### Simulation Output

| Function | Purpose |
|----------|---------|
| `add_tbv()` | Compute and store true breeding values in `ind_tbv` |
| `add_phenotype()` | Sample and store phenotype records in `ind_phenotype` |
| `add_ebv()` | Run BLUPF90 or parent average; store results in `ind_ebv` |
| `add_index()` | Compute weighted index from `ind_ebv` (or any table); store in `ind_index` |

### Selection Index

| Function | Purpose |
|----------|---------|
| `define_index()` | Register index weights in `index_meta` |
| `add_index()` | Calculate and store index values from EBVs or TBVs |

### Group / Litter Utilities

| Function | Purpose |
|----------|---------|
| `mutate_group_seq()` | Add a within-group sequence number column |
| `mutate_group_named()` | Add a column with named labels within groups |
| `mutate_group_concatenate()` | Concatenate values within groups into a string column |

### Replication & Archiving

| Function | Purpose |
|----------|---------|
| `archive_replicate()` | Save a replicate's database to an archive DuckDB |
| `summary_pop()` | Print a structured summary of the population state |

# Implementation Summary: `mutate_genome_meta()` and `define_chip()`

## Overview

Successfully implemented two new functions following the approved naming philosophy:

1. **`mutate_genome_meta()`** - General-purpose column mutation for genome_meta table
2. **`define_chip()`** - Convenient SNP chip definition helper

## Files Created

### Core Functions
- ✅ `R/mutate_genome_meta.R` - Main mutation function (mirrors `mutate_ind_meta()`)
- ✅ `R/define_chip.R` - Chip definition helper function
- ✅ `R/chip_helpers.R` - Internal helper functions for chip selection

### Tests
- ✅ `tests/testthat/test-mutate_genome_meta.R` - 10 test cases
- ✅ `tests/testthat/test-define_chip.R` - 15 test cases

### Documentation
- ✅ Comprehensive roxygen2 documentation for all functions
- ✅ Multiple usage examples in each function

## Next Steps

### 1. Update Package Documentation

Run this in R to regenerate NAMESPACE and documentation:

```r
devtools::document()
```

This will:
- Export `mutate_genome_meta()` and `define_chip()`
- Generate man pages
- Update NAMESPACE

### 2. Run Tests

```r
devtools::test()
```

Or test specific files:

```r
testthat::test_file("tests/testthat/test-mutate_genome_meta.R")
testthat::test_file("tests/testthat/test-define_chip.R")
```

### 3. Run Integration Test

```r
source("test_new_functions.R")
```

This will run a complete workflow demonstrating:
- Defining multiple chips
- Adding QTL effects to chip SNPs
- Integrating with `add_founders()` and `mutate_ind_meta()`
- Querying chip genotypes

### 4. Update DESIGN.md (Optional)

Lines 222-243 in DESIGN.md still reference the old `add_SNP_chip()` function. Update to use `define_chip()` instead:

**Old:**
```r
pop_A <- pop_A %>%
  add_SNP_chip(chip_name = "50K", n = 500)
```

**New:**
```r
pop_A <- pop_A %>%
  define_chip(name = "50k", n_snp = 500, method = "random")
```

## Function Overview

### `mutate_genome_meta(pop, ...)`

**Purpose**: Add or modify columns in genome_meta table

**Examples:**
```r
# Scalar (all loci same value)
pop %>% mutate_genome_meta(is_QTL = FALSE)

# Vector (different value per locus)
pop %>% mutate_genome_meta(
  maf = runif(n_loci, 0.01, 0.5),
  effect = rnorm(n_loci, 0, 1)
)

# Multiple fields
pop %>% mutate_genome_meta(
  is_50k = chip_indicator,
  is_HD = hd_indicator,
  gene_name = gene_names
)
```

**Reserved columns** (cannot modify):
- `locus_id`
- `locus_name`
- `chr`
- `chr_name`
- `pos_Mb`

### `define_chip(pop, name, ...)`

**Purpose**: Convenient SNP chip definition with multiple selection methods

**Selection Methods:**

1. **Random selection:**
```r
pop %>% define_chip(name = "50k", n_snp = 500, method = "random")
```

2. **Even spacing:**
```r
pop %>% define_chip(name = "HD", n_snp = 50000, method = "even")
```

3. **Chromosome-proportional:**
```r
pop %>% define_chip(name = "10k", n_snp = 10000, method = "chromosome_even")
```

4. **By locus IDs:**
```r
pop %>% define_chip(name = "custom", locus_ids = c(1, 10, 50, 100))
```

5. **By locus names:**
```r
chip_manifest <- c("Locus_1", "Locus_10", "Locus_50")
pop %>% define_chip(name = "custom", locus_names = chip_manifest)
```

6. **Custom column name:**
```r
pop %>% define_chip(name = "bovine_50k", n_snp = 500, col_name = "SNP_50k")
```

**Default column naming**: `is_{name}` (e.g., "50k" → "is_50k")

## Complete Workflow Example

```r
library(tidybreed)
library(dplyr)

# 1. Initialize genome
pop <- initialize_genome(
  pop_name = "test",
  n_loci = 1000,
  n_chr = 10,
  chr_len_Mb = 100,
  n_haplotypes = 100
)

# 2. Define SNP chip
pop <- pop %>%
  define_chip(name = "50k", n_snp = 500, method = "random")

# 3. Add QTL effects to chip SNPs
genome <- get_table(pop, "genome_meta") %>% collect()
n_loci <- nrow(genome)

effect_vec <- ifelse(
  genome$is_50k,                    # if on chip
  rnorm(n_loci, mean = 0, sd = 1),  # assign random effect
  0                                  # otherwise zero
)

pop <- pop %>% mutate_genome_meta(
  effect_50k = effect_vec,
  is_QTL_growth = genome$is_50k
)

# 4. Add founders
pop <- pop %>%
  add_founders(n_males = 10, n_females = 100, line_name = "A")

# 5. Mark individuals as genotyped
pop <- pop %>%
  mutate_ind_meta(genotyped_50k = TRUE, gen = 0)

# 6. Query chip genotypes
chip_loci <- genome %>%
  filter(is_50k == TRUE) %>%
  pull(locus_name)

chip_genotypes <- get_table(pop, "genome_genotype") %>%
  select(ind_id, all_of(chip_loci)) %>%
  collect()

# 7. View results
get_table(pop, "genome_meta") %>%
  filter(is_50k == TRUE) %>%
  select(locus_id, locus_name, chr, pos_Mb, effect_50k, is_QTL_growth) %>%
  head(10) %>%
  collect()
```

## Implementation Notes

### Design Decisions

1. **Separate functions**: `mutate_genome_meta()` for general use, `define_chip()` for convenience
2. **Consistent naming**: Mirrors `mutate_ind_meta()` behavior exactly
3. **Internal helpers**: Chip selection logic is in separate internal functions (not exported)
4. **Flexible column naming**: Default `is_{name}` but customizable with `col_name`

### Helper Functions (Internal)

Located in `R/chip_helpers.R`:
- `select_by_n_snp()` - Random, even, or chromosome-proportional selection
- `select_by_locus_ids()` - Selection by locus ID
- `select_by_locus_names()` - Selection by locus name

These are **not exported** - users should use `define_chip()` instead.

### Error Handling

Both functions include comprehensive error checking:

**`mutate_genome_meta()`:**
- Reserved column modification
- Vector length mismatch
- Invalid field names
- SQL keyword conflicts
- Unsupported data types

**`define_chip()`:**
- No selection method specified
- Multiple selection methods
- Invalid method name
- n_snp exceeds total loci
- Invalid locus IDs/names

### Testing

**25 total test cases:**
- 10 for `mutate_genome_meta()`
- 15 for `define_chip()`
- Includes integration test

All tests use in-memory databases (`:memory:`) for speed.

## Naming Philosophy

This implementation follows the approved naming philosophy:

- **`initialize_*`** - Create core database structure
- **`add_*`** - Add new entities (rows: individuals, progeny)
- **`mutate_*`** - Add/modify attributes (columns in any table)
- **`define_*`** - Helper functions (internally call `mutate_*`)

The old `add_SNP_chip()` design (from DESIGN.md) was **not** implemented, as it would be:
- Inconsistent with `mutate_ind_meta()`
- Redundant with `mutate_genome_meta()`
- Less flexible

## Performance Considerations

- Vector updates use temporary tables and JOIN for efficiency
- Scalar updates use simple UPDATE queries
- Chip selection with large `n_snp` may take a few seconds (chromosome_even method)
- All database operations are transactional

## Future Enhancements

Possible future additions:
1. **Additional spacing methods**:
   - `"chromosome_random"` - Equal count per chromosome, random within
   - `"ld_pruned"` - Prune by linkage disequilibrium
   
2. **Bulk chip definition**:
   - Load chip manifests from CSV/external files
   
3. **Chip querying helpers**:
   - `get_chip_snps(pop, chip_name)` - Get loci for a chip
   - `count_chips(pop)` - Count SNPs per chip

But the current implementation provides a solid, extensible foundation.

# Quick Start Guide: tidybreed

## Testing the Initial Implementation

### Step 1: Document and Load Package

In R, run:

```r
# Generate documentation
devtools::document()

# Load package
devtools::load_all()

# Or install and load
devtools::install()
library(tidybreed)
```

### Step 2: Test Basic Functionality

```r
# Initialize a simple genome
pop <- initialize_genome(
  pop_name = "test_pop",
  n_loci = 100,
  n_chr = 2,
  chr_len_Mb = 50
)

# View the object
print(pop)
# Should show:
# <tidybreed population>
# Population: test_pop
# Database: test_pop_tidybreed.duckdb
# Tables: genome_meta, genome_haplotype, genome_genotype
# Loci: 100
```

### Step 3: Explore the Genome Metadata

```r
# Get genome metadata table
genome <- get_table(pop, "genome_meta")

# Use dplyr verbs (lazy evaluation)
genome %>%
  filter(chr == 1) %>%
  select(locus_name, chr, pos_Mb) %>%
  collect()

# Get summary
genome %>%
  group_by(chr) %>%
  summarise(
    n_loci = n(),
    min_pos = min(pos_Mb),
    max_pos = max(pos_Mb)
  ) %>%
  collect()
```

### Step 4: Test with Different Parameters

```r
# Different chromosome lengths (like cattle genome)
pop_cattle <- initialize_genome(
  pop_name = "cattle",
  n_loci = 50000,
  n_chr = 30,
  chr_len_Mb = c(158, 137, 121, 120, 121, 119, 112, 113, 105, 104,
                 107, 91, 84, 84, 85, 81, 75, 66, 64, 72,
                 71, 61, 52, 62, 42, 51, 45, 46, 51, 42),
  db_path = "cattle_sim.duckdb"
)

print(pop_cattle)

# Custom locus names
pop_custom <- initialize_genome(
  pop_name = "custom",
  n_loci = 10,
  n_chr = 1,
  chr_len_Mb = 100,
  locus_names = paste0("rs", 1:10),
  db_path = ":memory:"  # in-memory for quick tests
)

get_table(pop_custom, "genome_meta") %>% collect()
```

### Step 5: Run Tests

```r
# Run all tests
devtools::test()

# Or run specific test file
devtools::test(filter = "initialize_genome")
```

### Step 6: Clean Up

```r
# Close database connections
close_pop(pop)
close_pop(pop_cattle)
close_pop(pop_custom)

# Remove database files if desired
file.remove("test_pop_tidybreed.duckdb")
file.remove("cattle_sim.duckdb")
```

## What's Working Now

✅ **Core S3 class** (`tidybreed_pop`)
- Wraps DuckDB connection
- Stores population metadata
- Print method for easy viewing

✅ **`initialize_genome()`**
- Creates file-based DuckDB database
- Generates three tables: `genome_meta`, `genome_haplotype`, `genome_genotype`
- Evenly distributes loci across chromosomes
- Supports custom locus names and chromosome lengths

✅ **`get_table()`**
- Returns lazy tibble (efficient querying)
- Works with all dplyr verbs
- No need to load entire table into memory

✅ **`close_pop()`**
- Safely closes database connection

## Next Steps to Implement

1. **`add_SNP_chip()`** - Add SNP chip indicators to genome_meta
2. **`add_founders()`** - Create initial individuals with genotypes/haplotypes
3. **`mutate_ind_meta()`** - Add custom fields to individual metadata
4. **`add_ind_meta()`** - Add values to custom fields

## File Structure

```
tidybreed/
├── R/
│   ├── tidybreed-package.R      # Package documentation
│   ├── tidybreed_pop.R          # S3 class definition
│   └── initialize_genome.R      # Genome initialization
├── tests/
│   └── testthat/
│       ├── testthat.R
│       └── test-initialize_genome.R
├── DESIGN.md                    # Design documentation (READ THIS!)
├── QUICKSTART.md               # This file
└── DESCRIPTION                  # Package metadata
```

## Common Issues

### Issue: Database already exists

**Error:** `Database file 'X_tidybreed.duckdb' already exists`

**Solution:**
```r
# Option 1: Use overwrite = TRUE
pop <- initialize_genome(..., overwrite = TRUE)

# Option 2: Delete file manually
file.remove("X_tidybreed.duckdb")
```

### Issue: Connection not valid

**Error:** `Database connection is no longer valid`

**Solution:** You may have closed the connection. Reinitialize the population.

### Issue: Table not found

**Error:** `Table 'X' not found`

**Solution:** Check available tables:
```r
pop$tables
```

## Design Philosophy

Read `DESIGN.md` for full details, but key points:

1. **Database-first**: All data in DuckDB, not R memory
2. **File-based**: Reduces memory, enables persistence
3. **Lazy evaluation**: Use dplyr without loading full tables
4. **Pipe-friendly**: All functions return modified `tidybreed_pop`
5. **Extensible**: Easy to add custom fields and chips

## Questions or Issues?

Check `DESIGN.md` section "Questions for Discussion" or open an issue on GitHub.

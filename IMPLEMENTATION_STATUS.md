# Implementation Status

**Date:** 2026-04-13
**Phase:** 1 - Core Infrastructure

## ✅ Completed

### 1. S3 Class Structure (`R/tidybreed_pop.R`)

Created the core `tidybreed_pop` S3 class with:

- **`new_tidybreed_pop()`** - Constructor for the population object
- **`validate_tidybreed_pop()`** - Validation function ensuring data integrity
- **`print.tidybreed_pop()`** - Custom print method showing:
  - Population name
  - Database path
  - Active tables
  - Number of individuals (if `ind_meta` exists)
  - Number of loci (if `genome_meta` exists)
- **`get_table()`** - Retrieve lazy tibble from database for dplyr operations
- **`close_pop()`** - Safely close database connection

**Key features:**
- Wraps DuckDB connection
- Stores metadata (n_loci, n_chr, chromosome lengths, etc.)
- Tracks active tables
- Provides clean interface for database access

### 2. Genome Initialization (`R/initialize_genome.R`)

Created **`initialize_genome()`** function that:

**Parameters:**
- `pop_name` - Population identifier
- `n_loci` - Total number of loci
- `n_chr` - Number of chromosomes
- `chr_len_Mb` - Chromosome lengths (single value or vector)
- `db_path` - Database file path (defaults to `{pop_name}_tidybreed.duckdb`)
- `locus_names` - Optional custom locus names
- `chr_names` - Optional custom chromosome names
- `overwrite` - Whether to overwrite existing database

**Creates three tables:**

1. **`genome_meta`** - Locus metadata
   - `locus_id` (INTEGER): 1, 2, 3, ...
   - `locus_name` (VARCHAR): "SNP_1", "SNP_2", ... (or custom)
   - `chr` (INTEGER): Chromosome number
   - `chr_name` (VARCHAR): Chromosome name
   - `pos_Mb` (DOUBLE): Position in megabases

2. **`genome_haplotype`** - Phased haplotypes (currently empty schema)
   - `ind_id` (VARCHAR)
   - `parent_origin` (INTEGER): 1 or 2
   - `locus_1`, `locus_2`, ..., `locus_n` (INTEGER): Alleles

3. **`genome_genotype`** - Genotypes 0/1/2 (currently empty schema)
   - `ind_id` (VARCHAR)
   - `locus_1`, `locus_2`, ..., `locus_n` (INTEGER): Genotypes

**Features:**
- Evenly distributes loci across chromosomes
- Positions loci evenly within each chromosome
- File-based database by default (memory-efficient)
- Prevents accidental overwriting
- Returns `tidybreed_pop` object ready for piping

### 3. Test Suite (`tests/testthat/test-initialize_genome.R`)

Created comprehensive tests:

- ✅ Population object creation
- ✅ Table existence and structure
- ✅ Correct locus assignment to chromosomes
- ✅ Different chromosome lengths
- ✅ Overwrite protection and permission
- ✅ Custom locus names
- ✅ Lazy tibble functionality with dplyr
- ✅ Print method output

### 4. Documentation

Created three documentation files:

1. **`DESIGN.md`** - Comprehensive design document covering:
   - Architecture philosophy
   - Database schema for all tables
   - Table structures and examples
   - Typical workflow examples
   - Implementation priority
   - Design decisions and rationale
   - Future considerations

2. **`QUICKSTART.md`** - Hands-on testing guide with:
   - Step-by-step testing instructions
   - Example code for different scenarios
   - Troubleshooting common issues
   - What's currently working

3. **`IMPLEMENTATION_STATUS.md`** - This file!

## 📋 Next Steps (Priority Order)

### Step 3: `add_SNP_chip()` Function

Create function to add SNP chip indicators to `genome_meta` table.

**API Design:**
```r
# Option 1: Specify number (evenly spaced)
pop %>% add_SNP_chip(chip_name = "50K", n = 50000)

# Option 2: Specify locus IDs
pop %>% add_SNP_chip(chip_name = "HD", loci = c(1, 5, 10, ...))

# Option 3: Logical vector
pop %>% add_SNP_chip(chip_name = "custom", loci = is_on_chip)
```

**Implementation needs:**
- Add column `is_{chip_name}` (BOOLEAN) to `genome_meta`
- Support multiple SNP chips
- Validate locus IDs exist
- Update population metadata

### Step 4: `add_founders()` Function

Create founder individuals with genotypes/haplotypes.

**API Design:**
```r
pop %>% add_founders(n_males = 10, n_females = 100)
```

**Implementation needs:**
- Create `ind_meta` table with required fields:
  - `pop_id`, `ind_id`, `sire_id`, `dam_id`
- Generate random haplotypes for founders
- Compute genotypes from haplotypes
- Populate `genome_haplotype` and `genome_genotype` tables
- ID naming convention (e.g., "A_F0_001")

### Step 5: Individual Metadata Functions

**`mutate_ind_meta()`** - Add new fields:
```r
pop %>% mutate_ind_meta(gen = 0, farm = "A", date_birth = Sys.Date())
```

**`add_ind_meta()`** - Set field values:
```r
# Single value for all
pop %>% add_ind_meta(field_name = "gen", value = 0)

# Vector of values
pop %>% add_ind_meta(field_name = "farm", values = c("A", "B", "C", ...))
```

**Implementation needs:**
- ALTER TABLE to add columns to `ind_meta`
- Type detection or specification
- UPDATE queries to set values
- Validation of vector lengths

## 🎯 Current State

The package now has:
- ✅ Working S3 class
- ✅ Database initialization
- ✅ Genome structure creation
- ✅ dplyr integration via lazy tibbles
- ✅ Comprehensive tests
- ✅ Design documentation

**You can now:**
1. Create populations with genomes
2. Query genome metadata with dplyr
3. Persist data to disk
4. Customize chromosome structures

**Still need to add:**
1. SNP chip annotations
2. Founder individuals
3. Custom metadata fields
4. Trait architecture
5. Phenotype generation
6. Mating and recombination

## 📦 Files Created

```
R/
├── tidybreed_pop.R           # S3 class (177 lines)
└── initialize_genome.R       # Genome initialization (185 lines)

tests/testthat/
└── test-initialize_genome.R  # Test suite (145 lines)

Documentation/
├── DESIGN.md                 # Design document (490 lines)
├── QUICKSTART.md            # Quick start guide (240 lines)
└── IMPLEMENTATION_STATUS.md  # This file

Modified/
└── R/tidybreed-package.R    # Updated imports
```

## 🧪 Testing Instructions

See `QUICKSTART.md` for detailed testing steps, but quick version:

```r
# In R console
devtools::load_all()

# Create test population
pop <- initialize_genome(
  pop_name = "test",
  n_loci = 100,
  n_chr = 2,
  chr_len_Mb = 50
)

# View it
print(pop)

# Query genome
get_table(pop, "genome_meta") %>%
  filter(chr == 1) %>%
  collect()

# Run tests
devtools::test()
```

## 💭 Design Questions to Consider

Before implementing next steps, we should decide:

1. **Founder haplotype generation:**
   - Random allele frequencies? (e.g., MAF = 0.5)
   - User-specified frequencies?
   - Read from real data?

2. **Individual ID format:**
   - `{pop_name}_{gen}_{seq}` (e.g., "A_0_001")?
   - User-specified prefix?
   - UUID for guaranteed uniqueness?

3. **SNP chip storage:**
   - Boolean columns `is_{chip_name}` in `genome_meta`? ✓ (current plan)
   - Separate junction table `chip_membership`?

4. **Missing data handling:**
   - NA for missing?
   - -9 (PLINK convention)?
   - NULL in database?

5. **Sex encoding:**
   - "M"/"F"?
   - "Male"/"Female"?
   - 1/2?
   - User configurable?

## 📝 Notes

- All table names use underscores (not hyphens) for SQL compatibility
- File-based database is default for memory efficiency
- Genotypes stored as 0/1/2 (can convert to -1/0/1 on export)
- Haplotypes store parent_origin as 1 (paternal) or 2 (maternal)
- Both genotypes AND haplotypes stored (disk cheap, compute expensive)

## 🚀 Ready for Review

The core infrastructure is complete and ready for testing. Please:

1. Review `DESIGN.md` for full architecture details
2. Test basic functionality using `QUICKSTART.md`
3. Run test suite
4. Provide feedback on design decisions
5. Confirm approach for next steps (add_SNP_chip, add_founders)

Once you confirm this foundation looks good, we can proceed with implementing the SNP chip and founder functions!

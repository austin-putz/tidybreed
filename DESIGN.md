# tidybreed Package Design

## Overview

`tidybreed` is an R package for breeding program simulation that uses:
- **Tidyverse-style API** with pipe operators
- **DuckDB backend** for memory-efficient storage
- **File-based database** to handle large simulations
- **Flexible, extensible architecture** for custom workflows

## Core Object: `tidybreed_pop`

The central object in tidybreed is the `tidybreed_pop` S3 class, which wraps a DuckDB connection and manages all simulation data.

### Object Structure

```r
tidybreed_pop <- list(
  db_conn = <DuckDB connection>,
  pop_name = "population_name",
  db_path = "path/to/database.duckdb",
  tables = c("genome_meta", "genome_haplotype", ...),
  metadata = list(n_loci, n_chr, chr_len_Mb, ...)
)
```

### Design Principles

1. **Database-first**: All data lives in DuckDB, not in memory
2. **Lazy evaluation**: Use `dplyr::tbl()` for efficient querying
3. **File-based by default**: Reduces memory footprint for large simulations
4. **Pipe-friendly**: All functions return the modified `tidybreed_pop` object
5. **Type-safe**: Explicit types for all table columns

## Database Schema

### Table Naming Convention

All tables use underscores (not hyphens) to avoid SQL issues:
- `genome_meta` ✓
- `genome_haplotype` ✓
- `genome_genotype` ✓
- `ind_meta` ✓
- `ind_phenotype` ✓
- `trait_meta` ✓

### Core Tables

#### 1. `genome_meta`

Metadata for all genomic loci.

**Required columns:**
```
- locus_id (INTEGER): Unique locus identifier (1, 2, 3, ...)
- locus_name (VARCHAR): Locus name (e.g., "Locus_1", "rs12345")
- chr (INTEGER): Chromosome number (1, 2, ..., n_chr)
- chr_name (VARCHAR): Chromosome name (for compatibility with different naming)
- pos_Mb (DOUBLE): Position in megabases
```

**User-added columns** (via `add_SNP_chip()` or similar):
```
- is_50K (BOOLEAN): TRUE if locus on 50K SNP chip
- is_HD (BOOLEAN): TRUE if locus on HD SNP chip
- is_QTL_trait1 (BOOLEAN): TRUE if locus is QTL for trait1
```

**Example:**
| locus_id | locus_name | chr | chr_name | pos_Mb | is_50K | is_HD |
|----------|------------|-----|----------|--------|--------|-------|
| 1        | Locus_1    | 1   | 1        | 0.5    | TRUE   | TRUE  |
| 2        | Locus_2    | 1   | 1        | 1.0    | TRUE   | TRUE  |
| 3        | Locus_3    | 1   | 1        | 1.5    | FALSE  | TRUE  |

#### 2. `genome_haplotype`

Phased haplotypes for each individual. Each individual has 2 rows (one per parent).

**Structure:**
```
- ind_id (VARCHAR): Individual ID
- parent_origin (INTEGER): 1 = paternal, 2 = maternal
- locus_1 (INTEGER): Allele at locus 1 (0 or 1)
- locus_2 (INTEGER): Allele at locus 2 (0 or 1)
- ...
- locus_n (INTEGER): Allele at locus n (0 or 1)
```

**Example (3 loci):**
| ind_id | parent_origin | locus_1 | locus_2 | locus_3 |
|--------|---------------|---------|---------|---------|
| A_001  | 1             | 0       | 1       | 0       |
| A_001  | 2             | 1       | 1       | 1       |
| A_002  | 1             | 0       | 0       | 1       |
| A_002  | 2             | 0       | 1       | 0       |

#### 3. `genome_genotype`

Genotypes (0/1/2 encoding) for each individual.

**Structure:**
```
- ind_id (VARCHAR): Individual ID
- locus_1 (INTEGER): Genotype at locus 1 (0, 1, or 2)
- locus_2 (INTEGER): Genotype at locus 2 (0, 1, or 2)
- ...
- locus_n (INTEGER): Genotype at locus n (0, 1, or 2)
```

**Example (3 loci):**
| ind_id | locus_1 | locus_2 | locus_3 |
|--------|---------|---------|---------|
| A_001  | 1       | 2       | 1       |
| A_002  | 0       | 1       | 1       |

**Note:** Genotypes are sum of haplotypes (0+1=1, 1+1=2, etc.)

#### 4. `ind_meta`

Individual metadata. Created when founders are added.

**Required columns:**
```
- pop_id (VARCHAR): Population ID
- ind_id (VARCHAR): Individual ID (primary key)
- sire_id (VARCHAR): Sire ID (parent 1)
- dam_id (VARCHAR): Dam ID (parent 2)
```

**User-added columns** (via `mutate_ind_meta()` or similar):
```
- gen (INTEGER): Generation number
- sex (VARCHAR): "M" or "F"
- farm (VARCHAR): Farm identifier
- date_birth (DATE): Birth date
- ... (any user-defined field)
```

#### 5. `ind_phenotype`

Phenotypic records in long format.

**Structure:**
```
- ind_id (VARCHAR): Individual ID
- trait (VARCHAR): Trait name
- value (DOUBLE): Phenotypic value
```

**Optional columns:**
```
- env (VARCHAR): Environment identifier
- rep (INTEGER): Replicate number
- date_measured (DATE): Measurement date
```

**Example:**
| ind_id | trait | value | env |
|--------|-------|-------|-----|
| A_001  | ADG   | 1.45  | F1  |
| A_001  | BF    | 12.3  | F1  |
| A_002  | ADG   | 1.38  | F1  |

#### 6. `trait_meta`

Trait architecture and genetic effects.

**Structure:**
```
- trait_name (VARCHAR): Trait identifier
- effect_name (VARCHAR): Effect identifier (e.g., "QTL_1", "locus_50")
- effect_type (VARCHAR): Type ("QTL", "polygenic", "environmental")
- level (VARCHAR): Level ("additive", "dominance", "epistasis")
- locus_id (INTEGER): Associated locus (if applicable)
- true_effect (DOUBLE): True genetic effect size
- distribution (VARCHAR): Distribution type ("normal", "uniform", etc.)
- mean (DOUBLE): Mean of distribution
- variance (DOUBLE): Variance of distribution
```

**Example:**
| trait_name | effect_name | effect_type | level    | locus_id | true_effect |
|------------|-------------|-------------|----------|----------|-------------|
| ADG        | QTL_1       | QTL         | additive | 50       | 0.25        |
| ADG        | QTL_2       | QTL         | additive | 150      | -0.18       |
| ADG        | residual    | environmental | NA     | NA       | 0.00        |

## Typical Workflow

### 1. Initialize Population

```r
library(tidybreed)

# Create population with genome
pop_A <- initialize_genome(
  pop_name = "A",
  n_loci = 1000,
  n_chr = 10,
  chr_len_Mb = 100
)

# Check structure
print(pop_A)
# <tidybreed population>
# Population: A
# Database: A_tidybreed.duckdb
# Tables: genome_meta, genome_haplotype, genome_genotype
# Loci: 1000
```

### 2. Add SNP Chips

```r
# Add 50K SNP chip (first 500 loci)
pop_A <- pop_A %>%
  add_SNP_chip(
    chip_name = "50K",
    n = 500  # evenly spaced
  )

# Or specify exact loci
pop_A <- pop_A %>%
  add_SNP_chip(
    chip_name = "HD",
    loci = c(1, 5, 10, 15, ...)  # specific loci
  )

# Or use logical vector
is_on_chip <- rep(c(TRUE, FALSE), each = 500)
pop_A <- pop_A %>%
  add_SNP_chip(
    chip_name = "custom",
    loci = is_on_chip
  )
```

### 3. Add Founders

```r
# Create founder animals
pop_A <- pop_A %>%
  add_founders(
    n_males = 10,
    n_females = 100
  )

# Add custom metadata fields
pop_A <- pop_A %>%
  mutate_ind_meta(
    gen = 0,
    farm = "A",
    date_birth = Sys.Date()
  )

# Or add field with values
pop_A <- pop_A %>%
  add_ind_meta(
    field_name = "farm",
    values = sample(c("A", "B", "C"), 110, replace = TRUE)
  )
```

### 4. Define Traits

```r
pop_A <- pop_A %>%
  add_trait(
    name = "ADG",
    n_qtl = 10,
    target_add_var = 1.0,
    res_var = 1.0
  )
```

### 5. Generate Phenotypes

```r
pop_A <- pop_A %>%
  add_phenotype(trait = "ADG")
```

### 6. Query and Filter

```r
# Get specific table
genome <- get_table(pop_A, "genome_meta")

# Use dplyr verbs
genome %>%
  filter(chr == 1, is_50K == TRUE) %>%
  select(locus_name, pos_Mb) %>%
  collect()

# Filter individuals (future functionality)
pop_A %>%
  filter(sex == "M", gen == 0) %>%
  select_animals()
```

## Implementation Priority

### Phase 1: Core Infrastructure (CURRENT)
- [x] S3 class `tidybreed_pop`
- [x] `initialize_genome()`
- [ ] `add_SNP_chip()`
- [ ] `get_table()` and print methods

### Phase 2: Population Management
- [ ] `add_founders()`
- [ ] `mutate_ind_meta()`
- [ ] `add_ind_meta()`
- [ ] Individual ID system

### Phase 3: Traits and Phenotypes
- [ ] `add_trait()`
- [ ] `add_phenotype()`
- [ ] `sample_phenotype()`

### Phase 4: Mating and Selection
- [ ] `mate()`
- [ ] `select_parents()`
- [ ] Recombination engine

### Phase 5: Advanced Features
- [ ] `filter()` method for populations
- [ ] `mutate()` method for populations
- [ ] Export functions (PLINK, VCF, etc.)
- [ ] Visualization functions

## Design Decisions

### Why DuckDB?

1. **Fast**: Columnar storage optimized for analytics
2. **Embedded**: No server setup required
3. **SQL**: Powerful querying via dplyr/dbplyr
4. **Scalable**: Handles datasets larger than RAM
5. **R-friendly**: Excellent integration with tidyverse

### Why File-Based?

- Large simulations (100K+ individuals, 50K+ SNPs) exceed memory limits
- File persistence enables:
  - Resumable simulations
  - Sharing results
  - Long-term storage
  - Parallel processing across sessions

### Why 0/1/2 Encoding?

- Standard in genomics (PLINK, VCF)
- Easy interpretation: count of minor/reference allele
- Efficient storage as integers
- Can convert to -1/0/1 or 0/0.5/1 for specific analyses

### Genotype vs Haplotype Storage

**Haplotypes:**
- 2 rows per individual (paternal/maternal)
- Required for: recombination, phased data export, parent-of-origin effects
- Storage: ~2x individuals × n_loci integers

**Genotypes:**
- 1 row per individual
- Required for: GWAS, genomic prediction, most analyses
- Storage: ~individuals × n_loci integers
- Can be computed from haplotypes, but stored for speed

**Decision:** Store both. Disk is cheap, recomputation is expensive.

## Future Considerations

1. **Sparse genotype storage**: Most simulations start with founders and fill in over time
2. **Compression**: DuckDB supports compression; may want to experiment
3. **Indexing**: Add indexes on frequently queried columns (ind_id, chr, etc.)
4. **Partitioning**: Consider partitioning by chromosome for very large datasets
5. **Parallelization**: DuckDB supports parallel queries; explore for mating/recombination

## Questions for Discussion

1. Should we store genotypes AND haplotypes, or compute genotypes on-the-fly?
2. How to handle missing data (NA vs -9 vs special encoding)?
3. Should SNP chips be stored as columns or in a separate junction table?
4. What's the best way to version/track simulation parameters?
5. How to handle pedigree structure? Separate table or part of ind_meta?

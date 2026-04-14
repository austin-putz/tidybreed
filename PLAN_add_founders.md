# Plan: Implement add_founders() Function

## Context

This implementation adds the `add_founders()` function to create founder individuals by sampling haplotypes from the `founder_haplotypes` table. This is the critical second step in the tidybreed workflow that connects genome initialization with population management.

**Why this change is needed:**
- Users need to create initial founder populations to start breeding simulations
- Founders are created by sampling from realistic haplotype pools (created by initialize_genome)
- Multiple populations can be added to the same database for multi-line selection studies
- This function bridges genome structure (initialize_genome) with individual management (mutate_ind_meta)

**Intended workflow:**
```r
pop <- initialize_genome(pop_name = "test", n_loci = 1000, n_chr = 10, n_haplotypes = 100) %>%
  add_founders(n_males = 10, n_females = 100, pop_name = "A") %>%
  mutate_ind_meta(gen = 0, farm = "A")
```

---

## Implementation Plan

### 1. Function Signature

```r
add_founders <- function(pop, n_males, n_females, pop_name)
```

**Parameters:**
- `pop` - tidybreed_pop object
- `n_males` - Number of male founders
- `n_females` - Number of female founders
- `pop_name` - Population identifier for ID generation

### 2. Core Algorithm

1. **Validate inputs**
   - Check pop is valid tidybreed_pop object
   - Validate n_males, n_females are positive integers (allow 0 but at least one must be > 0)
   - Validate pop_name format (letters, numbers, hyphens, underscores)
   - Check founder_haplotypes table exists

2. **Read founder haplotypes**
   - Load entire founder_haplotypes table into memory (once)
   - Extract locus columns and convert to matrix for efficiency
   - Validate table has data

3. **Determine ID sequence**
   - Check if ind_meta exists
   - Query max ID number for this population if exists
   - Start from 1 for new population, or continue from max+1

4. **Sample haplotypes**
   - For each founder: sample 2 haplotype indices with replacement
   - Use matrix indexing for efficient extraction
   - Store as n_founders × 2 matrix of indices

5. **Create ind_meta data**
   - Generate IDs: "{pop_name}-{seq}" (e.g., "A-1", "A-2")
   - Set parent_1 and parent_2 to NA (NULL in database)
   - Set population to pop_name
   - Set sex based on n_males/n_females

6. **Create genome_haplotype data**
   - 2 rows per individual (parent_origin 1 and 2)
   - Extract sampled haplotypes from matrix
   - Build data frame with ind_id, parent_origin, locus_1, locus_2, ...

7. **Create genome_genotype data**
   - 1 row per individual
   - Genotype = sum of two haplotypes (element-wise)
   - Build data frame with ind_id, locus_1, locus_2, ...

8. **Write to database**
   - ind_meta: CREATE if new, INSERT if exists
   - genome_haplotype: Append (table exists from initialize_genome)
   - genome_genotype: Append (table exists from initialize_genome)

9. **Update pop object**
   - Add "ind_meta" to pop$tables if new
   - Update pop$metadata$n_individuals
   - Return modified pop (invisibly)

### 3. Key SQL Patterns

**Find max ID for population:**
```sql
SELECT MAX(CAST(SUBSTRING(ind_id FROM POSITION('-' IN ind_id) + 1) AS INTEGER)) as max_num
FROM ind_meta
WHERE ind_id LIKE '{pop_name}-%'
```

**Data writes use DBI::dbWriteTable:**
```r
DBI::dbWriteTable(conn, "ind_meta", data_df, append = TRUE/FALSE)
```

### 4. Data Structures

**ind_meta:**
```
ind_id (VARCHAR): "A-1", "A-2", ..., "B-1"
parent_1 (VARCHAR): NULL for founders
parent_2 (VARCHAR): NULL for founders
population (VARCHAR): "A", "B", etc.
sex (VARCHAR): "M" or "F"
```

**genome_haplotype (2 rows per individual):**
```
ind_id (VARCHAR): Individual ID
parent_origin (INTEGER): 1 (paternal) or 2 (maternal)
locus_1 (INTEGER): 0 or 1
...
locus_n (INTEGER): 0 or 1
```

**genome_genotype (1 row per individual):**
```
ind_id (VARCHAR): Individual ID
locus_1 (INTEGER): 0, 1, or 2 (sum of haplotypes)
...
locus_n (INTEGER): 0, 1, or 2
```

### 5. Detailed Implementation Steps

#### 5.1 Input Validation
```r
# Validate pop object
stopifnot(inherits(pop, "tidybreed_pop"))
validate_tidybreed_pop(pop)

# Validate n_males and n_females
stopifnot(is.numeric(n_males), length(n_males) == 1, n_males >= 0)
stopifnot(is.numeric(n_females), length(n_females) == 1, n_females >= 0)
stopifnot((n_males + n_females) > 0)  # At least one founder

# Validate pop_name
stopifnot(is.character(pop_name), length(pop_name) == 1, nchar(pop_name) > 0)

# Validate field name format for pop_name
if (!grepl("^[a-zA-Z][a-zA-Z0-9_-]*$", pop_name)) {
  stop("pop_name must start with letter and contain only letters, numbers, underscores, or hyphens",
       call. = FALSE)
}

# Check founder_haplotypes table exists
if (!"founder_haplotypes" %in% pop$tables) {
  stop(
    "founder_haplotypes table does not exist. ",
    "Call initialize_genome() with n_haplotypes parameter to create founder haplotypes.",
    call. = FALSE
  )
}
```

#### 5.2 Read Founder Haplotypes
```r
# Read all founder haplotypes from database (once, into memory)
founder_haps_tbl <- get_table(pop, "founder_haplotypes") %>% 
  dplyr::collect()

# Get number of available haplotypes
n_haplotypes <- nrow(founder_haps_tbl)

if (n_haplotypes == 0) {
  stop("founder_haplotypes table is empty. Cannot sample haplotypes.", call. = FALSE)
}

# Get number of loci from founder_haplotypes columns
locus_cols <- grep("^locus_", colnames(founder_haps_tbl), value = TRUE)
n_loci <- length(locus_cols)

if (n_loci == 0) {
  stop("No locus columns found in founder_haplotypes table.", call. = FALSE)
}

# Convert to matrix for efficient access
hap_data_matrix <- as.matrix(founder_haps_tbl[, locus_cols])
```

#### 5.3 Determine ID Sequence
```r
# Check if ind_meta already exists
ind_meta_exists <- "ind_meta" %in% DBI::dbListTables(pop$db_conn)

# Determine starting ID number for this population
if (ind_meta_exists) {
  # Query max ID number for this population
  max_num_query <- paste0(
    "SELECT MAX(CAST(SUBSTRING(ind_id FROM POSITION('-' IN ind_id) + 1) AS INTEGER)) as max_num ",
    "FROM ind_meta ",
    "WHERE ind_id LIKE '", pop_name, "-%'"
  )
  
  max_num_result <- DBI::dbGetQuery(pop$db_conn, max_num_query)
  
  if (is.na(max_num_result$max_num)) {
    start_id <- 1
  } else {
    start_id <- max_num_result$max_num + 1
  }
} else {
  start_id <- 1
}
```

#### 5.4 Sample Haplotypes
```r
# Total number of founders
n_founders <- n_males + n_females

# Sample haplotypes: 2 per individual, with replacement
hap_indices <- matrix(
  sample(1:n_haplotypes, size = n_founders * 2, replace = TRUE),
  nrow = n_founders,
  ncol = 2
)
```

#### 5.5 Create ind_meta Data
```r
# Generate individual IDs
ind_ids <- paste0(pop_name, "-", seq(start_id, start_id + n_founders - 1))

# Create sex vector
sex_vector <- c(rep("M", n_males), rep("F", n_females))

# Create ind_meta data frame
ind_meta_df <- tibble::tibble(
  ind_id = ind_ids,
  parent_1 = NA_character_,  # NULL for founders
  parent_2 = NA_character_,  # NULL for founders
  population = pop_name,
  sex = sex_vector
)
```

#### 5.6 Create genome_haplotype Data
```r
# Build haplotype data frame (2 rows per individual)
hap_list <- vector("list", n_founders)

for (i in 1:n_founders) {
  # Get the two haplotype indices for this individual
  hap_idx_1 <- hap_indices[i, 1]
  hap_idx_2 <- hap_indices[i, 2]
  
  # Create 2-row data frame for this individual
  hap_list[[i]] <- tibble::tibble(
    ind_id = rep(ind_ids[i], 2),
    parent_origin = c(1L, 2L)
  )
  
  # Add locus columns
  for (j in 1:n_loci) {
    locus_name <- paste0("locus_", j)
    hap_list[[i]][[locus_name]] <- c(
      hap_data_matrix[hap_idx_1, j],
      hap_data_matrix[hap_idx_2, j]
    )
  }
}

# Combine all individuals into single data frame
genome_haplotype_df <- dplyr::bind_rows(hap_list)
```

#### 5.7 Create genome_genotype Data
```r
# Build genotype data frame (1 row per individual)
geno_df <- tibble::tibble(ind_id = ind_ids)

for (j in 1:n_loci) {
  locus_name <- paste0("locus_", j)
  
  # Extract the two haplotypes for this locus across all individuals
  hap1_vals <- hap_data_matrix[hap_indices[, 1], j]
  hap2_vals <- hap_data_matrix[hap_indices[, 2], j]
  
  # Sum to get genotype (0, 1, or 2)
  geno_df[[locus_name]] <- as.integer(hap1_vals + hap2_vals)
}

genome_genotype_df <- geno_df
```

#### 5.8 Write Data to Database
```r
# Write ind_meta table
if (ind_meta_exists) {
  DBI::dbWriteTable(pop$db_conn, "ind_meta", ind_meta_df, append = TRUE)
} else {
  DBI::dbWriteTable(pop$db_conn, "ind_meta", ind_meta_df, overwrite = FALSE)
  pop$tables <- c(pop$tables, "ind_meta")
}

# Append to genome_haplotype
DBI::dbWriteTable(pop$db_conn, "genome_haplotype", genome_haplotype_df, append = TRUE)

# Append to genome_genotype
DBI::dbWriteTable(pop$db_conn, "genome_genotype", genome_genotype_df, append = TRUE)
```

#### 5.9 Update and Return
```r
# Update metadata
if (is.null(pop$metadata$n_individuals)) {
  pop$metadata$n_individuals <- n_founders
} else {
  pop$metadata$n_individuals <- pop$metadata$n_individuals + n_founders
}

# Return modified pop object
invisible(pop)
```

### 6. Error Handling

**Error conditions:**
- Pop object invalid
- n_males/n_females not positive integers
- n_males + n_females == 0
- pop_name invalid format
- founder_haplotypes table doesn't exist
- founder_haplotypes table is empty
- No locus columns found

**Clear error messages:**
```r
"founder_haplotypes table does not exist. Call initialize_genome() with n_haplotypes parameter."
"pop_name must start with letter and contain only letters, numbers, underscores, or hyphens"
"At least one founder must be specified (n_males + n_females > 0)"
"founder_haplotypes table is empty. Cannot sample haplotypes."
"No locus columns found in founder_haplotypes table."
```

### 7. Testing Strategy

**Test file:** `tests/testthat/test-add_founders.R`

**Test categories:**

1. **Basic functionality**
   - Add founders to fresh database
   - Verify ind_meta created with correct columns
   - Verify genome_haplotype populated (2 rows per individual)
   - Verify genome_genotype populated (1 row per individual)
   - Verify genotypes = sum of haplotypes

2. **Multiple populations**
   - Add population A (n=10)
   - Add population B (n=10)
   - Verify IDs: "A-1", "A-2", ..., "B-1", "B-2"
   - Verify ind_meta has all individuals

3. **Sequential additions**
   - Add population A (n=10)
   - Add more population A (n=5)
   - Verify IDs continue: "A-11" through "A-15"

4. **Sex assignment**
   - Verify first n_males have sex="M"
   - Verify remaining have sex="F"

5. **Parent IDs**
   - Verify parent_1 and parent_2 are NULL (NA)

6. **Haplotype sampling**
   - Verify 2 haplotypes per individual
   - Verify sampling with replacement

7. **Error conditions**
   - Error if founder_haplotypes doesn't exist
   - Error if pop not valid tidybreed_pop
   - Error if n_males/n_females invalid
   - Error if pop_name invalid
   - Error if founder_haplotypes empty

8. **Integration**
   - Full pipeline: initialize_genome → add_founders → mutate_ind_meta

9. **Edge cases**
   - n_males = 0, n_females > 0
   - n_males > 0, n_females = 0
   - Large populations (1000+ founders)

### 8. Documentation

**Roxygen2 documentation:**
```r
#' Add founder individuals to population
#'
#' @description
#' Creates founder individuals by sampling haplotypes from the \code{founder_haplotypes}
#' table. Each founder receives two randomly sampled haplotypes (with replacement),
#' which are used to populate the \code{genome_haplotype} and \code{genome_genotype} tables.
#' 
#' The \code{ind_meta} table is created (if it doesn't exist) or appended to with
#' the new founders. Founder individuals have \code{NULL} for both parent IDs.
#'
#' @param pop A \code{tidybreed_pop} object
#' @param n_males Integer. Number of male founders to create
#' @param n_females Integer. Number of female founders to create
#' @param pop_name Character. Population identifier used for individual IDs.
#'   IDs are formatted as \code{"{pop_name}-{number}"} (e.g., "A-1", "A-2")
#'
#' @return The \code{tidybreed_pop} object (invisibly), modified in place
#'
#' @details
#' **Requirements:**
#' - The \code{founder_haplotypes} table must exist. Create it by calling
#'   \code{initialize_genome()} with the \code{n_haplotypes} parameter.
#'
#' **What it does:**
#' 1. Samples 2 haplotypes per founder from \code{founder_haplotypes} (with replacement)
#' 2. Creates/updates \code{ind_meta} table with founder metadata
#' 3. Populates \code{genome_haplotype} (2 rows per individual)
#' 4. Populates \code{genome_genotype} (1 row per individual, sum of haplotypes)
#'
#' **ID Format:**
#' - Individual IDs: \code{"{pop_name}-{number}"} (e.g., "A-1", "A-2", "B-1")
#' - Numbers are sequential within each population
#' - If founders already exist for a population, numbering continues from max ID
#'
#' **Multiple Populations:**
#' - Can be called multiple times to add different populations to same database
#' - Each population has independent ID numbering
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Initialize genome with founder haplotypes
#' pop <- initialize_genome(
#'   pop_name = "test",
#'   n_loci = 1000,
#'   n_chr = 10,
#'   chr_len_Mb = 100,
#'   n_haplotypes = 100
#' )
#'
#' # Add founders
#' pop <- pop %>%
#'   add_founders(n_males = 10, n_females = 100, pop_name = "A")
#'
#' # Add custom metadata
#' pop <- pop %>%
#'   mutate_ind_meta(
#'     gen = 0,
#'     farm = "FarmA",
#'     date_birth = Sys.Date()
#'   )
#'
#' # Add second population to same database
#' pop <- pop %>%
#'   add_founders(n_males = 5, n_females = 50, pop_name = "B")
#'
#' # View founders
#' get_table(pop, "ind_meta") %>% collect()
#' }
```

### 9. Critical Files

**New files:**
1. `R/add_founders.R` - Main implementation
2. `tests/testthat/test-add_founders.R` - Test suite
3. `man/add_founders.Rd` - Generated documentation

**Files to modify:**
4. `NAMESPACE` - Export add_founders
5. `DESIGN.md` - Update documentation
6. `IMPLEMENTATION_STATUS.md` - Mark complete

**Reference files:**
7. `R/tidybreed_pop.R` - S3 patterns
8. `R/initialize_genome.R` - Table creation patterns
9. `R/mutate_ind_meta.R` - Validation patterns

---

## Verification Plan

**Manual testing:**
```r
# Create population
pop <- initialize_genome(
  pop_name = "test",
  n_loci = 100,
  n_chr = 5,
  chr_len_Mb = 100,
  n_haplotypes = 50
) %>%
  add_founders(n_males = 10, n_females = 100, pop_name = "A")

# Verify ind_meta
ind_meta <- get_table(pop, "ind_meta") %>% collect()
print(head(ind_meta))
# Expected: A-1 through A-110, NULL parents, sex M/F

# Verify haplotypes (2 per individual)
haps <- get_table(pop, "genome_haplotype") %>% collect()
print(nrow(haps))  # Expected: 220 (110 × 2)

# Verify genotypes (1 per individual)
genos <- get_table(pop, "genome_genotype") %>% collect()
print(nrow(genos))  # Expected: 110

# Verify genotype calculation
for (i in 1:5) {
  ind_id <- ind_meta$ind_id[i]
  hap_rows <- haps %>% filter(ind_id == !!ind_id)
  geno_row <- genos %>% filter(ind_id == !!ind_id)
  
  # Check locus_1
  expect_equal(
    geno_row$locus_1,
    sum(hap_rows$locus_1)
  )
}
```

**Automated testing:**
1. Run `devtools::test()` - all tests pass
2. Test coverage for add_founders.R
3. Integration tests with complete workflow

---

## Performance Considerations

- Read founder_haplotypes once into memory
- Use matrix operations for haplotype extraction
- Build complete data frames in R, then batch write to database
- Avoid row-by-row database operations
- Pre-allocate lists for efficiency

---

## Next Steps After Implementation

1. Test with realistic scenarios (large populations, multiple populations)
2. Implement breeding functions (mate, select)
3. Add recombination logic for non-founder individuals
4. Consider adding founder import from external files

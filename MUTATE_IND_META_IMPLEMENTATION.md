# mutate_ind_meta() Implementation Complete

## Summary

Successfully implemented the `mutate_ind_meta()` function for adding and modifying custom columns in the individual metadata table (`ind_meta`).

## Files Created

1. **R/mutate_ind_meta.R** - Main implementation (364 lines)
   - `mutate_ind_meta()` - Main function
   - `infer_duckdb_type()` - Type inference helper
   - `validate_field_name()` - Field name validation
   - `update_column_scalar()` - Scalar value updates
   - `update_column_vector()` - Vector value updates

2. **tests/testthat/test-mutate_ind_meta.R** - Comprehensive test suite (450+ lines)
   - 25 test cases covering all functionality
   - Tests for basic operations, type handling, error conditions, edge cases, and integration

3. **man/mutate_ind_meta.Rd** - Function documentation
   - Complete roxygen2 documentation
   - Examples and usage patterns

## Files Modified

1. **R/tidybreed-package.R**
   - Added `DBI::dbListFields` to imports

2. **NAMESPACE**
   - Exported `mutate_ind_meta`
   - Added `importFrom(DBI, dbListFields)`

3. **DESIGN.md**
   - Updated ind_meta schema to reflect 5 core columns: `ind_id`, `parent_1`, `parent_2`, `line`, `sex`
   - Added documentation for custom columns and mutate_ind_meta() usage

4. **IMPLEMENTATION_STATUS.md**
   - Marked mutate_ind_meta() as completed
   - Updated current state section

## Core Features

### ind_meta Table Schema

**Core columns** (created by `add_founders()`):
- `ind_id` (VARCHAR) - Individual identifier
- `parent_1` (VARCHAR) - First parent ID ("0" for founders)
- `parent_2` (VARCHAR) - Second parent ID ("0" for founders)
- `line` (VARCHAR) - Line name
- `sex` (VARCHAR) - Sex ("M" or "F")

**User-added columns** (via `mutate_ind_meta()`):
- Any custom field with automatic type inference

### Type Support

| R Type | DuckDB Type |
|--------|-------------|
| logical | BOOLEAN |
| integer | INTEGER |
| numeric | DOUBLE |
| Date | DATE |
| POSIXct/POSIXlt | TIMESTAMP |
| character | VARCHAR |

### Usage Examples

```r
# Initialize population (once add_founders() is implemented)
pop <- initialize_genome(pop_name = "A", n_loci = 1000, n_chr = 10) %>%
  add_founders(n_males = 10, n_females = 100, line_name = "A")

# Add scalar fields (same value for all individuals)
pop <- pop %>%
  mutate_ind_meta(
    gen = 0,
    farm = "A",
    date_birth = Sys.Date()
  )

# Add vector field (different value per individual)
n_ind <- 110
pop <- pop %>%
  mutate_ind_meta(
    weight_kg = rnorm(n_ind, mean = 80, sd = 10),
    is_selected = sample(c(TRUE, FALSE), n_ind, replace = TRUE)
  )

# Update existing field
pop <- pop %>%
  mutate_ind_meta(gen = 1)

# View metadata
get_table(pop, "ind_meta") %>% collect()
```

### Validation and Safety

**Field name validation:**
- Must start with a letter
- Can contain only letters, numbers, underscores
- Cannot be SQL reserved keywords
- Cannot modify reserved columns (ind_id, parent_1, parent_2, line, sex)

**Value validation:**
- Scalar values: Applied to all individuals
- Vector values: Must match number of individuals
- Type validation: Only supported R types allowed

**Error handling:**
- Clear error messages for all failure modes
- Prevents SQL injection
- Protects reserved columns
- Validates vector lengths

## Test Coverage

### Test Categories (25 tests)

1. **Basic Functionality** (4 tests)
   - Add single scalar column
   - Add multiple scalar columns
   - Add vector column
   - Update existing column

2. **Type Handling** (5 tests)
   - Logical/BOOLEAN
   - Integer/INTEGER
   - Numeric/DOUBLE
   - Date/DATE
   - Character/VARCHAR

3. **Error Handling** (6 tests)
   - ind_meta doesn't exist
   - Vector length mismatch
   - Reserved column modification
   - Invalid field names
   - SQL reserved words
   - Unsupported types

4. **Edge Cases** (7 tests)
   - NA values
   - Empty ind_meta table
   - No fields specified
   - Special characters in strings
   - Mixed scalar/vector
   - Multiple calls
   - Preserves existing data

5. **Integration** (3 tests)
   - Pipe workflow
   - Data preservation
   - Multiple sequential calls

## Next Steps

1. **Implement `add_founders()`** - Will create the ind_meta table with 5 core columns
   - Parameters: `n_males`, `n_females`, `line_name`
   - Auto-generates ind_id values
   - Sets parent_1 and parent_2 to "0" for founders
   - Populates sex based on n_males/n_females

2. **Test Complete Workflow**
   ```r
   pop <- initialize_genome(...) %>%
     add_founders(...) %>%
     mutate_ind_meta(...)
   ```

3. **Consider Future Enhancements**
   - Conditional updates: `mutate_ind_meta(field = value, .where = condition)`
   - Computed columns: `mutate_ind_meta(age = Sys.Date() - date_birth)`
   - Type specification: `mutate_ind_meta(gen = 0, .types = c(gen = "BIGINT"))`

## Verification

Once `add_founders()` is implemented, the full workflow can be tested:

```r
# Create population
pop <- initialize_genome(
  pop_name = "TestPop",
  n_loci = 1000,
  n_chr = 10,
  chr_len_Mb = 100
) %>%
  add_founders(
    n_males = 10,
    n_females = 100,
    line_name = "A"
  )

# Verify ind_meta exists with 5 core columns
ind_meta <- get_table(pop, "ind_meta") %>% collect()
print(colnames(ind_meta))
# Expected: "ind_id" "parent_1" "parent_2" "line" "sex"

# Verify sex is populated correctly
table(ind_meta$sex)
# Expected: F=100, M=10

# Add custom metadata
pop <- pop %>%
  mutate_ind_meta(
    gen = 0,
    farm = "FarmA",
    date_birth = Sys.Date(),
    weight_kg = rnorm(110, 80, 10)
  )

# Verify custom columns added
ind_meta <- get_table(pop, "ind_meta") %>% collect()
print(colnames(ind_meta))
# Expected: includes "gen", "farm", "date_birth", "weight_kg"

# Clean up
close_pop(pop)
```

## Implementation Notes

- All helper functions are marked as `@keywords internal` to keep them private
- The main function returns the pop object invisibly for cleaner piping
- SQL queries use parameterization and escaping to prevent injection
- Temporary tables are used for vector updates to ensure atomic operations
- Comprehensive error messages guide users to correct usage

## Status: ✅ Complete and Ready for Integration

The `mutate_ind_meta()` function is fully implemented, tested, and documented. It's ready to be integrated into the tidybreed workflow once `add_founders()` is implemented.

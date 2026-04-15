# tidybreed 0.0.2 (2026-04-15)

## Bug Fixes

* `infer_duckdb_type()`: fixed incorrect VARCHAR inference when a typed vector
  starts with `NA` (e.g. `c(NA_real_, rnorm(n))`). Type is now determined from
  the R class of the full vector before inspecting individual elements, so
  `NA_real_` → `DOUBLE`, `NA_integer_` → `INTEGER`, `NA` → `BOOLEAN`. A bare
  all-`NA` vector (class `logical`) now warns with a message pointing users to
  typed NA constants.

# tidybreed 0.0.1 (2026-04-14)

## Breaking Changes

### Terminology: "population" → "line"

* `add_founders()` parameter renamed: `pop_name` → `line_name`
* `ind_meta` table column renamed: `population` → `line`

**Rationale:** "population" refers to one complete genome build from
`initialize_genome()`. "Line" refers to distinct genetic lines within that
population (e.g., line A selected for growth, line B for egg quality).

**Migration for existing databases:**
```r
conn <- DBI::dbConnect(duckdb::duckdb(), "your_database.duckdb")
DBI::dbExecute(conn, "ALTER TABLE ind_meta RENAME COLUMN population TO line")
DBI::dbDisconnect(conn)
```

## New Functions

* `mutate_genome_meta()` — add or update user columns in `genome_meta`;
  reserved columns (`locus_id`, `locus_name`, `chr`, `chr_name`, `pos_Mb`)
  are blocked
* `define_chip()` — convenience wrapper that marks loci as members of a named
  SNP chip by writing a `BOOLEAN` column `is_{chip_name}` to `genome_meta`;
  supports `"random"`, `"even"`, and `"chr_even"` selection methods

## Other Changes

* Removed stale design and planning documents (DESIGN.md, QUICKSTART.md,
  IMPLEMENTATION_STATUS.md, IMPLEMENTATION_SUMMARY.md, TODO_finalize.md,
  MUTATE_IND_META_IMPLEMENTATION.md)
* Moved manual test scripts to `tests/` directory
* Added CLAUDE.md with developer and AI context for the project

---

# tidybreed 0.0.0

* Initial package setup and framework

# tidybreed (development version)

## Breaking Changes (2026-04-14)

### Terminology Clarification: "population" → "line"

**BREAKING:** Renamed field and parameter to clarify terminology:
* `add_founders()` parameter renamed: `pop_name` → `line_name`
* `ind_meta` table column renamed: `population` → `line`

**Rationale:**
* **Population**: One complete genome build from `initialize_genome()` (one genetic system)
* **Line**: Different genetic lines within a population selected differently (e.g., different selection indices or breeding goals)

For example, you might have one population (genome) but create line A (selected for growth) and line B (selected for egg quality) within it.

### Migration Required

**For existing databases:**
```r
library(DBI)
conn <- dbConnect(duckdb::duckdb(), "your_database.duckdb")
dbExecute(conn, "ALTER TABLE ind_meta RENAME COLUMN population TO line")
dbDisconnect(conn)
```

**Update all code:**
```r
# OLD:
pop <- add_founders(pop, n_males = 10, n_females = 100, pop_name = "A")
ind_meta %>% filter(population == "A")
ind_meta$population

# NEW:
pop <- add_founders(pop, n_males = 10, n_females = 100, line_name = "A")
ind_meta %>% filter(line == "A")
ind_meta$line
```

**Updated reserved columns:**
* OLD: `ind_id`, `parent_1`, `parent_2`, `population`, `sex`
* NEW: `ind_id`, `parent_1`, `parent_2`, `line`, `sex`

---

## tidybreed 0.0.0.9000

* Initial package setup
* Added basic package framework

# tidybreed 0.3.0 (2026-04-21)

## New Functions

* `add_genotypes(pop, chip_name)` — marks a filtered subset of animals as
  genotyped on a named SNP chip by writing a `has_<chip_name>` BOOLEAN column
  to `ind_meta`. Follows the same `filter()` -> action pipe pattern as
  `add_phenotype()`. Operation is additive: animals already marked TRUE remain
  TRUE across multiple calls. Chip must exist in `genome_meta` (via
  `define_chip()`) before calling.
* `extract_genotypes(pop, chip_name)` — returns a tibble of genotypes
  (0/1/2 encoding) for animals marked as genotyped, restricted to chip loci.
  The returned set is the intersection of animals with `has_<chip_name> == TRUE`,
  any pending `filter()` predicates, and loci with `is_<chip_name> == TRUE`.
  Intended for use immediately before GBLUP/GWAS evaluation.

---

# tidybreed 0.2.2 (2026-04-20)

## Performance

- `add_founders()`: eliminated nested R loop (O(n_founders × n_loci) iterations)
  that caused multi-minute runtimes for large genomes. Haplotype and genotype
  frames are now built via vectorized matrix indexing and addition in C, reducing
  frame construction from minutes to ~0.18 s regardless of genome size. Also
  switched genome table writes to `duckdb_register` + `INSERT SELECT` for a
  further ~2× speedup on the DB write step. Typical runtimes: 0.38 s (1k loci)
  → 2.7 s (5k loci) → 8 s (10k loci) for 215 founders.

# tidybreed 0.2.1 (2026-04-20)

## API

- Renamed `name` argument to `trait_name` in `add_trait()` and `add_trait_simple()` for consistency
- Renamed `trait` argument to `trait_name` in `define_qtl()` and `set_qtl_effects()`
- Renamed `traits` argument to `trait_names` in `set_qtl_effects_multi()`
- Renamed `name` argument to `chip_name` in `define_chip()`

# tidybreed 0.2.0 (2026-04-20)

## New: Trait Architecture & Phenotype Simulation

Trait definition, QTL selection, additive-effect sampling (single and
correlated multi-trait), fixed and random covariates, phenotype generation
for continuous / count / binary / categorical traits, imprinting support,
and storage of true and estimated breeding values.

### New exported functions

* `add_trait()` — insert a row in `trait_meta` with target variance
  components, trait type (continuous / count / binary / categorical),
  expression rules (sex-limited, parent-of-origin), and index/economic
  weights. Creates the six trait-layer tables on first call.
* `define_qtl()` — mirror of `define_chip()` for QTL loci; writes an
  `is_QTL_{trait}` BOOLEAN column in `genome_meta`. Reuses the same six
  selection methods (by count + `random` / `even` / `chromosome_even`, by
  logical vector, by locus ids, by locus names).
* `set_qtl_effects()` — write the `add_{trait}` additive-effect column.
  Supports manual effects or sampled effects (`normal` / `gamma`) with
  automatic rescaling to hit `trait_meta$target_add_var`.
* `set_qtl_effects_multi()` — draw correlated additive effects across
  multiple traits from `MVN(0, G)` via `MASS::mvrnorm`; supports `"shared"`
  (pleiotropy) and `"union"` strategies.
* `set_residual_cov()` — store a residual covariance matrix `R` across
  traits in a new `trait_residual_cov` table. Consumed by `add_phenotype()`
  when multiple traits share the same filtered subset.
* `add_trait_covariate()` — append fixed or random covariate rows to a
  `trait_effects` table. Fixed-effect levels are serialised to a JSON-style
  VARCHAR; random effects store distribution + variance.
* `add_phenotype()` — the workhorse. Generates phenotypes for a subset of
  individuals for one or more traits and writes rows to `ind_phenotype`. Also
  computes and stores the underlying TBV in `ind_tbv`. Joint MVN residual
  draws when multiple traits share the subset and `R` is stored.
* `add_tbv()` — compute and store TBV without generating phenotype records.
* `add_ebv()` — ingest externally computed estimated breeding values into
  `ind_ebv`, tagged with a user-supplied model label.
* `add_trait_simple()` — one-shot wrapper chaining `add_trait()` +
  `define_qtl()` + `set_qtl_effects()`.

### filter() on tidybreed_pop

New S3 method `filter.tidybreed_pop()` stashes dplyr predicates on the
population object. The next `add_phenotype()` / `add_tbv()` call applies and
clears them. Multiple `filter()` calls stack with AND semantics.

```r
pop |>
  dplyr::filter(sex == "F", gen == 1L) |>
  add_phenotype("ADG")
```

### New tables

`trait_meta`, `trait_effects`, `trait_residual_cov`, `ind_phenotype`,
`ind_tbv`, `ind_ebv`. Also new columns written to `genome_meta`:
`is_QTL_{trait}` (BOOLEAN) and `add_{trait}` (DOUBLE) per trait.

### Other

* Added `MASS` to Imports for multivariate normal sampling. MASS is a
  Recommended R package that ships with base R distributions — no extra
  install burden.

---

# tidybreed 0.1.0 (2026-04-17)

## New Functions

* `add_offspring()` — core mating function. User supplies a `matings` tibble
  (one row per offspring) with required columns `id_parent_1`, `id_parent_2`,
  `sex`, and `line`. Gametes are produced via chromosomal crossover simulation
  (Haldane map: crossovers per chromosome ~ Poisson(chr_len_Mb / 100)).
  New `ind_meta`, `genome_haplotype`, and `genome_genotype` rows are written
  for all offspring. Animal-breeder aliases `id_sire` / `id_dam` are accepted.
  Any extra columns in `matings` (e.g. `gen = 2L`) are validated and written
  to `ind_meta`, with automatic `ALTER TABLE` if the column is new.
* `make_gamete()` — internal recombination helper (`R/recombination_helpers.R`).
  Not exported; used by `add_offspring()`.

---

# tidybreed 0.0.3 (2026-04-17)

## Breaking Changes

* `ind_meta`, `genome_haplotype`, `genome_genotype`: column `ind_id` renamed to
  `id_ind`; `parent_1` → `id_parent_1`; `parent_2` → `id_parent_2`. The `id_`
  prefix groups all identifier columns together for clarity.

**Migration for existing databases:**
```r
conn <- DBI::dbConnect(duckdb::duckdb(), "your_database.duckdb")
DBI::dbExecute(conn, "ALTER TABLE ind_meta RENAME COLUMN ind_id TO id_ind")
DBI::dbExecute(conn, "ALTER TABLE ind_meta RENAME COLUMN parent_1 TO id_parent_1")
DBI::dbExecute(conn, "ALTER TABLE ind_meta RENAME COLUMN parent_2 TO id_parent_2")
DBI::dbExecute(conn, "ALTER TABLE genome_haplotype RENAME COLUMN ind_id TO id_ind")
DBI::dbExecute(conn, "ALTER TABLE genome_genotype RENAME COLUMN ind_id TO id_ind")
DBI::dbDisconnect(conn)
```

## Bug Fixes

* `infer_duckdb_type()`: bare `NA` (untyped logical) now warns and returns
  `VARCHAR` instead of silently returning `BOOLEAN`
* `mutate_ind_meta()`: type check now runs before length check, so passing an
  unsupported type (e.g. a `list`) produces "Unsupported type" rather than a
  misleading length-mismatch error
* Fixed pre-existing test bugs: `expect_error(expr, NA)` inverted assertion in
  `test-add_founders.R`; null `db_conn` in `test-mutate_genome_meta.R` now uses
  a real in-memory connection

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

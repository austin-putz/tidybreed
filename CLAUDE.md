# tidybreed — Developer & AI Context

## What This Package Does

`tidybreed` simulates breeding programs. It provides a
pipe-friendly API where all genomic and individual data is stored in a
file-based DuckDB database (not in R memory). Users build a `tidybreed_pop`
object step by step, query tables with `dplyr`, and eventually run selection
and mating cycles.

## Development Status: Pre-1.0.0 — Break Freely

**The package has no users yet. Until `1.0.0`, there is NO backward-compatibility
obligation of any kind.** Make whatever breaking changes produce the best design:
rename or remove columns, functions, and arguments; change schemas; alter default
behavior. **Do not** write compatibility shims, deprecation aliases, or "legacy"
code paths, and **do not** preserve old behavior for its own sake. If a redesign is
cleaner, take it — a breaking change is not a cost to weigh, it is the expected mode
of work in this phase.

No one should be treated as a downstream user before `1.0.0`; the few people who
have seen the package have been told not to build against it. Before `1.0.0`,
breaking changes are preferred whenever they improve the long-term schema, API,
correctness, or implementation.

When a function, argument, table, or column is renamed before `1.0.0`, remove the
old name completely from the codebase. Do not keep deprecated wrappers, aliases,
manual-page examples, roxygen examples, tests, comments, or generated docs that
teach the old name unless the explicit task is to document a migration that still
exists in current code.

Leftover compatibility files are technical debt, not harmless history. If an old
entry point has been replaced, delete the old R file, remove it from `NAMESPACE`,
remove its manual page, update examples/tests/vignettes, and make the new path
the only documented path. Do not leave "deprecated for later" functions around
before `1.0.0`; there is no downstream compatibility contract to protect.

**The only reproducibility contract is forward-looking, not cross-version:**

1. **Same seed reproduces within the current code.** A given `set.seed()` (or base
   seed) must produce identical output on repeated runs of the *current*
   implementation. We do **not** care whether it matches any previous version's
   output — never write a test that compares against pre-change/"golden-from-old"
   output, and never contort a formula to stay "byte-identical to today."

   **"Identical" means bit-identical, not "within tolerance."** This binds more
   than the RNG stream: every value the simulation writes must be a function of
   the stored inputs alone, never of how many threads DuckDB happened to use.
   A parallel `SUM()` over more than two floating-point summands is the usual
   way this breaks — partial sums combine in whatever order the threads finish,
   which is not associative. The genome-effect evaluator therefore accumulates
   its one many-summand sum exactly (`GEV_ACC_TYPE` in
   `R/genome_effects_eval.R`); do not "optimize" that back to a plain `SUM()`.
   `tests/testthat/test-genome-effects-determinism.R` pins it with
   `expect_identical()`, which is the only assertion that can catch a
   regression here — `expect_equal()` passes on the broken code.
2. **R ↔ Rcpp parity.** Where the same algorithm exists in both R and C++, a given
   seed must produce identical output in both (this is a *within-current-code*
   guarantee, and the reason RNG choices like `dqrng` matter).

Before `1.0.0`, algorithms may change and seeded output may change across
versions. Within a single implementation, seeded runs should be reproducible.
After `1.0.0`, users must be able to reproduce simulations exactly for the same
package version, seed, inputs, and supported platform. Any algorithm implemented
in both R and Rcpp/C++ must have explicit parity tests showing identical results
for the same seed.

When `1.0.0` ships, this section is replaced by a normal semantic-versioning
compatibility policy. Until then: design first, break as needed.

## Design Principles

1. **Database-first** — all data lives in DuckDB, not R objects
2. **File-based by default** — enables simulations larger than RAM; supports
   replicates, resumable runs, and sharing easily
3. **Lazy evaluation** — use `dplyr::tbl()` / `get_table()` and filter before
   `collect()`-ing into R
4. **Pipe-friendly** — most exported functions accept a `tidybreed_pop` and
   return a `tidybreed_pop`; action functions (`add_phenotype`, `add_tbv`,
   `add_genotypes`, `extract_genotypes`, `define_chip`, `define_additive_effects`)
   accept a `tidybreed_table` from `get_table()` and return `tidybreed_pop`
5. **Type-safe** — all table columns have explicit DuckDB types; user-added
   columns are inferred via `infer_duckdb_type()`
6. **Disdain and intolerance for storing metadata** - storing data such as
   n_loci for loci count is silly when the function that needs that data
   can run a basic SQL query to pull out info and makes the program truly modular
7. **No implicit ordering** — never accept positional vectors (logical TRUE/FALSE
   or integer indices) to select rows from a database table. Row order is not
   guaranteed and changes silently when rows are added or reordered. All
   selection flows through named identifiers (`id_ind`, `locus_name`, `locus_id`
   as PK) or explicit SQL filter predicates via `get_table() |> filter()`.

If metadata is required, it needs to be considered carefully and likely 
should be stored in a table and never in the R object itself, for restarts
and for `restore_pop()` to be complete without having to recreate it or
given as options again. 

## Schema Design Bias

Schema choices should optimize for future breeding-program designs, especially
crossbreeding. Before adding or reshaping a table, ask whether the schema can
later support multiple lines, line-specific and population-wide effects,
reciprocal crosses, sex-specific recombination maps, non-additive effects, and
contributor-specific phenotype components without another fundamental rewrite.

Do not implement future biology before it is needed. It is enough to reserve
clean dimensions now when the schema would be painful to alter later. Prefer
long tables with explicit dimensions such as `line_name`, `sex`, `map_name`,
`contrast_name`, `component_name`, and `effect_name`. Use `NULL` deliberately for
shared/default behavior, such as population-wide genome effects or maps applying
to all lines.

## Naming Consistency Rules

1. **Function argument names must match the database column they populate exactly.**
   Never use a different name for the same thing (e.g. argument `line_name` → column
   `line_name`, not `line`).

2. **All primary numeric value columns follow the `{prefix}_value` pattern.**
   Examples: `pheno_value` in `ind_phenotype`, `tbv_value` in `ind_tbv`,
   `ebv_value` in `ind_ebv`, `cov_value` in `trait_var_comp`,
   `index_value` in `ind_index`.

3. **All name/label columns end in `_name`.**
   Examples: `trait_name`, `locus_name`, `chr_name`, `line_name`, `effect_name`,
   `index_name`. Never abbreviate to just `trait`, `line`, `chr`, etc. when
   used as a column or function parameter that maps to one of these columns.

4. **All ID foreign-key columns start with `id_`.**
   Examples: `id_ind`, `id_trait`, `id_ebv`, `id_tbv`. No `phenotype_id`-style names.

5. **No abbreviations in column names when the full word is unambiguous.**
   `index_weight` not `index_wt`; `trait_name_1`/`trait_name_2` not `trait_1`/`trait_2`.

## Two-Layer Phenotype Design (v0.31.0+)

The model is split into two distinct layers with a strict boundary between them:

**Genetic component layer** — managed by `define_trait()`:
- One row in `trait_meta` per underlying genetic quantity (e.g. `ADG_direct`, `ADG_social`, `WWD`, `WWM`)
- Has QTL effects in `genome_effects`, TBVs in `ind_tbv`, additive variance in `trait_var_comp`
- Arguments: `target_add_var`, `target_add_mean`, `description`, `units`
- No phenotype-level information at all — no mean, no residual, no type, no expressed_sex

**Observation layer** — managed by `define_phenotype()`:
- One row in `phenotype_meta` per observed phenotype individuals receive records for (e.g. `ADG`, `WW`, `mortality`)
- For simple traits, `phenotype_name` equals the `trait_name` of its single genetic component
- For composite traits (maternal, SGE), `phenotype_name` is new and one or more `trait_meta` rows feed into it via `phenotype_components`
- Arguments: `type`, `mean`, `expressed_sex`, `repeatable`, `min_value`, `max_value`, `prevalence`, `thresholds`, `cat_values`, `cat_names`, `store_liability`, `residual_var`, `components`, `formula_tbv`, `formula`, `missing_component_action`

**The rule**: if an argument describes the genetics (variance, QTL structure, parent-of-origin), it belongs in `define_trait()`. If it describes what observers record (mean, distribution, sex expression, residual noise, how to assemble from components), it belongs in `define_phenotype()`. Never put observation-layer arguments on `define_trait()` or genetic-layer arguments on `define_phenotype()`.

## Function Naming Convention

| Prefix        | Meaning                                        |
|---------------|------------------------------------------------|
| `open_`       | Opens or creates a population database/session |
| `restore_`    | Restores an existing population database       |
| `add_`        | Inserts simulation output rows                 |
| `define_`     | Writes model configuration / metadata          |
| `mutate_`     | Adds or updates columns in an existing table   |
| `extract_`    | Returns analysis/export data without changing simulation state |
| `remove_`     | Deletes selected rows                          |
| `archive_`    | Moves/stamps completed replicate data          |

**`add_*` vs `define_*` rule**: if the function writes rows that represent
simulation *output* (data produced by running the model), use `add_`. If the
function writes rows that configure *how* the model runs (parameters, weights,
effect definitions), use `define_`.

Examples: `add_founders()`, `add_phenotype()`, `add_tbv()`, `add_ebv()`,
`add_index()` — all write simulation output.  
`define_trait()`, `define_additive_effects()`, `define_effect_cov_matrix()`,
`define_chip()`, `define_index()` — all write model configuration.

## Schema and API Reference (lazy-loaded)

The per-table schema and the per-function reference live in two project skills,
loaded on demand rather than every session:

- **`tidybreed-schema`** (`.claude/skills/tidybreed-schema/SKILL.md`) — every
  table's columns, keys, reserved columns and invariants, the genome-effect
  views, and the `chr_inheritance` / `chr_recombination` resolution rules.
- **`tidybreed-api`** (`.claude/skills/tidybreed-api/SKILL.md`) — every
  implemented function: arguments, the three-stage `add_phenotype()` pipeline,
  the genome-effect evaluator, `schema()` maintenance, and custom-field
  forwarding.

**Load the matching skill before changing a table, a column, a writer, or an
exported function's behavior.** Update the skill in the same change, exactly as
this file was updated before.

## Hard Rules (always apply; detail in the skills above)

- **IDs:** integer primary keys are assigned by `next_int_id()`, never by DuckDB
  auto-increment. Individual ids never appear in SQL text — register them as a
  view (`resolve_subset_ids()`, `.ap_read_by_id()`, `R/contributor_tbv.R`).
- **RNG:** never call `dbWriteTable()` in the `add_phenotype()` path (it advances
  the RNG); sort before every RNG-consuming step; an individual without a record
  must not consume RNG or leave stochastic state. Nothing in `R/` touches
  `.Random.seed` — never add seed restoration to one function.
- **Failure contract (D7):** an `add_phenotype()` error leaves `ind_phenotype`
  and `phenotype_random_effects` untouched (the `add_tbv()` upsert is the one
  write that remains). `tests/testthat/test-add_phenotype_failure_contract.R`
  asserts it.
- **Genome effects:** never `SUM(genome_value)` over `genome_effect_loci` (an
  interaction term counts once per member). Row deletion from the three
  `genome_effect*` tables is refused — replace through
  `define_genome_effects(mode = ...)`. Do **not** re-add foreign keys *inside*
  that set: DuckDB 1.5.5 cannot delete parent and child rows in one transaction,
  so `validate_genome_effects()` checks orphans before every `COMMIT` instead
  (pinned in `tests/testthat/test-genome-effects-schema.R`).
- **One evaluator:** `add_tbv()` and `add_tgv()` share `R/genome_effects_eval.R`;
  never write a second implementation of the effect math. `add_tbv()` reads only
  reserved-owner, order-one `additive` terms.
- **No stored totals or derived values:** no `'total'` row in `ind_tgv` (use the
  `ind_tgv_total` view); `ad_terms()`'s implied mean `μ` is reported, never
  written to `phenotype_meta.mean`; no `replicate` column outside archives.
- **Covariance blocks** in `phenotype_var_comp` are declared whole, in one call,
  and locked once realized. There is no `force`.
- **Schema DDL:** `genome_meta.pos_bp` stays `BIGINT` — add columns with
  `ALTER TABLE`, never a table rewrite. `ind_phenotype`'s columns are all in the
  base `CREATE TABLE`. `ind_genotype` is an on-demand cache written only by
  `add_dosage()`.
- **Adding a table:** name it in both `.schema_table_order()` and the matching
  `.<group>_descriptions()` helper in `R/schema.R`;
  `tests/testthat/test-schema-print.R` checks that they agree.

## Roadmap

### Longer-Term

- `select_parents()` — selection index or truncation selection
- Export: PLINK `.bed/.bim/.fam`, VCF
- Visualization helpers
- Realized variance components from an arbitrary effect model
- Consolidating `ind_tbv` into `ind_tgv` (see `plans/consolidate_genetic_values.md`)

## Future Compiled Code Policy

When adding C++ code, prefer standard Rcpp/Rcpp Attributes and CRAN-compatible
dependencies that install cleanly on Linux, Windows, and macOS. Avoid required
system libraries, non-portable compiler flags, and architecture-specific code in
the main path.

Architecture-specific acceleration such as CUDA, Metal, or platform-specific
SIMD should be isolated behind optional backends with feature detection and an
R/Rcpp fallback. Do not make CUDA, Metal, or another accelerator required for
installation unless the package design deliberately splits those backends later.

## Versioning Policy

Use **three-part semantic versioning**: `MAJOR.MINOR.PATCH` (e.g. `0.0.1`).
Do **not** use the four-part devtools convention (`0.0.0.9000`).

**Pre-1.0.0 (current phase): version bumps are bookkeeping, not compatibility
promises.** Per "Development Status: Pre-1.0.0 — Break Freely" above, breaking
schema/API changes are normal development work. Document meaningful changes in
`NEWS.md`, but do not add compatibility layers. Reserve the `0 → 1` MAJOR bump
for the deliberate `1.0.0` stabilization. The table below is the policy that
takes effect **at and after 1.0.0**:

| Part  | Bump when… (≥ 1.0.0)                                 |
|-------|------------------------------------------------------|
| PATCH | Bug fixes, doc updates, minor internal changes       |
| MINOR | New exported functions or non-breaking feature additions |
| MAJOR | Breaking API changes                                 |

**Before every commit + push, update:**
1. `DESCRIPTION` — `Version:` field
2. `NEWS.md` — add an entry under the new version heading

## Design Rationale

**Why DuckDB?** Columnar, embedded (no server), SQL via dbplyr, handles
datasets larger than RAM, excellent R integration.

**Why 0/1/2 encoding?** Standard in genomics (PLINK, VCF). Easy to interpret
(count of alternate allele). Efficient integer storage.

**Why store both haplotypes and genotypes?** Haplotypes are required for
recombination and phased exports. Genotypes are required for GWAS and genomic
prediction. Computing genotypes on the fly during every query would be wasteful.

## Development Workflow

For focused changes, prefer `pkgload::load_all(".", quiet = TRUE)` plus targeted
`testthat::test_file()` calls over running the full suite first. Run broader
tests before declaring shared schema, RNG, or cross-module changes done.

Performance work must start from a small reproducible benchmark or profiling
script under `dev/benchmarks/`. Keep benchmarks deterministic, small enough to
run during development, and scalable enough to expose the intended bottleneck.

### Test Coverage (`covr`)

`covr` is in `Suggests` and is run **locally only** — there is no
`test-coverage.yaml` workflow and no Codecov integration or badge. Coverage is a
diagnostic for finding untested and dead code, not a metric to publish or chase.

Coverage measures which lines *execute* during the suite, not whether anything is
asserted about them. Treat a low number as reliable bad news and a high number as
weak good news.

**The invocation that works** (both deviations from the obvious call are load-bearing):

```r
cov <- covr::package_coverage(
  path = ".",
  type = "none",
  code = 'testthat::test_dir("tests/testthat", package = "tidybreed",
                             load_package = "installed", reporter = "summary")'
)
covr::report(cov)   # interactive HTML, uncovered lines in red
```

- **Keep `type = "none"` and scope explicitly to `tests/testthat`.** `covr`
  installs from source (`R CMD INSTALL` ignores `.Rbuildignore`), so anything
  parked in `tests/` gets executed even when `R CMD check` would skip it. The
  `.Rbuildignore` entry `^tests/test_[^/]*\.R$` still guards that slot — do not
  add unassertive demo scripts there again. The formal suite in `tests/testthat/`
  is the only test surface.
- **Do not use `testthat::test_local()`.** It calls `pkgload::load_all()`
  internally, which replaces `covr`'s instrumented package and silently reports
  **0% for every R file** while the C++ still reads ~99% (gcov instrumentation
  lives in the `.so` and survives). This fails silently as a plausible-looking
  result, not as an error. `test_dir(load_package = "installed")` uses the
  instrumented install.

A full instrumented run takes **~12 minutes** (the suite is DuckDB file-backed
and instrumentation adds 2–5×). Do not put it in a fast edit-test loop.

**Baseline at v0.63.0** — 80.7% across R code; `src/make_gametes.cpp` 98.9%.
Remaining gaps:

| File | Coverage | Why |
|------|----------|-----|
| `blupf90_helpers.R` | 22% | Needs the external BLUPF90 binary |
| `add_ebv.R` | 29% | Same — external solver dependency |
| `define_effect_cov_matrix.R` | 54% | Routing branches per `effect_name` |
| `tidybreed-package.R` | 67% | Mostly `.onLoad`/startup paths |
| `restore_pop.R` | 68% | Reopen-and-resume paths |
| `add_phenotype.R` | 71% | Composite/SGE and distribution branches |

The BLUPF90 paths are expected to stay low without a CI solver; do not chase
them. `add_tbv.R` was the standing concern at 37% and is now 99.6% — see
`tests/testthat/test-add_tbv_index.R`, which covers the `index_names` /
`weight_type` block, and `test-add_tbv.R`, which covers the line-precedence
crossbreeding join.

## Development Environment

### Running R Commands

The R executable path is platform-specific. When running R or Rscript via the
Bash tool, use the appropriate path:

**Windows (Hendrix Genetics AVD):**
- Check if working directory contains "Hendrix"
- R: `"/c/Program Files/R/R-4.5.1/bin/x64/R.exe"`
- Rscript: `"/c/Program Files/R/R-4.5.1/bin/x64/Rscript.exe"`

**Mac/Linux:**
- Use standard shell commands: `R` or `r` and `Rscript`

Example test command:
```bash
"/c/Program Files/R/R-4.5.1/bin/x64/Rscript.exe" -e "print('Hello from R')"
```

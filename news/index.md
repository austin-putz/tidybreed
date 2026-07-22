# Changelog

## tidybreed 0.58.1 (2026-07-22)

### Bug fixes

- Fixed a test that could only pass on a machine with the BLUPF90 suite
  installed.
  [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
  resolves the `renumf90` / `blupf90+` binaries *before* it creates the
  eval directory — a deliberate fail-fast order that avoids leaving
  empty eval folders behind when the suite is absent — so the assertion
  in `test-add_ebv.R` that the directory exists could never hold without
  the binaries on `PATH`. The test now skips in that case. Caught by the
  new R-CMD-check workflow on its first run; the suite had been passing
  locally only because BLUPF90 was installed.

## tidybreed 0.58.0 (2026-07-21)

### Documentation overhaul

The introduction vignette is now a real, working tour of the package
instead of a stub that pointed readers at an HTML file.
`vignettes/tidybreed-introduction.Rmd` walks the full pipeline — open a
population, define a genome and founder haplotypes, add founders, define
traits and phenotypes, assign QTL effects, compute TBVs and phenotypes,
add a genetically correlated second trait, produce offspring, define a
SNP chip, and build a selection index.

Every chunk executes at build time against an in-memory DuckDB database
(`db_name = ":memory:"`) with a fixed seed, so the vignette writes
nothing to disk, renders identical numbers on every build, and — because
it runs during `R CMD check` — fails loudly if the API changes
underneath it.

### pkgdown site and continuous integration

- New `_pkgdown.yml` publishes a documentation site to
  <https://austin-putz.github.io/tidybreed/>, with all 56 user-facing
  topics organised into 12 sections that mirror the package’s own
  workflow order and the genetic-layer / observation-layer split.
- New GitHub Actions workflows: `pkgdown.yaml` (builds and deploys the
  site) and `R-CMD-check.yaml` (Ubuntu and macOS).

### Bug fixes

- [`?schema`](https://austin-putz.github.io/tidybreed/reference/schema.md)
  showed the wrong documentation. A file-level roxygen block in
  `R/schema.R` used `@name schema` with `@keywords internal`, which
  collided with the exported
  [`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md)
  function and overwrote its topic — leaving `man/schema.Rd` with no
  `\usage`, `\arguments`, `\value`, or examples, and marked internal.
  [`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md)
  now documents itself.
- The internal `%||%` operator was marked `@keywords internal` rather
  than `@noRd`, generating a manual page whose `\name` contained an
  illegal `|` character. It is no longer documented.

### Housekeeping

- Orphaned Quarto assets moved from `docs/` to `tools/quarto/` so
  pkgdown can own `docs/`. The superseded `overview.qmd` and
  `tidybreed-quick-start.qmd` moved to `tools/quarto/legacy/`; both were
  written against the v0.9-era API and would otherwise have been
  published to the new site as if current.
- `.Rbuildignore` additions clear the “non-standard files at top level”
  check warning.

## tidybreed 0.57.0 (2026-07-15)

### Redesigned `print()` summary for a population

`print(pop)` (and typing `pop` at the console) now renders a grouped,
unicode-ruled summary instead of a flat table dump. It reuses the same
rule style as
[`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md)/[`summary()`](https://rdrr.io/r/base/summary.html)
and surfaces the data the package actually tracks:

- **Header** — population name, database location (basename, or
  `in-memory`), and live connection status (`[connected]` /
  `[disconnected]`).
- **Genome** — loci, chromosomes, physical length (Mb) and genetic
  length (cM), plus a founder haplotype-pool line when
  [`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)
  has run.
- **Model** — genetic-component traits, observed phenotypes, selection
  indices, and QTL count.
- **Individuals** — total, broken down by sex and (for multi-line
  populations) by line.
- **Records** — phenotype, TBV, EBV, and index-value row counts.

Every section is hidden when its underlying data does not yet exist, so
a freshly opened population collapses to just the header, database line,
and a pointer to
[`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md)
/
[`describe_table()`](https://austin-putz.github.io/tidybreed/reference/describe_table.md).
A closed connection prints `[disconnected]` without issuing any queries.
The old `Tables:` line (a raw comma-joined list of every table) has been
removed in favor of the `schema(pop)` pointer.

## tidybreed 0.56.0 (2026-07-08)

### Faster `ind_haplotype` write path (`add_founders()` + `add_offspring()`)

Profiling of the swine-scale founder build (40M `ind_haplotype` rows)
traced the cost to the bulk-insert PRIMARY KEY and the R-side long-frame
assembly. Two changes on the shared write core speed up both
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
and every per-generation
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
mating:

- **Dropped the `ind_haplotype`
  `PRIMARY KEY (id_ind, parent_origin, strand, locus_id)`.** A
  micro-bench showed the 4-column key’s index maintenance dominated
  insert cost (~2× faster to drop) and that dropping it also makes
  [`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md)’s
  `id_ind` DELETE ~74× faster (DuckDB prunes by zonemap, so a secondary
  index only adds overhead). Uniqueness is now guaranteed R-side —
  `id_ind` is `MAX(id)+1` from `ind_meta` (which keeps its own PK) and
  each individual emits each `(parent_origin, strand, locus_id)` tuple
  exactly once — and covered by a new test. **Schema change** (pre-1.0,
  no compatibility shim). `ind_genotype` keeps its PK (needed for
  `INSERT OR REPLACE` idempotency).
- **Replaced the `do.call(rbind, parts)` long-frame assembly with
  preallocated typed vectors** filled by slice, in both
  [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
  and the autosome path of
  [`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md).
  Output is byte-identical for a given seed; this only removes the O(n)
  `rbind` copies and per-part `data.frame` overhead.

## tidybreed 0.55.0 (2026-07-08)

### Install-time compiler guidance (source installs)

Since `v0.53.0` the package ships compiled C++17 code, so a source
install needs a C++ toolchain. This release makes that requirement
explicit and self-diagnosing.

- **README install section** now documents the exact toolchain per
  platform — **Windows** (Rtools45/44), **Ubuntu/Debian**
  (`r-base-dev` + `build-essential`), **Fedora/RHEL** (`gcc-c++` +
  `make`), and **macOS** (Xcode Command Line Tools via
  `xcode-select --install`) — plus how to force the pure-R kernel with
  `TIDYBREED_KERNEL = "r"`.
- **`configure` / `configure.win`** — new fail-safe probes. If no C++
  compiler is found they stop the install with a clear,
  platform-specific message (e.g. “no C++ compiler was found … macOS:
  xcode-select –install” or “Rtools is missing …
  <https://cran.r-project.org/bin/windows/Rtools/>”) instead of a
  cryptic build error. When a compiler is present they exit 0 and change
  nothing about the build (no flags, no `Makevars`).

## tidybreed 0.54.0 (2026-07-08)

### Swine vignette updated for the recombination refactor

- `vignettes/swine/swine-time-based-age-at-puberty-sex-semen.R` now
  exercises the v0.53.0 recombination API:
  `define_genome(cM_per_Mb = ...)` writes the default genetic map, a
  single [`set.seed()`](https://rdrr.io/r/base/Random.html) pins the
  whole replicate, and
  `add_offspring(seed = as.integer(cur_date), store_crossovers, batch_size)`
  makes each day’s matings byte-reproducible while surfacing the new
  `ind_crossover` table and memory-batching controls. Documented that
  the recombination kernel is compiled C++ by default
  (`TIDYBREED_KERNEL = "r"` forces the R reference).
- Correctness fix in the disabled cleanup block:
  `DELETE FROM ind_phenotype` keyed on `trait_name`; that table’s column
  is `phenotype_name`.

## tidybreed 0.53.0 (2026-07-07)

### C++ recombination kernel + R↔︎C++ parity (Stage 3)

The
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
autosome gamete kernel now runs in **compiled C++17** (`Rcpp` + the
`dqrng` C++ headers), driven by the same `dqrng` engine as the R
reference, so a given seed produces **byte-identical** output in R and
C++ (within-version R↔︎C++ parity; seeded output is unchanged from
0.52.0).

- **First compiled code in the package.** New
  `LinkingTo: Rcpp, dqrng, BH, sitmo` and `SystemRequirements: C++17`; a
  source install now needs a compiler (binary installs do not). The
  pure-R kernel is kept as the parity oracle and a no-compiler fallback.
- **Runtime kernel selector.** `getOption("tidybreed.kernel")` /
  `TIDYBREED_KERNEL` (`"auto"` default, `"r"` forces the R reference)
  picks the kernel; both paths produce identical results.
- **Gamete-flat, group-by-map seam.** The internal gamete generator was
  reshaped to packed integer arrays, `line_origin` integer codes, and
  long-native output (output-neutral), and now groups gametes by their
  producing parent’s resolved map so sex-specific / line-specific maps
  drive sire and dam gametes independently within one offspring.
- New dependency: **`Rcpp`** (Imports).

### dqrng per-gamete recombination + crossover storage (Stage 2)

[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
now generates gametes from **per-gamete `dqrng` sub-streams** instead of
base R’s global RNG. **This is an intentional seeded-output change**:
for a given seed, `ind_haplotype` offspring differ from 0.51.0.
Pre-1.0.0 there is no cross-version output guarantee — only forward
reproducibility (same seed reproduces on the current code).

- **`add_offspring(seed = NULL)`.** Each gamete draws from its own
  deterministic stream keyed on
  `(base seed, global offspring index, parent role)`, so output is
  **identical for any `batch_size`** (and, in a later stage, any thread
  count). `seed = NULL` draws one base seed from base R’s RNG, so an
  upstream [`set.seed()`](https://rdrr.io/r/base/Random.html) still
  reproduces the whole run; pass an integer to fix the gamete streams
  independently of surrounding base-R state. The resolved base seed is
  reported and attached as `attr(pop, "base_seed")` (not persisted to
  any table).
- **`add_offspring(store_crossovers = FALSE)`.** When `TRUE`, every
  crossover drawn during meiosis is written to the `ind_crossover` table
  (`id_ind`, `parent_origin`, `chr`, `chr_name`, `pos_cM`), including
  recombining special chromosomes. Off costs nothing. Absence of a row
  means that gamete’s chromosome did not recombine.
- **Uniform-only inversion Poisson sampler.** Crossover counts are drawn
  by a log-accumulation inversion sampler consuming only `dqrng`
  uniforms (no base `rpois`), the exact algorithm the future C++ kernel
  will mirror. Documented ceiling of 30 Morgans (3000 cM) per
  chromosome; errors above it.
- Special (sex-chromosome / organelle) inheritance also moved to
  `dqrng`, on an independent per-gamete `"special"` sub-stream.
- New dependency: **`dqrng`** (Imports).

## tidybreed 0.51.0 (2026-07-06)

### Haplotype write-path optimization (Stage 0 + Stage 1)

Internal performance work on the
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
/
[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
/
[`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)
write paths. **Output-neutral** — seeded `ind_haplotype` /
`founder_haplotypes` are byte-identical to before (verified across
autosome and sex-chromosome designs and every batch size).

- **Direct long, batched writes.** All three paths now insert
  `ind_haplotype` (and the founder pool) in long format directly,
  replacing the dense wide-matrix + `UNPIVOT` write. Peak memory is now
  **decoupled from the number of individuals** — a single high-fecundity
  mating streams to disk batch by batch. `add_offspring(200)` allocation
  churn dropped from ~43 GB to ~0.46 GB.
- **RAM-aware batching.**
  [`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
  and
  [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
  gain `batch_size` and `max_batch_mem` arguments; by default the batch
  size is auto-picked from detected available system memory (cross-OS:
  Linux, macOS Intel + Apple Silicon, Windows; conservative fallback,
  never errors). Any batch size yields byte-identical output.
- **One transaction per call.** Each
  [`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
  /
  [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
  / founder-pool write is wrapped in a single transaction — a mid-run
  failure rolls back cleanly instead of leaving a partially written
  generation.
- **New `ind_crossover` table** created (empty) by
  [`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)
  and registered across the schema, reserving crossover-event storage
  for a later stage; no rows are written yet.
- Extracted the
  [`make_gametes_batch()`](https://austin-putz.github.io/tidybreed/reference/make_gametes_batch.md)
  gamete seam (Stage 0) plus cheap structural wins
  (index-[`split()`](https://rdrr.io/r/base/split.html) parent lookup,
  preallocated special-row accumulator).

### Deprecated-code removal (pre-1.0.0 cleanup)

Per the pre-1.0.0 “break freely / no compatibility shims” policy, all
remaining deprecated wrappers and legacy code paths were removed:

- **Removed `initialize_genome()`** (the `.Deprecated` wrapper, since
  0.40.0). Use `open_pop() |> define_genome()` — and
  [`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)
  for the founder pool. `define_genome(cM_per_Mb = ...)` exposes the
  genetic-map rate the old wrapper could not.
- **Removed `set_qtl_effects_multi()`** (internal `.Deprecated`
  wrapper). Use
  [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
  with a character vector of trait names.
- **Removed `migrate_schema_meta()`** (legacy pre-v0.36.0 DB migration).
  Dead under the forward-only reproducibility contract; new databases
  always create `_schema_meta` via
  [`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md).
  Also dropped the now-orphaned internal `.genome_layer_descriptions()`
  helper.
- Cleaned all stale `initialize_genome()` references in roxygen
  examples, comments, and error strings; rewrote the loose
  `tests/test_*.R` manual scripts against the current long-format
  schema.

## tidybreed 0.50.0 (2026-07-06)

### Genome position coordinates & genetic map (breaking schema change)

Physical and genetic map positions are now stored separately, matching
the `genome_effects` precedent (sex/line/version dimensions become
**rows**, never a schema change). Pre-1.0.0: breaking changes ship as a
MINOR bump.

- **`genome_meta`: `pos_Mb DOUBLE` → `pos_bp BIGINT`** (physical
  position, base pairs, 1-based; VCF/PLINK convention, lossless export).
  `pos_bp` is created via an explicit typed `CREATE TABLE`, so user
  row-adds of an integer `pos_bp` widen to `BIGINT` automatically — no
  R-side conversion, and the type survives
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)/[`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md)/[`restore_pop()`](https://austin-putz.github.io/tidybreed/reference/restore_pop.md).
- **New `genome_map` table** holds the genetic map (`pos_cM DOUBLE`),
  keyed `(locus_id, sex, line_name, map_name)`; `sex`/`line_name` `NULL`
  = applies to both sexes / all lines; `map_name` defaults to
  `"default"`.
  [`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)
  writes a single default map; sex/line/version-specific maps are added
  as rows.
- **`define_genome(cM_per_Mb = 1.0)`** — new argument (scalar or
  length-`n_chr`, strictly positive) deriving
  `pos_cM = pos_bp/1e6 * cM_per_Mb`. Recombination and founder-LD are
  now driven by genetic distance (cM), not physical Mb.
- **`chr_meta.recombines` → `recombines_M` + `recombines_F`** (per-sex,
  matching `copy_mode_M`/`copy_mode_F`). `define_chr(recombines = ...)`
  remains the primary shorthand that sets both; pass
  `recombines_M`/`recombines_F` for single-sex achiasmy
  (e.g. *Drosophila* males).
- **Removed `genome_meta.introduced_gen`** (unused Stage-5 scaffolding).
- New internal helpers `resolve_genome_map()` (per-gamete sex/line map
  resolution with monotonicity + completeness validation) and
  `validate_genome_map()`.
- Founder-LD helpers and the recombination readers now key on the
  resolved map; the `founder_allele_freq` write no longer rewrites
  `genome_meta` (preserves the `BIGINT` type via `ALTER TABLE` +
  `UPDATE`).

Reproducibility is forward-only (same seed reproduces on the current
code); old seeded output is not preserved.

## tidybreed 0.49.1 (2026-07-03)

### Documentation audit

Full pass over every `man/` page’s roxygen source: verified
`@param`/`@return` against actual function signatures, executed every
runnable `@examples` block against a live in-memory pop, and added
examples for previously-undemonstrated key arguments (e.g. crossbreeding
TBV and `index_names`/`type = "both"` in
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md),
`economic_wts` in
[`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md),
sex-chromosome/organelle usage in
[`define_chr()`](https://austin-putz.github.io/tidybreed/reference/define_chr.md),
composite/SGE phenotype scenarios in
[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)).
Fixed pervasive stale references to the deprecated `initialize_genome()`
(superseded by `open_pop() |> define_genome()`), resolved all `Rd`
cross-reference warnings, and cleaned up numerous docs that had drifted
from the current schema/column names.

Two real (non-doc) bugs surfaced and fixed along the way:

- `TABLE_RESERVED_COLS[["ind_phenotype"]]` protected a nonexistent
  `trait_name` column instead of the actual `phenotype_name` column, so
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)/[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  custom-column writes did not correctly block overwriting
  `phenotype_name`.
- [`write_renum_par()`](https://austin-putz.github.io/tidybreed/reference/write_renum_par.md)
  (BLUPF90 export path) called a function, `load_effect_cov()`, that
  does not exist anywhere in the package — it now correctly calls
  [`load_phenotype_cov()`](https://austin-putz.github.io/tidybreed/reference/load_phenotype_cov.md)
  for the residual matrix and
  [`load_trait_cov()`](https://austin-putz.github.io/tidybreed/reference/load_trait_cov.md)
  for the additive genetic matrix.

## tidybreed 0.49.0 (2026-07-03)

### New features — sex chromosomes and organelles (Stage 4)

`chr_meta`/[`define_chr()`](https://austin-putz.github.io/tidybreed/reference/define_chr.md)
now support non-diploid inheritance rules: sex chromosomes (X/Y, Z/W,
X0/Z0) and organelles (mitochondria, plastids). Real polyploidy (ploidy
\> 2, uneven-ploidy crosses) remains out of scope — an `ind_meta.ploidy`
column is added for forward schema compatibility, but every individual
must be `ploidy = 2` in this version.

- **Schema**: `chr_meta.copies_M`/`copies_F` (absolute integers) renamed
  to `copy_mode_M`/`copy_mode_F` (`"full"`/`"half"`/`"none"`, relative
  to an individual’s own ploidy) — resolves the terminology clash
  flagged during Stage 4 planning. New `ind_meta.ploidy` column
  (`UTINYINT`, default `2`, reserved).
- **New
  [`define_chr()`](https://austin-putz.github.io/tidybreed/reference/define_chr.md)**
  (`R/define_chr.R`) — upserts a chromosome’s inheritance rule
  (`copy_mode_M`, `copy_mode_F`, `hemi_parent`, `recombines`) into
  `chr_meta`, with validation for the enum values and the
  `hemi_parent`/copy_mode consistency requirement.
- **[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)**
  gains a `ploidy` argument (must be `2`) and now writes `ind_haplotype`
  rows per chromosome according to `chr_meta`: both `parent_origin`
  slots for `"full"`, only the `hemi_parent`-designated slot for
  `"half"`, none for `"none"`. The founder haplotype-pool sampling draw
  itself is unchanged (still one flat
  [`sample()`](https://rdrr.io/r/base/sample.html) call), preserving
  byte-exact parity for existing diploid-autosome simulations.
- **[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)**
  now branches per chromosome: plain autosomes go through the original,
  unchanged recombination path; sex-linked/organelle chromosomes route
  through a new branch that resolves each contributing parent’s own copy
  count, recombines when the parent carries 2 copies (via the same
  [`make_gamete()`](https://austin-putz.github.io/tidybreed/reference/make_gamete.md),
  scoped to that one chromosome), and otherwise passes the parent’s
  single stored copy straight through via a new
  [`pass_through_gamete()`](https://austin-putz.github.io/tidybreed/reference/pass_through_gamete.md)
  helper (`R/recombination_helpers.R`) that consumes **zero** RNG draws
  when only one copy exists — so strictly hemizygous, non-recombining
  inheritance (patrilineal Y, matrilineal MT) never perturbs the random
  stream.
- **[`add_dosage()`](https://austin-putz.github.io/tidybreed/reference/add_dosage.md)/[`extract_genotypes()`](https://austin-putz.github.io/tidybreed/reference/extract_genotypes.md)**:
  the Stage 3 `assert_diploid_only()` guard (blocked *any* non-default
  `chr_meta` row) is replaced with
  [`assert_ploidy_2()`](https://austin-putz.github.io/tidybreed/reference/assert_ploidy_2.md)
  (`R/ploidy_helpers.R`), scoped to `ind_meta.ploidy` instead — their
  underlying `SUM(allele)` SQL was already correct for any row count per
  locus, so sex-linked/organelle chromosomes now work with no SQL
  changes, only the guard needed relaxing.
- **[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)**
  required no changes — its centered SQL already sums over however many
  haplotype rows exist per individual per locus.
- **`define_additive_effects(scale_to_target = TRUE)`** now errors
  clearly if any selected QTL locus is on a non-`"full"`-copy_mode
  chromosome, since the Falconer variance formula (`2*p*(1-p)*a^2`)
  assumes diploid/autosomal loci and would otherwise silently
  miscalibrate trait variance for sex-linked QTL.
- Parity preserved: `tests/testthat/test-parity.R`’s golden files are
  unchanged (no re-baseline) — confirms the diploid-autosome case is
  byte-identical before and after this stage’s restructuring.
- New tests: `test-define_chr.R`, `test-sex_chromosomes.R` (X/Y, Z/W,
  MT, X0, mixed autosome+sex-chromosome dosage/export), extended
  `test-ploidy_helpers.R`, `test-add_tbv.R`, and
  `test-define_additive_effects.R`.

## tidybreed 0.48.1 (2026-07-02)

### Hardening — dosage cache & export (Stage 3)

[`add_dosage()`](https://austin-putz.github.io/tidybreed/reference/add_dosage.md)
and
[`extract_genotypes()`](https://austin-putz.github.io/tidybreed/reference/extract_genotypes.md)
already shipped ahead of schedule in Stage 1; this release closes the
two gaps found when reviewing them against the Stage 3 exit criteria in
`plans/refactor_haplotype_stage_3.md`.

- Both functions now reject a `chip_name` whose constructed
  `is_<chip_name>`/ `has_<chip_name>` column identifier is unsafe
  (quotes, semicolons, spaces), closing a SQL-injection gap.
  Digit-leading, industry-standard chip names (e.g. `"50k"`, `"770k"`)
  continue to work unchanged, since the constructed identifier is always
  prefixed with a letter.
- All `IN (...)` list construction (`locus_names`, QTL/loci-filter
  names, individual IDs) now goes through a new internal
  [`sql_in_list()`](https://austin-putz.github.io/tidybreed/reference/sql_in_list.md)
  helper (`R/sql_utils.R`) that escapes embedded single quotes instead
  of interpolating raw values.
- New internal `assert_diploid_only()` guard (`R/ploidy_helpers.R`)
  makes both functions error clearly if `chr_meta` ever contains a
  non-diploid chromosome, rather than silently computing wrong dosage —
  forward defense ahead of Stage 4
  ([`define_chr()`](https://austin-putz.github.io/tidybreed/reference/define_chr.md)).
- Added test coverage for: cache-vs-direct extraction parity, partial
  `ind_genotype` cache population never affecting
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)/
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  results, `extract_genotypes(loci_tbl = ...)`, the above hardening, and
  the new diploid guard.

## tidybreed 0.48.0 (2026-07-02)

### New features — line-origin TBV (Stage 2)

[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
now uses `line_origin` (populated on every allele since Stage 1) to
compute correct crossbreeding additive TBV.

- Additive effects are matched to each haplotype allele by `line_origin`
  first; a population-wide effect (`genome_effects.line_name IS NULL`)
  is only used as a per-locus fallback when no line-specific row exists
  for that locus/line. This makes F1/F2/backcross TBV correct when
  different lines carry different QTL effects and/or base allele
  frequencies at the same locus.
- Centering is now folded per-allele into the summed SQL term
  (`(allele - base_allele_freq) * genome_value`), which reproduces the
  prior Falconer-centered result exactly when only population-wide
  effects exist (no behavior change for existing population-wide-only
  simulations — verified via the Stage 1 parity harness).
- Imprinted traits (`expressed_parent = "parent_1"`/`"parent_2"`)
  combine correctly with line-specific effects: the query still
  restricts to the expressed parent’s `parent_origin`.
- Internal-only `get_genotype_matrix()`/`get_haplotype_matrix()`
  helpers, superseded by the new SQL, have been removed.

### Dev tooling

- Added `dev/benchmarks/benchmark_haplotype_scale.R` — a standalone
  scale benchmark for the long-format haplotype storage (insert
  throughput for
  [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)/[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md),
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
  for both population-wide and line-specific effects, and
  [`extract_genotypes()`](https://austin-putz.github.io/tidybreed/reference/extract_genotypes.md)
  PIVOT export). Not run as part of `R CMD check`/`testthat`; see the
  script header for usage. This was the benchmark deferred from Stage 1
  (`plans/refactor_haplotype_stage_1.md`).

## tidybreed 0.47.0 (2026-07-02)

### Breaking changes — wide → long haplotype/genotype storage (Stage 1)

Haplotypes and genotypes are now stored in **long** format. This is
Stage 1 of the refactor in `plans/refactor_haplotype.md`; simulation
behavior is otherwise preserved (validated by a deterministic seeded
parity harness comparing haplotypes, dosages, TBVs, and exports
before/after).

- **Tables renamed and reshaped.** `genome_haplotype` → `ind_haplotype`
  (one row per individual × haplotype × locus:
  `id_ind, parent_origin, strand, line_origin, locus_id, locus_name, allele`;
  PK `(id_ind, parent_origin, strand, locus_id)`). `genome_genotype` →
  `ind_genotype` (long: `id_ind, locus_id, locus_name, dosage_value`; PK
  `(id_ind, locus_id)`). The old wide tables are removed.
- **`ind_genotype` is now on-demand.** It is no longer auto-populated by
  [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)/[`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md).
  Materialize dosage explicitly with the new
  [`add_dosage()`](https://austin-putz.github.io/tidybreed/reference/add_dosage.md).
- **`founder_haplotypes` is now long**
  (`line_name, haplotype_id, locus_name, allele`); the `hap_id` string
  column is replaced by an integer `haplotype_id` scoped per line.
- `genome_meta` gains an `introduced_gen` column (NULL for founding
  loci) and `locus_name` is validated to be unique.

### New features

- `add_dosage(tbl, chip_name, locus_names, overwrite_dosage)` —
  materializes 0/1/2 genotype dosages from `ind_haplotype` into the
  on-demand `ind_genotype` cache. Contrast with
  [`add_genotypes()`](https://austin-putz.github.io/tidybreed/reference/add_genotypes.md)
  (which only marks animals as physically genotyped). Idempotent via
  `INSERT OR REPLACE`.
- [`extract_genotypes()`](https://austin-putz.github.io/tidybreed/reference/extract_genotypes.md)
  gains a `loci_tbl` argument — a general locus filter via
  `get_table("genome_meta")` (e.g. autosomes only), unioned with
  `chip_name`/`effects_tbl`.
- `line_origin` is populated on every allele from the start (founders:
  their own line; offspring: inherited per-locus from the contributing
  parental segment during recombination). It is not yet *used* in TBV —
  that is Stage 2.
- New `chr_meta` table with default diploid-autosome rows
  (per-chromosome inheritance rules;
  [`define_chr()`](https://austin-putz.github.io/tidybreed/reference/define_chr.md)
  and non-default rules are Stage 4).

## tidybreed 0.46.2 (2026-07-01)

### Minor changes

- README “Global Options” section now documents all 10 `tidybreed.*`
  options (previously missing `tidybreed.archive_path`,
  `tidybreed.db_name_archive`, and `tidybreed.quiet`), with a full
  description table and archive path resolution order.
- The [`library(tidybreed)`](https://github.com/austin-putz/tidybreed)
  startup banner now prints all 10 options with their live current
  values instead of a hardcoded partial list of 5.

## tidybreed 0.46.1 (2026-06-30)

### Minor changes

- License changed from GPL-3 back to **MIT**. You are free to use,
  modify, and distribute tidybreed in commercial and proprietary
  projects provided the copyright notice is retained.

## tidybreed 0.46.0 (2026-06-26)

### Breaking changes

- `drop_rows()` renamed to
  [`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md)
  — `drop_` is not a defined prefix in the tidybreed naming convention.

### New features

- [`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md)
  now requires a prior
  [`filter()`](https://dplyr.tidyverse.org/reference/filter.html) by
  default. Calling without a filter stops with a clear error showing the
  row count and how to proceed. Pass `confirm_all = TRUE` to
  intentionally wipe an entire table.
- Bug fix: `ind_true_index` was missing from `IND_TABLE_ID_IND_COLS`, so
  `remove_rows(tables = "all")` did not clean up true index rows when
  removing individuals. Now included.

------------------------------------------------------------------------

## tidybreed 0.45.0 (2026-06-24)

### New functions

- [`mutate_group_seq()`](https://austin-putz.github.io/tidybreed/reference/mutate_group_seq.md)
  — assigns sequential integer group IDs (pen 1, 2, 3 …) to filtered
  rows. Supports fixed size, size range `c(min, max)`, and `n_groups`
  modes. Auto-continues from `MAX + 1` or restarts with `start`. Label
  reuse across unrelated rows never requires `overwrite = TRUE`.
- [`mutate_group_named()`](https://austin-putz.github.io/tidybreed/reference/mutate_group_named.md)
  — assigns named VARCHAR group labels by explicit `counts` or
  `proportions`. The `prop_method` argument selects largest-remainder
  (`"exact"`) or independent-draw (`"random"`) rounding.
  `underfill_action` handles the case where fewer animals are available
  than requested counts.
- [`mutate_group_concatenate()`](https://austin-putz.github.io/tidybreed/reference/mutate_group_concatenate.md)
  — creates a deterministic VARCHAR composite label by concatenating
  existing columns with a separator. `null_action` controls NULL
  propagation (`"propagate"` / `"skip"` / `"literal"`).

All three functions accept any `tidybreed_table` (not just `ind_meta`),
use sorted primary-key ordering before shuffling for reproducibility,
and share `overwrite`, `quiet`, and `seed` arguments following the
tidybreed conventions.

## tidybreed 0.44.1 (2026-06-24)

### Improvements

- [`archive_replicate()`](https://austin-putz.github.io/tidybreed/reference/archive_replicate.md):
  added detailed progress messages after each phase — archived tables
  with row counts, store_once tables copied vs. skipped, tables wiped
  from the working DB, and an explicit note that `tidybreed.replicate`
  was auto-incremented so users know not to manage the counter
  themselves.

## tidybreed 0.44.0 (2026-06-24)

### New functions

- [`archive_replicate()`](https://austin-putz.github.io/tidybreed/reference/archive_replicate.md)
  — accumulates multi-replicate simulation results in a separate archive
  DuckDB file without ballooning the working DB. Archives per-replicate
  tables (stamped with a `replicate` column), copies static metadata
  tables once, then deletes working-DB rows so the next replicate starts
  clean. Pipe-friendly (`pop |> archive_replicate()`). Controlled by
  three new global options: `tidybreed.replicate` (default `1L`),
  `tidybreed.archive_path` (directory), and `tidybreed.db_name_archive`
  (filename). Supports both in-process loops and HPC/SLURM array jobs
  (one archive file per job).

### Internal changes

- `"replicate"` added to `TABLE_RESERVED_COLS` for `ind_meta`,
  `ind_phenotype`, `ind_tbv`, `ind_ebv`, `ind_index`, and
  `ind_true_index` so
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  blocks this column name (it is reserved for
  [`archive_replicate()`](https://austin-putz.github.io/tidybreed/reference/archive_replicate.md)).

## tidybreed 0.43.1 (2026-06-23)

### Bug fixes

- [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md):
  directory creation for BLUPF90 runs now routes through
  [`.create_run_dir()`](https://austin-putz.github.io/tidybreed/reference/dot-create_run_dir.md)
  when `pop$run_dirs` is non-empty (i.e. when \[open_pop()\] was called
  with `tools`). Previously the eval directory was always created
  relative to the `run_dir` parameter (defaulting to the working
  directory), ignoring the managed folder structure entirely. In managed
  mode the directory now lands inside `pop$run_dirs[["blupf90"]]` with
  the standard `YYYYMMDD_HHMMSS_<6-char-random>` naming. If `"blupf90"`
  is not registered in `pop$run_dirs`,
  [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
  stops with a clear message explaining how to add it via
  `open_pop(..., tools = c("blupf90"))`. Passing `run_dir` or `eval_id`
  in managed mode now produces an informative warning that they are
  being ignored.

## tidybreed 0.43.0 (2026-06-22)

### Bug fixes

- [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md):
  added an explicit null/empty check on `mean` before building the
  insert tibble. Previously, passing `mean = NULL` (e.g. from a missing
  config key) caused `as.numeric(NULL)` → `numeric(0)`, which made
  [`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
  throw a silent length-mismatch error and left no row in
  `phenotype_meta`. Now stops immediately with a clear message pointing
  to the missing config value.

### Vignettes

- `age_at_puberty_sexed_semen.yaml`: added missing
  `wean_weight_mean: 6.0` under `general`; its absence was the root
  cause of WW not being registered in `phenotype_meta`.

## tidybreed 0.42.1 (2026-06-18)

### Documentation

- Regenerated `man/` docs to complete the v0.42.0 internal helper
  renames: deleted stale `ensure_effect_cov_table.Rd`,
  `get_effect_var.Rd`, `load_effect_cov.Rd`,
  `write_effect_cov_diagonal.Rd`; added `ensure_trait_var_comp.Rd`,
  `ensure_phenotype_var_comp.Rd`, `get_trait_var.Rd`,
  `get_phenotype_var.Rd`, `load_trait_cov.Rd`, `load_phenotype_cov.Rd`,
  `write_trait_var_diag.Rd`, `write_phenotype_var_diag.Rd`,
  `dot-migrate_var_comp_tables.Rd`.

### Vignettes

- Updated swine vignette (`swine-time-based-age-at-puberty-sex-semen.R`)
  to add `ADFI` as a seventh genetic trait and sixth phenotype. Updated
  additive genetic, residual, and pen variance/covariance matrices to
  include `ADFI`. Added symmetry checks
  ([`isSymmetric()`](https://rdrr.io/r/base/isSymmetric.html)) for all
  three matrices.

## tidybreed 0.42.0 (2026-06-09)

### Breaking changes

- **Schema rename: `trait_effect_cov` → `trait_var_comp`.** The table
  now stores only genetic-layer variance components (`gen_add`; future:
  `dominance`, `epistasis`). The primary key column is renamed
  `id_trait_var_comp`. Existing `.duckdb` files are migrated
  automatically by
  [`restore_pop()`](https://austin-putz.github.io/tidybreed/reference/restore_pop.md)
  via `ALTER TABLE RENAME`.

- **Schema rename: `phenotype_residual_cov` → `phenotype_var_comp` + new
  `effect_name` column.** The table now stores ALL phenotype-level
  variance components: residual noise (`effect_name = 'residual'`) and
  named random effects (e.g. `'hys'`, `'litter'`, `'pen'`). The primary
  key column is renamed `id_phenotype_var_comp`. Existing `.duckdb`
  files are migrated automatically by
  [`restore_pop()`](https://austin-putz.github.io/tidybreed/reference/restore_pop.md).

- **[`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md)
  routing change.** Variance for named random effects (HYS, litter, pen,
  etc.) is now written to `phenotype_var_comp` instead of
  `trait_var_comp`, correctly placing those components at the
  observation layer rather than the genetic layer.

- **[`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md)
  routing change.** Genetic effects (`gen_add`, `dominance`,
  `epistasis`) route to `trait_var_comp`; `"residual"` routes to
  [`define_residual_cov()`](https://austin-putz.github.io/tidybreed/reference/define_residual_cov.md)
  → `phenotype_var_comp`; any other name routes directly to
  `phenotype_var_comp` with that name as `effect_name`.

- **Internal helper renames** (affect code that called private helpers
  directly):

  | Old name | New name |
  |----|----|
  | `get_effect_var()` | [`get_trait_var()`](https://austin-putz.github.io/tidybreed/reference/get_trait_var.md) |
  | `load_effect_cov()` | [`load_trait_cov()`](https://austin-putz.github.io/tidybreed/reference/load_trait_cov.md) |
  | `write_effect_cov_diagonal()` | [`write_trait_var_diag()`](https://austin-putz.github.io/tidybreed/reference/write_trait_var_diag.md) |
  | `ensure_effect_cov_table()` | [`ensure_trait_var_comp()`](https://austin-putz.github.io/tidybreed/reference/ensure_trait_var_comp.md) |

  New helpers added:
  [`get_phenotype_var()`](https://austin-putz.github.io/tidybreed/reference/get_phenotype_var.md),
  [`load_phenotype_cov()`](https://austin-putz.github.io/tidybreed/reference/load_phenotype_cov.md),
  [`write_phenotype_var_diag()`](https://austin-putz.github.io/tidybreed/reference/write_phenotype_var_diag.md),
  [`ensure_phenotype_var_comp()`](https://austin-putz.github.io/tidybreed/reference/ensure_phenotype_var_comp.md).

- [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
  now accepts a `tidybreed_table` (from
  `get_table("founder_haplotypes") |> filter(...)`) as its first
  argument instead of a `tidybreed_pop`. Migrate by inserting
  `get_table("founder_haplotypes") |>` (and optionally `filter(...) |>`)
  before each
  [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
  call. The implicit `WHERE line_name = ?` / fallback-to-NULL filter is
  removed; callers now own subset selection explicitly via
  [`filter()`](https://dplyr.tidyverse.org/reference/filter.html).

  ``` r

  # Before
  pop <- add_founders(pop, n_males = 10, n_females = 50, line_name = "A")

  # After
  pop <- pop |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 10, n_females = 50, line_name = "A")

  # With explicit line filter
  pop <- pop |>
    get_table("founder_haplotypes") |>
    dplyr::filter(line_name == "Yorkshire") |>
    add_founders(n_males = 10, n_females = 50, line_name = "Yorkshire")
  ```

### Bug fixes

- Fixed
  [`compute_base_allele_freq()`](https://austin-putz.github.io/tidybreed/reference/compute_base_allele_freq.md)
  in
  [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
  to correctly select only locus columns from `founder_haplotypes`
  (previously included the `line_name` column added in v0.39.0, causing
  [`colMeans()`](https://rdrr.io/r/base/colSums.html) to fail with a
  non-numeric error when `base = "founder_haplotypes"`).

## tidybreed 0.41.0 (2026-06-09)

### New features

- Styled startup message printed on
  [`library(tidybreed)`](https://github.com/austin-putz/tidybreed).
  Shows package version, license/disclaimer, citation hint, help
  pointers, DuckDB backend version, a reminder that all data lives on
  disk in a `.duckdb` file, and the active default options. Suppress
  with `options(tidybreed.quiet = TRUE)` in `.Rprofile` or before
  loading the package.

- New option `tidybreed.quiet` (default `FALSE`) — set to `TRUE` to
  suppress the startup message entirely.

## tidybreed 0.40.0 (2026-06-08)

### Breaking changes

- `initialize_genome()` is deprecated. Replace with the new two-step
  API:

  ``` r

  pop <- open_pop(pop_name = "sim", db_name = "sim.duckdb") |>
    define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100)
  ```

### New features

- [`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md)
  — the new simulation entry point. Creates the DuckDB database, all
  core metadata tables, and the per-scenario folder structure (layers
  2–4). Reads six global options: `tidybreed.base_dir`,
  `tidybreed.output`, `tidybreed.scenario`, `tidybreed.tools`,
  `tidybreed.db_name`, `tidybreed.pop_name`. Folder creation is skipped
  for in-memory databases (`db_name = ":memory:"`).

- [`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)
  — pipe-friendly function that adds genome tables (`genome_meta`,
  `genome_haplotype`, `genome_genotype`) to an existing pop. Accepts
  `pop`, `n_loci`, `n_chr`, `chr_len_Mb`, `locus_names`, `chr_names`.

- `pop$run_dirs` — named character vector on every `tidybreed_pop`
  object mapping tool names to their pre-created layer-4 output
  directories. Always includes a `"base"` entry pointing to the layer-3
  scenario directory. `character(0)` for in-memory databases.

- `.create_run_dir(pop, tool)` (internal) — creates a timestamped
  `YYYYMMDD_HHMMSS_rrrrrr` run subdirectory inside
  `pop$run_dirs[[tool]]`. Errors immediately if `tool` is not registered
  in `pop$run_dirs`.

- [`restore_pop()`](https://austin-putz.github.io/tidybreed/reference/restore_pop.md)
  gains a `tools` argument (default:
  `getOption("tidybreed.tools", NULL)`) that reconstructs `pop$run_dirs`
  from the database location and declared tool names.

- Six global options are registered on package load via `.onLoad()`:
  `tidybreed.base_dir`, `tidybreed.output`, `tidybreed.scenario`,
  `tidybreed.tools`, `tidybreed.db_name`, `tidybreed.pop_name`.
  Pre-existing user options are preserved.

## tidybreed 0.39.0 (2026-06-07)

### New features

- [`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)
  gains a `line_name` argument. Supplying it tags each haplotype row in
  `founder_haplotypes` with a line label and prefixes `hap_id` values
  (`"LineA_hap_1"`, …). The function can now be called multiple times
  with distinct `line_name` values, enabling per-line LD structures
  (e.g. `"mosaic"` for one line, `"gaussian_copula"` for another).
  Calling twice with the same `line_name` (or twice without any
  `line_name`) still raises an error.

- [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
  now filters `founder_haplotypes` to only rows matching its own
  `line_name` argument. When no matching rows exist it falls back to
  rows with `line_name = NULL` (the shared pool from old-style single
  calls), preserving backward compatibility.

- `founder_haplotypes` table now includes a `line_name VARCHAR` column.
  Existing in-memory populations are unaffected; on-disk populations
  created before v0.39.0 that use
  [`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)
  will receive a migration-hint error and must be recreated.

## tidybreed 0.38.1 (2026-06-07)

### Documentation

- Updated swine vignette to demonstrate all six
  [`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)
  methods (`"uniform"`, `"fixed"`, `"beta"`, `"balding_nichols"`,
  `"mosaic"`, `"gaussian_copula"`), one per line (A–F).
- Added `.quarto/*` and `.DS_Store` to `.gitignore`.

## tidybreed 0.38.0 (2026-06-07)

### Breaking changes

- [`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)
  no longer accepts `fixed_allele_freq` or `allele_freq_dist`. Both are
  replaced by the new `method` dispatch argument.

  ``` r

  # Before:
  pop |> define_founder_haplotypes(n_haplotypes = 100, fixed_allele_freq = 0.5)
  pop |> define_founder_haplotypes(n_haplotypes = 100, allele_freq_dist = "uniform",
                                   min_allele_freq = 0.05, max_allele_freq = 0.95)

  # After:
  pop |> define_founder_haplotypes(n_haplotypes = 100, method = "fixed")
  pop |> define_founder_haplotypes(n_haplotypes = 100, method = "uniform",
                                   min_allele_freq = 0.05, max_allele_freq = 0.95)
  ```

### New features

- [`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)
  now supports six methods via the `method` argument:

  - `"uniform"` (default) — per-locus freq from Uniform(min, max);
    replaces `allele_freq_dist = "uniform"`
  - `"fixed"` — every locus gets the same `allele_freq` (default 0.5);
    replaces `fixed_allele_freq`
  - `"beta"` — per-locus freq from Beta(`beta_shape1`, `beta_shape2`);
    default `shape1 = shape2 = 0.5` (Jeffreys prior, biologically
    realistic MAF distribution)
  - `"balding_nichols"` — Beta parameterised by Wright’s `fst` and
    `mean_freq`; models allele frequency drift around a mean
  - `"mosaic"` — quick LD via haplotype block copying: `n_templates`
    template haplotypes are generated and new haplotypes are mosaics
    that switch templates at rate `switch_rate` per Mb; no external
    software required
  - `"gaussian_copula"` — fast LD via AR(1) latent normal thresholding;
    fully vectorised over haplotypes; inter-locus correlation decays as
    ρ = exp(−`decay_rate` × d_Mb)

- Passing a method-specific argument to the wrong method now raises an
  informative error naming both the argument and the method it belongs
  to.

- For LD methods (`"mosaic"`, `"gaussian_copula"`),
  `genome_meta$founder_allele_freq` stores the empirical column means of
  the generated haplotype matrix rather than the input parameters.

## tidybreed 0.37.0 (2026-06-07)

### Breaking changes

- `initialize_genome()` no longer accepts `n_haplotypes`,
  `fixed_allele_freq`, `allele_freq_dist`, `min_allele_freq`, or
  `max_allele_freq`. Founder haplotype generation is now a separate step
  via
  [`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md).

  ``` r

  # Before:
  pop <- initialize_genome(pop_name = "A", n_loci = 1000, n_chr = 10,
                           chr_len_Mb = 100, n_haplotypes = 100)

  # After:
  pop <- initialize_genome(pop_name = "A", n_loci = 1000, n_chr = 10,
                           chr_len_Mb = 100) |>
    define_founder_haplotypes(n_haplotypes = 100)
  ```

### New functions

- `define_founder_haplotypes(pop, n_haplotypes, fixed_allele_freq = NULL, allele_freq_dist = "uniform", min_allele_freq = 0.01, max_allele_freq = 0.99)`
  — generates a pool of founder haplotypes from per-locus Bernoulli
  distributions and writes them to the `founder_haplotypes` table. Also
  writes `founder_allele_freq` to `genome_meta`. Errors if
  `founder_haplotypes` already exists. Follows the `define_*`
  convention: configures genetic structure, not simulation output.

## tidybreed 0.36.4 (2026-06-06)

### Bug fixes

- [`define_schema_description()`](https://austin-putz.github.io/tidybreed/reference/define_schema_description.md)
  restored to returning `tidybreed_table` (v0.36.2 behaviour) so
  back-to-back chaining works without repeating
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md).
  [`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md)
  and
  [`describe_table()`](https://austin-putz.github.io/tidybreed/reference/describe_table.md)
  now also accept a `tidybreed_table` as their first argument, so
  calling `schema(pop)` no longer fails if `pop` was accidentally
  overwritten with the chain result.

## tidybreed 0.36.3 (2026-06-06)

### Bug fixes

- [`define_schema_description()`](https://austin-putz.github.io/tidybreed/reference/define_schema_description.md)
  reverted to returning `tidybreed_pop` (consistent with all other
  action functions). Returning `tidybreed_table` (v0.36.2) caused
  `schema(pop)` to fail when users assigned the result back to `pop`. To
  chain multiple column descriptions, repeat
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  between each call — each
  [`define_schema_description()`](https://austin-putz.github.io/tidybreed/reference/define_schema_description.md)
  returns `pop`, which
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  accepts.

## tidybreed 0.36.2 (2026-06-06)

### Bug fixes

- [`define_schema_description()`](https://austin-putz.github.io/tidybreed/reference/define_schema_description.md)
  now returns the `tidybreed_table` invisibly instead of
  `tidybreed_pop`, enabling back-to-back chained calls for multiple
  column descriptions on the same table. The DBI connection is a
  reference object so the original `pop` variable reflects all changes
  without reassignment.

## tidybreed 0.36.1 (2026-06-06)

### Breaking changes

- `add_schema_description()` renamed to
  [`define_schema_description()`](https://austin-putz.github.io/tidybreed/reference/define_schema_description.md)
  to match the package naming convention (`define_*` = model
  configuration / metadata).
- Signature changed from
  `(pop, table_name, description, column_name, notes)` to
  `(tbl, column_name = NULL, description, notes = NULL)` where `tbl` is
  a `tidybreed_table` from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md),
  consistent with
  [`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md)
  and all other action functions. Column name is now the first
  positional argument so
  `define_schema_description("sex", "Sex of individual")` works without
  named args.

## tidybreed 0.36.0 (2026-06-06)

### New features

- **In-database schema documentation** — descriptions for all tables and
  columns are now stored in a `_schema_meta` system table inside the
  DuckDB file. Descriptions persist across
  [`restore_pop()`](https://austin-putz.github.io/tidybreed/reference/restore_pop.md)
  restores and travel with the `.duckdb` file.
- `schema(pop)` — returns an aligned tibble of all user-visible tables
  with row counts, column counts, and descriptions.
- `describe_table(pop, "table_name")` — prints column names, types, and
  descriptions for any table; returns an invisible tibble.
- `add_schema_description(pop, table_name, description, column_name, notes)`
  — upserts a description for any table or column, including
  user-defined ones added via
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  or `...` in
  [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md).
- `migrate_schema_meta(pop)` — adds `_schema_meta` to databases created
  before v0.36.0 without recreating the database.
- [`print.tidybreed_pop()`](https://austin-putz.github.io/tidybreed/reference/print.tidybreed_pop.md)
  now shows a hint to use
  [`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md)
  and
  [`describe_table()`](https://austin-putz.github.io/tidybreed/reference/describe_table.md),
  and hides the system `_schema_meta` table from the table list.
- [`summary.tidybreed_pop()`](https://austin-putz.github.io/tidybreed/reference/summary.tidybreed_pop.md)
  now pulls table descriptions from `_schema_meta` rather than a
  hard-coded internal vector and shows them below each table header. The
  `_schema_meta` table is excluded from the summary output.
- [`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md)
  now auto-registers a column description in `_schema_meta` for the chip
  membership flag it creates.
- [`restore_pop()`](https://austin-putz.github.io/tidybreed/reference/restore_pop.md)
  now warns when opening a database that predates `_schema_meta` and
  suggests running `migrate_schema_meta()`.

### Internal

- `SYSTEM_TABLES` in `sql_utils.R` extended to include `_schema_meta`,
  `founder_haplotypes`, `phenotype_meta`, `phenotype_components`,
  `phenotype_residual_cov`, and `trait_random_effects`.

## tidybreed 0.35.1 (2026-06-06)

### Bug fixes

- [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  now validates `...` arguments against common misspellings of
  `phenotype_name` (e.g., `trait`, `phenotype`, `trait_name`,
  `pheno_name`). Previously, `phenotype_name = "ADG"` mistyped as
  `trait = "ADG"` would silently be treated as an extra column, leaving
  the selector NULL and adding **all** phenotypes instead of just ADG.
  Now raises a clear error directing users to the correct argument name.

## tidybreed 0.35.0 (2026-06-05)

### Breaking changes

- [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  no longer accepts plain R vectors (length \> 1) for setting per-row
  values. Positional vector assignment was unsafe: row order is not
  guaranteed and can silently drift when rows are added or reordered,
  causing wrong values to land on wrong rows with no error.

  Replace with an explicit tibble carrying the primary key column:

  ``` r

  # Before (unsafe):
  mutate_table(gen = c(1, 2, 3))

  # After (safe):
  mutate_table(gen = tibble::tibble(
    id_ind = c("A_1", "A_2", "A_3"),
    gen = c(1L, 2L, 3L)
  ))
  ```

  Scalar broadcasting is unchanged and unaffected.

## tidybreed 0.34.2 (2026-06-05)

### Breaking changes

- [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  argument renamed from `trait_name` to `phenotype_name`. The old name
  contradicted the two-layer design (traits = genetic components,
  phenotypes = observed records) and caused silent bugs: passing
  `phenotype_name = "ADG"` fell into `...` as an extra column, leaving
  the selector `NULL` and generating records for *all* phenotypes.
  Update all named calls: `trait_name = "ADG"` →
  `phenotype_name = "ADG"`. Positional calls (`add_phenotype("ADG")`)
  are unaffected.

## tidybreed 0.34.1 (2026-06-03)

### Breaking changes

- [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
  argument `trait_type` renamed to `type`. The old name conflicts with
  the genetic layer
  ([`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md)
  world); the observation layer should use the shorter, unambiguous
  name. Update all calls: `trait_type = "continuous"` →
  `type = "continuous"`.

- `phenotype_meta` DB column renamed from `trait_type` to `type` to
  match the argument name (per naming-consistency rule: argument names
  must match the column they populate). Existing databases must be
  recreated.

### Improvements

- [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
  man page now fully documents `formula_tbv`, `formula`, and
  `type = "derived_formula"`, with examples for continuous, count,
  categorical, maternal composite (both `components` data frame and
  `formula_tbv` shorthand), SGE, and derived-formula patterns.

- [`define_trait_simple()`](https://austin-putz.github.io/tidybreed/reference/define_trait_simple.md)
  argument `trait_type` likewise renamed to `type`.

- `sql_utils.R`: removed stale observation-layer columns (`trait_type`,
  `repeatable`, etc.) from `trait_meta` reserved cols (they moved to
  `phenotype_meta` in v0.31.0); added `phenotype_meta` reserved cols
  entry.

## tidybreed 0.34.0 (2026-05-24)

### New features

- [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
  gains two new arguments for formula-based phenotype specification,
  reducing boilerplate for common composite and derived patterns:

  - **`formula_tbv`** — a DSL string that composes the true breeding
    value from named traits using a small set of functions:

    - Bare symbol or `self(trait)` — the focal individual’s own TBV
    - `dam(trait)` — the dam’s TBV (NA for founders; handled by
      `missing_component_action`)
    - `sire(trait)` — the sire’s TBV
    - `group_sum(trait, col, table = "ind_meta")` — sum of group-mates’
      TBVs (self excluded)
    - `group_mean(trait, col, table = "ind_meta")` — mean of
      group-mates’ TBVs
    - Standard R arithmetic operators and a whitelist of math functions
      (`sqrt`, `log`, `exp`, `abs`, `round`, etc.)

    Example — maternal weaning weight:

    ``` r

    define_phenotype(pop, "WW", trait_type = "continuous",
                     mean = 230, residual_var = 180,
                     formula_tbv = "WWD + 0.5 * dam(WWM)")
    ```

    Example — SGE average daily gain:

    ``` r

    define_phenotype(pop, "ADG_obs", trait_type = "continuous",
                     mean = 850, residual_var = 300,
                     formula_tbv = "ADG_direct + group_sum(ADG_SGE, pen_id)")
    ```

    The existing `components` data frame path is retained unchanged for
    advanced cases (covariate weights, Legendre polynomial weights).

  - **`formula`** — an arithmetic expression over previously simulated
    `ind_phenotype` records, used with `trait_type = "derived_formula"`.
    No residual variance is drawn; the value is computed directly.
    Topological sort in
    [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
    ensures dependencies are evaluated first.

    Example — feed conversion ratio:

    ``` r

    define_phenotype(pop, "FCR", trait_type = "derived_formula",
                     formula = "ADFI / ADG")
    ```

    Supports chains (`FCR_pct = "FCR * 100"`), scalar coefficients
    (`"ADFI - 0.036*ADG - 0.0072*MBW"`), and whitelisted math functions
    (`"ADG / MBW^0.75"`). Division by zero and `sqrt(negative)` produce
    NA with a warning.

- [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  now topologically sorts phenotypes when any `derived_formula`
  phenotype is present, guaranteeing dependency order even when
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  is called with no phenotype name.

- `phenotype_meta` gains two new `VARCHAR` columns: `formula_tbv` and
  `formula`. Existing databases are migrated automatically on first use
  via
  [`ensure_trait_tables()`](https://austin-putz.github.io/tidybreed/reference/ensure_trait_tables.md).

- [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
  validates `formula_tbv` trait names against `trait_meta` at definition
  time, with close-match suggestions (using `agrep`) for unknown names.
  Group column existence is validated at
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  time with a clear error message and a hint to use the `table=`
  argument.

## tidybreed 0.33.0 (2026-05-23)

### Breaking changes

- `add_residual_cov()` renamed to
  [`define_residual_cov()`](https://austin-putz.github.io/tidybreed/reference/define_residual_cov.md)
  — residual covariance is model configuration, not simulation output.
- `define_effect_int()` renamed to
  [`define_effect_intercept()`](https://austin-putz.github.io/tidybreed/reference/define_effect_intercept.md)
  — abbreviation violated the no-abbreviation naming rule.
- `compute_derived()` renamed to
  [`mutate_derived()`](https://austin-putz.github.io/tidybreed/reference/mutate_derived.md)
  — aligns with the `mutate_` prefix for column-update operations.
- `delete_rows()` renamed to `drop_rows()` — `delete_` is not a defined
  prefix in the naming convention.
- `add_trait_covariate()` removed (was deprecated since v0.6.0).

## tidybreed 0.32.0 (2026-05-20)

### New features

- `contributor_type = "group"` is now fully implemented in composite
  phenotypes. The SGE (Social Genetic Effects / Bijma) model is
  supported: each focal individual’s social component is the sum (or
  mean, via `aggregation`) of pen-/cage-/litter-mates’ social TBVs. Self
  is always excluded from the group aggregation. Singletons (no
  group-mates) receive a social contribution of 0 and are not excluded.
  Relevant use-cases include pig ADG in nucleus pens, poultry cage
  mortality, litter effects, and any competition/cooperation model.

- [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
  gains a `missing_component_action` argument (`"skip"` by default, or
  `"error"`). This is stored in `phenotype_meta` and applied uniformly
  across **all** missing composite components (group assignment,
  dam/sire TBV, etc.) on every
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  call. `"skip"` returns no phenotype record for the affected individual
  and emits a warning with the count and up to 5 example IDs. `"error"`
  stops with an informative message.

### Schema changes

- `phenotype_meta` gains a
  `missing_component_action VARCHAR DEFAULT 'skip'` column. Existing
  databases are migrated automatically by
  [`ensure_trait_tables()`](https://austin-putz.github.io/tidybreed/reference/ensure_trait_tables.md).

## tidybreed 0.31.0 (2026-05-19)

### Breaking changes

- [`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md)
  no longer accepts phenotype-level arguments (`trait_type`,
  `residual_var`, `expressed_sex`, `mean`, `min_value`, `max_value`,
  `prevalence`, `thresholds`, `cat_values`, `cat_names`,
  `store_liability`, `index_weight`, `economic_value`, `recorded_on`).
  Pass these to the new
  [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
  function instead.

- `ind_phenotype.trait_name` is renamed to `phenotype_name`. SQL queries
  and direct table writes against `ind_phenotype` must use the new
  column name. Same rename applies to `trait_effects.phenotype_name` and
  `trait_random_effects.phenotype_name`.

- [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
  no longer applies `expressed_sex` filtering (that column was removed
  from `trait_meta`). Sex filtering is now handled by
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  via `phenotype_meta.expressed_sex`.

### New functions

- [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
  — registers an observed phenotype in the new `phenotype_meta` table.
  Accepts all phenotype-level metadata (mean, `expressed_sex`,
  `trait_type`, `repeatable`, categorical thresholds, prevalence, etc.)
  and an optional `residual_var` that writes to the new
  `phenotype_residual_cov` table.

- `add_residual_cov()` — writes conditional or unconditional residual
  (co)variance entries to `phenotype_residual_cov`, enabling
  heterogeneous residuals by sex, herd, or any `ind_meta` group.

### New database tables

- `phenotype_meta` — one row per observed phenotype (decoupled from
  `trait_meta`).
- `phenotype_components` — assembly rules for composite phenotypes
  (maternal, social, reaction-norm). Supports `contributor_type` values
  `"self"`, `"dam"`, and `"sire"`.
- `phenotype_residual_cov` — residual (co)variance store; replaces the
  `"residual"` rows formerly in `trait_effect_cov` for phenotype
  residuals.

### Changes

- `define_effect_cov_matrix(pop, "residual", R)` now routes to
  `phenotype_residual_cov` via `add_residual_cov()`.
- [`define_trait_simple()`](https://austin-putz.github.io/tidybreed/reference/define_trait_simple.md)
  chains through
  [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
  automatically.
- [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  reads all phenotype metadata from `phenotype_meta` and assembles
  composite TBVs from `phenotype_components` when present.
- [`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md)
  detects `phenotype_name` column in `ind_phenotype` automatically (no
  user action required).

## tidybreed 0.30.1 (2026-05-19)

### New

- Added `inst/CITATION` with two entries: an `@Manual` entry for the R
  package and an `@Misc` entry for the Quarto quick-start manual. Users
  can retrieve them with `citation("tidybreed")` or
  `toBibtex(citation("tidybreed"))`. A formal DOI will be assigned at
  v1.0.0.

## tidybreed 0.30.0 (2026-05-19)

### Breaking changes

- `set_qtl_effects_multi()` is removed from the public API. Pass a
  character vector to
  [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
  instead.

  ``` r

  # old
  tbl |> set_qtl_effects_multi(trait_names = c("ADG", "BW"), G = G)

  # new
  tbl |> define_additive_effects(c("ADG", "BW"), G = G)
  ```

  A deprecated wrapper remains and will forward calls with a warning,
  but will be fully removed in a future release.

### Enhancements

- [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
  now accepts `trait_name` as a character vector (length \>= 2) for
  correlated multi-trait effect sampling. New parameters:
  - `G` — additive-genetic (co)variance matrix; `NULL` reads from
    `trait_effect_cov`.
  - `method` — `"shared"` (default, all traits share the filtered QTL
    set) or `"union"` (per-trait QTL sets from existing `genome_effects`
    rows). Single-trait behaviour is unchanged.

## tidybreed 0.29.0 (2026-05-19)

### Breaking changes — rename `add_*` configuration functions to `define_*`

Eight functions renamed from `add_*` to `define_*` to match their
metadata/configuration semantics (pre-v1 API cleanup):

- `add_trait()` →
  [`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md)
- `add_additive_effects()` →
  [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
- `add_trait_simple()` →
  [`define_trait_simple()`](https://austin-putz.github.io/tidybreed/reference/define_trait_simple.md)
- `add_effect_fixed_class()` →
  [`define_effect_fixed_class()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_class.md)
- `add_effect_fixed_cov()` →
  [`define_effect_fixed_cov()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_cov.md)
- `add_effect_random()` →
  [`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md)
- `add_effect_cov_matrix()` →
  [`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md)
- `add_effect_int()` → `define_effect_int()`

The rule: `add_*` now exclusively means inserting simulation output rows
(individuals, phenotypes, TBVs, EBVs, index values). `define_*` means
writing model configuration that shapes how the simulation runs (trait
specs, effect definitions, variance matrices).

## tidybreed 0.28.0 (2026-05-18)

### Enhancements

- [`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md)
  is now fully generalized: any `tidybreed_table` with `id_ind`,
  `trait_name`, and a numeric value column can be passed as the first
  argument — not just `ind_ebv`. This enables phenotypic indexes
  (`ind_phenotype`), TBV-based indexes (`ind_tbv`), and user-defined
  tables.
  - New `value_col` argument (default `NULL`): auto-detected from the
    table name (`ind_ebv` → `"ebv_value"`, `ind_phenotype` →
    `"pheno_value"`, `ind_tbv` → `"tbv_value"`). Pass explicitly for
    unknown or user-defined tables.
  - A uniqueness check errors if any `(id_ind, trait_name)` combination
    has more than one row after filtering, preventing silent
    mis-computation from repeat phenotype records or multiple EBV
    evaluations.
  - A warning is issued when no filter is applied, reminding users to
    narrow to a single record per `(id_ind, trait_name)` before calling
    [`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md).

## tidybreed 0.27.1 (2026-05-18)

### Administrative

- License changed from MIT to **GPL-3**. Any modifications or derivative
  works distributed publicly must also be released under GPL-3,
  protecting the project’s open-source contributions.

## tidybreed 0.27.0 (2026-05-16)

### New features

- New
  [`define_table()`](https://austin-putz.github.io/tidybreed/reference/define_table.md)
  function. Creates a user-defined custom table inside the tidybreed
  DuckDB database and registers it in `pop$tables` so that
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  and
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  work on it immediately. Column names and types are declared with named
  typed-NA `...` arguments (same convention as
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  schema pre-declaration). Optional `primary_key` argument names a
  column as the `PRIMARY KEY`. Use
  [`DBI::dbAppendTable()`](https://dbi.r-dbi.org/reference/dbAppendTable.html)
  to insert rows after creation.

## tidybreed 0.26.0 (2026-05-16)

### New features

- [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
  gains three new arguments for computing true selection index values
  from TBVs:
  - `index_names`: character vector of named indices (from
    [`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md))
    for which a true index value should be computed. When `NULL`
    (default), only `ind_tbv` is written — existing behaviour is
    unchanged.
  - `type`: which weight column from `index_meta` to use — `"index"`
    (default, `index_weight`), `"economic"` (`economic_weight`), or
    `"both"` (writes two rows per individual distinguished by the
    `weight_type` column).
  - `overwrite_index`: when `FALSE` (default), individuals already
    present in `ind_true_index` for a given `(index_name, weight_type)`
    are skipped — avoids redundant recomputation across generations. Set
    to `TRUE` to recompute when index weights change.
- New table `ind_true_index` (`id_true_index`, `id_ind`, `index_name`,
  `weight_type`, `true_index_value`). Created automatically by
  [`ensure_trait_tables()`](https://austin-putz.github.io/tidybreed/reference/ensure_trait_tables.md)
  on the next call that triggers it (e.g., `add_trait()`,
  `initialize_genome()`). Existing databases are migrated transparently.

## tidybreed 0.25.0 (2026-05-16)

### New features

- [`extract_genotypes()`](https://austin-putz.github.io/tidybreed/reference/extract_genotypes.md)
  gains an `effects_tbl` argument. Pass a filtered `tidybreed_table`
  from `get_table(pop, "genome_effects")` to extract genotypes at QTL
  loci (selected by trait, effect type, effect size, line, or any other
  filter). `chip_name` and `effects_tbl` may be used together; locus
  sets are unioned and deduplicated. `chip_name` now has a `NULL`
  default (backward compatible — positional calls like
  `extract_genotypes(tbl, "50k")` are unchanged).

## tidybreed 0.24.1 (2026-05-16)

### Naming consistency fixes

- `index_meta.economic_value` renamed to `index_meta.economic_weight` to
  match the `index_weight` column and the `{prefix}_weight` convention.
  Existing databases are migrated automatically.
- [`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md)
  parameter `economic_values` renamed to `economic_wts` to match the
  `index_wts` parameter naming convention.
- `CLAUDE.md` updated: `index_meta` and `ind_index` table schemas
  documented;
  [`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md)
  /
  [`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md)
  added to the Implemented Functions section.

## tidybreed 0.24.0 (2026-05-16)

### New features

- `index_meta` gains an `economic_value DOUBLE` column. Existing
  databases are migrated automatically on the next call to
  [`ensure_trait_tables()`](https://austin-putz.github.io/tidybreed/reference/ensure_trait_tables.md).
- `add_trait()` now writes a global economic-value entry to `index_meta`
  (`index_name = NULL`) alongside its `trait_meta` row. When
  `overwrite = FALSE` (the default), an existing entry is preserved.
  When `overwrite = TRUE`, the `economic_value` is updated in place.
- [`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md)
  gains two new parameters:
  - `economic_values` — optional numeric vector (same length as
    `trait_names`) written to the `index_meta.economic_value` column.
    Some values may be `0`.
  - `overwrite = FALSE` — re-calling for an existing
    `(index_name, trait_name)` pair is now a no-op by default. Set
    `overwrite = TRUE` to update weights and economic values in place
    (the previous unconditional-upsert behaviour).

## tidybreed 0.23.0 (2026-05-15)

### Breaking changes — genome_effects table replaces genome_meta effect columns

- New `genome_effects` table stores all QTL effect data in long format.
  The dynamically-added wide columns `add_{trait}`, `is_QTL_{trait}`,
  and `base_allele_freq_{trait}` no longer exist in `genome_meta`.
  Existing databases created with prior versions are incompatible.

- `define_qtl()` is removed. Its role is absorbed by
  `add_additive_effects()`, which now accepts a filtered
  `tidybreed_table` from `get_table("genome_meta")` as its first
  argument. QTL selection and effect assignment happen in one step.

  **Migration:**

  ``` r

  # Old (v0.22.x):
  pop |> get_table("genome_meta") |> filter(...) |> define_qtl("T") |> add_additive_effects("T")

  # New (v0.23.0):
  pop |> get_table("genome_meta") |> filter(...) |> add_additive_effects("T")
  ```

- `add_additive_effects()` first argument is now always a
  `tidybreed_table` from `get_table("genome_meta")` (never a bare
  `tidybreed_pop`). The `base = "current_pop"` individual filter is now
  supplied via the new `base_tbl` argument instead of piping from an
  `ind_meta` table.

- `set_qtl_effects_multi()` first argument is now a `tidybreed_table`
  from `get_table("genome_meta")` (defining the shared/union loci set)
  instead of a `tidybreed_pop`. Add `line_name` and `base_tbl`
  arguments.

### Improvements

- [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  no longer duplicates TBV computation. It now calls
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
  internally and reads TBVs from `ind_tbv` for the phenotype model.
- `genome_effects` table supports a `line_name` column for future
  line-specific QTL effects (`NULL` = population-wide, the default).

## tidybreed 0.22.1 (2026-05-15)

### Breaking changes — rename `set_qtl_effects()` to `add_additive_effects()`

`set_qtl_effects()` has been renamed to `add_additive_effects()` to
align with the `add_` naming convention and to prepare for future
`add_dominance_effects()` and other effect-type functions. All callers
(tests, vignettes, documentation) have been updated.

**Migration:** replace every `set_qtl_effects(` with
`add_additive_effects(`. `set_qtl_effects_multi()` is unchanged.

## tidybreed 0.22.0 (2026-05-15)

### Breaking changes — eliminate implicit locus ordering

All positional locus selection methods have been removed. Loci are now
always selected via the
`get_table("genome_meta") |> filter(...) |> fn()` pattern, which is
order-safe and consistent with individual selection throughout the
package.

#### `define_chip()` — signature change

Old interface (removed):

``` r

define_chip(pop, chip_name, n, method, locus_tf, locus_ids, locus_names)
```

New interface — pipe a filtered `genome_meta` table:

``` r

pop |> get_table("genome_meta") |> filter(...) |> define_chip(chip_name)
```

The helper functions `select_by_n()`, `select_by_locus_ids()`, and
`select_by_locus_names()` in `chip_helpers.R` have been removed
entirely.

#### `define_qtl()` — signature change

Old interface (removed):

``` r

define_qtl(pop, trait_name, n, method, locus_tf, locus_ids, locus_names)
```

New interface — pipe a filtered `genome_meta` table; `trait_name` is now
optional:

``` r

pop |> get_table("genome_meta") |> filter(...) |> define_qtl("ADG")
pop |> get_table("genome_meta") |> filter(...) |> define_qtl(c("ADG", "BW"))
pop |> get_table("genome_meta") |> filter(...) |> define_qtl()  # all traits
```

Shared QTL pleiotropy pattern:

``` r

pop |> get_table("genome_meta") |> filter(is_QTL_ADG == TRUE) |> define_qtl("BW")
```

`define_qtl()` still returns `tidybreed_pop`, so it can pipe into
`set_qtl_effects()`.

#### `add_trait_simple()` — `qtl_method` argument removed

The `qtl_method` parameter no longer exists. QTL are always placed
randomly inside `add_trait_simple()`. For non-random placement, use
`define_qtl()` directly after filtering `genome_meta`.

#### `add_phenotype()` and `add_tbv()` — `trait_name` now optional

When omitted, all traits in `trait_meta` are used (in `id_trait` order).

``` r

pop |> get_table("ind_meta") |> add_phenotype()  # all traits
pop |> get_table("ind_meta") |> add_tbv()         # all traits
```

## tidybreed 0.21.0 (2026-05-15)

### Breaking changes — naming consistency

All column and argument names have been standardised. Scripts that
reference any of the old names must be updated.

#### Database column renames

| Table              | Old column | New column     |
|--------------------|------------|----------------|
| `ind_meta`         | `line`     | `line_name`    |
| `ind_phenotype`    | `value`    | `pheno_value`  |
| `ind_tbv`          | `tbv`      | `tbv_value`    |
| `ind_ebv`          | `ebv`      | `ebv_value`    |
| `trait_effect_cov` | `trait_1`  | `trait_name_1` |
| `trait_effect_cov` | `trait_2`  | `trait_name_2` |
| `trait_effect_cov` | `cov`      | `cov_value`    |
| `index_meta`       | `index_wt` | `index_weight` |

#### Public function parameter renames

- `add_phenotype(tbl, trait, ...)` — `trait` → `trait_name`
- `add_tbv(tbl, trait, ...)` — `trait` → `trait_name`
- `add_ebv(tbl, trait_name, ...)` — already correct; now consistent
  throughout
- Internal helpers (`build_data_file`, `parse_blupf90_solutions`,
  `update_covars_from_blupf90`, `compute_covariate_contribution`,
  `next_pheno_numbers`, `load_effect_cov`) — `trait` → `trait_name`

#### CLAUDE.md — new naming consistency rules

Five rules now documented in `CLAUDE.md`: 1. Function argument names
must match the database column they populate exactly. 2. All primary
numeric value columns follow the `{prefix}_value` pattern. 3. All
name/label columns end in `_name`. 4. All ID columns start with `id_`.
5. No abbreviations when the full word is unambiguous.

## tidybreed 0.20.1 (2026-05-14)

### Documentation & cleanup

- `NAMESPACE` now exports `compute_derived()` (was missing from v0.20.0
  commit).
- `man/compute_derived.Rd` added (roxygen docs for `compute_derived()`).
- `man/add_trait.Rd` and `man/add_phenotype.Rd` updated to reflect
  v0.19.0 binary → categorical unification (`cat_values`, `cat_names`,
  `store_liability` params; removed `"binary"` from `trait_type`
  choices).
- Added swine basic replication/generation YAML vignette script.

## tidybreed 0.20.0 (2026-05-14)

### New features

- `compute_derived()` — new action function that joins a filtered
  primary table with an optional secondary table, applies a
  user-supplied R function to produce a derived vector, and writes the
  result to one or more destination columns in any combination of
  tables. Reduces ~25-line dplyr workflows for computed date fields
  (e.g. `puberty_date = birth_date + AP_value`) to a single declarative
  pipe step. Accepts `write_to = c(table = "col")` syntax to write the
  same computed value under different column names in different tables.

## tidybreed 0.19.0 (2026-05-14)

### Breaking change

- `trait_type = "binary"` has been **removed**. Use
  `trait_type = "categorical"` with `prevalence = <p>` (for a 2-category
  trait) or `thresholds = c(<value>)` and add `cat_values = c(0, 1)` to
  recover 0/1 encoding. Binary was a strict subset of categorical — this
  eliminates the duplication.

### New features

- `add_trait()` gains three new parameters for categorical traits:

  - `cat_values` — numeric vector (length = n_categories) giving the
    phenotype value stored in `ind_phenotype` for each category
    (e.g. `c(0, 1)` for standard binary, or `c(1, 2, 3, 4)` for a
    4-level score). Defaults to `1, 2, ..., K` when `NULL`.
  - `cat_names` — character vector of human-readable labels per category
    (e.g. `c("Alive", "Dead")`), written to the reserved `cat_name`
    column in `ind_phenotype`. Must not contain commas.
  - `store_liability` — logical. When `TRUE`, the raw liability value on
    the continuous scale is written to the reserved `liability_value`
    column in `ind_phenotype` alongside the observed category value.

- `add_trait()` now accepts `prevalence` for `categorical` traits
  (previously only accepted for `binary`). A single prevalence value
  defines a 2-category threshold; the actual threshold is computed at
  phenotype time from
  `target_add_mean + qnorm(1 - prevalence) * sqrt(VA + VR)`.

- `add_trait()` validates that exactly one of `thresholds` or
  `prevalence` is supplied for categorical traits; supplying both is an
  error.

- Two reserved column names are added to `ind_phenotype`:
  `liability_value` (DOUBLE) and `cat_name` (VARCHAR). Both are added
  dynamically via `ALTER TABLE` on first use, so existing databases are
  not affected until a categorical trait with these options is
  phenotyped.

- `trait_meta` gains three new columns: `cat_values VARCHAR`,
  `cat_names VARCHAR`, `store_liability BOOLEAN`. Existing databases are
  migrated automatically on the next call to
  [`ensure_trait_tables()`](https://austin-putz.github.io/tidybreed/reference/ensure_trait_tables.md).

## tidybreed 0.18.3 (2026-05-14)

### Bug fix / Enhancement

- [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  now enforces the `repeatable` field from `trait_meta`. If
  `repeatable = FALSE` and individuals in the candidate set already have
  a phenotype record for that trait, they are silently skipped and a
  warning is issued reporting the number rejected and the number that
  will receive a new record. Individuals with no prior record are
  phenotyped as usual.

## tidybreed 0.18.2 (2026-05-14)

### Bug fix

- [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  binary trait now uses a **fixed theoretical threshold**
  (`target_add_mean + qnorm(1 - prevalence) * sqrt(VA + VR)`) rather
  than the empirical quantile of the current batch. The old approach
  mechanically forced exactly the target prevalence; the new approach
  correctly shows sampling variation around the target, matching the
  standard threshold model from quantitative genetics.
- Updated the existing testthat binary-trait test to use explicit
  absolute tolerance (`abs(rate - 0.1) <= 0.04`) instead of relative
  tolerance.
- Added `tests/test_binary_mortality.R` — a stand-alone demo script that
  sets up 10% mortality with VA = 1, VR = 9 (h² = 0.10) and verifies the
  observed rate is within 1.5 pp of the target.

## tidybreed 0.18.1 (2026-05-14)

### Documentation

- Removed all references to removed functions (`mutate_ind_meta()`,
  `mutate_genome_meta()`, `set_residual_cov()`) from roxygen examples,
  `@seealso` sections, vignettes, and test scripts; replaced with
  current API.
- Replaced `[add_trait_covariate()]` in `@seealso` of `add_trait()` and
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  with the current `add_effect_*()` functions.
- Tightened `@param` type documentation across all exported functions to
  explicitly state scalar vs. vector, valid choices, and R type
  suffixes.
- Updated vignette subtitle to v0.18.1 and corrected table name
  `ind_phenotype`.

## tidybreed 0.18.0 (2026-05-14)

### New features

- [`close_pop()`](https://austin-putz.github.io/tidybreed/reference/close_pop.md)
  gains a `results_dir` argument. When provided, the `.duckdb` file is
  moved to that directory after the connection is closed. The directory
  is created automatically if it does not exist. An error is raised if a
  file with the same name already exists at the destination. In-memory
  databases are silently ignored.

## tidybreed 0.17.0 (2026-05-14)

### Breaking changes

- `pop$metadata` has been removed from the `tidybreed_pop` S3 object.
  Any code that reads `pop$metadata$n_loci`, `pop$metadata$n_chr`,
  `pop$metadata$chr_len_Mb`, `pop$metadata$chr_names`, or
  `pop$metadata$n_individuals` will need to query the database directly
  (e.g. `SELECT COUNT(*) FROM genome_meta`).

### Improvements

- **Database is the single source of truth.** Genome shape (`n_loci`,
  chromosome lengths, etc.) is no longer cached on the pop object at
  initialization. All functions that need this information query
  `genome_meta` live, so they automatically reflect any changes made to
  the table.

- **Dynamic genomes are now possible.** Because no genome summary is
  cached, users can add rows to `genome_meta` and columns to
  `genome_haplotype` / `genome_genotype` after initialization (e.g. for
  novel mutations, gene editing, or adding newly discovered QTL).
  [`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
  and
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
  will automatically pick up the updated locus set without any metadata
  refresh step.

- [`build_chr_info()`](https://austin-putz.github.io/tidybreed/reference/build_chr_info.md)
  no longer takes a `chr_len_Mb` argument. Chromosome length for
  recombination is now derived as `MAX(pos_Mb)` per chromosome from the
  in-memory `genome_meta_df` slice — genomically equivalent since
  crossovers beyond the last locus position cannot affect observed
  alleles.

- [`restore_pop()`](https://austin-putz.github.io/tidybreed/reference/restore_pop.md)
  simplified: no metadata reconstruction block. Function now just opens
  the connection, verifies `genome_meta` is non-empty, infers
  `pop_name`, and returns the pop object.

## tidybreed 0.16.0 (2026-05-14)

### New features

- `restore_pop(db_path)` — reconnects to an existing `.duckdb` file and
  reconstructs a fully operational `tidybreed_pop` object. Allows
  simulations to be resumed after
  [`close_pop()`](https://austin-putz.github.io/tidybreed/reference/close_pop.md),
  a crash, or in a fresh R session without any stored metadata tables.
  All runtime metadata (`n_loci`, `n_chr`, `chr_len_Mb`, `chr_names`) is
  derived on the fly from the current state of `genome_meta`. An
  optional `pop_name` argument overrides the name inferred from the
  filename.

## tidybreed 0.15.2 (2026-05-13)

### New features

- [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md):
  new optional `phenotype` argument accepts a `tidybreed_table` from
  `get_table("ind_phenotype") |> filter(...)`. When supplied, only
  phenotype records matching the filter are included in the BLUPF90 data
  file; the pedigree (and which animals receive EBVs) is unaffected.
  This solves the “future phenotype” problem where records are sampled
  early (e.g., number weaned, age at puberty) but should be excluded
  from evaluations until actually observed. Ignored with a warning in
  `parent_avg` mode.

## tidybreed 0.15.1 (2026-05-13)

### Bug fixes

- [`print.tidybreed_table()`](https://austin-putz.github.io/tidybreed/reference/print.tidybreed_table.md):
  fixed crash when
  [`select()`](https://dplyr.tidyverse.org/reference/select.html) is
  called before printing. Previously, the print method derived column
  names from
  [`DBI::dbListFields()`](https://dbi.r-dbi.org/reference/dbListFields.html)
  (the full physical table schema) and tried to select them from `x$tbl`
  via
  [`dplyr::all_of()`](https://tidyselect.r-lib.org/reference/all_of.html).
  After a
  [`select()`](https://dplyr.tidyverse.org/reference/select.html) call
  narrows the lazy tbl to a subset of columns, those original column
  names no longer exist in `x$tbl`, causing
  [`dplyr::all_of()`](https://tidyselect.r-lib.org/reference/all_of.html)
  to error. The fix uses `colnames(x$tbl)` instead, which always
  reflects the current visible columns. The header field count now also
  matches the narrowed column set.

## tidybreed 0.15.0 (2026-05-13)

### Breaking changes

- `id_ind` format changed from `{line_name}-{n}` (e.g. `Libra-1020`) to
  `{line_name}_{n}` (e.g. `Libra_1020`). Underscore separator avoids SQL
  query issues and DML-string construction with glue/paste. **Existing
  `.ddb` files are incompatible — recreate databases with this
  version.**

- `line_name` no longer accepts hyphens. Only letters, digits, and
  underscores are allowed (regex `^[a-zA-Z][a-zA-Z0-9_]*$`).

- `ind_phenotype`: `id_record VARCHAR` replaced by
  `id_phenotype INTEGER`. The new column is a global auto-incrementing
  integer primary key rather than a per-trait formatted string like
  `"ADG-529"`.

- New integer primary-key columns added to all major tables:

  | Table              | New column                    |
  |--------------------|-------------------------------|
  | `ind_tbv`          | `id_tbv INTEGER`              |
  | `ind_ebv`          | `id_ebv INTEGER`              |
  | `ind_index`        | `id_index INTEGER`            |
  | `trait_meta`       | `id_trait INTEGER`            |
  | `trait_effect_cov` | `id_trait_effect_cov INTEGER` |
  | `index_meta`       | `id_index_name INTEGER`       |

  Each column holds a globally unique auto-incrementing integer, making
  it easy to reference specific rows without composing multi-column
  keys. Former composite primary keys are now enforced at the
  application level (unique constraint retained where relevant).

- New internal helper `next_int_id(conn, table, id_col)` in
  `R/sql_utils.R` returns `MAX(id_col) + 1` for any table, used by all
  insert functions.

## tidybreed 0.14.0 (2026-05-13)

### New features

- New `delete_rows(tbl, tables = NULL, dry_run = FALSE, verbose = TRUE)`
  — delete rows from one or more population tables using the
  `get_table() |> filter() |> delete_rows()` pipeline.
  - `tables = NULL` (default): delete from the current table only, using
    the full composite row key (e.g. filters on `ind_tbv` by
    `trait_name` delete only that trait’s rows, not all rows for those
    animals).
  - `tables = "all"`: delete matching individuals from every `ind_*`
    table plus `genome_haplotype` and `genome_genotype` that exist in
    the population (includes `ind_meta` for a complete “hard delete”).
  - `tables = c(...)`: delete from an explicit subset of tables; each
    must have an `id_ind` column.
  - `dry_run = TRUE`: reports what would be deleted without modifying
    the DB.
  - Uses a temp-table JOIN internally (consistent with
    [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md));
    safe for large ID sets.
- New internal constants in `R/sql_utils.R`:
  - `TABLE_ROW_KEYS`: composite row key columns per table, used for
    exact single-table deletion.
  - `IND_TABLE_ID_IND_COLS`: character vector of all tables with
    `id_ind`.

## tidybreed 0.13.2 (2026-05-12)

### Performance

- [`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
  /
  [`make_gamete()`](https://austin-putz.github.io/tidybreed/reference/make_gamete.md):
  ~3x speedup for gamete simulation.
  - Replaced the per-locus `for` / `while` loop in
    [`make_gamete()`](https://austin-putz.github.io/tidybreed/reference/make_gamete.md)
    with [`findInterval()`](https://rdrr.io/r/base/findInterval.html),
    which counts crossovers before each locus position in vectorized C
    code. Haplotype assignment becomes a single arithmetic expression:
    `(current_hap - 1 + n_toggles %% 2) %% 2 + 1`.
  - Added
    [`build_chr_info()`](https://austin-putz.github.io/tidybreed/reference/build_chr_info.md)
    helper that pre-computes per-chromosome locus indices, positions,
    and lengths once per
    [`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
    call. Previously, `genome_meta_df$chr == chr_id` (an O(n_loci) scan)
    was repeated for every chromosome × every gamete. The result is
    passed into each
    [`make_gamete()`](https://austin-putz.github.io/tidybreed/reference/make_gamete.md)
    call, eliminating redundant work.

## tidybreed 0.13.1 (2026-05-11)

### Bug fixes

- [`format_sql_value()`](https://austin-putz.github.io/tidybreed/reference/format_sql_value.md):
  added class-based guards for `Date` and `POSIXct`/`POSIXlt` values
  **before** the `db_type` string-matching branches. If a `Date` object
  loses its class attribute (e.g. via R’s `for`-loop iterator stripping
  class on some R builds),
  [`infer_duckdb_type()`](https://austin-putz.github.io/tidybreed/reference/infer_duckdb_type.md)
  could fall through to `"DOUBLE"`, causing `as.character(date)` =
  `"2026-02-26"` to be embedded **unquoted** in the UPDATE SQL. DuckDB
  then interpreted the bare expression as integer arithmetic
  (`2026 − 02 − 26 = 1998`) and threw
  `Conversion Error: Unimplemented type for cast (INTEGER -> DATE)`. The
  fix ensures Date and TIMESTAMP values always produce `DATE '...'` /
  `TIMESTAMP '...'` literals regardless of what `db_type` was inferred.

## tidybreed 0.13.0 (2026-05-10)

### Bug fixes

- [`infer_duckdb_type()`](https://austin-putz.github.io/tidybreed/reference/infer_duckdb_type.md):
  moved `inherits(value, "Date")` and
  `inherits(value, "POSIXct"/"POSIXlt")` checks **before**
  [`is.numeric()`](https://rdrr.io/r/base/numeric.html). `Date` objects
  are stored internally as doubles, so
  [`is.numeric()`](https://rdrr.io/r/base/numeric.html) could return
  `TRUE` for them in some configurations, causing the type to be
  inferred as `DOUBLE` instead of `DATE`. The symptom was
  `SET mate_date = 2026-02-26` (unquoted) in the UPDATE SQL, which
  DuckDB interpreted as integer arithmetic and threw
  `Conversion Error: Unimplemented type for cast (INTEGER -> DATE)`.

- [`format_sql_value()`](https://austin-putz.github.io/tidybreed/reference/format_sql_value.md):
  date and timestamp literals now use the explicit SQL type-keyword
  syntax (`DATE '2026-02-26'`, `TIMESTAMP '...'`) instead of bare quoted
  strings. This removes any remaining ambiguity about literal types
  regardless of column context.

## tidybreed 0.12.0 (2026-05-03)

### New features

- [`summary.tidybreed_pop()`](https://austin-putz.github.io/tidybreed/reference/summary.tidybreed_pop.md)
  — detailed per-table summary with frequency tables for low-cardinality
  columns, 5-number summaries for numeric columns, date ranges for
  date/timestamp columns, and safe handling of wide genome tables (locus
  columns are not SUMMARIZE’d). Dispatches via `summary(pop)`.
- [`print.tidybreed_summary()`](https://austin-putz.github.io/tidybreed/reference/print.tidybreed_summary.md)
  — tidyverse-style console display with box-drawing separators and
  aligned column output. Accepts optional `tables` and `max_values`
  parameters to control scope and display thresholds.

## tidybreed 0.11.3 (2026-05-02)

### Bug fixes

- [`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md):
  replaced [`stats::reshape()`](https://rdrr.io/r/stats/reshape.html) +
  [`sub()`](https://rdrr.io/r/base/grep.html) rename with a direct
  matrix-based EBV pivot.
  [`stats::reshape()`](https://rdrr.io/r/stats/reshape.html) behaves
  differently on tibbles (returned by
  [`dplyr::collect()`](https://dplyr.tidyverse.org/reference/compute.html)
  in the filter branch) vs plain data frames (returned by
  [`DBI::dbGetQuery()`](https://dbi.r-dbi.org/reference/dbGetQuery.html)
  in the no-filter branch), causing column name mangling that made
  `wide[[tr]]` return `NULL` and threw
  `vapply … values must be length 1, but FUN(X[[1]]) result is length 0`.
  The fix builds an `n_ind × n_traits` matrix using direct row lookup
  and [`match()`](https://rdrr.io/r/base/match.html), then computes
  index values via `%*%`. No behaviour change for the no-filter path;
  the filter path now works correctly for any number of traits.

## tidybreed 0.11.2 (2026-05-02)

### Bug fixes

- [`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md):
  renamed `replace_index` parameter to `overwrite_index` to prevent R’s
  partial argument matching from silently consuming user-defined extra
  column values named `rep`. Previously, passing `rep = 1L` as a custom
  column would be matched to `replace_index` (since `"rep"` is a unique
  prefix of `"replace_index"`), causing a `stopifnot(is.logical(...))`
  error. Now consistent with `overwrite_trait` in
  [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md).

## tidybreed 0.11.1 (2026-05-02)

### Bug fixes

- [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md):
  renamed `replace_trait` parameter to `overwrite_trait` to prevent R’s
  partial argument matching from silently consuming user-defined extra
  column values named `rep`. Previously, passing `rep = 1L` as a custom
  column would be matched to `replace_trait` (since `"rep"` is a unique
  prefix of `"replace_trait"`), causing the value to never reach `...`
  and leaving the `rep` column at its declared default (`NA`).

## tidybreed 0.11.0 (2026-05-01)

### New features

- New `define_index(pop, index_name, trait_names, index_wts, ...)`
  registers a named selection index in the new `index_meta` table. Not
  all population traits need to be indexed — specify only the traits
  with non-zero economic weights. Supports extra user-defined columns
  via `...` (same pattern as other `add_*()` / `define_*()` functions).
  Re-calling with the same `(index_name, trait_name)` updates the
  weights in place (upsert).

- New
  `add_index(tbl, index_name, replace_index = FALSE, delete_all = FALSE, ...)`
  computes a selection index by multiplying EBVs by the weights defined
  in
  [`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md)
  and appends results to the new `ind_index` table. Requires piping
  through `get_table("ind_ebv")` (the only `add_*()` function that takes
  `ind_ebv` as its starting table, not `ind_meta`). Each run increments
  an `index_number` column per individual, mirroring `ind_ebv`’s
  `eval_number` pattern. `replace_index = TRUE` clears prior runs for
  the named index; `delete_all = TRUE` clears the entire `ind_index`
  table. Issues a warning when no filter is applied (auto-selects the
  latest `eval_number` per individual per trait). Errors immediately if
  any individual is missing an EBV for any required index trait.

- Two new database tables created by `initialize_genome()`:

  - `index_meta (index_name, trait_name, index_wt)` — index definitions;
    supports user-defined extra columns via `define_index(...)`
  - `ind_index (id_ind, index_name, index_number, index_value)` —
    computed index values; supports user-defined extra columns via
    `add_index(...)`

## tidybreed 0.10.0 (2026-05-01)

### Breaking changes

- `ind_tbv` no longer has a `date_calc` column. The column was removed
  from the schema, the
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
  function signature (`date_calc` parameter dropped), and the internal
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  TBV write path. Users who need a date can add a custom column via
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md).

- `ind_ebv` no longer has a `date_calc` column. The
  [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
  `date_calc` parameter has been removed. The
  [`parse_blupf90_solutions()`](https://austin-putz.github.io/tidybreed/reference/parse_blupf90_solutions.md)
  internal function no longer accepts `date_calc`; its last argument is
  now `eval_nums` (named integer vector).

- `ind_ebv` primary key changed from `(id_ind, trait_name, model)` to
  `(id_ind, trait_name, model, eval_number)`. Existing databases are
  automatically migrated: `date_calc` is dropped and
  `eval_number INTEGER` is added (defaulting to 1 for any pre-existing
  rows).

### New features

- `ind_ebv` gains an `eval_number INTEGER` column. Each call to
  [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
  for a given trait increments the counter by 1 (global per
  `trait_name`, across all models). The latest evaluation can be found
  with `ORDER BY eval_number DESC LIMIT 1`.

- [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
  gains two new logical parameters:

  - `replace_trait = TRUE` — deletes all `ind_ebv` rows for the traits
    being added before inserting; new rows receive `eval_number = 1`.
  - `delete_all = TRUE` — deletes **all** rows in `ind_ebv` before
    inserting; new rows receive `eval_number = 1`. Takes precedence over
    `replace_trait`.

- `add_ebv(..., parent_avg = TRUE)` now queries parent EBVs using
  `MAX(eval_number)` per parent per `(trait_name, model)`. A warning is
  printed when any parent has more than one evaluation for a given
  trait, noting which `eval_number` is used.

## tidybreed 0.9.6 (2026-05-01)

### Bug fixes

- [`upsert_ind_tbv()`](https://austin-putz.github.io/tidybreed/reference/upsert_ind_tbv.md)
  / `upsert_ind_ebv()`: replaced the DELETE + `dbWriteTable` pattern
  with a DuckDB-native `INSERT … ON CONFLICT DO UPDATE SET` UPSERT.
  Previously, calling
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
  (or
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md))
  after a prior
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
  that had written user-defined extra columns (e.g. `rep`) would
  silently reset those columns to `NULL` in existing rows for other
  traits or on re-runs without the extra-column argument. The new UPSERT
  only updates columns that are explicitly present in the incoming data
  frame, leaving all other columns in existing rows untouched.

## tidybreed 0.9.5 (2026-05-01)

### Bug fixes

- [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md):
  removed the `existing_ids` skip guard that prevented re-computing TBVs
  for individuals already in `ind_tbv`. The guard was redundant with the
  DELETE + INSERT logic inside
  [`upsert_ind_tbv()`](https://austin-putz.github.io/tidybreed/reference/upsert_ind_tbv.md)
  and blocked custom column updates (e.g. `rep`) when looping over
  replicates.

## tidybreed 0.9.4 (2026-04-29)

### New features

- [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  gains `.set_default` parameter to create SQL DEFAULT constraints on
  new columns. When `.set_default = TRUE`, future INSERT operations from
  [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md),
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
  and other `add_*()` functions automatically use the default value when
  the column is not explicitly specified. This enables easy schema
  pre-declaration and consistent metadata across individuals.

## tidybreed 0.9.2 (2026-04-29)

### Documentation

- Added `vignettes/tidybreed-introduction.Rmd` — a placeholder
  introduction vignette directing users to the QMD-based quick-start in
  `vignettes/qmd/` until the package reaches v1.0.0.

## tidybreed 0.9.1 (2026-04-28)

### Internal

- `genome_haplotype` and `genome_genotype` locus columns now use
  `UTINYINT` instead of `INTEGER`, reducing per-locus memory from 4
  bytes to 1 byte (4× reduction). `UTINYINT` range (0–255) covers
  biallelic haplotypes (0/1), diploid genotypes (0/1/2), and polyploid
  up to 8n or 16n.

## tidybreed 0.9.0 (2026-04-27)

### New Features

- **Eager table initialization**: `initialize_genome()` now creates all
  core database tables upfront (`ind_meta`, `ind_phenotype`, `ind_tbv`,
  `ind_ebv`, `trait_meta`, `trait_effects`, `trait_effect_cov`,
  `trait_random_effects`). Users can call
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  and
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  on any table immediately after `initialize_genome()`, before any data
  has been added.

- **Custom field forwarding in `add_*` functions**:
  [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md),
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md),
  and
  [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
  now accept `...` for custom column values that are written atomically
  with the new rows. Column types are inferred from the R type (`0L` →
  INTEGER, `0` → DOUBLE, `"text"` → VARCHAR, `TRUE` → BOOLEAN). Scalars
  are broadcast; vectors must match the inserted row count. Reserved
  column names are blocked.

  ``` r

  # Founders with custom columns in a single call
  pop <- pop |>
    add_founders(n_males = 10, n_females = 100, line_name = "A",
                 gen = 0L, farm = "Iowa")
  ```

- **[`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  on empty tables**: calling
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  on an empty table now creates the column schema (via
  `ALTER TABLE ADD COLUMN`) instead of warning and returning early. This
  enables pre-declaring typed column schemas using typed NAs before any
  rows exist:

  ``` r

  pop <- pop |>
    get_table("ind_meta") |>
    mutate_table(gen = NA_integer_, farm = NA_character_)
  ```

- **[`prepare_extra_cols()`](https://austin-putz.github.io/tidybreed/reference/prepare_extra_cols.md)
  internal helper** (in `R/sql_utils.R`): shared validation +
  type-inference + ALTER TABLE logic used by all `add_*` functions with
  `...` support.

- **`TABLE_RESERVED_COLS`** and **`TABLE_PRIMARY_KEYS`** in
  `R/sql_utils.R` expanded to cover `ind_tbv`, `ind_ebv`, `trait_meta`,
  `trait_effects`, and `trait_effect_cov`.

------------------------------------------------------------------------

## tidybreed 0.8.2 (2026-04-24)

### Bug Fixes

- **[`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
  /
  [`parse_blupf90_solutions()`](https://austin-putz.github.io/tidybreed/reference/parse_blupf90_solutions.md)**:
  fixed two bugs when reading back blupf90+ solutions with
  `OPTION origID` active:
  - The parser was looking for a file named `solutions`, but blupf90+
    writes `solutions.orig` when `OPTION origID` is present. The file
    path is now correct.
  - `original_id` values that look numeric (e.g. `10`, `20`) were parsed
    as integers by `read.table`, causing the `%in% all_ped_ids`
    comparison against a character vector to silently return zero
    matches. The column is now coerced to character immediately after
    parsing.
- Added `tests/testthat/test-parse_blupf90_solutions.R` covering:
  correct parsing of the aligned/header format, numeric-looking IDs,
  missing-file error, and no-match warning.

------------------------------------------------------------------------

## tidybreed 0.8.1 (2026-04-24)

### Bug Fixes

- **[`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
  /
  [`write_renum_par()`](https://austin-putz.github.io/tidybreed/reference/write_renum_par.md)**:
  fixed three errors in the generated `renum.par` parameter file for
  BLUPF90:
  - The intercept (mu) effect was declared as `cross` (class variable)
    instead of `cov` (covariate). The mu column is always 1 for every
    row, so it must be `cov`.
  - The animal random effect was missing its own `EFFECT` block before
    `RANDOM animal`. Without it, `renumf90` was associating the last
    fixed-effect `EFFECT` block (e.g. sex at column 3) with the animal
    random effect, assigning the wrong column number instead of column 2
    (`id_ind`).
  - As a consequence of the missing animal `EFFECT` block, the
    fixed-effect column for sex was not written correctly. With the
    animal `EFFECT` block restored at column 2, the sex `EFFECT` block
    at column 3 now lands correctly.

------------------------------------------------------------------------

## tidybreed 0.8.0 (2026-04-24)

### New Features

- **[`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)**
  — new function to populate the `ind_ebv` table with estimated breeding
  values. Two modes:
  - `software = "blupf90"`: runs `renumf90` + `blupf90+` from PATH.
    Builds pedigree, data, and (optionally) genotype files automatically
    from the database. Parameter file (`renum.par`) is auto-generated
    from `trait_effects` and `trait_effect_cov`. All input/output files
    are written to a timestamped subfolder of `run_dir` for full
    reproducibility. Supports BLUP and ssGBLUP (when `chip_name` is
    supplied). After the run, solutions are parsed and EBVs are stored
    in `ind_ebv`.
  - `parent_avg = TRUE`: computes EBVs as the simple average of parent
    EBVs already in `ind_ebv` for the given `model`. Returns `NA` (with
    a warning) for animals whose parents lack EBVs.
- [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
  accepts a `tidybreed_table` (from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md) +
  optional
  [`filter()`](https://dplyr.tidyverse.org/reference/filter.html)) as
  its first argument, following the same calling convention as
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  and
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md).

------------------------------------------------------------------------

## tidybreed 0.7.2 (2026-04-24)

### Breaking Changes

- **`ind_phenotype` schema simplified**: removed `env`, `rep`, and
  `date_measured` columns. These are replaced by a single
  `pheno_number INTEGER` column that auto-increments per individual ×
  trait (1 = first record, 2 = second, etc.). Users can add their own
  columns via `get_table("ind_phenotype") |> mutate_table(...)`.

- **[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  parameters removed**: `env`, `rep`, and `date_measured` arguments have
  been dropped from the function signature. Remove them from any
  existing calls.

- **`ind_phenotype` now registered for
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)**:
  primary key (`id_record`) and reserved columns (`id_record`, `id_ind`,
  `trait_name`, `value`, `pheno_number`) are registered, enabling
  filtered column additions via
  `get_table("ind_phenotype") |> filter(...) |> mutate_table(my_col = ...)`.

------------------------------------------------------------------------

## tidybreed 0.7.1 (2026-04-24)

### Bug Fixes

- [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
  now checks for existing TBV records before computing. Individuals that
  already have a TBV for a requested trait are skipped with an
  informative message rather than silently overwritten.

## tidybreed 0.7.0 (2026-04-24)

### Breaking Changes

- **`add_effect_cov_matrix(pop, effect_name, cov_matrix, trait_names, tol)`**
  — new unified function for storing variance/covariance matrices. Use
  `effect_name = "gen_add"` for additive genetic (co)variances and
  `effect_name = "residual"` for residual (co)variances. Any named
  random effect (e.g. `"litter"`, `"dam"`) is also supported. Replaces
  `set_residual_cov()` and `set_random_effect_cov()`.

- **`set_residual_cov()` removed** — use
  `add_effect_cov_matrix(pop, "residual", R)` instead.

- **`set_random_effect_cov()` removed** — use
  `add_effect_cov_matrix(pop, effect_name, R)` instead.

- **`trait_meta` columns removed**: `target_add_var` and `residual_var`
  are no longer stored in `trait_meta`. Supply them via
  `add_trait(target_add_var = ...)` / `add_trait(residual_var = ...)`
  (values are written to `trait_effect_cov`) or call
  `add_effect_cov_matrix()` directly.

- **`trait_effects` column removed**: `variance` column removed from
  `trait_effects`. Random effect variances are now stored exclusively in
  `trait_effect_cov`.

- **`set_qtl_effects_multi()` — `G` parameter is now optional** (default
  `NULL`). If omitted, the additive genetic covariance matrix is read
  from `trait_effect_cov` (stored via
  `add_effect_cov_matrix("gen_add", ...)`).

- **`add_effect_random()` — `variance` parameter is now optional**
  (default `NULL`). If a diagonal entry for `(effect_name, trait_name)`
  already exists in `trait_effect_cov`, it is used automatically. Only
  required when no stored value exists.

### New Tables

- **`trait_effect_cov`** — unified variance/covariance table replacing
  `trait_residual_cov` and `trait_random_effect_cov`. Schema:
  `(effect_name VARCHAR, trait_1 VARCHAR, trait_2 VARCHAR, cov DOUBLE)`.
  Both `(i,j)` and `(j,i)` pairs stored for symmetric lookup.

### Bug Fixes / Internal Changes

- All writes to `trait_effect_cov` use
  [`DBI::dbExecute()`](https://dbi.r-dbi.org/reference/dbExecute.html)
  with raw SQL instead of
  [`DBI::dbWriteTable()`](https://dbi.r-dbi.org/reference/dbWriteTable.html)
  to avoid consuming R’s RNG state (DuckDB’s `dbWriteTable` internally
  touches the RNG, which would shift simulation reproducibility).

------------------------------------------------------------------------

## tidybreed 0.6.0 (2026-04-23)

### New Features

- **`add_effect_int(pop, trait_name, mean)`** — sets the intercept
  (`target_add_mean`) for a trait; convenience alternative to setting it
  at `add_trait()` time.

- **`add_effect_fixed_class(pop, trait_name, effect_name, source_column, levels, source_table, overwrite)`**
  — discrete fixed effect. Maps levels of a grouping column to numeric
  shifts. Errors (not silent 0) at phenotyping time if an individual’s
  level is not in `levels`. Replaces the `effect_class = "fixed"` path
  of `add_trait_covariate()`.

- **`add_effect_fixed_cov(pop, trait_name, effect_name, source_column, slope, center, source_table, overwrite)`**
  — continuous covariate regression term. Contribution =
  `slope * (x - center)`. When `center = NULL` the mean of
  `source_column` is computed from the current `source_table` and stored
  for reproducibility.

- **`add_effect_random(pop, trait_name, effect_name, source_column, variance, distribution, source_table, overwrite)`**
  — random group effect. Drawn values are now persisted in a new
  `trait_random_effects` table so the same group receives the same shift
  on repeated calls to
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
  without requiring a fixed `seed`.

- **`set_random_effect_cov(pop, effect_name, traits, R, tol, overwrite)`**
  — stores a covariance matrix enabling joint MVN draws of a named
  random effect across multiple traits. Analogous to
  `set_residual_cov()`.

- All `add_effect_*()` functions accept a `source_table` parameter
  (default `"ind_meta"`). Any database table with an `id_ind` column can
  serve as the source for effect levels, enabling future
  repeated-measures and multi-table workflows.

### Schema Changes

- `trait_effects` gains three nullable columns: `source_table VARCHAR`,
  `slope DOUBLE`, `center DOUBLE`. Existing rows (written by
  `add_trait_covariate()`) are backward-compatible — `source_table`
  defaults to `"ind_meta"` when `NULL`.

- Two new tables: `trait_random_effects` (stores drawn group-level
  random effect values) and `trait_random_effect_cov` (covariance
  structure for correlated random effects across traits).

### Deprecations

- `add_trait_covariate()` now emits a deprecation warning. It continues
  to work for backward compatibility but users should migrate to
  `add_effect_fixed_class()` or `add_effect_random()`.

------------------------------------------------------------------------

## tidybreed 0.5.0 (2026-04-23)

### Breaking Changes

- `add_trait()`: parameter `mean` renamed to `target_add_mean`. The
  column in `trait_meta` is also renamed. Existing serialized databases
  with the old column name will not load correctly — re-initialize when
  upgrading.

- `set_qtl_effects()` and `set_qtl_effects_multi()`: first argument
  renamed from `pop` to `x` and now accepts either a `tidybreed_pop` or
  a `tidybreed_table` (for `base = "current_pop"` to specify which
  individuals define the base).

### New Features

- `set_qtl_effects()` / `set_qtl_effects_multi()`: new `base` parameter
  (`"founder_haplotypes"` default, or `"current_pop"`). Controls which
  allele frequencies are used for effect scaling and TBV centering.

- Effect scaling now uses the **Falconer formula**
  (`V_A = Σ 2·p·(1−p)·α²`) instead of the realized `var(G %*% alpha)`.
  This gives a theoretically grounded `target_add_var` guarantee rather
  than an empirical one.

- `set_qtl_effects()` writes a `base_allele_freq_{trait}` column to
  `genome_meta` recording which allele frequencies were used. This
  column is read by
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
  and
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  to center TBVs: `TBV_i = (G_i − 2·p_base) · α`, ensuring `E[TBV] ≈ 0`
  for the base population.

- `base = "founder_haplotypes"` computes allele frequencies from the
  actual rows in the `founder_haplotypes` table (the pool used to sample
  founders), not from the stored theoretical `founder_allele_freq`
  column.

------------------------------------------------------------------------

## tidybreed 0.4.3 (2026-04-23)

### Internal / Cleanup

- Action functions
  ([`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md),
  [`add_genotypes()`](https://austin-putz.github.io/tidybreed/reference/add_genotypes.md),
  [`extract_genotypes()`](https://austin-putz.github.io/tidybreed/reference/extract_genotypes.md))
  fully migrated to the `tidybreed_table`-first calling convention
  announced in 0.4.1. Removed the internal `resolve_pending_filter()`
  helper; each function now collects and extracts `id_ind` directly from
  the `tidybreed_table`.

- [`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md)
  and `define_qtl()` now call
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  directly instead of relying on the removed `mutate_genome_meta()`.

- [`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
  now uses
  [`validate_sql_identifier()`](https://austin-putz.github.io/tidybreed/reference/validate_sql_identifier.md)
  for extra column validation, consistent with the rest of the package.

- Source files for deleted functions (`R/mutate_genome_meta.R`,
  `R/mutate_ind_meta.R`, `R/add_ebv.R`) and their associated man pages
  and tests are now physically removed from the repository.

- Expanded roxygen documentation: new `man/` pages for many exported and
  internal functions, and updated
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  /
  [`infer_duckdb_type()`](https://austin-putz.github.io/tidybreed/reference/infer_duckdb_type.md)
  docs.

- Test suite updated throughout to use the new
  `get_table() |> filter() |> action()` pattern.

## tidybreed 0.4.2 (2026-04-23)

### New Features

- [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  chains now support the full `slice_*` family:
  [`slice_max()`](https://dplyr.tidyverse.org/reference/slice.html),
  [`slice_min()`](https://dplyr.tidyverse.org/reference/slice.html),
  [`slice_head()`](https://dplyr.tidyverse.org/reference/slice.html),
  [`slice_tail()`](https://dplyr.tidyverse.org/reference/slice.html),
  and
  [`slice_sample()`](https://dplyr.tidyverse.org/reference/slice.html).
  These can be used directly on a `tidybreed_table` before passing to
  action functions such as
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md).

## tidybreed 0.4.1 (2026-04-23)

### Breaking Changes

- [`filter()`](https://dplyr.tidyverse.org/reference/filter.html) can no
  longer be called directly on a `tidybreed_pop` object. You must now
  call
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  first to identify which table to filter. This applies to
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md),
  [`add_genotypes()`](https://austin-putz.github.io/tidybreed/reference/add_genotypes.md),
  and
  [`extract_genotypes()`](https://austin-putz.github.io/tidybreed/reference/extract_genotypes.md).

  ``` r

  # Old (no longer works):
  pop |> dplyr::filter(sex == "F") |> add_phenotype("ADG")

  # New (required):
  pop |> get_table("ind_meta") |> dplyr::filter(sex == "F") |> add_phenotype("ADG")

  # Any table with id_ind now works:
  pop |> get_table("ind_phenotype") |> dplyr::filter(value > 500) |> add_phenotype("ADG2")
  ```

- [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
  [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md),
  [`add_genotypes()`](https://austin-putz.github.io/tidybreed/reference/add_genotypes.md),
  and
  [`extract_genotypes()`](https://austin-putz.github.io/tidybreed/reference/extract_genotypes.md)
  now require a `tidybreed_table` as the first argument (from
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)),
  not a `tidybreed_pop`.

- [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
  has been removed. A proper evaluation runner will replace it in a
  future version.

## tidybreed 0.4.0 (2026-04-23)

### Breaking Changes

- `mutate_ind_meta()` and `mutate_genome_meta()` have been **removed**.
  Use the new generic
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  instead:

  ``` r

  # Old:
  pop <- mutate_ind_meta(pop, gen = 1L)
  # New:
  pop <- pop |> get_table("ind_meta") |> mutate_table(gen = 1L)
  ```

### New Functions

- `mutate_table(tbl_obj, ...)` — generic column add/update for any
  table. Chain after
  [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  (and optionally
  [`filter()`](https://dplyr.tidyverse.org/reference/filter.html)) to
  add new columns or update existing ones. Returns `pop` invisibly.
  Supports scalar and vector values, type inference, reserved-column
  blocking, and informative messages.

  ``` r

  # All rows:
  pop <- pop |> get_table("ind_meta") |> mutate_table(gen = 1L)

  # Filtered rows only (females get NULL for a new column):
  pop <- pop |>
    get_table("ind_meta") |>
    filter(sex == "M") |>
    mutate_table(gen = 2L)
  ```

### Changes

- [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  now returns a `tidybreed_table` S3 object instead of a raw
  `tbl_duckdb_connection`. Backward-compatible via
  `collect.tidybreed_table`, `filter.tidybreed_table`,
  `select.tidybreed_table`, `arrange.tidybreed_table`,
  `pull.tidybreed_table`, and `count.tidybreed_table` S3 methods — all
  existing `get_table(...) |> dplyr::filter(...) |> dplyr::collect()`
  patterns continue to work unchanged.
- [`infer_duckdb_type()`](https://austin-putz.github.io/tidybreed/reference/infer_duckdb_type.md)
  moved to `R/sql_utils.R` and is now a shared internal utility
  available to all mutation helpers.
- [`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md),
  `define_qtl()`, `set_qtl_effects()`, and `set_qtl_effects_multi()`
  updated to call
  [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  internally.

------------------------------------------------------------------------

## tidybreed 0.3.0 (2026-04-21)

### New Functions

- `add_genotypes(pop, chip_name)` — marks a filtered subset of animals
  as genotyped on a named SNP chip by writing a `has_<chip_name>`
  BOOLEAN column to `ind_meta`. Follows the same
  [`filter()`](https://dplyr.tidyverse.org/reference/filter.html) -\>
  action pipe pattern as
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md).
  Operation is additive: animals already marked TRUE remain TRUE across
  multiple calls. Chip must exist in `genome_meta` (via
  [`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md))
  before calling.
- `extract_genotypes(pop, chip_name)` — returns a tibble of genotypes
  (0/1/2 encoding) for animals marked as genotyped, restricted to chip
  loci. The returned set is the intersection of animals with
  `has_<chip_name> == TRUE`, any pending
  [`filter()`](https://dplyr.tidyverse.org/reference/filter.html)
  predicates, and loci with `is_<chip_name> == TRUE`. Intended for use
  immediately before GBLUP/GWAS evaluation.

------------------------------------------------------------------------

## tidybreed 0.2.2 (2026-04-20)

### Performance

- [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md):
  eliminated nested R loop (O(n_founders × n_loci) iterations) that
  caused multi-minute runtimes for large genomes. Haplotype and genotype
  frames are now built via vectorized matrix indexing and addition in C,
  reducing frame construction from minutes to ~0.18 s regardless of
  genome size. Also switched genome table writes to `duckdb_register` +
  `INSERT SELECT` for a further ~2× speedup on the DB write step.
  Typical runtimes: 0.38 s (1k loci) → 2.7 s (5k loci) → 8 s (10k loci)
  for 215 founders.

## tidybreed 0.2.1 (2026-04-20)

### API

- Renamed `name` argument to `trait_name` in `add_trait()` and
  `add_trait_simple()` for consistency
- Renamed `trait` argument to `trait_name` in `define_qtl()` and
  `set_qtl_effects()`
- Renamed `traits` argument to `trait_names` in
  `set_qtl_effects_multi()`
- Renamed `name` argument to `chip_name` in
  [`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md)

## tidybreed 0.2.0 (2026-04-20)

### New: Trait Architecture & Phenotype Simulation

Trait definition, QTL selection, additive-effect sampling (single and
correlated multi-trait), fixed and random covariates, phenotype
generation for continuous / count / binary / categorical traits,
imprinting support, and storage of true and estimated breeding values.

#### New exported functions

- `add_trait()` — insert a row in `trait_meta` with target variance
  components, trait type (continuous / count / binary / categorical),
  expression rules (sex-limited, parent-of-origin), and index/economic
  weights. Creates the six trait-layer tables on first call.
- `define_qtl()` — mirror of
  [`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md)
  for QTL loci; writes an `is_QTL_{trait}` BOOLEAN column in
  `genome_meta`. Reuses the same six selection methods (by count +
  `random` / `even` / `chromosome_even`, by logical vector, by locus
  ids, by locus names).
- `set_qtl_effects()` — write the `add_{trait}` additive-effect column.
  Supports manual effects or sampled effects (`normal` / `gamma`) with
  automatic rescaling to hit `trait_meta$target_add_var`.
- `set_qtl_effects_multi()` — draw correlated additive effects across
  multiple traits from `MVN(0, G)` via
  [`MASS::mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html); supports
  `"shared"` (pleiotropy) and `"union"` strategies.
- `set_residual_cov()` — store a residual covariance matrix `R` across
  traits in a new `trait_residual_cov` table. Consumed by
  [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  when multiple traits share the same filtered subset.
- `add_trait_covariate()` — append fixed or random covariate rows to a
  `trait_effects` table. Fixed-effect levels are serialised to a
  JSON-style VARCHAR; random effects store distribution + variance.
- [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  — the workhorse. Generates phenotypes for a subset of individuals for
  one or more traits and writes rows to `ind_phenotype`. Also computes
  and stores the underlying TBV in `ind_tbv`. Joint MVN residual draws
  when multiple traits share the subset and `R` is stored.
- [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
  — compute and store TBV without generating phenotype records.
- [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
  — ingest externally computed estimated breeding values into `ind_ebv`,
  tagged with a user-supplied model label.
- `add_trait_simple()` — one-shot wrapper chaining `add_trait()` +
  `define_qtl()` + `set_qtl_effects()`.

#### filter() on tidybreed_pop

New S3 method `filter.tidybreed_pop()` stashes dplyr predicates on the
population object. The next
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
/
[`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
call applies and clears them. Multiple
[`filter()`](https://dplyr.tidyverse.org/reference/filter.html) calls
stack with AND semantics.

``` r

pop |>
  dplyr::filter(sex == "F", gen == 1L) |>
  add_phenotype("ADG")
```

#### New tables

`trait_meta`, `trait_effects`, `trait_residual_cov`, `ind_phenotype`,
`ind_tbv`, `ind_ebv`. Also new columns written to `genome_meta`:
`is_QTL_{trait}` (BOOLEAN) and `add_{trait}` (DOUBLE) per trait.

#### Other

- Added `MASS` to Imports for multivariate normal sampling. MASS is a
  Recommended R package that ships with base R distributions — no extra
  install burden.

------------------------------------------------------------------------

## tidybreed 0.1.0 (2026-04-17)

### New Functions

- [`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
  — core mating function. User supplies a `matings` tibble (one row per
  offspring) with required columns `id_parent_1`, `id_parent_2`, `sex`,
  and `line`. Gametes are produced via chromosomal crossover simulation
  (Haldane map: crossovers per chromosome ~ Poisson(chr_len_Mb / 100)).
  New `ind_meta`, `genome_haplotype`, and `genome_genotype` rows are
  written for all offspring. Animal-breeder aliases `id_sire` / `id_dam`
  are accepted. Any extra columns in `matings` (e.g. `gen = 2L`) are
  validated and written to `ind_meta`, with automatic `ALTER TABLE` if
  the column is new.
- [`make_gamete()`](https://austin-putz.github.io/tidybreed/reference/make_gamete.md)
  — internal recombination helper (`R/recombination_helpers.R`). Not
  exported; used by
  [`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md).

------------------------------------------------------------------------

## tidybreed 0.0.3 (2026-04-17)

### Breaking Changes

- `ind_meta`, `genome_haplotype`, `genome_genotype`: column `ind_id`
  renamed to `id_ind`; `parent_1` → `id_parent_1`; `parent_2` →
  `id_parent_2`. The `id_` prefix groups all identifier columns together
  for clarity.

**Migration for existing databases:**

``` r

conn <- DBI::dbConnect(duckdb::duckdb(), "your_database.duckdb")
DBI::dbExecute(conn, "ALTER TABLE ind_meta RENAME COLUMN ind_id TO id_ind")
DBI::dbExecute(conn, "ALTER TABLE ind_meta RENAME COLUMN parent_1 TO id_parent_1")
DBI::dbExecute(conn, "ALTER TABLE ind_meta RENAME COLUMN parent_2 TO id_parent_2")
DBI::dbExecute(conn, "ALTER TABLE genome_haplotype RENAME COLUMN ind_id TO id_ind")
DBI::dbExecute(conn, "ALTER TABLE genome_genotype RENAME COLUMN ind_id TO id_ind")
DBI::dbDisconnect(conn)
```

### Bug Fixes

- [`infer_duckdb_type()`](https://austin-putz.github.io/tidybreed/reference/infer_duckdb_type.md):
  bare `NA` (untyped logical) now warns and returns `VARCHAR` instead of
  silently returning `BOOLEAN`
- `mutate_ind_meta()`: type check now runs before length check, so
  passing an unsupported type (e.g. a `list`) produces “Unsupported
  type” rather than a misleading length-mismatch error
- Fixed pre-existing test bugs: `expect_error(expr, NA)` inverted
  assertion in `test-add_founders.R`; null `db_conn` in
  `test-mutate_genome_meta.R` now uses a real in-memory connection

## tidybreed 0.0.2 (2026-04-15)

### Bug Fixes

- [`infer_duckdb_type()`](https://austin-putz.github.io/tidybreed/reference/infer_duckdb_type.md):
  fixed incorrect VARCHAR inference when a typed vector starts with `NA`
  (e.g. `c(NA_real_, rnorm(n))`). Type is now determined from the R
  class of the full vector before inspecting individual elements, so
  `NA_real_` → `DOUBLE`, `NA_integer_` → `INTEGER`, `NA` → `BOOLEAN`. A
  bare all-`NA` vector (class `logical`) now warns with a message
  pointing users to typed NA constants.

## tidybreed 0.0.1 (2026-04-14)

### Breaking Changes

#### Terminology: “population” → “line”

- [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
  parameter renamed: `pop_name` → `line_name`
- `ind_meta` table column renamed: `population` → `line`

**Rationale:** “population” refers to one complete genome build from
`initialize_genome()`. “Line” refers to distinct genetic lines within
that population (e.g., line A selected for growth, line B for egg
quality).

**Migration for existing databases:**

``` r

conn <- DBI::dbConnect(duckdb::duckdb(), "your_database.duckdb")
DBI::dbExecute(conn, "ALTER TABLE ind_meta RENAME COLUMN population TO line")
DBI::dbDisconnect(conn)
```

### New Functions

- `mutate_genome_meta()` — add or update user columns in `genome_meta`;
  reserved columns (`locus_id`, `locus_name`, `chr`, `chr_name`,
  `pos_Mb`) are blocked
- [`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md)
  — convenience wrapper that marks loci as members of a named SNP chip
  by writing a `BOOLEAN` column `is_{chip_name}` to `genome_meta`;
  supports `"random"`, `"even"`, and `"chr_even"` selection methods

### Other Changes

- Removed stale design and planning documents (DESIGN.md, QUICKSTART.md,
  IMPLEMENTATION_STATUS.md, IMPLEMENTATION_SUMMARY.md, TODO_finalize.md,
  MUTATE_IND_META_IMPLEMENTATION.md)
- Moved manual test scripts to `tests/` directory
- Added CLAUDE.md with developer and AI context for the project

------------------------------------------------------------------------

## tidybreed 0.0.0

- Initial package setup and framework

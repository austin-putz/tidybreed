# Genome Position Coordinates & Genetic Map

**Status**: ✅ IMPLEMENTED (v0.50.0, branch `feat/position-coordinates-genome-map`).
All 4 steps landed; testthat suite green (0 fail / 1227 pass, incl. new
`test-genome_map.R`); `R CMD check` clean of new issues. Revised after two review
rounds + BIGINT-mechanism resolved (Option B).
**Relationship**: **Prerequisite for** `plans/optimize_add_offspring.md`. That plan's
Morgans-native recombination kernel reads the genetic map produced here. This is a
**separate, self-contained project**: it changes the physical/genetic coordinate
storage and every current reader/writer, lands and is verified on its own, and
does **not** touch the speed/C++/RNG work. Do this first.

---

## Pre-1.0.0: no cross-version compatibility

Per CLAUDE.md → "Development Status: Pre-1.0.0 — Break Freely", this migration owes
**nothing** to the current implementation's output. We are **not** preserving
byte-identical founder/offspring haplotypes across the change, not deriving formulas
to match old values, and not capturing "golden-from-old" fixtures. The only
reproducibility requirement is **forward**: a given seed reproduces on repeated runs
of the *new* code. This simplifies several critique items (see below): the genetic
position is derived straight from the sampled `genome_meta.pos_bp`, with no
old-`pos_Mb` shadow value.

## Resolution of the review critique (2026-07-05)

The critique in `plans/position_coordinates_critique.md` was reviewed against code.
Where a point existed only to preserve old behavior, it is **dropped** under the
pre-1.0.0 stance; the structural points are kept:

- **#1 (moot under pre-1.0.0):** the critique objected that `pos_cM` derived from
  rounded `pos_bp` wouldn't match the old unrounded `pos_Mb`. We don't need the
  match. **`pos_bp` is the single sampled physical coordinate; `pos_cM` derives
  directly from it** (`pos_cM = pos_bp/1e6 · cM_per_Mb`). No `pos_Mb_old` shadow —
  the draft's original formula stands (Step 1).
- **#2 (accepted — but only the BIGINT part):** `founder_haplotype_helpers.R:187-189`
  does `SELECT * → dbWriteTable(overwrite=TRUE)`, which would demote `pos_bp BIGINT`.
  Replaced with `ALTER TABLE ADD COLUMN` + `UPDATE` **for type preservation**. The
  critique's "must stay RNG-neutral to match the historical writer" reason is
  **void** (we don't match old output); we only require the new path be
  self-consistent and deterministic. Step 2.
- **#3/#12 (accepted):** single `resolve_genome_map()` reader + `validate_genome_map()`
  validator with explicit NULL-key normalization (see "Map helpers").
- **#4 (accepted, scoped):** validate `pos_cM` monotonicity on the **resolved**
  map, in the helper. The critique's "require complete overrides per locus" is
  deferred to the future `define_genetic_map()` writer (v1 writes only the default map).
- **#5 (accepted — your call):** add a `map_name` dimension now (default
  `"default"`), so versions are rows not schema.
- **#6 (accepted):** `genome_map` stores **both** `locus_id` and `locus_name`;
  internal joins/ordering key on `locus_id` (mirrors `ind_haplotype`).
- **#7 (accepted — your call): the `recombines_M`/`_F` split is IN this migration**
  (no longer "recommended/vetoable").
- **#8 (accepted):** the recombination reader fetches parent `ind_meta.line_name`
  alongside sex/ploidy, even though v1 only exercises the default map.
- **#9 (accepted):** `pos_bp` non-null + strictly-increasing for *simulated*
  genomes; nullable only via a future import writer; assert no duplicate bp after rounding.
- **#10 (accepted — your call): `introduced_gen` removal stays in this plan**, with
  an enumerated touchpoint audit (see "Removing `introduced_gen`").
- **#11 (dropped under pre-1.0.0):** no golden fixtures captured from old `HEAD`.
  Replaced by a **forward** reproducibility test (same seed reproduces on the new
  code) — see Testing.
- **#13/#14/#15 (accepted):** document cM chr-length semantics, validate `cM_per_Mb`,
  own the doc/test churn.

### Second review round (tightening) — all accepted

A follow-up critique reviewed this revision. Every point was folded in:

- **R1: residual compat wording removed** — Step 3 no longer says "value-neutral" or
  compares to "old `pos_Mb`"; it states only that the reader consumes the same
  newly-written default map rows.
- **R2 + R14: `CLAUDE.md` is updated dead last**, only after code + tests are green,
  documenting what actually landed — never as a design target.
- **R3: `introduced_gen` rationale trimmed** to a short bullet + touchpoint audit;
  the Stage-5 event-modeling design is out of scope here.
- **R4: explicit `pos_bp` failure mode** — validate before writing; hard error when a
  chromosome is too short to place all loci as positive, strictly-increasing
  integers; no silent `pmax()`/jitter repair.
- **R5: `cM_per_Mb > 0` reconciled with flat maps** — `define_genome()` always writes
  positive `pos_cM`; organelle "no recombination" comes from `recombines_* = FALSE`,
  not a flat cM map.
- **R6: `id_genome_map` uses `next_int_id()`** (the package's next-ID helper, as in
  `genome_effects`), not DB auto-increment.
- **R7: strict map-helper contracts** — see "Map helpers".
- **R8: full logical key everywhere** — `(sex=NULL, line_name=NULL, map_name="default")`.
- **R9: recombination map keyed on the parent's `ind_meta.line_name`** — never
  offspring line or allele `line_origin`.
- **R11: `recombines_M/_F` split is unconditional** — "if decision #7" phrasing
  removed; `define_chr(recombines=)` is the primary shorthand, not a compat alias.
- **R12: nullable `pos_bp` (import) is out of first-stage scope** — `define_genome()`
  enforces non-null.
- **R13: no `pos_cM` in `genome_meta`** — physical in `genome_meta`, genetic in
  `genome_map`, join via `locus_id`.

### Post-implementation critique (2026-07-06, Codex) — all addressed

A review of the landed code found 8 items; all fixed in this branch:

- **C1 (High): per-gamete map resolution was not actually wired** — `add_offspring()`
  resolved one default map for all gametes. **Fixed:** it now caches a descriptor
  set per unique `(parent_sex, parent_line)` and each gamete uses its producing
  parent's resolved map. Behavior-neutral in v1; a sex/line `genome_map` override
  now changes recombination (new regression test).
- **C2 (Med): `CLAUDE.md` attributed `cM_per_Mb` to deprecated `initialize_genome()`**
  — **Fixed:** section reframed to `open_pop() |> define_genome()` as the real
  surface; `initialize_genome()` noted as a deprecated wrapper without `cM_per_Mb`.
- **C3 (Med): stale recombination docs** — **Fixed:** `add_offspring()` roxygen now
  says `Poisson(chr_len_cM/100)`; `pass_through_gamete()` doc uses `recombines_M/_F`;
  the "BYTE-IDENTICAL to pre-Stage-4" comment removed.
- **C4/C5 (Med): validator not strict enough** — **Fixed:** `validate_genome_map()`
  now also checks `pos_cM` non-null/finite and `id_genome_map` uniqueness; new tests
  cover invalid sex, NULL `map_name`, non-finite `pos_cM`, duplicate id, line
  precedence, missing-locus resolve error, and `genome_map` schema-metadata.
- **C6 (Low): `chr_len_Mb` unvalidated** — **Fixed:** now required finite + positive.
- **C7 (Low): double `define_genome()` partial clobber** — **Fixed:** preflight error
  if `genome_meta` is already non-empty.
- **C8 (Low): crossover-length semantics** — documented in
  `plans/optimize_add_offspring.md` (the `max(pos_cM)` marker-span convention and
  its `ind_crossover` implication).

Post-fix: full testthat suite green (0 fail / 1236 pass).

---

## Why

`genome_meta.pos_Mb DOUBLE` conflates two genuinely different genomics
coordinates, and it stores the genetic one in the wrong *shape*:

- **Physical position** (base pairs) is a single, integer, per-locus property →
  belongs in `genome_meta`.
- **Genetic position** (centiMorgans) is **multi-dimensional** — it varies by
  **sex** (sex-differential maps; achiasmy in *Drosophila* males / Lepidoptera
  females) and by **line/population** (breeds have different maps), and may have
  **versions**. That is the same shape as QTL effects, which this package
  already stores in a **separate long table** (`genome_effects`), *not* as
  columns in `genome_meta`. The genetic map follows the same precedent.

So: physical stays in `genome_meta`; the genetic map moves to its own long table
`genome_map`. Adding a sex-specific, line-specific, or versioned map then becomes
**adding rows, never a schema change** — the whole point.

## Locked decisions

1. **`genome_meta`: `pos_Mb` → `pos_bp`** (physical only). No `pos_cM` column.
2. **`pos_bp` type = `BIGINT`, 1-based.** Same 8 bytes as the old `DOUBLE`;
   bulletproof for any genome size; 1-based matches VCF `POS` / PLINK `.bim`
   col 4 for lossless export. (`UINTEGER` = a later 4-byte optimization.)
   **Mechanism (Option B, resolved — see "Resolved: `pos_bp` BIGINT mechanism"):**
   `genome_meta` is built with an explicit typed `CREATE OR REPLACE TABLE
   (... pos_bp BIGINT)` + `duckdb_register`/`INSERT` (matching the existing
   RNG-safe chr_meta write), **not** `dbWriteTable(overwrite=TRUE)` (which
   infers `DOUBLE`/`INTEGER` from R). The type is a **column** property: once
   declared, DuckDB widens every inserted value to `BIGINT` — user row-adds of a
   plain R `integer` need **no R-side conversion**; only *recreating* the table
   (`overwrite=TRUE`, which the package no longer does to `genome_meta`) can demote it.
3. **New `genome_map` table** holds `pos_cM DOUBLE`, keyed
   `(locus_id, sex, line_name, map_name)`; **`sex IS NULL` = applies to both sexes**,
   **`line_name IS NULL` = applies to all lines** — mirroring `genome_effects`'s
   `line_name IS NULL` = population-wide convention. Stores both `locus_id` (internal
   join/order key) and `locus_name` (denormalized, like `ind_haplotype`).
   `map_name` defaults to `"default"` so map versions are **rows, not schema** (critique #5/#6).
4. **Recombination + founder-LD both key on the genetic map (`genome_map.pos_cM`).**
   Everything distance-driven uses cM, not physical bp.
5. **`define_genome(cM_per_Mb = 1.0)`** (scalar or length-`n_chr`). **`pos_cM`
   derives directly from the sampled `genome_meta.pos_bp`**:
   `pos_cM = pos_bp/1e6 · cM_per_Mb`. `pos_bp` is the single physical source of
   truth; there is no old-`pos_Mb` shadow value and no cross-version parity goal
   (pre-1.0.0). `cM_per_Mb` is validated: numeric, finite, non-missing, `> 0`,
   scalar or length-`n_chr`. **Because `cM_per_Mb > 0` and `pos_bp` is strictly
   increasing, every `define_genome()` map has strictly-increasing `pos_cM` — it
   never produces a flat (`pos_cM = 0`) map** (R5). Whole-chromosome "no
   recombination" (organelles, achiasmy) is expressed by `recombines_* = FALSE`
   (decision #7), not by a flat cM map. A genuinely flat/non-linear map is a future
   `define_genetic_map()` concern, out of scope here.
6. **Map resolution precedence** for a gamete from parent sex `S`, line `L`
   (first match per locus wins), mirroring `genome_effects`'s line-fallback:
   `(sex=S, line=L)` → `(sex=S, line NULL)` → `(sex NULL, line=L)` →
   `(sex NULL, line NULL)`, all within `map_name = "default"`. v1 writes only that
   default map → a single shared map.
7. **LOCKED (accepted in review): `chr_meta.recombines` → `recombines_M` +
   `recombines_F`.** `chr_meta` already carries `copy_mode_M`/`copy_mode_F`, so a
   single `recombines` boolean is the odd one out. Making it per-sex (a) restores
   consistency and (b) gives a **compact** whole-chromosome achiasmy switch
   (e.g. *Drosophila* males: `recombines_M = FALSE` on every chromosome) without
   needing a degenerate sex-specific `genome_map`. `define_chr(recombines = ...)`
   is the **primary shorthand** that sets both columns — a first-class arg in the
   new API, not a compatibility alias (R11).
8. **Remove `introduced_gen` from `genome_meta`.** Unused (nothing reads it; set
   to `NA` for every locus), forces the rejected "generation" concept, and a
   per-locus scalar cannot represent provenance — gene-editing (CRISPR) applies
   the same change to a **set** of individuals at once, and provenance is
   derivable anyway. See "Removing `introduced_gen`".

---

## Schema

### `genome_meta` (modified)

| Column         | Type    | Notes                                          |
|----------------|---------|------------------------------------------------|
| locus_id       | INTEGER | unchanged (PK)                                 |
| locus_name     | VARCHAR | unchanged                                      |
| chr            | INTEGER | unchanged                                      |
| chr_name       | VARCHAR | unchanged                                      |
| **pos_bp**     | BIGINT  | physical position, base pairs, **1-based**     |

**`introduced_gen` is removed** (decision #8) — it was unused Stage-5 scaffolding,
forced the "generation" concept, and (per the CRISPR case) a per-locus scalar
cannot model provenance anyway. See "Removing `introduced_gen`" below.

**Reserved** (`TABLE_RESERVED_COLS$genome_meta`): `locus_id`, `locus_name`,
`chr`, `chr_name`, `pos_bp`. (`pos_Mb`/`pos_cM`/`introduced_gen` are **not** in
`genome_meta`.)

#### Removing `introduced_gen` (decision #8)

Unused Stage-5 scaffolding — `NA` for every locus, nothing reads it. It is both
**derivable** (founding loci = those in `founder_haplotypes`) and the wrong shape
for provenance, so it is a **pure deletion** here, not a redesign. Provenance/event
modeling is a future Stage-5 concern, explicitly out of scope for this coordinate plan.

**Touchpoint audit (R3 — one deletion sweep).** Run `rg introduced_gen`:

- `R/define_genome.R` — stop writing the column.
- `R/schema.R` — remove from the `genome_meta` schema definition.
- `R/sql_utils.R` — remove from `TABLE_RESERVED_COLS$genome_meta`.
- Tests snapshotting `genome_meta` columns (`test-define_genome.R`,
  restore/archive/schema-summary) — update expected column sets.
- Confirm `restore_pop()` / `archive_replicate()` / `summary_pop()` don't reference it.
- `CLAUDE.md` — remove the row + "Reserved" mention **(done last, with all other
  doc updates — see Step 4).**

### `genome_map` (new) — the genetic map

One row per (locus × sex × line) that has a defined genetic position. Populated
by `define_genome()` (a single default map — `(sex=NULL, line_name=NULL,
map_name="default")` — in v1) and, later, by a `define_genetic_map()`-style
function for sex/line-specific maps.

| Column        | Type    | Notes                                                          |
|---------------|---------|----------------------------------------------------------------|
| id_genome_map | INTEGER | Surrogate PK, assigned via `next_int_id(conn, "genome_map", "id_genome_map")` — the package's next-ID helper (as in `genome_effects`/`trait_meta`), **not** DB auto-increment (R6) |
| locus_id      | INTEGER | FK to `genome_meta.locus_id` — **internal join/order key** (critique #6) |
| locus_name    | VARCHAR | FK to `genome_meta.locus_name` — denormalized, like `ind_haplotype` |
| sex           | VARCHAR | `NULL` = both sexes; `'M'`/`'F'` = sex-specific map            |
| line_name     | VARCHAR | `NULL` = all lines; set for line-specific maps                 |
| map_name      | VARCHAR | Map version/identity; default `"default"` (critique #5)        |
| pos_cM        | DOUBLE  | genetic-map position, centiMorgans                             |

**Logical key** `(locus_id, sex, line_name, map_name)` (nullable `sex`/`line_name`,
enforced in R by `validate_genome_map()` with explicit NULL normalization — same
approach as `founder_haplotypes`/`index_meta`, critique #3). Hot joins and the
resolved map order by `locus_id`.
**Reserved**: all columns (managed by map-writing functions).
Register `genome_map` in `TABLE_RESERVED_COLS`, `SYSTEM_TABLES`, schema metadata,
and `pop$tables` (critique #3).

Why a table, not `pos_cM_M`/`pos_cM_F` columns: adding a map dimension (a new
sex, a line, a version) is **rows**, never `ALTER TABLE` — identical to how
`genome_effects` scales across trait/line/effect-type.

### `chr_meta` (modified — decision #7, locked)

`recombines BOOLEAN` → `recombines_M BOOLEAN` + `recombines_F BOOLEAN` (both
default `TRUE`), matching the existing `copy_mode_M`/`copy_mode_F` pattern. This is
part of the schema migration, not optional.

---

## Edge cases & future-proofing — does this cover everything?

The goal is: no schema refactor when these arrive. Coverage of the scenarios
raised in review:

| Scenario | Handled? | How — **without** a future schema change |
|---|---|---|
| **Per-chromosome recombination rate** | ✅ | `cM_per_Mb` per-chr, or arbitrary `genome_map.pos_cM`; chr genetic length = its cM span |
| **Within-chromosome hotspots** | ✅ (engine/gen only) | `genome_map.pos_cM` can be **non-linear** in `pos_bp`; only the *generator* needs enhancing, schema already supports it |
| **Pseudoautosomal region (partial recombination on sex chr)** | ✅ (engine/gen only) | non-PAR loci at equal `pos_cM` ⇒ zero genetic length ⇒ no crossovers there; PAR loci spread in cM. Needs a non-linear map + cM-space crossover draws (the Morgans kernel) |
| **Sex-differential maps** | ✅ | `sex='M'`/`'F'` rows in `genome_map` |
| **Achiasmy (Drosophila ♂, Lepidoptera ♀ — no recombination in one sex)** | ✅ | `chr_meta.recombines_M/_F = FALSE` (compact), or a flat sex-specific `genome_map` |
| **Line/breed-specific maps** | ✅ | `line_name` rows in `genome_map` |
| **Polyploid plants (auto/allo)** | ✅ (position side) | position is **ploidy-independent** — unchanged by `genome_map`/`pos_bp`; polyploid *meiosis* mechanics live in `ind_meta.ploidy`/`chr_meta`/the engine, not here. `BIGINT` covers giant plant genomes |
| **Mitochondria / plastids** | ✅ | `chr_meta.recombines_* = FALSE` + uniparental via `hemi_parent` suppresses recombination regardless of `pos_cM` (organelle loci still get positive `pos_cM` from `define_genome()`, which is simply never used); circular genome → standard linear coordinate |
| **Haplodiploidy (bees/ants/wasps)** | ⚠️ verify | representable via all-chromosome `copy_mode_M = "half"` (haploid males); **add a test** to confirm, no schema change expected |
| **Holocentric chromosomes (Lepidoptera, sedges)** | ✅ | segregation-mechanics concern, position-agnostic |
| **Crossover interference (Kosambi/gamma)** | ✅ (engine) | `pos_cM` is map-function-agnostic; interference is a **recombination-engine** model change (Poisson→gamma), not a schema change |
| **Importing real maps (unknown/duplicate bp)** | 🔜 future writer, **not v1** | a future import/map writer may allow `pos_bp` NULL for unplaced markers; **`define_genome()` in this stage always writes non-null, strictly-increasing `pos_bp`** (R12) |

**Known limitation (not a schema issue): true X–Y PAR *exchange*.** X and Y are
separate chromosomes in `chr_meta`, so recombination *between* the X and Y in the
PAR (they pair only there during ♂ meiosis) is not modeled — each recombines only
with its own homolog. Very niche (few PAR loci); a meiosis-model limitation to
document, not a coordinate-schema one.

**Position conventions to pin now (cheap, avoids export/import pain later):**
- `pos_bp` is **1-based** (VCF/PLINK convention).
- **This stage: `pos_bp` is non-null and per-chromosome strictly-increasing** for
  every genome `define_genome()` creates (assert on write; hard error otherwise —
  see Step 1). Nullable `pos_bp` for imported maps with unplaced markers is a
  future writer's concern, not part of this first-stage contract (R12).
- `pos_cM` **monotonic non-decreasing** within a chromosome (crossover math
  assumes it) — assert on the resolved map (see `resolve_genome_map()`).

---

## Step 1 — `define_genome()` writes `pos_bp` + the default `genome_map`

New arg `cM_per_Mb = 1.0` (scalar or length-`n_chr`, validated per decision #5).
Keep the current interior even-spacing. **`pos_bp` is sampled first as the single
physical coordinate; `pos_cM` derives from it** (no old-`pos_Mb` shadow, no
cross-version parity — pre-1.0.0), per chromosome `i`:

```r
chr_len_bp_i <- round(chr_len_Mb[i] * 1e6)                   # physical span, bp
pos_bp_i     <- round(seq(0, chr_len_bp_i,                   # sampled physical coordinate, bp
                          length.out = n_chr_loci + 2)[2:(n_chr_loci + 1)])
pos_cM_i     <- (pos_bp_i / 1e6) * cM_per_Mb[i]              # genetic position, cM (from pos_bp)
```

- **Create `genome_meta` with explicit typed columns (Option B — replaces the
  `dbWriteTable(overwrite = TRUE)` at `R/define_genome.R:124`):**

  ```r
  DBI::dbExecute(db_conn, "CREATE OR REPLACE TABLE genome_meta (
    locus_id INTEGER, locus_name VARCHAR, chr INTEGER, chr_name VARCHAR, pos_bp BIGINT)")
  # genome_meta_df has pos_bp as an R double (from round()); inserting into the
  # pre-typed BIGINT column widens it losslessly (verified: 3e9 stores exactly).
  # Use duckdb_register + INSERT (NOT dbWriteTable), matching the existing chr_meta
  # write at define_genome.R:159-164 — dbWriteTable advances R's RNG via a random
  # temp name and would shift the seeded draw sequence.
  duckdb::duckdb_register(db_conn, "__tmp_genome_meta", as.data.frame(genome_meta_df))
  DBI::dbExecute(db_conn,
    "INSERT INTO genome_meta SELECT locus_id, locus_name, chr, chr_name, pos_bp FROM __tmp_genome_meta")
  duckdb::duckdb_unregister(db_conn, "__tmp_genome_meta")
  ```

  Build `genome_meta_df` **without** `introduced_gen` (decision #8) — columns
  `locus_id, locus_name, chr, chr_name, pos_bp` only.
- **Create `genome_map` with explicit typed columns** (same approach; `pos_cM
  DOUBLE`, nullable `sex`/`line_name`/`map_name` VARCHAR), then write one row per
  locus: `id_genome_map` (via `next_int_id`), `locus_id`, `locus_name`,
  `sex = NULL`, `line_name = NULL`, `map_name = "default"`, `pos_cM = pos_cM_i`.
  Run `validate_genome_map()` after the write.
- Add the per-sex `chr_meta` `recombines_M`/`recombines_F` columns (decision #7)
  as part of genome setup.
- **Validate `pos_bp` *before* writing any table (R4).** Per chromosome, the
  rounded `pos_bp` must be all-positive, strictly increasing, and duplicate-free.
  Very short chromosomes or very dense marker counts make `round(seq(...))` collide
  or hit `0`. On failure, **error clearly** — e.g. *"chromosome {chr}: physical span
  {chr_len_bp} bp is too short to place {n_chr_loci} loci as positive,
  strictly-increasing integer positions; increase `chr_len_Mb[i]` or reduce loci."*
  **Do not silently repair** with `pmax()`/jitter (no such rule is adopted or
  tested). Minimum requirement for a simulated genome: `chr_len_bp_i` large enough
  that even spacing yields `n_chr_loci` distinct positive integer positions.
- Then assert `pos_cM` non-decreasing per chromosome (guaranteed here since
  `pos_bp` is strictly increasing and `cM_per_Mb > 0`).
- `chr_len_Mb` **stays** the physical-length argument; only `pos_cM` derives from
  the rate.

**Deliverable:** `define_genome()` writes `genome_meta.pos_bp` + a single-map
`genome_map`; schema descriptions + reserved cols updated; `genome_map`
registered in `pop$tables`; unit test on generated values.

## Step 2 — founder-LD helpers key on `genome_map` (`d_cM`); fix the `genome_meta` write

The mosaic and gaussian-copula LD models use adjacent-locus distance `d_Mb`
(`R/founder_haplotype_helpers.R:53-73`, `:123-152`):
`p_switch = 1 - exp(-switch_rate · d_Mb)`, `rho = exp(-decay_rate · d_Mb)`.

- Change `d_Mb → d_cM`, distances from the resolved default map via
  `resolve_genome_map(conn)` (LD is a population property → the shared default map).
  Traversal **order** unchanged.
- `switch_rate`/`decay_rate` meaning shifts "per Mb" → "**per cM**" (at the default
  1 cM/Mb the numbers coincide). Update the `@param` docs and the
  `requires ... pos_Mb` guards (`R/define_founder_haplotypes.R:312-314`,
  `:326-328`) to reference the map.

**Fix the `genome_meta` rewrite (critique #2 — type preservation).** The current
`founder_allele_freq` write is
`SELECT * FROM genome_meta → dbWriteTable(overwrite = TRUE)`
(`R/founder_haplotype_helpers.R:187-189`), which round-trips through an R data
frame and would **demote `pos_bp BIGINT`**. Replace with:

```r
# add column once (idempotent), then UPDATE in place — never rewrite the table
if (!"founder_allele_freq" %in% cols) dbExecute(conn,
  "ALTER TABLE genome_meta ADD COLUMN founder_allele_freq DOUBLE")
# UPDATE via a registered temp relation keyed on locus_id
```

The old code used `dbWriteTable` partly to match a historical RNG advance
(`R/founder_haplotype_helpers.R:207-208`); under pre-1.0.0 **we don't preserve that
old stream** — the only requirement is that the new path is deterministic and
self-consistent (same seed reproduces going forward). Free to change it; just keep
`pos_bp` typed `BIGINT`.

**Deliverable:** founder-LD helpers read `d_cM`; the `genome_meta` write no longer
rewrites the table (add a test asserting `pos_bp` DuckDB type stays `BIGINT` after
`define_founder_haplotypes()`); params documented per-cM; forward reproducibility
test (same seed → identical founder haplotypes on repeated runs).

### Map helpers (R3/R7/R12) — build these in Step 1, before any reader uses them

These are the single reliable source of positions for all downstream code, so their
contracts are strict (R7 — this is where most downstream bugs are prevented).

```r
resolve_genome_map(conn, sex = NULL, line_name = NULL, map_name = "default")
```

Returns **exactly one row per `genome_meta.locus_id`** — columns
`locus_id, locus_name, chr, chr_name, pos_bp, pos_cM` — **ordered by `locus_id`**.
Within the chosen `map_name`, applies the decision-#6 per-locus precedence
`(sex=S, line=L) → (sex=S, NULL) → (NULL, line=L) → (NULL, NULL)`. It **errors** if:

- any `locus_id` in `genome_meta` is missing from the resolved map;
- `pos_cM` is not monotonic non-decreasing within any chromosome **after** fallback
  resolution (validate the *resolved* map, not just raw stored rows — critique #4).

Founder LD calls it `(sex=NULL, line=NULL, map_name="default")`; recombination calls
it with the parent's sex/line. Single map source shared by founder LD, R
recombination, and (later) the optimized kernel.

```r
validate_genome_map(conn)   # run after every map write
```

Structural integrity of the stored rows (R7):

- **logical-key uniqueness** `(locus_id, sex, line_name, map_name)` with explicit
  NULL normalization (SQL `NULL ≠ NULL` won't catch duplicate
  `(locus, NULL, NULL, "default")` rows on its own);
- every `locus_id` and `locus_name` **agrees with `genome_meta`** — no orphan loci,
  and one `locus_id` is never paired with two different `locus_name`s;
- every row has a valid `sex` (`NULL`/`'M'`/`'F'`) and a **non-missing `map_name`**.

## Step 3 — recombination readers use `genome_map`

Reader swap so the package works after `pos_Mb` is gone. The readers now consume the
default map rows written in Step 1; they no longer reference `genome_meta.pos_Mb`
(which no longer exists):

- `build_chr_info()` (`R/recombination_helpers.R`): reads `pos_cM` from a resolved
  map; `chr_pos`/`chr_len` are cM. **The seam is fully wired (post-critique):**
  `add_offspring()` resolves + caches one descriptor set per unique
  `(parent_sex, parent_line)` among the parents and each gamete uses its
  **producing parent's** resolved map. The chromosome classification (autosome vs
  special, compact index remap) is position-independent and computed once. In v1
  only the default map exists so every gamete resolves to it (behavior-neutral,
  guarded by the forward-repro test); adding a sex/line-specific `genome_map` row
  now changes recombination for the matching gametes with no further code change
  (regression test: `test-genome_map.R` "add_offspring() uses the gamete-producing
  parent's resolved map").
- **Parent `line_name` source (R9):** `add_offspring()` fetches `ind_meta.line_name`
  alongside sex/ploidy and passes parent sex **and** parent line to the resolver.
  `L` = the **parent's** `ind_meta.line_name` (a gamete's map is its parent's) —
  **never** the offspring line or the allele's `line_origin`.
- `make_gamete()` (`:50-80`): `lambda = chr_len / 100` already yields Morgans, and
  `cM / 100 = Morgans` **by definition** — so the `/100` becomes *exactly* correct
  rather than an implicit 1-cM/Mb assumption. Crossovers/`findInterval` in cM.
- `add_offspring.R:217`: source positions from `genome_map`, not `pos_Mb`.
- Special-chromosome branch in `add_offspring` reads `recombines_M`/`recombines_F`
  (decision #7 is locked — no conditional) per the gamete's parent sex.
- **`R/chr_meta_helpers.R`** also reads `recombines`: `is_plain_autosome()`
  (`:99-102`) and `get_chr_meta_map()` (`:33-38`). `is_plain_autosome()` must read
  the **sex-specific** column for the gamete-producing parent, not a single
  `recombines`. This file is a required touchpoint for the split.

Full sex/line-aware resolution is only *exercised* once non-default map rows exist
(future / optimization plan); the **schema and read seam ship now**.

**Deliverable:** recombination reads `genome_map`; forward reproducibility test
(same seed → identical offspring haplotypes on repeated runs of the new code).

## Step 4 — `define_chr()`, tests, then docs **last**

- `define_chr()` (decision #7): accept `recombines_M`/`recombines_F`, with
  `recombines` as the primary shorthand setting both; default rows updated.
- Update tests referencing `pos_Mb`; add the position + map + resolution tests (see
  Testing).
- `DESCRIPTION` version bump + `NEWS.md` (breaking `genome_meta`/`chr_meta` change).
- **Touchpoint completeness gate.** The `rg` sweep is authoritative — the lists in
  this plan are the *known* set, but run and clear all three before merging:
  - `rg pos_Mb` — 7 R files (`define_genome`, `founder_haplotype_helpers`,
    `define_founder_haplotypes`, `recombination_helpers`, `add_offspring`,
    `schema`, `sql_utils`) **and** the tests
    (`test-define_genome.R`, `test-restore_pop.R`, `test-mutate_table.R`, plus the
    manual scripts `tests/test_new_functions.R`, `tests/test_dynamic_genome.R`).
  - `rg '\brecombines\b'` — includes `chr_meta_helpers.R`, `define_chr.R`,
    `recombination_helpers.R`, `add_offspring.R`, `schema.R`, `sql_utils.R`,
    `define_genome.R`, and tests (`test-define_chr.R`, `test-long-schema.R`,
    `test-sex_chromosomes.R`).
  - `rg introduced_gen` — per the deletion sweep above.
- **`CLAUDE.md` and user-facing docs are updated DEAD LAST (R2), only after code +
  tests are green**, and describe **only what actually landed** — never aspirational
  design. Update: `genome_meta` (`pos_Mb`→`pos_bp`, `introduced_gen` removed), the
  new `genome_map` table, `chr_meta` per-sex `recombines`. Teach the split
  explicitly (R13): **`genome_meta` = physical bp + chip/locus metadata; `genome_map`
  = genetic cM; join via `locus_id`** when both are needed — so the old conflation
  can't reappear under a new column name.

---

## Testing strategy

No cross-version parity tests (pre-1.0.0 — we don't compare against old `HEAD`).

1. **Forward reproducibility (key).** Fixed seed, run the *new* code twice →
   identical `genome_map`, founder haplotypes (both methods), and offspring
   haplotypes. Same-seed determinism of the current implementation, not a match to
   any prior version.
2. **Position generation.** `pos_bp` integer, strictly increasing per chromosome,
   1-based, no duplicate after rounding; `pos_cM` non-decreasing;
   `pos_cM = pos_bp/1e6 · rate` (derived straight from the sampled `pos_bp`).
3. **Map resolution.** With the default map + one `sex='M'` override row, a ♂ gamete
   resolves to the `M` row and a ♀ gamete to the default (decision #6 precedence);
   `resolve_genome_map()` returns exactly one row per locus, `locus_id`-ordered,
   monotonic `pos_cM` per chromosome, and **errors on a missing locus** (R7).
4. **Sex-specific recombination.** `recombines_M = FALSE` on all chromosomes ⇒
   ♂ gametes are pass-through (achiasmy); ♀ recombine normally.
5. **Non-default rate.** `cM_per_Mb = 2` scales genetic distances → both
   recombination and LD scale accordingly (sanity).
6. **Type preservation.** `pos_bp` DuckDB type stays `BIGINT` after
   `define_founder_haplotypes()` (guards critique #2's write fix).
7. **Old columns gone (R14).** Assert the exact `genome_meta` column set — no
   `pos_Mb`, no `pos_cM`, no `introduced_gen`; confirm `genome_map` and per-sex
   `chr_meta.recombines_M/_F` exist and are registered in schema metadata / `SYSTEM_TABLES`.
8. **`pos_bp` failure mode (R4).** A chromosome too short for its loci raises the
   clear error, and no tables are written (validation happens before any write).
9. **Map-helper strictness (R7).** `validate_genome_map()` catches a duplicate
   default-key row, an orphan `locus_id`, and a `locus_id`↔`locus_name` mismatch;
   `resolve_genome_map()` errors on a locus missing from the resolved map.
10. **Haplodiploidy smoke test.** all-chr `copy_mode_M = "half"` → functional
    haploid males.
11. **Existing suite** green; `R CMD check`.

---

## Explicitly out of scope (engine or later, **no schema change** when built)

- Non-linear / hotspot map **generation**, PAR maps, sex/line-specific map
  **authoring** function (`define_genetic_map()`), map versions — schema is ready.
- Crossover interference model (Poisson → gamma) — engine.
- Polyploid *meiosis* mechanics — engine + `ind_meta.ploidy`.
- PLINK/VCF export — `pos_bp` merely makes it lossless.
- Morgans-native kernel / C++ / `dqrng` / batching — `plans/optimize_add_offspring.md`.
- `UINTEGER` storage optimization for `pos_bp`.

---

## Implementation order

| Step | What | Breaking? | Status |
|---|---|---|---|
| 1 | `define_genome()` builds `genome_meta` via explicit `CREATE OR REPLACE TABLE (... pos_bp BIGINT)` + register/INSERT (Option B), and the default `genome_map` (explicit typed CREATE; `pos_cM = pos_bp/1e6·rate`; `locus_id`, `map_name="default"`); `cM_per_Mb` validated; `resolve_genome_map()` + `validate_genome_map()`; schema/reserved/system-table registration; remove `introduced_gen` (touchpoint audit) | schema | ✅ Done |
| 2 | founder-LD helpers key on `d_cM` via `resolve_genome_map()`; **fix `genome_meta` write (ALTER+UPDATE, preserves `BIGINT`)** + type test | schema | ✅ Done |
| 3 | recombination readers (`build_chr_info`, `make_gamete`, `add_offspring`, `chr_meta_helpers`) read resolved map; fetch parent sex + `line_name` | reader swap to resolved map | ✅ Done |
| 4 | `chr_meta` per-sex `recombines_M`/`_F` + `define_chr()` shorthand; docs (`CLAUDE.md`), tests + forward-reproducibility test, version bump | schema (chr_meta) | ✅ Done |

Steps 1–4 land together as one migration (no broken intermediate state). Only
**after** this is merged and green does `plans/optimize_add_offspring.md` begin.

---

## Resolved: `pos_bp` BIGINT mechanism (Option B)

**Decision (chosen by the user): Option B — explicit typed `CREATE OR REPLACE
TABLE`.** `define_genome()` creates `genome_meta` (and `genome_map`) with explicit
column types, then appends rows — replacing the current
`dbWriteTable(overwrite = TRUE)` (`R/define_genome.R:124`), which infers types from R
and would produce `DOUBLE` (from `round()`) or, with `as.integer()`, a 32-bit
`INTEGER` that overflows above ~2.1 Gb. No new dependency; explicit control of every
column type.

*Rejected:* Option A (`bit64::integer64`) works but adds a dependency and
`integer64` coercion caveats; Option C (store `pos_bp` as `DOUBLE`) abandons
decision #2's clean-integer / lossless-export rationale.

### The type is column-level and "sticky" — verified against DuckDB

Once `pos_bp` is declared `BIGINT`, the type belongs to the **column**, not to each
inserted value. Every write path DuckDB uses casts the incoming value to `BIGINT`
(all rows below were run and confirmed):

| Operation | R source type | `pos_bp` after |
|---|---|---|
| `CREATE ... pos_bp BIGINT` then `dbWriteTable(append = TRUE)` | `double` (from `round()`) | **BIGINT** — `3e9` stored exactly |
| `duckdb_register` + `INSERT ... SELECT` | `double` | **BIGINT** |
| user `INSERT`/append of a plain `integer` | `integer` | **BIGINT** (auto-widened) |
| `define_chip()` `ALTER ADD COLUMN` / `mutate_table()` `UPDATE` | — | **BIGINT** (table not recreated) |
| `restore_pop()` (reopens the DuckDB *file*) | — | **BIGINT** (types persist in the file) |
| `dbWriteTable(overwrite = TRUE)` ⚠️ | any | **demoted** (recreates the table) |

**Answer to "what if a user adds rows with an integer `pos_bp`?"** DuckDB widens it
to `BIGINT` automatically on insert — **no R-side conversion or special exception is
needed.** Even a value beyond the 32-bit `INTEGER` ceiling stores correctly, because
the destination column is 64-bit.

**The only way to demote the type is to *recreate* the table** (`overwrite = TRUE` /
`CREATE TABLE AS` from R). After this migration the package never does that to
`genome_meta`: the two existing `overwrite` writers are both being removed (Step 1 →
explicit `CREATE`; Step 2 → `ALTER + UPDATE`), and every "change anything" tool
(`define_chip`, `mutate_table`, user appends, `restore_pop`) uses
`ALTER`/`UPDATE`/`INSERT`/file-reopen — all type-preserving. A user could still
demote it by deliberately collecting the whole table and writing it back with
`overwrite = TRUE`, but that is a destructive full-table replace — true of every
typed column in every table, not specific to `pos_bp`.

### Scope note

Option B *enables* adding real DB-level constraints later (`UNIQUE(locus_name)` /
`PRIMARY KEY(locus_id)`), since the table is no longer recreated on every write — but
this migration **keeps** the existing R-side `locus_name` uniqueness validation
(`R/define_genome.R:85-93`) and does **not** add DB constraints. Separate, optional
follow-up.

All reviewed items are resolved; there are no open blockers.

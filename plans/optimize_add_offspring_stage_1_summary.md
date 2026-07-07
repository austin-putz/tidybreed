# Optimize the haplotype write path — Stage 1 summary (direct long write + batching)

**Date:** 2026-07-07
**Branch:** `feat/position-coordinates-genome-map`
**Scope executed:** Stage 1 of
[`plans/optimize_add_offspring.md`](optimize_add_offspring.md) — convert **all three**
haplotype write paths from the wide-matrix + `UNPIVOT` write to **direct long, batched,
one-transaction-per-call** inserts, and create/register the `ind_crossover` schema.
Allele values still come from the base-R `make_gamete()` seam (Stage 0), so Stage 1 is
**output-neutral** for `ind_haplotype` / `founder_haplotypes`. Stages 2–4 (dqrng R kernel,
C++ kernel, parallel) were **not** started.

**Status: complete and green.** Full suite **1229 pass / 0 fail / 0 error** (9 pre-existing
warnings, unrelated). No `DESCRIPTION`/`NEWS` bump and no commit yet (reserved for when
requested). No new package dependency — memory detection is hand-rolled.

---

## What changed

### New shared helpers — `R/haplotype_write_helpers.R`
- **`.write_long_haplotypes(conn, long_frame)`** — the single write core for
  `add_offspring()` and `add_founders()`. Registers a long frame
  `(id_ind, parent_origin, strand, line_origin, locus_id, allele)` and does one
  `INSERT ... SELECT` joining `genome_meta` for `locus_name`. **No `UNPIVOT`.** Uses
  `duckdb_register` + `INSERT` (RNG-neutral, zero-copy) — Appender API rejected (marginal
  gain, C++ boundary; documented in the helper).
- **`detect_available_memory()`** — hand-rolled cross-OS available-RAM probe
  (`/proc/meminfo` MemAvailable on Linux; `sysctl hw.pagesize` + `vm_stat` on macOS
  Intel+ARM; `wmic`/PowerShell on Windows). **Never errors**; returns `NA` on failure so
  callers fall back to a conservative fixed budget.
- **`resolve_batch_size()`** — precedence: explicit `batch_size` > `max_batch_mem` >
  RAM-aware auto-pick (~25% of available) > conservative 512 MB fallback. Clamped to
  `[1, n_total]`. `bytes_per_offspring_row = 200`, calibrated from Stage-1 measurements
  (each offspring emits `2 x n_loci` long rows at ~90 B/row peak incl. transient copy).
- **`parse_mem_size()`** — `512e6` / `"512MB"` / `"2GB"` → bytes.

### `add_offspring()` — `R/add_offspring.R`
New args `batch_size = NULL, max_batch_mem = NULL`. Offspring are generated + written in
batches of `B` in `matings` row order inside **one transaction**. Each batch's wide seam
matrices are reshaped to long parallel vectors in R (bounded by `B`, no dense
`n_offspring` frame, no `UNPIVOT`) and written via `.write_long_haplotypes()`. **RNG
discipline:** all gamete draws happen in the loop; the single `dbWriteTable(ind_meta)`
(which advances R's RNG) is issued **after** the loop, matching the pre-Stage-1 stream
position exactly — preserving this call's *and* downstream seeded output. `matings`
normalization stays at the top (non-foreclosing for a future count-based spec).

### `add_founders()` — `R/add_founders.R`
New args `batch_size = NULL, max_batch_mem = NULL`. The single `sample()` of all haplotype
indices stays up front (the **founder RNG trap** — batch the write, never the draw); the
selected `(2 x n_founders) x n_loci` wide matrix and per-chr `UNPIVOT` are gone. Rows are
emitted directly in long form from the pre-drawn indices against the pool, honoring
`copy_mode_M/F` / `hemi_parent`, in one transaction.

### `define_founder_haplotypes()` pool write — `R/founder_haplotype_helpers.R`
The full `n_hap x n_loci` `fh` frame is no longer materialized. The pool is streamed to
`founder_haplotypes` in haplotype batches inside one transaction. **RNG-preserving trick:**
`dbWriteTable`'s RNG advance is fixed and size-independent (verified empirically), so the
**first** batch uses `dbWriteTable` (creates the table + the one historical advance) and
the rest use `register`+`INSERT` (neutral) — exactly one advance per call, downstream
stream unchanged.

### `ind_crossover` schema (created here; rows in Stage 2)
Empty table created in `define_genome()` and registered in `R/sql_utils.R`
(`TABLE_RESERVED_COLS`, `TABLE_PRIMARY_KEYS`, `TABLE_ROW_KEYS`, `SYSTEM_TABLES`,
`IND_TABLE_ID_IND_COLS`), described in `R/schema.R`, and added to `pop$tables`. Survives
`restore_pop()`. No `store_crossovers` argument yet (the base-R path emits no events; the
arg + writes land with the Stage-2 kernel).

---

## Verification

- **Output-neutrality (byte-identical vs pre-Stage-1 golden):** founder pool, founder
  `ind_haplotype`, and full-pipeline offspring `ind_haplotype` all IDENTICAL, for both the
  autosome design and a `define_chr()` sex-chromosome (half/full) design — including a
  **forced-`B=3` multi-batch run of all three paths** (pool, founders, offspring).
- **Batch invariance:** `batch_size ∈ {1, 3, 5, 7, all, default}` → identical
  `ind_haplotype` (autosome + sex chr). Locked in `tests/testthat/test-haplotype_write_helpers.R`.
- **Transaction atomicity:** a mid-write failure injected into batch 2 rolls back the whole
  generation — `ind_meta` and `ind_haplotype` unchanged. (testthat + manual.)
- **Helpers:** `parse_mem_size`, `detect_available_memory` (positive/`NA`, never errors),
  `resolve_batch_size` (precedence + clamp) unit-tested.
- **Full suite:** `testthat::test_local()` → **1229 pass / 0 fail / 0 error**.

### Measured wins (5000 loci, 10 chr)

| metric | pre-Stage-0 | Stage-0 | **Stage-1** |
|---|---|---|---|
| `add_offspring(200)` mem_alloc | 43.2 GB | 1.55 GB | **464 MB** |
| `add_offspring(2000)` peak, default | — | — | 1867 MB |
| `add_offspring(2000)` peak, `batch_size=100` | — | — | **329 MB** |
| `add_offspring(8000)` peak, `batch_size=100` | — | — | **183 MB** (4x offspring, peak bounded) |

Peak memory is now **decoupled from total offspring** — a single high-fecundity mating
streams to disk batch by batch. The RAM-aware default sizes `B` so peak ≈ 25% of available
memory (e.g. `B ≈ 107` at 50k loci for a 1 GB budget → ~950 MB peak), protecting 8–16 GB
machines automatically.

---

## Files touched

**Package code:**
- `R/haplotype_write_helpers.R` — **new**: helpers above.
- `R/add_offspring.R`, `R/add_founders.R` — batched long write + args + one transaction.
- `R/founder_haplotype_helpers.R` — batched pool write, RNG-preserving.
- `R/define_genome.R`, `R/sql_utils.R`, `R/schema.R` — `ind_crossover` create + registration.
- `tests/testthat/test-haplotype_write_helpers.R` — **new**: helpers, batch invariance,
  `ind_crossover` schema, atomicity.
- `man/*` + `NAMESPACE` — regenerated (also removed orphaned `.Rd`/exports for
  already-deleted `initialize_genome`/`migrate_schema_meta`/`set_qtl_effects_multi`).

**No new dependency** (`duckdb`/`DBI`/`stats` already imported).

---

## Open items

1. **Commit Stage 1** with `DESCRIPTION` version + `NEWS.md` entry when ready (deferred).
2. **Optional:** refresh the line-range buckets in
   `dev/benchmarks/profile_haplotype_paths.R` (its per-step labels reference pre-Stage-1
   line numbers; a focused Stage-1 bench was used instead).
3. **Regroup gate:** confirm the memory-ceiling + write-time win before starting Stage 2
   (dqrng R reference kernel — the first *intentional* seeded-output change).

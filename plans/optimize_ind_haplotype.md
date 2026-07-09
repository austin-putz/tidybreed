# Optimize the `ind_haplotype` write path (PRIMARY-KEY + R-side assembly)

**Status**: Problem statement + solution options for review. **No code changed.**
Profiled 2026-07-08 on macOS (Opus session). This is the follow-up to
`plans/optimize_add_offspring.md` Stages 0–3, which already put `add_founders()`
and `add_offspring()` on the shared batched long-write path
(`resolve_batch_size()` loop → `.write_long_haplotypes()`). This doc isolates the
**remaining** bottleneck that profiling surfaced: the cost of the
`ind_haplotype` **PRIMARY KEY** during bulk inserts, plus the R-side frame
assembly.

**Scope note**: `.write_long_haplotypes()` and the `ind_haplotype` table are
**shared** by `add_founders()` *and* `add_offspring()`. Every fix here speeds up
both — the founder build and the per-generation matings in the date loop.

---

## 1. Symptom

The swine vignette (`vignettes/swine/…`) is "super super slow". `add_founders()`
was suspected of being un-optimized. It is **not** — it is on the optimized bulk
path — but at swine scale it still takes ~2 minutes for one founder call.

Swine scale: `n_loci = 10000`, `n_chr = 18`, 1000-haplotype line-A pool,
**2000 founders** (400M + 1600F) → **40,000,000 `ind_haplotype` rows**
(2 homologs × 10000 loci × 2000 individuals).

---

## 2. Profiling evidence

### 2a. `add_founders()` end-to-end (`Rprof`, swine scale)

```
add_founders wall-clock:      116.70 s   (40M rows written; 342,753 rows/s)

by self-time:
  rapi_execute (DuckDB exec)   62.25 s  62.4%   <- INSERT in .write_long_haplotypes()
  rbind                        15.99 s  16.0%   <- do.call(rbind, parts), ~16 GB churned
  mapply                       10.24 s  10.3%   <- inside frame build / register coercion
  data.frame                    3.26 s   3.3%
  gc / paste / match / …        ~few s
by total-time:
  .write_long_haplotypes       57.18 s  57%     (DuckDB insert)
  do.call(rbind, …)            29.01 s  29%     (R-side long-frame assembly)
```

**≈53% DuckDB insert, ≈25% R-side `rbind`/frame assembly.**

### 2b. Attribution of the DuckDB insert (isolated micro-bench, 10M rows)

| Variant | Time |
|---|---|
| **A. Current** — INSERT…SELECT **JOIN** genome_meta, **PK on** | 11.33 s |
| **B. No join** (locus_name precomputed in R), **PK on** | 14.72 s |
| **C. No join + NO primary key** | **3.03 s** |
| D. `dbAppendTable`, no PK | 3.20 s |
| R-side: build `locus_name` column in R (`gm$locus_name[locus_id]`) | 0.35 s |

**Conclusion**: the **4-column PRIMARY KEY `(id_ind, parent_origin, strand,
locus_id)`** is the dominant insert cost. Dropping it: **~11–15 s → ~3 s
(~4–5×)**. The `genome_meta` **join is a non-issue** (A vs B is within noise;
removing it did not help). Extrapolated to the real 40M-row write, the PK is
costing on the order of **~45 s**.

### 2b′. Integrated bench — PK vs no-index vs `id_ind` index (settles §6.1)

Run 2026-07-08 (macOS, Opus session), `scratchpad/bench_pk_vs_index.R` +
`bench_delete_confirm.R`. 10M rows (n_loci=5000, n_ind=1000), inserted in **20
incremental batches through the production write path** (`duckdb_register` +
`INSERT…SELECT JOIN genome_meta`), then a `remove_rows()`-style
`DELETE … WHERE id_ind IN (…)` of 10 individuals (100k rows), 3 reps on
early/late/scattered id sets:

| Config | Insert (10M) | vs PK | DELETE 10 ind (100k rows) |
|---|---|---|---|
| **A. current 4-col PK** | 25.0 s | 1.00× | ~0.74 s |
| **B. no index (1a)** | **12.0 s** | **2.08×** | **~0.010 s** |
| **C. `id_ind` index (1b)** | 16.6 s | 1.50× | ~0.155 s |

**Two decisive findings:**

1. **1a dominates 1b on *both* axes.** No index is faster to insert (2.08× vs
   1.50×) *and* faster to delete than the `id_ind` index. The index is strictly
   dominated — it buys nothing.
2. **The plan's assumption that an `id_ind` index would speed up `remove_rows()`
   was wrong.** DuckDB deletes via a vectorized scan with min/max **zonemap**
   pruning (rows land clustered by insertion order, which mirrors sequential
   founder numbering), so a secondary index doesn't accelerate the scan — it only
   adds maintenance cost. No index (1a) deletes ~10 ind in **~10 ms**; the
   `id_ind` index is ~15× slower at ~155 ms; the 4-col PK is ~74× slower at
   ~740 ms. So dropping the PK actually makes `remove_rows()` *faster*, not
   slower — the opposite of the concern raised under Lever 1.

*Calibration note*: the insert speedup here is **~2×**, not the ~4–5× of the
isolated §2b bench, because the per-batch `register` + `JOIN` is a fixed cost
that dilutes the PK saving at 20 small batches. Production uses RAM-driven batch
sizing (`resolve_batch_size()`); larger/fewer batches push the realized speedup
back up toward §2b's ~4×. Treat §4's projections as the optimistic end (few big
batches) and ~2× as the conservative end (many small batches).

### 2c. Why the PK is safe to reconsider

- `ind_haplotype` is created with the PK at `R/define_genome.R:188–193`.
- **The PK is redundant, not load-bearing.** Uniqueness of `ind_haplotype` rows
  is already guaranteed *upstream*, independent of this table's own constraint:
  - `id_ind` is unique across individuals because it is assigned as
    `MAX(id)+1` **from `ind_meta`** (`R/add_founders.R:196–200`), and `ind_meta`
    itself carries a PK on `id_ind` (`R/open_pop.R:265`). New founder/offspring
    IDs therefore cannot collide with existing ones.
  - Within one individual, the `(parent_origin, strand, locus_id)` tuples are
    emitted deterministically, exactly once each, by the write path.
  - So the 4-column tuple is unique by construction; the DB PK only *re-checks*
    what R already guarantees.
- **Correction to an earlier draft**: it is **not** true that "nothing deletes
  `ind_haplotype`." `remove_rows()` issues `DELETE FROM ind_haplotype WHERE
  id_ind IN (...)` — both in explicit-table mode and in `tables = "all"` mode,
  which sweeps every `ind_*` table (`R/remove_rows.R:248–253`). This does **not**
  block dropping the PK (a `DELETE` is correct without one), but it (a) means the
  safety argument must rest on the upstream guarantee above, not on "append-only",
  and (b) has a performance consequence — see Lever 1's *Cost/risk*.
- No **DB-level foreign key** references `ind_haplotype` (the schema's "FK to …"
  columns are logical/R-enforced, not `FOREIGN KEY` constraints). Confirm with a
  `grep -rn "FOREIGN KEY" R/` before pulling the PK — a DuckDB FK target requires
  a PK/unique, and that would be the one thing that hard-blocks 1a.
- DuckDB is columnar: reads/filters/joins on `ind_haplotype` do **not** require
  the PK for performance (scans lean on min/max zonemaps, not the ART index).
- Contrast `ind_genotype`: it **does** use `INSERT OR REPLACE` (idempotent dosage
  cache) and genuinely needs its PK — **do not touch that one.**

---

## 3. The levers

### Lever 1 — the `ind_haplotype` PRIMARY KEY  (~53% of runtime; the big win)

The PK's ART-index maintenance is what makes the insert super-linear as the table
grows. Options:

- **1a. Drop the PK entirely.** Remove `PRIMARY KEY (…)` from the
  `CREATE TABLE ind_haplotype` DDL in `R/define_genome.R`. Insert ~4–5× faster;
  benefits `add_offspring()` too. Uniqueness remains guaranteed upstream (§2c).

  *Edit surface (smaller than an earlier draft claimed)*: only **two** places.
  - The `CREATE TABLE ind_haplotype` DDL at `R/define_genome.R:188–193`.
  - The human-readable PK mention in the schema doc string at `R/schema.R:102`.

  *Not affected — do not list these as work*:
  - `restore_pop()` recreates **no** DDL. It just `dbConnect`s to the existing
    file (`R/restore_pop.R:65–66`); the table keeps whatever shape it was
    created with. Nothing to change.
  - `archive_replicate()` derives archive columns **generically** from
    `information_schema` and never replicates a PK (`R/archive_replicate.R:363–371`).
    Archive copies of `ind_haplotype` already have no PK today — which is also
    live proof the table reads/writes fine without one.

  *Cost/risk (genuine, accept explicitly)*:
  - **Loses the loud-failure backstop.** Today a double-run of `add_founders()`
    or a buggy re-insert errors immediately on PK violation. Drop it and the same
    mistake silently duplicates up to 40M rows, surfacing later as doubled TBVs.
    In this "break freely" pre-1.0 phase — where cells get re-run constantly —
    this is the backstop most likely to actually be missed. Partial mitigation:
    a debug-mode / test-only uniqueness assertion on the built population (catches
    it in the suite, **not** in an interactive session).
  - **`remove_rows()` actually gets *faster*, not slower** (measured — §2b′).
    The earlier worry that an unindexed `DELETE … WHERE id_ind IN (...)`
    (`R/remove_rows.R:248–253`) would be slow proved false: DuckDB deletes via a
    zonemap-pruned vectorized scan, so no-index delete of 10 individuals runs in
    ~10 ms vs ~740 ms under the current PK (~74× faster). This is a *benefit* of
    1a, not a cost.
- **1b. Demote to a non-unique index on `id_ind`. ~~Candidate~~ REJECTED by the
  §2b′ bench.** The idea was that a single-column `id_ind` index would recover
  most of the insert win while keeping `remove_rows()` fast (every read/delete
  path filters by `id_ind`: `add_tbv` `R/add_tbv.R:210`, `add_dosage`
  `R/add_dosage.R:128`, `extract_genotypes` `R/extract_genotypes.R:250`,
  `remove_rows`). The measurement killed it: the index is **slower to insert than
  1a** (1.50× vs 2.08×) **and slower to delete than 1a** (~155 ms vs ~10 ms),
  because DuckDB scans lean on zonemaps, not secondary indexes — the index buys
  no read/delete speedup and only adds maintenance cost. **Strictly dominated by
  1a.** Keep this bullet only as the record of why "just add an index" is the
  wrong instinct here.
- **1c. Drop+recreate the PK around each bulk load.** `ALTER TABLE … DROP/ADD
  CONSTRAINT` before/after the batch loop. **Rejected** for a persistent,
  incrementally-appended table: rebuilding the index over the *entire* growing
  table on every `add_founders`/`add_offspring` call is worse than 1a and can be
  worse than the status quo at large table sizes.

**Recommendation: 1a (drop the PK, no replacement index) — settled by the §2b′
bench.** 1a is the fastest on insert *and* on `remove_rows()` DELETE; 1b is
strictly dominated; 1c is rejected on principle. The only thing 1a gives up is
the loud-failure backstop (a correctness-ergonomics cost, not a performance one),
mitigated by the test-only uniqueness assertion in §7.

### Lever 2 — R-side long-frame assembly  (~25%; ~29 s)

`add_founders()` builds a `parts` list (one `data.frame` per chr ×
parent_origin) and calls `do.call(rbind, parts)` per batch — O(n) copies, ~16 GB
churn, and the `mapply`/`data.frame`/`Make.row.names` overhead. We know the exact
total row count for a batch up front, so:

- Preallocate typed vectors (`id_ind`, `parent_origin`, `strand`, `line_origin`,
  `locus_id`, `locus_name`, `allele`) at full batch length and fill by slice, then
  build one `data.frame`/register once — **no `rbind`**. Or use a single
  vectorized `rep()`-based construction across all chromosomes at once.
- **Byte-identical output, no schema change.** Expected to reclaim ~20 s and
  likely dissolve most of the 10 s `mapply`.
- Applies identically to `add_offspring()`'s long-frame build.

### Lever 3 — fold `locus_name` into R, drop the join  (small, cleanup)

Precompute `locus_name` in R (`gm$locus_name[locus_id]`, 0.35 s) so
`.write_long_haplotypes()` can `INSERT … SELECT *` (or `dbAppendTable`) with **no
`genome_meta` join**. Negligible speed win on its own (the join isn't the
bottleneck), but it simplifies the write and pairs naturally with Lever 2. Keep
only if it doesn't complicate the helper's contract with `add_offspring()`.

---

## 4. Expected outcome

Rough, at swine scale (40M rows):

| | Insert | R build | Total |
|---|---|---|---|
| Current | ~57 s | ~29 s | **~117 s** |
| Lever 1 only (drop PK) | ~12–15 s | ~29 s | ~45–50 s |
| Lever 1 + 2 | ~12–15 s | ~5–8 s | **~25–35 s** |

⇒ a plausible **3–5× overall** speedup for `add_founders()`, with `add_offspring()`
inheriting the same wins.

---

## 5. Cross-cutting note (not this function's job)

In the *full* swine script the single largest waste is upstream of the write path:
`define_founder_haplotypes()` builds **six** pools (lines A–F, ~60M rows) while
only line A seeds founders (~50M rows never used as founders). That is a **script**
choice (the other lines are there to demo generation methods / future
crossbreeding), not an `add_founders()` defect. Flag to the user; out of scope for
the function fix.

---

## 6. Open decisions (settle before implementing)

1. **PK handling — SETTLED: 1a (drop, no replacement index).** The §2b′
   head-to-head bench decided it: 1a is fastest on both insert (2.08× vs the
   PK; 1b only 1.50×) and DELETE (~10 ms vs 1b's ~155 ms vs the PK's ~740 ms);
   1b is strictly dominated and rejected. The only remaining decision is whether
   to accept the lost loud-failure backstop (yes — mitigated by the §7 test-only
   uniqueness assertion). **No `restore_pop()`/`archive_replicate()` DDL to
   update** — neither recreates the table (see Lever 1 *Edit surface*); the only
   edits are the `define_genome.R` DDL and the `schema.R` doc string.
2. **Founder sampling option** (raised separately by the user; design-only here) —
   add `replace = TRUE/FALSE` (or `allow_self`/similar) to `add_founders()`.
   Semantics to pin down, because pure "without replacement" is impossible when
   `2 × n_founders > n_haplotypes` (swine: 4000 draws from a 1000 pool):
   - **(a) no self-homozygosity**: a single founder cannot draw the same pool
     haplotype for both of its copies (still with-replacement across founders).
   - **(b) globally unique**: each pool haplotype used by at most one founder —
     only viable when the pool ≥ 2× founders; otherwise error.
   Likely default stays `(with replacement)` = current behavior; option adds (a)
   and/or (b). Note `add_founders()` currently uses base-R `sample()` (no `seed`
   arg), unlike `add_offspring()`'s `dqrng` streams — decide whether to align.

---

## 7. Verification plan (for when implemented)

- **Reproducibility**: same `set.seed()` (base-R sampling) → identical
  `ind_haplotype` contents before vs after each lever (Levers 2 & 3 are
  output-neutral; Lever 1 changes only the table's constraint, not row contents).
  **Hard requirement — sort before comparing.** With no PK, on-disk row order
  follows insertion order and may differ from the PK-ordered layout, so the
  comparison MUST `ORDER BY id_ind, parent_origin, strand, locus_id` (or sort the
  collected frame) on both sides. A naive row-order-sensitive diff will fail
  spuriously even though the data is identical.
- **Correctness backstop for Lever 1**: a test that asserts
  `(id_ind, parent_origin, strand, locus_id)` uniqueness on a built population
  (R-side check replacing the dropped DB constraint). Note this only guards the
  test suite — an interactive double-insert would no longer error at write time
  (see Lever 1 *Cost/risk*).
- **FK precheck for Lever 1**: `grep -rn "FOREIGN KEY" R/` to confirm nothing
  declares a DB-level FK targeting `ind_haplotype` (such a constraint requires a
  PK/unique target and would hard-block 1a).
- **Perf**: re-run `scratchpad/prof_add_founders.R` (swine scale) — expect
  ~117 s → ~25–35 s; confirm `add_offspring()` in a small date loop also drops.
- **Regression**: full `testthat` suite green on both kernels; `R CMD check` 0
  errors; parity golden unchanged (write-path changes are output-neutral).

---

## 8. Explicitly out of scope

- Changing `ind_genotype`'s PK (it needs it — `INSERT OR REPLACE` cache).
- The genome-effects / TBV / dosage steps (separate hotspots; the 7× `add_tbv`
  and `add_genotypes(9k)` in the vignette are their own investigation).
- Any change to the C++ recombination kernel (Stage 3) — it is not on the
  founder path and not implicated here.
- Parallelism (`RcppParallel`) — that is `optimize_add_offspring.md` Stage 4.
```

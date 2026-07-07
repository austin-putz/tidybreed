# Optimize `add_offspring()` (and `add_founders()` write path)

**Status**: Draft for review — revised after critique round 2, then revised again
(2026-07-06) to make **`add_founders()` and `define_founder_haplotypes()`
first-class** in this optimization (they, not just `add_offspring()`, are the
slow founder path) and to **refresh all stale `R/*.R` line references** against
the current working tree (`add_offspring.R` is now 612 lines post-coordinate
migration; the old refs were off by ~120 lines).
**Depends on**: `plans/position_coordinates.md` (the `pos_Mb` → `pos_bp` +
`genome_map` schema migration). **That prerequisite is now IMPLEMENTED and green
in the working tree** (`feat/position-coordinates-genome-map`): `genome_meta.pos_bp`
(`BIGINT`), the long `genome_map` table, `resolve_genome_map()`/`validate_genome_map()`,
and the per-`(sex,line)` gamete-map cache in `add_offspring()` all exist. This plan's
Morgans-native kernel reads that resolved map.
**History**: this is the sole plan of record for offspring optimization. An earlier
exploratory doc (`plans/speed_up_offspring.md`) surveyed the option space (Option A
vs B, wide+UNPIVOT vs long, Rcpp packaging) and has been **deleted** — its decisions
are locked below and the one piece of standalone rationale worth keeping (the
CUDA/GPU rejection) is inlined in "Explicitly out of scope".

---

## Critique round 2 (2026-07-06) — resolutions

The second Codex review (`plans/optimize_add_offspring_critique.md`) raised 16
points. **All 16 are accepted; none required code changes** (Stage 0+ has not started —
these are all plan refinements for future work). The plan below is revised accordingly.
Point-by-point, with where each is now addressed:

| # | Critique | Resolution | Where |
|---|---|---|---|
| 1 | Kernel signature can't represent per-parent/sex/line maps | **Group-by-map calls** — one kernel call per resolved `(sex, line, map_name)` group; preserves the map seam | "Multiple resolved maps" (Stage 0) |
| 2 | Long-output kernel arrives before long-write/batching | **Stages reordered** — long write + batching becomes Stage 1, before any long-output kernel | "Stage sequence"; Stages 1–3 |
| 3 | Batching must move earlier | Folded into Stage 1 (prerequisite for the dqrng kernel) | Stage 1 (1b) |
| 4 | RNG/seed API unspecified | **`add_offspring(seed = NULL)`**; `NULL` draws one base seed from base-R RNG (upstream `set.seed()` still reproduces) | "Seed API" |
| 5 | dqrng stream scheme open but foundational | Key tuple **pinned**; **pre-Stage-3 uniform-parity spike** gates the C++ kernel | "dqrng stream-derivation"; Spike row |
| 6 | Kernel input layout underspecified for C++ | **Packed parent array** (homolog-major), `parent_*_idx` = packed-list indices | "Kernel input memory layout" |
| 7 | Special-chr dqrng/crossover path underspecified | dqrng `"special"` sub-stream for row choice + recombining special chr; **wrapper-level** R↔C++ parity | "Special-chromosome RNG under dqrng" |
| 8 | `add_founders()` batching RNG trap | **Pre-draw all indices once**, batch only materialization/write | Stage 1 (founder RNG trap) |
| 9 | `ind_crossover` should follow conventions | `next_int_id()`, created by `define_genome()`, registered in all schema tables | `ind_crossover` section; Stage 1 (1c) |
| 10 | Long writes should be transaction-scoped | **One transaction per `add_offspring()` call** | Stage 1 (1d) |
| 11 | Poisson sampler numerical/param bounds | **Log-accumulation** (not raw product); documented max-`lambda` + guard; stress test | Poisson sampler note; test 5 |
| 12 | Crossover-consistency test needs controlled parents | Fixture with homologs distinct at **every** marker (or homolog trace) | Testing strategy (test 6c) |
| 13 | Stale status / line references | Prereq marked **implemented**; **line refs genuinely refreshed in the round-3 revision** (round-2 claimed this but the refs were still ~120 lines stale after the coordinate migration grew `add_offspring.R` to 612 lines — now corrected against the working tree) | throughout; hot-spots; Phase 0 tables; Stage 0/1 |
| 14 | Benchmark deps / deliverables | `profvis`/`bench` dev-only, `requireNamespace()`-guarded + skip; new per-step instrumentation | Phase 0 Tools |
| 15 | Compiler-fallback wording optimistic | Source install **requires a compiler** post-Stage-3; R path is runtime/debug fallback + oracle, **not** a no-compiler install | Stage 3 (compiler requirement) |
| 16 | Output ordering should be canonical | Canonical order **pinned**; execution/parent-grouping must not change it or `o` | "Canonical output order" |

---

## Decisions locked in this round (constraints, not open questions)

These were decided in discussion and are **not** to be re-litigated during
implementation:

1. **No wide matrices — anywhere in the write path.** The generator does
   **not** build an `n_offspring × n_loci` allele matrix and then `UNPIVOT` it
   into long form. It samples a gamete from the sire and a gamete from the dam
   per offspring and emits **long-format rows directly** (one row per offspring
   × parent_origin × locus). Wide format is a legacy artifact of the current
   code being deleted — it costs ~1 GB per matrix at target scale plus a
   `UNPIVOT` pass, to build something we immediately un-build. **The only
   legitimate wide representation in the package is the user-facing
   `extract_genotypes()`** (`R/extract_genotypes.R` — one row per individual,
   one column per locus, for `colMeans()` allele frequencies / PCA / export),
   computed on demand from the long store. Wide never appears in
   `add_offspring()` / `add_founders()`.

2. **RNG model = per-gamete deterministic sub-streams via `dqrng`.** *(Changed
   from the earlier "match base R's global stream" idea.)* Each gamete —
   identified by `(offspring_index, parent_role ∈ {1,2})` — draws from its own
   independent, deterministically-seeded counter-based RNG sub-stream derived
   from one base seed. Consequences:
   - **Reproducible regardless of core count or scheduling** — gamete *i* always
     uses stream *i*, so single-threaded and any-threaded runs are identical.
   - **Enables across-individual parallelism with no rewrite** (see the
     Parallelization section).
   - Uses a counter-based generator (`dqrng`: threefry/pcg) that has
     **identical implementations in R and C++**, so R↔C++ parity is preserved.
   - **Does not match base R's `set.seed()` + `rpois`/`runif` stream** — that is
     an accepted, deliberate trade for reproducible parallelism (pre-v1, no
     legacy stream to preserve).

3. **RNG contract = R↔C++ parity, via `dqrng` in both.** The pure-R reference
   and the C++ kernel must produce **byte-identical output for a given base
   seed**, both drawing from the same `dqrng` per-gamete streams. The prohibition
   is on RNG **engines**: no `std::mt19937`, no `rand()`, no ad-hoc `<random>`
   engine as the source of randomness — the engine is always `dqrng`. A
   `<random>`/boost *distribution* (e.g. `std::poisson_distribution`) fed by the
   `dqrng` engine would be permitted only if proven to match the R-side sampler —
   but decision #8 sidesteps this by using a uniform-only inversion sampler
   instead (see the Poisson sampler note). Base R's `R::rpois`/`unif_rand` are
   also out now — superseded by `dqrng` per decision #2.

4. **Scale target = large.** 10k+ offspring, 50k+ loci. Memory batching is a
   first-class requirement — a dense intermediate at this scale (~1 GB+) would
   OOM before CPU time matters.

5. **No `sample()` in the canonical algorithm.** The start-homolog choice is
   `u < 0.5` on a `dqrng` uniform draw on both sides. Avoids R 3.6's
   `R_unif_index` `sample.int` internals entirely. See the Poisson sampler note
   (decision #8) in the canonical-algorithm section for the crossover-count draw.

6. **C++ = one entry point per batch, internally decomposed.** A single
   `// [[Rcpp::export]]` `make_gametes_batch()` crosses the R↔C++ boundary
   **once per batch**; modularity lives in non-exported C++ helpers inside it.
   Rationale table in the "C++ kernel granularity" section.

7. **Parallelize across individuals, never per chromosome.** Rationale in the
   Parallelization section. The per-gamete RNG streams (#2) are what make this
   reproducible.

8. **Poisson via `dqrng`-uniform inversion; keep the position-based model.**
   *(Resolves open question #1.)* Crossover counts are drawn by a hand-written
   inversion/Knuth sampler consuming only `dqrng` uniforms — identical in R and
   C++ — so we retain exact genetic (cM) crossover positions (needed for #9) with
   no fragile cross-language Poisson matching. The uniform-only Markov
   reformulation is rejected because it has no continuous crossover location.

9. **Optional crossover-event storage, opt-in (default off).** `add_offspring()`
   gains `store_crossovers = FALSE`. When `TRUE`, each crossover drawn during
   meiosis is written as one row to a new `ind_crossover` table (see Optional
   feature). When `FALSE`, **zero cost** — the kernel does not emit the buffer.
   Rationale: crossovers are sparse (~1 per chromosome), so the table is ~0.06%
   the size of `ind_haplotype` — cheap when on, free when off.

10. **Physical `pos_bp` + a `genome_map` table + Morgans-native kernel.**
    `genome_meta.pos_Mb` → `pos_bp` (physical, `BIGINT`); the **genetic map moves
    to its own long table `genome_map(locus_name, sex, line_name, pos_cM)`**
    (mirrors `genome_effects`, so sex/line/version maps are rows, not schema
    changes). The recombination **kernel consumes genetic positions/lengths in
    Morgans** (`pos_cM / 100`), which the caller resolves **per gamete** for its
    `(parent sex, line)` (decision-#6 precedence in the migration plan) — the
    engine never sees bp/Mb. **The schema migration is a separate prerequisite
    project, fully specified in `plans/position_coordinates.md`** — merged and
    green before this plan starts. This plan assumes `genome_map` exists.

---

## The founder path — `define_founder_haplotypes()` + `add_founders()` (first-class here)

The **founder path is the other half of the slow write problem** and is in scope
for this plan, not an afterthought. It has three functions and a clear data flow;
**all three share the same wide-matrix → long-write hot spot** that
`add_offspring()` has, and none of the founder path involves meiosis.

### Data flow (grounded in code)

```text
define_founder_haplotypes()                add_founders()                      add_offspring()
  ── generates wide haplotype_matrix         ── reads founder_haplotypes (long)   ── reads parent ind_haplotype (long)
     (n_hap × n_loci) in R                      back into R                          ── recombines (gamete kernel)
  ── .write_founder_haplotypes()             ── reconstructs wide hap_data_matrix  ── writes ind_haplotype (long)
     materializes the FULL long frame           (n_hap × n_loci)
     (n_hap × n_loci rows) in R, then         ── sample() 2 haps/founder, w/ repl
     dbWriteTable → founder_haplotypes        ── builds wide hap_matrix + hap_wide
        (long pool: line_name,                ── per-chr UNPIVOT → ind_haplotype
         haplotype_id, locus_name, allele)       (long)
```

- **`define_founder_haplotypes()`** writes the **`founder_haplotypes` pool** (long:
  `line_name, haplotype_id, locus_name, allele`). It generates a dense wide
  `haplotype_matrix` (`R/define_founder_haplotypes.R:271-329`) and hands it to
  **`.write_founder_haplotypes()`** (`R/founder_haplotype_helpers.R:183-227`),
  which builds the **entire** long frame in R via `rep()` (one row per
  haplotype × locus) and `dbWriteTable`s it in one shot — **no batching, full
  materialization**. At 5k haplotypes × 50k loci that is a 250M-row R data.frame
  on top of the 1 GB wide matrix. This is a real hot spot and is **in scope**.
- **`add_founders()`** consumes that pool: it reads it back (`R/add_founders.R:155`),
  **reconstructs a second dense wide matrix** `hap_data_matrix`
  (`R/add_founders.R:178-180`), draws all haplotype indices in **one** `sample()`
  (`R/add_founders.R:200-205`, with-replacement — this is the population pool draw),
  builds a **third dense wide matrix** `hap_matrix` + `hap_wide`
  (`R/add_founders.R:244-257`), and writes `ind_haplotype` via the same per-chr
  `UNPIVOT` (`R/add_founders.R:279-308`) that `add_offspring()` uses.

**So the founder path round-trips a dense `n × n_loci` matrix through R three
times** (generate → write-long → read-long → reconstruct → select-wide → UNPIVOT),
and `add_founders()` holds **two** such matrices (`hap_data_matrix` +
`hap_matrix`) simultaneously. That is why founders are slow — same root cause as
`add_offspring()`'s wide intermediate (hot spot #3), doubled.

### No recombination for founders (confirmed) — founders are a *pool draw*

`add_founders()` has **no meiosis**: a founder's two haplotypes are **whole
haplotypes sampled with replacement** from the pool — the pool *is* the population
to draw from, and a founder simply gets two draws from it. There is **no
crossover, no genetic map, no per-gamete Poisson/dqrng stream** in founder
creation. So the gamete-generation machinery below (Stages 2–3: dqrng kernel,
Morgans map, inversion Poisson) **does not apply to founders** and must not be
forced onto them. A founder "gamete" is just a copied haplotype row.

### Shared core vs. alias — DECISION: share the write/assembly core, do **not** alias

**Do not alias `add_founders()` to `add_offspring()`.** They are distinct
user-facing functions with distinct signatures (`add_founders(get_table(
"founder_haplotypes"), n_males, n_females, line_name)` vs `add_offspring(pop,
matings)`) and distinct semantics (**sample-with-replacement pool draw** vs
**meiosis**). Aliasing would conflate two genuinely different operations and
confuse users; it also couples their argument surfaces forever. Rejected.

**Do share an internal core.** Both ultimately emit the *identical* long shape —
`(id_ind, parent_origin, strand = 1, line_origin, locus_id, locus_name, allele)`
— and both currently pay the wide→`UNPIVOT` tax to get there. Factor the common
part into **one seam** used by both:

> **The swappable seam is the "allele source"; everything downstream is shared.**
> - **`.write_long_haplotypes()`** (Stage 1): the batched, transaction-scoped
>   long insert into `ind_haplotype`. **Shared verbatim** by both functions.
> - The **"selected allele/`line_origin` block → long output vectors"**
>   materialization. **Shared** — it is the same packing for both.
> - The **allele source** differs and is the *only* difference:
>   - **founder source** — for each founder, copy sampled whole-haplotype rows
>     into the two `parent_origin` slots (`line_origin` = the founder's own line
>     for every locus). Zero crossovers, no map, no dqrng.
>   - **gamete source** — for each offspring, run the recombination kernel on the
>     two parents' haplotypes (Stages 2–3).

Concretely: express the founder source as the **degenerate zero-recombination
case** of the same long-materialization + write path. It feeds
`.write_long_haplotypes()` exactly as the gamete kernel does. This gives founders
`add_offspring()`-grade write efficiency (no dense wide matrix, batched, one
transaction) **without** dragging meiosis into founder creation. If Phase 0 shows
the founder gather/pack itself (not just the write) dominates at scale, the same
C++ long-materialization helper the kernel uses (the allele-gather + pack step,
**minus** the crossover draw) can serve the founder source too — a "no-crossover
gather" entry into the shared packer, decided by measurement, not up front.

**Scope caveat — the shared core is the autosome long write, not special
chromosomes.** The two functions handle non-autosomes **differently** and this must
**not** be unified blindly: `add_founders()` applies per-chromosome `copy_mode`
(`"full"`/`"half"`/`"none"`) **inline** in its write loop (`R/add_founders.R:279-308`),
writing only the applicable `(sex, parent_origin)` slots for a sampled haplotype;
`add_offspring()` routes special chromosomes through its **separate** `special_rows`
path (`R/add_offspring.R:447-503`, `:594-603`) with pass-through/recombining logic.
So the shared `.write_long_haplotypes()` covers the **autosome** long frame + the
generic "long frame → `ind_haplotype` INSERT" mechanism; each caller keeps (and
batches) its **own** special-chromosome emission and appends those long rows to the
same insert. Keep these seams distinct.

### What each function gets from this plan (efficiency parity)

| Function | Stage 1 (long write + batching + txn) | dqrng kernel (Stage 2) | C++ kernel (Stage 3) | Parallel (Stage 4) |
|---|---|---|---|---|
| `add_offspring()` | ✅ shared `.write_long_haplotypes()` | ✅ recombination | ✅ recombination | ✅ across offspring |
| `add_founders()` | ✅ **same** helper + batching + pre-draw RNG discipline (critique #8) | ✗ no meiosis | ▵ optional "no-crossover gather" only if Phase 0 flags the pack step | ▵ trivially parallel pool draw if ever needed |
| `define_founder_haplotypes()` | ✅ stream the pool write in batches (stop full-frame `fh` materialization in `.write_founder_haplotypes()`) | ✗ | ✗ | — |

The end state: **all three** founder/offspring write paths are dense-matrix-free,
batched, and transaction-scoped after Stage 1 — the founder path stops being the
simulation's slowest step. `add_offspring()` alone additionally gains the
Morgans-native recombination kernel (Stages 2–4).

---

## Current hot spots (grounded in code)

Ranked by expected ROI at the large-scale target. To be **confirmed by Phase 0
profiling** before committing effort.

1. **Per-offspring R gamete loop** — `add_offspring()` calls `make_gamete()`
   twice per offspring in a plain `for` loop (`R/add_offspring.R:427-503`, autosome
   calls at `:439`/`:443`); each call loops over chromosomes doing `rpois` /
   `sample` / `runif` / `findInterval` / `matrix[cbind(...)]`
   (`R/recombination_helpers.R:50-80`). For 10k offspring × ~30 chr × 2 parents
   that is ~600k R-level chromosome iterations of tight numeric work — the classic
   wrong-tool-for-R cost.

2. **Parent-matrix construction is O(parents × rows)** — the loop at
   `R/add_offspring.R:320-345` does `parent_haps_raw[parent_haps_raw$id_ind == pid, ]`
   (`:321`), a full `data.frame` scan **per parent**. Replace with one `split()` pass.

3. **Wide intermediate + `UNPIVOT`** — `hap_sire_mat`/`hap_dam_mat` +
   `line_origin` matrices (`R/add_offspring.R:416-419`) are dense
   `n_offspring × n_loci` integer/character matrices held in R (~1 GB each at
   target scale), assembled into wide frames (`:513-522`, `:565-574`) then
   unpivoted (`:577-588`). Removed entirely by decision #1.

4. **`special_rows` accumulation** — a fresh `data.frame(...)` per (offspring ×
   special chromosome) (`R/add_offspring.R:492`) then `do.call(rbind, ...)`
   (`:595`). Only matters when sex chromosomes/organelles are configured;
   preallocate.

5. **Bulk `INSERT` throughput** — once the generator is fast and wide is gone,
   the long insert of `n_offspring × n_loci × 2` rows (`R/add_offspring.R:577-588`)
   becomes the floor. Check `ind_haplotype` PK/constraint enforcement cost
   (`R/schema.R`) and consider DuckDB's Appender API vs `register` + `INSERT SELECT`.

6. **Founder path: three dense `n × n_loci` matrices + a full long
   materialization** (see "The founder path" above) — `define_founder_haplotypes()`
   builds the pool wide matrix (`R/define_founder_haplotypes.R:271-329`) then
   materializes the **whole** long frame in R (`R/founder_haplotype_helpers.R:206-225`);
   `add_founders()` reconstructs the pool matrix (`R/add_founders.R:178-180`) and
   builds the selected wide matrix (`:244-257`) before the per-chr `UNPIVOT`
   (`:279-308`). Same wide-intermediate root cause as #3, across the founder path.
   **Confirm in Phase 0 whether the founder path or the offspring path is the
   larger share of total simulation wall time** — the user reports founder
   creation as a top cost.

---

## Phase 0 — Profiling (do this FIRST, before any change)

Goal: a per-step time+memory breakdown so we optimize what actually costs, and a
baseline to measure every later stage against. **No optimization work starts
until Phase 0 numbers exist.**

### Tools

- **`profvis`** — line-level flame graph of a full `add_offspring()` call. This
  is the primary tool: it attributes time to individual source lines, so we see
  whether the loop, the parent-matrix build, or the DB write dominates.
- **`bench::mark()`** (preferred over `microbenchmark` — reports memory
  allocation + GC too) — for tight before/after comparison of individual
  extracted functions (e.g. the gamete kernel in isolation).
- **`Rprofmem()` / `bench`'s `mem_alloc`** — to catch the allocation churn from
  per-offspring vector/matrix allocation and the dense wide matrices.
- Existing harness: `dev/benchmarks/benchmark_haplotype_scale.R` (has `small` =
  5k loci / 500 founders always-on, `large` = 50k loci / 5k founders gated
  behind `TIDYBREED_BENCH_LARGE=1`).

**Dependency handling (critique #14).** `profvis` and `bench` are **dev-only** — they
do **not** become package `Suggests` (the benchmark is a `dev/` script, not part of the
package build). The benchmark script must therefore `requireNamespace()`-guard them and
**skip that measurement gracefully with a clear message** if absent, rather than error.
The existing harness (`benchmark_haplotype_scale.R:59-98`) is **broad storage
benchmarking**, not the per-step profile this plan needs — Phase 0 requires a **new
timing wrapper / instrumentation mode** that measures the named steps below separately.

### Instrumentation — split total wall time into named steps

Add the new per-step timing wrapper (above) to time these steps **separately** (wrap
each in its own timer, report as a table):

**`add_offspring()`:**

| Step | Code today | Why measure |
|---|---|---|
| Validation + ID assignment | `R/add_offspring.R:91-258` | Usually cheap; confirm it is |
| Gamete-map resolution + cache | `R/add_offspring.R:260-300` | Per-`(sex,line)` `resolve_genome_map()`; cheap if cached |
| Parent haplotype load (SQL) | `R/add_offspring.R:302-309` | I/O; scales with #parents × #loci |
| Parent matrix build (R loop) | `R/add_offspring.R:320-405` | O(P×N) scan (`:321`) — suspected hot spot #2 |
| **Gamete generation** | `R/add_offspring.R:427-503` | Suspected hot spot #1 |
| `line_origin` gather | `R/add_offspring.R:441,445` | Matrix `cbind` indexing per offspring |
| Long-row assembly | (new — replaces wide build `:513-522`) | New critical path under decision #1 |
| DuckDB write / INSERT | `R/add_offspring.R:559-603` | Suspected floor after kernel is fast |

**Founder path (profile these too — hot spot #6):**

| Step | Code today | Why measure |
|---|---|---|
| Pool generation (wide matrix) | `R/define_founder_haplotypes.R:271-329` | LD methods are O(n_hap × n_loci) |
| Pool long materialization + write | `R/founder_haplotype_helpers.R:206-225` | Full long frame in R — suspected founder hot spot |
| Pool read-back + matrix rebuild | `R/add_founders.R:155-180` | Second dense matrix from long |
| Founder sample() | `R/add_founders.R:200-205` | Single call; expected cheap |
| Selected wide build + `UNPIVOT` write | `R/add_founders.R:244-308` | Third dense matrix + per-chr UNPIVOT |

### Scenarios to profile (must stress the real shape)

Simple 1:1 pairings hide the per-offspring-loop cost. Profile at least:

- **Founder creation** (the founder path, hot spot #6): `define_founder_haplotypes()`
  then `add_founders()` at a large pool (e.g. 5k haplotypes → 10k founders × 50k loci).
  Measure both functions' per-step tables (above) so we know whether the founder
  path or the offspring path is the bigger share of total wall time — the user
  reports founders as a top cost, so this is not optional.
- **Litter/hatch shape**: few parents, many offspring each (e.g. 20 sires × 100
  dams, 1 offspring each → 2000; and 1 sire × 2000 offspring for the extreme).
- **Target scale**: 10k offspring, 50k loci, ~30 chr (env-gated `large`).
- With and without a `chr_meta`-configured sex chromosome (exercises the
  special-chromosome branch).

### Deliverable

A committed baseline table (per-step time + peak memory, at `small` and `large`,
for the founder path **and** each mating design) as a fixture/comment in the
benchmark script, plus a short written conclusion: **which steps dominate across
both the founder and offspring paths**, so Stages 1–3 are ordered by measured ROI
rather than assumption.

---

## Target architecture (after all stages)

```text
add_offspring(pop, matings, store_crossovers = FALSE)
  ├─ validate + resolve parent sex/ploidy + offspring IDs          (unchanged)
  ├─ load parent haplotypes once (SQL)                             (unchanged query)
  ├─ build parent haplotype/line_origin structures  ── split(), O(N) once
  ├─ FOR each batch of B offspring:                 ── memory bound (Stage 1)
  │     ├─ gametes <- make_gametes_batch(..., store_crossovers)  ── R ref OR C++ kernel
  │     │      returns LONG parallel vectors (autosomes):
  │     │        parent_origin, locus_idx, allele, line_origin_code
  │     │      + optional ragged crossover buffer (when store_crossovers)
  │     ├─ append special-chromosome rows (R path)   ── same long shape
  │     ├─ caller adds id_ind; maps locus_idx→locus_id/locus_name;
  │     │   sets strand = 1L
  │     ├─ register long frame + INSERT ... SELECT   ── no wide, no UNPIVOT
  │     ├─ if store_crossovers: INSERT crossover rows → ind_crossover
  │     └─ free batch
  └─ (special chromosomes handled in same long emit; see below)
```

The generator is the single swappable seam. It returns **kernel-local** columns
(`locus_idx` = 1..n_autosome_loci, `line_origin_code` = integer dictionary code);
the caller attaches `id_ind`, resolves `locus_idx → locus_id`/`locus_name` via
`genome_meta`, sets `strand = 1L`, and decodes `line_origin_code → line_origin`.
Everything below the generator (DB write) is shared with `add_founders()`.

**Naming (used consistently throughout this plan):** kernel inputs/outputs use
`parent_1_idx`/`parent_2_idx` (the canonical `parent_1`/`parent_2` schema names,
*not* the `sire`/`dam` alias), `parent_origin ∈ {1,2}` (1 = parent_1, 2 =
parent_2 — matches `ind_haplotype.parent_origin`), `locus_idx` (kernel-local),
and `line_origin_code`. DB column names are unchanged from the schema in
`CLAUDE.md`.

---

## Data path & memory model — where the bytes actually live (R ↔ C++ ↔ DuckDB)

This pins the single most important implementation question: **the C++ kernel does
NOT read or write DuckDB tables directly. It never opens a connection, never sees a
table. It is a pure numeric function over in-memory arrays.** R owns all DB I/O; C++
owns all computation. Crucially, "going through R" here is **not** a
serialization/copy tax — it is zero-copy pointer passing. The whole data path per
batch:

```text
  DuckDB ind_haplotype/founder_haplotypes
     │   (1) ONE SQL read per call: DBI::dbGetQuery(...)  ── parents << offspring, cheap
     ▼
  R integer vectors / packed parent matrix           ── R-allocated memory
     │   (2) passed to make_gametes_batch() — Rcpp hands C++ the RAW POINTER
     │       (INTEGER(x) / IntegerVector::begin()) to R's memory. NO COPY.
     ▼
  C++ kernel: loops over loci/chromosomes on int*    ── the only heavy compute
     │   (3) writes results into PREALLOCATED R output vectors via raw pointer. NO COPY.
     ▼
  R long output vectors (parent_origin, locus_idx, allele, line_origin_code)
     │   (4) duckdb_register(conn, "__tmp", frame) — DuckDB's scan reads these R
     │       column vectors IN PLACE; INSERT ... SELECT appends to ind_haplotype.
     ▼
  DuckDB ind_haplotype   ── ONE write into table storage (unavoidable for any insert)
```

**The only two data movements are the unavoidable ones:** read parents out of DuckDB
once (step 1), write new rows into DuckDB once (step 4). Everything between is pointer
passing. Three distinctions that resolve the common confusion:

1. **R the interpreter is not in the hot path.** R makes a handful of orchestration
   calls per *batch*; the per-locus loops are 100% C++. R's interpreter overhead is
   paid ~once per batch, never per locus.
2. **R the allocator is not a copy layer.** Inputs and outputs live in R-allocated
   memory; Rcpp gives C++ a raw `int*` into it, so C++ reads inputs and fills outputs
   **in place**. No marshalling, no SEXP-per-value conversion for the numeric columns.
3. **The R→DuckDB handoff is zero-copy on read.** `duckdb_register()` exposes the R
   output vectors as a virtual table; DuckDB's dataframe scan reads the R integer/double
   columns **directly** during `INSERT ... SELECT` (no serialization, no temp file, no
   second buffer). The INSERT then writes into the table's column store — that write is
   intrinsic to *any* insert and is not overhead we introduced. This is why
   `line_origin` crosses the C++ boundary as an **integer code** (Stage 1 encoding): so
   the heavy columns are all integer and stay on the fast, conversion-free scan path;
   the string decode happens once, in SQL, at insert.

### Why C++ does not write to DuckDB directly (and why that's optimal, not a compromise)

DuckDB has a C++ `Appender`, and because the `duckdb` R package *is* libduckdb in the
same process, one *could* reach into its connection external-pointer and append from
our Rcpp code. **Rejected**, deliberately:

- It requires `LinkingTo: duckdb` and an **exact ABI match** with the user's installed
  `duckdb` package — version-fragile, breaks on `duckdb` updates, CRAN-hostile.
- The win is **marginal**: it would save allocating the per-batch R output vectors and
  nothing else. Steps (1) and (4) — the actual DB reads/writes — cost the same whether
  driven from R or C++, because the data still has to leave/enter DuckDB's storage.
- Batching (Stage 1) already bounds the R output-vector memory, and `duckdb_register()`
  already makes the read zero-copy, so the thing C++-direct would optimize is already
  small. Not worth trading the portable `duckdb`-package abstraction for it.

**The "Appender API vs register+INSERT SELECT" A/B (Stage 1 deliverable) is between two
R-side DuckDB write APIs** — `duckdb`'s appender vs `duckdb_register`+`INSERT SELECT` —
**not** C++-direct-to-DuckDB. C++-direct is out of scope (see Explicitly out of scope).

### What this means for how the Rcpp is written

`make_gametes_batch()` is a `// [[Rcpp::export]]` function whose **entire signature is
plain numeric arrays** (packed parent alleles, chr position/length vectors in Morgans,
parent index vectors, base seed, flags) and whose **return is a list of plain numeric
vectors** (the long output). No `tidybreed_pop`, no `DBI`, no `data.frame`, no
connection handle ever enters C++ — the same seam constraint Stage 0 already imposes.
That is exactly what makes it (a) portable, (b) unit-testable in isolation against the R
reference for parity, and (c) zero-copy at both boundaries. The DuckDB read, the
`id_ind`/`locus_id`/`line_origin` resolution, and the INSERT all stay in R, wrapped in
the one per-call transaction.

---

## Position coordinates & the genetic map (decision #10)

**This was a separate prerequisite project — see `plans/position_coordinates.md`
for the full specification — and it is now IMPLEMENTED and green in the working
tree.** It moved physical position to `genome_meta.pos_bp` (`BIGINT`) and the genetic
map to the long table `genome_map(id_genome_map, locus_id, locus_name, sex, line_name,
map_name, pos_cM)`, updated `define_genome()`, `define_founder_haplotypes()`, and the
recombination readers, and added `resolve_genome_map()`/`validate_genome_map()` plus the
per-`(sex,line)` gamete-map cache in `add_offspring()`. Swapping the readers to the
resolved default map was value-neutral (same positions); pre-1.0.0 there is **no**
cross-version output guarantee — only forward reproducibility (same seed reproduces on
the new code). This plan builds directly on that merged work.

What this plan needs from it, and nothing more:

- The `genome_map(locus_name, sex, line_name, pos_cM)` table exists (the genetic
  map), plus `genome_meta.pos_bp` (physical, unused by the engine).
- **Kernel is Morgans-native, map resolved per gamete:** the caller resolves the
  map for each gamete's `(parent sex, line)` (decision-#6 precedence), converts
  `pos_cM → Morgans` (`/100`), and passes `chr_pos_M`/`chr_len_M` to
  `make_gametes_batch`. The engine never sees bp/Mb. **Design consequence:** sire
  and dam gametes may use *different* maps (sex-specific / achiasmy), so the
  kernel takes the map layout as an argument rather than assuming one shared map.
  **This caller seam already exists** in the current R `add_offspring()` (post the
  coordinate migration): it caches a resolved descriptor set per unique
  `(parent_sex, parent_line)` and each gamete uses its producing parent's map. This
  plan replaces the *per-gamete R loop* with the batch kernel but keeps that
  per-`(sex,line)` resolution/caching. In v1 all gametes still resolve to the single
  default map, so no map differs yet — but the resolution is genuinely exercised.
- Crossover positions are drawn and stored in cM (`ind_crossover.pos_cM`).

---

## The canonical recombination algorithm (the parity contract)

Same Haldane model as today. The **unit of randomness is the gamete** =
`(offspring_index o, parent_role r ∈ {1,2})`. Each gamete draws from its **own
`dqrng` sub-stream**, seeded deterministically from `(base_seed, o, r)` (e.g.
`dqset.seed(base_seed)` then a per-gamete jump/stream id, or a hashed seed). This
is what decouples reproducibility from iteration order and thread count:
iterating by parent, by offspring, or across N cores all yield identical output,
because gamete `(o, r)` always consumes stream `(o, r)` from its start.

### Seed API — the public reproducibility contract (critique #4)

`add_offspring()` currently uses base R's global RNG via the caller's `set.seed()`.
The per-gamete `dqrng` streams need a **base seed**. **Decision (public API):**

- **`add_offspring(pop, matings, ..., seed = NULL)`.**
- **`seed = NULL` (default):** draw one integer base seed from the **current base-R
  RNG** (a single `sample.int(.Machine$integer.max, 1L)`-style draw). Consequence: an
  upstream `set.seed(k)` still fully determines the run — the same `set.seed(k)` before
  the same `add_offspring()` call reproduces, satisfying "same global seed reproduces on
  current code." This keeps one seeding idiom for users across `add_founders()` (base R)
  and `add_offspring()` (dqrng).
- **`seed = <int>`:** use it directly as the base seed, bypassing the global RNG draw —
  for callers who want a gamete stream independent of surrounding base-R state.
- **Audit / storage:** by default the resolved base seed is **not** persisted to a
  table (CLAUDE.md's "disdain for storing metadata"); it can be surfaced to the caller
  (e.g. a message or invisible return) for note-taking. **Open question #2 (below):** if
  `store_crossovers = TRUE`, do we persist the base seed alongside the crossover rows for
  exact replay? Recommendation: yes — a single cheap `base_seed` scalar/column makes an
  `ind_crossover` run self-describing; decide at Stage 2 (where the seed API lands).
- **`add_founders()`** keeps base R `sample()` for founder haplotype selection (no
  meiosis, no dqrng) — see the founder RNG note in Stage 1.

### dqrng stream-derivation scheme — pin before Stage 2 (critique #5)

The rest of the design assumes independent, collision-free per-gamete streams; this must
be **resolved before the C++ kernel**, not during it. **Decision + spike:**

- **Stream key tuple:** `(base_seed, global_offspring_index o, parent_origin r,
  stream_kind)` where `stream_kind ∈ {"autosome", "special"}` (the two must not share a
  counter — see Special chromosomes). `o` is the **original `matings` row index**, never
  a batch- or group-local index.
- **Derivation — RESOLVED by the spike (2026-07-07), simpler than the original
  hash64 proposal.** `dqrng::dqset.seed()` accepts only a **scalar or 2-int seed and a
  scalar stream** (a length-≥3 vector errors `"out-of-range seed"`), so the tuple is
  folded to a single **stream id** `sid = ((o*2 + (r-1))*2 + (kind-1))`
  (`kind`: 1 = autosome, 2 = special; `r` = `parent_origin` ∈ {1,2}). No hand-rolled
  `splitmix64`/`hash64` is needed — **dqrng's own `convert_seed()` does the 64-bit
  folding**, so feeding the identical integer vector through it in R and C++ yields a
  bit-identical generator seed by construction. **Adopted (Scheme A):**
  - **R reference (Stage 2):** `dqrng::dqRNGkind("Xoroshiro128++"); dqrng::dqset.seed(
    c(base_seed, sid)); dqrng::dqrunif(n)`.
  - **C++ kernel (Stage 3):** `dqrng::convert_seed<uint64_t>(c(base_seed, sid))` →
    `dqrng::random_64bit_wrapper<dqrng::xoroshiro128plusplus>` → `uniform01()`
    (`= (x >> 11) * 2^-53`, which the spike confirmed `dqrunif` uses for `[0,1)`).
  - **Bound:** `sid` fits signed-32-bit for `o < ~5e8` (>> realistic `n_offspring`);
    if ever exceeded, split `sid` into a 2-int stream — `convert_seed` handles it.
  - Scheme B (`dqset.seed(base_seed, sid)` scalar stream → C++ `seed(base)+jump(sid)`)
    was **also** proven bit-identical and is the fallback; Scheme A is preferred (no
    jump semantics, everything through `convert_seed`).
- **Pre-implementation spike — DONE and GREEN** (`dev/spikes/dqrng_parity/`,
  dqrng 0.4.1): Scheme A and Scheme B both `ALL IDENTICAL` across 6 tuples × {1,5,50}
  draws; distinct keys give distinct streams (`ALL DISTINCT`); `dqrunif` matches C++
  `uniform01()`. This retires the foundational R↔C++ RNG risk; the kernel (Stage 3) is
  built on Scheme A.

**Within one gamete's stream**, the draw order across `C` chromosomes (ascending
`chr`) is fixed and identical in R and C++ — this per-gamete order is the parity
contract the golden test enforces. Everything here is in **Morgans** (genetic
distance from `pos_cM / 100`, decision #10). Per chromosome `c` of genetic length
`len_M_c` Morgans, the Poisson rate is simply `lambda_c = len_M_c` (expected
crossovers = genetic length in Morgans):

1. **Crossover count** `n_cross ~ Poisson(len_M_c)`, drawn by an
   **inversion/Knuth sampler on `dqrng` uniforms** (see below).
2. **Start homolog** `start = (u < 0.5) ? 1 : 2`, one `dqrng` uniform `u`.
3. **Crossover positions** `n_cross` uniforms scaled to `U(0, len_M_c)` Morgans,
   then **sorted** (sorting consumes no RNG). Drawing uniform in **genetic**
   distance (not physical Mb) is the correct Haldane model and is what the
   genetic map buys us. These cM positions are the real crossover locations —
   retained for the deterministic step and **optionally persisted** to
   `ind_crossover` (see Optional feature).

Then, **no more RNG** — allele/homolog assignment is deterministic:
`findInterval(locus_cM, sorted_crossovers_cM)` parity → homolog index per locus →
gather allele from the parent's two homologs and `line_origin` from the same
homolog. Deterministic part is vectorizable in R and a tight loop in C++.

> **Chromosome genetic length `len_M_c` — pin the convention (post-migration).**
> The current `build_chr_info()` uses `len_cM = max(pos_cM)` for the chromosome
> (the last **marker's** genetic position), so `lambda_c = max(pos_cM)/100`.
> Crossovers beyond the last marker do not change any stored haplotype, so this is
> the correct **observed-marker inheritance** length for the allele output. It is
> **not** the full biological chromosome length unless the last locus sits at the
> chromosome end. Consequence for `ind_crossover` (decision #9): with this
> convention the stored crossovers are exactly those that fall within the marker
> span — i.e. the crossovers that can affect represented loci, **not** all
> biological crossovers (nothing is drawn in the terminal `[max(pos_cM),
> chr_end_cM]` interval). If the crossover-storage feature is meant to reflect all
> biological recombination, this plan must add an explicit chromosome-end genetic
> length (a per-chromosome terminal cM, sourced from `chr_meta` or a map endpoint)
> and draw `lambda_c` from it. **Default decision unless changed:** keep the
> `max(pos_cM)` marker-span convention and document that `ind_crossover` stores
> marker-span crossovers only.

> **Poisson sampler — RESOLVED (open question #1), refined for numerical safety
> (critique #11).** We keep the **Poisson-count + uniform-position** model (not the
> uniform-only Markov reformulation), because storing crossover **locations (in cM)**
> requires continuous positions, which the Markov switch-per-locus model does not
> produce. Parity is preserved by drawing the Poisson count with a **fixed
> inversion/Knuth algorithm that consumes only `dqrng` uniforms**. **Use
> log-accumulation, not the raw-product form:** accumulate `s += log(u_k)` and count
> steps while `s > -lambda_c` (equivalently `k` such that `sum_{j<=k} log(u_j) <
> -lambda_c`). This is algebraically identical to "multiply uniforms until the product
> drops below `exp(-lambda_c)`" but **does not underflow** when `exp(-lambda_c)` is
> denormal (large `lambda_c`). Implemented identically in the R reference and the C++
> kernel, it is a uniform-only sampler (so `dqrng`'s identical R↔C++ uniform stream
> guarantees parity) **and** yields the exact count, after which the `n_cross` position
> uniforms are drawn. This gives exact genetic (cM) positions *and* removes the fragile
> "native Poisson must match across languages" risk.
>
> **Parameter bounds & guard.** Cost is `O(lambda_c)` uniforms per chromosome. For
> livestock-scale maps `lambda_c` ≈ 1 (chromosome length in Morgans), so this is
> negligible — but `cM_per_Mb` is user-specified and the package aims to be
> species-flexible, so `lambda_c` can be large. **Decision:** document a supported
> `lambda_c` ceiling for the exact inversion path (target: correct and fast to
> `lambda_c ≈ 30`, i.e. a 3000 cM chromosome), and add a guard that errors (or, if we
> later choose, switches to a parity-safe alternative — PTRS/transformed rejection
> proven bit-identical R↔C++) above it rather than silently spinning `O(lambda_c)`.
> Both R reference and C++ kernel use the **same** hand-written log-accumulation
> sampler (not base `rpois`, not `std::poisson_distribution`) so they match. Tests must
> cover `lambda = 0`, small (`~0.3`), typical (`~1`), and a stress value near the
> ceiling — see Testing strategy.

### Special chromosomes (sex chromosomes / organelles)

Handled on the **R path** (not the C++ kernel in v1) and emitted into the **same
long output**. Non-recombining / single-copy inheritance (`pass_through_gamete`,
`R/recombination_helpers.R:107`, row choice `sample.int(k, 1)` at `:109`)
consumes RNG only when `k == 2` (one uniform `< 0.5` row choice) — keep that rule.

**Stream independence (important under per-gamete streams):** because autosomes
are drawn in C++ and special chromosomes in R, the two must **not** share one
RNG stream position — otherwise the R side would have to know where the C++
kernel left the counter (fragile across the boundary). Instead, key the two on
**independent** per-gamete sub-streams: autosomes use `(base_seed, o, r,
"autosome")`, special chromosomes use `(base_seed, o, r, "special")`. Both are
deterministic and reproducible; neither depends on the other's draw count. This
replaces the old "special strictly after autosome" ordering discipline
(`R/add_offspring.R:447`, comment at `:218`), which existed only because both shared base R's
single global stream — a constraint the per-gamete stream model removes. These
chromosomes are low-locus-count; correctness first, speed second.

**Special-chromosome RNG under dqrng (critique #7).** After Stage 1 the special path
must **also** move to `dqrng` (not base `sample.int()`/`runif()`), so it draws from the
`"special"` sub-stream and inherits thread-/batch-invariance. Concretely:

- **`pass_through_gamete()` row choice** (`k == 2`, non-recombining special chr): the
  homolog pick becomes `u < 0.5` on **one `dqrng` uniform** from the `"special"` stream,
  not `sample.int(2, 1)`. `k == 1` still draws nothing (strictly hemizygous inheritance
  consumes no RNG — keep that rule).
- **Recombining special chromosome** (`k == 2 && recombines_*`, e.g. the
  pseudoautosomal-style case): use the **same canonical algorithm** as autosomes — the
  log-accumulation Poisson count + uniform cM positions + `findInterval` parity — but run
  it in the shared R code on the `"special"` stream, keyed per gamete.
- **Crossover storage:** recombining special chromosomes emit `ind_crossover` rows too
  when `store_crossovers = TRUE`, with their crossover cM positions from the `"special"`
  stream. Non-recombining special chromosomes emit none (nothing is drawn).
- **Parity scope (important):** the C++ kernel is **autosome-only**. "R↔C++ parity on a
  special-chromosome config" therefore means **wrapper-level equality** — autosomes come
  from C++, special chromosomes from the shared R code, and the *assembled*
  `ind_haplotype` (and `ind_crossover`) match between the all-R reference run and the
  C++-autosome + R-special run. It does **not** mean the C++ kernel itself handles
  special chromosomes. The parity test fixture must include at least one configured sex
  chromosome to exercise this seam.

---

## Optional feature — crossover-event storage (`ind_crossover`)

Opt-in via `add_offspring(..., store_crossovers = FALSE)` (decision #9). When
`TRUE`, every crossover drawn during meiosis is recorded as one row. This lets
users verify recombination: count crossovers per gamete, check their genetic
(cM) distribution, and plot them.

**Why it is cheap:** the positions are *already drawn* (Block 3) — we emit
instead of discard. Crossovers are sparse (Poisson mean = chromosome genetic
length in Morgans ≈ 1 per 100 cM), so at the 10k-offspring × 30-chr target the
table is ~600k rows vs `ind_haplotype`'s ~1e9 (**~0.06%**). Off ⇒ the kernel
skips the buffer entirely (zero cost). On ⇒ a second small `INSERT`.

### `ind_crossover`

Long format, one row **per crossover event**. Absence of a row for a
`(id_ind, parent_origin, chr)` means that gamete's chromosome did not recombine
(no zero-rows stored — consistent with the "don't store derivable/absent data"
principle). Includes crossovers on recombining special chromosomes
(`k == 2 && recombines`) when `store_crossovers = TRUE`.

| Column        | Type     | Notes                                                        |
|---------------|----------|-------------------------------------------------------------|
| id_crossover  | INTEGER  | Primary key, assigned R-side via `next_int_id(conn, "ind_crossover", "id_crossover")` — **not** DB auto-increment (matches `id_genome_map` precedent) |
| id_ind        | VARCHAR  | FK to `ind_meta.id_ind` — the offspring                     |
| parent_origin | UTINYINT | 1 = event in parent_1's (sire's) meiosis, 2 = parent_2's (dam's); matches `ind_haplotype.parent_origin` |
| chr           | INTEGER  | Chromosome number                                           |
| chr_name      | VARCHAR  | Chromosome name                                             |
| pos_cM        | DOUBLE   | Crossover **genetic** location in cM — the exact drawn position (decision #10) |

**Reserved**: all columns (managed by `add_offspring()`).

**Schema registration (critique #9).** `ind_crossover` follows the same package
conventions as every other system table — it is **not** ad-hoc. Stage 1 must:

- create the empty table in **`define_genome()`** (the genome-setup surface), so
  `restore_pop()` sees it and it exists before the first opt-in insert. (Note:
  `initialize_genome()` was **removed** in v0.51.0 — `define_genome()` is the only
  genome-setup entry point.)
- register it in `R/sql_utils.R` — `TABLE_RESERVED_COLS` (all columns reserved),
  `TABLE_PRIMARY_KEYS` (`id_crossover`), `TABLE_ROW_KEYS`, and `SYSTEM_TABLES`;
- add its description to the schema metadata in `R/schema.R`;
- include it in `pop$tables` so `get_table("ind_crossover")` works.

Assign `id_crossover` per batch via `next_int_id()` (`COALESCE(MAX(id),0)+1` R-side),
consistent with `id_genome_map` from the coordinate migration.

**On the unit:** crossovers are drawn uniform in *genetic* distance, so `pos_cM`
is the exact, native location. A **physical `pos_bp`** is *derivable* by linearly
interpolating between the flanking loci's `pos_bp`/`pos_cM` — under a constant
rate it's a simple scale; under a variable-rate map it needs the local
interpolation. Left out of the stored table for now (avoid storing derivable
data); revisit if users want physical crossover coordinates materialized.

Aggregation examples (the user's stated use cases):

```sql
-- crossovers per gamete per chromosome
SELECT id_ind, parent_origin, chr, COUNT(*) AS n_crossovers
FROM ind_crossover GROUP BY id_ind, parent_origin, chr;

-- distribution of crossover positions on chr 1 (for plotting)
SELECT pos_cM FROM ind_crossover WHERE chr = 1;
```

**Kernel impact:** the C++ kernel gains a **ragged second output** (variable
crossover count per gamete) — long format handles this naturally, and it is
gated on the `store_crossovers` flag so the default path allocates nothing for
it. Under per-gamete streams the crossover rows are as parallel-safe as the
alleles (each gamete owns its events).

---

## C++ kernel granularity — one entry point, internal helpers

**Decision: a single `// [[Rcpp::export]] make_gametes_batch()` per batch,
decomposed internally into non-exported C++ helpers** (`draw_crossovers()`,
`assign_homologs()`, `gather_alleles_and_lo()`). Modularity lives *inside* C++,
where it is free; the R↔C++ boundary is crossed exactly once per batch.

| Concern | Many exported C++ fns (called from an R loop) | One entry point (internal C++ helpers) |
|---|---|---|
| Boundary crossings | One per call — per-chr/per-offspring calls reintroduce the R-loop overhead we are removing | **Once per batch**; SEXP marshalling + RNG-engine setup paid once |
| Intermediates (crossover positions, homolog vectors) | Become R objects → allocation + GC churn, large at 50k loci | Stay on the C++ stack/heap → cache-friendly, no SEXP alloc, no GC |
| RNG streams | Multiple engine enter/exit; easier to mis-order draws | Streams created/consumed in one place → cleanest parity guarantee |
| Memory peak | Forces intermediate/wide materialization across boundary | Writes straight into preallocated long output vectors |
| Testability | Each piece callable from R | Whole kernel via the R↔C++ parity test; helpers via a small C++ harness if needed |

What stays in R (I/O-bound, DuckDB's job — not ported): SQL parent load, the
`split()` structuring of parent haplotypes, and the long `INSERT`. What the one
kernel does: consume per-gamete streams → produce the long parallel output
vectors. The kernel operates on **one bounded batch** (Stage 1 caps it at `B`
offspring), so "one function" never implies "whole population in memory."

---

## Parallelization — across individuals, never per chromosome

**Decision: parallelize across individuals (gametes); do not parallelize
chromosomes within a gamete.** Four independent reasons:

1. **Task count / scalability.** ~10k+ offspring vs ~30 chromosomes.
   Across-individual → thousands of independent tasks, scales to any core count.
   Per-chromosome → capped at ~30 tasks.
2. **Load balancing.** Chromosomes differ in length/locus-count and draw a
   random crossover count; per-chr tasks are uneven. A whole gamete/individual is
   a balanced fixed chunk.
3. **Output locality / false sharing.** Each individual owns a contiguous block
   of the long output (`2 × n_loci` rows) → disjoint per-core regions.
   Per-chr-across-individuals interleaves writes into one individual's range →
   false sharing.
4. **Per-task overhead.** A single chromosome is trivial work; dispatch overhead
   would swamp it. A gamete/individual amortizes it.

**Why it is reproducible:** the per-gamete `dqrng` streams (decision #2) mean
gamete `(o, r)` always uses stream `(o, r)`, independent of which core runs it or
in what order — so output is identical for any thread count. This is the whole
reason we adopted per-gamete streams now rather than base R's global,
thread-unsafe stream.

**Parent data under parallelism (the "don't re-grab the same parent" point):**
each unique parent's haplotype + `line_origin` matrices are loaded **once** into
**read-only, shared** structures (the `split()` fix in Stage 0). A sire mated to
500 dams is *read* 500 times, never mutated → no lock, no re-fetch, zero
contention across cores. Iterating/grouping generation **by parent** also keeps
that parent's matrix hot in cache (a 50k-locus sire matrix ≈ 400 KB → better L3
reuse than interleaving sire/dam per offspring). Parent-grouping (data access)
and per-gamete streams (randomness) are orthogonal and coexist cleanly: group by
parent for locality while each gamete still draws from its own `(o, r)` stream.

**Toolchain — DECISION: `RcppParallel`, not OpenMP (install portability).**
`dqrng` exposes its engine to C++ via `dqrng.h`; drive the threading with
**`RcppParallel::parallelFor` over the offspring range**. Single- and
multi-threaded runs must pass the *same* golden test (thread-count invariance,
test 2). This is **in scope** but is the **last** stage, scheduled only after the
single-threaded kernel and write path are correct and parity-verified (Stage 4
after Stage 3).

> **Why `RcppParallel` and not OpenMP — this is a portability decision, not a
> preference.** Raw **OpenMP is the classic "won't install on macOS" trap**:
> Apple's default clang (Xcode Command Line Tools) ships **no `libomp`**, so a
> `devtools::install_github()` source build with `$(SHLIB_OPENMP_*FLAGS)` **fails
> on a stock Mac** unless the user has separately installed `libomp` — exactly the
> friction we must avoid for a GitHub-distributed package. **`RcppParallel` bundles
> Intel TBB** and handles platform detection, so it compiles cleanly on
> macOS/Windows/Linux with only the standard R toolchain — **no system library, no
> OpenMP**. It is the ecosystem-standard, CRAN-safe choice for portable multicore.
>
> **Zero-dependency alternative:** `std::thread` (C++11, already required by the
> toolchain) partitioning the flat offspring range needs **no extra package at
> all** — viable here because the kernel is a pure numeric function that never calls
> the R API off-thread (the one hard rule for off-main-thread code). If
> `RcppParallel` ever proves unavailable on a target platform, fall back to
> `std::thread`; the single-threaded Stage-3 kernel is always the ultimate fallback.
> **OpenMP is rejected as a hard requirement** on macOS-install grounds (consistent
> with CLAUDE.md's compiled-code policy: no required system libraries).

### Install portability (why Stage 4 does not hurt GitHub installs)

The whole point of layering parallelism last is that **the package is fully
functional and universally installable long before any threading dependency
enters**:

- **Stages 0–2:** pure R + `dqrng` (an ordinary CRAN R package). Installs
  everywhere, no compiler beyond what R itself needs.
- **Stage 3 (single-threaded C++):** deps `Rcpp`, `dqrng`, `BH`, `sitmo` — **all
  header-only or standard CRAN packages, no system libraries, no OpenMP**. Compiles
  from GitHub source with only the toolchain every source-install user already has
  (Rtools / Xcode CLT / build-essential). This is the same install surface as any
  Rcpp package.
- **Stage 4 (parallel):** adds `RcppParallel` (bundled TBB, no system lib) **or**
  `std::thread` (nothing). Still the big-three, still no OpenMP.

So the only way installs break is requiring OpenMP or a system library — both
**forbidden** here. Binary channels (r-universe / CRAN) need no compiler at all. The
R reference generator (Stage 2) also remains as a runtime fallback if the compiled
object is ever unavailable.

---

## Stage sequence (revised — critique #2, #3)

Execution order, authoritative:

| Stage | What | Kernel | Write path |
|---|---|---|---|
| **0** | Extract gamete seam + cheap wins (`split()`, preallocate) | base-R `make_gamete()` | **wide (unchanged)** |
| **1** | Direct long write + batching + `ind_crossover` schema | base-R `make_gamete()` | **long, batched** |
| **2** | dqrng R reference kernel (intentional stream change) | dqrng R, long output | long, batched |
| **3** | C++17 kernel + R↔C++ parity | dqrng C++ | long, batched |
| **4** | Parallelize across individuals | dqrng C++ (threaded) | long, batched |

The key reorder from the first draft: **the long-write + batching path (Stage 1) lands
before any long-output kernel is wired into `add_offspring()`.** This gives a stable,
bounded memory ceiling *before* the dqrng kernel starts allocating long vectors for
large batches, and it means the dqrng kernel (Stage 2) drops into an already-long write
path rather than forcing an awkward re-widen. Stage 0 keeps the existing wide write
untouched (pure structural refactor); Stage 1 converts the write path to long while
still driving allele values from base-R `make_gamete()` (so only the write *shape*
changes, not the values); Stage 2 swaps in the dqrng kernel (the one intentional
RNG-stream change).

> **Phase 0 finding (2026-07-07) — reprioritization.** The committed baseline
> (`dev/benchmarks/profile_haplotype_paths.R`;
> `plans/optimize_add_offspring_stage_0_summary.md`) showed **gamete generation is
> only ~2–6% of `add_offspring()` time**; the dominant costs are the parent scan
> (now fixed in Stage 0) and the **wide→`UNPIVOT` DB write** (34–70% of
> `add_offspring()`, ~69% of `add_founders()`, ~96% of the founder-pool write),
> which also drove ~43 GB of allocation churn for 200 offspring. **Consequence:
> Stage 1 (long/batched write) is the highest-ROI next step and the priority; the
> dqrng/C++/parallel kernel (Stages 2–4) is a *later, smaller* optimization of a
> <5% slice.** The parity spike being done keeps Stages 2–3 unblocked whenever we
> want them, but they no longer lead. **Stage 0 (profile + seam + `split()` +
> preallocation) is complete and green** (output-neutral: parity golden + seam-
> equivalence proof; full suite 1229 pass / 0 fail; O(P×N) scan ~12× faster, churn
> ~28× lower).

---

## Stage 0 — Extract the gamete seam + cheap wins (structural, output-neutral)

Prerequisite for everything. Pull gamete generation out of `add_offspring()`
into a standalone, dependency-free function — **no `tidybreed_pop`, no DBI, no
data.frame** in its signature — so the R and C++ versions are drop-in
swappable and independently testable/benchmarkable. **In Stage 0 the seam wraps
today's base-R `make_gamete()` and still feeds the current wide write path**
(`UNPIVOT` unchanged) — this stage is a pure structural refactor, output-neutral by
construction.

**Target kernel interface (finalized in Stages 2–3).** The signature below is the
*eventual* dqrng/long-output kernel. Stage 0 only extracts the seam and can keep a
base-R, wide-feeding shape internally; the long-output signature, `base_seed`, and
`store_crossovers` arrive with the dqrng reference kernel (Stage 2). It is shown here
because the whole design (multi-map grouping, memory layout, canonical order below)
targets this interface:

```r
make_gametes_batch(
  parent_allele,          # per-parent 2 x n_autosome_loci integer (homolog x locus)
  parent_line_origin_code,# per-parent 2 x n_autosome_loci integer line_origin codes
  parent_1_idx,           # integer, length n_offspring: parent_1 (sire) index per offspring
  parent_2_idx,           # integer, length n_offspring: parent_2 (dam) index per offspring
  chr_start, chr_end,     # integer per chr: locus-index range into 1..n_autosome_loci
  chr_pos_M, chr_len_M,   # double per chr: locus genetic positions + length, in MORGANS (pos_cM/100)
  base_seed,              # integer: per-gamete streams derived from (base_seed, o, r)
  store_crossovers = FALSE# logical: emit the ragged crossover buffer (decision #9)
) -> list(
  # haplotype output (always):
  parent_origin    = <int, length 2 * n_offspring * n_autosome_loci>,  # 1=parent_1, 2=parent_2
  locus_idx        = <int, ...>,   # 1..n_autosome_loci; caller maps to locus_id/locus_name
  allele           = <int, ...>,
  line_origin_code = <int, ...>,   # dictionary code; caller decodes to line_origin VARCHAR
  # crossover output (only when store_crossovers = TRUE; ragged, ~0.06% of the above):
  xover_gamete_o   = <int, length total_crossovers>,  # offspring index o per event
  xover_parent_origin = <int, ...>,                   # 1=parent_1, 2=parent_2
  xover_chr        = <int, ...>,                       # chromosome (autosomes here)
  xover_pos_cM     = <double, ...>                     # crossover genetic location in cM (Morgans*100)
)
```

The kernel is **Morgans-native** (decision #10): the caller resolves each gamete's
map for its `(parent sex, line)` from `genome_map`, converts `pos_cM → Morgans`,
and passes `chr_pos_M`/`chr_len_M`; the kernel never sees bp/Mb. Crossover positions
come back in **cM** for `ind_crossover`.

### Multiple resolved maps — group-by-map calls (critique #1)

A gamete's map is its **producing parent's** `(sex, line)` map, so a sire gamete and
a dam gamete — or two parents on different lines — can carry **different** maps
(sex-specific recombination, achiasmy, line-specific maps). The signature above has a
**single** `chr_pos_M`/`chr_len_M`, i.e. it represents **one** map. The chosen design
resolves that:

> **Decision (locks critique #1): group gametes by resolved map key and call the
> kernel once per `(sex, line, map_name)` group.** The caller already caches a
> resolved descriptor set per unique `(parent_sex, parent_line)` (the coordinate
> migration built this — `R/add_offspring.R`). Stage 2+ groups the batch's gametes by
> their resolved-map key and issues one `make_gametes_batch()` call per group, each
> with that group's `chr_pos_M`/`chr_len_M`. The per-gamete `dqrng` stream key
> (below) carries the **global** offspring index `o`, so splitting the batch across
> map-groups does **not** change any gamete's random draws — output is independent of
> how gametes are partitioned into calls.

Rejected alternatives, and why: a **map-indexed batch** (pass `gamete_map_idx` +
packed per-map descriptors in one call) keeps a single call but complicates the C++
input layout for no throughput gain at realistic map-group counts (≤ a few); a
**v1-only shared-map kernel** was rejected outright — it would silently regress the
sex/line map seam the migration just established. In v1 all gametes still resolve to
the single `(NULL, NULL)` default map, so there is exactly one group and one call —
but the grouping is real, exercised, and future-proof.

### Kernel input memory layout (critique #6)

`add_offspring()` has many unique parents; the kernel must receive them in a layout
that is cheap in R and tight in C++. **Decision:**

- **Packed parent array.** Build one contiguous integer buffer holding, for each
  unique parent in a fixed order, its `2 × n_autosome_loci` allele block (homolog-major:
  homolog 1 row then homolog 2 row), plus a parallel `line_origin_code` buffer of the
  same shape. In R this is a single `matrix`/array assembled once from the `split()`
  parent structure (Stage 0); in C++ it is a flat `int*` with stride `n_autosome_loci`.
- **Parent index vectors.** `parent_1_idx` / `parent_2_idx` (length `n_offspring`) are
  **0-based (C++) / 1-based (R) row indices into the packed parent list**, not
  `id_ind`. The caller owns the `id_ind → packed-index` map.
- **Special-chromosome-only parents** (a parent that contributes only to the R special
  path — e.g. supplies Y but no autosome variation) still appear in the packed autosome
  list with their real autosome haplotypes; nothing special is needed because every
  diploid parent has autosomes. Parents are packed once per `add_offspring()` call and
  shared read-only across all offspring and threads (see Parallelization).

This layout lets the R reference and the C++ kernel take **byte-identical** inputs, so
the parity test compares like with like.

### Canonical output order (critique #16)

The long output vectors have a **fixed, documented order**, independent of execution
order:

```text
offspring in original `matings` row order
  -> parent_origin 1 (parent_1/sire) then parent_origin 2 (parent_2/dam)
    -> chromosome ascending
      -> locus ascending (locus_idx)
```

**Execution order may differ from output order.** Parent-grouping for cache locality
(Parallelization) and group-by-map calls (above) reorder *computation*, but the caller
writes rows in the canonical order above, and every gamete's `dqrng` stream is keyed on
its **global** `o` (original `matings` row), so neither grouping changes output bytes or
the stream index. Tests assert this order so DB inserts, batching-invariance, and R↔C++
parity all compare deterministically.

When `store_crossovers = FALSE` the four `xover_*` vectors are length-0 (or
absent) and nothing is allocated for them.

**Scope of the kernel:** autosomes only. `chr_start`/`chr_end` assume each
chromosome's loci are **contiguous** in the compacted autosome ordering. This
holds today by construction — `define_genome()` assigns loci in contiguous
blocks (`R/define_genome.R:133-134`: `rep(seq_len(n_chr), times = loci_per_chr)`
with `locus_id = seq_len(n_loci)` at `:173`), and dropping whole special-chromosome blocks
preserves contiguity of the rest. **Caveat:** any loci added after founding
(Stage-5 novel mutations / gene edits) may append out-of-order `locus_id`s and
break this, so the kernel should either (a) assert contiguity on entry, or (b)
accept an explicit per-locus `chr_of_locus` vector instead of start/end ranges. Prefer (b)
if it costs nothing — it is future-proof. Special chromosomes (sex chromosomes /
organelles) are **not** in the kernel in v1 — they stay on the R path (low locus
count, "correctness first") and are appended to the same long output by the
caller.

Also in Stage 0 (cheap, measured against Phase 0):

- Replace the O(P×N) parent scan (the `parent_haps_raw[... == pid, ]` subset at
  `R/add_offspring.R:321`, inside the `unique_parents` loop `:320-405`) with a
  single `split(parent_haps_raw, parent_haps_raw$id_ind)`.
- Preallocate the `special_rows` structure instead of `data.frame`-per-iteration.

**Stage 0 is a pure structural refactor** — it keeps today's `make_gamete()`
internals, base-R RNG, **and the current wide write path**, so seeded `ind_haplotype`
output is unchanged *by construction* (seam extraction, `split()`, and preallocation
are RNG-neutral; nothing about the write shape changes yet). This is a
**refactor-safety property** (Stage-0 output equals the immediately-prior code on this
branch), **not** a cross-version compatibility contract — per CLAUDE.md's pre-1.0.0
stance we don't preserve or test against any prior version's output. A transient
before/after diff during development is enough; **no committed golden-vs-old fixture is
required**. The write path becomes long in **Stage 1** (values still base-R, so still
neutral); the intentional RNG-stream change lands in **Stage 2** (dqrng). From then on
the reproducibility contract is **R↔C++ parity + forward same-seed determinism +
distributional checks**, never a match to base R's old stream.

---

## Stage 1 — Direct long write + batching (base-R kernel, shared with `add_founders()`)

Convert the write path to long format and add batching **before** any dqrng long-output
kernel is wired into `add_offspring()` (critique #2, #3). Allele values still come from
base-R `make_gamete()` via the Stage-0 seam, so only the write *shape* changes — this
stage remains output-neutral for `ind_haplotype`. This ordering pins the memory ceiling
first, so the dqrng kernel (Stage 2) drops into an already-long, already-batched path.

### 1a. Direct long insert (no wide, no UNPIVOT)

Shared internal helper, e.g. `.write_long_haplotypes(conn, long_frame,
line_origin_decode)`, used by **both** `add_offspring()` and `add_founders()`
(this is the shared "allele source" seam from "The founder path" section — the
write core is identical, only the source differs):

- **Seam shape in Stage 1 (important — the seam is still WIDE here).** The Stage-0
  `make_gametes_batch()` returns the four **wide** matrices
  (`sire_allele`/`dam_allele`/`sire_lo`/`dam_lo`, each `B × n_autosome_loci`); the
  long-native seam signature is a **Stage 2** change (see Stage 0 / the target
  architecture). So in Stage 1 the caller **reshapes each batch's wide matrices to
  long parallel vectors in R** — a cheap `as.vector()`/`rep()` reshape **bounded by
  the batch `B`**, *not* a dense full-`n_offspring` frame and *not* a SQL `UNPIVOT` —
  then hands `(parent_origin, locus_idx/locus_id, allele, line_origin)` to the helper.
  Stage 2 makes the seam emit those long vectors directly and this reshape disappears.
- The helper registers the long parallel vectors as a temp Arrow/DuckDB relation.
  Caller has already attached `id_ind`, `strand = 1L`, and either resolved
  `locus_idx → locus_id` or left the join to the helper.
- `INSERT INTO ind_haplotype (id_ind, parent_origin, strand, line_origin,
  locus_id, locus_name, allele) SELECT ... JOIN genome_meta` — a plain join
  resolving `locus_id`/`locus_name` (and, once the kernel emits codes, decoding
  `line_origin_code → line_origin`), **no `UNPIVOT`** (the current UNPIVOT+INSERT at
  `R/add_offspring.R:~594-620` — post-Stage-0 line numbers — and
  `R/add_founders.R:279-308` both go away). `strand` is always `1L` in this
  version (ploidy 2).
- Convert `add_founders()` to build its rows in long form too (it currently
  builds the wide `hap_matrix` + `hap_wide` at `R/add_founders.R:244-257`);
  founders have no meiosis and just materialize **sampled whole-haplotype** rows
  in long form (the founder allele source). At the large-scale target
  `add_founders()` also reconstructs the dense pool `hap_data_matrix`
  (`R/add_founders.R:178-180`) and the selected `(2 · n_founders) × n_loci` matrix
  (`R/add_founders.R:244-246`) — **two** dense matrices — so it must use the
  **same batching** (1b) to bound memory, not only `add_offspring()`. Materialize
  each batch's long rows **directly from the pre-drawn sample indices against the
  pool** (see RNG trap), skipping both dense matrices entirely.

  **Founder RNG trap (critique #8).** `add_founders()` currently draws **all**
  haplotype indices in **one** `sample()` call (`R/add_founders.R:200-205`). If founder
  batching sampled per batch, the base-R RNG output would change. To keep founder
  long-write **equivalence** (test 8), the refactor must **pre-draw all haplotype
  indices in the single existing `sample()` call up front**, then batch only the
  *materialization + write* of the selected haplotypes — never the sampling. (Batching
  the write, not the draw, is what makes founder output byte-identical to the pre-refactor
  code for a fixed seed.) At target scale, also avoid building the whole dense selected
  matrix before writing — materialize each batch's long rows directly from the
  pre-drawn indices, so peak memory is bounded by `B`, not `n_founders`.

  **`define_founder_haplotypes()` pool write (the third dense materialization).**
  `.write_founder_haplotypes()` (`R/founder_haplotype_helpers.R:206-225`) currently
  builds the **entire** long `founder_haplotypes` frame (`n_hap × n_loci` rows) in R
  via `rep()` before `dbWriteTable`. At target pool scale this is a 250M-row R
  data.frame on top of the wide `haplotype_matrix`. Stage 1 must **batch the pool
  write** the same way — stream the wide `haplotype_matrix` to the long
  `founder_haplotypes` table in column/haplotype batches (or via a DuckDB
  `register` + `INSERT ... SELECT` over a wide temp relation, analogous to the
  `ind_haplotype` path but keyed on `haplotype_id`/`locus_name`) so peak memory is
  bounded, not `O(n_hap × n_loci)`. `founder_haplotypes` has its **own** schema
  (`line_name, haplotype_id, locus_name, allele`), so it does **not** reuse
  `.write_long_haplotypes()` verbatim, but it uses the same batched-register-and-insert
  mechanism. RNG note: the pool generator draws are unchanged — only the write is
  batched — so seeded pools stay identical (same discipline as the founder RNG trap).

### 1b. Batching for memory

Process `matings` in batches of `B` offspring (generate → write → free → `gc()`).
Bounds peak memory at `~B × n_loci` long rows **regardless of total `n_offspring`**.
**Batch index must preserve the original `matings` row order (and the global
offspring index `o`)** so that (i) in Stage 1 the base-R draw order is unchanged by
batching, and (ii) from Stage 2 on, gamete `(o, r)` uses stream `(o, r)` regardless of
which batch it lands in. Either way any `B` (including `B = n_offspring`) yields
identical output — no ordering constraint at batch boundaries.

**`B` is exposed AND RAM-aware by default (resolves Open Q#3).** Decided:

- **`add_offspring(..., batch_size = NULL, max_batch_mem = NULL)`.** `batch_size`
  is an explicit offspring-per-batch override; `max_batch_mem` is an explicit
  per-batch memory budget (bytes/`"512MB"`). Either wins over the default when set.
- **Default (both `NULL`): auto-pick `B` to not overrun RAM.** Detect *available*
  system memory, target a conservative fraction of it (e.g. ~25%), and set
  `B = target_bytes / (n_loci × bytes_per_offspring_row)` where the per-row constant
  is calibrated from the Phase-0 numbers. This protects the 8–16 GB laptop
  automatically and lets a 256 GB HPC node use large batches without configuration.
- **Cross-OS memory detection (build correct once — main platforms only).** Support
  **Windows, macOS (Intel + Apple Silicon), and mainstream Ubuntu/Linux**; do **not**
  chase exotic OSes. Prefer a small, portable path: parse `/proc/meminfo`
  (`MemAvailable`) on Linux, `sysctl hw.memsize` + `vm_stat` on macOS,
  `GlobalMemoryStatusEx`/`wmic OS get FreePhysicalMemory` on Windows — or lean on a
  maintained helper (`memuse`/`benchmarkme`) if it covers all three cleanly. **If
  detection fails, fall back to a fixed conservative default** (target ~512 MB/batch)
  — never error, never assume unlimited RAM.
- **Speed/RAM trade-off:** smaller `B` → lower peak, slightly more per-batch overhead
  (transaction, register/unregister, INSERT round-trips); the penalty is mild over a
  broad range and only bites at very small `B`. The long-write itself (1a) and integer
  `line_origin` codes are win-win (lower memory *and* faster); `B` is the only genuine
  dial.

### 1b-note. High-fecundity / low-memory contract (aquaculture, poultry)

The binding RAM quantity is **`n_offspring_in_one_call × n_loci`**, not the whole
database (which lives on disk in DuckDB). Because offspring are independent given
their parents, batching operates on **offspring rows** and bounds peak memory the
same way whether those offspring come from many pairs or from **one** mega-mating
(e.g. an aquaculture spawn of hundreds of thousands, or a poultry line). **So a single
huge mating is handled by batching alone — no separate mechanism is needed**, and a
16 GB machine can generate an arbitrarily large single-mating cohort by looping
batches, streaming each to disk. The parent side is loaded once and is tiny (2 parents
× 2 haplotypes × `n_loci`), so it does not scale with offspring.

Reference bounds (order-of-magnitude, anchored on the Phase-0 measurement of 418 MB
peak for 500 × 5k; **confirm with the gated large run**):

- **RAM, current wide path (pre-Stage-1):** peak ≈ `~50–150 B × n_offspring × n_loci`.
  On 16 GB the worry line is `product ≳ 5×10⁷` (e.g. 1,000 offspring × 50k loci).
- **RAM, post-Stage-1 (batched, long, integer codes):** peak ≈ `~64 B × B × n_loci`,
  independent of total offspring. 8 GB → `B ≈ 600` @ 50k loci; 16 GB → `B ≈ 1,250`;
  256 GB → `B ≈ 20k` or unbatched.
- **Disk (DuckDB, cumulative):** `ind_haplotype` is 2 × `n_loci` rows/individual,
  heavily compressed (alleles bit-pack; `id_ind`/`line_origin` dictionary+RLE;
  `locus_id` delta) to ~1–4 B/row, so ~400 GB usable holds **~1–2 M individuals at
  50k loci**. Disk is not the binding constraint. (A later schema question: dropping
  the denormalized `locus_name` from `ind_haplotype` in favor of a join would cut disk
  further — out of scope here.)

**Deferred, deliberately (non-foreclosing): a count-based `matings` spec.** For an
extreme single mating the *user-supplied* `matings` tibble is itself `O(n_offspring)`.
An optional `n_offspring` count column (one row expands to N offspring internally)
would avoid materializing millions of identical rows. **Not built now** — it adds a
second input convention and per-row extra-column semantics we don't want to confuse
users with yet. It is **purely additive, needs no schema/DB change, and stays
non-breaking whenever added.** Stage 1 must not foreclose it: keep `matings`
normalization/expansion at the **top** of `add_offspring()` (before the batch loop) so
a future count column feeds the same batched path with no refactor.

### `line_origin` encoding (kernel boundary, Stage 2+)

`line_origin` is `VARCHAR` but has bounded cardinality (number of founding
lines). Once the kernel boundary exists (Stage 2+), pass it as **integer codes**
(dictionary built in R from the parents' `line_origin` values), gather codes in C++,
and decode back to strings (or a DuckDB `ENUM`) in the write helper — avoids moving
strings through the kernel and repeating them per allele row. Stage 1's base-R path
carries `line_origin` strings directly (no kernel boundary yet).

### 1c. `ind_crossover` schema (created here; rows written from Stage 2)

Create the empty `ind_crossover` table in `define_genome()` and register it across the
schema tables **now** (critique #9 — full list in the `ind_crossover` section above), so
`restore_pop()` sees it and it exists before the first insert. The **row writes** wait
for the dqrng kernel's `xover_*` buffer (Stage 2) — Stage 1's base-R path does not emit
crossover events. Once wired, the opt-in write registers the batch's `xover_*` buffer
and `INSERT`s to `ind_crossover` in the same batched loop, assigning `id_crossover` via
`next_int_id()` and joining `genome_meta` only to attach `chr_name`. Off ⇒ no
`ind_crossover` writes (the empty table still exists, holding zero rows).

### 1d. Transaction scope (critique #10)

Wrap each `add_offspring()` call in **one transaction** (`BEGIN … COMMIT`) covering all
its writes — `ind_meta`, every batch's `ind_haplotype` insert, and (when opted in)
`ind_crossover`. With batching, a mid-run failure could otherwise leave a **partial
generation** (some batches written, some not). One transaction per call makes the whole
generation atomic — it lands completely or rolls back cleanly, and `restore_pop()` never
sees a half-written generation. (Per-batch transactions were considered and rejected:
they reintroduce the partial-generation hazard for negligible gain, since the batch loop
already bounds memory.) `add_founders()` gets the same one-call transaction wrapper,
and `define_founder_haplotypes()` wraps its (now batched) pool write in one
transaction too.

**Deliverable:** shared `.write_long_haplotypes()` helper; **all three** write paths
on the long, batched path — `add_offspring()`, `add_founders()` (with the pre-draw
RNG discipline), and `define_founder_haplotypes()`'s pool write (no full-frame
materialization); one transaction per call each; `ind_crossover` created + registered
in schema; before/after write-time + memory benchmark vs the Phase 0 wide baseline for
**both** the offspring and founder paths; a documented decision on Appender API vs
`register`+`INSERT SELECT` from a direct A/B.

---

## Stage 2 — Pure-R reference generator (`dqrng`, long output)

Implement the canonical per-gamete-stream algorithm in R inside
`make_gametes_batch()`, using **`dqrng`** (not base `runif`/`rpois`):

- Seed a `dqrng` sub-stream per gamete `(o, r)` from the base seed.
- Within each gamete's stream, draw the crossover count via the
  **inversion/Knuth sampler on `dqrng` uniforms** (decision #8), then the
  `n_cross` position uniforms, in the fixed per-chromosome order.
- Do the deterministic `findInterval` + parity + gather **vectorized across the
  batch** where possible; emit long parallel vectors directly.
- When `store_crossovers = TRUE`, also emit the ragged `xover_*` buffer (the cM
  positions already in hand). Gated so the default path allocates nothing.

This is the **reference implementation** — the C++ parity oracle and the
no-compiler fallback. The inversion Poisson sampler here is the **exact
algorithm the C++ kernel must mirror** (do not use base `rpois`), so it is the
authoritative spec for the sampler.

**Deliverable:** `make_gametes_batch()` (R, `dqrng`-based, inversion Poisson,
optional crossover buffer) wired into `add_offspring()`; Phase 0 benchmark
re-run.

---

## Stage 3 — C++17 kernel (Rcpp + `dqrng`) + parity test

Port `make_gametes_batch()` to `src/make_gametes.cpp`, following the canonical
per-gamete order exactly, driving randomness from the **`dqrng` engine via
`dqrng.h`** (`LinkingTo: dqrng, BH, sitmo`) — the same generator family the R
reference uses, so streams match. The Poisson count uses the **same
log-accumulation inversion sampler** as the Stage-2 R reference (decision #8) —
**not** `std::poisson_distribution` or any `<random>`/`rand()` source — so parity
is guaranteed by the shared uniform stream. Emit the ragged `xover_*` buffer when
`store_crossovers = TRUE`.

**Prerequisite:** the dqrng R↔C++ uniform-parity spike (see the stream-derivation
section) must pass **before** this stage — the kernel is built on the proven scheme,
not around an unverified one.

- Work on raw `int*` / `Rcpp::IntegerVector`; no per-locus allocation, no
  `cbind` index matrices; write straight into preallocated long output vectors.
- Packaging: `LinkingTo: Rcpp, dqrng, BH, sitmo`, `Imports: Rcpp, dqrng` in
  `DESCRIPTION`; `useDynLib(tidybreed, .registration = TRUE)` +
  `importFrom(Rcpp, sourceCpp)` in `NAMESPACE`; run
  `Rcpp::compileAttributes()`; add `src/Makevars`/`Makevars.win` if needed.
- **Compiler requirement (critique #15).** Once the package ships C++ sources +
  `LinkingTo`, a **source install requires a working compiler** — this is normal for
  a compiled R package and is the accepted state after this stage. The R reference is
  **not** a no-compiler install path; it is kept as a **runtime/debug fallback and the
  parity oracle** (an internal selector prefers the compiled kernel and falls back to R
  if the shared object is unavailable at runtime). Binary installs (CRAN/r-universe)
  need no compiler as usual. We are **not** building an optional-compile `configure`
  path — that complexity is not worth it now.

**Parity test (the core deliverable):**

```r
base_seed <- S
r_out   <- make_gametes_batch_R(args,   seed = base_seed, store_crossovers = TRUE)
cpp_out <- make_gametes_batch_cpp(args, seed = base_seed, store_crossovers = TRUE)
# compares haplotype output: parent_origin, locus_idx, allele, line_origin_code
# AND crossover output:      xover_gamete_o, xover_parent_origin, xover_chr, xover_pos_cM
testthat::expect_identical(r_out, cpp_out)
```

Run across several base seeds, several `n_chr`/`chr_len` configs, a
special-chromosome config, **and single- vs multi-threaded** (must be identical —
that is the point of per-gamete streams). Testing with `store_crossovers = TRUE`
also validates the inversion Poisson sampler R↔C++ (the crossover positions
expose any count mismatch immediately). This satisfies the R↔C++ RNG contract
(decisions #2, #3).

**Deliverable:** C++ kernel; parity test green; `R CMD check` passes with and
without the compiled code; Phase 0 benchmark quantifies the win (expect a large
multiple over Stage 2 on the kernel itself).

*(The long-format direct write + batching + `ind_crossover` schema — formerly Stage 3
— now lands in **Stage 1**, before any long-output kernel; the opt-in crossover **row
write** is wired here once the kernel emits the `xover_*` buffer.)*

---

## Testing strategy

0. **Stage 0/1 refactor-safety (transient, within-branch)** — fixed seed,
   `ind_haplotype` contents identical to the immediately-prior code on the branch
   through **both** the Stage-0 seam extraction and the Stage-1 wide→long write
   conversion (both keep base-R `make_gamete()`, so allele values are unchanged; only
   the write shape moves). This is a development-time before/after diff, **not** a
   committed golden-vs-old fixture, and it is dropped the moment **Stage 2**
   intentionally changes the stream (pre-1.0.0 — no cross-version contract).
1. **R↔C++ parity** (Stage 3) — identical output for a fixed base seed, multiple
   configs, **including at least one configured sex chromosome** (exercises the
   wrapper-level parity described in Special chromosomes — C++ autosomes + R special).
   The core RNG contract.
2. **Thread-count invariance** (Stage 4) — same base seed, 1 core vs N cores →
   identical `ind_haplotype`. The payoff of per-gamete streams.
3. **Batching invariance** — same base seed, `B = n_offspring` vs `B` small →
   identical `ind_haplotype` (from Stage 1 for base-R order-preservation, and from
   Stage 2 for full stream-keyed independence).
4. **Correctness / distribution** — recombination-fraction sanity: observed
   between-adjacent-loci switch rate ≈ Haldane `0.5*(1-exp(-2d/100))` over many
   offspring; allele frequencies preserved from parents in expectation.
5. **Poisson inversion sampler** (Stage 2) — direct unit tests of the
   log-accumulation sampler at `lambda = 0` (always 0 crossovers, no RNG consumed
   beyond the guard), small (`~0.3`), typical (`~1`), and a **stress `lambda` near the
   documented ceiling** (confirms no underflow and `O(lambda)` termination). Mean and
   variance ≈ `lambda`; R and C++ counts identical for the same stream.
6. **Crossover storage** (`store_crossovers = TRUE`) — (a) crossover **count**
   per gamete per chromosome ≈ Poisson(genetic length in Morgans) in mean and
   variance; (b) crossover **positions** ≈ uniform on `(0, len_cM)`; (c)
   **consistency, with a controlled fixture (critique #12):** use a fixture where the
   two parental homologs are encoded **differently at every marker** (or expose an
   internal homolog trace from the generator), so every homolog switch is visible in
   `ind_haplotype` — then assert the recorded `ind_crossover` rows reproduce exactly
   those switch points. Do **not** rely on random founder alleles, under which a switch
   between identical-allele/identical-`line_origin` homologs is invisible and the test
   would silently under-check.
7. **`store_crossovers` has no side effects** — `ind_haplotype` contents are
   byte-identical with the flag off vs on for the same seed (storing crossovers
   must not perturb the allele stream).
8. **Founder-path long-write equivalence** — for a fixed seed: (a)
   `define_founder_haplotypes()`'s `founder_haplotypes` pool is byte-identical after
   batching the pool write (only the write shape changes, generator draws unchanged);
   and (b) `add_founders()`'s `ind_haplotype` allele values are unchanged after moving
   to the long/batched path — both require the founder-RNG pre-draw discipline (Stage 1
   founder note). Also assert the founder path builds **no** dense `n × n_loci` R
   matrix (peak memory bounded by `B`), the point of the refactor.
9. **Transaction atomicity** (Stage 1) — an injected mid-generation failure leaves
   `ind_meta`/`ind_haplotype`/`ind_crossover` with **no** rows from the aborted call
   (rollback), never a partial generation.
10. **Existing suite** — full `tests/testthat/` green; `R CMD check`.

---

## Explicitly out of scope

- **Wide matrices / `UNPIVOT`** in the write path — removed by decision. The only
  wide representation is the user-facing `extract_genotypes()`.
- **C++ reading/writing DuckDB tables directly** (linking `duckdb` into the package
  and using the C++ `Appender` / reaching into the R package's connection pointer) —
  rejected for ABI fragility and CRAN-hostility, for a marginal gain the batched
  zero-copy `register`+`INSERT` design already neutralizes. See "Data path & memory
  model". The C++ kernel is a pure in-memory numeric function; R owns all DB I/O.
- **CUDA/GPU** — meiosis simulation is **branchy** (interval search per locus) and
  **memory-bound** (scatter writes into long-format rows), not the dense
  linear-algebra shape GPUs accelerate. It would also fragment the dev/test story:
  CLAUDE.md targets both a Mac and a Windows AVD, and **Apple Silicon has no CUDA
  path at all**. Rejected entirely so it isn't re-litigated. *(Rationale salvaged
  from the retired `plans/speed_up_offspring.md`.)*
- **Real polyploidy** — tracked in `plans/refactor_haplotype.md`.
- **Schema changes** to `ind_haplotype` / `ind_meta`.

*(Multithreading is now **in scope** — Stage 4 — because the per-gamete `dqrng`
stream model makes it reproducible. It is no longer a deferred stretch.)*

---

## Open questions for review

- ~~**Poisson-parity resolution**~~ — **RESOLVED** (decision #8): Poisson via a
  `dqrng`-uniform inversion/Knuth sampler, keeping the position-based model. The
  crossover-storage feature (decision #9) drove this — it needs exact genetic
  (cM) positions, which the uniform-only Markov model cannot produce.

1. ~~**Per-gamete stream seeding scheme**~~ — **RESOLVED by the spike** (2026-07-07,
   `dev/spikes/dqrng_parity/`): fold the tuple to `sid = ((o*2 + (r-1))*2 + (kind-1))`
   and seed with the **2-int seed `c(base_seed, sid)` through dqrng's own
   `convert_seed`** (Scheme A) — no hand-rolled hash. Proven bit-identical R↔C++ (both
   Scheme A and the seed+stream/jump Scheme B). See the stream-derivation section.
2. **Base-seed storage / audit** (critique #4) — when `store_crossovers = TRUE` (or
   always), persist the resolved `base_seed` so a run is exactly replayable?
   Recommendation: yes, a single cheap `base_seed` scalar/column with the crossover
   rows; decide at Stage 2. Default off otherwise (CLAUDE.md's "disdain for metadata").
3. ~~**Batch size `B`**~~ — **RESOLVED (2026-07-07):** expose `batch_size` and
   `max_batch_mem` args **and** default to a RAM-aware auto-pick (detect available
   memory on Windows/macOS-Intel+ARM/Ubuntu, target a conservative fraction, fixed
   fallback if detection fails). Batching bounds peak independent of total offspring,
   which also covers a single high-fecundity mega-mating (aqua/poultry) with no
   separate mechanism. See Stage 1 (1b) + the high-fecundity/low-memory contract.
4. **`line_origin` as DuckDB `ENUM`** vs decoded `VARCHAR` on insert — decide
   from the Stage 1 write benchmark.
5. **Appender API vs register+INSERT** — decide from the Stage 1 A/B, not a priori.
6. **Keep the R reference path permanently** as a fallback, or drop it once the
   kernel is stable? (Pre-v1, low cost to keep; it is also the parity oracle.)

---

## Implementation order

| Stage | What | New dep? | Status |
|---|---|---|---|
| **Pre** | **`plans/position_coordinates.md`** — `pos_Mb` → `pos_bp` + `genome_map`/`pos_cM` migration | None | ✅ **Implemented + green** (working tree, `feat/position-coordinates-genome-map`; 0 fail / 1236 pass) |
| **Spike** | dqrng R↔C++ uniform-parity check for `(base_seed, o, r, stream_kind)` (blocks Stage 3) | dqrng | ✅ **Done + green** (`dev/spikes/dqrng_parity/`; Scheme A adopted — 2-int seed `c(base_seed, sid)` via `convert_seed`, `xoroshiro128plusplus.uniform01()`) |
| 0 (profile) | Extend harness: per-step timing + memory (skip gracefully if `profvis`/`bench` absent), realistic mating designs; committed baseline table | None | ✅ **Done + green** (`dev/benchmarks/profile_haplotype_paths.R`; baseline + conclusion committed) |
| 0 (seam) | Extract `make_gametes_batch()` seam wrapping base-R `make_gamete()`, **wide write kept**; `split()` fix (index-split); preallocate `special_rows`; output-neutrality check | None | ✅ **Done + green** (parity golden + seam-equivalence; 1229 pass / 0 fail; O(P×N) ~12× faster, churn ~28× lower) |
| 1 | **Direct long write + batching** across **all three** write paths — `add_offspring()`, `add_founders()` (founder pre-draw discipline, kills both dense matrices), **and `define_founder_haplotypes()`'s pool write** (batched, no full-frame materialization) — `line_origin` path, `ind_crossover` schema + registration, **one transaction per call**, **exposed `batch_size`/`max_batch_mem` with a RAM-aware auto-default (cross-OS: Win/macOS-Intel+ARM/Ubuntu, conservative fallback)**, `matings` expansion kept at the top (non-foreclosing count-spec) — still base-R `make_gamete()` internally | None | ✅ Done (2026-07-07) — output-neutral (byte-identical golden, autosome + sex-chr, multi-batch); batch invariance + txn atomicity + 1257 tests green; `add_offspring(200)` mem_alloc 1.55GB→464MB, peak decoupled from offspring count. See [Stage-1 summary](optimize_add_offspring_stage_1_summary.md). |
| 2 | Pure-R reference generator, per-gamete `dqrng` streams, **log-accumulation** inversion Poisson sampler, long output, optional crossover buffer, seed API (intentional stream change) | dqrng | ⬜ Not started |
| 3 | C++17 kernel (Rcpp + `dqrng.h`), one entry point, autosome-only, same log-accumulation Poisson, group-by-map calls, R↔C++ parity test (incl. crossovers + sex-chr wrapper parity) | Rcpp, dqrng, BH, sitmo | ⬜ Not started |
| 4 | Parallelize across individuals (**`RcppParallel::parallelFor`** over offspring — **not** OpenMP, for macOS install portability; `std::thread` as zero-dep fallback); thread-count invariance test | RcppParallel (header-only TBB, no system lib) | ⬜ In scope (after 3) |

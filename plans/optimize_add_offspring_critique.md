# Critique: `plans/optimize_add_offspring.md`

Review date: 2026-07-06

Scope: plan review for the next implementation stage. I reviewed the current
plan against the current code, especially `add_offspring()`, `add_founders()`,
`recombination_helpers.R`, `genome_map_helpers.R`, `DESCRIPTION`, and the
benchmark harness. I did not change implementation code.

## Overall Assessment

The plan has the right high-level direction:

- move the hot gamete work out of the per-offspring R loop;
- use genetic map positions in Morgans, not physical bp/Mb;
- avoid wide matrices and `UNPIVOT` in the hot write path;
- batch to bound memory;
- use deterministic per-gamete streams so batching and parallelism can be
  reproducible;
- keep a pure-R reference and require R-to-C++ parity.

The coordinate prerequisite now appears to satisfy the plan's key assumption:
current `add_offspring()` resolves and caches a descriptor set per unique
`(parent_sex, parent_line)`, and sire/dam/special-chromosome gametes use the
gamete-producing parent's map (`R/add_offspring.R:222-292`,
`R/add_offspring.R:432-445`, `R/add_offspring.R:476-478`).

However, the optimization plan is **not implementation-ready as written**. The
main problem is that the proposed kernel interface and stage order do not yet
match the plan's own requirements: per-gamete maps, long output, batching,
special-chromosome RNG, and a C++ fallback/build story need to be made precise
before implementation starts.

## Highest-Priority Concerns

### 1. The proposed kernel signature cannot represent per-parent/sex/line maps

The plan correctly states that the caller resolves the map per gamete and that
sire and dam gametes may use different maps (`plans/optimize_add_offspring.md:
254-265`). But the proposed `make_gametes_batch()` signature has only one
`chr_pos_M` / `chr_len_M` input (`plans/optimize_add_offspring.md:489-497`).
That signature can represent a single shared map, not a batch where parent_1 and
parent_2, or different parent lines, use different map descriptors.

This is the most important interface gap. The current R implementation now has
the right cache structure, so the optimized kernel must preserve it.

Recommended fix: choose one of these designs and write it into the plan:

- **Group-by-map calls:** split gametes by resolved map key and call the kernel
  once per `(sex, line, map_name)` group. Simpler input shape, more calls.
- **Map-indexed batch:** pass `gamete_map_idx` for each gamete plus packed map
  descriptors (`map_chr_start`, `map_chr_end`, `map_chr_pos_M`, `map_chr_len_M`)
  keyed by map id. More complex, but keeps one batch call.
- **Temporary v1 shared-map scope:** explicitly state that the first optimized
  kernel only supports the default map and errors if non-default map rows exist.
  This is the least future-proof option and conflicts with the stated goal, so I
  would avoid it.

Without this change, the implementation will quietly regress the sex/line map
seam that the coordinate migration just established.

### 2. Stage ordering is internally inconsistent: long-output kernel arrives before long-write/batching

Decision #1 says no wide matrices anywhere in the write path
(`plans/optimize_add_offspring.md:20-31`). The target architecture has the
kernel emit long vectors directly (`plans/optimize_add_offspring.md:205-218`).
But the implementation order puts direct long insert and batching in Stage 3
(`plans/optimize_add_offspring.md:628-686`), after Stage 0/1/2 have already
introduced a long-output kernel.

That leaves an implementation gap: if Stage 0/1 emits long rows but the current
write path still expects wide matrices and `UNPIVOT`, the code either has to
re-widen the long output or partially implement Stage 3 early. Re-widening
violates the plan's core decision and would hide performance/memory behavior.

Recommended fix: reorder the stages:

1. Phase 0: profile.
2. Stage 0: extract current base-R gamete logic behind a seam but keep the
   current wide write path, **or** implement the direct long write helper first.
3. Stage 1: direct long write + batching, still using the current base-R
   `make_gamete()` internally.
4. Stage 2: dqrng R reference kernel.
5. Stage 3: C++ kernel.
6. Stage 4: parallelism.

The safer path is to implement direct long write + batching before any
long-output kernel is wired into `add_offspring()`. That gives a stable memory
ceiling before allocating long vectors for large batches.

### 3. Batching must move earlier, before any full long-output kernel is wired into production

The scale target is 10k+ offspring and 50k+ loci (`plans/optimize_add_offspring.md:
59-61`). A full-generation long output is about:

```text
2 * n_offspring * n_loci rows
```

At 10k x 50k, that is 1 billion haplotype rows before accounting for vectors,
line-origin codes, crossovers, or DB staging. The plan currently places batching
in Stage 3 (`plans/optimize_add_offspring.md:654-662`), after the Stage 1
long-output R reference is "wired into `add_offspring()`"
(`plans/optimize_add_offspring.md:577-579`).

That is too late. A no-wide long-output design only works at target scale if
batching is part of the first production path that emits long vectors.

Recommended fix: make batching a prerequisite for Stage 1 and Stage 2, not a
later write-path improvement. The batch index must preserve original mating row
index `o` for per-gamete streams so batch size and ordering do not affect output.

### 4. The RNG API is not specified: how does the user provide `base_seed`?

The plan's reproducibility contract depends on a `base_seed`
(`plans/optimize_add_offspring.md:274-275`, `plans/optimize_add_offspring.md:
496`). Current `add_offspring()` uses the global base R RNG via `set.seed()`;
there is no `seed` or `base_seed` argument.

The plan needs to decide the user-facing API before implementation:

- Add `seed = NULL` to `add_offspring()` and use it as the base seed?
- If `seed = NULL`, sample a base seed from the current R RNG, or use a fixed
  default?
- Should repeated calls with the same global `set.seed()` but no explicit
  `seed` reproduce?
- Is the base seed stored anywhere for audit, especially if `store_crossovers`
  is on?
- Does `add_founders()` keep base R `sample()` or also move to a seed argument
  / dqrng?

This decision is part of the public reproducibility contract and should not be
left as a Stage 1 implementation detail.

### 5. The dqrng stream-seeding scheme is still open but is foundational

The plan correctly lists the stream scheme as an open question
(`plans/optimize_add_offspring.md:747-751`), but the rest of the design assumes
it is solved. This should be resolved before coding, not during Stage 1.

The scheme must define:

- exact key tuple: `(base_seed, global_offspring_index, parent_origin,
  stream_kind)` where `stream_kind` is `"autosome"` or `"special"`;
- whether `global_offspring_index` is the original row in `matings`, not a batch
  or grouped index;
- how keys are converted to dqrng streams in both R and C++;
- collision guarantees if hashing is used;
- whether the R and C++ APIs can produce bit-identical uniforms from the same
  key without relying on undocumented behavior.

Recommended pre-implementation task: a tiny dqrng spike that proves R and C++
can generate identical uniform sequences for several `(seed, offspring, parent,
stream_kind)` keys before building the gamete kernel around it.

## Medium-Priority Concerns

### 6. The kernel input layout is underspecified for C++

The proposed signature says:

```r
parent_allele           # per-parent 2 x n_autosome_loci integer
parent_line_origin_code # per-parent 2 x n_autosome_loci integer
```

But `add_offspring()` has many unique parents. The plan needs a concrete memory
layout:

- 3D array `[parent, homolog, locus]`?
- flat integer vector with offsets?
- list of matrices for the R reference and packed matrix for C++?
- parent index ordering and whether `parent_1_idx` / `parent_2_idx` index this
  packed parent list;
- how missing/special-chromosome-only parents are represented.

C++ performance and correctness depend on this layout. It also affects whether
the R reference and C++ kernel can be compared directly.

### 7. Special chromosomes are partly specified, but their dqrng and crossover output path needs more detail

The plan keeps special chromosomes on the R path (`plans/optimize_add_offspring.md:
335-352`), which is reasonable. But after Stage 1 the special path cannot keep
using base `sample.int()` if `store_crossovers = TRUE` must not perturb allele
output and if thread/batch invariance is required.

The plan should specify:

- a dqrng-based replacement for `pass_through_gamete()` row choice when `k == 2`;
- a dqrng-based special-chromosome recombination path for `k == 2 &&
  recombines_*`;
- whether special-chromosome crossover rows are included in
  `ind_crossover` and how their stream key is derived;
- parity testing scope: the C++ kernel is autosome-only, so "R↔C++ parity on a
  special-chromosome config" must mean wrapper-level equality where autosomes
  are C++ and special chromosomes are shared R code, not that the C++ kernel
  handles special chromosomes.

### 8. `add_founders()` batching has an RNG trap

The plan says `add_founders()` skips the gamete generator but should move to the
same long write path and batching (`plans/optimize_add_offspring.md:646-652`).
Current `add_founders()` draws all haplotype indices in one `sample()` call
(`R/add_founders.R:200-205`), then builds a dense selected haplotype matrix
(`R/add_founders.R:244-257`).

If founder batching samples per batch, base-R RNG output changes. That may be
acceptable pre-1.0.0, but the plan also says founder long-write equivalence
should hold (`plans/optimize_add_offspring.md:716-717`).

Recommended fix: state the intended RNG behavior explicitly:

- If founder output must remain equivalent for the long-write refactor,
  pre-draw all haplotype indices once, then batch materialization/write only.
- If founder sampling may change, say so and replace equivalence with
  forward same-seed determinism plus distribution checks.

Also, `hap_data_matrix` itself is `n_haplotypes x n_loci`; at target scale this
may be large enough that direct long materialization should avoid building the
whole selected dense matrix even before writing.

### 9. `ind_crossover` schema should follow package ID/schema conventions

The plan says `id_crossover` is "auto-incrementing"
(`plans/optimize_add_offspring.md:379`). Current package convention for new
integer IDs is usually `next_int_id()`, and the coordinate plan explicitly
avoided DB auto-increment for `id_genome_map`.

Recommended fix:

- specify `id_crossover` assignment via `next_int_id(conn, "ind_crossover",
  "id_crossover")`;
- add `ind_crossover` to `TABLE_RESERVED_COLS`, `TABLE_PRIMARY_KEYS`,
  `TABLE_ROW_KEYS`, `SYSTEM_TABLES`, schema metadata, and `pop$tables`;
- say it is created by `define_genome()`, not `initialize_genome()` only
  (`plans/optimize_add_offspring.md:677-679` currently mentions
  `initialize_genome()`, which is deprecated).

### 10. Long writes should be transaction-scoped

Both current `add_offspring()` and the plan write multiple tables:

- `ind_meta`;
- `ind_haplotype`;
- optional `ind_crossover`;
- maybe schema/table creation the first time crossovers are enabled.

With batching, a failure midway could leave a partial generation. The plan does
not discuss transactions. This matters more once a batch loop can write
hundreds of millions of rows.

Recommended fix: define transaction boundaries. Prefer one transaction per
`add_offspring()` call if feasible; otherwise one transaction per batch with a
clear statement that completed earlier batches are committed. For correctness,
the first option is cleaner.

### 11. The Poisson inversion sampler needs numerical and parameter bounds

The Knuth sampler is fine for lambda near 1, which is the expected chromosome
scale. But the plan describes it as multiplying uniforms until product <
`exp(-lambda)` (`plans/optimize_add_offspring.md:323-332`).

For very large maps, `exp(-lambda)` underflows and the loop cost is `O(lambda)`.
That may be irrelevant in livestock-style maps, but the package aims to be
species-flexible and `cM_per_Mb` can be user-specified.

Recommended fix:

- implement the sampler using log accumulation rather than raw product;
- document an expected/max supported lambda for the exact inversion path;
- add a guard or alternate parity-safe algorithm for unusually large lambda;
- include tests around lambda 0, small lambda, typical lambda, and a larger
  stress lambda.

### 12. Testing "recorded crossovers reproduce haplotype switches" needs a controlled parent design

The testing section says `ind_crossover` rows and `ind_haplotype` switches should
agree (`plans/optimize_add_offspring.md:707-711`). That is only directly visible
if the two parental homologs are distinguishable at the relevant loci. If both
homologs carry the same allele or same `line_origin` over an interval, a homolog
switch can be invisible in `ind_haplotype`.

Recommended fix: use a controlled test fixture where homolog 1 and homolog 2 are
encoded differently at every marker, or expose/compare an internal homolog trace
from the generator in tests. Do not rely on random founder alleles to infer
switch points.

## Lower-Priority Plan Cleanups

### 13. Update stale status and line references

The implementation-order table still says the position-coordinate pre-stage is
"Not started" (`plans/optimize_add_offspring.md:766`), but that migration is now
implemented in the working tree. The line references to current code also need a
refresh after the coordinate fixes; for example, current `add_offspring()`
already has the per-map cache and the hot-loop line numbers have shifted.

### 14. Benchmark dependencies and deliverables need tightening

The plan names `profvis` and `bench`, but `DESCRIPTION` currently suggests only
`testthat`, `knitr`, and `rmarkdown`. Since the benchmark is a dev script, they
do not necessarily need to become package `Suggests`, but the plan should say
whether the benchmark script checks for them and skips gracefully.

The existing benchmark harness is broad storage benchmarking, not the step-level
profile the plan requires (`dev/benchmarks/benchmark_haplotype_scale.R:59-98`).
Phase 0 should explicitly require a new timing wrapper or instrumentation mode
that measures validation, parent load, parent matrix build, gamete generation,
line-origin gather, long-row assembly, and DB write separately.

### 15. Fallback wording around compiled code is optimistic

The plan says the internal selector can call the compiled kernel when available,
else the R reference, and that this "guarantees source installs without a
compiler still work" (`plans/optimize_add_offspring.md:596-600`). Once the
package contains C++ sources and `LinkingTo` dependencies, source installs
generally still require a working compiler unless the build is deliberately
configured to make compilation optional.

Recommended fix: either:

- require a compiler for source installs after Stage 2, while keeping the R path
  as a runtime/debug fallback; or
- design an explicit optional-compile configure path. That is extra complexity
  and probably not worth it right now.

### 16. Output ordering should be canonical

The long output vectors need a documented order for testing, deterministic DB
inserts, and R/C++ parity. For example:

```text
offspring original row order -> parent_origin 1 then 2 -> chromosome order -> locus order
```

If parent grouping is used for cache locality, it must not change the public
row order or the `o` index used for stream keys. The plan hints at grouping by
parent for locality, but should pin the distinction between execution order and
output order.

## Readiness Recommendation

Do not start implementation directly from the current plan. It is good enough
for technical discussion, but it should be revised before coding.

Minimum revisions before implementation:

1. Define the kernel input shape for multiple parents and multiple resolved
   maps.
2. Move batching/direct-long-write earlier or explicitly split Stage 0 into a
   wide-output seam extraction followed by long-write conversion.
3. Decide the `add_offspring()` seed/base-seed API and the exact dqrng stream
   derivation scheme.
4. Specify special-chromosome RNG/crossover behavior under dqrng.
5. Define `ind_crossover` schema registration and ID assignment conventions.
6. Add transaction boundaries and canonical output ordering.

Once those are pinned, the plan should be ready for professor review and then
implementation in small, testable stages.

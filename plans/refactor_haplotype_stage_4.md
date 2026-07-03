# Stage 4 — Sex Chromosomes & Organelles (COMPLETE)

**Status**: Done (tidybreed 0.49.0, 2026-07-03)
**Parent plan**: [refactor_haplotype.md](refactor_haplotype.md) (Stage 4 of 5)
**Scope decision**: sex chromosomes and organelles only (X/Y, Z/W, X0/Z0, MT, plastids). Real
polyploidy (ploidy > 2, uneven-ploidy crosses like 2N × 4N → 3N) is explicitly out of scope for
this pass — deferred to a future stage once a meiosis-pairing model is designed (autopolyploid vs.
allopolyploid pairing remains unresolved, per Open Questions in the parent plan).

---

## Outcome

- **Full test suite: 0 failures, 0 errors, 1192 passing** (up from 1055 at Stage 3) — 137 new
  tests across 2 new files (`test-define_chr.R`, `test-sex_chromosomes.R`) and 4 extended files.
- **Parity preserved**: `tests/testthat/test-parity.R`'s golden files are byte-identical, **no
  re-baseline** — direct evidence that restructuring `add_founders()`/`add_offspring()` to branch
  per chromosome did not change the diploid-autosome RNG sequence or output at all.
- Exploration during planning found the codebase was **fully diploid-hardcoded with zero
  `chr_meta` awareness** in `add_founders()`/`add_offspring()` — even the narrow X/Y-only scope
  required restructuring both functions to operate per-chromosome-group rather than as one flat
  whole-genome matrix.

---

## What shipped

### Schema (`R/define_genome.R`, `R/open_pop.R`)

- **`chr_meta.copies_M`/`copies_F`** (absolute `UTINYINT`, `0`/`1`/`2`) → **`copy_mode_M`/
  `copy_mode_F`** (`VARCHAR`, `"full"`/`"half"`/`"none"`, **relative** to an individual's own
  ploidy). This resolves the terminology clash identified during planning (see
  `refactor_haplotype.md`, "Ploidy vs. sex-linkage" section, dated 2026-07-03): the same field was
  previously used to mean both "sex-linked copy reduction" (0–2) and "polyploid baseline count"
  (4/6/8), which conflict in any species that is both polyploid and has a hemizygous sex
  chromosome. `hemi_parent`/`recombines` unchanged.
- **`ind_meta.ploidy`** — new `UTINYINT NOT NULL DEFAULT 2`, reserved column. Declared explicitly
  at `add_founders()` time (must be `2`); computed at `add_offspring()` time as
  `sire_ploidy %/% 2 + dam_ploidy %/% 2` (always `2` in this stage, but computed rather than
  hardcoded so this code doesn't need revisiting when real polyploidy ships).
- Greenfield rename, no migration — consistent with the project's established "rename and break"
  convention (Stage 1).

### `define_chr()` — new function (`R/define_chr.R`)

Upserts a chromosome's inheritance rule by `chr_name` (`ON CONFLICT` pattern, modeled on
`define_index()`, not `define_chip()`'s filtered-tidybreed_table pattern, since there's no
existing table being subset). Validates: `chr_name` exists in `genome_meta`; `copy_mode_M`/
`copy_mode_F` are one of the three valid values; `hemi_parent` is `NULL` iff both copy_modes are
`"full"`, and one of `"parent_1"`/`"parent_2"` otherwise. `overwrite = TRUE` by default (deliberate
deviation from `define_index()`'s `overwrite = FALSE` — re-calling to fix a typo'd `hemi_parent` is
a normal edit workflow, not append-only simulation output). RNG-neutral write (`duckdb_register()`
+ `dbExecute`), matching the Stage 1 invariant.

### `add_founders()` restructuring (`R/add_founders.R`)

New `ploidy = 2L` argument (errors if not `2`). The single flat founder-sampling draw
(`sample(1:n_haplotypes, size = n_founders * 2, replace = TRUE)`) is **untouched** — this was the
top identified risk during planning, since moving it inside a per-chromosome loop would have
changed the RNG draw count/sequence for every existing simulation, even ones that never touch
`chr_meta`. Only the **write step** became chromosome-aware: one UNPIVOT+INSERT per chromosome,
using `chr_meta.copy_mode_M`/`copy_mode_F` (for the founder's sex) to decide which
`parent_origin` slot(s) get the already-drawn haplotype's alleles for that chromosome's loci —
both slots for `"full"`, only the `hemi_parent`-designated slot for `"half"`, neither for `"none"`.

### `add_offspring()` restructuring (`R/add_offspring.R`)

Two compounding generalizations:

1. **Parent-haplotype loading** no longer hard-asserts exactly `2 × n_loci` rows per parent.
   Chromosomes are classified once per call as **plain autosome** (`copy_mode_M = copy_mode_F =
   "full"`, `recombines = TRUE` — chr_meta's default row) or **special** (anything else). Each
   parent's autosome haplotypes are loaded into a **compact** `2 × n_autosome_loci` matrix (special
   loci excluded entirely, not just zeroed) using a `full_to_local` position remap so
   `autosome_chr_info` (derived from `build_chr_info()`) can be re-expressed with local indices
   without touching `make_gamete()` itself. Special-chromosome haplotypes are loaded separately
   into compact `k × n_chr_loci` matrices per (parent, chromosome), `k` resolved from that
   specific parent's own sex via the new shared `resolve_chr_copy_count()` helper
   (`R/chr_meta_helpers.R`).

   (An earlier draft kept the autosome matrix at full `n_loci` width with special columns left at
   0. This looked simpler but was wrong: `make_gamete()`'s returned `homolog` vector defaults to
   `0` at untouched positions, and `0` is R's "drop this element" sentinel for matrix-style
   `m[cbind(r, c)]` indexing — so `line_origin` lookups silently returned a shorter-than-expected
   vector whenever a genome had more than one autosome alongside a special chromosome, throwing
   "number of items to replace is not a multiple of replacement length." Caught by the X0 test
   during Stage 4's own test-writing pass — see "What the X0 test caught" below.)

2. **Gamete generation** splits into an **autosome path** (byte-identical to pre-Stage-4 code —
   same single `make_gamete()` call per parent, now against the compact matrix/local-index
   descriptor pair, which is trivially identical to the original for the all-default case) and a
   **special-chromosome path**, executed strictly *after* the autosome path for each offspring
   (never reordered ahead of it, so special-chromosome RNG draws can never perturb autosome
   draws). For each special chromosome: resolve the offspring's own copy_mode; `"none"` writes
   nothing; `"half"` resolves the contributing parent via `hemi_parent`, then recombines across
   that parent's own two stored copies (via `make_gamete()`, scoped to that one chromosome's local
   descriptor) if the parent has 2 copies and `recombines = TRUE`, otherwise passes the parent's
   single stored copy through unchanged via the new `pass_through_gamete()`; `"full"` does the same
   resolution independently for both mating parents.

### `pass_through_gamete()` — new helper (`R/recombination_helpers.R`)

Handles non-recombining and single-copy inheritance (Y, W, MT, most organelles) without touching
`make_gamete()` (whose whole contract is "simulate crossovers" — reusing it for a chromosome that
never recombines would mean carrying spurious `rpois()`/`runif()` calls). When the contributing
parent has only 1 stored copy, the result is passed through with **zero** random draws — patrilineal
Y / matrilineal MT inheritance should never consume RNG state, and now provably doesn't (see
`test-sex_chromosomes.R`'s direct unit tests of this function).

### Guard relaxation (`R/ploidy_helpers.R`, `R/add_dosage.R`, `R/extract_genotypes.R`)

`assert_diploid_only()` (blocked *any* non-default `chr_meta` row) → `assert_ploidy_2()` (checks
`ind_meta.ploidy` for the candidate individual set). Confirmed by reading the code directly that
`add_dosage()`'s `SUM(h.allele) ... GROUP BY (id_ind, locus_id, ...)` and `extract_genotypes()`'s
zero-initialized dense matrix already produce correct results for any row count per locus — the
old guard was blocking a case that didn't actually need blocking. Only genuine non-diploid organism
ploidy remains guarded against.

### `add_tbv()` and `define_additive_effects.R`

`add_tbv()` needed **no changes** — its Stage 2 centered SQL already sums over however many
haplotype rows exist per individual per locus, with no hardcoded row-count assumption. Verified by
a new hand-calculated hemizygous-locus TBV regression test rather than by code inspection alone.

`rescale_effects_to_target()` (`R/define_additive_effects.R`) has a real, previously-unflagged gap:
its Falconer variance formula (`2 * p * (1-p) * a^2`) assumes every QTL is diploid/autosomal, and
would silently overstate variance contribution for a hemizygous locus (whose correct genic variance
is `p*(1-p)*a^2`). Rather than generalizing the formula (an approximation that gets murky for
mixed-sex populations), `define_additive_effects()` now calls a new `assert_qtl_autosomal()`
(`R/chr_meta_helpers.R`) before rescaling, erroring clearly if `scale_to_target = TRUE` and any
selected locus is on a non-`"full"`-copy_mode chromosome. Fail loud, consistent with this project's
existing `assert_diploid_only()`-style forward-defense pattern.

---

## What the X0 test caught

The X0 test (a hemizygous chromosome with no partner chromosome defined, alongside **two** plain
autosomes rather than one) failed during test-writing with `number of items to replace is not a
multiple of replacement length` — a real bug in the first implementation draft's "full-width matrix
with unfilled special columns" design (see "add_offspring() restructuring" above for the root
cause). This is exactly the kind of edge case the parent plan's own reasoning anticipated needing
explicit test coverage for ("X0/Z0 systems... flagged here so no future code assumes sex-determination
systems always pair two named chromosomes") — it surfaced a genuine implementation defect, not just
a documentation gap. Fixed by making the autosome path fully compact (local-indexed) rather than
full-width-with-holes; re-verified against the full parity suite afterward.

---

## Deviations from the Stage 4 plan (intentional)

- **Autosome matrices are compact (local-indexed), not full-width-with-holes.** The parent plan's
  design sketch (and this stage's own initial implementation attempt) described keeping parent
  haplotype matrices at full `n_loci` width for simplicity. This turned out to be incorrect for any
  genome with more than one plain autosome alongside a special chromosome (see above) — fixed by
  fully compacting both the autosome and special-chromosome matrices, with a `full_to_local` index
  remap for the autosome path's `chr_info` descriptors.
- **RNG unit test targets `pass_through_gamete()` directly, not an integration-level comparison.**
  An early test draft tried to prove "MT inheritance draws zero RNG values" by comparing two
  differently-sized genomes' autosome results — this was flawed because `define_genome()` splits
  loci evenly across `n_chr`, so varying total locus count while holding chromosome count fixed
  changes *every* chromosome's locus count and position layout, not just the one being varied.
  Replaced with a direct, unambiguous unit test of `pass_through_gamete()`'s RNG behavior instead
  (`k == 1` → zero draws; `k == 2` → exactly one `sample.int(2, 1)` draw).

---

## What's next (not done in Stage 4)

- **Stage 5** — sparse mutation storage (`add_mutation()`; `introduced_gen` column already
  present in `genome_meta`).
- **Real polyploidy** (ploidy > 2, uneven-ploidy crosses) — explicitly deferred; needs a resolved
  meiosis-pairing model (autopolyploid vs. allopolyploid) before `ind_meta.ploidy` can take values
  other than `2`. See Open Questions in `refactor_haplotype.md`.
- **Haplodiploidy** (bees, ants, wasps) and **biparental plastid inheritance** — documented
  non-goals per the parent plan; do not fit the `chr_meta` per-chromosome model at all (see
  "Haplodiploidy is architecturally different" in `refactor_haplotype.md`).

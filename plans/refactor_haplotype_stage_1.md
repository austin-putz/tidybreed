# Stage 1 — Wide → Long Haplotype/Genotype Storage (COMPLETE)

**Status**: Done (tidybreed 0.47.0, 2026-07-02)
**Parent plan**: [refactor_haplotype.md](refactor_haplotype.md) (Stage 1 of 5)
**Scope decisions**: greenfield (no migration built), benchmark deferred.

Stage 1 changed the *storage shape* of haplotypes/genotypes from wide to long
while preserving simulation behavior exactly. Stages 2–5 build on this and are
still open (see "What's next").

---

## Outcome

- **Full test suite: 0 failures, 0 errors, ~1005 passing.**
- **Parity preserved**: a deterministic seeded 2-line → F1 → F2 simulation
  produces semantically identical haplotypes, dosages, TBVs, and exported
  genotype matrices before vs. after the refactor.
- **`grep` clean**: no `genome_haplotype` / `genome_genotype` references in `R/`.

---

## What shipped

### Schema
- `genome_haplotype` → **`ind_haplotype`** (long): one row per
  (individual × haplotype × locus) —
  `id_ind, parent_origin, strand, line_origin, locus_id, locus_name, allele`;
  PK `(id_ind, parent_origin, strand, locus_id)`. `strand` always 1 (diploid).
- `genome_genotype` → **`ind_genotype`** (long, **on-demand cache**):
  `id_ind, locus_id, locus_name, dosage_value`; PK `(id_ind, locus_id)`.
  No longer auto-populated — filled only by `add_dosage()`.
- **`founder_haplotypes`** reshaped to long:
  `line_name, haplotype_id, locus_name, allele` (the old `hap_id` string is now
  an integer `haplotype_id` scoped per line).
- **`chr_meta`** new table with default diploid-autosome rows
  (`copies_M=2, copies_F=2, hemi_parent=NULL, recombines=TRUE`).
  `define_chr()` and non-default rules are Stage 4.
- `genome_meta` gains `introduced_gen` (NULL for founding loci; Stage 5 mutation
  marker); `locus_name` validated unique (R-side, survives the overwrite pattern
  in `.write_founder_haplotypes()`).
- The old wide tables are **removed** (greenfield — no compatibility layer).

### Behavior
- **`line_origin` populated on every allele from the start.**
  `make_gamete()` (`R/recombination_helpers.R`) now returns
  `list(allele, homolog)` — `homolog` is the per-locus contributing parental
  homolog, derived from values it already computed (`current_hap` /
  `hap_indices`), so the random-draw sequence is unchanged. Founders get their
  own line; offspring inherit `line_origin` per-locus from the segment that
  contributed each allele (correct across F1/F2/backcross). **TBV does not yet
  *use* `line_origin`** — that is Stage 2.
- **New `add_dosage(tbl, chip_name, locus_names, overwrite_dosage)**` — computes
  0/1/2 dosages from `ind_haplotype` into `ind_genotype`. Idempotent via
  `INSERT OR REPLACE`. Docs contrast it with `add_genotypes()` (which only marks
  animals as physically genotyped).
- **`extract_genotypes(loci_tbl=)`** new argument — a general locus filter via
  `get_table("genome_meta")` (e.g. autosomes only), unioned with
  `chip_name` / `effects_tbl`. Output columns kept as `locus_<id>` (preserves the
  existing contract; sourced from `ind_haplotype` via aggregation).
- Reader helpers now source from long: `get_genotype_matrix()`,
  `compute_base_allele_freq()` (both `founder_haplotypes` and `current_pop`
  bases) in `R/define_additive_effects.R`, `get_haplotype_matrix()` in
  `R/phenotype_helpers.R`, and the BLUPF90 genotype export in
  `R/blupf90_helpers.R`.

---

## The critical invariant: DB writes must not perturb the genetic RNG stream

Discovered mid-refactor and load-bearing for parity:

- **`dbAppendTable()` and `dbWriteTable()` advance R's `.Random.seed`** (they
  generate an internal random temp-view name). The leak amount is *constant*
  regardless of frame shape.
- **`duckdb_register()` + `INSERT ... SELECT`** and **`dbExecute("INSERT ...")`**
  are **RNG-neutral**. `dbGetQuery()` (reads) are also neutral.

Consequences applied throughout:
- All new/changed writes (chr_meta seed, `ind_haplotype` inserts) use the
  RNG-neutral primitives, so the genetic draw sequence is untouched and the
  wide-captured golden still matches the long code — no golden re-capture needed.
- `founder_haplotypes` kept using `dbWriteTable` (constant leak preserved) so
  seeded tests that sample right after `define_founder_haplotypes()` (e.g.
  `slice_sample` in `test-define_additive_effects.R`) don't shift.
- The parity harness re-seeds before every stochastic step, so incidental write
  leaks between generations never affect the golden.

**If you touch a writer in later stages: never introduce a `dbAppendTable`/
`dbWriteTable` on the pre-sampling path unless you intend the RNG shift.**

---

## How writes work now (pattern for later stages)

`add_founders()` / `add_offspring()` still build the **wide** allele matrix in R
exactly as before, then let DuckDB do the wide→long explosion via `UNPIVOT`
(so R never materializes the 2·n_ind·n_loci long frame):

```sql
INSERT INTO ind_haplotype (...)
SELECT u.id_ind, u.parent_origin, 1, <line_origin>, gm.locus_id, gm.locus_name, u.allele
FROM (UNPIVOT __tmp_hap ON COLUMNS(* EXCLUDE (id_ind, parent_origin))
      INTO NAME locus_col VALUE allele) u
JOIN genome_meta gm ON u.locus_col = 'locus_' || gm.locus_id
```

`add_offspring()` additionally builds a parallel wide `line_origin` frame and
joins the two UNPIVOTs on `(id_ind, parent_origin, locus_col)` so each allele
carries its inherited line.

---

## Verification harness (keep this — Stages 2–5 rely on it)

- `tests/testthat/helper-parity.R` — `run_parity_sim()` plus format-agnostic
  `canonicalize_haplotypes()` / `canonicalize_dosage()` / `canonicalize_export()`
  helpers.
- `tests/testthat/test-parity.R` — captures golden RDS on first run, compares on
  every subsequent run. **To intentionally re-baseline** (after an approved
  behavior change, e.g. Stage 2 changes TBV), delete
  `tests/testthat/parity_golden/` and re-run.
- `tests/testthat/test-long-schema.R` — long DDL / `chr_meta` / `introduced_gen`
  checkpoints.
- `tests/testthat/test-long-writers.R` — `ind_haplotype` structure + `line_origin`
  inheritance (F1 one line per parent_origin; F2 mosaic of both).
- `tests/testthat/test-add_dosage.R` — dosage correctness, idempotence, chip /
  locus_names filters, `overwrite_dosage`.

Run the suite: `Rscript -e 'devtools::test()'` (Mac/Linux) or the Windows Rscript
path in CLAUDE.md.

---

## Deviations from the Stage 1 plan (intentional)

- **`add_tbv` was not rewritten to centered `GROUP BY` SQL.** It still uses the
  matrix helpers, now sourced from `ind_haplotype`, producing identical results
  with lower risk. The centered-SQL rewrite belongs in **Stage 2**, where
  `add_tbv` is rewritten for `line_origin` anyway — doing it twice was wasteful.
- **`import_founder_haplotypes()` was not added** (it was in the parent plan's
  file table but not the Stage 1 sub-plan). Additive convenience; defer.
- A now-unreachable wide-fallback branch remains in `helper-parity.R` — harmless
  (documents the transition; the harness always reads long).

---

## What's next (not done in Stage 1)

- **Stage 2 — line-origin TBV**: make `add_tbv()` *use* `line_origin` for
  crossbreeding additive TBV (centered SQL for population-wide + line-specific
  effects with fallback). `line_origin` is already populated, so F2 tests can
  reuse Stage 1 data. When TBV output changes, re-baseline the parity golden.
- **Stage 3** — dosage cache / export hardening beyond current `add_dosage()`.
- **Stage 4** — `define_chr()`, sex chromosomes, organelles, polyploidy
  (`chr_meta` shape is ready).
- **Stage 5** — sparse mutation storage (`add_mutation()`; `introduced_gen`
  column already present).
- **`migrate_haplotype_to_long()`** — only if pre-0.47 databases ever need
  preserving (greenfield today).
- **Scale benchmark** (5000 × 50k) — recommended before large runs;
  `summary_pop` and the matrix readers now build wide frames from the long table.

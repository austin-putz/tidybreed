# Doubled haploids — inbred founders and `add_doubled_haploids()`

**Status:** design proposal, nothing implemented. **Created:** 2026-09-19.
**Baseline:** v0.68.4 on `feat/genome-effects-v49`.

---

## Summary of decisions

1. **Nothing changes in `define_founder_haplotypes()`.** A founder pool is a bag of
   *unpaired* haplotypes; "doubled" is a statement about how the two slots of an
   individual get filled, and the pool has no slots. The request "simulate doubled
   haploids from founders" is a pairing rule and belongs in `add_founders()`.
2. **`add_founders()` gains `inbred = FALSE`.** When `TRUE`, each founder draws one
   pool haplotype and carries it in both `parent_origin` slots. Not called "doubled
   haploid" because nothing was haploid and nothing was doubled — no gamete was made.
3. **A new exported `add_doubled_haploids()`** performs the real DH operation: one
   meiosis from a heterozygous parent (typically an F1), the gamete written to both
   slots. It mirrors `add_offspring()` — `pop` + a one-row-per-individual tibble.
4. **Not one merged function.** The two inputs (pool vs individuals) differ in
   meiosis (none vs one), pedigree (`NULL` vs `parent = parent`), RNG (one base-R
   `sample()` vs per-gamete `dqrng` streams), and ploidy source (declared vs
   computed). A function dispatching on the class of its first argument would carry
   two mutually exclusive argument sets and spend its docs explaining which apply
   when. One transmission mechanism per exported verb; one engine underneath.
5. **The sharing happens below the exported surface.** `add_offspring()`'s §7–§12
   (parent resolution → kernel → writes) is refactored into an internal progeny
   engine parameterised by an *assembly rule*; `add_offspring()` and
   `add_doubled_haploids()` become thin callers. This is also the seam a future
   `add_clones()` (zero gametes) plugs into.
6. **No schema change.** No new table, no new column. `ind_meta`, `ind_haplotype`,
   `ind_crossover` hold everything as they are today.

---

## Terminology (so the docs say one thing)

| Term | Meaning here | Where it lands |
|---|---|---|
| **Inbred founder** | A founder whose two haplotypes are the same pool haplotype (F = 1 by construction, no pedigree). | `add_founders(inbred = TRUE)` |
| **Doubled haploid (DH)** | An individual produced by one meiosis of a single parent, with the resulting gamete duplicated. Fully homozygous; a mosaic of the parent's two haplotypes. | `add_doubled_haploids()` |
| **Selfing** | Two independent gametes from the same parent. Already supported: `add_offspring()` with `id_parent_1 == id_parent_2` (there is no sex check on parents). | `add_offspring()` — document it |

"Doubled haploid", not "double haploid": the standard term in plant breeding is
*doubled* haploid (DH). Function name follows the no-abbreviation rule:
`add_doubled_haploids()`, not `add_dh()`.

---

## Part A — `add_founders(inbred = TRUE)`

### Behaviour

Today (`R/add_founders.R:211`):

```r
hap_indices <- matrix(
  sample(1:n_haplotypes, size = n_founders * 2, replace = TRUE),
  nrow = n_founders, ncol = 2)
```

Proposed:

```r
if (inbred) {
  # One pool haplotype per founder, WITHOUT replacement: inbred founders are
  # meant to be distinct lines. Both slots carry the same pool row.
  if (n_founders > n_haplotypes) stop(...)
  one <- sample(n_haplotypes, size = n_founders, replace = FALSE)
  hap_indices <- cbind(one, one)
} else {
  hap_indices <- matrix(sample(1:n_haplotypes, n_founders * 2, replace = TRUE),
                        nrow = n_founders, ncol = 2)   # unchanged, byte-identical
}
```

Everything after step 4 is untouched: the batched long write already reads
`hap_indices[, 1]` for slot 1 and `hap_indices[, 2]` for slot 2, `line_origin` is the
founder's line either way, `chr_inheritance` routing (which slots are emitted per
chromosome for the founder's sex) is independent of which pool row fills them.

### Decisions

- **Without replacement.** With replacement, two "inbred lines" could be identical
  individuals, which is a clone, not a second line. `inbred = TRUE` therefore requires
  `n_males + n_females <= n_haplotypes` (after any `filter()` on the pool) and errors
  otherwise with a message naming both counts. Users who want duplicates can call
  `add_founders()` twice.
- **RNG.** The `inbred = FALSE` path is byte-identical to today (same `sample()`
  call, same stream consumption). `inbred = TRUE` is a new stream — no
  golden-output concern (pre-1.0.0).
- **Hemizygous/absent chromosomes** need no special handling: a male founder's X is
  emitted from slot 2 only, and slot 2 reads the same pool row as slot 1 would have.
  The per-parent `copy count > 1` error stays.
- **`ploidy`** stays `2`. An inbred diploid is still a diploid.

### Variance consequence (document, do not "fix")

`define_additive_effects(scale_to_target = TRUE)` sets
`V_A = Σ_j n_eligible,j · p_j q_j a_j²` — the HWE / outbred genic variance. Among
fully inbred founders every locus is homozygous, so the realised additive variance
**among the founders themselves** is `Σ 2·(2 p_j q_j a_j²)` — **twice the target**
(F = 1 inflates genic variance by `1 + F`). This is the correct biology, not a bug:
the target describes the random-mating population the inbred lines were derived
from, and it is what re-appears after one generation of crossing among them.

Actions:
- `@param inbred` roxygen carries this note.
- `define_additive_effects()` `scale_to_target` docs gain one sentence: "realised
  among individuals with inbreeding `F` is `(1 + F)` times this quantity".
- No warning at run time — `add_founders()` runs before effects exist, and
  `define_additive_effects()` does not know how founders were paired without a
  haplotype scan. Not worth a scan for a documented property.

### Files

| File | Change |
|---|---|
| `R/add_founders.R` | New formal `inbred = FALSE` (after `ploidy`, before `...`); `.check_flag()` validation; the `hap_indices` branch above; roxygen `@param inbred`, a `@details` paragraph, one example; the completion message appends `" (inbred)"` when `TRUE` |
| `R/define_additive_effects.R` | One-sentence `(1 + F)` note under `scale_to_target` |
| `man/add_founders.Rd`, `man/define_additive_effects.Rd` | regenerate |
| `tests/testthat/test-add_founders.R` | see Tests |
| `NEWS.md`, `DESCRIPTION` | bump |

---

## Part B — `add_doubled_haploids()`

### Signature

```r
add_doubled_haploids(pop, parents, seed = NULL, store_crossovers = FALSE,
                     batch_size = NULL, max_batch_mem = NULL)
```

`parents` is a tibble/data.frame, **one row per doubled haploid**:

| Column | Required | Notes |
|---|---|---|
| `id_parent` | yes | The single parent; must exist in `ind_meta`. Written to **both** `id_parent_1` and `id_parent_2`. No `id_sire`/`id_dam` aliases — there is no sire/dam here. |
| `sex` | yes | `"M"`/`"F"` — same constraint as `add_offspring()` (hermaphrodite support is a separate question; this plan does not touch it) |
| `line_name` | yes | Line of the new individual; drives `{line_name}_{n}` IDs, same as `add_offspring()` |
| extra columns | no | Forwarded to `ind_meta` via `prepare_extra_cols()`, exactly like `matings` |

`n_per_parent` is deliberately **not** an argument: one row per individual is how
`matings` already works (`rep()` the parent id). Two conventions for "how many" would
be worse than one.

Returns `pop` invisibly with `attr(pop, "base_seed")`, and prints
`"Added N doubled haploids (base_seed = S)"`, matching `add_offspring()`.

### Naming note

`id_parent` is not a column of any table, which strains naming rule 1 ("argument
names match the column they populate"). The alternative — require `id_parent_1` and
silently copy it to `id_parent_2` — makes the user type a lie. `id_parent` is the
honest name for "the one parent"; the roxygen says explicitly that it fills both
pedigree columns.

### Semantics

For DH row `i` with parent `P`, offspring sex `S`, offspring line `L`:

1. **One gamete** is drawn from `P` using `P`'s resolved map (`P`'s sex, `P`'s
   line), on the stream `(o = i, parent_origin = 1, kind)`. Autosomes go through
   `make_gametes_batch()` (C++ or R reference, unchanged); special chromosomes go
   through the existing R special path with `r = 1` only.
2. The gamete's `(locus_id, allele, line_origin)` rows are written **twice**:
   `parent_origin = 1` and `parent_origin = 2`, both `strand = 1`. The duplication is
   an R-side row copy before `.write_long_haplotypes()`; **no kernel change**, so
   R↔C++ parity is inherited, not re-proven.
3. `ind_meta`: `id_parent_1 = id_parent_2 = P`, `sex = S`, `line_name = L`,
   `ploidy = 2` (computed as `2 × (ploidy_P %/% 2)`, asserted `== 2L`, same style as
   `add_offspring()`).
4. `ind_crossover` (when `store_crossovers = TRUE`): the gamete's crossovers are
   recorded **once**, under `parent_origin = 1`. One meiosis happened; writing them
   under both slots would double-count in every downstream crossover summary. The
   roxygen says so, and says that the absence of `parent_origin = 2` rows for a DH is
   expected.
5. `line_origin` flows through the gamete, so a DH of an A×B F1 is a homozygous
   mosaic of A- and B-labelled segments. `add_tbv()`/`add_tgv()` need no change:
   the evaluator reads labels per copy, and crossbred/line-specific effects and
   dominance/epistasis at homozygous states all evaluate correctly as-is.

### Chromosome rules

DH is only defined where the offspring inherits exactly one copy from each parent
slot and the parent has two copies to recombine. Resolution, per DH row:

- **Offspring side** (`resolve` at `S`, `L`): every chromosome must have
  `from_parent_1 == 1 && from_parent_2 == 1`. A hemizygous (`0,1`/`1,0`) or absent
  (`0,0`) chromosome for the offspring's sex → **hard error** naming the chromosome
  and sex ("doubled-haploid derivation is not defined for a chromosome the offspring
  carries in fewer than two copies").
- **Parent side**: the parent must hold `k = 2` copies (its own sex/line rule). `k
  = 1` (e.g. a male parent's X when the offspring is female) → hard error, same
  family as `add_offspring()`'s existing `k_p` checks.
- **Recombination** (`resolve` at `P`'s sex, `P`'s line): a chromosome that is
  `1,1` for both sexes but non-recombining in `P`'s sex (achiasmy) is classified
  "special" by `is_plain_autosome()` and passes through **one randomly chosen**
  homolog via `pass_through_gamete()` (one `dqrunif`, as today). That is the
  correct DH result for an achiasmatic parent — the whole chromosome is one
  parental homolog, doubled.

Net: in the common plant case (all autosomes, both sexes recombine) the check is a
no-op and everything goes through the fast kernel path.

### RNG design

- `seed` semantics identical to `add_offspring()`: `NULL` draws one base seed from
  base-R (so `set.seed()` upstream reproduces the call); an integer is used as-is;
  the resolved value is messaged and attached as an attribute, never persisted.
- Stream id: `.gamete_stream_id(o, parent_origin = 1L, kind)` with `o` = the
  **global** `parents` row index. `parent_origin = 2` streams are never seeded (no
  draws happen for the copy). Output is therefore independent of `batch_size`, map
  grouping, and whether special chromosomes exist — the same invariants
  `add_offspring()` already tests.
- The int32 bound on `o` is the same `(.Machine$integer.max - 3L) %/% 4L`; reuse the
  guard verbatim (it lives in the engine after Part C).
- `add_doubled_haploids()` and `add_offspring()` are separate calls with separate
  base seeds, so using `r = 1` in both is not a collision.

### What does **not** change

- Schema: none. `schema()`/`describe_table()`/`_schema_meta`: untouched (no new
  table). `archive_replicate()`: untouched (`ind_haplotype`/`ind_crossover` are
  already `reset_only`).
- Kernels: `make_gametes_batch_r()`, `make_gametes_batch_cpp()`,
  `.draw_chr_recombination()`, `pass_through_gamete()`: untouched.
- `add_dosage()`, `add_tbv()`, `add_tgv()`, `add_phenotype()`: untouched; DH rows
  look like any other diploid to them.

### Known limitation to document (not to solve)

Pedigree-based inbreeding for a DH: with `id_parent_1 == id_parent_2 == P`, the
tabular A-matrix yields `F_DH = 0.5 · (1 + F_P)`, not `1`. Genomic F is exact.
Users running pedigree BLUP on DH lines should know this; it is a property of
pedigree algebra, not of tidybreed. Storing an `is_dh` flag in `ind_meta` would be
metadata (CLAUDE.md principle 6); users who want it pass `dh = TRUE` as an extra
column in `parents`. The roxygen `@details` carries this paragraph.

---

## Part C — the shared progeny engine (refactor before Part B)

`add_offspring()` is 823 lines and the most intricate file in the package. Copying
its §7–§12 into a second function would be the exact "two implementations of the
same algorithm" CLAUDE.md warns against. Instead, extract an internal engine first,
as a **pure refactor** (byte-identical output for `add_offspring()`, verified by the
existing seeded tests), then add DH as a second caller.

### What is already parent-role-agnostic in `add_offspring()`

Reading `R/add_offspring.R`, the following sections depend only on "a set of parent
ids" and "a list of gametes", not on there being two distinct parents:

| § | Content | Role-agnostic? |
|---|---|---|
| 7 | Validate parent ids exist; fetch sex/ploidy/line | yes, given `unique_parents` |
| 8 | Chromosome classification; `build_gamete_info()`; per-`(sex,line)` map cache | yes |
| 9 | Load parent haplotypes; compact autosome + special matrices per parent | yes |
| 10 | Offspring ID numbering per line | yes, given `line_name` vector |
| 11 | Packed kernel inputs (`parent_allele`, `parent_lo_code`, `lo_levels`, `chr_arrays_by_key`) | yes |
| 12 | Seed resolution, batch loop, transaction, `ind_meta` append | **mostly** — see below |

### The five places that know about two parents

1. `unique_parents <- unique(c(matings$id_parent_1, matings$id_parent_2))`
2. `gam_o / gam_origin / gam_parent / gam_key` construction (2 gametes per
   offspring, interleaved `1,2,1,2,…`).
3. The special-chromosome loop `for (r in c(1L, 2L))` and its
   `pid <- if (r == 1L) sire else dam`.
4. `offspring_ploidy <- ploidy[p1] %/% 2L + ploidy[p2] %/% 2L`.
5. `ind_meta_new` (`id_parent_1`, `id_parent_2` columns) and the row-duplication
   step that DH needs and biparental does not.

### Proposed engine

```r
.write_progeny(pop, plan, assembly = c("biparental", "doubled"),
               seed, store_crossovers, batch_size, max_batch_mem, verb)
```

`plan` is a normalised tibble with columns `id_parent_1`, `id_parent_2`, `sex`,
`line_name`, plus extra user columns. Both callers normalise into it:

- `add_offspring()`: alias handling (`id_sire`/`id_dam`), then pass through.
- `add_doubled_haploids()`: `id_parent_1 <- id_parent_2 <- parents$id_parent`.

`assembly` controls the five points above:

| | `"biparental"` | `"doubled"` |
|---|---|---|
| Gametes per offspring | 2, origins `(1, 2)`, parents `(p1, p2)` | 1, origin `1`, parent `p1` |
| Special-chromosome `r` loop | `c(1L, 2L)` | `1L` |
| Chromosome precondition | existing (`cnt_r > 1` error) | additionally: offspring `1,1` on every chr; parent `k == 2` |
| After kernel/special output | write as-is | duplicate every row with `parent_origin = 2`; crossovers not duplicated |
| `ploidy` | `p1 %/% 2 + p2 %/% 2` | `2 * (p1 %/% 2)` |
| Message verb | `"offspring"` | `"doubled haploids"` |

The row duplication for `"doubled"` is one step after `long_frame` is assembled and
before `.write_long_haplotypes()`:

```r
if (assembly == "doubled") {
  dup <- long_frame; dup$parent_origin <- 2L
  long_frame <- rbind(long_frame, dup)
}
```

(In the preallocated autosome path it is cleaner to allocate `2 × n_auto` up front
and fill the second half with `ai_po = 2L` — same bytes, no `rbind`.)

### Refactor rules

- **Byte-identical for `add_offspring()`.** Every existing seeded test in
  `test-add_offspring.R`, `test-make_gametes_parity.R`, `test-parity.R` must pass
  unchanged. The `dbWriteTable(ind_meta)` call stays last so base-R RNG position is
  unchanged (the comment at §12 explains why).
- The engine is `@keywords internal`, not exported. `add_offspring()`'s roxygen is
  the user-facing documentation of the recombination model; `add_doubled_haploids()`
  `@inherits`/cross-references it rather than restating it.
- File placement: `R/progeny_engine.R` (new) holds `.write_progeny()` and the
  private helpers it pulls out of `add_offspring.R`. `add_offspring.R` shrinks to
  validation + normalisation + one call. `add_doubled_haploids.R` (new) is the same
  shape.

### Why do Part C at all rather than a `dh` flag inside `add_offspring()`

A `matings$dh` column was considered and rejected: a mating has two parents and a DH
has one, so users would write `id_parent_1 = id_parent_2 = "F1_3", dh = TRUE` to
get a DH — a leaky abstraction that also invites the wrong reading "selfing with a
flag". The engine gives the same code reuse with an honest surface.

---

## Implementation order

| Stage | Work | Verifies |
|---|---|---|
| **0** | Freeze decisions above (this document). | — |
| **1** | Part A: `add_founders(inbred = TRUE)` + tests + docs. Independent of everything else; ship first. | `test-add_founders.R` |
| **2** | Part C: extract `.write_progeny()` from `add_offspring()`, `assembly = "biparental"` only. **No behaviour change.** | full existing suite; seeded golden-within-current-code checks in `test-add_offspring.R` |
| **3** | Part B: `assembly = "doubled"` + `add_doubled_haploids()` + tests + docs. | new `test-add_doubled_haploids.R`; parity test |
| **4** | Docs sweep: `add_offspring()` `@details` gets a "Selfing" bullet (`id_parent_1 == id_parent_2` works today); `_pkgdown.yml`; vignette paragraph; `CLAUDE.md`; `NEWS.md`. | `devtools::document()`, `pkgdown::check_pkgdown()` |

Stages 1 and 2 can be separate commits/PRs. Stage 3 depends on 2.

---

## Files to alter — checklist

### New
- [ ] `R/progeny_engine.R` — `.write_progeny()` and extracted helpers (Stage 2)
- [ ] `R/add_doubled_haploids.R` — exported wrapper (Stage 3)
- [ ] `tests/testthat/test-add_doubled_haploids.R` (Stage 3)
- [ ] `man/add_doubled_haploids.Rd` (generated)

### Modified
- [ ] `R/add_founders.R` — `inbred` (Stage 1)
- [ ] `R/add_offspring.R` — becomes a thin caller (Stage 2); `@details` selfing bullet (Stage 4)
- [ ] `R/define_additive_effects.R` — `(1 + F)` sentence under `scale_to_target` (Stage 1)
- [ ] `NAMESPACE` — `export(add_doubled_haploids)` (Stage 3, via roxygen)
- [ ] `_pkgdown.yml` — add `add_doubled_haploids` under **Individuals** after `add_offspring` (Stage 3)
- [ ] `vignettes/tidybreed-introduction.Rmd` — short "Inbred founders and doubled haploids" subsection, or a pointer (Stage 4)
- [ ] `CLAUDE.md` — new subsection under Implemented Functions for `add_doubled_haploids()`; `add_founders()` entry mentions `inbred`; the `ind_haplotype` table note on row counts unchanged (Stage 4)
- [ ] `dev/package_summary/package_summary.md` — regenerate (Stage 4)
- [ ] `NEWS.md`, `DESCRIPTION` — one bump per shipped stage
- [ ] `tests/testthat/test-add_founders.R` (Stage 1)
- [ ] `tests/testthat/test-add_offspring.R` — add an explicit selfing test if none exists (Stage 2/4)

### Untouched (confirm, do not edit)
- `R/define_founder_haplotypes.R`, `R/recombination_helpers.R` (kernels),
  `src/make_gametes.cpp`, `R/define_genome.R`, `R/schema.R`,
  `R/archive_replicate.R`, `R/add_tbv.R`, `R/add_tgv.R`, `R/genome_effects_eval.R`

---

## Tests

### Stage 1 — `add_founders(inbred = TRUE)`
- Every founder is homozygous at every locus: for each `id_ind`, the
  `parent_origin = 1` and `= 2` allele vectors are identical (`ind_haplotype`
  self-join, `COUNT(*) WHERE a1.allele <> a2.allele == 0`).
- Each founder's haplotype equals *some* pool haplotype (join back to
  `founder_haplotypes` on all loci).
- Without replacement: `n_founders` distinct pool rows used; `n_founders >
  n_haplotypes` errors with a message naming both numbers.
- With a filtered pool (`filter(line_name == "A")`), the bound is the filtered count.
- `inbred = FALSE` under `set.seed()` is byte-identical to the current output
  (assert against a run of the *same* code path, not an old golden).
- Hemizygous chromosome (`define_chromosome("X", offspring_sex = "M", 0, 1)`):
  male inbred founders have one X row-set, female two, both slots identical where
  both exist.
- `inbred` validation: non-flag values rejected via `.check_flag()`.
- Sanity on the variance note: with `define_additive_effects(scale_to_target =
  TRUE)` and `base = "founder_haplotypes"`, `var(tbv_value)` among many inbred
  founders is ≈ `2 × target_add_var` (loose tolerance; documents the behaviour).

### Stage 2 — engine extraction (no behaviour change)
- All existing `add_offspring()` tests pass unchanged.
- Add, if missing: a seeded `add_offspring()` run captured **before** the refactor
  in the same session (`local_test_pop` helper), compared byte-for-byte to the run
  **after** — done once during development, then replaced by the ordinary
  "same seed → same output within current code" test.

### Stage 3 — `add_doubled_haploids()`
- **Homozygosity**: every DH is homozygous at every locus (same self-join as Stage 1).
- **Mosaic of parent**: at every locus the DH allele equals the parent's
  `parent_origin = 1` or `= 2` allele (never a third value); with `store_crossovers
  = TRUE`, the switch points in the DH's haplotype coincide with the recorded
  `pos_cM` (between adjacent loci).
- **line_origin**: for an A×B F1 parent, DH `line_origin` values are a subset of
  `{"A","B"}` and both appear across a chromosome with ≥1 crossover.
- **Pedigree**: `id_parent_1 == id_parent_2 == id_parent`; `ploidy == 2`;
  `line_name`/`sex` as supplied; extra columns forwarded (typed, `prepare_extra_cols`).
- **Crossovers recorded once**: `ind_crossover` rows for a DH all have
  `parent_origin = 1`; none have `2`.
- **Batch invariance**: `batch_size = 1`, `= 7`, `= nrow(parents)` under the same
  `seed` produce identical `ind_haplotype` and `ind_crossover`.
- **R↔C++ parity**: `options(tidybreed.kernel = "r")` vs C++ under the same seed →
  identical output (the engine calls `make_gametes_batch()`, so this is one test,
  not a new kernel contract).
- **`seed = NULL` under `set.seed()`** reproduces; `attr(pop, "base_seed")` is set
  and matches the message.
- **Chromosome-rule errors**: female DH from a male parent with an X defined as
  `(M: 0,1)` → parent `k = 1` error; male DH where Y is `(M: 1,0)` → offspring
  `1,0` error; achiasmatic-parent chromosome (`recombines = FALSE` for the parent's
  sex) → DH carries exactly one parental homolog for that chromosome, doubled.
- **Validation**: missing `id_parent`; parent not in `ind_meta`; `id_sire` given
  (rejected, with a message pointing to `id_parent`); bad `sex`; bad `line_name`;
  reserved extra column names; zero rows.
- **Downstream**: `add_tbv()` on DH individuals equals `2 × Σ (allele − p) a` over
  the single haplotype (hand-computed for a tiny genome); `add_tgv()` with a
  dominance term evaluates only homozygous states (dominance contribution is the
  homozygous-state value, never the heterozygote's).
- **Selfing regression** (Stage 4 docs claim): `add_offspring()` with
  `id_parent_1 == id_parent_2` produces two *independent* gametes (some loci
  heterozygous in a heterozygous parent), distinguishing it from DH.

---

## Documentation

### `add_founders()` roxygen
- `@param inbred` — definition, without-replacement rule and its bound, the
  `(1 + F)` variance note, and the pointer: "for doubled haploids derived by meiosis
  from existing individuals see `add_doubled_haploids()`".
- One example: inbred maize-style base population.

### `add_doubled_haploids()` roxygen
- `@description`: one meiosis, gamete duplicated, homozygous mosaic.
- `@param parents`: the table above, explicit that `id_parent` fills both pedigree
  columns; no aliases.
- `@details`: chromosome-rule preconditions; crossovers recorded once under
  `parent_origin = 1`; RNG/`seed` (`@inheritParams add_offspring` for `seed`,
  `store_crossovers`, `batch_size`, `max_batch_mem`); the pedigree-F caveat;
  comparison table Inbred founder / DH / Selfing (from *Terminology*).
- `@seealso add_offspring(), add_founders()`.
- Example: F1 from two inbred lines → 200 DH lines → `add_tbv()`.

### `add_offspring()` roxygen
- `@details` "Mating design flexibility" gains: "**Selfing**: use the same id for
  `id_parent_1` and `id_parent_2`; two independent gametes are drawn. For a
  doubled haploid (one gamete, duplicated) use `add_doubled_haploids()`."

### `CLAUDE.md`
- Function Naming table unchanged (`add_` fits: simulation output).
- New "### `add_doubled_haploids()`" subsection under Implemented Functions,
  after `add_offspring()`'s (there is no `add_offspring()` subsection today — add a
  short one alongside, or fold both into one "Progeny" subsection that also names
  `.write_progeny()` and the assembly rule).
- `add_founders()` subsection: mention `inbred`.

---

## Open questions (decide before Stage 3)

1. **Argument name for the DH design table.** `parents` is proposed (parallel to
   `matings`). Alternatives: `plan`, `derivations`. `parents` reads slightly like a
   character vector; the roxygen's first sentence should dispel that.
2. **Should `add_doubled_haploids()` accept a `parents` row whose parent is itself
   fully homozygous?** Biologically pointless (the DH equals the parent) but not
   wrong. Proposed: allow silently; no haplotype scan to detect it.
3. **`store_crossovers` for DH** — record under `parent_origin = 1` (proposed) vs a
   sentinel. `1` is truthful (that stream produced it) and needs no schema note.
4. **Hermaphrodite `sex`.** Out of scope; DH does not make the `M`/`F` requirement
   worse than `add_offspring()` already does. Track separately.
5. **Future `add_clones()`.** `assembly = "clonal"` (zero gametes, copy both parent
   slots, no RNG) drops into the same engine. Not built now; the engine's `assembly`
   argument is the reserved dimension, at zero cost.

---

## Non-goals

- No change to `define_founder_haplotypes()` (see Summary §1).
- No `is_dh` / `origin_type` column in `ind_meta`.
- No merged founders-or-parents function.
- No haploid individuals (`ploidy = 1`) — a DH is stored as the diploid it is.
- No polyploid DH (`ploidy > 2` remains unsupported everywhere).

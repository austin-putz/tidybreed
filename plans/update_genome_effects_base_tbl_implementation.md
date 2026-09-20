# `base_tbl` — implementation summary

**Implements:** `plans/update_genome_effects_base_tbl.md` (v5, approved).
**Shipped as:** v0.69.0 on `feat/genome-effects-v49`, 2026-09-20, in the three
commits the plan's §3.7 called for plus one docs commit.
**Result:** every decision in the plan is implemented as written. Four small
deviations, all recorded below with the reason; none changes the surface.

| Commit | Plan step | What |
|---|---|---|
| `48ebce2` | §3.7 step 1 | `extract_allele_freq()`, `.validate_base_tbl()`, `.founder_base_empty_error()`; tests 1–9 |
| `6a7d3a5` | step 2 | `define_additive_effects()` rewired; `base`, `base_line_name`, `compute_base_allele_freq()` deleted; 73 call sites migrated incl. the six vignette calls; tests 10–17 |
| `845d1cd` | step 3 | `base_freq` on `.ge_build()`, `base_tbl` on `define_genome_effects()`; tests 18–22 |
| (docs) | steps 4–5 | roxygen family sentence on both writers, `CLAUDE.md`, `NEWS.md`, `DESCRIPTION` 0.69.0, `plans/TODO.md` |

---

## What was built

### `R/extract_allele_freq.R` (new, exported)

Exactly the §3.2 sketch. Three shapes dispatched on `table_name`; the user's
filter rendered with `dbplyr::sql_render()` and embedded as a subquery; the
individual path a `JOIN (SELECT DISTINCT id_ind FROM (...))`; a final `LEFT
JOIN genome_meta` so an absent locus is `NA`; an all-`NA` result is an error
(founder shape → the "Available: …" diagnostic, which names an unnamed pool as
such); every non-`NA` value checked finite and in `[0, 1]`; output typed and
`locus_id`-ordered. Never warns.

`.validate_base_tbl(base_tbl, pop = NULL, arg)` — class, `validate_tidybreed_pop()`,
same-connection when `pop` is given, projected columns per shape
(`founder_haplotypes`: `locus_name`, `allele`; `ind_haplotype`: `locus_id`,
`allele`; else `id_ind`). The public helper calls it with `pop = NULL`.

### `define_additive_effects()`

Signature now `(tbl, trait_name, effects, distribution, G, method, base_tbl,
line_name, parent_origin, scale_to_target, seed)`. Three new internals:

- `.dae_resolve_base(pop, base_tbl, line_name)` → `list(p_base, label)`.
  Default path through `.dae_default_base()`, explicit path through
  `.validate_base_tbl()`; both end in `extract_allele_freq()`.
- `.dae_default_base(pop, line_name)` — the §3.3 resolver verbatim: named pool
  → shared pool → `.founder_base_empty_error()`. Wahlund warning only on the
  population-wide branch, counting `COALESCE(line_name, '')`.
- `.dae_require_base_at(p_base, written_tf, locus_names, trait)` — refuses
  `NA` at any locus the call will write, naming up to five loci. Runs before
  either `rescale_effects_to_target()` call site; per trait in the multi-trait
  path on `qtl_tf_mat[, t] & !is.na(effects_mat[, t])`.

The completion message now reads `base: founder_haplotypes [1 filter]` /
`base: ind_meta` instead of the old enum value.

### `define_genome_effects()`

`base_tbl = NULL` added between `origin` and `require_complete`. Validated
whenever supplied; queried only when an additive/dominance row has a missing
centre, with the column normalised first (the §3.4 amendment). `.ge_build()`
gained `base_freq = NULL`; the fill sits between members assembly and
`.ge_check_member_fields()`, and a temporary `fill_failed` column lets the
validator append *"— and base_tbl has no allele copies at this locus"* to its
existing per-row message. The column is dropped before anything is written.

---

## Deviations from the plan

1. **Base resolution runs after argument validation, not before both paths.**
   The first cut resolved the base once at the top of `define_additive_effects()`,
   which made the Wahlund warning fire *before* the multi-trait mixed-origin
   error — `gate 44` in `test-genome-effects-writer.R` then saw a warning it
   never had. Moved into each path at the point the old
   `compute_base_allele_freq()` call sat, via `.dae_resolve_base()`. Argument
   errors are reported before the base is touched or warned about. The plan's
   §3.3 sketch showed the resolution at the top; this is the correct placement.
2. **`.founder_base_empty_error()` lives in `R/extract_allele_freq.R`**, not
   `define_additive_effects.R`, because the helper needs it and the helper's
   file is the lower layer. Same signature as the plan's v5 fix: `(conn, what)`.
3. **The "Did you call `define_founder_haplotypes()`?" hint stayed**, on
   `.dae_default_base()` rather than `get_table()`, and now suggests
   `base_tbl = get_table(pop, "ind_meta")` as the alternative. The plan left
   this open ("or dropped if `get_table()`'s message is judged sufficient");
   the hint is more useful than the generic table-missing message.
4. **Test 19's "no query" check is a semantic one, not a call counter.** The
   plan suggested mocking `extract_allele_freq()` or counting `DBI` calls. A
   sharper check needs neither: pass a base whose *query would fail* (no copies
   at any locus) with every centre explicit — if the base were queried the call
   would error; it does not. The same base with a missing centre does error.

Everything else — the three shapes, the default precedence, the warning
boundary, the `NA`/`0` contract, the QTL-coverage error, the fill rules, the
same-connection check, the pool-level (never per-locus) fallback, the vignette
migration rule — is exactly as the v5 plan specifies.

---

## Tests

| Plan test | Where | Status |
|---|---|---|
| 1–9 helper contract | `test-extract_allele_freq.R` (12 blocks incl. the NULL-pool diagnostic and the founder/copies pair) | pass |
| 10–15 default resolution, warning boundary, validation, QTL gap, seed | `test-define_additive_effects.R` | pass |
| 16 crossbreeding end to end (common + two line defaults, F1 + backcross, independent per-copy TBV) | `test-add_tbv.R` "crossbreeding end to end with default bases" | pass |
| 17 `line_origin` vs `line_name` on a fixed F1 fixture | `test-extract_allele_freq.R` (last block) | pass |
| 18–21 writer fill | `test-genome-effects-writer.R` | pass |
| 22 generator ≡ writer, centres omitted, common and composed scope | `test-genome-effects-writer.R` "generator == writer" | pass |
| swine vignette | run through the last `define_additive_effects()` call (line 1520) on the full 10k-locus / 2000-founder config with `load_all()`; all six migrated calls and the WWD/WWM correlated call ran with `base: ind_meta`; exit 0 | pass |

Migrated and green: `test-add_tbv.R`, `test-genome-effects-eval.R`,
`test-genome-effects-writer.R` (`gew_base_A()` helper replaces ten
`base_line_name = "A"`), `test-define_additive_effects.R`. Twelve
`suppressWarnings()` wrappers around pooled-default calls were replaced with an
explicit `base_tbl = get_table(pop, "founder_haplotypes")` — the Q11 "tell"
is gone from the suite.

Full suite at the docs commit: **2864 pass, 0 fail, 9 warnings, 1 skip** — the
warnings are pre-existing (`formula_tbv` scalar-constant notices and the like) in
files unrelated to this change. One incidental pooled default in
`test-print-pop.R` was made an explicit `base_tbl` in the docs commit.

---

## Call-site inventory, as executed

`grep -rn 'base = "\|base_line_name\|compute_base_allele_freq' R/ tests/testthat/ vignettes/`
returns nothing after commit 2.

| File | Before | After |
|---|---|---|
| `vignettes/swine/…sex-semen.R` | 6 × `base = "current_pop"` + 1 comment | `base_tbl = get_table(pop, "ind_meta")` × 6; comment rewritten |
| `tests/testthat/test-add_tbv.R` | 2 × `base = "current_pop", base_tbl = x`; 1 × `base = "current_pop"`; 3 × `suppressWarnings(default)`; 1 test title naming the deleted function | `base_tbl = x`; `base_tbl = get_table(pop, "ind_meta")`; explicit pooled base; title reworded |
| `tests/testthat/test-define_additive_effects.R` | `base = "founder_haplotypes"`, `base = "current_pop", base_tbl =`, whole `base_line_name` section | dropped / `base_tbl =` / section rewritten as tests 10–15 |
| `tests/testthat/test-genome-effects-writer.R` | 10 × `base_line_name = "A"`; 1 × `tidybreed:::compute_base_allele_freq()` | `base_tbl = gew_base_A(pop)`; `extract_allele_freq(gew_base_A(pop))` |
| `tests/testthat/test-genome-effects-eval.R` | 1 × `base = "founder_haplotypes"`; 9 × `suppressWarnings(default)` | explicit pooled base |
| `R/define_founder_haplotypes.R` roxygen | pointed multi-line users at `base = "current_pop"` | points at the per-line default and `line_origin` |
| `R/define_additive_effects.R` roxygen | `base` paragraph, two `@param`s, three examples | rewritten; family sentence added |

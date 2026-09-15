# Genome Effects v4.9 — overall review after Phases A–E

**Reviewed:** 2026-09-13, on `feat/genome-effects-v49` at `c165e49` (0.68.0).
**Against:** `plans/update_genome_effects_v4.md` (v4.9) and the five phase
notes `update_genome_effects_phase_A.md` … `_phase_E.md`.
**Verdict:** the plan is complete. Every phase row, every gate (1–54), and every
★-marked correction in the plan has a corresponding implementation and test.
The suite is green (see §Verification). What follows is (1) what was fixed
during this review because it was plainly wrong, and (2) the items that are
*recorded but not scheduled* — things you may want to change or do better, none
of which blocks a merge.

---

## 1. Fixed in this review (obvious errors)

All small; all documentation or test-adjacent, none a behaviour change.

| # | Where | What was wrong | Fix |
|---|---|---|---|
| 1 | `R/schema.R:117,134,136` | Three column descriptions for `genome_effect_members.id_genome_effect` and `genome_effect_member_origins.{id_genome_effect, member_slot}` still said "FK to genome_effects" / "FK to genome_effect_members". Phase C dropped both keys; `describe_table()` was teaching a constraint that does not exist | Reworded: "matching … (R-enforced by `validate_genome_effects()`; no SQL FK)" |
| 2 | `plans/update_genome_effects_v4.md` §Schema DDL | The DDL block still declared `FOREIGN KEY (id_genome_effect) REFERENCES genome_effects` and the composite origin → member key — the two keys the plan's own header says Phase C removed. §Validation's "Declared in SQL" list named them too | DDL block and the list now carry the ★ Phase C note instead of the keys |
| 3 | `plans/update_genome_effects_v4.md` §Output contract | A stray edit left "`archive_replicate()` Written by **`add_tgv()`**, beside…" as one garbled sentence | Split back into two sentences |
| 4 | `CLAUDE.md` §Schema Design Bias | Listed `genome_effect_type` — deleted in 0.65.0 — as an example of a preferred long-table dimension | Now `contrast_name`, `component_name` |
| 5 | `vignettes/swine/swine-time-based-age-at-puberty-sex-semen.R:910` | `get_table("genome_effects") \|> count(locus_name)` — an executable line that errors, because the term table has no `locus_name`. Phase C fixed the two `expressed_parent` calls in this file (1406, 1439) and Phase E's grep for legacy *names* could not catch this one, since `genome_effects` is still a live table name | Reads `genome_effect_loci` |
| 6 | `README.md` §Database Tables | `genome_effects` described as "1 per (locus × trait × effect type) — Additive QTL effect sizes"; no row for the two child tables, the views, or `ind_tgv` | Five rows describing the term/member/origin shape, the two views, and `ind_tgv` |
| 7 | `R/genome_effects_eval.R`, `R/genome_effects_helpers.R` | Six `[.internal_fn()]` roxygen links inside `@noRd` blocks pointing at other undocumented internals — six "Could not resolve link" warnings on every `roxygenise()`, which is how a future real warning gets ignored | Code spans instead of links; `roxygenise()` is now silent and `man/` unchanged |

---

## 2. Phase-by-phase verification

Checked by reading the code, not the notes.

| Phase | Plan row | Verified |
|---|---|---|
| **A** | 19 fixtures, two evaluators, no schema change | `helper-genome-effects.R` (`gefx_*`), `test-genome-effects-fixtures.R`. The `terms`-format deviation was carried to Phase C as the round-trip test (gate 53, `test-genome-effects-writer.R:194`) |
| **B** | `genome_meta` PK; `open_pop.R` DDL deleted; three tables in `define_genome()`; `GENOME_TABLES` = 10 + `GENOME_VIEWS`; `ind_tgv` in `define_trait()`'s lazy block; 3 views; registries | `define_genome.R:4-12, 296, 455-546`; `open_pop.R:162` (comment only — no DDL); `define_trait.R:230-238, 357`; `sql_utils.R:85-115, 153-156, 183, 211-235, 294-310`; `schema.R:100-190`; `archive_replicate.R:110-114`. Both Phase B plan corrections (no `replicate` on `ind_tgv`; `TABLE_NO_ROW_DELETE` not `TABLE_ROW_KEYS`) are what is built |
| **C** | `define_genome_effects()` with the long `terms` format + four modes + `origin` + `require_complete` + `allow_reserved_owner`; `ad_terms()`, `genotype_terms()`; `define_additive_effects()` on the writer with `parent_origin`, `replace_scope`, origin-aware `scale_to_target`, parent-only warning; `expressed_parent` gone at all 14 sites; intra-set FKs dropped | Signatures match the plan exactly (`define_genome_effects.R:133-143`, `genome_effect_terms_builders.R:47, 190`, `define_additive_effects.R:178-190`). `expressed_parent` survives only as absence-assertions and "what `parent_origin` replaced" prose. `validate_sql_identifier(effect_owner)` at `:150`. `.ge_duplicate_hint()` present (Q2) |
| **D** | One evaluator; `add_tgv()`; `add_tbv()` as a filtered call; per-**member** freedom (`"*"` sentinel); preflight with `tidybreed.label_vector_warn/max`; zero-copy state materialized against genotype-member loci; `.gev_warn_tbv_stale()` (Q1); benchmark script | `genome_effects_eval.R` (all `.gev_*` listed), `add_tgv.R`, `add_tbv.R:141-173` reads `.gev_read_model` → `.gev_reserved_additive` → `.gev_evaluate` — no second implementation. `dev/benchmarks/benchmark_tgv_scale.R` exists |
| **E** | `restore_pop()` guard (both signals + the 0.68.0 rename, connection released); `extract_genotypes()` on `genome_effect_loci` (nominal check **and** locus resolution); `genome_effect_loci.genome_value`; `n_qtl` → `n_causal`; `phenotype_components.genome_effect_types` → `component_names` | `restore_pop.R:111-139`; `extract_genotypes.R:126-131, 210-221`; `genome_effects_helpers.R:69-79`; `tidybreed_pop.R:163-170`; `open_pop.R:325` |

**Gate coverage.** Every gate 1–54 has a test. The ones the Phase D table attributes to
"✅ Phase C" without a `gate NN` label in the file are covered under descriptive
names: 11 → `writer.R:314`, 12 → `:314, :880`, 13 → `:553`, 19 → `:464`, 20 → `:573`,
24 → `schema.R:121` + `writer.R:837`, 36 → `writer.R:753`, 37 → `writer.R:927`.

**Docs.** `NAMESPACE` exports `define_genome_effects`, `ad_terms`, `genotype_terms`,
`add_tgv`; all six `man/` pages present; `roxygenise()` produces no diff.
`NEWS.md` has 0.65.0 → 0.68.0 entries, one per phase. `CLAUDE.md` schema
sections match the DDL (spot-checked all three effect tables, `ind_tgv`, the
views, `chr_*`, `phenotype_components.component_names`).

---

## 3. Recorded, not scheduled — for your review

These are the items the phase notes and the plan explicitly *record without
fixing*. I have not changed any of them. Each ends with a **Recommendation**
naming the fix I would make, where it goes, and the test that pins it. Ordered
roughly by how likely they are to bite a user.

| # | Item | Recommendation | Size | When |
|---|---|---|---|---|
| 3.1 | Mismatched-centre trap | Write-time warning in `define_genome_effects()` | ~30 lines + 2 tests | **Before merge** |
| 3.2 | Realized variance / orthogonality diagnostic | New `extract_genetic_variance()`; schedule it **ahead of** consolidation | New file, ~200 lines | Next plan |
| 3.3 | Dead generic variant under two reciprocals | Write-time warning from the variant map | ~25 lines + 1 test | Before merge, same commit as 3.1 |
| 3.4 | `require_complete` is scope-blind | Leave; document the limitation in the roxygen | 3 lines of roxygen | Now |
| 3.5 | Ploidy proof is line-agnostic | Leave; one comment at the call site | 2 lines | Now |
| 3.6 | `remove_rows()` guard ordering | Move the `TABLE_NO_ROW_DELETE` check to the top | ~10 lines + 1 test | Before merge |
| 3.7 | `ind_tbv` ↔ `ind_tgv` redundancy | Leave; it is the consolidation plan's whole job | — | Next plan |
| 3.8 | `package_summary.md` untracked and stale | Delete it | `git rm` | Before merge |

### 3.1 The mismatched-centre trap (Phase D §Tracked for later; plan Q1 last row)

`define_additive_effects(base = "current_pop")` centres the additive term at the
realized frequency; a hand-written `ad_terms(p = 0.30)` dominance term centres
wherever the user typed. Both writes are legal (different families, they sum),
and the additive coefficient quietly stops being an average effect. Today the
only net is `add_tbv()`'s warning — right check, wrong moment.

Three candidate fixes recorded in the Phase D notes, none built:
1. Write-time warning in `define_genome_effects()` when a `dominance` member
   lands at a locus whose stored additive variant has a different
   `center_value`.
2. Realized-orthogonality reporting (needs 3.2).
3. `ad_terms()` inheriting `p` from the stored additive term when `p` is omitted.

**Recommendation — do (1) now, before merge; (2) arrives with 3.2; skip (3).**

- *Where:* `define_genome_effects.R:161-164` already has both halves in hand —
  `built` (the candidate rows) and `model <- .ge_read_model(conn)` (the stored
  rows). Add one helper, `.ge_warn_centre_mismatch(built, model, labels)`, called
  right after `.ge_resolve_deletes()` and before `.ge_commit()`. Pattern to copy:
  `.dae_warn_parent_only()` in `define_additive_effects.R:338`, which is the same
  shape (a write that is legal but confusable, reported with both values named).
- *Rule:* for every candidate **order-one `dominance`** member at `locus_id` L,
  look up stored **order-one `additive`** members at L for the same trait, any
  owner, whose origin scope matches or contains the dominance member's scope
  (common scope always matches). If any has `center_value` differing by more than
  `1e-8`, warn once per locus: *"dominance term_id 'X' at Locus_10 is centred at
  0.30 but the additive term (owner 'generated_additive_tbv') is centred at
  0.3333; the Cockerham contrast is orthogonal only at the additive term's centre,
  so `tbv_value` for 'ADG' will stop being the average effect. Pass p = 0.3333, or
  re-run define_additive_effects() with the same p."* Check the reverse direction
  too (an additive write landing on a stored dominance term), same message with
  the roles swapped.
- *Why a warning and not an error:* different centres across **scope variants**
  are legal and intended — a line-A variant centres at line A's frequency — which
  is exactly why `center_value` is outside the family signature. The check has to
  be scope-aware (containment, not locus-only) or it fires falsely on the
  common/A/B pattern of gate 35.
- *Why not (3):* `ad_terms()` is a pure builder with no `pop`; giving it a
  database read changes its contract for one convenience, and it does not help
  the `data.frame`-typed path at all.
- *Tests:* one asserting the warning on the Phase D 4-locus reproduction
  (`base = "current_pop"` then `ad_terms(p = 0.30, coding = "cockerham")`), one
  asserting silence on gate 35's common/A/B sequence and on a Cockerham write
  with the matching `p`.

### 3.2 Realized variance / orthogonality diagnostic (§Future limitations 6)

The plan calls this "the first thing to build after this plan". Storage for
dominance and epistasis ships with no path from a user's target `V_D` to
coefficients and no way to read back what variance a stored model has. It is
also the only *true* (population-measured, not `center_value`-declared) check
for 3.1.

**Recommendation — schedule it as the very next plan, ahead of
`consolidate_genetic_values.md`.** It has no schema dependency on consolidation,
it turns 3.1's warning into a measurement, and it is the thing that makes
dominance/epistasis storage *usable* rather than merely storable.

- *Shape:* `extract_genetic_variance(tbl, trait_name = NULL, by = NULL)` —
  `extract_` prefix per the naming table (returns analysis data, changes no
  state); takes a `tidybreed_table` of individuals like `add_tgv()`. Returns a
  tibble `(trait_name, component_name, var_value, n_ind)` from the realized
  `ind_tgv` rows, plus an `orthogonality` attribute or second tibble giving
  `E[x_A x_D]` per locus where both contrasts exist. Base-population variance
  from stored effects + founder frequencies (the LE-assuming formula the plan's
  §Out of scope names) is a second function or a `method =` switch, not the
  default, because the realized number is the honest one.
- *Test:* on a purebred Cockerham model in HWE, `V_A` from the function agrees
  with `sum(2pq a²)` to sampling error and `E[x_A x_D] ≈ 0`; on the 3.1
  reproduction it is measurably non-zero.

### 3.3 Both reciprocals defined ⇒ the generic `{A,B}` variant is unreachable (Phase A finding 3)

Two reciprocal dominance variants `{A:1@p1, B:1@p2}` and `{A:1@p2, B:1@p1}`
partition the `{A,B}` label space, so a generic `{A:1, B:1}` fallback beneath
them can never be selected. Correct by the rules, but a user who writes all
three has written one term that does nothing, silently. Phase A named it as a
candidate for "the same class of warning as gate 50". Nothing was built and
nothing documents it in `define_genome_effects()`'s roxygen.

**Recommendation — write-time warning, same commit as 3.1.** The machinery
already exists: `.gev_variant_map()` (`genome_effects_eval.R:369`) resolves every
`(family, label-vector)` to a winning `id_genome_effect`. After
`.ge_resolve_deletes()` in the writer, build the post-write model (stored minus
`drop` plus `built`), run the map over the label alphabet, and any term in a
scoped family that wins **no** label-vector is dead. Warn once, naming the user's
`term_id` and the sibling variants that shadow it. Restrict to families with ≥ 2
variants so the common single-variant case costs nothing; the label alphabet
comes from `ind_haplotype`, so on an empty population (effects defined before
founders) skip the check rather than warn on an empty alphabet.

- *Also:* add one sentence to the `@details` of `define_genome_effects()` under
  the reciprocal example: *"If both reciprocals are written, a generic `{A, B}`
  variant beneath them can never be selected; the writer warns."*
- *Test:* Phase A's F04 fixture with the second reciprocal added — warns and
  names the generic term; F04 as written (one reciprocal) stays silent.

### 3.4 `require_complete` is scope-blind (Phase C correction 4; plan §Writer API)

Coverage is checked over the terms *in the call*, grouped by
`(trait, owner, ordered locus set)`. A surface completed across two calls at
two scopes is not checked as a whole.

**Recommendation — leave the behaviour; state it in the roxygen.** No gate
covers it, the meaning of a per-scope grid under partial containment is
genuinely unsettled, and nobody has written a scoped indicator surface. Add to
`@param require_complete`: *"Coverage is checked over the terms in this call
only, and is scope-blind: a surface written at one origin scope is not checked
against the same surface at another."* Three lines, no code.

### 3.5 Dominance ploidy proof is line-agnostic (Phase B finding 6)

`.ge_validate_dominance_ploidy()` and the writer both call
`resolve_chr_inheritance(conn, sex)` with no `line_name`
(`genome_effects_helpers.R:349`, `define_genome_effects.R:367`), so a
line-specific karyotype rule in `chr_inheritance` would not be consulted. No
such rule can be written today — `define_chromosome()` has no `line_name`
argument — so this is a latent coupling, not a bug.

**Recommendation — leave it; add one comment at each of the two call sites:**
`# line_name = NULL: line-specific chr_inheritance rules do not exist yet;
revisit here if define_chromosome() gains line_name.` The right time to change
the code is when `define_chromosome(line_name = )` lands, and the comment is what
makes that change find its way here.

### 3.6 `remove_rows()` guard ordering (Phase B findings 4 and 9)

`remove_rows()` checks "no filter applied" (`remove_rows.R:200`) and "filter
matched 0 rows" (`:218`) *before* `TABLE_NO_ROW_DELETE` (`:228`), so an
unfiltered `get_table("genome_effects") |> remove_rows()` — or one with
`confirm_all = TRUE` — reports the generic "would delete ALL rows" message, and
the `confirm_all = TRUE` form even suggests a command that will then fail for a
different reason. Either message stops the user, so nothing is deleted.

**Recommendation — move the `TABLE_NO_ROW_DELETE` check above guard B, before
merge.** It is a table-level property and belongs before any row-level
reasoning: a table that refuses row deletion refuses it regardless of the filter.
Ten lines: hoist the `if (table_name %in% names(TABLE_NO_ROW_DELETE)) stop(...)`
block to immediately after `table_name` is known, and delete it from section D.
Behaviour change is confined to the *message* on refused tables (views, the three
effect tables, and whatever else is registered); no table that could be deleted
from before becomes undeletable.

- *Test:* unfiltered and `confirm_all = TRUE` calls on `genome_effects` and on
  `genome_effect_loci` each error with the `define_genome_effects` /
  derived-view reason, not the "no filter" text. Phase B's existing filtered-call
  tests (`test-genome-effects-schema.R:182-190`) stay as they are.

### 3.7 `ind_tbv` ↔ `ind_tgv` redundancy (plan §Known, time-boxed redundancy)

`ind_tbv.tbv_value` = the reserved owner's breeding value;
`ind_tgv.'order1_additive'` = *every* additive-structured term including
`custom` ones. Identical unless someone hand-writes an order-one additive term
under a non-reserved owner, at which point they differ and `add_phenotype()`
reads the former.

**Recommendation — leave it; this is the consolidation plan's entire job.**
`add_tbv()` already warns on exactly the case where they diverge
(`.gev_warn_tbv_stale()`, first row of the Q1 table), so the trap is not silent
in the meantime. Do not patch `add_phenotype()` piecemeal — Task 3 of
`plans/consolidate_genetic_values.md` moves it onto `ind_tgv_total` in one step.

### 3.8 Two cosmetic observations

- **Plan audit item 3** still reads "composite FOREIGN KEY … Deletes do not
  cascade, so replacement must delete children before parents." That was true of
  DuckDB when written and is now moot (no intra-set FKs).
  **Recommendation — leave as-is.** It is the historical audit against
  `4e9a72b`, the ★ correction at the top explains the change, and the DDL block
  below it is now consistent (fix #2 above). Rewriting audit history to match
  the outcome is the thing the plan's own "corrected in place, ★-marked"
  convention avoids.
- **`package_summary.md`** in the repo root is untracked, stamped 0.65.0 /
  2026-09-10, and describes the schema as of Phase B.
  **Recommendation — `git rm`-equivalent: delete it before merge.** It
  duplicates `CLAUDE.md` (which *is* the maintained package summary) at a stale
  version, and anything that lands on `main` untracked will be re-created stale
  again. If a one-page external summary is wanted, generate it from `CLAUDE.md`
  at release time under `dev/`, not the root.

---

## 4. Verification

Run on 2026-09-13 after the fixes in §1:

- `roxygen2::roxygenise(".")` — **0 warnings**, no diff in `man/` or `NAMESPACE`.
- Legacy-name grep (`genome_effect_type\b|genome_effect_types|base_allele_freq\b|expressed_parent|effect_set_name`) over `R/ tests/ vignettes/ man/ README.md CLAUDE.md`, minus `compute_base_allele_freq()`: every remaining hit is one of the deliberate categories Phase E itemized (absence-assertions, "what `parent_origin` replaced" prose, the `restore_pop()` guard's own message, fixtures that build the old shape to be refused).
- Full suite — see the result appended below.

```
testthat::test_dir("tests/testthat")  (pkgload::load_all, SilentReporter, full tally)
RESULT: pass 2772  fail 0  error 0  warn 10  skip 1
```

Identical to the Phase E figure. The 10 warnings are the pre-existing set Phase D
and E itemized (Wahlund, monomorphic map loci, scalar `formula_tbv`, `hap_id`
tibble); the skip is the `add_founders` large-line performance guard.

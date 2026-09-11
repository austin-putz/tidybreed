# Phase C — the writers

**Status: complete.** Branch `feat/genome-effects-v49`, version `0.66.0`.
Companion to `plans/update_genome_effects_v4.md` (v4.9), which was corrected in
place where this phase found it wrong — every such correction is flagged below
and marked ★ in the plan itself.

Phase C builds everything that *writes* genome effects. Nothing here evaluates
one; the evaluator is Phase D, and `add_tbv()` therefore still reads a column
shape that no longer exists.

---

## What shipped

| File | What |
|---|---|
| `R/define_genome_effects.R` (new, ~860 lines) | The general writer, the `terms`/`origin` parser, copy-count inference, `require_complete`, the four replacement modes, and the single-transaction commit |
| `R/genome_effect_terms_builders.R` (new, ~250 lines) | `ad_terms()` and `genotype_terms()` — the two `terms` constructors |
| `R/define_additive_effects.R` (rewritten) | Rebuilt on the general writer; new `parent_origin`; origin-aware `scale_to_target`; the parent-only re-run warning; `method = "union"` reads the new tables |
| `R/genome_effects_helpers.R` | `labels` on both validators, `.ge_scope_label()`, `.ge_duplicate_hint()` |
| `R/define_genome.R` | The two intra-set foreign keys removed, with the reason in the DDL comment |
| `R/define_trait.R`, `R/add_tbv.R`, `R/schema.R`, `R/sql_utils.R` | `expressed_parent` deleted |
| `tests/testthat/test-genome-effects-writer.R` (new) | 183 assertions |
| `tests/testthat/test-define_additive_effects.R` | Retargeted onto a test-local flat view (see **Verification**) |
| `tests/testthat/test-genome-effects-schema.R` | FK tests retargeted; the DuckDB delete-in-transaction behaviour pinned |
| `tests/testthat/test-add_tbv.R`, `helper-parity.R`, `test-define_trait.R` | Imprinting retargeted onto origin rows (gate 40) |

Three exported functions are new: `define_genome_effects()`, `ad_terms()`,
`genotype_terms()`.

---

## The `terms` format, as built

A long data frame, one row per (term × locus). Accepted columns are exactly
`term_id`, `genome_value`, `effect_name`, `locus_name`, `contrast_name`,
`center_value`, `copy_count_value`, `dosage_value`; anything else is **rejected,
never silently dropped**. `term_id` is user-facing only and is replaced by
`next_int_id()`; a single-term call may omit it.

Scope lives in `origin`: `NULL`, a named scalar list, or a data frame keyed by
`locus_name`. The scalar list produces one origin row per member — all an
`additive` member may have, and usually *not* enough for a genotype member,
whose multiset must sum to the state's copy count. The validator says which,
naming the locus.

**Every message about malformed input names the `term_id` the user typed.** This
took a small change to Phase B's validators: both now take an optional
`labels` map from `id_genome_effect` to a display label. Post-insert whole-set
validation passes it too, so a conflict between a new term and a stored one
names the new one by the user's label and the stored one by its real id — which
is correct, because the stored one *is* visible to the user, in
`genome_effect_terms`.

### What the format buys, asserted

- constant `genome_value` and `effect_name` within a `term_id`
- each locus at most once per term
- `contrast_name` ↔ state/centre agreement, checked before any SQL
- origin row-local rules (`'exact'` needs a line, `'any'` needs a parent, …)
  checked in R **as well as** by the SQL `CHECK`, so the user sees a sentence
  about their own call rather than a DuckDB expression
- silent coercions refused: a character `genome_value` (which `as.numeric()`
  would turn into `NA`) and a fractional `dosage_value` (which `as.integer()`
  would truncate into a *different genotype state*)

---

## Corrections to `update_genome_effects_v4.md`

### 1. The intra-set foreign keys had to go — a plan error, and a Phase B one

The plan specified, and Phase B declared and verified:

```
genome_effect_members.id_genome_effect          -> genome_effects
genome_effect_member_origins.(id, member_slot)  -> genome_effect_members
```

Phase B probed these against real **inserts**. It did not probe a **delete
inside an explicit transaction**, which is what every `replace_*` mode does.
DuckDB 1.5.5 refuses it:

```
Constraint Error: Violates foreign key constraint because key "id: 1"
is still referenced by a foreign key in a different table.
```

Probed to exhaustion before changing anything (`fk_probe*.R` in the scratchpad):

| Case | Result |
|---|---|
| single-column FK, delete child then parent, **in a transaction** | FAIL |
| single-column FK, same sequence, **autocommit** | OK |
| single-column FK, delete parent first | FAIL |
| composite FK, delete child then parent, in a transaction | FAIL |
| composite FK, same sequence, autocommit | OK |
| composite FK, whole-table `DELETE` on the child first | FAIL |
| insert then delete the same row in one transaction | OK |

The FK index retains entries for rows deleted in the current uncommitted
transaction. It is not an ordering problem and not a composite-key problem; it
is a transaction-visibility limitation, and autocommit is the only escape.

Autocommit is not an option here. Deleting the old rows outside the transaction
that inserts the new ones means a failed insert leaves the trait with no effects
at all — and a half-replaced effect model is not a weaker version of the
requested model, it is a **different genetic model**, silently.

**Resolution.** The two keys are dropped. The integrity they would buy is
enforced by `validate_genome_effects()`, which already reported orphans in both
directions and runs inside every write transaction before `COMMIT`. This is a
real weakening in exactly one respect: a raw `DBI::dbExecute()` can now create an
orphan that survives until the next write. That route bypasses every other guard
in the package too — `remove_rows()` refuses these tables, every column is
reserved, `mutate_table()` is blocked, and `define_genome_effects()` is the only
writer.

The `locus_id → genome_meta` key **stays**. `genome_meta` rows are never deleted,
so it never sits in the failing position, and it is the one relationship R cannot
cheaply re-derive. It is also the FK gate 47 tests.

`test-genome-effects-schema.R` gained two tests: one asserting the R-level orphan
detection in all three directions, and one **pinning the DuckDB behaviour itself**
so a version that fixes it is noticed rather than leaving the workaround standing
forever.

### 2. Gate 42 was unsatisfiable as written

> All four `(line_name, parent_origin)` combinations round-trip … and
> `replace_scope` on one leaves the other three untouched.

Two of the four **cannot coexist in one fallback family**. `('exact' A, parent
ANY)` and `('any', parent 1)` overlap without either containing the other — a
line-A paternal copy matches both — so rule 3 of §Origin resolution rejects the
pair at write time. Correctly: no variant could be selected for that copy.

This is a property of the containment lattice, not a defect, and the plan
describes it two sections above the gate. The gate is met with **one locus per
scope**, which puts the four in four families. `replace_scope` keys on
`(trait, owner, scope)` and ignores locus, so its isolation is still exactly what
is being tested. A second test asserts the rejection in **both** orders, so the
incompatibility is pinned rather than merely worked around.

### 3. Q2's combined-equivalent message was specified but never scheduled

The plan says a duplicate caused only by a different `center_value` should be
rejected **with the combined equivalent in the message**, and Q2 recommends it —
but no gate covers it, so it would have shipped as a bare duplicate error.
Implemented as `.ge_duplicate_hint()`:
`a₁(g − c₁) + a₂(g − c₂) = (a₁ + a₂)(g − c′)`, `c′ = (a₁c₁ + a₂c₂)/(a₁ + a₂)`.
Restricted to order-one `additive` — the Cockerham dominance contrast is not
linear in its centre, so no such combination exists for it — and the cancelling
case (`a₁ + a₂ = 0`) reports that the combination is a constant, which this model
has no place for.

### 4. `require_complete` is scope-blind

Coverage is checked over the terms in the call, grouped by
`(trait, owner, ordered locus set)`. A surface written at one origin scope is not
checked against a surface at another. Completing a scoped surface is a coherent
thing to want; nothing prevents it, it simply is not verified. Recorded in the
plan rather than fixed: no gate covers it, and the meaning of a per-scope grid
under partial containment is not settled.

---

## Implementation issues found in review

1. **`'any'` reached SQL instead of R.** A scope that constrains nothing hit the
   row-local `CHECK` mid-transaction, which aborts the transaction and reports a
   DuckDB expression. All origin row-local rules are now checked in R first; the
   `CHECK` remains as the backstop that makes the rule true of the *table*, not
   just of the writer.

2. **Named `parent_origin` lost its names.** `as.integer(c(ADG = 1L))` drops
   `names()`, so `v[[nm]]` threw `subscript out of bounds` — the named form,
   which is the one that makes per-trait origins ergonomic, was entirely broken.
   Caught by the gate 44 test, not by any earlier smoke run.

3. **`.ad_rbind_fill()` failed on a zero-row frame.** `d[[col]] <- NA` on a
   0-row data frame errors, so `ad_terms(a = 0.4, d = 0)` — the single most
   likely call — died. `rep(NA, nrow(d))`.

4. **`ge_seed_term()` was borrowed across test files.** The new test file called
   a function defined in `test-genome-effects-schema.R`. It works in a full-suite
   run (alphabetical order) and fails on `test_file()`. Inlined.

5. **A fixture was mis-classified in the round-trip.** F14 (dominance at a
   hemizygous locus) carries `expect_eval_error`, not `expect_invalid`, because
   Phase A had no `chr_inheritance` table and could only pin it at the evaluator.
   With a database it is refused at the **writer**, which is where the fixture's
   own note said it belonged. The round-trip now excludes it and asserts the
   write-time refusal explicitly — Phase C strengthens a Phase A fixture rather
   than skipping it.

6. **`add_tbv()` was left with dead scaffolding.** Deleting the
   `expressed_parent` branch left `parent_filter <- ""` as a constant, a
   `meta_rows` reorder with no reader, and an `m <- meta_rows[...]` row lookup
   used by nothing. All removed. Pre-1.0 policy: no vestiges.

7. **`GE_RESERVED_OWNERS` repeated `GE_ADDITIVE_OWNER`'s literal.** Two
   constants holding the same string is a silent-drift hazard; the reserved list
   is now derived from the owner constant.

8. **Silent coercions would have stored plausible wrong values.** A character
   `genome_value` became `NA` under `as.numeric()`; a fractional `dosage_value`
   became a **different genotype state** under `as.integer()`. Both now error.
   Neither was in any gate; both were found by reading the parser back against
   what R's coercion rules actually do.

9. **`origin_slot` did not sort NULLs last as specified.** The plan pins the
   canonical order as `(line_match_type, line_name, parent_origin, copy_count)`,
   **NULLs last**; the first implementation used `""` for a missing line, which
   sorts first. Unreachable as a tie-breaker — a NULL line name occurs only under
   `unknown`/`any`, which already differ on the primary key — but the code should
   say what the rule is, not what the current lattice happens to allow.

---

## Deleting `trait_meta.expressed_parent`

All fourteen sites, per the plan's inventory:

| Site | Done |
|---|---|
| `define_trait.R` — DDL, argument, `match.arg`, written value | ✅ |
| `add_tbv.R` — the `SELECT`, the `parent_filter` branch, both prose blocks | ✅ |
| `sql_utils.R` `TABLE_RESERVED_COLS`, `schema.R` column description | ✅ |
| `man/define_trait.Rd`, `man/add_tbv.Rd` | ✅ regenerated |
| `helper-parity.R`, `test-add_tbv.R`, `test-define_trait.R` | ✅ retargeted |
| `CLAUDE.md` — four sites | ✅ |
| `vignettes/tidybreed-introduction.Rmd:112` | ✅ |
| `vignettes/swine/…-sex-semen.R:1406, 1439` — two **executable** calls | ✅ |

`test-define_trait.R`'s test was not deleted but **inverted**: it now asserts the
column and the argument are both absent, so a re-introduction is caught rather
than silently tolerated.

`test-add_tbv.R`'s imprinting fixture (gate 40) keeps its full shape — Duroc and
Landrace variants with distinct coefficients and per-line centres, an F1, paternal
expression — and now asserts the stored scope explicitly: two
`('exact', line, parent 1, copy_count 1)` rows, one per line. That assertion is
the point. A population-wide `('any', 1)`-only fixture would pass with the line
dimension silently dropped, which is exactly the bug that stamping `'any'`
unconditionally would introduce: both line-specific calls would land on one
identical scope and the second would replace the first.

The one surviving mention of the old name is in `define_additive_effects()`'s
`@param parent_origin`, which says what it replaces. That is migration
orientation for the one person migrating, not documentation that teaches the old
name.

---

## Gate status

| Gate | Status |
|---|---|
| 34 — general-writer defaults survive rerunning the generator | ✅ three reruns; `custom` and `my_model` untouched, values unchanged |
| 35 — successive common / A / B preserve every variant | ✅ three variants, **one** family; re-running the common call drops loci absent from the new set |
| 41 — `'any'` constraints | ✅ no-parent `'any'` refused by the writer **and** by direct insert; `'any'` on a genotype member refused in R |
| 42 — wrapper scope matrix | ✅ restated (see above); all four mappings round-trip, `replace_scope` isolates, and the incomparable pair is refused in both orders |
| 43 — origin-aware variance target | ✅ unparented and parent-qualified models both realize `V = 1.0` to 1e-8; the latter would have landed at 0.5 before |
| 44 — multi-trait origin contract | ✅ scalar, positional-vector and named forms; mixed + `G` rejected with the covariance reason, nothing written |
| 50 — parent-only re-run | ✅ warns and leaves both variants; common/A/B and a reciprocal `(exact A, 1)`/`(exact A, 2)` pair both silent |
| 53 — `terms` round-trips | ✅ all three worked examples; every valid Phase A fixture; malformed calls named by the user's `term_id` |

Gates 1–33, 36–40, 45, 51–52 belong to Phase D (evaluation) and 46–49 to Phases
B and E.

---

## Verification

| File | Result |
|---|---|
| `test-genome-effects-writer.R` (new) | **183 pass, 0 fail, 0 warn** |
| `test-genome-effects-schema.R` (Phase B) | 99 pass, 0 fail — after the FK retarget |
| `test-genome-effects-fixtures.R` (Phase A) | 248 pass, 0 fail — unchanged |
| `test-define_trait.R` | 26 pass, 0 fail |
| `test-define_genome.R` | 80 pass, 0 fail |
| `test-define_additive_effects.R` | 40 pass, **2 errors**, both `add_tbv()` |

### `test-define_additive_effects.R` was retargeted, not deferred

519 lines written against the flat `genome_effects` shape. `add_tbv()` is Phase
D's problem; **this file is Phase C's**, because it tests the function this phase
rewrote. It now installs a test-local view, `gen_add_flat`, reconstructing
`(locus_name, line_name, base_allele_freq, genome_value)` from the term / member
/ origin tables, restricted to the generated owner. The package ships no such
view and should not — the point of the new schema is that a term is not a locus
— but every assertion in that file concerns *generated additive* effects, which
are exactly the order-one, single-origin case the flat shape described correctly.
The assertions and their reasoning (Falconer variance, Wahlund, per-line
centring) are untouched.

The two remaining errors in it are the two tests that call `add_tbv()`.

### Full suite

```
2,194 pass · 0 fail · 144 error · 10 warn
```

Phase B ended at 1,887 / 3 / 180. Every number moved the right way: +307 passes,
failures to zero, errors down 36.

**Every one of the 144 traced to the old table shape.** Verified individually,
not assumed:

| Cause | Count | Owner |
|---|---|---|
| `Binder Error: Referenced column "genome_effect_type" not found` | 136 | Phase D (`add_tbv()`, directly or via `add_phenotype()`) |
| `extract_genotypes()`'s `effects_tbl` path wants `locus_name` | 5 | Phase E (`extract_genotypes.R:118-127`) |
| Test-local SQL against the old flat shape in `test-add_tbv.R` / `test-sql_injection_hardening.R` | 3 | Phase D / E, with their functions |
| `object 'pop_mean' not found` | 1 | An `on.exit` cascade from an `add_phenotype()` failure earlier in the same test |

By file: `test-mutate_derived` 25, `test-formula_phenotype` 17,
`test-add_tbv_index` 17, `test-remove_rows` 16, `test-phenotype_composite` 15,
`test-add_phenotype` 12, `test-schema-registries` 11, `test-add_tbv` 10,
`test-extract_genotypes` 5, `test-define_effect_fixed_cov` 5, and eight files
with one or two each. Each of those files reaches `add_tbv()` through its own
setup; none of them is testing a Phase C surface.

Nothing failed for any reason Phase C introduced, and no test *failed* at all —
every one of the 144 is an error raised by a query against a column that no
longer exists.

The 10 warnings are all pre-existing and unrelated: pooled base allele
frequencies (a deliberate Wahlund warning the tests provoke on purpose),
monomorphic loci in a map test, a scalar constant in a `formula_tbv`, and one
`hap_id` tibble warning in `add_founders`.

**A reporting note for future phases.** The first two full-suite runs of this
phase reported "10 failures" and named a single file. That was
`SummaryReporter`'s default `max_reports = 10` truncating the list, not the
result — `options(testthat.progress.max_fails = ...)` does **not** raise it,
because that option belongs to `ProgressReporter`. Pass a constructed reporter
(`SummaryReporter$new(max_reports = 5000L)`) or the silent reporter and tally the
data frame. The Phase B report's "3 fail, 180 error" came from a run that did not
hit the cap, so the comparison above is sound.

---

## Next

**Phase D** — one evaluator built to §Evaluation strategy (label alphabet,
resolved variant map, member reduction including the synthesized zero-copy state,
family partitioner, containment resolver, label-vector preflight, term
evaluator), `add_tgv()` writing `ind_tgv`, and `add_tbv()` as a **thin filtered
call into that same evaluator** — reserved owner, order-one `additive` only — not
a second implementation. Gates 1–39, 40, 45, 51–52. This is the phase that turns
the current red green.

**Phase E** — delete the last old-shape remnants, the `restore_pop()` guard for
pre-change files, and `extract_genotypes()`'s `effects_tbl` path (`:118-127`,
plus the `genome_effect_type` mention in its roxygen at `:33`). Gate 49.

Merge to `main` when Phase D is green.

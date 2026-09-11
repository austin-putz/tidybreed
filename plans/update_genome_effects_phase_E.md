# Phase E — the last readers, and the end of the old shape

*Companion to `plans/update_genome_effects_v4.md` (v4.9). Phase A:
`update_genome_effects_phase_A.md`. Phase B: `..._phase_B.md`. Phase C:
`..._phase_C.md`. Phase D: `..._phase_D.md`.*

Phases B–D replaced the storage and rebuilt everything that computes from it.
Phase E is the sweep: one reader still assumed one row per locus, one test
helper still aliased a deleted column name back into existence, one column still
named a vocabulary deleted in 0.65.0, and a database written before the change
had no way to say so. After this phase **no legacy genome-effect shape
remains anywhere in the package** — not in `R/`, not in the test helpers, not in
the schema descriptions, not in the vignettes.

The suite goes from **2746 pass / 0 fail / 7 error** (end of Phase D) to
**2772 pass / 0 fail / 0 error**. The 7 errors were one root cause and are gone,
and the failure they were hiding (below, issue 7) is fixed too.

---

## What shipped

| File | Status | What it is |
|---|---|---|
| `R/extract_genotypes.R` | rewritten reader | `effects_tbl` takes `get_table("genome_effect_loci")`; the locus set is read straight off `locus_id` |
| `R/restore_pop.R` | new guard | A pre-0.65.0 or pre-0.68.0 file stops here, with the connection released |
| `R/genome_effects_helpers.R` | view widened | `genome_effect_loci` gains `genome_value` |
| `R/tidybreed_pop.R` | rename | `n_qtl` → `n_causal`; the query already read the view |
| `R/open_pop.R`, `R/define_trait.R`, `R/define_phenotype.R` | column rename | `phenotype_components.genome_effect_types` → `component_names` |
| `R/sql_utils.R`, `R/schema.R` | registries | Reserved columns and descriptions follow both changes |
| `tests/testthat/test-restore_pop.R` | +3 tests | Gate 49, both halves, plus the 0.68.0 rename |
| `tests/testthat/test-genome-effects-schema.R` | +1 test | The coefficient is on every member row, and is not split |
| `tests/testthat/test-extract_genotypes.R` | +1 test, 5 retargeted | The large-effect filter, end to end |
| `tests/testthat/test-sql_injection_hardening.R` | retargeted | `loci_tbl` is now the only name-driven `IN` list |
| `tests/testthat/helper-genome-effects-db.R` | alias dropped | `gen_add_flat` exposes `center_value`, not `base_allele_freq` |
| `tests/testthat/test-genome-effects-writer.R` | join qualified | The two views now share three columns; gate 35's join had to say which one it meant |
| `tests/testthat/parity_golden/tbv.rds` | re-captured | Phase C's origin-aware rescale, finally visible — see issue 6 |
| `vignettes/`, `CLAUDE.md`, `NEWS.md`, `man/` | docs | The old shape is not taught anywhere |

---

## The three moves

### 1. `extract_genotypes()` reads the locus view

The term table has neither `locus_name` nor `locus_id` — a coefficient spans one
or more loci, and which loci is a property of `genome_effect_members`. So the
locus grain exists only in `genome_effect_loci`, and that is what `effects_tbl`
now takes. The function reads `locus_id` off the collected view and unions it
into the locus set; the old `locus_name` → `genome_meta` → `locus_id` round trip
is deleted, not replaced.

The semantics widen in one way worth stating: a **multi-locus term contributes
every one of its loci**. Filtering `trait_name == "ADG"` on a trait with an A×A
term returns both of that term's loci, which is the right answer to "which loci
are causal for ADG" and a change from the old table, where an interaction could
not be expressed at all.

### 2. `restore_pop()` stops a pre-change file — gate 49

```
'.../old.duckdb' carries the pre-v0.65.0 'genome_effects' shape (one row per
locus). Genome effects are now stored as terms: one coefficient in
'genome_effects', its loci in 'genome_effect_members', and their scope in
'genome_effect_member_origins'. There is no in-place migration -- rebuild the
population with define_genome() and re-declare the effects with
define_additive_effects() / define_genome_effects().
```

Detected on two independent signals, because a half-migrated file fails either
way: `genome_effects` missing `effect_owner`, or either child table absent. The
check runs before `new_tidybreed_pop()`, so nothing is constructed from a file
that cannot be read.

★ **The rename in move 3 creates the same hazard one release later**, and gets
the same treatment. A file written by 0.65.0–0.67.0 has the term/member genome
effects but `phenotype_components.genome_effect_types`, and would fail inside
`define_phenotype(components = )` as a DuckDB append error naming a column the
user never typed. Both checks now go through one `stop_stale()` helper, so the
guard reads as what it is — *this file was written by a version whose schema
this code cannot read* — rather than as two unrelated special cases.

### 3. The last column naming a deleted vocabulary

`phenotype_components.genome_effect_types` (default `'additive'`) was written by
`define_phenotype()` and read by nothing. Its name pointed at
`genome_effects.genome_effect_type`, deleted in 0.65.0. It is now
`component_names`, default `'order1_additive'` — the `ind_tgv.component_name`
vocabulary that exists today and that it will actually select from when
non-additive genetic values reach the phenotype layer.

The column stays reserved and stays unread. **Activating** it is explicitly out
of this plan (§Implementation order, "Out of this plan, contract defined here"),
and removing it would throw away a dimension the schema deliberately reserves.
Renaming it is neither of those: it stops the schema teaching a vocabulary that
no longer exists.

---

## Corrections to the plan

### 1. ★ The locus view could not express "large-effect QTL"

The plan's §Views fixes `genome_effect_loci` at seven columns, and
`genome_value` is not one of them. But `extract_genotypes()`'s own documented
example filtered `abs(genome_value) > 0.15`, and `get_table()` reads **exactly
one relation** — there is no join available to a `tidybreed_table`. Moving the
reader onto the view as specified would therefore have silently deleted a
documented workflow: the coefficient is a term attribute, the locus a member
attribute, and before this change no relation carried both.

`genome_value` is added to the view, repeated on every member row of a term.
The description carries the hazard that repetition creates:

> The whole term's coefficient, from `genome_effects`, repeated on every member
> row. Filter on it to select large-effect loci; never `SUM` it, because an
> interaction term would be counted once per member.

`trait_name` and `effect_owner` were already denormalized onto this view for the
same reason, so this is the existing pattern rather than a new one.

### 2. ★ The Phase E row named the nominal check, not the reader

The plan says Phase E moves "`extract_genotypes()`'s nominal `table_name` check
(`:124-127`)". That is one of two sites. The **locus resolution** (`:203-219`)
also read `locus_name` off the collected table — a column `genome_effects` no
longer has. Fixing only the nominal check would have left the function accepting
the right relation and then failing on the wrong column. Both moved.

### 3. ★ The `tidybreed_pop.R:159` item was already done

The plan's Phase E row also lists "move QTL extraction (`tidybreed_pop.R:159`)
to the locus view". The print method has read
`COUNT(DISTINCT locus_id) FROM genome_effect_loci` since the view existed; only
the internal variable name still said `n_qtl`. The row is stale rather than
wrong — recorded so the ★ audit does not read as unfinished work.

---

## Implementation issues found in review

1. **The injection test rewrote a column that no longer exists.** It did
   `UPDATE genome_effects SET locus_name = ...` across three tables to keep the
   denormalized name consistent. Under the term schema that statement errors
   before any assertion runs. Worse, the test's *point* — that `sql_in_list()`
   escapes a quote-embedded locus name — no longer applies to the effects path
   at all, because that path now passes integers. Retargeted to `loci_tbl`,
   which is the only remaining `extract_genotypes()` path that builds an `IN`
   list out of names, with a companion assertion that the effects path returns
   the same columns and so cannot carry the quote into SQL.

2. **Widening the view made a join ambiguous, and a test caught it.** Gate 35
   in `test-genome-effects-writer.R` ran
   `SELECT locus_name, scope_description, genome_value FROM genome_effect_terms
   JOIN genome_effect_loci USING (id_genome_effect)`. `genome_effect_terms`
   already carried `genome_value`; adding it to `genome_effect_loci` made the
   reference ambiguous and the query error. The two views now share three
   columns — `trait_name`, `effect_owner`, `genome_value` — so any join between
   them has to qualify them. Fixed by aliasing both sides. Worth stating rather
   than quietly patching: denormalizing onto a view buys a filter and costs the
   unqualified join, and the existing `trait_name` / `effect_owner` overlap
   meant that price was already being paid.

3. **The refusal has to disconnect before it stops.** `restore_pop()` opens the
   file first; a bare `stop()` would leave a live DuckDB handle on it, and every
   later `dbConnect()` to that path fails with a lock error that has nothing to
   do with the real problem. The existing `ind_meta` guard already gets this
   right and was the model. The test reopens the file after the error to pin it.

4. **`base_allele_freq` survived as a test-only alias.** `ge_flat_view()`
   aliased `m.center_value AS base_allele_freq` so Phase-A-era assertions could
   stay unchanged. That is a compatibility shim in everything but name, and
   pre-1.0 policy has no room for one. Alias dropped, ~20 call sites renamed.
   The one surviving mention is `expect_false("base_allele_freq_ADG" %in%
   genome_cols)`, which asserts the *absence* of the old `genome_meta` column
   and has to name it.

5. **A schema description claimed a join that no longer exists.**
   `ind_haplotype.locus_name` was described as "denormalized for direct joins to
   `genome_effects`". The effect tables key on `locus_id`; the column now
   earns its place as an export/query convenience, and the description says so.
   Same fix in `CLAUDE.md` and in a `test-define_genome.R` comment.

6. **`n_qtl` in the print method.** The query counted causal loci; the variable
   still said QTL, which is precisely the conflation §Causal-locus terminology
   exists to end.

7. **A Phase D failure was masked by a Phase D error.** `test-parity.R`
   errored inside `run_parity_sim()` at the `extract_genotypes()` call, so the
   golden-artifact comparison it exists for never ran. With the error fixed, the
   comparison ran and failed: the **IMP** (paternal-only) TBVs differ from the
   July golden by **exactly sqrt(2) at every individual**, while ADG,
   haplotypes, dosage and the exported matrix are bit-identical.

   That factor is not drift, it is the fix. Phase C made `scale_to_target`
   origin-aware — `V_A = sum_j n_eligible,j * p_j q_j a_j^2` with
   `n_eligible = 1` for a parent-qualified copy instead of 2 — so the **old
   golden recorded a paternal-only trait carrying half its requested additive
   variance**. `tbv.rds` was re-captured (the other four goldens were left
   alone); the property itself is pinned independently by gate 43, so the golden
   is a regression net rather than the specification. The file's header now
   records which artifact was re-captured, when, and why.

   The general lesson is about the shape of the Phase D report, not about
   parity: **an error earlier in a test file hides every assertion after it**,
   so "7 errors, 0 failures" was never 7 known problems — it was 7 known
   problems plus an unknown number of unreached assertions.

8. **The vignette's re-run note was stale since Phase C.** It claimed
   `define_additive_effects()` "always overwrites its rows in `genome_effects`".
   Since 0.66.0 it replaces the variant at the *same scope* and leaves other
   scopes standing — which is the whole point of `replace_scope`, and the thing
   that makes successive common / line-A / line-B calls compose.

---

## Gate status

| Gate | Result |
|---|---|
| 49 — pre-change database | ✅ `restore_pop()` stops with the term/member message; the connection is released so the file reopens. Three cases: the old one-row-per-locus table, a term table whose children are gone, and (★ beyond the gate) a 0.65.0–0.67.0 file carrying `genome_effect_types` |
| 48 — fresh vs restored parity | ✅ still met with the widened view; `genome_effect_loci` is registered in all three schema lists and its reserved-column list now includes `genome_value` |
| "No legacy columns remain" | ✅ 19 grep hits over `R/`, `tests/`, `vignettes/`, `man/`, every one deliberate — roxygen explaining what `parent_origin` replaced, absence-assertions, the guard's own message, and the fixtures that build the old shape in order to be refused. Itemized under §Verification |

---

## Verification

**Full suite** — `testthat::test_dir("tests/testthat")` against
`pkgload::load_all()`:

```
pass 2772   fail 0   err 0   warn 11   skip 1
```

The one skip is the pre-existing `add_founders handles large lines efficiently`
performance guard. Of the 11 warnings, 10 are the same pre-existing set Phase D
reported; the 11th was incidental to my own new fixture — it restored a
genome-less database, so `restore_pop()`'s unrelated "No 'genome_meta' table
found" warning fired before the guard under test. The fixture now defines a
genome, which also makes the test say what it means (the genome-effects guard
passes, *then* the second check fires), and `test-restore_pop.R` re-runs at
**27 pass / 0 fail / 0 error / 0 warning**.

**Per file, after every change:**

| File | Result |
|---|---|
| `test-extract_genotypes.R` | 31 pass (was 29 with 5 errors before the phase) |
| `test-restore_pop.R` | 27 pass (was 22) |
| `test-genome-effects-schema.R` | 100 pass (was 98) |
| `test-genome-effects-writer.R` | 183 pass |
| `test-sql_injection_hardening.R` | 18 pass |
| `test-define_additive_effects.R` | 44 pass |
| `test-add_tbv.R` | 22 pass |
| `test-parity.R` | 6 pass, 0 skip (capture run: 4 pass, 1 skip) |
| `test-phenotype_composite.R`, `test-define_phenotype.R` | 35 / 38 pass — the `component_names` rename |
| `test-schema-print.R`, `test-schema-registries.R`, `test-open_pop.R` | 56 / 110 / 26 pass — registries agree after both schema changes |

**The parity drift, measured rather than assumed.** Before re-capturing, the
fresh run was compared to the July golden individual by individual:

```
ADG  ratio  min 1.000000  max 1.000000
IMP  ratio  min 1.414214  max 1.414214      sqrt(2) = 1.414214
haplotypes identical: TRUE
dosage     identical: TRUE
export_vals identical: TRUE
export_order identical: TRUE
```

Every ADG value and every non-TBV artifact is bit-identical; every IMP value is
off by exactly `sqrt(2)`. That is the signature of `n_eligible` going from 2 to
1 in `V_A = sum_j n_eligible,j * p_j q_j a_j^2`, and nothing else.

**No legacy name remains.** `grep -rn "genome_effect_type\|genome_effect_types\|
base_allele_freq\|expressed_parent"` over `R/`, `tests/`, `vignettes/` and
`man/`, minus `compute_base_allele_freq()` (a live function), returns 19 lines,
every one deliberate:

| Lines | What |
|---|---|
| 4 | roxygen (3 in `R/`, 1 generated in `man/`) explaining what `parent_origin` *replaced* |
| 5 | `expressed_parent` absence-assertions in `test-define_trait.R`, plus one fixture note |
| 3 | `base_allele_freq` assertions: `genome_meta` has no `base_allele_freq_ADG`; `terms` rejects the column by name |
| 3 | the `restore_pop()` guard — a comment and the message that has to name the old column |
| 4 | the gate-49 / 0.68.0 fixtures, which have to *build* the old shape in order to be refused |

---

## Deliberately not done

- **Activating `component_names`.** Out of this plan by its own text. The column
  is renamed, not wired.
- **A migration path for pre-0.65.0 files.** There is no shape-preserving
  mapping from one row per locus to terms with scoped members in the general
  case, and pre-1.0 there is nothing to preserve. The guard says so in one
  sentence rather than half-migrating.
- **An additive-QTL / epistatic-only breakdown.** The plan puts it in
  `describe_table()` or a summary helper, not the print header, and nothing yet
  needs it.
- **`ind_tbv` / `ind_tgv` consolidation.** `plans/consolidate_genetic_values.md`.

---

## Next

The v4 plan is complete. The branch `feat/genome-effects-v49` carries Phases
B–E and is ready to merge to `main`.

The two items this plan recorded and did not schedule stay open:

- **The mismatched-centre trap** (`..._phase_D.md`, §Tracked for later) — a
  hand-written `ad_terms(p = ...)` centred somewhere other than the stored
  additive term's realized frequency makes `tbv_value` stop being the model's
  breeding value, and the writer cannot reject it.
- **Realized variance reporting** (§Future limitations item 6) — the only check
  that measures `E[x_A x_D]` against the real population instead of assuming
  HWE at a declared `p`. It is also the fix that would turn the item above from
  a warning into a measurement.

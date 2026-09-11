# Phase D — the evaluator

*Companion to `plans/update_genome_effects_v4.md` (v4.9). Phase A:
`update_genome_effects_phase_A.md`. Phase B: `..._phase_B.md`. Phase C:
`..._phase_C.md`.*

Phase D builds the one piece the previous three phases were storage for: a
function from stored rows to a number. `add_tgv()` writes every component of it
to `ind_tgv`; `add_tbv()` is the same evaluator filtered to the reserved owner
and to order-one `additive` terms — **not** a second implementation.

The suite goes from **2194 pass / 144 error** (end of Phase C) to
**2746 pass / 0 fail / 7 error**. The 7 are one root cause, in
`extract_genotypes()`, and are Phase E's scheduled work.

---

## What shipped

| File | Status | What it is |
|---|---|---|
| `R/genome_effects_eval.R` | **new**, ~600 lines | The evaluator: model reader, family partitioner, label alphabets, slot-freedom analysis, resolved variant map, containment resolver, preflight, the single evaluation statement, and the shared entry-point plumbing |
| `R/add_tgv.R` | **new**, ~185 lines | `add_tgv()` and its transactional idempotent write |
| `R/add_tbv.R` | rewritten | The TBV math replaced by a filtered evaluator call; the index block untouched |
| `R/add_phenotype.R` | edited | The genome-effect precondition reads the new shape |
| `R/genome_effects_helpers.R` | edited | `.ge_family_key()` → vectorized `.ge_family_keys()` |
| `R/define_additive_effects.R` | edited | One call site of the above |
| `R/schema.R` | edited | One column description naming a deleted concept |
| `tests/testthat/test-genome-effects-eval.R` | **new**, 209 assertions | Gates 1–10, 14, 15–18, 21–23, 25–33, 38, 39, 45, 51, 52, and Q1 |
| `tests/testthat/helper-genome-effects-db.R` | **new** | `ge_flat_view()`, deduplicated out of `test-define_additive_effects.R` |
| `tests/testthat/test-add_tbv.R` | edited | The oracle rewritten against the term shape; gate 40 unchanged and green |
| `dev/benchmarks/benchmark_tgv_scale.R` | **new** | Gate 51's wall-clock half |
| `vignettes/tidybreed-introduction.Rmd` | edited | Causal-locus prose and the `genome_effect_terms` view |
| `CLAUDE.md`, `NEWS.md`, `DESCRIPTION` | edited | 0.67.0 |

## The shape of the evaluation

The semantic definition is per **evaluation tuple** — one unit from each member,
one variant selected per tuple by predicate containment. It is never executed
that way. From §Evaluation strategy:

1. An origin predicate reads a copy's `(line_origin, parent_origin)` **label**
   and nothing else, so the winning variant is a function of the label.
2. Tuples group by label-vector, and the inner sum factors inside each group.

So resolution runs once over `#families × ∏_m |labels_m|` rows and never during
evaluation, and the evaluation itself is **five SQL statements** for an additive
model, whatever the population size.

```
read model (3 queries) → families → label alphabets (1 query)
   → slot freedom → preflight → resolved variant map (R, no query)
   → one statement: member reduction ⋈ map → product → sum
```

The member reduction is the part that makes it set-based: each member reduces to
one value per `(id_ind, id_genome_effect, member_slot, label)` **before** members
are combined, so no row ever carries a Cartesian pairing of raw haplotype rows.
`product()` does the multiplication; `HAVING COUNT(*) = n_members` is what makes
a partially matched label-vector contribute nothing rather than a truncated
product.

### The zero-copy state is materialized

An individual with no `ind_haplotype` row at a locus has no row to reduce, so the
genotype grain left-joins against the `(individual × genotype-member locus)` list
and synthesizes `copy_count = 0`. The plan left open what that list should be;
**it is every locus a genotype member names**, not the loci `chr_inheritance`
says the individual should carry — a `(copy_count 0, dosage 0)` indicator is a
"carries no copy here" effect and must fire whatever the reason for the absence.
That is also what the Phase A fixtures do.

### Dominance at a non-diploid state

The reduction emits `NULL` rather than a number, and the aggregate counts those
`NULL`s per selected tuple. A tuple whose predicate does not match is never
selected, so a dominance member legitimately scoped to a multiset of two never
trips it; a *selected* dominance member at a copy count other than 2 stops with
the individuals named. That reproduces the Phase A naive evaluator's `stop()`
exactly (gate 23 / fixture F14).

---

## Corrections to the plan

### 1. ★ A high-order **unscoped** term tripped the resource guard

The plan says the fast path is "for a *common/unscoped* family only. When no
variant in a family carries any origin row, the sum factors exactly into
`∏ over members of (Σ over units)`." My first working implementation enumerated
label-vectors for every family regardless, and the benchmark caught it
immediately: a 50-locus unscoped dominance term reported

```
Label-vector enumeration for family '...' would be 1125899906842624 rows,
above the hard cap of 1000000.
```

for a model whose evaluation is a single `GROUP BY`.

The fix is sharper than the plan's wording, and the sharper form is the one
implemented: freedom is a property of a **member**, not of a family. If no
variant in the family scopes member *m*, then `variant(L)` does not depend on
`L_m`, so

```
Σ over L_m  Σ over units carrying L_m   =   Σ over all units
```

and that member carries the sentinel label `"*"` — one reduced row instead of
one per label. The plan's family-level fast path is the special case where every
member is free. This also handles a **partially** scoped term (one member scoped,
another not), which the family-level rule could not express at all.

`"*"` cannot collide with a real label because every real label contains an `@`.

**Plan updated**: §Contribution's fast-path paragraph now states the per-member
rule, ★-marked.

### 2. ★ The preflight's subject had to follow

`.gev_preflight()` counts `∏_m |labels_m|` with free members contributing 1.
Without that, the guard fires on models it should wave through — which is worse
than not having it, because the message tells the user to reduce the order of a
term that costs nothing.

### 3. ★ `.ge_family_key()` was quadratic, and so was the validator

The plan's family signature is defined per term, and Phase B/C implemented it
that way: `.ge_family_key(trait, owner, members_of_one_term)`. Every caller —
the Phase C validator, the parent-only warning, and now the evaluator — wanted
keys for a **whole model**, and each was filtering the member frame once per
term. On a 500-QTL three-line model that is 1500 filters over a 1500-row frame
inside `validate_genome_effects()`, which runs inside every write transaction.

Replaced by `.ge_family_keys(terms, members)`, one `tapply`. Reading a 1500-term
model dropped from 0.28 s to 0.02 s, and the effect on the writer's validation
path is the same order. All four call sites updated; the scalar form is gone
rather than kept as a wrapper.

This is a latent Phase B/C defect that Phase D found, not a Phase D regression.

### 4. Gate 51's statement-count assertion needed a mechanism

The plan asks that "the statement count is **independent of the number of
individuals**". Counting statements needs a hook. `trace()` on the S4 generic
`DBI::dbGetQuery` counts nothing — dispatch goes to the duckdb method — so the
test traces `dbGetQuery`/`dbExecute` in `asNamespace("duckdb")` with the
signature `c("duckdb_connection", "character")`, and keeps the counter in
`globalenv()` because a traced function in a namespace resolves free variables
through the namespace's parents, which end there.

It is joined by a second, mechanism-free assertion: `.gev_sql()` is a pure
function of four registered frame names, so an id list **cannot** be inlined into
it — which is what a per-individual loop looks like from the outside.

---

## Implementation issues found in review

1. **`tapply()` returns an array, not a vector.** Three expectations compared
   `unname(parts[ids])` against a plain numeric and failed on `dim`. Test-side.

2. **`get(x, envir = e)$n <- 0L` is not a valid assignment target.** The
   statement counter's tracer expanded to a non-language object and errored on
   first use. Replaced with `assign("n", ..., envir = e)` inside `local()`.

3. **The dominance error listed individuals in hash order.** Non-deterministic
   message text is a flaky test waiting to happen; sorted.

4. **`.gev_write_tgv()` unregistered its temp frames before rolling back.**
   `on.exit` handlers run in the order added, and the unregister handler was
   added first, so a failed insert would drop the registered views and *then*
   attempt `ROLLBACK`. Rewritten with `tryCatch`, which makes the ordering
   explicit instead of positional, and both frames are now registered before
   `BEGIN`.

5. **`ind_tgv`'s replacement had to be by trait, not by component.** An upsert
   keyed on `(id_ind, trait_name, component_name)` satisfies gate 22 and is still
   wrong: drop the dominance terms from a model, re-run, and the stale
   `'order1_dominance'` row survives and keeps being summed by `ind_tgv_total`.
   The delete is by `(individual, trait)`. A test pins it.

6. **`.gev_read_model()` returned a differently-shaped frame when empty** —
   no `effect_order`, `component_name` or `family_key`. Harmless today because
   every caller checks `nrow` first; a trap for the next caller. The empty
   return now carries the same columns.

7. **Three redundant model reads per trait.** `.gev_require_terms()` re-read the
   whole model from SQL to count rows, once per trait, and then `.gev_evaluate()`
   read it again. It now takes the model the caller already has; `add_tbv()` and
   `add_tgv()` read once. This is what brings `.gev_evaluate()` to five
   statements.

8. **`MIN(mp.id_genome_effect)` in the aggregate.** Correct — every row of a
   `map_id` group carries the same id — but it reads as though a choice were
   being made. Both columns moved into the `GROUP BY`.

9. **`add_phenotype()`'s precondition required a population-wide term.** The old
   check was `genome_effect_type = 'additive' AND line_name IS NULL`. Carried
   over literally, that would reject a purebred multi-line design whose only
   terms are line-specific — which evaluates perfectly well. The requirement is
   now "a term `add_tbv()` can read": order one, `additive`, reserved owner.

10. **`add_tbv()`'s missing-match message contained a garbled clause** —
    *"…exclude them from the individuals carry, or exclude them from the
    subset"*, a pre-existing copy-paste in `add_tbv.R:233-234`. Rewritten, and
    it now also names a line-scoped term as a cause, not only a parent-scoped
    one.

---

## Performance

`dev/benchmarks/benchmark_tgv_scale.R`, 2000 loci / 500 QTL / two lines, MacBook:

| Model shape | n = 200 | n = 800 |
|---|---|---|
| common (one unscoped variant per locus) | 0.21 s (1.05 ms/ind) | 0.31 s (0.39 ms/ind) |
| lines (common + A + B) | 0.53 s (2.64 ms/ind) | 0.63 s (0.78 ms/ind) |
| dominance (additive + 50 dominance terms) | 0.18 s (0.88 ms/ind) | 0.23 s (0.28 ms/ind) |

Per-individual cost **falls** with population size: total time is nearly flat,
which is the signature of a fixed model cost plus a set-based scan. The three
optimizations that got it there, in order of effect:

1. **`.ge_family_keys()` vectorized** — the single largest win.
2. **Resolution cached on a family signature that excludes `locus_id`.** A
   500-QTL model with common/A/B variants is 500 copies of *one* resolution
   problem. Solving it once took the `lines` shape from 2.85 s to 1.9 s; the key
   is built from raw origin rows, so a cache hit costs no predicate
   construction either.
3. **Column-wise splits instead of data-frame subsetting in the per-family
   loop**, and one `data.frame()` at the end instead of 500 plus `rbind`.

Net: 2.85 s → 0.53 s on the `lines` shape, and the `common` shape is faster than
the engine it replaces — the correlated `NOT EXISTS` line subquery is gone,
replaced by a join against a precomputed map, exactly as the plan predicted.

---

## Gate status

| Gate | Status | Where |
|---|---|---|
| 1–10 | ✅ | Every Phase A fixture, through the real tables, against `gefx_eval_naive()` — for every individual, not only the named ones |
| 11–13 | ✅ Phase C | Writer gates |
| 14 | ✅ | Injected dominance *and* injected interaction each move `ind_tgv_total` |
| 15–18 | ✅ | Fixtures F19 (three families at one locus), F09 (copy counts), F11 (partial specificity) |
| 19–20 | ✅ Phase C | |
| 21 | ✅ | All four components present, one row per (individual, component), components sum to the derived total |
| 22 | ✅ | Idempotent, plus the stale-component test |
| 23 | ✅ | F14 stops; an autosomal dominance term never does |
| 24 | ✅ Phase C | |
| 25 | ✅ | Functional vs Cockerham differ by exactly `μ`; repeated `ad_terms()` appends accumulate per locus |
| 26–31 | ✅ | Fixtures F02 (parent-qualified beats generic), F08/F12 (reciprocals), F09 (copy counts 0/1/2 at dosage 0), F13 + a direct zero-copy test |
| 32–33 | ✅ | Custom surface and custom interaction move `tgv_value`, leave `tbv_value` identical |
| 34–36 | ✅ Phase C | |
| 37 | ✅ Phase C | Write-time ploidy proof |
| 38 | ✅ | Warns at the configured threshold, stops above the cap, silent at defaults — **and** an unscoped high-order term does neither |
| 39 | ✅ | `order1_other`, summing into the total |
| 40 | ✅ | `test-add_tbv.R` imprinting fixture: Duroc/Landrace variants, per-line centres, F1 offspring, paternal expression, value against the independent per-copy oracle |
| 41–44 | ✅ Phase C | |
| 45 | ✅ | Paternal X term contributes 0 beside an autosomal term; alone it errors for males and evaluates for females |
| 46–48 | ✅ Phase B | |
| 49 | ⏳ Phase E | `restore_pop()` guard |
| 50 | ✅ Phase C | |
| 51 | ✅ | Statement count equal for 8 and 400 individuals (5 statements); no id can reach the SQL text; benchmark script for wall-clock |
| 52 | ✅ | The resolved map for 8 and 400 individuals is identical; a mixed-line multi-scope model agrees value-for-value with an independent per-copy computation |
| 53–54 | ✅ Phase C | |

---

## Verification

```
2746 pass · 0 fail · 7 error · 10 warn
```

Phase C ended at 2194 / 0 / 144. The 10 warnings are pre-existing and unrelated.

**All 7 remaining errors are one root cause and are Phase E's scheduled work**:
`extract_genotypes(effects_tbl = )` reads `genome_effects.locus_name`, which now
lives only in the `genome_effect_loci` view. Six are `test-extract_genotypes.R`;
the seventh is `test-parity.R`, whose fixture calls the same path. The plan
assigns this to Phase E ("move QTL extraction and `extract_genotypes()`'s
nominal `table_name` check to the locus view"), and it is left alone rather than
folded into Phase D.

Per-file, verified directly:

| File | Result |
|---|---|
| `test-genome-effects-eval.R` | 209 / 0 |
| `test-genome-effects-writer.R` | 183 / 0 |
| `test-genome-effects-schema.R` | 98 / 0 |
| `test-genome-effects-fixtures.R` | 248 / 0 |
| `test-add_tbv.R` | 22 / 0 |
| `test-add_tbv_index.R` | 47 / 0 |
| `test-define_additive_effects.R` | 44 / 0 |
| `test-add_phenotype.R` | 32 / 0 |

---

## Deliberately not done

- **`ind_tbv` is not consolidated into `ind_tgv`.** The plan's §Known,
  time-boxed redundancy is unchanged: `ind_tbv.tbv_value` is the reserved
  owner's breeding value, `ind_tgv`'s `'order1_additive'` sums every
  additive-structured term including custom ones. They differ the moment someone
  hand-writes an additive term under `custom`. `add_phenotype()` reads `ind_tbv`
  for simple, composite and SGE assembly; rewiring it is Task 3 of
  `plans/consolidate_genetic_values.md`.
- **No realized-variance reporting.** §Future limitations item 6 calls it "the
  first thing to build after this plan"; it is still after this plan.
- **`phenotype_components.genome_effect_types`** is written, never read, and
  names a vocabulary (`genome_effect_type`) that no longer exists. Its
  *description* is corrected; the column itself belongs to Phase E's "no legacy
  columns remain".

## Q1, decided during review

**Should `add_tbv()` warn when a trait has terms outside the reserved owner?**
Decided: **yes, but on a sharper condition than the plan's options offered.**

None of the plan's three options is right, because "the trait has non-additive
terms" is not the condition under which `tbv_value` stops being the model's
breeding value. The real condition is whether some term contributes to the
additive component or shifts the coefficients `add_tbv()` reads:

| Non-reserved term | Warns | Why |
|---|---|---|
| order-one `additive` | yes | It is part of A and is skipped — owners always sum |
| `indicator` surface | yes | Raw functional coding: the stored `a` is no longer `α = a + d(q − p)` |
| interaction (≥ 2 members) | yes | Additive projection depends on other loci and on LD |
| order-one `dominance` at the additive term's centre | **no** | HWE-orthogonal; `tbv_value` exact |
| order-one `dominance` centred elsewhere | yes | Orthogonality is a property of the centring, not the contrast name |

The exception is the point. A Cockerham dominance model is the common way to
write dominance in this package, and under it the additive coefficient *is* the
average effect — `tbv_value` is exact and there is nothing to say. Warning there
would train the user to ignore the message in exactly the case where it matters.

`.gev_warn_tbv_stale()`, one warning per trait, no behaviour change: a test
asserts `ind_tbv` is byte-identical with and without the warning. Seven tests
cover the table above.

This also drove a small refactor. `add_tbv()` now reads the **whole** stored
model once and narrows it in R through `.gev_reserved_additive()` — one named
definition of "what counts as a breeding-value coefficient" — because the check
needs the terms the function does *not* read. `.gev_read_model()`'s
`order1_additive_only` parameter is gone rather than kept alongside it.

### ★ Tracked for later — the mismatched-centre trap

The last row of that table is a case **nobody had written down before**, found
while implementing the check rather than while designing it. It deserves
following up on its own, because a user can walk into it without doing anything
that looks wrong:

```r
# centre = the realized frequency in the current population, say 0.41
pop |> get_table("genome_meta") |> define_additive_effects("ADG", base = "current_pop")

# centre = 0.30, because that is what the user typed
pop |> define_genome_effects("ADG", ad_terms("Locus_10", a = 0, d = 0.5, p = 0.30,
                                             coding = "cockerham"))
```

Both calls are legal, both are reasonable in isolation, and the result is a
dominance contrast that is **not** orthogonal to the additive contrast sitting at
the same locus. Reproduced on a 4-locus fixture: the additive member stores
`center_value = 0.3333` (the realized frequency) and the dominance member
`0.3000` (the typed `p`), both writes are accepted, and Phase D's warning fires
at `add_tbv()` — so the stored additive coefficient quietly stops being the
average effect. Nothing in the writer objects: `center_value` is deliberately
excluded from the family signature (a line-specific variant legitimately centres
against its own line's frequency), so the two terms are different families and
simply sum.

Phase D's warning catches it **at `add_tbv()` time**, which is the right safety
net but the wrong moment — by then the model is stored and the user is reading a
number. Three things are worth considering later, in rough order of appeal:

1. **Warn at write time**, in `define_genome_effects()`, when a `dominance`
   member lands at a locus whose existing additive variant carries a different
   `center_value`. That is where the mistake is made, and the writer already has
   both rows in hand. It cannot be an *error* — different centres across scope
   variants are legal and intended — so it has to be a warning with the two
   values named.
2. **Report realized orthogonality** as part of the variance-reporting work
   (§Future limitations item 6), which can measure `E[x_A x_D]` against the
   actual population instead of assuming HWE at a declared `p`. That is the only
   check that is *true* rather than merely consistent.
3. **Let `ad_terms()` read `p` from the stored additive term** when the locus
   already has one, instead of requiring the user to retype it. Removes the
   opportunity rather than detecting the mistake, but only helps the `ad_terms()`
   path.

None is scheduled. Recorded here so the decision is not re-derived from scratch,
and so the one-line table row does not have to carry the whole argument.

## Next

**Phase E.** Delete the last old-shape remnants; `restore_pop()` guard for
pre-change files (gate 49); `extract_genotypes()`'s `effects_tbl` path and the
`genome_effect_type` mention in its roxygen, moved to `genome_effect_loci`. That
turns the remaining 7 errors green and is the last phase before this branch
merges.

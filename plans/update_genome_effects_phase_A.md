# Genome Effects — Phase A results

**Spec:** `plans/update_genome_effects_v4.md` (v4.9). **Branch:** `feat/genome-effects-v49`.
**Status:** complete. Phase A gate met; **no change to the v4.9 schema was required.**
**Date:** 2026-09-10.

Phase A writes no package code and creates no tables. It exists to answer two
questions before the DDL is committed:

1. Is every case in the origin truth table **representable** in the three proposed
   tables, and accepted or rejected by the stated rules?
2. Is its value **derivable from the stored rows alone**, and does it equal a number
   computed by hand?

Both are now answered yes, for 19 fixtures.

---

## What shipped

| File | Contents |
|---|---|
| `tests/testthat/helper-genome-effects.R` | Genotype context (9 individuals, 6 loci); row constructors for the three tables; the validator (row-local `CHECK`s + cross-row rules); **two independent evaluators**; the 19-fixture registry with hand-computed expected values |
| `tests/testthat/test-genome-effects-fixtures.R` | 248 assertions over the registry |

Nothing here is throwaway. The registry is the fixture set Phase B's validator tests
and Phase D's evaluator tests consume, and the naive evaluator is the oracle Phase D's
SQL is checked against — the same pattern as the committed oracle at
`tests/testthat/test-add_tbv.R:16-40`.

**Two evaluators, deliberately.** `gefx_eval_naive()` implements the semantic
definition literally: Cartesian product over evaluation units, one variant selected per
tuple. `gefx_eval_grouped()` implements the label-vector factorization from
§Evaluation strategy. Every fixture asserts
**naive == hand-computed == grouped**, so v4.9's central performance claim is proven
before any SQL is written against it.

**Loci are named, not id'd.** `locus_name` → `locus_id` resolution is the writer's job
in Phase C and has no bearing on the semantics under test.

---

## The genotype context

Nine individuals; an absent copy is an absent row.

| Individual | What it is for |
|---|---|
| `pure_A` | both copies line A |
| `f1_AB` | sire A, dam B |
| `f1_BA` | the reciprocal — sire B, dam A |
| `unk_A` | one line-A copy, one copy of unknown founding line (`NULL`) |
| `male_X` / `male_X0` | hemizygous at `LX` (one copy, from the dam), allele 1 / allele 0 |
| `fem_XX0` | two `LX` copies, both allele 0 — dosage 0 at copy count 2 |
| `fem_noY` | carries no `LY` copy at all — the synthesized zero-copy state |

Centring constants: `p_A = 0.4`, `p_B = 0.6`, common `p = 0.5`.

---

## Fixtures and hand computations

Contrast definitions used throughout: additive `x_A = Σ over eligible copies (allele − c)`;
dominance `x_D(0) = −2p²`, `x_D(1) = 2pq`, `x_D(2) = −2q²` with `c = p`;
indicator `x_I = 1` iff `(copy_count, dosage)` equals the member's state.

### Valued fixtures

**F01 — common vs A-specific additive @L1.** Common `a=2.0, c=0.5`; `(exact A)` `a=3.0, c=0.4`.
Per-copy fallback; the A-specific scope is a strict subset of the common one.

- `f1_AB`: A copy `3.0(1−0.4)=1.8` + B copy `2.0(1−0.5)=1.0` → **2.8**
- `pure_A`: `3.0(1−0.4) + 3.0(0−0.4) = 1.8 − 1.2` → **0.6**
- `unk_A`: `1.8` + NULL copy on the common variant `2.0(0−0.5) = −1.0` → **0.8**

**F02 — generic A vs paternal A @L2.** `(exact A, ANY)` `a=1.0`; `(exact A, parent 1)` `a=5.0`; `c=0.4`.
This is the tie v4.1's row-count rule got wrong. Non-A copies match nothing and score 0.

- `pure_A`: paternal `5.0(1−0.4)=3.0` + maternal generic `1.0(1−0.4)=0.6` → **3.6**
- `f1_AB`: paternal A allele 0 → `5.0(0−0.4)` → **−2.0**
- `f1_BA`: paternal B unmatched; maternal A generic `1.0(0−0.4)` → **−0.4**

**F03 — common vs A/B dominance @L1.** Common `d=4.0`; `{A:1, B:1}` `d=10.0`; `c=0.5`.

- `f1_AB`: multiset `{A,B}` matches, strictly inside common; dosage 2 → `10.0 × (−0.5)` → **−5.0**
- `pure_A`: multiset `{A,A}` disjoint from `{A,B}` → common; dosage 1 → `4.0 × 0.5` → **2.0**

**F04 — generic A/B vs one reciprocal @L3.** `{A:1,B:1}` `d=2.0`; `{A:1@p1, B:1@p2}` `d=7.0`; `c=0.5`; dosage 1 so `x_D = 0.5` in both directions.

- `f1_AB`: reciprocal matches and is strictly inside generic → `7.0 × 0.5` → **3.5**
- `f1_BA`: reciprocal needs A from the sire → falls back to generic → `2.0 × 0.5` → **1.0**
- `pure_A`, `unk_A`: neither multiset matches → **0**

**F06 — common vs origin-specific A×A @L1,L2.** Common `k=1.0, c=0.5`; both members `(exact A)` `k=4.0, c=0.4`.
The fixture that shows the sum does **not** factor into a plain product of per-member sums.

`f1_AB`, four label-vectors:
`(A,A) → 4.0(1−0.4)(0−0.4) = −0.96` · `(A,B) → 1.0(0.5)(0.5) = 0.25` ·
`(B,A) → 1.0(0.5)(−0.5) = −0.25` · `(B,B) → 1.0(0.5)(0.5) = 0.25` → **−0.71**
(a naive all-common factorization would give `1.0 × 1.0 × 0.0 = 0.0`).

`pure_A`: every label is A, so the scoped variant is selected throughout →
`4.0 × [(1−0.4)+(0−0.4)] × [(1−0.4)+(1−0.4)] = 4.0 × 0.2 × 1.2` → **0.96**

**F07 — partial specificity at one member of two @L1,L4.** Common `k=1.0`; scoped `k=3.0`
with `(exact A)` on member 1 only. All four alleles are 1 and `c=0.5`, so every contrast
value is 0.5.
`(A,A) → 0.75` · `(A,B) → 0.75` · `(B,A) → 0.25` · `(B,B) → 0.25` → **2.0**
All-common would be 1.0 and all-scoped 3.0, so the number separates all three readings.

**F08 — two disjoint specifics @L1.** `(exact A, parent 1)` `a=2.0`; `(exact B, parent 2)` `a=5.0`; `c=0.5`.

- `f1_AB`: `2.0(0.5) + 5.0(0.5)` → **3.5** — both apply, to different copies
- `f1_BA`: sire B fails variant 1 (wrong line) and variant 2 (wrong parent); dam A likewise → **0**
- `pure_A`: paternal copy only → **1.0**

**F09 — indicator states at `LX`.** Three terms: `(0,0) → 10.0`, `(1,0) → 20.0`, `(2,0) → 30.0`.
Three *families* (copy count is in the signature), so they sum rather than compete.

`fem_noY` **10.0** · `male_X0` **20.0** · `fem_XX0` **30.0** · `male_X` (1,1) **0** · `pure_A` (2,2) **0**.
Under dosage alone all three of the first row are "dosage 0".

**F10 — imprinting as `('any', parent 1)` @L2.** `a=2.0, c=0.5`.
`pure_A` **1.0** · `f1_AB` **−1.0** · `f1_BA` **1.0** — numerically the deleted
`expressed_parent = "parent_1"` filter, but per locus and per effect owner.

**F11 — line-specific imprinting @L1.** `(exact A, parent 1)` `a=1.0, c=0.4`;
`(exact B, parent 1)` `a=4.0, c=0.6`. Disjoint.
`f1_AB` → `1.0(1−0.4)` = **0.6** · `f1_BA` → `4.0(1−0.6)` = **1.6**.
Stamping `'any'` instead of composing line with parent would collapse both variants onto
one scope and lose the per-line centre — the wrapper bug v4.7 fixed, now pinned by a number.

**F12 — `{A, NULL}` with an `unknown` variant @L1.** Common `a=2.0, c=0.5`;
`(exact A)` `a=3.0, c=0.4`; `(unknown)` `a=7.0, c=0.5`.
`unk_A`: `1.8` + `7.0(0−0.5) = −3.5` → **−1.7**. `exact` and `unknown` are disjoint;
both sit strictly inside common.

**F13 — hemizygous centring @LX.** Common `a=2.0, c=0.5`.
`male_X` = `2.0(1−0.5)` = **1.0** (not `dosage − 2c` = 0) · `male_X0` **−1.0** ·
`pure_A` **2.0** · `fem_noY` **0** (no eligible copy).

**F19 — three families at one locus all fire.** Additive main `2.0` + dominance main `4.0`
+ A×A interaction `1.0` over L1,L4, for `f1_AB`:
`2.0 × 1.0` + `4.0 × (−0.5)` + `1.0 × 1.0 × 1.0` = `2.0 − 2.0 + 1.0` → **1.0**.
Origin precedence never makes them suppress one another.

### Rejection fixtures

| # | Case | Enforced by |
|---|---|---|
| F05 | `(exact A, ANY)` vs `(ANY, parent 1)` in one family — overlap, neither contains | R validator |
| F14 | `dominance` evaluated at copy count 1 | write-time in the real system (`chr_inheritance`); pinned at the evaluator here |
| F15 | duplicate family + scope identity | R validator |
| F16 | mixed centring at one locus (Q2) — collides as a duplicate | R validator; combined equivalent `a=3.0, c'=0.4` asserted |
| F17 | `('any', parent NULL)` — the common scope written the long way | row-local `CHECK` |
| F18 | `'any'` on a genotype member | cross-table (R) |

---

## What Phase A confirms

1. **The v4.9 schema needs no change.** No fixture required a column, table or
   `line_match_type` value that v4.9 does not already have. The storage design survived
   first contact with the full origin truth table.
2. **The label-vector factorization is correct.** Naive per-tuple and grouped
   label-vector evaluation agree to `1e-12` on every valued fixture × every individual.
   This was an unverified claim when v4.9 was written; it is now the thing 13 fixtures
   test. Phase D can build the SQL on it.
3. **Predicate containment is decidable and behaves as the worked cases state** —
   including the three cases v4.1 got wrong (parent-qualified beats generic; reciprocals
   are disjoint; line-only vs parent-only is incomparable and rejected).
4. **The per-tuple-zero / trait-level-error split is implementable.** A tuple that
   matches nothing contributes 0 inline (F02, F08); an individual for whom every term is
   unmatched (F08 `f1_BA` = 0 with no copy matched anywhere) is a separate check across
   families, not something the evaluator does inline.
5. **`center_value` belongs to the selected variant.** F01 and F12 both mix centres
   across scopes of one family, so a member's reduction cannot be precomputed
   independently of which variant won. v4.9's artifact table already keys member
   reduction by `id_genome_effect` — that is load-bearing, not incidental.

---

## What Phase A exposed (carry into B–D)

**1. Exact-multiset matching must be a full bijection search, not greedy.** *(Phase C)*
A demand set mixing a parent-qualified row with an ANY-parent row — `{A:1@p1, A:1}` — can
be satisfiable while greedy consumption fails it, because the ANY row can eat the copy the
qualified row needed. The Phase A validator and evaluator both use exhaustive matching
(`gefx_bijection()`), and Phase C's validator must do the same. Sizes are ≤ 2 items at
diploidy, so the cost is nil.

**2. Genotype containment needs the same search.** *(Phase C)* "P's parent assignments
refine Q's" is a bijection question for the identical reason (`gefx_refines()`).

**3. When both reciprocals are defined, the generic `{A,B}` variant is unreachable.**
*(documentation, possibly a writer warning)* Two reciprocals partition the `{A,B}` label
space, so a generic fallback beneath them can never be selected. Not a bug, but a user who
writes all three has written one term that does nothing. F04 was built with **one**
reciprocal precisely so the fallback path is actually exercised. Candidate for the same
class of warning as the parent-only re-run warning (gate 50).

**4. The zero-copy indicator is sharper than the plan's Y-in-female framing.** *(Phase D)*
`(copy_count 0, dosage 0)` matches **every** individual with no row at that locus — in F09,
`f1_AB`, `f1_BA` and `unk_A` all score 10.0 at the `LX` `(0,0)` term. That is semantically
correct: a `(0,0)` indicator is a "carries no copy here" effect. It does raise a question
the plan does not answer: **what does the zero-copy left join join against — every locus in
`genome_meta`, or the loci the individual is expected to carry under `chr_inheritance`?**
The two differ for an individual whose chromosome is absent for their sex versus one whose
rows are simply missing. Recommend joining against the full locus list (the fixture
behaviour), and treating a missing row that `chr_inheritance` says should exist as a data
error surfaced elsewhere. **Open for Phase D.**

**5. `product()` is not enough on its own for higher-order terms.** *(Phase D)* Because the
selected variant varies by label-vector (F06, F07), the SQL cannot be "reduce every member,
then multiply". It must be "reduce every member **per label**, join on the label-vector,
then multiply within the selected variant". v4.9 says this; F06 is the number that fails if
an implementation forgets it.

---

## Gate status

| Phase A gate | Status |
|---|---|
| Every fixture representable | ✅ 19/19, no schema change |
| Value derivable from stored rows alone | ✅ 13 valued fixtures, hand-computed, all agree |
| Containment order settled before DDL | ✅ additive and genotype lattices tested directly |
| Multi-locus fallback settled before DDL | ✅ F06, F07 |
| Each fixture expressed in the `terms` format | ⚠️ **deferred to Phase C** — see below |

**The one deviation.** v4.9's Phase A row asks that each fixture also be written out in the
`terms` input format, as an early usability check on that format. Phase A instead expresses
fixtures at the **stored-row** level (`gefx_term()` / `gefx_member()` / `gefx_origin()`),
because that is what proves representability and derivability — the two things the gate
actually turns on. Writing them a second time in `terms` form only becomes meaningful when
there is a writer to canonicalize them, and the natural test is then a **round-trip**:
`terms` → stored rows → compare against these fixtures. That is gate 53, and it lands in
Phase C. The registry is built to be that comparison target.

---

## Verification

```
testthat::test_file("tests/testthat/test-genome-effects-fixtures.R")
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 248 ]
```

Mutation-checked: inverting the specificity selection in `gefx_select_variant()` (pick the
least specific matching variant instead of the most specific) fails the suite, so the
fixtures constrain behaviour rather than merely executing it. `test-add_tbv.R`,
`test-schema-registries.R` and `test-schema-print.R` still pass — the new helper adds only
`gefx_*` names and touches nothing else.

---

## Next

**Phase B** — `genome_meta` gains its `PRIMARY KEY`; the `open_pop.R:286` `genome_effects`
DDL is deleted **in the same commit** that adds the three effect tables to `GENOME_TABLES`;
effect tables move into `define_genome()`; 4 tables + 24 registry entries; SQL constraints;
containment checker; R validator; the three views registered in all three schema lists.
Gates 46–48 and 54.

Phase B is the first phase that leaves the package **red** — `add_tbv.R` and
`define_additive_effects.R` still read the old `genome_effects` shape until Phase D. B, C
and D stay on this branch and only merge once D is green.

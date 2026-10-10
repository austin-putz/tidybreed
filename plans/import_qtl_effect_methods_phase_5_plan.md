# Import QTL-effect methods — Step 5 plan (Part C)

**Spec:**
- `plans/import_qtl_effect_methods.md`: §1.3 (the non-additive method), §3 (the mapping),
  §4.3, §5 (owners and replacement), §6C (targets), §7 (Part A, whose machinery is reused),
  §9 (Part C), Q5, Q8, Q12, Q14, Q17, Q21 and §10 Step 5.
- Gates C1–C15 and C17–C20 (§11). C16 shipped in step 4.
- Source: `$SRC = /Users/austinputz/Claude/simulate_qtl_effects`, commit `318e54f`;
  `non-additive/R/qtl_effects_nonadd.R` (531 lines) and
  `non-additive/tests/test_qtl_effects_nonadd.R` (18 tests).

**Versions:** 0.75.3 (5a), 0.76.0 (5b), 0.76.2 (5c; 0.75.4 and 0.76.1 are the review follow-ups of 5a and 5b). The main plan gives step 5 one version
(0.76.0). This plan splits it in three (decision D1).
**Status:** planned 2026-10-06, revised the same day after the
[Codex review](import_qtl_effect_methods_phase_5_codex_review.md) and again after its
re-review. **All decisions D1–D8 were made by the user (2026-10-06 to 2026-10-09), each
as recommended; D5 adds a rank `message()` for singular targets (5b.7 item 8).** **5a done
2026-10-09 (0.75.3); 5b done 2026-10-09 (0.76.0, review fixes 0.76.1); 5c done 2026-10-09
(0.76.2) — step 5 complete**; results and deviations in
[import_qtl_effect_methods_phase_5.md](import_qtl_effect_methods_phase_5.md).
**Starting point:** 0.75.2 (`0276271`).

**Changes after the Codex review (all eight findings and the smaller points accepted):**

| # | Finding | Where it changed |
|---|---|---|
| 1 | D5's "a singular `G_A` is below the floor anyway" is false (at `p = 0.5`, `C = 0`); zero non-additive blocks need no coupling | D5 rewritten as an implementation restriction with a zero-block exception; 5b.4 dispatch; G1 |
| 2 | `solve_dd_mean()` fails on feasible linear and degenerate cases, and expects a dense `M` | 5b.4 stage 2 (robust solver on anchor objects); G4 |
| 3 | `.noia_terms()` drops zeros and errors on an all-zero model, so additive-only row identity needs a shared *storage* path too | 5b.4 dispatch, 5b.5; C4 gates with zero and rank-one targets |
| 4 | `inbreeding_depression` subset policy, a pre-draw zero-`G_D` check, and no closeness promise for k ≥ 2 | D4, 5b.2, 5b.3, `.na_calibrate()` interface; G5 |
| 5 | Floor-boundary round-off; the source's scale test cannot catch a raw-scale floor; the QR rank convention differs | 5b.4 floor from residual coefficients; unit-transform gate G6; rank convention kept (SVD, cached) |
| 6 | Exact-HWE fixtures still differ from the genic sum by `n/(n−1)` | 5c.1 prevalence (a), (b) |
| 7 | One `n(2m+r)` guard for every call breaks additive-only parity; knowable A×A shortages should fail before any draw | 5b.2 steps 5–7; G7 |
| 8 | "First bad term" is not enough: the order is term, then rule, then member; no current test pins it | D7, 5a.2, new `test-genome-effects-writer-order.R` |
| — | Target precedence across three blocks; rollback mid-target-write; replacement benchmarks; forward-conversion profiling; the delivered-degree summary; oracle isolation | 5b.2 (D8), C11, 5a.3, 5b.7, 5b.9 |

**Changes after the Codex re-review (all five findings and the smaller points accepted;
each was reproduced with the review's probes before editing):**

| # | Finding | Where it changed |
|---|---|---|
| R1 | The stable quadratic loses both roots at `B = 0` (`sign(0) = 0`); reflection leaves a preferred mean already on `ID = 0` at zero variance; the branch scale depends on `sd` | 5b.4 `.na_solve_dd_mean()`: solve in `x = μ/sd`, `sign₊`, `q = 0` branch, finite **positive** `V_D` for every candidate, explicit displacement; G4 |
| R2 | The floor's coordinates were ambiguous (original-unit vs standardised residual) | 5b.4 `.na_additive_stage()`: every quantity named in standardised coordinates, original units for messages only; G6 |
| R3 | Positive `G_D` with both degree parameters zero is a knowable refusal; step 3 used `m` / `n` before step 5 resolved them | 5b.2 reordered (loci and base before counts; guard before collection); new step-1 refusal; G7 / C20 |
| R4 | G7 asked a successful sampling call to preserve `.Random.seed` | G7: preservation for refusals only; the success case matches `define_additive_effects()`'s rows and post-call RNG state; lowered limit in the fixture |
| R5 | The oracle list omitted `nonadd_decompose()` and the source's `.na_*` helpers | 5b.9 oracle isolation: transitive dependency set, explicit parent environment, smoke test |
| — | D8 mixed-source and scope gates; the non-empty-model argument; the floor is conditional on the sampled architecture; D3 gates per block | 5b.9 G8, 5b.5, 5b.4 / 5b.7 / 5b.8, D3 |

Step 5 adds the generator `define_genome_effects()`. It samples additive, dominance and
additive-by-additive (A×A) effects and **calibrates** them so that the stored targets
`G_A`, `G_D`, `G_AA` are delivered exactly under a named anchor. It writes them under the
reserved owner `"generated"`. It is the non-additive sibling of `define_additive_effects()`:
the same pipe subject, target rules, base resolution and owner. Everything it needs to
build on exists:

- the congruence (`R/qtl_congruence.R`, Part A);
- the target resolver and the atomic target write (`.dae_resolve_target()`,
  `.tvc_write_block()`, `.ge_commit(before_commit = )`);
- the NOIA conversion and term builders (`.noia_to_stored()`, `.noia_terms()`,
  `aa_terms()`, step 4, made sparse-safe in 0.75.2);
- the measuring instrument (`extract_genetic_variance()`, step 4);
- the one evaluator, and phenotypes that read the total (step 3).

The one missing piece found in step 4 is **writer speed**. Writing 19,900 A×A pairs took
109–196 s, and profiling shows why (5a). The generator writes through the same engine, so
5a comes first.

---

## Decisions (all decided; see the table)

Each decision has a recommendation. Where the main plan already decided something, it is
not reopened here, except where checking it against the code or the source showed a
problem (D3, D4). Codex agreed with every recommendation (D4 and D5 only after the
corrections now in their text). D8 is new from the review.

| | Decision | Recommendation |
|---|---|---|
| D1 | Three commits (5a / 5b / 5c) | yes — **decided 2026-10-06: three commits** |
| D2 | Draw order A, D, pairs, A×A | yes — **decided 2026-10-06: reordered (A, D, pairs, A×A)** |
| D3 | Loci fixed in the base keep their drawn effects (both generators) | (a), no zeroing — **decided 2026-10-06: (a)** |
| D4 | `inbreeding_depression`: exact for one trait, reported for several | (a) — **decided 2026-10-09: (a)** |
| D5 | Non-zero D / A×A targets require a positive-definite `G_A` (an implementation limit) | accept for this release; generalised solver later — **decided 2026-10-09: (a), accept the limit; singular targets stay valid but get a rank `message()` (5b.7 item 8)** |
| D6 | Diagnostics per block and base kind | the table — **decided 2026-10-09: as proposed** |
| D7 | Writer refactor keeps rows, messages and error order | yes — **decided 2026-10-09: yes** |
| D8 | Passed matrices and `trait_var_comp_tbl` together, per block | (a), allowed per block — **decided 2026-10-09: (a)** |

**D1. Three commits.** **Recommended**, rather than step 5 as one 0.76.0 commit.

| Sub-step | Version | Content | Why separate |
|---|---|---|---|
| 5a | 0.75.3 | Writer speed. Internal only, no API change. | Its gate is "every existing writer test still passes", plus a benchmark. Mixing it into 5b would hide a writer regression inside a new feature. |
| 5b | 0.76.0 | `define_genome_effects()`, the calibration internals, the owner rules, gates C1–C12, C14, C15, C17, C18. | The generator and its algebra. |
| 5c | 0.76.2 | End-to-end: phenotypes (C13), prevalence, `remove_generated_effects()`, the extractor round trips (C19, C20), the vignette. | These cross into phenotypes and the extractor, and the vignette describes the finished set of four paths. |

As in step 4, the user can choose to run them in one pass (each commit with the full suite
green) or to review between them.

**D2. Draw order: additive, dominance, pairs, A×A.** §9.2 step 3 orders the draws
"additive, then random pairs, then dominance, then A×A". **Recommended change:** draw the
random pairs **after** the dominance degrees, immediately before the A×A architecture.
Each block's draws then come after those of every block that is always present, so:

- adding or removing an A×A block changes neither the additive nor the dominance draw;
- adding a dominance block changes neither the additive draw (C17 already requires this)
  nor which loci the additive architecture uses.

The §9.2 order changes the dominance draw whenever an A×A block is added. It costs nothing to
avoid, and it makes "add epistasis to an A + D model" a controlled comparison under one
seed. Gate C17 gains the dominance half.

"Unchanged" means the **architecture draws** (`B_a`, the degrees `h`), not the final
coefficients: stage 3 changes the additive coefficients whenever the coupling `C` changes.
The draw step is its own internal function (`.dge_draw()`), and C17 tests it directly
under one seed, with both supplied and random pairs (Codex review).

**D3. Effects at loci (and pairs) that do not segregate in the base.** The source sets
the effects of *ineligible* loci to 0. A locus is ineligible when the anchor gives it no
variance: `2pq = 0` under `"genic"`, a constant dosage under `"realised"`. The same applies
to D and A×A eligibility (`qtl_effects_nonadd.R:205–216`; source test 6, "fixed loci
carry no effect"; test 13, all-heterozygous loci). **Part A does not.**
`define_additive_effects()` writes a drawn, calibrated effect at every selected locus,
including those fixed in the base. A fixed locus's effect contributes nothing to the
anchor's variance, and it is centred at `p = 0` or `1`, so it contributes nothing in the
base either. But it is real in any other population where the locus segregates: another
line, when the base is a reference line (§6C, "other lines get whatever their frequencies
give").

C4 requires an additive-only `define_genome_effects()` call to write the same rows as
`define_additive_effects()`. So the two generators must make the same choice.

| Option | Effect |
|---|---|
| **(a) Neither generator zeroes** (**recommended**) | QTL are QTL whether or not they segregate in the base. The calibration is exact on the base, and other populations get whatever their frequencies give, which §6C already documents. Source tests 6 and 13 are ported as "contributes zero variance under the anchor" instead of "carries no effect". No change to Part A. |
| (b) Both zero, with no RNG change (zeroing happens after the draw) | Matches the source. Changes Part A's written rows (pre-1.0, allowed): a locus fixed in the reference line carries nothing in other lines, so crossbred line differences lose those QTL. |
| (c) Part C zeroes D and A×A only; additive as Part A | C4 still holds (it is additive-only), but the two blocks then follow different rules in one model. |

(a) keeps one rule across both generators and does not throw away architecture that matters
to crossbreeding. Under (a), a pair with a fixed member still enters the additive stage
through `e·c` (step 4's B12 finding), and the calibration accounts for it exactly.
Gates under (a): the written coefficients at a fixed locus are non-zero; the relevant
centred block has zero variance at the base; and a second population in which the locus
segregates carries that locus's effect (Codex review). The zero-variance assertion is per
block and per locus or pair (the fixed locus's own `A` / `D` contrast columns, or the pairs
containing it in `A×A`), **never** the whole model's variance: a pair with a fixed member
still induces its partner's additive effect through `e·c` (Codex re-review).

**D4. `inbreeding_depression` is exact for one trait only, not for "diagonal `G_D`".** The
main plan (§1.3, §9.1, gate C7) says the inbreeding-depression target is "exact for `k = 1`
or diagonal `G_D`". **It is not exact for k ≥ 2.** The source solves the dominance-degree
mean per trait, and then the congruence `B_d ← B_d A` mixes the columns. With diagonal `G` the
polar-factor congruence gives `A = C^{-1/2} G^{1/2}`, which is diagonal only when the drawn
`C = B_dᵀ M_D B_d` is. Checked on the source itself (`sim_qtl_effects_nonadd`, n = 400,
m = 60, `G_A = diag(2, 3)`, `G_D = diag(0.5, 0.8)`, `inbr_depr = c(1.2, -0.4)`): delivered
`G_D` is exact, but the delivered depression is **1.206 and −0.418** (0.5% and 4.6% off). With
k = 1 it is 1.2 exactly.

| Option | Effect |
|---|---|
| **(a) Accept k ≥ 2; exact only for k = 1** (**recommended**) | Each requested trait's degree mean is solved as if it stood alone; the joint congruence then mixes the columns, so for k ≥ 2 the delivered depression is **approximate, with no closeness guarantee** (the small misses above are one example, not a bound). A `message()` reports requested and delivered depression by trait name, and says why they differ. Nothing is stored (§9.1), so no stored value is false. |
| (b) Refuse `inbreeding_depression` for k ≥ 2 | Honest but loses a useful approximate control. |
| (c) Make it exact for diagonal `G_D` | Calibrate each trait's `d` column alone: a scalar rescale, which preserves the per-trait ratio. That delivers the variances, but the zero dominance covariances only if the drawn columns are `M_D`-orthogonal, which they are not. So it trades one inexact quantity for another. Rejected. |

Under (a), §1.3, §9.1 and C7 in the main plan are corrected: "exact for k = 1; for k ≥ 2
the delivered value is reported".

**Subsets** (Codex review, finding 4). `inbreeding_depression` is a named numeric vector:
distinct, non-missing names, each a call trait, aligned by name (never by position). Only
the named traits' degree means are solved. The other traits keep their drawn degrees
`dominance_degree_mean + dominance_degree_sd · z`. A named trait whose resolved
`G_D[t, t]` is 0 is an error **before any draw**.

**D5. Non-zero dominance or A×A targets require a positive-definite `G_A` — an
implementation limit, not a mathematical one** (revised after the Codex review, finding 1).
The additive stage solves `(B_a T + C)ᵀ M_A (B_a T + C) = G_A` for the k×k factor `T`, and
the ported solver needs `P = B_aᵀ M_A B_a` of full rank k. The source draws `B_a` i.i.d.
normal, so `P` is full rank almost surely. tidybreed must draw it with
`.draw_additive_architecture()`, as `MVN(0, G_A)` rows for k ≥ 2 (the price of C4). When
`G_A` is singular (a genetic correlation of ±1, or a zero variance), that draw is
rank-deficient and the solver cannot run.

The first draft said such a target was "below the floor anyway". **That is false.** The
floor `R − QᵀP⁻¹Q` is positive *semi*definite, not definite. At `p = 0.5` the coupling
`C = diag(b) B_d + E_c` is exactly 0 (`b = q − p = 0`, `c = 2p − 1 = 0`), and Codex's
two-locus example (`B_a = [1 1; −1 −1]`, `B_d = I`) attains a singular `G_A` with full-rank
non-zero dominance exactly.

| Case | Rule |
|---|---|
| No D / A×A block, or every present D / A×A target is the **zero matrix** | No coupling to solve. The call takes the additive-only route (5b.4: Part A's draw, `.qtl_calibrate()` and storage), so any feasible `G_A`, singular included, is accepted. The zero blocks are still resolved, written when passed, and drawn (so D2's draw contract holds), and they produce no terms. |
| A non-zero D or A×A target, and `G_A` positive definite on its correlation scale | The quadratic solver. |
| A non-zero D or A×A target, and `G_A` singular | **Refused before any draw**, worded as a limit of this implementation: "this release calibrates dominance or epistasis only with a positive-definite `G_A` (the sampled additive architecture must have full rank); a correlation of ±1 or a zero additive variance is accepted for additive-only models". It must not claim the target is impossible. |

**Recommended:** accept the restriction for this release. A reduced-rank solver (solve in
the rank-r span of the architecture; feasible iff `G̃ ⪰ 0` and `rank(G̃) ≤ r`) gets its own
design if someone needs it. Rejected: drawing `B_a` i.i.d. whenever a non-additive block is
present. That breaks C4's shared sampler and C17's "adding a block never changes the
additive draw".

**D6. Diagnostics per block (`warn_bounds`, §7.4).** Part A compares the delivered additive
covariance with what another population sees. For Part C, **recommended**, per base kind:

| Base / anchor | Comparison | Report |
|---|---|---|
| `"realised"` | each block's genic limit at the base `p` (closed forms, step 4's `.egv_genic` algebra) | warning per block outside `warn_bounds` |
| `"genic"`, base selects individuals | each block **realised** on the base individuals (step 4's projection, `.egv_realised` algebra on the collected dosages) | warning per block outside `warn_bounds` |
| `"genic"`, founder pool | Part A's additive pool expectation; dominance and A×A are **not compared** (a pool expectation for them needs the haplotypes' multi-locus LD, which nothing computes) | message, as Part A (Q22) |
| `"genic"`, `ind_haplotype` copies | none, as Part A | — |

Nothing is stored (Q7). A size-limited comparison is skipped with a message, as in Part A.
The A×A part of the realised comparison is accumulated in pair chunks (step 4's
`.egv_aa_values()`), so a genic call never allocates an `n × r` pair matrix just to warn.
As in Part A, the comparison is computed **before** the commit and reported after it. The
founder-pool line is labelled "pool expectation (additive only); dominance and A×A not
compared" (Codex review).

**D7. The writer optimisation keeps every row, every message and the error order.** 5a
rewrites loops, not rules. **Recommended.** Messages naming the user's `term_id` stay
word-for-word, and existing test files are not edited. "The first bad term" is not enough
(Codex review, finding 8). Today `.ge_build()` walks term by term, applying its rules in a
fixed order (one coefficient, one label, no locus twice), so term 1's label error beats
term 2's coefficient error. A vectorised version must report the **lexicographically first
violation by (term position in user order, rule position, member position)**. It must not
report the first term that fails whichever rule happens to run first.
`.ge_validate_frames()` returns an ordered **vector** of violations; the whole vector,
including its outer phase order and family-pair order, is preserved, not just its first
element. NA handling and the original value order inside a message stay as they are.

**D8. Passed matrices and `trait_var_comp_tbl` together** (new, Codex review "smaller
improvements"). Part A refuses `G` with `trait_var_comp_tbl`; it has one block, so that is
"one block from two sources". With three blocks the rule must be stated per block.

| Option | Effect |
|---|---|
| **(a) Per block** (**recommended**) | Each block comes from exactly one source: a passed matrix (`G_A`, `G_D`, `G_AA`) or the resolved rows (`trait_var_comp_tbl`, or the stored `line_name IS NULL` rows when it is `NULL`). A block that is passed **and** present in the resolved rows is an error, as is a passed block already stored anywhere in the table (§6C, checked on the whole table so a filter cannot hide it). E.g. `G_D = …` plus `trait_var_comp_tbl` filtered to the stored additive rows is allowed. |
| (b) As Part A: any passed matrix with `trait_var_comp_tbl` is an error | Simpler; forces the user to store blocks first with `define_effect_cov_matrix()` to mix sources. |

(a) is what the default already does (`trait_var_comp_tbl = NULL` reads the stored rows
while `G_D` is passed), so making the explicit table behave the same is the consistent
choice. Part A keeps its one-block refusal. Part A's "target dependents" refusal
(generated scopes it would *keep* that fall back to a new target) is **not** inherited:
`define_genome_effects()` deletes every generated scope of its traits.

---

## 5a — Writer speed (0.75.3) — *done; see the results file*

**As built.** As specified below, plus the same one-pass grouping in
`.gev_target_kind()` / `.gev_term_line()` / `.gev_term_parent()`, the three
`define_additive_effects()` scope helpers and the evaluator's `.gev_variant_map()` /
`.gev_preflight()`, which had the same per-id pattern and which 5b and 5c reach with large
models (the last two were found by the extractor benchmark below). The contract test landed in the 5a commit (D1 fixes three
commits) after passing against unmodified 0.75.2. `validate_genome_effects()` needed no
change of its own: after `.ge_validate_frames()` it is linear.

### 5a.1 The measurement

Profile of `define_genome_effect_terms()` writing 4,000 `aa_terms()` pairs (300 loci,
`Rprof`, 2026-10-06; 10.0 s total):

Shares are **inclusive** (a function's time includes its callees), so they overlap and do
not sum to 100%.

| Function | Inclusive share | Cause |
|---|---|---|
| `.ge_build()` | 75% (6.0 s alone) | Per-term `terms[idx == k, ]` subsetting, **twice** (term rows and members): `O(terms × rows)`, so quadratic. `[.data.frame` is 52% of the whole profile. |
| `.ge_validate_frames()` (from `.ge_commit()`) | 23% | Per-term loops over `terms$id_genome_effect` and per-member loops (`R/genome_effects_helpers.R:456–474`). |
| `validate_genome_effects()` | 12% | Whole-table SQL plus per-family R checks (`:521–527`). |

That is why 1,000 / 2,000 / 4,000 pairs took 2.2 / 5.0 / 11.5 s and 19,900 took up to 196 s.
`.ge_resolve_deletes(mode = "replace_scope")` has the same pattern
(`model$members[model$members$id_genome_effect == id, ]` per id), which matters to
`define_additive_effects()` re-runs on large models. `.stored_to_functional()` (step 4)
loops per term with named-vector indexing, and it is the 13 s of the step-4 genic
benchmark at 20,000 terms.

### 5a.2 Changes

- **`.ge_build()`**: one `split()` of the rows by term (or an `order()` plus run lengths),
  never a per-term logical subset.
  - Per-term scalar checks become vectorised: one `genome_value` per term, one
    `effect_name` per term, no locus twice in a term. Use `tapply`/`duplicated` on
    `(term, value)` keys, then report the first offending term in user order.
  - Members are built by one `order(idx, locus_id)` and `ave()` for `member_slot`.
  - The error messages and their order are unchanged (D7: term, then rule, then member).
- **`.ge_validate_frames()`**: the same treatment, keyed per term and per member, with
  the returned violation vector unchanged in content and order.
- **`validate_genome_effects()`**: measure after the two above; vectorise only what still
  shows in the profile. It must keep validating the **whole table** before `COMMIT`
  (CLAUDE.md, genome effects). Do not narrow it to the new rows.
- **`.ge_resolve_deletes(replace_scope)`**: split members and origins once.
- **`.stored_to_functional()`**: vectorise per contrast shape (additive, dominance, the
  three indicator states, A×A) with `rowsum()` into the locus index. The cancellation rule
  of 0.75.2 (`.cancelled()`) needs per-coefficient `sum|contrib|` and count, which
  `rowsum()` gives as well.

No change to what is written, to the transaction, or to `.ge_commit()`'s whole-table
validation.

### 5a.3 Gates (5a)

- **A new error-order contract test, written and committed against the 0.75.2 code before
  any refactor**: `tests/testthat/test-genome-effects-writer-order.R`. Hand-authored
  expectations of current behaviour (not golden output):
  - two terms malformed by the same rule;
  - two terms malformed by different rules, where the later rule fails on the earlier term;
  - user `term_id`s whose lexical order differs from input order;
  - several `.ge_validate_frames()` violations at once, asserting the whole message.

  It must pass before and after 5a. It is what makes the "report the last bad term"
  mutation a real gate (Codex review, finding 8).
- **Every existing writer, builder, evaluator and extractor test passes unchanged.** No test
  file is edited in 5a: `test-genome-effects-*.R`, `test-define_genome_effect_terms*.R`,
  `test-genome-effect-terms-builders.R`, `test-extract_genetic_variance.R`. That is the
  correctness gate for a refactor that must not change behaviour.
- **Equivalence during development (not committed).** A scratch script builds random term
  sets (all three contrasts, indicator surfaces, origins, multi-owner) with the old and new
  `.ge_build()` and asserts `identical()` frames and identical error messages on a set of
  bad inputs. CLAUDE.md forbids committed golden-from-old tests, so the script is a
  development check only, and its result is recorded in the results file. It also compares
  the old and new `.stored_to_functional()`:
  - cancellation cases (0.75.2's `.cancelled()`);
  - unusual locus names;
  - multi-owner models.

  The vectorised grouping must be deterministic on repeat. Old floating-point bits are not
  a compatibility promise, but agreement within 1e-12 is expected.
- **Benchmark** `dev/benchmarks/benchmark_genome_effect_writer.R`: write 1,000 / 4,000 /
  16,000 / 64,000 A×A pairs, plus all 124,750 pairs of 500 loci. Record the time per pair at
  each size. **Target:** at most 2× growth in time per pair from 4,000 to 64,000 pairs (near
  linear), and the 124,750-pair write in under 60 s. Re-run the step-4 extractor benchmark
  at the planned 124,750 pairs (it was cut to 19,900 because of the writer).
  - Also time **replacement into a populated model**: re-writing one trait's 16,000 pairs
    when the database already holds two other traits' models, plus scoped and multi-owner
    terms. Whole-table validation cost depends on what is already stored.
  - Record the hardware, DuckDB thread count, trait count, setup vs write time, and peak R
    memory next to every number.

### 5a.4 Docs (5a)

NEWS (internal performance; the numbers), DESCRIPTION 0.75.3, results file
`plans/import_qtl_effect_methods_phase_5.md` (5a section). No skill change unless an
internal named in the skills changes its contract (none should).

---

## 5b — `define_genome_effects()` (0.76.0) — *done 2026-10-09; see the results file*

### 5b.1 Signature (§9.1, unchanged except D4/D5)

```r
define_genome_effects(tbl, trait_name,
  G_A = NULL, G_D = NULL, G_AA = NULL,
  trait_var_comp_tbl = NULL,
  pairs   = NULL,                         # data frame locus_name_1, locus_name_2
  n_pairs = NULL,                         # random matching; NULL = floor(m / 2)
  anchor  = c("genic", "realised"),
  dominance_degree_mean = 0.19,
  dominance_degree_sd   = 0.097,
  inbreeding_depression = NULL,           # named numeric, by trait_name
  base_tbl    = NULL,
  warn_bounds = c(0.8, 1.25))
```

- `tbl` is `get_table(pop, "genome_meta") |> filter(...)`, as `define_additive_effects()`.
  The filtered loci are the QTL of every trait (`"shared"` semantics; there is no
  `method`).
- No `seed` (Q12), no `line_name` / `parent_origin` (common scope only, Q14), no
  `distribution` (the additive architecture is `.draw_additive_architecture(mask,
  "normal", G_A)`), no `effects` (Q21). The roxygen says which generator to use when:
  `define_additive_effects()` for scoped or crossbred additive models, this one when
  dominance or A×A may be added.
- Returns the `tidybreed_pop` invisibly.

### 5b.2 Validation and resolution (§9.2 step 1). No RNG, no write.

In this order, so every **knowable** refusal leaves `.Random.seed` and the database
unchanged (C20, as A22). That holds both when `.Random.seed` already exists and when it does
not. A failure that depends on the draw itself (the chosen pairs' rank, the floor, the
dominance-mean solver, the final verification) comes after the draw. It leaves the database
unchanged but has consumed the draws; nothing in `R/` restores `.Random.seed` (CLAUDE.md).

1. `tbl` is a `genome_meta` table; traits exist and are distinct; scalars and
   `warn_bounds` are valid.
   - `dominance_degree_sd ≥ 0`.
   - `inbreeding_depression` (D4): finite, with distinct, non-missing names, each a call
     trait; aligned by name. It needs `dominance_degree_sd > 0`.
2. **Targets** (§6C, D8), for `effect_name ∈ {additive, dominance, additive_by_additive}`,
   resolved **once, per block**:

   | Block source | Rule |
   |---|---|
   | passed matrix (`G_A` / `G_D` / `G_AA`) | dimnames checked, never overwritten; refused if any row for that `effect_name` × call traits × `line_name IS NULL` exists anywhere in `trait_var_comp` (whole table, so a filter cannot hide it); written in the commit |
   | `trait_var_comp_tbl` given | its rows for the call traits, per `effect_name`; a block also passed is an error |
   | `trait_var_comp_tbl = NULL` | the stored `line_name IS NULL` rows, per `effect_name`, for blocks not passed |

   - Then, per resolved block: one complete symmetric k×k block, no row linking a call
     trait to an outside trait (unless the explicit table chose that), and one candidate set.
     This is the existing `.dae_resolve_target()` logic, generalised to a per-block resolver
     that both generators share. Part A calls it for `additive` only and keeps its
     one-block "`G` and `trait_var_comp_tbl` together" refusal, and its non-additive refusal
     reads the same resolved rows.
   - The per-block resolver keeps `.dae_resolve_target()`'s **scope check** for every block
     (A, D and A×A): explicit `line_name` rows cannot describe the common-scope terms this
     call generates, so a `trait_var_comp_tbl` selecting them is refused (Codex re-review).
   - Part A's target-dependents check is not used: this call deletes every generated scope
     of its traits.
   - Each block is validated by `.qtl_target_std()` (finite, symmetric, PSD on the
     correlation scale).
   - An absent `additive` block is an error.
   - `inbreeding_depression` without a resolved dominance block is an error, and so is a
     named trait whose resolved `G_D[t, t]` is 0 (Codex review, finding 4).
   - **D5:** classify the call. Additive-only route when no D / A×A block exists or every
     present one is the zero matrix. Otherwise `G_A` must be positive definite on its
     correlation scale, or the implementation-limit error.
   - **Degree parameters that cannot carry dominance** (Codex re-review, finding 3).
     A non-zero resolved `G_D` with `dominance_degree_mean = 0` **and**
     `dominance_degree_sd = 0` is refused: every draw then gives `B_d = 0`, which no
     calibration can scale to a non-zero target. The message calls it an
     architecture-parameter restriction (set a non-zero mean or sd), not an infeasible
     target. The same parameters with an absent or zero `G_D` are valid: there is no
     dominance to carry (the call is additive-only unless `G_AA` is non-zero).
3. **Owner content and replacement count** (§5).
   - Read the trait's `generated` terms. Count those this call will delete (all of them,
     `mode = "replace_owner"`) and how many are line-scoped, for the message.
   - Custom owners are untouched.
4. **Loci and base, lightweight** (reordered after the Codex re-review, finding 3: the
   pair and cohort checks below need `m`, the locus keys and `n`). `genome_order`, the
   sorted locus keys and their count `m`, the mask, `assert_qtl_autosomal()` (C10), and
   `.dae_resolve_base(pop, base_tbl, NULL)` with the same Wahlund warning as Part A.
   - A selected QTL with no copies in the base is an error, as in Part A.
   - Under `"realised"`: `.dae_check_realised()` (individuals only), and the base's
     individual count `n` by a `COUNT` query. **No dosages are collected yet.**
5. **Pairs** (§9.1, Q8), whenever an A×A block is resolved, zero or not. The rank
   checks below apply to a non-zero block only.
   - Argument-only checks (shape, column names, types; `pairs` / `n_pairs` without a
     resolved A×A block, or both together) could run in step 1; they are grouped here
     for reading, and nothing between changes RNG or the database.
   - Checks needing the resolved keys: unknown loci, loci outside the filter, self-pairs
     and repeats in either order are errors naming the keys.
   - `n_pairs` is an integer in `1..floor(m/2)`; a larger value errors with the maximum
     and points to `pairs`.
   - **Knowable shortages, before any draw** (Codex review, finding 7):
     - `m < 2` with random pairing;
     - `rank(G_AA) > r`, where `r` is the number of pairs (supplied, `n_pairs`, or
       `floor(m/2)`): r pair columns cannot carry a rank above r under any anchor;
     - under `"realised"`, `rank(G_c) > n − 1` for any block (a centred cohort of n has at
       most n − 1 directions).
6. **Dosage guard, then collection** (`"realised"` only). The guard counts **the designs
   this route actually keeps** (Codex review, finding 7):
   - additive-only route: Part A's guard, `n × m`, with Part A's message;
   - otherwise, `n × (m_A + m_D + r)` cells, where `m_D = m` with a non-zero D block, else
     0, and `r` counts only with a non-zero A×A block.

   The error names `n`, each block's column count, the product and the limit, and says
   the limit bounds the retained design arrays, not peak memory (sweeps and pair products
   allocate temporaries). Only then `.collect_dosages()` and the design arrays.
7. **Anchor feasibility per block**, before the draw, after the designs exist:
   `rank(G_c) ≤ rank(M_c)` for A and D (Part A's first rank error, A5 (a)), and for A×A
   when the pairs are supplied.
   - **Rank convention** (Codex review, finding 5): the existing one, singular values
     squared against `max · rel_tol` on the eigenvalue scale (`.qtl_anchor_design()$rank`).
     No QR shortcut: `qr(diag(1, 1e-6), tol = 1e-10)` gives rank 2 where the convention
     gives 1.
   - The resolved rank is cached on the anchor object, so `.qtl_calibrate()` does not
     repeat the decomposition.
   - With random pairs, the A×A anchor rank is checked after the pair draw. The message
     says that the chosen pair design is rank-deficient under the named anchor, and gives
     both ranks.

### 5b.3 Draw (D2 order). The only RNG use.

`.dge_draw(mask, G_A, has_d, has_aa, pairs_spec, dd)` returns every architecture, so C17
can test the draws directly:

1. `B_a ← .draw_additive_architecture(mask, "normal", G_A)` (C4).
2. If D (zero or not): `z ← matrix(rnorm(m·k), m, k)`, `h ← mean + sd · z`,
   `B_d ← h · |B_a|` (source line 223; `|B_a|` is the architecture, before stage 3).
   `rnorm(n, μ, σ)` computes `μ + σ · norm_rand()`, so this is the source's stream. `z` is
   kept and passed to the calibrator, which therefore never divides by `|B_a|` to
   recover degrees (Codex review, finding 4).
3. If A×A (zero or not) and `pairs = NULL`: QTL sorted by `locus_id`, one `sample()`
   permutation, paired off consecutively, the first `n_pairs` kept. Each pair is put in
   `aa_terms()`'s C-locale order, then pairs are sorted by `(locus_id_1, locus_id_2)`
   (§9.1). Supplied pairs go through the same canonicalisation, without RNG.
4. If A×A: `B_aa ← matrix(rnorm(r·k), r, k)`.

A zero block is still drawn, so a zero block and a non-zero block of the same kind consume
the same stream (D5's table).

### 5b.4 Calibration — `R/genome_effects_calibration.R` (§9.4)

The source's `sim_qtl_effects_nonadd()` body becomes `.na_calibrate(p, design, G, B_a,
z, B_aa, pairs, anchor, dd, inbreeding_depression)`. It is pure: every architecture,
including the standard-normal degree deviations `z`, is supplied; nothing is drawn, nothing
is written. Only the **non-additive route** (D5) reaches it.

**Anchors as objects, not matrices** (step-4 lesson: no `m×m` or `r×r` dense anchors).
Extend `R/qtl_congruence.R`'s anchor objects with a cross term:

- `.qtl_anchor_diag(w)` gains `cross(B1, B2) = crossprod(B1 * w, B2)`.
- `.qtl_anchor_design(Xc, den)` gains `cross(B1, B2) = crossprod(Xc %*% B1, Xc %*% B2) / den`.

Per block:

| Block | `"genic"` | `"realised"` (design, denominator `n − 1`) |
|---|---|---|
| A | `diag(2pq)` | `Z_A = X − 2p` |
| D | `diag((2pq)²)` | `Z_D = W_c − Z_A diag(b)`, observed `b` (0 at monomorphic loci) |
| A×A | `diag(4 p_k q_k p_l q_l)` | `Z_AA`, the centred products of the pairs' `Z_A` columns |

with `p` from the base, and `b = q − p` (genic) or the observed within-locus regression
(realised).

**Stages** (§1.3, highest order first):

1. **A×A:** `B_aa ← .qtl_calibrate(B_aa, .qtl_target_std(G_AA), M_AA)`.
2. **Dominance:** `B_d = (mean + sd · z) · |B_a|`. For each trait named in
   `inbreeding_depression`, its column becomes `μ_t |B_a| + sd · z|B_a|` with `μ_t` from
   `.na_solve_dd_mean()`; the other columns keep the drawn degrees (D4). Then
   `B_d ← .qtl_calibrate(B_d, .qtl_target_std(G_D), M_D)`.
3. **Additive:** `C = diag(b) B_d + E_c`, where `E_c` accumulates `c_l · e_r` onto each pair's
   partner (vectorised with `rowsum()`, never a loop over pairs). Then the new
   `.na_additive_stage(B_a, C, G_A, M_A)`.

**`.na_solve_dd_mean()` — the source's equation, a robust solver** (Codex review, finding
2; **not** a verbatim port). The source always applies the quadratic formula and divides
by `2A`. It fails on feasible inputs: all coefficients zero (one genic locus, where
`ID/√V_D = sign(d)`, with `ρ = ±1`) gives `missing value where TRUE/FALSE needed`. `A = 0`
with `B ≠ 0` gives roots `NaN, Inf` where the answer is `μ = 0.1` (Codex's probes). The
port keeps the equation and the root policy, and replaces the arithmetic:

Revised again after the Codex re-review (finding 1): the first revision's stable formula
lost both roots at `B = 0`, its reflection could land on `ID = 0` with zero variance, and
its coefficient scale depended on `sd`.

- **Solve in `x = μ / sd`** (`inbreeding_depression` requires `sd > 0`). With
  `u = |B_a[, t]|`, `v = z_t ∘ u` and `w = 2pq`, the degrees are `d = sd (x u + v)`, so
  `ID = sd (x w'u + w'v)` and `V_D = sd² (x² u'Mu + 2x u'Mv + v'Mv)`. The three quadratic
  forms come from `M_D.cov(cbind(u, v))`, a 2×2 from the anchor object; there is no dense
  `M`. The requested ratio `ρ = ID / √V_D` does not depend on `sd`, and neither do the
  coefficients below: `sd²` factors out of `ρ² V_D − ID² = sd² (A x² + B x + C)`, with
  `A = ρ² u'Mu − (w'u)²`, `B = 2 (ρ² u'Mv − w'u · w'v)`, `C = ρ² v'Mv − (w'v)²`.
  The degree mean is `μ = sd · x`.
- **Rounding budgets**, both relative to the operands, not to their cancelled differences:
  - a coefficient is ≈ 0 when `|coef| ≤ tol · s`, with
    `s = ρ² (u'Mu + 2|u'Mv| + v'Mv) + (|w'u| + |w'v|)²`;
  - a negative discriminant `B² − 4AC` is residue (clamped to 0) when
    `≥ −tol · (B² + 4|A||C|)`, and otherwise an error naming the Cauchy–Schwarz bound.
- **Branches:**
  - **all three ≈ 0** (the ratio holds wherever `V_D > 0` and the sign fits; e.g. one
    genic locus, where `ID/√V_D = sign(d)`): candidates in a fixed order, the first
    accepted wins: the preferred `x_p = dominance_degree_mean / sd`; its reflection
    `2 x₀ − x_p` across the root `x₀ = −w'v / w'u` of `ID`; then `x₀ + σ`, one `sd` past
    the root on the side whose `ID` has the sign of `ρ`. The last covers a preferred mean
    already on the root (Codex's `u = 1`, `v = −1.9`, `sd = 0.1`, mean 0.19, where the
    reflection is the same point and `V_D = 0`). With `w'u ≈ 0`, `ID` is constant: the
    wrong sign is an error; with `ρ = 0` the candidates are `x_p`, `x_p + 1`, `x_p − 1`
    (a non-zero `V_D` vanishes at one `x` at most);
  - **`A ≈ 0`, `B` not:** the one linear root `−C / B`;
  - **`A ≈ 0`, `B ≈ 0`, `C` not:** the ratio is a limit no finite mean attains; error,
    naming the bound;
  - **otherwise** the quadratic, with `sign₊(B) = 1` for `B ≥ 0` (never R's `sign()`,
    which is 0 at 0 and lost both roots of Codex's two-locus case: `p = (0.5, 0.5)`,
    `u = (1, 1)`, `v = (1, −1)`, `ρ = 1`, roots `±0.1` in `μ`). `q = −(B + sign₊(B)√disc)/2`;
    if `q = 0` (only when `B = 0` and `disc = 0`, so `C = 0`), the one repeated root is
    `x = 0`; otherwise the roots are `q / A` and `C / q`.
- **Acceptance**, the same for every branch and every candidate:
  - `x` finite;
  - `ID(x)` has the sign of `ρ` (`ρ = 0`: `|ID| ≤ tol` relative to its operands);
  - `V_D(x)` finite and **positive**, beyond `tol · sd² (x² u'Mu + 2|x||u'Mv| + v'Mv)`. A
    zero dominance variance cannot be calibrated to a positive `G_D[t, t]`, whatever its
    `ID`; this applies to `ρ = 0` too.

  Among accepted quadratic roots, the largest, as the source does. No accepted candidate is
  an error saying which condition failed (sign, or zero dominance variance).
- **Verify** `|ID(μ)/√V_D(μ) − ρ| ≤ tol · max(1, |ρ|)` before accepting it. A miss is an
  error, never a silently wrong mean.

**One congruence (§9.4).** Stages 1 and 2 use Part A's `.qtl_calibrate()`, not the source's
`congruence()`. That gives the correlation-scale rank and PSD judgement (a trait in small units
is never truncated), and verification against the target as stored at
`QTL_CALIBRATION_TOL`. The source's `congruence()` and `.qtl_congruence()` have the same
algebra (§9.4: "port it once").

**`.na_additive_stage()`** ports `solve_additive_stage()` with three changes:

- **Correlation scale, with every quantity named in its coordinates** (Codex re-review,
  finding 2). With `S = diag(√diag G_A)` (invertible on this route by D5), the stage works
  entirely in standardised coordinates; nothing is rescaled twice:

  ```text
  R_A      = S⁻¹ G_A S⁻¹
  B_s      = B_a S⁻¹ ;  C_s = C S⁻¹
  P_s      = M.cov(B_s) ;  Q_s = M.cross(B_s, C_s)        # anchor objects, no dense M
  C_res_s  = C_s − B_s solve(P_s, Q_s)                    # coupling the span cannot cancel
  floor_s  = M.cov(C_res_s)                               # already on the correlation scale
  G̃_s      = R_A − floor_s
  floor    = S floor_s S                                  # original units: messages only
  ```

  The solved standardised coefficients are rescaled to original units **once**, at the end.
  `floor_s` equals the source's `R − QᵀP⁻¹Q` in these coordinates, but it is PSD by
  construction, without subtracting two large covariances (Codex review, finding 5).
  Computing the residual in original units would also be valid; this plan does not.
- **Floor test with a rounding budget.** The source judged `min eig(G̃)` against
  `max |eig(G̃)|`. At an exactly feasible target `G̃ ≈ 0`, so round-off looked decisively
  negative: it refused `G_A = 12.3077497977446`, the independently computed floor, at
  −1.78e-15. Here a negative eigenvalue of `G̃_s` is residue when
  `≥ −tol · max(1, ‖floor_s‖)`, a budget relative to the operands, not to their cancelled
  difference. Residue is clamped to 0, and the final delivered-target check (below) still
  applies.
- **The floor is conditional on the sampled architecture** (Codex re-review). Its diagonal
  is the least additive variance reachable **within the drawn additive span, given the
  calibrated dominance and A×A effects**: not a biological minimum fixed by `G_D` / `G_AA`
  alone. A different draw gives a different floor, and with `C = 0` (e.g. every base
  `p = 0.5` under `"genic"`) it is zero however large `G_D` is. Errors, messages and the roxygen say
  "for this sampled architecture".
- **Errors:**
  - Below the floor (C5): names the floor's diagonal in original units (`diag(floor)`,
    the minimum additive variance per trait **for this sampled architecture**), the
    smallest eigenvalue of `G̃_s`, and "raise `G_A` or lower `G_D` / `G_AA`".
  - Rank-deficient `P` (source test 10): unreachable after D5 except numerically, and then
    it names the architecture.
  - Every error names the anchor: a target feasible under `"genic"` can be below a
    realised cohort's floor.

The `T` choice that maximises `tr(T)` (the larger root for k = 1) is kept.

**The additive-only route (C4, D5).** With no D / A×A block, or only zero ones, the
generator takes **Part A's whole additive path, shared as one helper both generators call**
(Codex review, finding 3):

- the anchor and its guard;
- `.qtl_calibrate(B_a, std_A, M_A)`;
- **Part A's writer frames** (`.dae_build()` per trait, stacked with `.dae_stack()`).

It does not call the quadratic with `C = 0`. For k ≥ 2 the quadratic's orthogonal choice
(maximise `tr T`) and the congruence's polar factor give different, equally exact `B`.
Storage matters as much as calibration: `.dae_build()` writes a zero-valued term at every
selected locus of a zero-variance trait, where `.noia_terms()` would drop it, and `.noia_terms()`
errors on an all-zero model. Only the one shared path makes C4 (b) row-identical, zeros
included. The zero D / A×A blocks are resolved, written when passed, and drawn (5b.3), and
they produce no terms.

**Verification.** After stage 3, every present block's delivered covariance
(`M_A.cov(B_alpha)`, `M_D.cov(B_d)`, `M_AA.cov(B_aa)`) is checked with
`.qtl_target_error() ≤ QTL_CALIBRATION_TOL` against the target as stored. A miss is an
error before anything is written, as in Part A.

### 5b.5 Storage (§9.2 steps 4–5)

- **Additive-only route:** Part A's frames, above.
- **Non-additive route**, per trait: `stat ← .noia_to_stored(a = B_a[, t], d = B_d[, t],
  pairs(e = B_aa[, t]), p)`, then `.noia_terms(stat)`. These are Cockerham `ad_terms()` +
  `aa_terms()` rows at `p = p_base`, keyed by locus name. A trait on this route never
  builds an empty model. The argument is about the model, not about α (Codex re-review): D5
  makes `G_A` positive definite, so every trait's delivered additive variance under the
  anchor is positive, and the verified model therefore has non-constant genetic values. If
  every stored α, d and e of a trait were zero, every genotype would have the same value and
  the delivered additive component would be zero. (The stored `stat$alpha` is the
  HWE-referenced re-derivation, which under `"realised"` differs from the calibrated
  `B_alpha` and need not itself be non-zero.) `.noia_terms()`'s all-zero error stays an
  internal assertion. `.noia_to_stored()` loops over
  pairs; profile it in 5b at the benchmark size and vectorise with `rowsum()` only if it
  shows (Codex review).
  - Under `"genic"`, `stat$alpha` equals the calibrated `B_alpha` (`b = q − p` in both).
    The test asserts this at 1e-12.
  - Under `"realised"`, `.noia_to_stored()` re-derives α with `b = q − p` (§4.3.1). The
    total is exact, and the stored `additive` / `dominance` split is HWE-referenced (the
    4.134-vs-4 effect of §4.2). `extract_genetic_variance(anchor = "realised")` on the base
    individuals recovers the targets exactly (C20), because it re-projects with the
    observed `b`.
- On the non-additive route, the builders drop exact zeros (`ad_terms()` α = d = 0 rows,
  `aa_terms()` e = 0 pairs). Under D3 (a) a fixed locus keeps its non-zero drawn effects,
  so only true zeros go: a **zero D or A×A** target column (e.g. `G_D = diag(c(1, 0))`)
  writes no dominance terms for that trait. The zero-row policy for the **additive** block
  is Part A's (above).
- One `.ge_commit()`. Delete every `generated` term of every call trait (`replace_owner`);
  insert all traits' rows; `before_commit` writes each **passed** target block with
  `.tvc_write_block()` (up to three). One transaction (C11).
- `.ge_build()` is called per trait and stacked with `.dae_stack()`, as in Part A. Every
  route builds at least one term per call trait, so no empty build reaches
  `.dae_stack()` / `.ge_commit()`. If a later change makes an empty model possible, it needs
  a typed empty build and a delete-plus-target commit first.

### 5b.6 Owner rules elsewhere (§5)

- **`define_additive_effects()` refuses** when a call trait's `generated` terms include
  anything but order-one `additive` terms. The check comes after argument validation and
  before target resolution, any write, and `set.seed()`. The error names the fix:

  ```r
  get_table(pop, "genome_meta") |> filter(...) |>
    define_genome_effects("ADG", trait_var_comp_tbl = get_table(pop, "trait_var_comp") |>
      filter(effect_name == "additive", is.na(line_name)))
  ```

  It also says the later `define_additive_effects()` call needs the same filter while the
  `dominance` / `additive_by_additive` targets remain stored (C8 (b)).
- `define_additive_effects()`'s stored-non-additive-target error (§6C) gains
  `define_genome_effects()` as a third fix (A16).
- `.ge_resolve_deletes()`'s `replace_trait` refusal message names both generators.
- `remove_generated_effects()`'s roxygen names `define_genome_effects()` models (removed
  whole by `remove_generated_effects(pop, trait)`).
- **The rank note (5b.7 item 8) is one shared helper** (`.qtl_rank_note()`), also called by
  `define_additive_effects()` for its additive target and by `define_effect_cov_matrix()`
  when it stores a block, so a singular target is reported wherever it enters. No change to
  what is accepted or refused.
- **Grep gate R1 (step-5 form)**, run before 5b starts: nothing in `R/`, `man/`, `tests/`
  or `vignettes/` names `define_genome_effects()` (only the API skill's forward reference,
  which is updated in 5b).

### 5b.7 Messages (§9.2 step 6). Reported, not stored.

1. Replaced: "Replaced N generated terms of trait T (L of them line-scoped: crossbred
   variants from define_additive_effects())". Only when N > 0, and always naming L when
   L > 0 (C8 (a)).
2. Random pairs: the §9.1 text with the drawn count and `floor(m / 2)`.
3. Delivered per block under the anchor ("exact under the genic anchor"), as Part A's
   closing message, with covariance text per present block.
4. The additive floor's diagonal in original units, and how far `G_A` is above it,
   labelled as conditional on the sampled architecture.
5. `inbreeding_depression`: requested vs delivered by trait name, with "approximate for
   k ≥ 2, no closeness guarantee" (D4). Without the argument, the implied depression `Σ 2pq d` per trait is still reported
   (§9.2 step 6). Sign convention: the drop in the mean per unit of inbreeding F, positive
   = depression, as the source's `inbr_depr`.
6. A one-line functional summary per trait: the number of QTL with non-zero `a`, with
   non-zero `d`, and the number of pairs, plus the **sampling parameters used**
   (`dominance_degree_mean`, `dominance_degree_sd`, and each solved `μ_t`). No "delivered
   dominance degree" is reported. The degrees scale the pre-calibration `|B_a|`, and stage
   3 changes the additive effects, so `d / |a_final|` (with its zero denominators) would
   not be the parameter the user set (Codex review).
7. The D6 diagnostics.
8. **Rank note for singular targets** (user decision 2026-10-09, with D5: transparency,
   not a `warning()`). For every resolved block that is singular on its correlation scale
   (the same rank convention as validation), one `message()` per block naming the block,
   its rank and the reason, most specific first:
   - traits with **zero variance** in that block, by name ("FCR has zero dominance
     variance");
   - trait pairs with genetic correlation **±1** (within the rank tolerance), by name;
   - otherwise "a linear dependency among traits …", naming the traits in the null space.

   It is informational: singular targets are valid and, on the additive-only route,
   delivered exactly. It catches a typed 1 meant as 0.99 without refusing a deliberate
   design. A zero block (D5's exception) gets no note; it is reported as absent-by-zero in
   message 3.

### 5b.8 Roxygen content

- **Title/description:** samples and calibrates A (+ D, + A×A) to stored targets, writes
  under `"generated"`, replaces the trait's whole generated model.
- **When the result is exact:** every present block, under the named anchor, to
  `QTL_CALIBRATION_TOL` on the correlation scale; checked, never assumed.
- **The additive floor**, in words: dominance and epistasis already imply additive
  variance; `G_A` below it is refused, and the error names the minimum. The floor belongs
  to the sampled architecture (another seed gives another floor; with every base
  `p = 0.5` under `"genic"` it is zero), not to `G_D` / `G_AA` alone.
- **Targets** (§6C): passed vs stored, `trait_var_comp_tbl`, absent vs zero block,
  never overwritten. The additive-only route, and "re-run to add a block later".
- **Pairs:** supplied (hubs allowed) vs random matching (each locus once, AlphaSimR); the
  `n_pairs` default.
- **The anchor**, and that `"realised"` stores HWE-referenced contrasts (§4.3.1): the total
  is exact; `extract_genetic_variance(anchor = "realised")` on the same individuals gives
  the targets back.
- **Crossbreeding** (§9.1): common scope still gives heterosis through frequency
  differences; the default base pools the founder table, so the targets hold at **pooled**
  frequencies (Wahlund warning); line-specific non-additive effects are Q14, not this
  release.
- **What a re-run replaces** (§5), including line-scoped additive variants, and that
  `define_additive_effects()` then refuses until an additive-only re-run.
- **Generated means calibrated** (Q21): exact coefficients go through
  `define_genome_effect_terms()` with `ad_terms()` / `aa_terms()` (functional coding, Q17);
  the §9.3 worked example.
- **Inbreeding depression:** exact for one trait, approximate and reported for several (D4).
- **The meaning of `additive` under LD** (step-4 follow-up): calibration is exact for the
  anchor; a realised cohort's `additive` row is a contrast component, not a joint
  regression.

### 5b.9 Gates (5b) — `tests/testthat/test-define_genome_effects.R` and `tests/testthat/test-genome-effects-calibration.R`

**Internals (`test-genome-effects-calibration.R`, no database):** the source suite ported
(§9.4). The source's RNG happens inside `sim_qtl_effects_nonadd()`, so a test helper
reproduces its draws and passes them to `.na_calibrate()` as supplied architectures (the
source's test 18 path).

**Oracle isolation** (Codex review; dependency list completed after the re-review,
finding 5). `tests/testthat/helper-nonadd-generator-oracle.R` copies the **transitive
dependency set** of the source entry points the tests call, with the commit (`318e54f`) and
file recorded:

- entry points: `sim_qtl_effects_nonadd()`, `nonadd_covariates()`, `congruence()`,
  `solve_additive_stage()`, `solve_dd_mean()`, `zeng_appendix_A()`;
- their dependencies: `nonadd_decompose()` (called unconditionally for the founder
  diagnostics) and the source's `.na_*` helpers (`.na_max_abs()`,
  `.na_validate_numeric_matrix()`, `.na_validate_scalar()`, `.na_validate_exact_dim()`,
  `.na_psd_eigen()`, `.na_target()`, `.na_relative_spectrum()`, `.na_spectrum_text()`,
  `.na_trace_ratio()`). The copy is checked against `codetools::findGlobals()` of the
  entry points when it is made, not by hand.

They are evaluated in a dedicated environment whose parent is `baseenv()`, with the base
and recommended-package functions they use (`stats`, `MASS` if any) reached explicitly,
never the package namespace or `.GlobalEnv`. So the package's own `.na_*` internals cannot
shadow them, and they cannot silently find a package function. A smoke test runs a small
A + D + A×A `sim_qtl_effects_nonadd()` call without sourcing the checkout, and asserts that
`environment()` of every source-owned function it reaches is the oracle environment.

Production intentionally diverges: correlation-scale congruence, D3 eligibility, the
robust solver. So the comparisons are on what must agree — delivered covariances,
feasibility, genotype-value identities — not on coefficients where the two may
legitimately differ. Hand-derived cases cover the source's own defects (the solver, the
floor boundary). The supplied-architecture and hand-derived cases are the substantive
correctness gates; the oracle is a cross-check.

| Source test | Port |
|---|---|
| 1 exactness, three blocks, both anchors, k = 1, 2, 3 | as is (1e-10) |
| 2 k = 1 genic vs Zeng Appendix A | as C6 (1e-12; `zeng_appendix_A()` copied as oracle) |
| 3 below the floor | as C5 |
| 4 no D, no A×A reduces to the congruence | as C4 (a) |
| 5 inbreeding depression k = 1 | as C7, plus the D4 k = 2 "reported, not exact" case, plus G4 |
| 6 fixed loci | restated under D3 (a): zero variance under the anchor (or as written, if D3 (b)) |
| 7, 8 package-style scaling, SS_NOIA | **not ported**: comparison methods of the source's paper, not part of the generator |
| 9 reduced-rank architecture | as is, through `.qtl_calibrate()` |
| 10 rank-deficient `B_a` | as the D5 refusal at the public level; the internal error kept |
| 11 scale invariance, 24 orders of magnitude | as is, kept as the global-scale sweep; **not** sufficient for per-trait units (it scales every block by one scalar, which the raw-scale source already passes), so see G6 |
| 12 randomised sweep | as is, a few seeds |
| 13 all-heterozygous loci | restated under D3 |
| 14 zero-rank targets | as is |
| 15 names | as is |
| 16 mismatch warning | as the D6 diagnostics, at the public level |
| 17 edge sizes | as is (m = 1 with A + D; one pair) |
| 18 supplied architectures | as is (the internal API is the supplied-architecture path) |

**Public (`test-define_genome_effects.R`):**

- C1–C3, C9, C14, C15, C17, C18 as §11 states, with C17 extended by D2: the dominance draw
  is unchanged by adding A×A.
- C4 (a) internals; (b) public row identity with `define_additive_effects()`, k = 1, 2,
  both anchors, after dropping `id_*` columns and ordering by `(trait_name, locus_id)`.
  Include `G_A = 0`, `G_A = diag(c(1, 0))` and a feasible rank-one target, comparing rows
  **including the zero-valued terms**; replacement of an existing generated model; target
  persistence; and custom owners untouched (Codex review, finding 3). Also the D5
  exception: the same targets plus explicit zero `G_D` / `G_AA` give the same additive rows
  and no dominance or A×A terms.
- C5 floor error with the floor's diagonal in the message; C6; C7 per D4.
- C8 (a)–(d), including the message counts and the refusal in `define_additive_effects()`.
- C10 non-diploid / non-autosomal refused.
- C11 atomicity: trait 2 below the floor leaves trait 1's prior model and `trait_var_comp`
  intact. **Injected failure** inside the transaction (a mocked `.tvc_write_block()` that
  fails on the second of two passed blocks, after the first is inserted) restores all
  three genome-effect tables and `trait_var_comp`, including prior scoped rows and custom
  owners (Codex review).
- C12 (a) exact genotype-count fixture; (b) individual values against the test's own
  computation.
- **C20** (refusals): every knowable validation error leaves `.Random.seed` and the
  database unchanged. Use `expect_identical()` on `.Random.seed` before and after, with and
  without a pre-existing seed; without one, assert `.Random.seed` is still absent
  (`exists(".Random.seed", globalenv())`), never read an undefined variable. A
  draw-dependent failure leaves the database unchanged (and is not expected to preserve
  the seed). The size guard's message is in G7. Includes the zero-degree refusal (5b.2
  step 2): positive `G_D` with `dominance_degree_mean = 0` and `dominance_degree_sd = 0`
  refuses, both with and without a seed; explicit zero `G_D` with the same parameters
  succeeds.
- **New G1 (D5):**
  - a singular `G_A` (correlation 1) with no non-additive block succeeds;
  - the same `G_A` with explicit zero `G_D` and `G_AA` succeeds (additive-only route);
  - the same `G_A` with a non-zero `G_D` errors before any draw, with the
    **implementation-limit** wording (the test asserts the message does not say
    "infeasible" or "impossible");
  - Codex's hand-derived feasible case (`p = 0.5`, `B_a = [1 1; −1 −1]`, `B_d = I`) is
    pinned as such an intentional refusal, with a comment that it is attainable.
- **G1, rank note:** `expect_message()` names the zero-variance trait for
  `G_D = diag(c(0.4, 0))`, the trait pair for a correlation of 1, and the dependent traits
  for a three-trait rank-two `G_A` (`T3 = T1 + T2`); a positive-definite target emits no
  note; the same notes from `define_additive_effects()` and `define_effect_cov_matrix()`.
- **New G2 (determinism, CLAUDE.md):** `expect_identical()` stored rows across two seeded
  runs, under `SET threads = 1` and `8`, and after `restore_pop()` with the same
  `base_tbl` filter (as A14).
- **New G4 (dominance-mean solver, Codex finding 2 and re-review finding 1):**
  - Codex's first two probes: all-zero coefficients with `ρ = ±1`, and `A = 0` with answer
    `μ = 0.1`;
  - the `B = 0` two-locus genic case (`p = (0.5, 0.5)`, `u = (1, 1)`, `v = (1, −1)`,
    `M_D = 0.25 I`, `sd = 0.1`): `ρ = 1` gives `μ = 0.1`, `ρ = −1` gives `μ = −0.1`;
  - the degenerate case with the preferred mean on `ID = 0` (`u = 1`, `v = −1.9`,
    `sd = 0.1`, mean 0.19): `ρ = 1` gives a mean above 0.19 (0.29 under the fixed
    displacement), `ρ = −1` one below (0.09);
  - a repeated root at `x = 0` (`B = 0`, `disc = 0`) with positive `V_D` there;
  - `ρ = 0` where the only root has `V_D = 0`: refused as zero dominance variance;
  - `ρ = 0`, negative `ρ`, and `ρ` just inside and just outside the attainable bound;
  - an unattainable one-locus request;
  - an `sd` sweep (1e-8 … 1e8) of one supplied architecture: the same branch and the same
    `x = μ / sd` at every `sd`.

  Each asserts finite coefficients, finite **positive** `V_D` and the delivered signed
  ratio, not just a successful return. (All but the sweep's exact values were checked
  against a prototype of the spec when the plan was revised.)
  Plus a public one-locus A + D call with `inbreeding_depression`.
- **New G5 (subsets, Codex finding 4):**
  - one requested trait in a two-trait model: the other column keeps its drawn degrees;
  - reversed name order gives identical rows;
  - duplicate or unknown names are refused;
  - a requested trait with `G_D[t, t] = 0` is refused before any draw;
  - every deterministic refusal preserves `.Random.seed`, both when it existed beforehand
    and when it did not.
- **New G6 (per-trait units, Codex finding 5):**
  - two-trait fixtures under `S = diag(1e-6, 1e6)`, with targets and supplied architectures
    transformed consistently: per-trait correlation-scale errors, the small-unit trait
    preserved, and the same feasibility decision as the untransformed problem;
  - targets **on** the floor (computed independently from residual coefficients), just
    above it, and meaningfully below it: accept, accept, refuse;
  - a one-trait case with a non-unit variance (`G_A = 4`) and a supplied coupling whose
    hand-computed original-unit floor is 1: the standardised residual target is
    `1 − 1/4 = 0.75` (a twice-standardised floor would give 0.9375); assert the reported
    floor (1, original units) and the delivered covariance (Codex re-review, finding 2).

  This is the gate that kills the "floor test on the raw scale" mutation.
- **New G7 (resources and knowable shortages, Codex finding 7 and re-review finding 4):**
  - the fixture lowers `QTL_REALISED_MAX_CELLS` with
    `testthat::local_mocked_bindings(QTL_REALISED_MAX_CELLS = …, .package = "tidybreed")`
    (checked: it rebinds the non-function constant and restores it on exit); no test
    allocates ~2e7 cells, and memory is measured in the benchmark only;
  - an additive-only realised call just under the (lowered) `n × m` limit **succeeds**, and
    is compared with `define_additive_effects()` from the same starting seed and inputs:
    identical rows **and identical post-call `.Random.seed`** (successful sampling advances
    RNG, by the same amount in both);
  - the deterministic refusals, each preserving the database and `.Random.seed` (and its
    absence when no seed existed): A + D + A×A over the limit at the new count, naming
    every block's column count; default pairing with one QTL; a rank-two `G_AA` with one
    supplied pair; a realised `rank(G_c) > n − 1`.
- **New G8 (D8 sources and scope, Codex re-review):**
  - passed `G_D` + `trait_var_comp_tbl` selecting only the stored additive rows: allowed;
  - passed `G_A` + stored `G_D` / `G_AA` (default `NULL` table): allowed;
  - a passed block also present in the explicit selection: refused;
  - a passed block already stored but filtered out of `trait_var_comp_tbl`: refused (the
    whole-table check);
  - an intentionally filtered-out stored block (e.g. stored `G_AA` excluded): that block is
    absent for the call, and its stored rows are left untouched;
  - explicit `line_name` rows for any of A, D or A×A: refused by the scope check.
- **New G3 (storage identity):** under `"genic"`, `.noia_to_stored()`'s α equals the
  calibrated `B_alpha` (1e-12), and the stored terms evaluated by `add_tgv()` equal
  functional `g` minus `μ` per individual (C2 is the all-genotype version).

### 5b.10 Docs (5b)

- **Skills:**
  - `tidybreed-api`: a new `define_genome_effects()` section; the `define_additive_effects()`
    refusal; the calibration internals.
  - `tidybreed-schema`: no DDL; the owner paragraph mentions both generators.
- **`CLAUDE.md`:**
  - the Q21 hard rule names both generators (it already says "generators");
  - the two-layer paragraph says targets also enter through `G_A` / `G_D` / `G_AA`;
  - the function table's examples gain `define_genome_effects()`.
- **Package docs:**
  - `_pkgdown.yml`: next to `define_additive_effects`.
  - `package_summary.md`: 46 → 47 exports.
  - `README.md` if it lists generators.
  - NEWS, DESCRIPTION 0.76.0.
- **Plans:** main plan corrections from D2/D3/D4/D5 (§1.3, §9.1, §9.2 step 3, C7, C17).

---

## 5c — End to end and the vignette (0.76.2)

### 5c.1 Gates — mostly `tests/testthat/test-define_genome_effects-integration.R`

- **C13 (phenotypes):** one trait, `G_D` large relative to a feasible `G_A`, no A×A.
  - `add_phenotype()` phenotypic variance minus the residual tracks `V_A + V_D` (not `V_A`).
  - Use a large enough cohort and a tolerance derived from the sampling variance of a
    variance; or, deterministically: the phenotype minus its residual equals
    `ind_tgv_total` per individual, and `var(ind_tgv_total)` equals the extractor's realised
    `total`.
  - Prefer the deterministic form; keep one coarse stochastic check.
- **Prevalence (step-3 review item):** a **generated** A + D + A×A model on a `prevalence`
  phenotype.
  - (a) At the generation base, under `"genic"` with an exact-HWE + LE genotype fixture, the
    extractor's genic blocks equal the stored targets. `between_components` is 0. The
    threshold's summed diagonal equals the realised total **as a population variance**,
    `(n − 1)/n · V_realised`. The extractor keeps its `n − 1` sample convention, so the
    realised sample total is `n/(n − 1)` times the genic sum even at exact equilibrium
    (Codex review, finding 6: `X = (0, 1, 1, 2)`, target 1, sample variance 4/3). C19's
    comparison of two sample covariances needs no factor.
  - (b) Off equilibrium (an LD or inbred cohort), the realised `total` differs from the
    summed targets by `between_components`, by drift in each block, and by the
    `n/(n − 1)` factor. The extractor shows the first two (the sum-plus-cross-covariance
    identity of B3).
  - The threshold still uses the stored sum. This documents the assumption; it is not an
    error.
  - (c) The one-parent-scope rule: a `define_genome_effects()` model is common-scope, so it
    never trips the rule. A test pins that a generated D / A×A model plus a parent-scoped
    additive variant from `define_additive_effects()` cannot exist. The refusal of 5b.6
    prevents it.
- **`remove_generated_effects()`** on a real `define_genome_effects()` model removes every
  term (additive, dominance, A×A), leaves targets and custom terms, and a later
  `define_additive_effects()` then works without the 5b.6 refusal.
- **C19:** `anchor = "genic"` model; `extract_genetic_variance(anchor = "genic", base_tbl =
  <the generation base>)` blocks equal the targets to 1e-10. On an exact-HWE + LE fixture
  at the stored centres, the `ind_tgv` `additive` covariance equals the realised `additive`
  block (1e-10). One fixture away from those conditions shows the documented difference.
- **C20 (round trip):** `anchor = "realised"` model; `extract_genetic_variance(anchor =
  "realised")` on the same base individuals returns `G_A`, `G_D`, `G_AA` (1e-10).
  `warn_bounds` fires on an inbred panel; `NULL` silences it.
- **Generator benchmark** `dev/benchmarks/benchmark_define_genome_effects.R`: 2,000 base
  individuals, 1,000 QTL, k = 2, A + D + A×A with `floor(m/2)` random pairs, both anchors.
  Time and peak memory. Plus one hub design with 20,000 supplied pairs under `"genic"`.

### 5c.2 The vignette — `vignettes/genetic-models.Rmd`

Short. The four paths of the Codex review's "recommended first-release contract":

1. **The writer:** `define_genome_effect_terms()` with `ad_terms()` / `aa_terms()` /
   `genotype_terms()` for known coefficients (§9.3 worked example, functional coding).
2. **The additive generator:** `define_additive_effects()`, exact `G`, and lines.
3. **The genome generator:** `define_genome_effects()` with A + D + A×A, and the floor.
4. **The measuring instrument:** `extract_genetic_variance()` and the target join
   (`inner_join` / `anti_join`).

The scope promise (main plan, step 5): the result is exact for a feasible `G` under the named
anchor and QTL set; a line call calibrates its own variant only; `"union"` does not hit
non-zero off-diagonals; one call has one `parent_origin` scope; the extractor measures the
population the user means.

Runs in seconds: small genome, `:memory:` database, `eval = TRUE` chunks. Listed in
`_pkgdown.yml` articles.

### 5c.3 Docs (5c)

NEWS, DESCRIPTION 0.76.2, results file (5c), main plan "As built (5)", the API skill's
cross-references to the vignette.

---

## Out of scope (stays later)

- Line-specific and parent-of-origin non-additive effects (Q14). The "common variant
  required" rule comes with them.
- Architecture samplers and `marginal =` (Q5); supplied architectures at the public level.
- A whole-model Cockerham builder (Q17 (b)).
- `extract_breeding_values()` and its estimator choice under LD (main plan §4.3.3, as
  qualified in 0.75.2).
- Removing `seed` from `define_additive_effects()` (Q12's package-wide question).
- Evaluation parameters vs generation targets (Q20).

## Risks and things to watch

- **C4 row identity is brittle by design.** It needs the same draw call, the same
  `.qtl_calibrate()` call on the same `B0` rows in the same order, the same `p`, and the same
  term columns (`effect_name` `NA`, `center_value` `p`, no indicator fields). Any
  "harmless" reordering in either generator breaks it. Keep the additive-only dispatch as
  one shared helper both generators call.
- **The realised dosage guard** counts only the designs a route keeps (5b.2 step 6): the
  additive-only route keeps Part A's `n × m`. It bounds the retained arrays, not peak
  memory; the benchmark records the peak.
- **SVD cost of the realised A×A rank check:** `svd(Z_AA)` is `O(n · r · min(n, r))`. At
  n = 2,000 and r = 10,000 that is seconds to tens of seconds. The rank convention is fixed
  (5b.2 step 7), and QR and Gram-matrix shortcuts disagree with it on near dependencies. So
  the SVD stays, computed once and cached on the anchor object. A faster method is adopted
  only with a proof (and a gate) that it agrees on exact and near dependencies.
- **Floor feasibility under `"realised"` depends on the base cohort's LD.** A target
  feasible under `"genic"` can be below the realised floor. The error must name the anchor.
- **5a must not change behaviour.** The equivalence script, plus every unedited test file,
  is the guard. Do not tune `validate_genome_effects()` by validating less.
- **Prevalence on a generated A + D + A×A model** trusts the summed targets. Gate 5c (b)
  shows where that differs from the cohort; the roxygen of `define_phenotype(prevalence = )`
  already says so (step 3).

## Verification (each sub-step)

1. `devtools::document()` and `pkgdown::check_pkgdown()` clean.
2. Focused files first, then the full suite with `NOT_CRAN=true`: 0 failures. The warning
   count stays at the 0.75.2 baseline (6, in the same five files) unless a new warning is
   intended and tested.
3. **Mutation checks** (scripted, file restored in `finally`), each of which must fail a
   gate:

   | Sub-step | Mutation | Must fail |
   |---|---|---|
   | 5a | Report the last bad term instead of the first | `test-genome-effects-writer-order.R` |
   | 5a | Run each rule across all terms before the next rule | `test-genome-effects-writer-order.R` (different-rule case) |
   | 5a | Drop the duplicate-locus-in-term check | the writer tests |
   | 5b | Drop `E_c` from `C` | C1, C2, C3 |
   | 5b | Use observed `b` in storage under `"realised"` | C2 (total no longer `g − μ`) |
   | 5b | Floor test on the raw scale | G6 (per-trait units) |
   | 5b | Floor as `R − QᵀP⁻¹Q` with the source's relative-to-difference tolerance | G6 (on the floor) |
   | 5b | Verbatim `solve_dd_mean()` | G4 |
   | 5b | Additive-only storage through `.noia_terms()` | C4 (zero targets) |
   | 5b | Solve unrequested traits' degree means too | G5 |
   | 5b | D5 check without the zero-block exception | G1 |
   | 5b | Draw pairs before dominance | C17 (D2 half) |
   | 5b | Additive-only through the quadratic instead of `.qtl_calibrate()` | C4 (b) at k = 2 |
   | 5b | Skip the D5 check | G1 (non-zero `G_D` case) |
   | 5b | `define_additive_effects()` content refusal removed | C8 (b) |
   | 5c | Extractor re-projection with `b = q − p` under `"realised"` | C20 |

4. Benchmarks recorded in `plans/import_qtl_effect_methods_phase_5.md`.
5. Commit and push each sub-step with DESCRIPTION and NEWS bumped, directly to `main`.

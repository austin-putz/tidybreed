# Import QTL-effect methods — Step 4 plan (Part B)

**Spec:**
- `plans/import_qtl_effect_methods.md`: §3 (the mapping), §4.3, §8 (Part B), §9.3 (`aa_terms()`
  and the fixed column set), Q9, Q13, Q16, Q17 and §10 Step 4.
- Gates B1–B12 and C16 (§11).

**Versions:** 0.75.0 (4a), 0.75.1 (4b). The main plan gives step 4 one version (0.75.0). This
plan splits it in two (decision D6).
**Status:** planned 2026-10-05, revised the same day after the
[Codex review](import_qtl_effect_methods_phase_4_codex_review.md); not started.
**Starting point:** 0.74.5 (`7c6e2b7`).

**Changes after the Codex review (all eight findings accepted):**

| # | Finding | Where it changed |
|---|---|---|
| 1 | Block rows follow the canonical decomposition, not stored term kinds | 4b.4 "Block availability", B16 |
| 2 | D5's `x`-joined `term_id`s collide (`A`,`B` vs `AxB`) | D5, 4a.1, C16 |
| 3 | A monomorphic pair member still induces an additive effect | 4b.2 rule 5, B12 (main plan too) |
| 4 | The guard bounds `n×m` only; pairs and dense anchors are unbounded | 4b.2 rule 4, 4b.4, B17, benchmark |
| 5 | Default genic base must be the resolved individuals, not `tbl` reused | 4b.1, B18 |
| 6 | `.gev_resolve_traits(NULL)` returns every trait, not traits with terms | 4b.1 |
| 7 | "All generated" is not enough for like-for-like; realised targets are not genic | 4b.5 |
| 8 | "between = total − blocks" is algebraically equal, so it is no mutation | Verification 4 |

The D1–D7 recommendations stand, with D5's encoding replaced and qualifications added to D1
and D7. The main plan's §8 (monomorphic sentence, block-row rule, like-for-like sentence),
§9.3 caveat, Q16 rationale and gates B8/B12 were corrected in the same revision.

Step 4 adds the measuring instrument: `extract_genetic_variance()`. It reports the genetic
covariance blocks of a cohort in `trait_var_comp`'s shape, so that "did the model deliver
its target?" is one `inner_join()`. It also adds the builder and conversion pieces that
both this step and step 5 need: `aa_terms()`, one fixed column set for every term builder,
and the NOIA conversion pair. Nothing in this step writes to the database. There is no DDL
and no stored-string change, so files written by 0.74.x stay readable.

## Decisions for the user (answer before 4a starts)

Each decision has a recommendation. Where the main plan already decided something, it is not
reopened here. These are the points it left open, or ones that reading the code showed.

**D1. How the output names its reference population** (step-4 bullet 1, Codex review item 3).
Every estimate means something different depending on anchor and base. `"genic"` is an
expectation at `base_tbl`'s frequencies. `"realised"` measures the selected cohort. And
`base_tbl = NULL` follows the cohort as it drifts.

| Option | Shape |
|---|---|
| **(a) An `anchor` column on every row, plus one `message()` per call naming the population** ← recommended | `anchor` is `"genic"` or `"realised"`. The message reads, for example, "realised covariances of 240 selected individuals" or "genic expectation at the allele frequencies of `founder_haplotypes` (filtered)". The join keys are unchanged |
| (b) A `reference` text column on every row | Self-describing, but wordy, and the same text repeats on every row |
| (c) Message only | The tibble no longer says what it is once it is saved or passed on |

(a) keeps the tibble honest when it travels, and keeps the join on
`(effect_name, trait_name_1, trait_name_2)`.

**Limitation (Codex review).** `anchor` makes a saved tibble self-describing about the
*estimator*, not about the reference population: it does not record which filtered cohort
or base produced the rows. The message must distinguish cohort-derived frequencies
(`base_tbl = NULL`) from an explicit copy/pool base, and the roxygen says that only the
message carries that. No `reference` column and no stored metadata are added unless saved
provenance becomes a requirement.

**D2. Widen case 1 to every single-locus diploid `indicator`.** §8 accepts only the
heterozygote indicator `(copy_count_value = 2, dosage_value = 1)`. But `1[g = 2]` and
`1[g = 0]` are also exact `(a, d)` plus a constant: `1[g=2] = ((g-1) + 1 - 1[g=1]) / 2`, and
`1[g=0] = 1 - 1[g=1] - 1[g=2]`. **Recommended: accept all three.** The algebra is exact and
costs one line in `.stored_to_functional()`. Without it, a hand-written one-locus genotype
surface (`genotype_terms()` on one locus) would drop to case 3 for no reason. Multi-locus
indicator surfaces stay case 3.

**D3. Case 3 keeps §8's rule** (covered terms projected on NOIA, everything else in
`unpartitioned`), with one consequence to confirm. In a model with common-scope `dominance`
and line-scoped `additive` (a plausible crossbred model before Q14), the scoped additive
family goes to `unpartitioned` (see 4b.3, the family rule). Projecting the dominance terms
on NOIA puts their `(q - p) d` part into the **`additive`** block. So that block is then "the
additive value of the covered terms", not the breeding value of the whole model. That is
correct, and the `decomposition = "partial"` label already says so. **Recommended: keep it,
and say it in the roxygen.** The alternative, reporting the stored `ind_tgv` components
without re-projection whenever a model is partial, would give a coding-dependent answer, which
§8 rules out.

**D4. A test oracle for gate B2.** B2 compares with the source project's
`nonadd_decompose()`. A testthat file cannot `source()` `$SRC`.
**Recommended: copy `nonadd_covariates()` and `nonadd_decompose()` verbatim** into
`tests/testthat/helper-nonadd-oracle.R`. That is about 60 lines of pure functions, with the
source commit (`8f8a97c`) in a header comment. B2 then uses fixed random `(B_a, B_d, B_aa)`.
Calibration does not affect the decomposition, so B2 does not need the source's generator.
The alternative, `skip_if_not(dir.exists($SRC))`, would never run in CI or on any other
machine.

**D5. `genotype_terms()` `term_id`s collide under `rbind()`.** Today it numbers terms
`1, 2, …`. Two surfaces bound together silently share `term_id`s. The result is either a
malformed term or an error far from the cause. C16 promises that the builders bind in any
combination.

~~`term_id = "<locus_1>x<locus_2>…_<row>"`~~ — **rejected by the Codex review.** Locus names
are arbitrary nonempty strings, so the surfaces on `c("A", "B")` and on the single locus
`"AxB"` both get `"AxB_1"`. With equal coefficients the writer accepts the bound frame as
**one three-member indicator term**, silently turning a sum into a product (reproduced
against the current writer). `aa_terms()` pairs `(A, BxC)` and `(AxB, C)` collide the same
way.

**Recommended:** one shared internal encoder, `.term_id(builder, locus_name, suffix)`, used
by every builder. It is a builder prefix plus **length-prefixed** locus names plus a suffix:
`"geno:1:A|1:B#1"`, `"aa:1:A|3:BxC"`, `"ad:1:A#a"`, `"ad:1:A#d"`. Length-prefixing makes the
locus list decode uniquely whatever characters the names contain, and the prefix keeps
builders apart. Locus order is canonical (C-locale radix) where order carries no meaning
(`aa_terms()`; `genotype_terms()` already orders its loci). The ids are deterministic, so
the builders stay pure. No random ids and no session counters.

What this promises (stated in the roxygen): builder outputs never share a `term_id` unless
they describe the same term on the same loci. Overlapping definitions on the same loci can
still be refused by the writer's family rules. That is the writer's policy, not an id
collision. Every builder returns `term_id` as character.

**D6. Two commits.** **Recommended**, rather than step 4 as one 0.75.0 commit.

| Sub-step | Version | Content |
|---|---|---|
| 4a | 0.75.0 | `aa_terms()`, the fixed builder column set (C16), D5, `.stored_to_functional()` / `.noia_to_stored()` |
| 4b | 0.75.1 | `extract_genetic_variance()` and gates B1–B12 |

4b reads 4a's conversion, and its B2 / B10 fixtures use `aa_terms()`. Each commit has the full
suite green and waits for your review, as in step 3.

**D7. `aa_terms()` drops pairs with `e = 0`**, and errors if every `e` is zero. That matches
`ad_terms()`, which drops zero `a` / `d` rows. **Recommended.** The alternative is to keep
zero-valued terms. They evaluate to nothing, but they would still enter the model and make
those loci causal for the trait. Lengths, names, finite coefficients and frequencies, and
duplicate pairs are validated **before** zero pairs are dropped, so a zero `e` never hides
malformed input.

---

## 4a — Builders and the NOIA conversion (0.75.0)

### 4a.1 The fixed column set (§9.3, C16)

Every term builder (`ad_terms()`, `aa_terms()`, `genotype_terms()`) returns exactly these
columns, in this order, with typed `NA` where a column does not apply:

| Column | Type | `ad_terms` additive | `ad_terms` functional `d` | `ad_terms` Cockerham `d` | `aa_terms` | `genotype_terms` |
|---|---|---|---|---|---|---|
| `term_id` | character | `ad:…#a` | `ad:…#d` | `ad:…#d` | `aa:…` | `geno:…#row` (D5) |
| `locus_name` | character | ✓ | ✓ | ✓ | ✓ (2 rows) | ✓ |
| `contrast_name` | character | `additive` | `indicator` | `dominance` | `additive` | `indicator` |
| `center_value` | double | 0.5 / p | NA | p | 0.5 / p_k | NA |
| `copy_count_value` | integer | NA | 2 | NA | NA | ✓ |
| `dosage_value` | integer | NA | 1 | NA | NA | ✓ |
| `genome_value` | double | a / α | d | d | e | value |
| `effect_name` | character | ✓ | ✓ | ✓ | ✓ | ✓ |

- One internal constructor, `.terms_frame(...)`, builds every frame. `.ad_rbind_fill()`,
  which fills with logical `NA`, is deleted.
- The writer already accepts this shape: `.ge_check_member_fields()` requires `NA` centres
  on `indicator` rows and `NA` states on non-indicator rows, and
  `.ge_infer_copy_counts()` fills an omitted indicator copy count where it is unambiguous
  (Codex probe: additive + functional dominance wrote cleanly). Test the fixed typed frame
  against the writer first; relax nothing in the writer unless that test fails.

### 4a.2 `aa_terms()` (§9.3)

`aa_terms(locus_name_1, locus_name_2, e, p_1, p_2, coding = c("functional", "cockerham"),
effect_name = NULL, report = TRUE)`, in `R/genome_effect_terms_builders.R`, exported. It is
implemented exactly as §9.3 specifies:
- Recycling and NA checks reuse `.ad_recycle()`.
- `p_1` / `p_2` are required under both codings.
- A pair whose two loci are the same is refused.
- A pair repeated after canonicalisation is refused.
- Within each pair the two loci are put in canonical C-locale order
  (`order(method = "radix")`), and each `p` moves with its locus.
- `e` is never rescaled. A pair with `e = 0` is dropped, and an all-zero `e` is an error,
  as `ad_terms()` does for `a` and `d` (D7).
- Under functional coding, `report` prints each pair's
  `e (2 p_1 - 1)(2 p_2 - 1)` and their running total, using `ad_terms()`'s message
  wording ("reported only — written to no table").
- The roxygen carries §9.3's caveat: Cockerham `ad_terms()` takes α, and functional `a` is
  not α once pairs exist (Q17).

### 4a.3 The NOIA conversion pair (Q13)

Internal, in `R/genome_effect_terms_builders.R`, next to `ad_terms()`. Both work on one
trait's term set and are pure, with no database access.

- `.stored_to_functional(terms)` takes the evaluator's term frame (as
  `.gev_read_model()` returns it: terms plus members, every member in the common scope)
  and returns functional `(a, d, e)`. `a` and `d` are per locus, and `e` is per canonical
  pair, all keyed by locus identity (`locus_id`, never by position; a locus that appears
  only in a pair or an indicator still gets its `a` / `d` entry). It also returns the dropped constant
  `kappa`, with one sign convention: **stored term = functional term + kappa**. Per term,
  with centre `c` (functional coding uses `g − 1` and `1[g = 1]`):

  | Stored term | Adds to functional | kappa |
  |---|---|---|
  | `additive`, value v | a += v | v(1 − 2c) |
  | `dominance`, value v | d += v; a −= v(1 − 2c) | −v(c² + (1 − c)²) |
  | `indicator (2,1)`, value v | d += v | 0 |
  | `indicator (2,2)`, value v | a += v/2; d −= v/2 (D2) | v/2 |
  | `indicator (2,0)`, value v | a −= v/2; d −= v/2 (D2) | v/2 |
  | `additive × additive`, centres c_k, c_l, value v | e_kl += v; a_k += v(1 − 2c_l); a_l += v(1 − 2c_k) | v(1 − 2c_k)(1 − 2c_l) |

  (All six rows re-derived by hand while answering the Codex review; they agree with its
  table.)

  The dominance row follows from the evaluator's contrast
  (`R/genome_effects_eval.R:794`, `−2c², 2c(1−c), −2(1−c)²`) and the identity
  `1[g=1] = 2pq + (q−p)(g−2p) + x_D` in §3. Every row is checked by a test that evaluates
  both sides on all genotypes (gate N1).
- `.noia_to_stored(a, d, e, pairs, p)` is the forward map of §3. It returns
  α = a + (q−p)d + Σ_l e_jl(2p_l − 1), plus d and e at centres p, plus μ, with the
  convention **functional model = converted statistical model + μ**. A matching round trip
  therefore has `kappa = −μ`, and N2 asserts exactly that sign. It returns **coefficient
  data** (per-locus α, d with centres; per-pair e with centres), not a writer frame. A
  small internal `.noia_terms()` turns that into the fixed builder frame (through
  `ad_terms(coding = "cockerham")` + `aa_terms(coding = "cockerham")`). N2 reads the
  written frame back through `.gev_read_model()` before calling `.stored_to_functional()`,
  so the interface between the two is exercised, not assumed. Step 5's generator is its
  main user. 4a ships it with a round-trip test, because Q13 puts the pair in step 4.
- `ad_terms()` is **not** rewired to call them. Its functional branch computes only the
  single-locus μ, which is `.noia_to_stored()`'s μ with no pairs. The test asserts that the
  two agree. Q13's "`ad_terms()` calls it" was written for a Cockerham `(a, d) → α` path,
  which does not exist (Q17 (a)). Record this as a plan change.

### 4a.4 Gates (4a)

| Gate | Test file | Checks |
|---|---|---|
| C16 | `test-genome-effect-terms-builders.R` (new, or extend the existing builder tests) | §11 C16 as written. Also: `rbind()` of all builder pairs in both orders and both codings gives identical column names and types; two `genotype_terms()` surfaces over different loci bind and write as separate terms. **Collision gates (D5):** `genotype_terms()` on `c("A","B")` + on `"AxB"`; `aa_terms()` pairs `(A, BxC)` + `(AxB, C)`; locus names containing `:`, `|`, `#`; and one locus set used by every builder. Each asserts the **written** term count, each term's member set, and the evaluated `ind_tgv_total` equal to the sum of the separately written surfaces — not merely that the builder ids look distinct |
| N1 | same | Each row of the 4a.3 table: on every genotype in {0,1,2} (pairs: all nine states of {0,1,2}²), the stored term's evaluated value equals its functional form plus `kappa`, to 1e-12. Centres arbitrary and ≠ 0.5, unequal within a pair; two owners in one fixture. Evaluated through `define_genome_effect_terms()` + `add_tgv()` on a small panel, so the check is against the real evaluator, not a re-derivation |
| N2 | same | Round trip: `.stored_to_functional(.noia_to_stored(a, d, e, p))` returns `(a, d, e)` to 1e-12, at random `p`. Writing `.noia_to_stored()` terms and functional `ad_terms()` + `aa_terms()` terms for the same model gives `ind_tgv_total` values that differ by μ alone (1e-10), and `.noia_to_stored()`'s μ equals the source formula in §3 |
| N3 | same | `aa_terms()` on Cockerham coding writes centres p_k, p_l, and the evaluated interaction equals `e (g_k − 2p_k)(g_l − 2p_l)` |

### 4a.5 Docs (4a)

NEWS 0.75.0. Add `aa_terms` to `_pkgdown.yml` after `ad_terms`. Update the API skill's
builder section (the fixed column set, `aa_terms()`, D5) and `package_summary.md`
(45 functions). Fix §9.3's worked example if `rbind()` behaves differently after 4a.
The `aa_terms()` / `ad_terms()` roxygen caveat uses the corrected §9.3 wording: feeding
functional `a` as Cockerham α changes the **genotypic values**, not just the component
split. Example: `a = 0`, `d = 1`, `p = 0.3` — the heterozygote indicator is not the centred
dominance contrast alone; it needs an additive coefficient of `q − p = 0.4`.

---

## 4b — `extract_genetic_variance()` (0.75.1)

### 4b.1 Signature and contract (§8, unchanged)

```r
extract_genetic_variance(tbl, trait_name = NULL, base_tbl = NULL,
                         anchor = c("realised", "genic"))
```

- `tbl` selects the cohort through `resolve_subset_ids(tbl, all_if_null = TRUE)`. Any table
  with `id_ind` works, as for the other action functions.
- `trait_name = NULL` means every trait **with active stored terms**, in `id_trait` order.
  `.gev_resolve_traits(conn, NULL)` does **not** do this: it returns every trait in
  `trait_meta` (Codex review, confirmed in the code). The extractor filters explicitly and
  leaves the shared resolver unchanged, since `add_tgv()` relies on its default. An
  explicitly named trait with no terms errors with its name; a model with no terms at all
  errors too.
- `base_tbl` is genic only. A non-`NULL` `base_tbl` with `"realised"` is an error.
  An explicit `base_tbl` goes through `.validate_base_tbl()`, and its frequencies come from
  `extract_allele_freq()` with that function's existing **allele-copy** semantics (so
  `ind_haplotype |> filter(parent_origin == 1L)` means paternal copies, on purpose).
- `base_tbl = NULL` means the **whole-individual** frequencies of the resolved cohort. It
  never reuses `tbl` as a frequency selection: when `tbl` is `ind_haplotype`, `tbl`'s
  filtered rows pick individuals for `resolve_subset_ids()` but would pick *copies* for
  `extract_allele_freq()` (probe: 0.50 vs 0.35 on the same cohort). The ids are resolved
  once and the frequencies are computed from all copies of those individuals, through a
  registered id view (ids never in SQL text). `extract_allele_freq()` is not changed.
- Every retained covered locus needs a finite, nonmissing frequency. `extract_allele_freq()`
  returns `NA` where the base has no copies, which is not frequency 0 or 1; that is an
  error naming the loci.
- Read-only. Values come from `.gev_evaluate()`, which writes nothing.
- The output is a tibble with columns `effect_name, trait_name_1, trait_name_2, cov_value,
  n_ind, decomposition, anchor` (D1). Trait pairs form a **full square** in both orientations,
  as `.tvc_write_block()` stores targets, so the join matches every stored row. Rows are
  sorted by effect (fixed order: `additive, dominance, additive_by_additive, unpartitioned,
  between_components, total`) and then by trait order.
- The effect name for the A×A block is **`additive_by_additive`** (the §6B vocabulary, as in
  `trait_var_comp`), not the `ind_tgv` component name `interaction`. That is what makes the
  join work. State it in the roxygen.

### 4b.2 Cohort rules (§8)

1. Fewer than two selected individuals: error.
2. **Realised:** evaluate the whole model once with `.gev_evaluate()`. If a selected
   individual is absent for any of the call's traits, the call errors, giving a count per
   trait and asking for a narrower `tbl`. `n_ind` is one number for the result.
3. **Genic:** no evaluation is needed. `n_ind` is the number of selected individuals, which
   matter only when `base_tbl = NULL`, through their frequencies.
4. **Resources.** Port the algebra, not the source's dense anchors: no `m × m` locus
   covariance, no `r × r` pair covariance, no `diag(weights)` (genic forms use the weight
   vectors directly, as `nonadd_decompose()` does for its genic summaries). The arrays
   that remain are:
   - the `n × m` dosage matrix of covered loci, under the existing
     `QTL_REALISED_MAX_CELLS` guard (error names n, m and the limit);
   - A×A values, computed in **deterministic pair chunks** (fixed chunk size, pairs in
     canonical order), accumulating an `n × k` value matrix, so no `n × r` matrix exists.
     The chunk size is chosen so one chunk stays under the same cell limit;
   - `n × k` value matrices per block and `k × k` covariances, which are small.

   All knowable sizes are checked after the model and cohort are resolved and **before**
   any evaluation or allocation. `.dae_collect_dosages()` is generalised into a shared
   helper (`.collect_dosages(conn, id_ind, locus_ids, caller)`, in
   `R/genome_effects_helpers.R`), and both callers use it. The generator's existing size
   and completeness contract is kept unchanged, and the error text keeps naming the
   caller's own fix. When a trait has no covered loci, no dosages are collected at all.
5. **Monomorphic loci** (Codex review, finding 3). A locus with zero dosage variance in the
   cohort gets regression b_j = 0, and the **centred A×A column** of any pair containing it
   is zero, so the A×A block gets nothing from that pair. The pair's **induced additive
   effect is kept**: canonicalisation adds `e·c` to the partner's α before zero-variance
   columns are handled. With locus 2 fixed at dosage 2, `e(g1 − 1)(g2 − 1) = e(g1 − 1)`, a
   real additive effect at locus 1 (fixed at 0 gives `−e`; fixed heterozygote gives `c = 0`
   and no induced effect). The output is finite, and B12 compares with a reduced model that
   carries that induced main effect.

### 4b.3 Dispatch (§8 cases, per trait)

Read the model with `.gev_read_model()`. Then classify each **family**
(`genome_effect_terms.family_key`: trait, owner and member states).

- **The family rule (new; not in §8).** Scope variants of one family **compete**: the most
  specific matching scope wins. So a family is covered or uncovered **as a whole**. If any
  variant of a family has a line or origin scope, every variant of it, the common one
  included, goes to the uncovered set. Splitting a common variant from its scoped siblings
  and evaluating them separately would count a line-A copy twice. Gate B13 pins this.
- **Covered term shapes** (D2): every member is in the common scope, and every member locus
  resolves to `chr_inheritance` `1,1` for both sexes. The accepted shapes are order-one
  `additive`, order-one `dominance`, order-one diploid `indicator` (any of the three
  states), and two-member `additive × additive`.
- **Case 1, `decomposition = "full"`:** every family is covered.
- **Case 2, `"additive_only"`:** every term is order-one `additive`, and at least one family
  is scoped.
- **Case 3, `"partial"`:** anything else.

**Owners.** Families are classified per owner first, using the evaluator's actual
`family_key` (owner and member states included; centres are deliberately not part of it),
so each owner's variant competition is preserved. Only then are the **covered**
coefficients of all owners summed during canonicalisation. Pooling before classifying could
change which terms apply. Two custom owners with additive and dominance terms at the same
loci still make one case-1 model (B4).

**Partial models.** The uncovered values are aligned to the full cohort by `id_ind`. An
individual that matches a covered family but no uncovered family gets an uncovered value
of 0, not a missing value. The "individual absent" error of 4b.2 rule 2 applies to the
**whole trait model**, never to the uncovered submodel alone.

### 4b.4 Computation

**Case 1 and the covered part of case 3:** `.stored_to_functional()` on the covered terms,
then the source's NOIA projection, ported from `nonadd_covariates()` / `nonadd_decompose()`.
- **Realised:** p comes from the cohort (`colMeans(X)/2`), and b is the observed regression
  (0 at monomorphic loci). Then α = a + b·d + Σ e·c and the value matrices
  BV = Z_A α, DD = Z_D d, AA = Z_AA e, as in the source, with AA accumulated in pair
  chunks (4b.2 rule 4). Observed covariances use the sample divisor `n − 1`; the
  within-locus regression may use population-form moments, since the divisor cancels.
- **Genic:** p comes from `base_tbl`, b = q − p, and the blocks are the closed forms
  `α' diag(2pq) α`, `d' diag((2pq)²) d` and `e' diag(4 p_k q_k p_l q_l) e`.
  `total` is their sum. There is no `between_components` row and no X.

**Case 2:** the `additive` and `total` rows are covariances of the evaluated `additive`
component and of the evaluated total. These are equal up to rounding, since every term is
additive. The roxygen calls this the **evaluated additive variance**, and the D1 message says
so.

**Case 3, uncovered part:** evaluate the uncovered terms alone (`.gev_evaluate(model = sub)`)
to get the value U. `unpartitioned` = Cov(U_t1, U_t2).

**Realised totals and `between_components`.**
- `total` = Cov(g_t1, g_t2) of the **evaluated** total, i.e. what `ind_tgv_total` holds.
- `between_components` = Σ over ordered pairs b ≠ b' of Cov(b_t1, b'_t2). Here b ranges over
  the reported blocks and `unpartitioned`, every orientation included.
- It is computed directly from the value matrices, **not** as total minus the blocks, so B3's
  sum-to-total check is a real check of the algebra. The block value matrices sum to the
  centred evaluated total to rounding, by the §3 identities.
- An internal assertion compares them, per trait, with a mixed absolute/relative tolerance
  on the **centred** values (a large raw genetic offset must not by itself widen the
  allowed error). If they disagree, the call errors and names the trait, rather than
  returning numbers that do not add up. That error would mean the canonicalisation is wrong.

**Block availability** (Codex review, finding 1; replaces "a block row appears only when
the model has terms of that kind"). Stored `component_name` / `.gev_target_kind()` describe
term *structure*; the report is a *statistical* decomposition, so availability is decided
**after canonicalisation**, never from stored term kinds. Counter-example to the old rule: a
one-locus `genotype_terms()` surface `c(0, 1, 2)` is exactly dosage, stores only
indicators, and has its whole variance in the additive block (probe: additive = total =
0.6778). Dominance and A×A terms also induce additive coefficients (`b·d`, `e·c`).

- **`additive`**: always reported for a trait with any covered term (value may be 0).
- **`dominance`**: reported when the canonical model has **dominance support**: some
  locus with canonical `d_j ≠ 0` after all covered owners are summed. This is decided on
  the canonical coefficients, not on stored contrast names, so it is coding-invariant: the
  dosage-only surface above has `d = 0` exactly (the D2 rows only halve) and reports no
  dominance row, just like `ad_terms(a = 1)`. A model without support gets no `dominance`
  row, so a target the generator was told to leave out (§6C) still falls into the
  `anti_join()` (B8's intent kept).
- **`additive_by_additive`**: reported when some canonical pair has `e ≠ 0`.
- A supported block whose value is 0 in this cohort is reported as 0. "Missing row" means
  "the model has no such block", never "measured zero".
- **Per trait pair:** an off-diagonal block row `(b, t1, t2)` appears only when `b` is
  available for **both** traits. One trait's dominance terms never make another trait's
  absent dominance block look measured. Cov(D_t1, ·_t2) for a trait t2 without dominance is
  a cross-block term and goes into `between_components`.
- `unpartitioned` appears for a trait with any uncovered family (case 3); `total` always.

**Off-diagonal rows when the traits are in different cases:** `decomposition` is the less
complete of the two (`partial` > `additive_only` > `full`). `cov_value` is computed the same
way regardless.

**Determinism.** X is ordered by `(id_ind, locus_id)`, terms by `id_genome_effect`, and every
reduction runs in R (`stats::cov`, `crossprod`). The output is a function of the stored rows
alone. B9 asserts `expect_identical()`, including across DuckDB thread counts in the style
of `test-genome-effects-determinism.R`.

### 4b.5 Roxygen content (§8 + step-4 bullets)

- The two anchors, the three cases, `decomposition`, and D1's message.
- `base_tbl`'s two differences from the generators: its default is the cohort, and it is
  genic-only.
- `unpartitioned` is not `residual`.
- **When target vs measured is like-for-like** (Codex review, finding 7). All terms
  `generated` is necessary but not sufficient. The anchor, the reference selection, the
  measured model and the calibrated scope must also match. A generated line variant
  competes with common fallback terms, so a mixed cohort measures the active combination,
  not that variant alone; a restored or later-edited founder pool need not be the
  generation-time base. With custom owners a mismatch is expected. Mismatch is never an
  error.
- **Worked target-vs-measured examples** (step-4 bullet 1). Both filter `trait_var_comp` by
  `line_name` (`is.na(line_name)` for the population-wide block, or a line), drop
  `line_name`, and use `inner_join()` and `anti_join()`:
  - genic: a fully covered, common-scope model calibrated with `anchor = "genic"`,
    measured with `anchor = "genic", base_tbl = <the generation base>`;
  - realised: a model calibrated `"realised"`, measured on its actual calibration
    individuals.
  For a scoped model the example keeps the "evaluated additive variance" wording and
  explains why a mismatch with one line's target can be meaningful.
- **Stored diagonals are generation targets under their calibration contract**, not
  universally "the genic total": the additive generator already supports realised
  calibration, and a scoped target calibrates a variant, not the whole cohort.
- **Why the prevalence rule can fail** (step-4 bullet 3). `define_phenotype(prevalence = )`
  sums the stored diagonals. On a cohort with LD or departures from HWE, or drifted
  frequencies, the realised `total` differs by `between_components` and by the drift in
  each block; this function shows that gap. Even a matching variance does not guarantee the
  prevalence: the threshold also assumes a near-normal liability, which a skewed
  finite-locus genetic model need not give.
- The interpretation of each estimate: genic on a selected cohort is a projection at
  `base_tbl`'s frequencies, realised measures that cohort, and `base_tbl = NULL` drifts with
  the cohort.

### 4b.6 Gates (4b), `tests/testthat/test-extract_genetic_variance.R`

B1–B12 as written in §11, with these specifics:

- **B2:** use the D4 oracle. Write the model with Cockerham `ad_terms()` + `aa_terms()` from
  fixed random `(B_a, B_d, B_aa)` converted with `.noia_to_stored()`. Evaluate it on base
  individuals made by `add_founders()` from an LD panel (`qtl_mkhap`). Compare the oracle's
  `real_*`, `cov_A_D`, `cov_A_AA` and the test's own Cov(D, AA) at 1e-10. The oracle takes
  the cohort's dosage matrix through `extract_genotypes()`.
- **B2** checks each **block** independently against the oracle; **B3** checks the
  cross-block accounting against value matrices the test builds itself. Add hand-derived
  fixtures (a two-locus case worked on paper) beside the oracle, since agreement with the
  same source algebra is not independent proof of its interpretation. The oracle stays
  independent of the production conversion and projection helpers.
- **B3:** use an inbred or repulsion-LD cohort with k = 2 and a deliberately **asymmetric**
  two-trait fixture (Cov(b_t1, b'_t2) ≠ Cov(b'_t1, b_t2)), with a D–A×A covariance and an
  `unpartitioned` part both nonzero. Assert the full ordered sum, and that the blocks alone
  do **not** sum to `total` (difference > 1e-6), so the test cannot pass on an HWE/LE
  panel by accident.
- **B12** (rewritten, Codex finding 3): a locus fixed at 0, 1 and 2 in turn, in a pair.
  The output is finite and equals a **reduced model** that drops the pair and adds its
  induced main effect `e·c` to the partner (and drops the constant). Fixed at 2 and at 0
  the additive block is nonzero for a pair-only model (probe: 0.6778); fixed at 1 it is
  not.
- **B6 and B4:** build the line-scoped F1 through `add_offspring()` from two lines with
  common + line-A + line-B generated variants.
- **B10:** write the same model twice (functional coding, then Cockerham coding) and get the
  same report under both anchors.

Three gates are new in this plan:

| Gate | Checks |
|---|---|
| B13 | Family rule: a model with common + line-A additive variants **and** a common dominance term is case 3. Its `total` equals `var()` of `ind_tgv_total` computed in the test, and the additive families appear in `unpartitioned`. With the rule broken (common variant projected separately) `total` would still match but the blocks would not sum, so the test asserts the B3 identity too |
| B14 | D2: a one-locus `genotype_terms()` surface (all three states) is case 1, and its report equals the same model written with `ad_terms(coding = "functional")` to 1e-10 |
| B15 | D1: every row's `anchor` matches the call. The message names the cohort size (realised) or the `base_tbl` table (genic), and says whether the frequencies are the cohort's or an explicit copy/pool base |
| B16 | Block availability (Codex finding 1), report **contents** asserted, not only totals and coding invariance: dosage-only one-locus indicator surface (additive row = total, no dominance row); heterozygote-only model with unequal genotype frequencies (additive and dominance rows); functional pair-only model (additive and A×A rows); Cockerham pair-only model measured after frequency drift; two traits with different block support (no off-diagonal dominance row; the cross term lands in `between_components`) |
| B17 | Resources (Codex finding 4): a many-pairs model whose `n × r` would exceed the cell limit runs in chunks and matches the oracle; chunked A×A values equal the dense form and every chunk size to 1e-12 (floating-point grouping differs, so not `expect_identical()` across chunk sizes; a given chunk size is bit-identical on repeat, and the chunk size is a function of `n` alone); a model whose dosage matrix exceeds the limit errors before any evaluation, naming n, m and the limit |
| B18 | Default base (Codex finding 5): the same cohort selected through `ind_meta`, through paternal `ind_haplotype` rows, and through repeated phenotype records gives identical genic output with `base_tbl = NULL`; an explicit copy-filtered `base_tbl` gives the appropriately different result; a base with no copies at a required locus errors naming it. Also: `trait_name = NULL` skips a trait with no terms; naming one errors |

### 4b.7 Docs (4b)

NEWS 0.75.1. Add `extract_genetic_variance` to `_pkgdown.yml` with the `extract_*` functions.
- API skill: a new section with the signature, cases, family rule, output, and D1.
- `package_summary.md`: 46 functions.
- CLAUDE.md design principle 4: add `extract_genetic_variance` to the functions that take a
  `tidybreed_table`.
- The schema skill needs nothing, since there is no table change.
- The main plan: an "As built" paragraph under Step 4, and an update to §8 for the family
  rule, D1, D2, D5 and the `additive_by_additive` name, and §11 for B13–B18. Q13 gets the
  `ad_terms()` note from 4a.3.

---

## Out of scope (stays later)

- `extract_breeding_values()`, the generation-t breeding value (§4.3.3, §8). It needs its own
  reference-population spec.
- `extract_functional_effects()` (Q9, low priority).
- `extract_genetic_covariance()`, the block × block breakdown (Q16 (b)).
- `nonadditive_terms()`, the manual Cockerham path with pairs (Q17 (b)).
- Genic decomposition of scoped models (§8: it needs per-line and per-origin frequencies).

## Risks and things to watch

- **The writer and `NA` columns (4a.1).** If `define_genome_effect_terms()` rejects a fixed
  column set with `NA`s in irrelevant columns, the fix goes in the writer. Do not have the
  builders drop columns, because C16 depends on one shape.
- **D5 changes `genotype_terms()` `term_id`s.** Existing tests that compare written terms by
  `term_id` will change. `term_id` is never stored, so only builder-output assertions move.
  Census them at the start of 4a.
- **Realised totals on large cohorts.** The evaluator runs on every selected individual. That
  is the same cost as `add_tgv()`, plus one R covariance. The X guard bounds only the covered
  loci. Add one benchmark in `dev/benchmarks/benchmark_extract_genetic_variance.R`
  (2,000 individuals, 500 loci, A + D + A×A, **2 traits, 1,000 pairs**, stated in the
  script) and a second, many-pairs case (all 124,750 pairs of the 500 loci, where a literal
  port would need ~2 GB for one pair matrix). Record both times and peak memory in the
  phase summary.
- **Exactness of the 1e-10 gates** depends on the scale of the coefficients. Keep fixtures
  O(1). If a gate is tolerance-sensitive, assert a relative tolerance and say so in the test.
- **The sum-to-total assertion (4b.4) could fire on legitimate models** if the evaluator's
  DECIMAL accumulation differs from R's double sums by more than the tolerance. Measure on B3
  first. If needed, scale the tolerance to `max(abs(g))`.

## Verification (each sub-step)

1. `pkgload::load_all()`, then `testthat::test_file()` on the touched files.
2. Full suite with `NOT_CRAN=true`, in the background: 0 failures. The warning count is
   compared with 0.74.5 (6), and any new warning is explained.
3. `devtools::document()` and `pkgdown::check_pkgdown()` run clean.
4. Mutation spot-checks:
   - 4a: flip the sign of the dominance row's `a −= v(1 − 2c)` → N1 must fail.
   - 4b: project a common variant separately from its scoped siblings → B13 must fail.
   - 4b: omit the D–A×A covariance, omit the `unpartitioned` cross terms, or use only one
     off-diagonal orientation → B3 (asymmetric two-trait fixture) must fail. ("total minus
     blocks" is **not** a mutation: it is algebraically equal to the ordered sum when the
     identities hold, so an honest test may pass on it. The direct computation stays as the
     implementation requirement.)
   - 4b: drop a monomorphic pair without carrying its `e·c` → B12 must fail.
   - 4b: decide block availability from stored contrast names → B16 must fail.
5. After each sub-step: write the main plan's "As built" paragraph, then commit and push,
   then **wait for your review**. Results go in `import_qtl_effect_methods_phase_4.md`
   (one file for all of step 4).

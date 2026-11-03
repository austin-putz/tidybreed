# Plan: check policy — refuse by default, `allow =` to bypass

Status: draft for review (user + Codex). Nothing implemented.
Decisions recorded 2026-10-13 (see "Check policy" in `plans/TODO.md`); this plan
turns them into a mechanism, an inventory and a phased implementation.

## 1. Goal

Users may model anything: unusual genomes, chromosome rules that differ between
lines, any custom genome effect. The package's job is to stop or warn loudly by
default when an input breaks a modelling assumption, and to let a user who
understands the issue say so explicitly and carry on. From then on the package
computes whatever is well-defined and accepts what the user did, without
second-guessing them in later calls.

## 2. Policy

### 2.1 Four kinds of check

| Kind | What it protects | Bypassable? |
|---|---|---|
| **0. Input shape** | Malformed arguments: wrong type or length, `NA`, unknown names, a missing table or column | No. Nothing to bypass to |
| **1. Integrity and safety** | Database consistency (orphans, ids, the `genome_effect*` deletion rule), bit-identical determinism, failure/rollback contracts, the closed formula grammar (stored formulas are `eval()`ed), SQL built only from registered ids, the 0.76.10 `define_chromosome()` guard | Never |
| **2. Not implemented** | Cases with no code path (ploidy > 2, `2,0` transmission) | Only once a behavior is defined; then it becomes kind 3 |
| **3. Modelling assumption** | The math is well-defined either way, but the user may get something they did not mean | Yes, by name |

### 2.2 When a kind-3 check may be bypassed

**A check may be bypassed only when every later function has a defined,
documented behavior for the state the check let through.** The bypass and that
behavior ship together. Otherwise the bypass moves the crash to a later call.

Example: a dominance term on a chromosome that line H carries once.

- The writer's `allow = "non_diploid_dominance"` stores the term.
- `add_tgv()` meets H animals with one copy. By default it stops. With its own
  `allow = "non_diploid_dominance"` it skips that term's contribution for
  those animals and gives a message with the counts.

### 2.3 Default level

- **Refuse** when the result would be wrong without anything saying so:
  calibration, a missing genotype state, a variance target that would not be
  met.
- **Warn** when the result is literally what the user asked for, only unusual.
  Example: additive effects on a hemizygous chromosome, which are counted per
  copy, correctly.

### 2.4 Later code re-checks, nothing is stored

A relaxation is not recorded in the database. A function that relies on an
assumption checks it again itself, and needs its own `allow =` to proceed.

Example: `define_phenotype(prevalence = )` relies on every `"generated"` term
being calibrated (Q21). If a generator ran with the diploid check relaxed,
`define_phenotype()` re-checks the copy counts of those terms' loci and refuses
unless it is given the same `allow` value.

This keeps the database free of policy flags, in line with the project's
preference for computing over storing metadata. It also keeps each function
responsible for its own assumptions, and keeps scripts self-describing: every
risk a script takes is visible in its calls.

## 3. Mechanism

### 3.1 The argument

- **Name and type:** `allow = character()`, a character vector of relaxation
  names. It has the same name in every function that has relaxable checks.
- **Position:** after the function's model arguments and before `seed` where
  there is one, so `seed` stays last.
- **Unknown names are refused** (kind 0), listing the names that function
  accepts. A misspelt `allow` must never be ignored.
- **Unused names:** a valid name that the call never needed is silently
  accepted. A script can then pass the same `allow` to every call in a
  pipeline.

### 3.2 One registry

`R/check_policy.R` (new) holds `ALLOW_NAMES`, a data frame with:

- `allow_name`;
- the default level (`refuse` or `warn`);
- a one-line description;
- the functions that accept it.

From this one table:

- **`.check_allow(allow, fn)`** validates the argument against `fn`'s names.
- **One documentation page, `?allow`,** lists every name, what it relaxes and
  what later functions then do. It is a topic page, not a function. Each
  function's `@param allow` lists its own names and links there.
- **A consistency test** (like `test-schema-print.R`) checks three things:
  every registry name is used by code, every name the code uses is in the
  registry, and every function's documented names match the registry.

### 3.3 Raising and relaxing

Two internal helpers, used by every kind-3 check:

```r
.refuse_unless(cond_ok, allow, "non_diploid_dominance", msg)
#   ok          -> nothing
#   not ok, not allowed -> stop() with class "tidybreed_check_error";
#                  `allow_name` is a field. The message ends with
#                  'To accept this, pass allow = "non_diploid_dominance".'
#   not ok, allowed     -> one message, "Allowed non_diploid_dominance: <counts>",
#                  respecting the tidybreed.quiet option

.warn_unless(cond_ok, allow, "hemizygous_additive", msg)
#   not ok, not allowed -> warning() with class "tidybreed_check_warning",
#                  ending the same way ("To silence this, pass ...")
#   not ok, allowed     -> silent (acknowledged)
```

Classed conditions let users and tests catch a specific check without matching
message text.

### 3.4 Invariants every relaxation keeps

- **No RNG and no writes before checks.** All checks, relaxed or not, still run
  before the first RNG draw and the first write. A relaxed call consumes the
  RNG exactly as the same call would if the check had passed naturally. So
  relaxing never changes seeded output of the valid part of a call.
- **Atomic.** A refused call writes nothing (existing contracts: D7,
  `expect_refusal()` in the genome-effect tests).
- **Kind 1 never gets an `allow` name.** The consistency test asserts that no
  integrity check routes through `.refuse_unless()`.

### 3.5 Existing opt-out arguments

The package already has arguments that act like opt-outs. Each is either kept
because it means something different, folded into `allow`, or renamed. See the
inventory, §4.4.

## 4. Inventory

Compiled from every `stop()` and `warning()` in `R/` at 0.76.12 by four
read-only agents, one per area. The raw tables, with `file:line` for each
check, are in [check_policy_inventory.md](check_policy_inventory.md).

### 4.1 Counts

| Area | Kind 0 | Kind 1 | Kind 2 | Kind 3 | Warnings |
|---|---|---|---|---|---|
| A. Genome effects, generators, `add_tgv()` | ~165 | ~30 (16 of them Q21) | 12 | 7 | 11 |
| B. Phenotype, covariance, index, EBV | 134 | 23 | 26 | 22 (5 already have an argument) | 19 |
| C. Genome, chromosomes, breeding | 80 | 28 | 11 | 3 (+5 proposed new checks) | 2 |
| D. Infrastructure, tables, groups | 62 | 10 | 13 | 15 (8 refusals, 7 warnings) | 7 more, informational |

About 630 refusals in all. Roughly 70% are kind 0 and stay as they are. The
policy work is the ~50 kind-3 checks, plus fixing the warnings that are at the
wrong level.

**Q21 refusals stay kind 1.** These keep "generated means calibrated" true.
Most of them are not assumptions at all; they are targets with no
mathematical solution:

- the anchor or architecture rank is below the target rank;
- the additive floor is missed;
- the inbreeding-depression ratio cannot be reached;
- calibration misses its target beyond tolerance.

No `allow` name can produce a number that does not exist.

### 4.2 Proposed `allow` names

One name per situation. A name can be shared across functions where the
situation is the same. "New" marks a check that does not exist today.

**Genome, chromosomes and genome effects**

| `allow` name | Default | Where (writer → re-checkers) | Behavior when allowed |
|---|---|---|---|
| `non_diploid_dominance` | refuse | `define_genome_effect_terms()` (`.ge_validate_dominance_ploidy`), `define_chromosome()` re-check → `add_tgv()`, `add_phenotype()`, `extract_genetic_variance()` | The dominance contribution is 0 for an animal carrying ≠ 2 copies at that locus. Message with counts of animals and terms. |
| `incomplete_genotype_surface` | **refuse** (today it passes silently unless `require_complete = TRUE`) | `define_genome_effect_terms()` | Unlisted states score 0, as the evaluator already does. `require_complete` is removed. |
| `unreachable_genotype_state` (new) | warn | `define_genome_effect_terms()` | The term is stored and never fires. Today this case is silent. |
| `no_matching_copy` | refuse | `add_tgv()`, `add_phenotype()`, `extract_genetic_variance()` | The animal's genetic value is 0, written as explicit 0 rows so it reads as 0, not as a missing component. |
| `trait_without_effects` | refuse | `add_tgv()`, `add_phenotype()` (also covers stage 1's "phenotype without genetic effects") | Same 0 rows; a purely environmental phenotype becomes possible. |
| `ambiguous_chr_rule` | refuse | `define_chromosome()` (shadowing check) | The resolver's precedence decides: sex-specific beats line-specific. |
| `long_genetic_map` | refuse | `define_genome()` (new, earlier) → `add_offspring()` (R sampler **and** C++ kernel) | Chromosomes over 3000 cM are drawn exactly, just slowly. The name must reach both kernels, or R ↔ C++ parity breaks. |
| `excess_templates` | **warn** (refuses today) | `define_founder_haplotypes()` | Already well-defined; unused templates cost only time. |
| `parent_sex_role` (new) | refuse when a sex-chromosome rule depends on parent roles; otherwise no check | `add_offspring()` | Transmit by parent slot, as now. Today a female in `id_parent_1` gives a cryptic error at `add_offspring.R:740`. |
| `founder_pool_other_line` (new) | warn | `add_founders()` | Sample what the filter returned, as now. The docs claim it samples only its own line's pool; that gets fixed. |

**Phenotype and evaluation**

| `allow` name | Default | Where | Behavior when allowed |
|---|---|---|---|
| `prevalence_from_stored_target` | refuse | `add_phenotype()` (5 checks: user terms, non-standard terms, two parent scopes, conditional residual, gamma effect) | The threshold is placed from the stored variances even though they do not describe this liability. The realized prevalence will differ. |
| `redefine_realized_block` | refuse | `define_phenotype()`, `define_residual_cov()`, `define_effect_cov_matrix()`, `define_effect_random()` | New records condition on the stored draws under the new matrix, which the draw code already supports. **Also closes a hole:** `define_effect_random(overwrite = TRUE)` deletes one member's draws and so bypasses the lock today. |
| `non_finite_formula_value` | refuse | `formula_tgv` (refuses today) **and** derived `formula` (warns and writes NA today) | The value is treated as missing and routed through `missing_component_action`. One rule for both formula kinds. |
| `unlisted_class_level` | refuse | `add_phenotype()` fixed-class effects | Shift 0 (the reference level); the next line of code already does this. |
| `index_missing_values` | refuse | `add_index()` | `index_value = NA` for those animals. |
| `structural_additive_index` (warning today) | refuse | `add_tgv(index_names =)` | A true index on `additive` for traits with hand-written indicator or interaction terms. Computed as asked. |

**Data manipulation and infrastructure**

| `allow` name | Default | Where | Behavior when allowed |
|---|---|---|---|
| `empty_selection` | warn | `remove_rows()`, `mutate_table()`, `mutate_group_*()` | No-op; the name silences the warning. |
| `unknown_ids` | refuse | `mutate_table()` with a tibble | Unmatched ids are skipped. Today a mistyped id is dropped silently. |
| `untyped_na_column` | **refuse** (warns today) | `mutate_table()`, `define_table()`, extra columns on `add_*()` | The column is created as VARCHAR, as now. Today later numeric writes are stored as text with no further warning. |
| `unassigned_rows` | warn | `mutate_group_seq()`, `mutate_group_named()` | Rows left NULL. |
| `empty_groups` | refuse | `mutate_group_seq()` | Fewer groups than asked. |
| `delete_all_rows` | refuse | `remove_rows()` (replaces `confirm_all`) | Empties the table, but only where no other table depends on it (see §4.3). |
| `replace_file` | refuse | `close_pop(results_dir =)`, and `open_pop()` if Q3 says so | The destination is replaced. |

**Not proposed:**

- `ignore_stored_target`: filtering `trait_var_comp_tbl` already does this,
  and it is explicit.
- `replace_stored_target`: "a stored target is never overwritten" is a design
  rule, not an assumption.
- `qtl_outside_base`: dropping QTL silently changes their count. It stays a
  refusal.

**Waiting on code (kind 2 until then):**

- **`non_psd_covariance`.** The draw code re-checks positive
  semi-definiteness, so a bypass first needs a deterministic projection to the
  nearest valid matrix, applied identically in `add_phenotype()` and in the
  BLUPF90 R. The writer and reader also use different tolerances today, which
  should be aligned whatever we decide.
- **`non_autosomal_qtl`.** See Q1.

### 4.3 Integrity gaps and bugs found (no `allow`; fix as kind 1)

These are not policy questions. The inventory found them, and each is either
a bug or a missing kind-1 check. "Verified" means checked by hand.

| # | Problem | Where | Proposed fix |
|---|---|---|---|
| G1 | **Line names containing `_` break id numbering (verified).** A second `add_founders(line_name = "A_B")` and any `add_offspring()` into `A_B` fail with a DuckDB conversion error. Line `A` also collides with `A_B_*` ids. | `add_founders.R:197`, `add_offspring.R:463` | See Q4. |
| G2 | **Making a BLUPF90 run directory consumes the R RNG (verified).** `sample()` names the folder, so `add_ebv()` shifts every later seeded draw. | `open_pop.R:233` | Name it from the time stamp plus a counter, with no RNG. |
| G3 | **`remove_rows()` can delete configuration and parent rows (verified table list):** `chr_inheritance`, `chr_recombination`, `genome_map`, `genome_meta`, `founder_haplotypes`, `trait_meta`, `phenotype_meta`, `ind_meta`. No validator runs, so it bypasses the 0.76.10 guard and leaves orphans and unresolvable rules. | `sql_utils.R:174` (`TABLE_ROW_KEYS`) | Refuse configuration tables (they have writers with `overwrite`/`mode`). For parent tables, refuse when dependants exist, or cascade through the existing multi-table mode. |
| G4 | **`define_chromosome()` never re-runs the genome-effect checks**, so "effects first, then rules" escapes them, even for line-agnostic X/Y rules. | `define_chromosome.R:298` | Re-validate before `COMMIT` (part of phase 2). |
| G5 | **The three ploidy checks read only line-agnostic rules.** | `assert_qtl_autosomal`, `.ge_validate_dominance_ploidy`, `.ge_reachable_copy_counts` | Resolve over every (sex, line) context: diploid = two copies everywhere, reachable counts = union (phase 2). |
| G6 | **`add_offspring()` pastes parent ids into SQL text**, against the hard rule; an id containing `'` breaks the query. | `add_offspring.R:234`, `:365` | Register them as a view, like `resolve_subset_ids()`. |
| G7 | **`add_founders()` with a partial pool fills unselected loci with allele 0** without saying so. | `add_founders.R:188` | Refuse unless the pool covers every locus. |
| G8 | **A derived `formula` symbol with no records gives silent NA records.** The docs promise a hard error, which happens only when *no* symbol has records. | `formula_helpers.R:537` | Refuse any unknown symbol at `add_phenotype()` time. |
| G9 | **`define_chip()` adds the column and fills it outside a transaction.** A failed fill leaves the column, and every retry is then refused. | `define_chip.R:98`, `:103` | One transaction. |
| G10 | **`archive_replicate()` writes several tables outside a transaction**, and nothing stops the same replicate being archived twice (duplicate rows). | `archive_replicate.R:188` | Transaction, plus refuse an existing replicate number. |
| G11 | **`slice_sample()` on a `tidybreed_table` uses DuckDB's own `random()`,** so `set.seed()` does not control it. Filtered writers also ignore the slice. (Reported, not verified.) | `mutate_table.R:251` | Refuse `slice_*` on `tidybreed_table`, or record it reproducibly. |
| G12 | **Legacy paths that CLAUDE.md says to delete:** the "old format" `founder_haplotypes` branch, the pre-v0.36 `_schema_meta` warning in `restore_pop()`, and the unreachable `assert_ploidy_2` re-check. | `define_founder_haplotypes.R:281`, `restore_pop.R:264`, `ploidy_helpers.R:31` | Delete. Keep one named "ploidy > 2 not implemented" refusal at the writers. `add_offspring.R:567` is a bare `stopifnot` today. |

The infrastructure agent also flagged the `seed = NULL → set.seed(seed)`
pattern as touching `.Random.seed`. That is the established, accepted pattern
(an argument the user passes, not seed restoration), so it is not a gap.

### 4.4 Existing arguments that act like opt-outs

| Argument | Today | Proposal |
|---|---|---|
| `missing_component_action`, `null_class_action`, `condition_change_action`, `underfill_action`, `null_action`, `include_leftover`, `repeatable` | Choose between defined behaviors, mostly stored per phenotype or effect | **Keep.** They are model configuration, not acknowledgement of a risk. |
| `overwrite` (`define_trait`, `define_phenotype`, `define_effect_*`, `define_table`, `mutate_group_*`, `define_index`, `define_chromosome`) | Replace an existing definition | **Keep as the verb for replacing a definition** (Q2), but make it consistent:<br>• `define_index(overwrite = FALSE)` is a silent no-op today; it should refuse.<br>• `mutate_table()` replaces values with only a warning, while `mutate_group_*()` refuses.<br>• `define_phenotype(overwrite = TRUE)` has no lock once records exist (cf. `redefine_realized_block`). |
| `confirm_all` (`remove_rows`) | Opt-out of the whole-table refusal | **Fold into** `allow = "delete_all_rows"`. |
| `require_complete` (`define_genome_effect_terms`) | Opt-**in** strictness | **Remove.** Refuse by default; `allow = "incomplete_genotype_surface"`. |
| `clean = TRUE` (`open_pop`) | Silently deletes an existing database, by default | Q3. |
| `warn_bounds = NULL` (generators) | Silences a calibration diagnostic | **Keep** as a tuning parameter. It sets bounds; it is not a check. |
| `options(tidybreed.label_vector_max / _warn)` | Resource guard on the evaluator | **Keep as options** (Q5). Resource limits do not change results, only whether the computation is attempted. |
| `allow_reserved_owner` (internal `.ge_write_terms()`) | Lets the generators write owner `"generated"` | **Rename** (e.g. `.generator_owner`) so it cannot be confused with `allow =`. Internal only. |

### 4.5 Warnings at the wrong level

- **Should refuse (an argument is silently ignored, or a failure goes
  unreported):**
  - `define_additive_effects()` `distribution` ignored for multi-trait calls (`:360`);
  - `add_ebv()` `phenotype` filter ignored (`:184`) or matching nothing (`:191`);
  - BLUPF90 solutions matching no ids (`blupf90_helpers.R:567`);
  - unknown table names in `summary()` (`summary_pop.R:48`);
  - `pop$tables` out of sync with the database (`tidybreed_pop.R:68`).
- **Should be a message (a normal state):**
  - `restore_pop()` on an empty or undefined genome (`:95`, `:100`);
  - `add_ebv()` using the latest of several evaluations (`:330`).
- **Wrong shape:** `add_ebv()` warns once per animal × trait for a missing
  parent EBV (`:369`). It should be one warning with a count.
- **Should refuse before closing:** `close_pop()` warns about an in-memory
  `results_dir` only after it has already disconnected (`:296`).
- **Noise:** `add_index()` "no filter applied" (`:192`) duplicates the refusal
  at `:215`.

## 5. Phases

Each phase is one or more commits, with its own short plan section or file
when you want one. The full suite and `R CMD check` run before each commit,
and you review between phases.

**Phase 1: integrity fixes (no `allow`).** G1–G3 and G6–G12 from §4.3, plus
the warning-level fixes in §4.5. All are independent of the mechanism. Most
are small and could be split into 2–3 commits by area.

**Phase 2: mechanism, and its first use on chromosomes and genome effects.**
- `R/check_policy.R`: the `ALLOW_NAMES` registry, `.check_allow()`,
  `.refuse_unless()`, `.warn_unless()`. Generic argument helpers stay in
  `R/arg_checks.R`.
- The `?allow` topic page and the registry consistency test.
- G4 and G5: copy counts resolved over every context, and
  `define_chromosome()` re-validation.
- Names: `non_diploid_dominance` (writer, `add_tgv()`, `add_phenotype()`
  passthrough, `extract_genetic_variance()`), `incomplete_genotype_surface`
  (removing `require_complete`), `unreachable_genotype_state`,
  `no_matching_copy`, `trait_without_effects`, `ambiguous_chr_rule`.
- Rename `allow_reserved_owner`.
- Tests for each name:
  - the default refuses with the classed condition, leaving the database and
    RNG untouched;
  - `allow` proceeds;
  - the defined downstream behavior holds;
  - a re-checker without its own `allow` refuses;
  - both configuration orders (rules then effects, effects then rules).

**Phase 3: breeding and founders.** `long_genetic_map` (R and C++, with a
parity test), `excess_templates`, `parent_sex_role`,
`founder_pool_other_line`.

**Phase 4: phenotype path.** `prevalence_from_stored_target`,
`redefine_realized_block` (and closing the `overwrite` hole),
`non_finite_formula_value`, `unlisted_class_level`, `index_missing_values`,
`structural_additive_index`.

**Phase 5: data manipulation.** `empty_selection`, `unknown_ids`,
`untyped_na_column`, `unassigned_rows`, `empty_groups`, `delete_all_rows`,
`replace_file`, and the `overwrite` consistency fixes (§4.4).

**Later, each needing its own design:** `non_psd_covariance` (a
deterministic projection), and `non_autosomal_qtl` if Q1 chooses real
hemizygous calibration.

Phases 3–5 do not depend on each other and can go in any order after phase 2.

## 6. Questions for review

**Q1. Generated QTL on chromosomes that are not diploid everywhere.**
`assert_qtl_autosomal()` is kind 2 today, because the generators can only
calibrate with 2pq. You said users should be able to "ignore the chromosomes
and just use them as is". Options:
- **(a)** `allow = "non_autosomal_qtl"`: the generator calibrates as if every
  animal were diploid, and says so. Calibration-dependent code re-checks the
  loci's copy counts and refuses without its own `allow`. The main case is
  `define_phenotype(prevalence = )`, and BLUPF90's G prior should get the same
  re-check. This turns Q21's "generated means calibrated" from trust into a
  re-check, which is what §2.4 asks for anyway.
- **(b)** implement real hemizygous calibration (pq a² per copy, weighted by
  sex), then the check disappears. Larger, and the variance then differs by
  sex, which needs its own decision.
- **Recommendation: (a) now, (b) only if someone needs calibrated sex-linked
  variance.**

**Q2. `overwrite` vs `allow`.**
Keep `overwrite` as the verb for "replace this definition", and use `allow` only
for modelling assumptions? Replacing a definition is something the user asks
for, not a risk they accept. The exception is when replacing changes realized
data, which gets an `allow` name (`redefine_realized_block`). The
inconsistencies in §4.4 get fixed either way. **Recommendation: keep
`overwrite`.**

**Q3. `open_pop(clean = TRUE)` deletes an existing database by default.**
It is the most destructive default in the package, but re-running a script
from the top depends on it. Options:
- (a) keep it as it is;
- (b) default `clean = FALSE`, refusing an existing file;
- (c) keep deleting, but print a message naming the file.

**Recommendation: (c).** Re-running a script is the normal workflow and
`restore_pop()` exists for reopening, so a refusal would mostly get in the way.

**Q4. Line names with `_`.**
Options:
- (a) forbid `_` in line names, a kind-0 refusal in every writer that takes
  `line_name`;
- (b) parse the number after the **last** `_`, with an exact prefix match.

**Recommendation: (b).** Ids stay `{line}_{n}`, and names such as `Duroc_US`
work.

**Q5. Resource guards stay as `options()`.**
These are the evaluator's label-vector limits and the dosage-matrix cell
limits. They don't change results, only whether a computation is attempted,
so they sit outside the "visible in the call" rule. **Recommendation: yes.**

**Q6. `extract_genotypes()` / `add_dosage()` on absent copies.**
A female's Y, or a missing MT copy, reads as dosage 0 today. That is the same
value as homozygous reference. Options:
- (a) return `NA` for absent copies by default;
- (b) keep 0 with `allow = "non_diploid_dosage"`, warning by default.

**Recommendation: (a).** `NA` is the true answer, and downstream tools treat
it as missing.

**Q7. Name of the documentation page.** `?allow` is short, but generic.
`?tidybreed_allow` is unambiguous. **Recommendation: `?tidybreed_allow`,**
with each function's `@param allow` linking to it.

**Q8. Warn-level names silence the warning.** For a warn-level check, passing
its name means "I know, stop warning me". Agreed?

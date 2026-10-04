# Import QTL-effect methods — Step 3 plan

**Spec:**
- `plans/import_qtl_effect_methods.md`: §6 (P1), §6A (P2), §6B (value half), Q18, Q21 and
  §10 Step 3.
- Gates T1–T9 and PH1–PH8 (§11).

**Versions:** 0.74.0 (3a), 0.74.1 (3b), 0.74.2 (3c).
**Status:** 3a built 2026-10-04 (0.74.0; as-built notes in the main plan's Step 3). 3b and 3c planned.
**Starting point:** 0.73.2 (`be43687`).

Step 3 makes `ind_tgv` the only table of true genetic values. Every phenotype sees the
total genetic value. The prevalence threshold can trust a stored target, because
`generated` comes to mean "calibrated". After this step, dominance and A×A terms written
by Part C (step 5) reach phenotypes.

## Decisions taken while planning (2026-10-04)

1. **`add_tgv(component_name = )`**, not `component =` (naming rule 1: the argument
   filters `ind_tgv.component_name`). The default is `"additive"`, because selection
   indices are on breeding values. `"total"` is also allowed. The formula DSL keeps the
   named-only `component =` that Q18 decided, because it is DSL syntax, not a column
   argument.
2. **`ind_true_index` gains `component_name VARCHAR NOT NULL`.** Its row key becomes
   `(id_ind, index_name, weight_type, component_name)`, in `TABLE_ROW_KEYS`, the DDL and
   `.sm_col()`. Without it an additive index and a total index share a key, and one
   silently replaces the other.
3. **Three commits.** Each is pushed with the full suite green, and each waits for the
   user's review:

| Sub-step | Version | Content |
|---|---|---|
| 3a | 0.74.0 | P1 consolidation + P2 readers + the value-name half of §6B + the active-block prevalence rule |
| 3b | 0.74.1 | Q21: `effects` / `scale_to_target` removed, the owner rule, the `define_effect_cov_matrix()` refusal, the reworded "already stored" error |
| 3c | 0.74.2 | Q18: `formula_tbv` → `formula_tgv`, the DSL's named-only `component =` / `table =`, identifier validation |

The order matters:
- 3a can keep the old `effects` path, because those terms are still `generated`, and
  `add_tgv()` reads every owner.
- 3b's owner rule needs `effects` gone first. Otherwise `generated` proves nothing.
- 3c is a rename plus parser work, independent of the others.

## Census (0.73.2; recount at the start of each sub-step)

| Pattern | R/ | tests/ | Elsewhere |
|---|---|---|---|
| `add_tbv(` calls in tests | — | **108** in 13 files (49 `test-add_tbv_index.R`, 23 `test-add_tbv.R`, 16 `test-genome-effects-eval.R`) | — |
| `ind_tbv` | 11 files | 25 files | 2 vignettes, README, CLAUDE.md, both skills, `package_summary.md`, `dev/` |
| `tbv_value` / `id_tbv` | 8 / 4 files | 20 files | swine vignette (13) |
| `order1_*` | 5 files | 2 files | `man/add_tgv`, `man/define_phenotype`, schema skill, swine vignette |
| `effects =` / `scale_to_target =` on `define_additive_effects()` | `add_tbv.R` example | **93** distinct calls in 8 files | swine vignette 6×, README 2×, `dev/benchmarks` 4× |
| `formula_tbv` | 9 files | 4 files (36 in `test-formula_phenotype.R`) | swine vignette, README, CLAUDE.md, skills |

`tools/quarto/legacy/` and `docs/` are frozen or rendered output. Both are excluded from
the gates, as step 2 excluded them: legacy still has `target_add_var`.

---

## 3a — Consolidation and P2 (0.74.0)

### 3a.1 Value-component rename (§6B, value half)
- `.gev_component()` (`R/genome_effects_eval.R:78-84`): `order1_additive` →
  `additive`, `order1_dominance` → `dominance`, `order1_other` → `indicator`. The
  `interaction` name is unchanged.
- The `component_name` description (`R/schema.R:544-545`), the `add_tgv.R` roxygen and
  comments, `test-genome-effects-eval.R` and `test-genome-effects-schema.R`.
- Add a closed set `TGV_COMPONENT_NAMES <- c("additive", "dominance", "indicator",
  "interaction")`, next to `GENETIC_EFFECT_NAMES`. It is an R-validated constant, not an
  SQL `CHECK` (§6.2). The evaluator asserts that it only produces these names.

### 3a.2 `add_tgv()` inherits from `add_tbv()`
New signature:
`add_tgv(tbl, trait_name = NULL, index_names = NULL, weight_type = c("index",
"economic", "both"), component_name = "additive", overwrite_index = FALSE, ...)`.
- `...` is scalar custom-field forwarding into `ind_tgv` through `prepare_extra_cols()`.
  Validation is copied from `R/add_tbv.R:147-156`.
- The true-index block is ported from `R/add_tbv.R:194-325`:
  - it reads `ind_tgv` rows for `component_name`, or `ind_tgv_total` for `"total"`;
  - a component with no row for an individual × trait contributes 0, but only if the
    individual has some `ind_tgv` row for that trait. An individual with no row at all
    is the existing "no TBV rows" error, reworded;
  - it writes `component_name` on every row;
  - `overwrite_index` deletes by `(id_ind, index_name, weight_type, component_name)`.
- `component_name` is validated against `c(TGV_COMPONENT_NAMES, "total")`.
- `ind_true_index`'s value is still `tgv_mat %*% wt_vec` in R. This is a fixed-order
  product over at most a few traits, not a SQL `SUM()`, so it stays deterministic.
- New roxygen sentence (§6.2): "the breeding value is `component_name = 'additive'` for
  generated effects (statistical coding at one base `p`). For hand-written functional
  terms it is the functional additive effect." The same sentence goes on
  `define_genome_effect_terms()` and `ad_terms()`.

### 3a.3 Delete `ind_tbv` and `add_tbv()`
- **Files:** delete `R/add_tbv.R` and `man/add_tbv.Rd` / `upsert_ind_tbv.Rd`. Remove the
  `NAMESPACE` export and the `_pkgdown.yml:81` entry.
- **DDL:** remove `CREATE TABLE ind_tbv` from `R/define_trait.R:206-214`, and rewrite
  the comment at :216-222 that compares `ind_tgv` with it.
- **`ind_true_index` DDL** (`R/define_trait.R:269-276`): add
  `component_name VARCHAR NOT NULL`.
- **Registries** (`R/sql_utils.R`):
  - drop `ind_tbv` from `TABLE_RESERVED_COLS:113`, `TABLE_PRIMARY_KEYS:155`,
    `TABLE_ROW_KEYS:182`, `SYSTEM_TABLES:299` and `IND_TABLE_ID_IND_COLS:310`;
  - add `component_name` to `ind_true_index` in `TABLE_RESERVED_COLS:141` and
    `TABLE_ROW_KEYS:186`.
- **`R/schema.R`:**
  - remove the `ind_tbv` block (:524-534) and its entry in `.schema_table_order()`
    (:780-782);
  - fix the `trait_meta.trait_name` text (:322);
  - fix the `ind_true_index` texts (:612, :622, now "from `ind_tgv`"), and add a
    `component_name` `.sm_col()`;
  - fix the `.schema_table_order()` roxygen (:755).
- **Evaluator helpers** (`R/genome_effects_eval.R`): delete `.gev_reserved_additive()`
  (:170-191), the `tbv` branch of `.gev_require_terms()` (:808-816) and
  `.gev_warn_tbv_stale()` (:822-890). Fix the header comment (:5) and the roxygen at
  :627.
- **Other readers:**
  - `R/add_index.R:118-122`: the map drops `ind_tbv` and gains `ind_tgv = "tgv_value"`.
    Its existing "more than one value per individual × trait" error is the T6 safety
    net. Update the roxygen at :13-71 with the `filter(component_name == "additive")`
    example;
  - `R/archive_replicate.R:110-119`: drop `ind_tbv` from `store_and_reset`;
  - `R/tidybreed_pop.R:202`: the print summary counts `ind_tgv` as "TGV";
  - roxygen mentions in `remove_rows.R:132,171`, `define_index.R:45`,
    `add_dosage.R:12`, `sql_utils.R:362,473`, `define_additive_effects.R:8,576` and
    `define_trait.R:5`.
- **`restore_pop()`** (pattern `stop_stale()`, `R/restore_pop.R:106-118`). Refuse a
  file that still has:
  - an `ind_tbv` table;
  - an `ind_true_index` with no `component_name`;
  - any stored `order1_*` value in `ind_tgv.component_name` or
    `phenotype_components.component_names`.
  The comments at :109 and :150 also stop naming `add_tbv()`.

### 3a.4 P2 — every phenotype path reads `ind_tgv` (§6A)
- **Stage 1** (`R/add_phenotype_stages.R:285-367`, `.ap_materialize_tbvs()`): the three
  `add_tbv()` calls (:298, :332, :360) become `add_tgv()`. The internal name stays until
  3c renames it.
- **One reader replaces `.tbv_by_id()`.** It is
  `.tgv_by_id(conn, trait_name, ids, components)` in `R/contributor_tbv.R`.
  - `components = "total"` reads the `ind_tgv_total` view.
  - A listed set reads `CAST(SUM(CAST(CASE WHEN component_name IN (…) THEN tgv_value
    ELSE 0 END AS GEV_ACC_TYPE)) AS DOUBLE)`, grouped by `id_ind` over the trait's
    rows.
  - An individual with no row is `NA`, and `missing_component_action` then applies. A
    present individual with no row of a listed component gets 0.
  - It registers ids as a view, never writing them into SQL (CLAUDE.md).
- **Simple phenotypes** (`:482-487`) use `.tgv_by_id(..., "total")`.
- **Composite** (`R/add_phenotype.R:326-370`, `.assemble_composite_tbv()`) reads
  `comp$component_names` (split on `,`) for self, dam, sire and group.
- **Group mates** (`.group_mate_tbv()`, `R/contributor_tbv.R:121-148`) take
  `components`. The inner join source is `ind_tgv_total`, or the filtered exact sum over
  `ind_tgv`. It keeps `GEV_ACC_TYPE`, the integer `COUNT`, the
  `.gev_accumulator_error` handler, and the rule that a focal with no mates gets 0.
- **The formula DSL** (`R/formula_helpers.R:305-321`, `.build_tbv_env()`) reads
  `"total"` for every reference. 3c adds `component =`.
- **`ind_tgv_total` sums deterministically** (`R/genome_effects_helpers.R:87-93`).
  *As built:* `list_sum(list(tgv_value ORDER BY component_name))`, not the
  `GEV_ACC_TYPE` cast first planned: that round trip moved ~10% of one-component
  totals by one ulp (see Risks). Views are rebuilt by `define_trait()`. `restore_pop()` needs no
  migration, because a pre-0.74 file is refused for `ind_tbv` anyway.
- **`component_names` default** becomes `'total'`:
  - in the DDL (`R/open_pop.R:325`) and `R/define_phenotype.R:550`;
  - the roxygen (:78-79) drops "Reserved";
  - the `.sm_col()` text (`R/schema.R:422-423`) is rewritten;
  - `define_phenotype()` validates each listed name against
    `c(TGV_COMPONENT_NAMES, "total")`, refuses `"total"` mixed with other names, and
    refuses duplicates.
- **The simple-phenotype precheck** (`R/add_phenotype_stages.R:148-163`) becomes "the
  trait has at least one term, any owner" (`.gev_read_model(conn, t)`). The message
  names `define_additive_effects()` or `define_genome_effect_terms()`, and no longer
  names `formula_tbv`.
- **The prevalence threshold uses the active-block rule** (§6A).
  - A new helper, `.ap_prevalence_genetic_var(conn, trait)`, sorts the trait's terms
    (all scopes, any owner) into kinds:
    - an order-1 `additive` term is `additive`;
    - an order-1 `dominance` term is `dominance`;
    - an order-2 term with two additive members is `additive_by_additive`;
    - anything else (indicator, A×D, order ≥ 3) is unsupported.
  - It sums the `line_name IS NULL` diagonal of each present kind from `trait_var_comp`.
  - It errors, naming `define_phenotype(thresholds = )`, if:
    - the model has an unsupported kind; or
    - a present kind has no stored diagonal.
  - It is used by both `.ap_check_prevalence()` (`:1152-1170`, the PLAN stage, before
    any RNG) and `.ap_liability_records()` (`:1194-1207`), so the two cannot disagree.
  - 3b adds the owner rule to the same helper.
- **The D7 contract** is unchanged in substance. The one surviving write becomes the
  `add_tgv()` replace in `.gev_write_tgv()`. It is retargeted in CLAUDE.md and
  `test-add_phenotype_failure_contract.R` (:21-34, :55-62, :199).

### 3a.5 Tests (3a)
- **Mechanical changes:**
  - `add_tbv(` → `add_tgv(` (108 sites).
  - Reads of `ind_tbv` / `tbv_value` → `ind_tgv |> filter(component_name ==
    "additive")` / `tgv_value`. A test-only helper `tgv_additive(pop, trait)` in
    `helper-pop.R` keeps the diffs small.
  - `type =` → `weight_type =` in `test-add_tbv_index.R`.
- **File renames with `git mv`:**
  - `test-add_tbv.R` → `test-add_tgv_breeding_value.R`, keeping the first-principles
    oracle at :16-40 (T2);
  - `test-add_tbv_index.R` → `test-add_tgv_index.R`.
- **Determinism:** `test-group-contributor-determinism.R` is retargeted to `ind_tgv`
  (PH5 form: total and a listed component).
- **New gates:**
  - T1: part of PH8's grep, run after 3c.
  - T3, T4: mixed model through `define_genome_effect_terms()` + `ad_terms(coding =
    "cockerham", d = )`, one A×A term from raw rows, and one indicator from
    `genotype_terms()`.
  - T5, T6 and T7 (additive, total, and both indices coexisting, with
    `component_name` stored).
  - T8 (archive stamps; `remove_rows()` on `ind_tgv`; `restore_pop()` refuses
    `ind_tbv`, `order1_*`, and an old `ind_true_index`).
  - T9.
  - PH1, PH3, PH4 and PH5: dominance fixtures come from the writer, because no generator
    writes dominance before step 5.
  - PH6.
  - PH7 parts that need no owner rule:
    - A + D target sum. The generated dominance terms are planted through the internal
      `.ge_write_terms(allow_reserved_owner = TRUE)`, test-only, as the
      generator==writer test already does.
    - A stored A×A target with no A×A terms is left out.
    - A model with an unsupported kind errors.
    - Composite with `prevalence` errors.
    - Explicit `thresholds` works.
    - Each error leaves `ind_phenotype` and `phenotype_random_effects` unchanged.
- **Parity goldens:** `helper-parity.R` calls `add_tbv()` once, and `add_tgv()` draws no
  RNG, so the goldens must still pass unchanged. If they move, stop and investigate.

### 3a.6 Docs (3a)
- **CLAUDE.md:** the D7 sentence, "One evaluator" (it drops "`add_tbv()` reads only
  reserved-owner…"), design principle 4's action list, the naming-rule-2 example
  (`tbv_value`), the Roadmap line (consolidation is done), and the Two-Layer bullet
  "TBVs in `ind_tbv`".
- **Both skills:**
  - the `add_tbv` section is removed, and `add_tgv` is rewritten with the index args;
  - in the schema skill, `ind_tbv` is removed, `ind_true_index` gains
    `component_name`, and the `component_name` vocabulary and the
    `component_names` default are updated.
- **Other docs:** `README.md`, `package_summary.md`, `dev/package_summary/*`,
  `dev/benchmarks/*` (3 files), the introduction vignette and the swine vignette script
  (15 `ind_tbv`, 12 `add_tbv`, 13 `tbv_value`).
- **Bookkeeping:** NEWS 0.74.0 (breaking: older databases are not readable), DESCRIPTION,
  and `devtools::document()`.

---

## 3b — Q21: generated means calibrated (0.74.1)

### 3b.1 `define_additive_effects()` loses `effects` and `scale_to_target`
The dead code to remove (`R/define_additive_effects.R`):
- the checks (:313-316, :331-333, :434-440);
- the A19 refusal (:319-326);
- the manual branch (:471-477) and its "not calibrated" message (:546-547);
- the `need` modes `"none"` / `"sigma"` (:369-371, and `.dae_resolve_target()` :850-867);
  every `need == "target"` test becomes unconditional;
- the unscaled union warning (:414-418);
- the `.dae_check_realised()` first branch and its two parameters (:816-822);
- the roxygen at :14-15, :92-95, :163-165, :206-207 and :639-642.

New roxygen paragraph: the generator always samples and calibrates. Exact coefficients
(GWAS estimates, a QTL map) go through `define_genome_effect_terms()` + `ad_terms()`
under a user owner, and then the prevalence threshold needs `thresholds =`.

Error strings that name the removed arguments are reworded to point at the writer:
- `R/chr_meta_helpers.R:357,389-395` (sex-linked QTL);
- `R/add_phenotype_stages.R:1161-1165`.

### 3b.2 The owner rule (Q21 (a))
- `.ap_prevalence_genetic_var()` also errors when **any** active term of the trait has
  `effect_owner != 'generated'`, even with a stored target. The message names
  `thresholds =` and gives the reason: the target is not known to describe hand-written
  terms.
- The `define_phenotype(prevalence = )` roxygen replaces the 0.72.3 caveat
  (`R/define_phenotype.R:32-44`) with the rule.
- The `mean` roxygen (:24) states that `mean` is an intercept (§6A "Mean"), with the
  recipe for setting a realised base mean.

### 3b.3 `define_effect_cov_matrix()` refuses a block under generated terms
- This is only in the exported function, never in `.tvc_write_block()`, so
  `define_additive_effects(G = )` keeps working (§7.1).
- It refuses an `effect_name` in `GENETIC_EFFECT_NAMES` when any trait of the block has
  `generated` terms of that kind at the block's scope:
  - a `line_name = NULL` block is checked against terms whose origin `line_name` is NULL
    (common or parent-only);
  - a `line_name = "C"` block is checked against line-C-scoped terms.
- A line-C target written before line-C effects are generated still works (A17).
- The error gives the sequence: `remove_rows()` the old block, then
  `define_additive_effects(G = )`, which writes target and terms in one transaction.
- Step 2's "already stored" error in `.tvc_write_block()` is reworded to give the same
  sequence, not only `remove_rows()`.

### 3b.4 Migrating the 93 call sites (tests), by category
| Category | Sites | Migration |
|---|---|---|
| (a) Known coefficients for a TGV oracle (`test-add_tgv_breeding_value.R` 20, `test-genome-effects-eval.R` 17, `test-extract_allele_freq.R` 2, `test-add_tgv_index.R` 1 on Y loci) | ~40 | Use the test helper `with_additive_terms(pop, trait, loci, a, base_tbl, origin = NULL, owner = "custom")`. It wraps `define_genome_effect_terms(data.frame(locus_name, contrast_name = "additive", genome_value = a), base_tbl = )`, so the writer fills `center_value = p_base` and values match the generator's centring. Assertions do not change, because `add_tgv()` reads every owner. |
| (d) Filler effects while testing generator behaviour: centring, line/base precedence, `base_tbl` validation and injection, half bases, `parent_origin` validation, writer gates 34/42 | ~40 | Keep calling the generator, sampled and calibrated, with a stored target (`with_additive_target()` or `G =`). Assert on what the test is about (`center_value`, scope rows, errors), not on `genome_value == effects`. Where a value really must be known, use (a). |
| (b) Unscaled per-trait QTL sets for `method = "union"` (A7, R5) and `test-add_phenotype.R:320` | 6 | Add the test helper `plant_generated_additive()`, which calls `.ge_write_terms(allow_reserved_owner = TRUE)`, test-only, as the precedent at `test-genome-effects-writer.R:1199`. The prevalence fixture uses custom terms and now expects the owner-rule error. |
| (c) Tests of the removed arguments (A19 :466/468, :479, A22 :524, `test-define_additive_effects.R:316`) | 5 | Delete them. Add one test that `effects =` and `scale_to_target =` are "unused argument" errors. A22's RNG check keeps another refusal. |
| `scale_to_target = TRUE` / `effects = NULL` no-ops | 3 | Drop the argument. |

Outside `tests/`:
- swine vignette (6×) and README (2×): drop `scale_to_target = TRUE`. Also fix README
  :760's stale `base = "current_pop"`;
- the old `add_tbv.R` roxygen example: it is deleted in 3a, so its line-specific
  manual-effects example is rewritten in `add_tgv()`'s roxygen with the writer;
- `dev/benchmarks/benchmark_tgv_scale.R:63-67` and `benchmark_phenotype_scale.R:108`:
  move to the writer;
- the skill prose (`tidybreed-api` :319, :392).

### 3b.5 New gates (3b)
- The rest of PH7:
  - the Codex finding-1 reproduction: target 1, `ad_terms()` effects of 10 at ten
    loci, and `prevalence = 0.1`, which now errors;
  - stored target + any non-`generated` term errors;
  - after generation, `define_effect_cov_matrix()` on the trait's `additive` block
    errors, including after `remove_rows()` of the old target;
  - `remove_rows()` + `define_additive_effects(G = )` succeeds;
  - a line-C target before line-C generation succeeds.
- PH7's Q23 regression: the exported writer cannot write `generated`. It is already
  pinned at 0.73.2; keep it.
- "Unused argument" for `effects` / `scale_to_target`.

### 3b.6 Docs (3b)
- the skills (step list, writer route for exact coefficients);
- the main plan §7.1's "Signature from step 3" note is now true; mark it as built;
- gate A19 is marked moot;
- NEWS 0.74.1 (breaking) and DESCRIPTION.

---

## 3c — Q18: `formula_tgv` and the DSL arguments (0.74.2)

### 3c.1 Rename (no alias)
- The argument, and the `phenotype_meta` column (`R/open_pop.R:287-306`).
- `sql_utils.R:123` and `schema.R`.
- Internals: `.FORMULA_TBV_DSL_FUNS`, `.eval_formula_tbv()`,
  `.validate_formula_tbv()`, `.walk_formula_tbv_ast()`, `.substitute_tbv_ast()`,
  `.build_tbv_env()`, `.ap_materialize_tbvs()`, `.assemble_composite_tbv()` and the
  `tbv_kind` values. The `.tbv_*` placeholders become `.tgv_*`.
- Use `git mv R/contributor_tbv.R R/contributor_tgv.R` (and its Rd topic).
- The roxygen wording "composite TBV" becomes "composite genetic value".
- `restore_pop()` refuses `phenotype_meta.formula_tbv`.

### 3c.2 Parser (`R/formula_helpers.R:156-292`)
- Each of `self` / `dam` / `sire` takes exactly one positional trait, plus an optional
  named `component =` (a string in `TGV_COMPONENT_NAMES` or `"total"`, default
  `"total"`).
- `group_sum` / `group_mean` take positional `trait, col`, plus optional named
  `component =` and `table =` (default `"ind_meta"`). The positional third argument is
  removed.
- Anything else errors in `define_phenotype()`, naming the call. This covers an unknown
  named argument, an extra positional argument, or a non-string or non-symbol value.
- A ref becomes `list(trait, type, col, table, component, placeholder)`.
- `.substitute_tgv_ast()` matches on **every** ref field. Today it ignores `table`, so
  two group terms differing only in table collide. That is a latent bug, fixed here, and
  a test covers it.
- `col` and `table` pass `validate_sql_identifier()` (as `define_phenotype.R:517,523`).
  `define_phenotype()` checks that the table exists and has the column. The current
  "validated at `add_phenotype()` time" message goes.
- `.build_tgv_env()` passes `component` to `.tgv_by_id()` / `.group_mate_tgv()`.
- The `formula_tgv` roxygen documents `component =` and `table =`.

### 3c.3 Gates (3c)
- **PH2 in full:**
  - the default total;
  - `dam(WWM, component = "additive")`;
  - `component = "bogus"` errors;
  - `table =` reads the named table;
  - an unknown table or column, a non-identifier, and a positional third argument each
    error before any write.
- The duplicate-ref fix.
- **PH8 grep:** `ind_tbv`, `add_tbv`, `tbv_value`, `id_tbv`, `formula_tbv`, `order1_`,
  `.gev_reserved_additive` and `.gev_warn_tbv_stale` appear nowhere in `R/`, `tests/`,
  `man/`, `vignettes/`, `dev/`, `NAMESPACE`, `CLAUDE.md`, `README.md`, `_pkgdown.yml`,
  `package_summary.md` or the skills. Past NEWS, closed plans, `tools/quarto/legacy/`
  and `docs/` are excepted.
- 36 `formula_tbv` uses in `test-formula_phenotype.R` are renamed. The existing `table =`
  test (:379) is kept in its named form.
- **Bookkeeping:** NEWS 0.74.2, DESCRIPTION, docs.

---

## Risks and things to watch
- **Phenotype values change only where the model is not additive-only.** Every existing
  phenotype test uses additive-only `generated` models, where total = additive. PH1
  asserts `expect_identical()` there. Any other change in an existing phenotype test is
  a bug.
- **The exact view changes `ind_tgv_total` bits.** With one component per row the
  DECIMAL round trip is exact for |x| < 1e20 at 18 decimals. A value with more than 18
  fractional digits rounds at ~1e-18, which is below double precision for |x| ≥ ~1e-2
  but not for tiny values. Check that parity goldens (1e-8) and `expect_identical`
  tests still pass. If the round trip moves a small additive-only total, read the single
  row directly when the trait has one component, or accept and document it. Measure
  before choosing.
- **`add_tgv()` in stage 1 now evaluates every component** of every contributor trait.
  It uses the same single evaluator pass. Run `dev/benchmarks/benchmark_phenotype_scale.R`
  before and after 3a, and report the result.
- **`.gev_require_contribution()`** errors exactly as before (a female with Y-only
  terms). No change.
- **PH7's A + D case plants `generated` dominance terms through the internal engine**
  until step 5 writes real ones. That is test-only, and it is noted in the test.

## Verification (each sub-step)
1. `pkgload::load_all()`, then `testthat::test_file()` on the touched files.
2. Full suite with `NOT_CRAN=true`, in the background: 0 failures. The warning count is
   compared with 0.73.2 (11), and any new warning is explained.
3. `devtools::document()` runs clean, and `pkgdown::check_pkgdown()` too.
4. The introduction vignette, purled and sourced under `load_all()`; the swine script
   sourced.
5. A mutation spot-check per sub-step:
   - 3a: revert the exact view to `SUM()` → PH5 must fail;
   - 3b: drop the owner check → PH7's Codex reproduction must fail;
   - 3c: drop the `table` field from ref matching → the duplicate-ref test must fail.
6. After each sub-step: the main plan's "As built" paragraph, then commit + push, then
   **wait for the user's review**. Add the sub-step's section to
   `import_qtl_effect_methods_phase_3.md` (results, one file for all of step 3).

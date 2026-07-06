# Fix Plan: `add_ebv()` — Schema Inconsistencies & Maternal Model Support

> Revised after `plans/fix_add_ebv_critique.md`. All corrections in that
> critique were verified against the current code and accepted. Scope
> decisions made in response to the critique are recorded inline as
> **Decision:** notes.

## Background

`add_ebv()` and its helpers in `blupf90_helpers.R` were written before the
two-layer phenotype design (v0.31.0) was fully settled. As a result they
contain **two** schema mismatches that cause hard crashes (not one — see Bug 2
and Bug 2b). Additionally, `add_ebv()` has no support for composite (maternal)
phenotypes — users currently have to manually build the BLUPF90 parameter
file, run the binaries, and insert results themselves (as seen in the swine
vignette).

The source of truth for composite phenotype structure for the BLUPF90 path in
this plan is `phenotype_meta.formula_tbv` (e.g. `"WWD + dam(WWM)"`), parsed via
the existing `.walk_formula_tbv_ast()`.

`phenotype_components` (written by `define_phenotype(..., components = ...)`)
is **not** unused — `add_phenotype()` reads it directly and
`.assemble_composite_tbv()` (`R/add_phenotype.R`) is a full implementation over
it (self/dam/sire/group contributors, weights, covariates, SGE aggregation).

**Decision:** `add_ebv()` in this pass only resolves components via
`formula_tbv`. Phenotypes defined only via `phenotype_components` (no
`formula_tbv`) are explicitly unsupported in `add_ebv()` for now and produce a
clear error telling the user to add a `formula_tbv` equivalent. Do not describe
`phenotype_components` as unused anywhere in code or docs — it is live and
owned by `add_phenotype()`.

---

## Scope for This Pass

**Decision:** ship a narrow, verifiable first version rather than general
composite support:

- At most **one** phenotype per `add_ebv(..., software = "blupf90")` call may
  be composite (maternal). Mixing a composite phenotype with other phenotypes
  (composite or simple) in the same BLUPF90 call is out of scope and errors.
- Composite phenotypes are resolved only via `formula_tbv`, and only the
  `self` + `dam` contributor roles are supported. Any `formula_tbv` containing
  `sire`, `group_sum`, `group_mean`, a scalar constant, or more than one `dam`
  role errors with a clear message — it does not silently degrade to a partial
  or wrong model.
- `formula` (derived, e.g. FCR-style ratio) phenotypes are out of scope for
  both BLUPF90 and parent-average modes in this pass — `add_ebv()` errors on
  them. (`derived_formula` phenotypes have no independent residual variance by
  design and need their own residual/R strategy; see critique #7.)
- Parent-average (`parent_avg = TRUE`) mode gets the same self+dam-only,
  formula_tbv-only restriction, for consistency and because it shares the
  resolver with the BLUPF90 path.

This intentionally defers: `sire`/`group` contributor roles in `add_ebv()`,
multi-phenotype BLUPF90 calls that mix composite and simple traits, and
`phenotype_components`/`formula`-based composites. A follow-up plan can
generalize once the `renum.par`/`solutions.orig` mapping for the single
self+dam case is validated end-to-end.

**Decision (back-compat):** the `add_ebv()` argument is renamed to
`phenotype_name`, but `trait_name` is kept as a **permanent silent alias** (no
deprecation warning) — both names route to the same argument indefinitely, to
avoid breaking existing scripts (including `plans/run_blupf90_manual.md`).

---

## Bugs to Fix

### Bug 1 — `load_effect_cov()` ✅ ALREADY FIXED

`write_renum_par()` already calls `load_phenotype_cov()` (residual) and
`load_trait_cov()` (gen_add). No action needed.

---

### Bug 2 — `ind_phenotype.trait_name` column doesn't exist (hard crash)

**File:** `R/blupf90_helpers.R:119`, also lines 123, 127, 129, 133

```r
# CURRENT (broken)
"SELECT id_ind, trait_name, AVG(pheno_value) FROM ind_phenotype
 WHERE trait_name IN (...)"
```

`ind_phenotype` stores `phenotype_name`, not `trait_name` (confirmed in
`R/define_trait.R`, `ind_phenotype` DDL). This query fails at runtime with a
DuckDB column-not-found error.

**Fix:** change every reference to `trait_name` in the `ind_phenotype` query to
`phenotype_name`. The pivot logic at lines 127/129/133 must also use
`phenotype_name`. For the maternal model, the query uses the `phenotype_name`
("WW"), not the component trait names ("WWD", "WWM") — there is only one
observation column in the data file.

---

### Bug 2b — `trait_effects.trait_name` column doesn't exist (second hard crash)

**File:** `R/blupf90_helpers.R:77-84` (`build_data_file()`)

```r
# CURRENT (broken)
effects_df <- DBI::dbGetQuery(
  pop$db_conn,
  paste0("SELECT trait_name, effect_name, effect_class, source_column, source_table ",
         "FROM trait_effects ",
         "WHERE trait_name IN (", ...
```

`trait_effects`'s real schema uses `phenotype_name` as its key column
(`R/define_trait.R:182-196`; PK is `(phenotype_name, effect_name)`), and effect
helpers (`define_effect_fixed_class()`, `define_effect_fixed_cov()`) write rows
keyed by `phenotype_name`. This query also crashes with a column-not-found
error, independently of Bug 2.

**Fix:** change the query to select `phenotype_name`, not `trait_name`, and
filter by the **observed phenotype name(s)** being evaluated (e.g. `"WW"`),
not the component trait names (`"WWD"`, `"WWM"`) — fixed effects are attached
to the observed phenotype model, not to individual genetic components. In
`write_renum_par()`, per-observation fixed-effect column presence must be
checked against `effects_df$phenotype_name`, not component trait names.

---

### Bug 3 — G matrix passed wrong names for maternal model

**File:** `R/blupf90_helpers.R:282-283`

```r
# CURRENT — function names are correct but argument is wrong for maternal
R_mat <- load_phenotype_cov(pop, "residual", trait)   # trait = "WW" ✓
G_mat <- load_trait_cov(pop,     "gen_add",  trait)   # trait = "WW" ✗
```

`load_phenotype_cov(pop, "residual", "WW")` is correct — residuals live in
`phenotype_var_comp` keyed by phenotype name.

`load_trait_cov(pop, "gen_add", "WW")` returns `NULL` for a maternal phenotype
because `"WW"` is not in `trait_var_comp` — only `"WWD"` and `"WWM"` are.
The G matrix must use the component trait names extracted from `formula_tbv`.

**Fix:** `write_renum_par()` needs two separate arguments — `phenotype_name`
(for R matrix) and `component_trait_names` (for G matrix):

```r
R_mat <- load_phenotype_cov(pop, "residual", phenotype_name)
G_mat <- load_trait_cov(pop,     "gen_add",  component_trait_names)
```

For simple traits, `component_trait_names` equals `phenotype_name` and both
calls work identically to today.

---

### Bug 4 — `add_ebv()` parameter `trait_name` is ambiguous; validates wrong table

**File:** `R/add_ebv.R:132`, `R/add_ebv.R:152-208`

The function accepts `trait_name` and validates against `trait_meta`. But
`ind_phenotype` is keyed by `phenotype_name`, and for composite phenotypes the
user's entry point is the phenotype name (e.g. `"WW"`), not the underlying
genetic trait names (`"WWD"`, `"WWM"`). For simple traits the two are identical,
so this silently works; for composite traits it breaks.

**Fix:** rename the parameter to `phenotype_name` throughout `add_ebv()`, and
accept `trait_name` as a permanent silent alias (see Scope decision above —
both resolve to the same internal variable; supplying both is an error).
Validate the supplied names against `phenotype_meta.phenotype_name` instead of
`trait_meta.trait_name`. Component trait names are then derived internally by
parsing `formula_tbv` (see Feature section).

The stored column in `ind_ebv` remains `trait_name` (genetic component) — EBVs
are stored as `"WWD"` and `"WWM"`, not `"WW"`.

---

### Bug 5 — `@param update_covars` references wrong table name ✅ ALREADY FIXED

The doc at `R/add_ebv.R:58` already says `trait_var_comp`. No action needed.

---

### Bug 6 — Empty phenotype filter silently falls back to all records

**File:** `R/blupf90_helpers.R` (`build_data_file()`), `R/add_ebv.R`

`add_ebv()` warns when a `phenotype` filter matches no records
(`R/add_ebv.R` — `pheno_ids` branch), but `build_data_file()` only applies
`pheno_ids` when `length(pheno_ids) > 0`. If the filter matches zero rows,
`pheno_clause` becomes empty and **all** phenotype records are included —
silently contradicting the filter the user supplied.

**Fix:** if `phenotype` is supplied and matches zero rows, stop before
entering BLUPF90 mode (upgrade the existing warning to an error in that
specific case) rather than passing an empty filter through.

---

## Feature: Maternal Model Support in BLUPF90

### Two distinct name-sets — do not conflate them

Per critique #2, `add_ebv()` must track two separate vectors instead of one
`trait`:

- `phenotype_name` — the observed phenotype(s). Drives `data.txt` observation
  column(s), `trait_effects` fixed-effect lookups, and the R (residual)
  matrix. Exactly one value in this pass (see Scope).
- `component_trait_names` — the genetic component trait(s) resolved from
  `formula_tbv` (e.g. `c("WWD", "WWM")`). Drives the G matrix, `ind_ebv`
  storage, and `eval_nums`/delete logic. For a simple trait,
  `component_trait_names == phenotype_name` (length 1).

In `write_renum_par()`, `n_traits` for the **data file / R matrix** side must
be based on `phenotype_name` (one BLUPF90 "trait" column per observed
phenotype), while the G matrix dimension is based on `component_trait_names`.
These are not interchangeable, and (per Scope) this pass only needs to support
`length(phenotype_name) == 1` with `length(component_trait_names) %in% c(1, 2)`
(simple, or self+dam maternal).

### Parsing `formula_tbv` — reuse existing `.walk_formula_tbv_ast()`

**Do not write a new helper.** `R/formula_helpers.R` already has a full AST
walker (`.walk_formula_tbv_ast()`) that returns exactly what is needed: a list
of `trait_refs`, each with `$trait` (the trait name) and `$type` (`"self"`,
`"dam"`, `"sire"`, `"group_sum"`, `"group_mean"`).

After fetching `formula_tbv` from `phenotype_meta`, call:

```r
walk_res        <- .walk_formula_tbv_ast(parse(text = formula_tbv)[[1]])
role_counts     <- table(vapply(walk_res$trait_refs, `[[`, character(1), "type"))
unsupported     <- setdiff(names(role_counts), c("self", "dam"))

if (length(unsupported) > 0 || isTRUE(walk_res$has_scalar_constant))
  stop("add_ebv(software = 'blupf90') only supports formula_tbv expressions ",
       "using self + dam roles. Phenotype '", phenotype_name,
       "' uses unsupported role(s): ", paste(unsupported, collapse = ", "),
       call. = FALSE)

direct_traits    <- unique(vapply(
  Filter(function(r) r$type == "self", walk_res$trait_refs),
  `[[`, character(1), "trait"))
maternal_traits  <- unique(vapply(
  Filter(function(r) r$type == "dam",  walk_res$trait_refs),
  `[[`, character(1), "trait"))

if (length(direct_traits) > 1 || length(maternal_traits) > 1)
  stop("add_ebv(software = 'blupf90') supports at most one direct (self) and ",
       "one maternal (dam) component trait per phenotype.", call. = FALSE)

is_maternal      <- length(maternal_traits) == 1
component_trait_names <- unique(c(direct_traits, maternal_traits))
```

For a simple phenotype with no `formula_tbv`, `component_trait_names <-
phenotype_name` and `is_maternal <- FALSE`.

For a phenotype with `phenotype_meta.formula` set instead of `formula_tbv`
(derived/ratio phenotypes), or with only `phenotype_components` rows and no
`formula_tbv`, `add_ebv()` errors per the Scope decision above — these are not
resolved in this pass.

### Changes to `add_ebv()` main function

1. **Rename** `trait_name` → `phenotype_name` in signature and all internal
   references; keep `trait_name` as a silent alias argument (Scope decision).
2. **Validate** supplied names against `phenotype_meta.phenotype_name` (not
   `trait_meta`). Error if `length(phenotype_name) > 1` and any of them
   resolves to a composite (see Scope — one composite phenotype at a time;
   multiple simple phenotypes in one call remain fine since they degenerate to
   today's behavior).
3. **Parse `formula_tbv`** early to resolve `component_trait_names` and
   `is_maternal`, applying the self+dam-only validation above.
4. **`eval_nums` computation** must loop over `component_trait_names` (e.g.
   `"WWD"`, `"WWM"`), not `phenotype_name` (`"WW"`), because `ind_ebv` is
   keyed by `trait_name` (genetic component).
5. **`overwrite_trait` / `delete_all`** delete logic also uses
   `component_trait_names` when touching `ind_ebv`.
6. **Pass `phenotype_name` and `component_trait_names`** as separate arguments
   into `ebv_blupf90()` and `ebv_parent_avg()`.
7. **Zero-row phenotype filter** (Bug 6): stop before dispatching to either
   mode if `phenotype` was supplied and matched zero rows.

### Changes to `ebv_blupf90()` signature

Currently receives a single `trait_name` argument. Must be updated to receive
both `phenotype_name` (scalar in this pass) and `component_trait_names`
(length 1 or 2) as distinct arguments and route them to the correct downstream
helpers.

### Changes to `build_data_file()`

- Fix Bug 2: use `phenotype_name` in the `ind_phenotype` query (one observation
  column for "WW").
- Fix Bug 2b: use `phenotype_name` in the `trait_effects` query, filtered by
  observed phenotype name(s), not component trait names.
- For a maternal model, add a dam ID column fetched from `ind_meta.id_parent_2`.
  Missing dam IDs (founders) get `"0"` (BLUPF90 missing indicator).
- Use `sql_in_list()` (`R/sql_utils.R`) instead of ad hoc
  `paste0("'", x, "'", collapse = ", ")` for all `IN (...)` clauses touched in
  this pass, and `gsub("'", "''", x)` (or `format_sql_value()`) for any other
  string literal construction. This is not a full SQL-hardening pass — just do
  not add new unescaped string interpolation while touching this file.

Column layout for maternal model:
```text
col 1: mu (intercept)
col 2: id_ind (animal ID — for direct random effect)
col 3: [fixed effects if any, looked up via phenotype_name]
col 4: WW observation
col 5: id_parent_2 (dam ID — for maternal random effect)
```

### Changes to `write_renum_par()`

- Accept both `phenotype_name` (scalar, drives R matrix / TRAITS / data
  column) and `component_trait_names` (drives G matrix) as arguments.
- R matrix: `load_phenotype_cov(pop, "residual", phenotype_name)`
- G matrix: `load_trait_cov(pop, "gen_add", component_trait_names)`
- For the maternal model, add an `OPTIONAL mat` line to the `RANDOM animal`
  block and add the dam ID column to `FIELDS_PASSED TO OUTPUT` tracking in
  `col_map`.
- `animal_effect_num` stays `1 + n_fixed_effs + 1`. For maternal, the maternal
  effect is `animal_effect_num + 1`.

**Prerequisite investigation (not yet resolved by this plan — per critique
item 3):** before implementing the maternal branch, produce a synthetic
`renum.par` for a small self+dam example and confirm against `renumf90`
documentation/behavior:

- The exact `OPTIONAL mat` line syntax, including which column holds the dam
  ID.
- Whether a separate `EFFECT` block for the dam ID column is required, or
  whether `OPTIONAL mat` alone modifies the existing `RANDOM animal` block.
- Whether `FIELDS_PASSED TO OUTPUT` must include the dam column, and whether
  that changes `solutions.orig` layout.
- How `renumf90` numbers direct vs. maternal effects in `solutions.orig` (to
  confirm the `animal_effect_num + 1` assumption).

This should be validated against a real `renumf90` run (or authoritative
BLUPF90 docs) before `write_renum_par()`'s maternal branch is written, not
inferred from the direct-effect case.

### Changes to `parse_blupf90_solutions()`

For a maternal model, BLUPF90 outputs two sets of animal solutions:
- `effect_num == animal_effect_num` → direct EBVs → stored as `trait_name = component_trait_names[direct]` (e.g. `"WWD"`)
- `effect_num == animal_effect_num + 1` → maternal EBVs → stored as `trait_name = component_trait_names[maternal]` (e.g. `"WWM"`)

`eval_nums` carries entries for both component trait names (resolved before
this call).

For simple phenotypes, only `effect_num == animal_effect_num` is parsed —
unchanged.

### Changes to `ebv_parent_avg()` — self+dam composite support only

**File:** `R/add_ebv.R:295` (not `blupf90_helpers.R` — that file reference in
the original plan draft was stale).

The current implementation receives `trait_name` and looks it up directly in
`ind_ebv`. For a composite phenotype like "WW" (stored under "WWD"/"WWM"),
this returns nothing and silently produces all-NA EBVs.

**Fix:** resolve `component_trait_names` before entering `ebv_parent_avg()`
(same resolution done for the BLUPF90 path, same self+dam-only restriction —
see Scope) and pass them in. The function then:

1. Looks up parent EBVs for each component trait separately in `ind_ebv`.
2. Computes parent average for each component independently.
3. Returns one EBV row per (animal × component trait).

| Phenotype | `formula_tbv` | PA operates on |
|---|---|---|
| `"WW"` | `"WWD + dam(WWM)"` | "WWD", "WWM" via `.walk_formula_tbv_ast()` |
| `"ADG"` | none (simple) | "ADG" directly |

Ratio/derived phenotypes (e.g. FCR via `formula = "ADFI / ADG"`) and
`sire`/`group_sum`/`group_mean` roles are out of scope for `ebv_parent_avg()`
in this pass too (see Scope) — `add_ebv(parent_avg = TRUE)` errors on them with
the same message as the BLUPF90 path.

---

## Summary of All Changes

| File | Change |
|---|---|
| `R/add_ebv.R` | Rename `trait_name` → `phenotype_name` (keep `trait_name` as permanent silent alias); validate against `phenotype_meta`; parse `formula_tbv` to resolve `component_trait_names` with self+dam-only validation; enforce at-most-one-composite-phenotype-per-call; pass both `phenotype_name` and `component_trait_names` downstream; route `eval_nums`/delete logic to component trait names; error before dispatch on zero-row `phenotype` filter |
| `R/blupf90_helpers.R:ebv_blupf90()` | Accept `phenotype_name` + `component_trait_names` as distinct args; route each to correct downstream helper |
| `R/blupf90_helpers.R:build_data_file()` | Fix `ind_phenotype` query (Bug 2) and `trait_effects` query (Bug 2b): `trait_name` → `phenotype_name`; fixed effects lookup uses `phenotype_name`; for maternal model add dam ID column; use `sql_in_list()`/escaping for new/touched queries |
| `R/blupf90_helpers.R:write_renum_par()` | Accept `phenotype_name` + `component_trait_names`; R matrix via `load_phenotype_cov(phenotype_name)`; G matrix via `load_trait_cov(component_trait_names)`; for maternal add `OPTIONAL mat` block (syntax to be confirmed per prerequisite investigation above) |
| `R/blupf90_helpers.R:parse_blupf90_solutions()` | For maternal model parse two effect blocks, map to component trait names |
| `R/add_ebv.R:ebv_parent_avg()` | Accept `component_trait_names` (self+dam only); look up PA for each component separately; return one row per animal × component |
| `man/add_ebv.Rd`, `man/build_data_file.Rd`, `man/write_renum_par.Rd`, `man/write_meta_file.Rd`, `man/parse_blupf90_solutions.Rd`, `man/update_covars_from_blupf90.Rd` | Regenerate via roxygen after any signature change to the functions they document |

---

## What Does NOT Change

- `ind_ebv` schema — `trait_name` column stays as-is; EBVs are always stored
  under genetic component names (`"WWD"`, `"WWM"`), never under phenotype name
  (`"WW"`)
- `load_trait_cov()` and `load_phenotype_cov()` — no changes to these functions
- `.walk_formula_tbv_ast()` and `.extract_all_symbols()` in `formula_helpers.R`
  — reused as-is, no changes
- Simple phenotype path — all existing call sites using `add_ebv("ADG")`,
  `add_ebv("AP")` etc. continue to work without modification; `component_trait_names`
  just equals the phenotype name and every code path degenerates to current
  behaviour

## Explicitly Out of Scope for This Pass (see Scope decision)

- `phenotype_components`-based composites in `add_ebv()` (still fully
  supported in `add_phenotype()` — unrelated)
- `sire`, `group_sum`, `group_mean` contributor roles in `add_ebv()`
- `formula` (derived/ratio, e.g. FCR) phenotypes in `add_ebv()`
- Multiple composite phenotypes, or mixed composite + simple phenotypes, in a
  single `add_ebv()` call
- A deprecation warning for `trait_name` (kept as a silent permanent alias
  instead)

All of the above error clearly rather than silently producing a wrong or
partial model. A follow-up plan can lift these restrictions once the
self+dam `renum.par`/`solutions.orig` mapping is validated in practice.

---

## Test Coverage

Tests that don't require BLUPF90 binaries to be installed:

- `build_data_file()` uses `phenotype_name` and no longer crashes on
  `ind_phenotype` (Bug 2) or `trait_effects` (Bug 2b).
- Fixed effects are loaded from `trait_effects.phenotype_name`, filtered by
  observed phenotype name.
- A zero-row `phenotype` filter errors instead of silently including all
  records (Bug 6).
- Component resolver: simple phenotype, `formula_tbv` self+dam maternal,
  `formula_tbv` with an unsupported role (`sire`/`group_sum`/`group_mean`) →
  clear error, `phenotype_components`-only phenotype → clear error, `formula`
  (derived) phenotype → clear error.
- `ebv_parent_avg()` returns component rows (both "WWD" and "WWM") for a
  maternal phenotype.
- `write_renum_par()` emits the expected maternal `renum.par` for a small
  synthetic case (once the prerequisite `OPTIONAL mat` investigation above is
  done).
- `parse_blupf90_solutions()` maps synthetic direct and maternal solution rows
  to the correct component trait names.
- `add_ebv(trait_name = "ADG", ...)` (named, old argument name) still works
  identically to `add_ebv(phenotype_name = "ADG", ...)`.
- `add_ebv()` errors clearly when asked to evaluate two composite phenotypes,
  or a composite + simple phenotype, in one BLUPF90 call.

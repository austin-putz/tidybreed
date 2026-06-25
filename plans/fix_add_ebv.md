# Fix Plan: `add_ebv()` — Schema Inconsistencies & Maternal Model Support

## Background

`add_ebv()` and its helpers in `blupf90_helpers.R` were written before the
two-layer phenotype design (v0.31.0) was fully settled. As a result they
contain several hard crashes and schema mismatches. Additionally, `add_ebv()`
has no support for composite (maternal) phenotypes — users currently have to
manually build the BLUPF90 parameter file, run the binaries, and insert results
themselves (as seen in the swine vignette).

The source of truth for composite phenotype structure is
`phenotype_meta.formula_tbv` (e.g. `"WWD + dam(WWM)"`). `phenotype_components`
is an unused artifact and should not be relied on.

---

## Bugs to Fix

### Bug 1 — `load_effect_cov()` doesn't exist (hard crash)

**File:** `R/blupf90_helpers.R:250-251`

```r
# CURRENT (broken)
R_mat <- load_effect_cov(pop, "residual", trait)
G_mat <- load_effect_cov(pop, "gen_add",  trait)
```

`load_effect_cov()` is not defined anywhere in the package. The correct
functions are `load_trait_cov()` and `load_phenotype_cov()`, which already exist
in `R/define_effect_cov_matrix.R`.

**Fix:** replace both calls (see Bug 3 for correct routing).

---

### Bug 2 — `ind_phenotype.trait_name` column doesn't exist (hard crash)

**File:** `R/blupf90_helpers.R:113`, also lines 127 and 129

```r
# CURRENT (broken)
"SELECT id_ind, trait_name, AVG(pheno_value) FROM ind_phenotype
 WHERE trait_name IN (...)"
```

`ind_phenotype` stores `phenotype_name`, not `trait_name` (confirmed in
`R/define_trait.R:214`). This query fails at runtime with a DuckDB
column-not-found error.

**Fix:** change every reference to `trait_name` in the `ind_phenotype` query to
`phenotype_name`. The pivot logic at lines 127/129 must also use `phenotype_name`.

---

### Bug 3 — Residual variance looked up from wrong table

**File:** `R/blupf90_helpers.R:250`

Even if Bug 1 were fixed by renaming to `load_trait_cov()`, the residual call
would still be wrong. Residual variances are stored in `phenotype_var_comp`
(routed there by `define_effect_cov_matrix(..., "residual", ...)`), not in
`trait_var_comp`. The column names also differ
(`phenotype_name_1`/`phenotype_name_2` vs `trait_name_1`/`trait_name_2`), so
`load_trait_cov()` would always return `NULL` for `"residual"` and the function
would stop with "Residual covariance matrix not found."

**Fix:** after parsing `formula_tbv` (see Feature below), route each matrix to
the correct function and table:

```r
# G matrix — component trait names from formula_tbv → trait_var_comp
G_mat <- load_trait_cov(pop, "gen_add", component_trait_names)

# R matrix — phenotype name → phenotype_var_comp
R_mat <- load_phenotype_cov(pop, "residual", phenotype_name)
```

For simple traits where `phenotype_name == trait_name` and no `formula_tbv` is
set, `component_trait_names` is just `phenotype_name` and both lookups still
work correctly.

---

### Bug 4 — `add_ebv()` parameter `trait_name` is ambiguous; validation checks wrong table

**File:** `R/add_ebv.R:19`, `R/add_ebv.R:96`, `R/add_ebv.R:116`,
`R/add_ebv.R:165-172`

The function currently accepts `trait_name` and validates it against
`trait_meta`. But `ind_phenotype` is keyed by `phenotype_name`, and for
composite phenotypes the user's entry point is the phenotype name (e.g. `"WW"`),
not the underlying genetic trait names (`"WWD"`, `"WWM"`). For simple traits the
two are identical, so this silently works; for composite traits it breaks.

**Fix:** rename the parameter from `trait_name` to `phenotype_name` throughout
`add_ebv()`. Validate the supplied names against `phenotype_meta.phenotype_name`
instead of `trait_meta.trait_name`. Component trait names are then derived
internally by parsing `formula_tbv`.

The stored column in `ind_ebv` remains `trait_name` (genetic component) — EBVs
are stored as `"WWD"` and `"WWM"`, not `"WW"`.

---

### Bug 5 — `@param update_covars` references wrong table name (doc bug)

**File:** `R/add_ebv.R:50`

> "write estimated variance components back to `trait_effect_cov`"

The table is `trait_var_comp`. `update_covars_from_blupf90()` is a stub anyway,
but the doc pointer is wrong.

**Fix:** change `trait_effect_cov` → `trait_var_comp` in the roxygen comment.

---

## Feature: Maternal Model Support in BLUPF90

### Design

`add_ebv()` will parse `phenotype_meta.formula_tbv` minimally — enough to detect
model type and extract component trait names. No full DSL evaluation is needed;
two regex checks cover all current cases:

```r
formula_tbv <- # fetched from phenotype_meta for this phenotype

direct_traits   <- # plain names:     "WWD" from "WWD + dam(WWM)"
maternal_traits <- # dam(...) names:  "WWM" from "WWD + dam(WWM)"

is_maternal <- length(maternal_traits) > 0
is_simple   <- !is_maternal  # no dam() or sire() terms
```

A helper (e.g. `.parse_formula_tbv()`) should be extracted to a shared location
so `add_phenotype()` and `add_ebv()` can reuse the same parsing logic.

### Changes to `build_data_file()`

For a maternal model, the data file needs an extra column: the dam ID of each
animal. This must be fetched from `ind_meta.id_parent_2`.

Column layout for maternal model:
```
col 1: mu (intercept)
col 2: id_ind (animal ID — for direct random effect)
col 3: [fixed effects if any]
col 4: WW observation
col 5: id_parent_2 (dam ID — for maternal random effect)
```

Missing dam IDs (founders, animals with unknown dam) get `"0"` (BLUPF90
missing indicator).

### Changes to `write_renum_par()`

For a maternal model, the BLUPF90 parameter file needs:

1. **G matrix**: 2×2 from `trait_var_comp` for `c("WWD", "WWM")` via
   `load_trait_cov()`:
   ```
   G = [var_WWD    cov_WWD_WWM]
       [cov_WWD_WWM var_WWM   ]
   ```

2. **R matrix**: scalar from `phenotype_var_comp` for `"WW"` via
   `load_phenotype_cov()`.

3. **RANDOM animal block** with `OPTIONAL mat` to activate the maternal effect:
   ```
   EFFECT
     <id_col per trait> cross alpha
   RANDOM
     animal
   OPTIONAL
     mat
   FILE
     pedigree.txt
   ...
   (CO)VARIANCES
     <G_WW 2x2>
   ```

For simple traits the `OPTIONAL mat` line is omitted and the G matrix is 1×1
(or n×n for multi-trait).

### Changes to `parse_blupf90_solutions()`

For a maternal model, BLUPF90 outputs two sets of animal solutions:
- `effect_num == animal_effect_num` → direct EBVs → stored as `trait_name = "WWD"`
- `effect_num == animal_effect_num + 1` → maternal EBVs → stored as `trait_name = "WWM"`

`eval_nums` must carry entries for both `"WWD"` and `"WWM"` (not `"WW"`).

For simple traits, only `effect_num == animal_effect_num` is parsed (current
behaviour unchanged).

### Changes to `add_ebv()` main function

The `eval_nums` computation and the delete-old-rows logic currently operate on
`trait_name`. For a maternal phenotype like `"WW"`, the internal trait names
`"WWD"` and `"WWM"` must be used for those lookups against `ind_ebv`, not the
phenotype name `"WW"`.

After parsing `formula_tbv`, resolve `component_trait_names` early and use those
names everywhere `ind_ebv` is touched.

---

## Summary of All Changes

| File | Change |
|---|---|
| `R/add_ebv.R` | Rename parameter `trait_name` → `phenotype_name`; validate against `phenotype_meta`; parse `formula_tbv` to get component trait names; route `eval_nums` / delete logic to component trait names |
| `R/add_ebv.R:50` | Fix doc: `trait_effect_cov` → `trait_var_comp` |
| `R/blupf90_helpers.R:build_data_file()` | Fix `ind_phenotype` query: `trait_name` → `phenotype_name`; for maternal model add dam ID column |
| `R/blupf90_helpers.R:write_renum_par()` | Replace `load_effect_cov()` with `load_trait_cov()` (G) and `load_phenotype_cov()` (R); for maternal model add `OPTIONAL mat` block |
| `R/blupf90_helpers.R:parse_blupf90_solutions()` | For maternal model parse two effect blocks and map to `"WWD"` / `"WWM"` |
| `R/blupf90_helpers.R` (new helper) | Add `.parse_formula_tbv(formula_tbv)` returning list of `direct_traits` and `maternal_traits` |
| `man/add_ebv.Rd` | Regenerate after roxygen changes |

---

## What Does NOT Change

- `ind_ebv` schema — `trait_name` column stays as-is; EBVs are always stored
  under genetic component names (`"WWD"`, `"WWM"`), never under phenotype name
  (`"WW"`)
- `load_trait_cov()` and `load_phenotype_cov()` — no changes to these functions
- `phenotype_components` — remains unused; no changes
- Simple trait path — all existing call sites using `add_ebv("ADG")`,
  `add_ebv("AP")` etc. continue to work without modification

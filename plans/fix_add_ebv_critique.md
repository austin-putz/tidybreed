# Critique: `plans/fix_add_ebv.md`

## Overall Assessment

The plan correctly identifies the main public/API mismatch: `add_ebv()` is still written around genetic `trait_name`, while observed records now live under `phenotype_name`. It also correctly catches the `ind_phenotype.trait_name` hard crash in `build_data_file()`.

However, the plan is not yet safe to implement as written. The biggest problems are:

1. It misses a second schema crash in `trait_effects`.
2. It treats observed phenotypes and genetic component traits as interchangeable in several BLUPF90 places where they are not.
3. It understates the model-design work needed for a maternal animal model.
4. It would break named `trait_name = ...` callers unless an alias/deprecation path is kept.
5. It incorrectly says `phenotype_components` is unused.

## High Severity Issues

### 1. The plan misses the `trait_effects.trait_name` schema crash

The plan focuses on `ind_phenotype.trait_name`, but `build_data_file()` also queries `trait_effects.trait_name`:

- Current code: `R/blupf90_helpers.R:77-84`
- Actual schema: `R/define_trait.R:181-196` defines `trait_effects(phenotype_name, effect_name, ...)`
- Effect helpers insert by `phenotype_name`: `R/define_effect_fixed_class.R` and `R/define_effect_fixed_cov.R`

The plan's proposed fix is also backwards: it says fixed effects should be looked up using `component_trait_names` (`plans/fix_add_ebv.md:164-165`). That is wrong for the current schema. Fixed effects are attached to the observed phenotype model, so a maternal phenotype `"WW"` should look up fixed effects where `trait_effects.phenotype_name = 'WW'`, not `"WWD"` or `"WWM"`.

Corrective plan change:

- Change the query to select `phenotype_name`, not `trait_name`.
- Filter fixed effects by observed phenotype names.
- In `write_renum_par()`, per-observation fixed effect presence should be checked against `effects_df$phenotype_name`, not component trait names.

### 2. The BLUPF90 trait count is under-specified and likely wrong for maternal models

The plan says the maternal data file has one observation column (`"WW"`) but the G matrix uses two component traits (`"WWD"`, `"WWM"`). That distinction is correct. The plan does not then specify which vector drives:

- `n_traits`
- `TRAITS`
- per-trait fixed EFFECT columns
- residual matrix dimension
- `trait_num` mapping in `solutions.orig`

In the current `write_renum_par()`, `n_traits <- length(trait)`, `TRAITS` uses `col_map[trait]`, and fixed/random EFFECT columns are repeated over `trait` (`R/blupf90_helpers.R:276-346`). If `trait` is replaced by `component_trait_names`, the parameter file would try to model two observed traits while `data.txt` has one observation column. If `trait` remains `phenotype_name`, the G matrix needs a separate path that does not assume G dimension equals the number of `TRAITS` columns.

Corrective plan change:

- Introduce separate concepts explicitly:
  - `phenotype_names`: observed columns in `data.txt` and names for R.
  - `component_trait_names`: genetic effects in G and names for `ind_ebv`.
  - `phenotype_to_components`: mapping from each observed phenotype to its component traits and contributor roles.
- For the first implementation, strongly consider supporting only one maternal phenotype per BLUPF90 run and erroring on mixed simple+maternal multi-phenotype calls until the parameter-file design is validated.

### 3. Maternal BLUPF90 syntax and effect numbering are not sufficiently specified

The plan says to add an `OPTIONAL mat` line and assume `animal_effect_num + 1` is the maternal effect (`plans/fix_add_ebv.md:181-184`, `186-194`). That is too thin for an implementation plan.

The implementation needs to specify:

- The exact `OPTIONAL mat` line syntax, including the dam ID column number.
- Whether a separate `EFFECT` block for the dam ID column is required or whether `OPTIONAL mat` modifies the existing `RANDOM animal` block.
- Whether `FIELDS_PASSED TO OUTPUT` must include the dam column, and whether that affects `solutions.orig`.
- How `renumf90` numbers direct and maternal effects in `solutions.orig`.
- A synthetic or real expected `renum.par` example.

Without this, `parse_blupf90_solutions()` may map the wrong rows or no rows.

### 4. The plan incorrectly says `phenotype_components` is unused

The plan states that `phenotype_components` is an unused artifact and should not be relied on (`plans/fix_add_ebv.md:12-14`, `246`). That conflicts with the codebase:

- `define_phenotype(..., components = ...)` still documents and writes `phenotype_components` (`R/define_phenotype.R:75-86`, `427-470`).
- `add_phenotype()` actively reads `phenotype_components` (`R/add_phenotype.R:211-226`).
- `.assemble_composite_tbv()` is a full implementation over `phenotype_components` (`R/add_phenotype.R:881-1008`).

Corrective plan change:

- Either support both `formula_tbv` and `phenotype_components` when resolving EBV component traits, or explicitly mark `components`-based composites as unsupported in `add_ebv()` with a clear error.
- Do not describe `phenotype_components` as unused unless a separate migration/removal plan exists.

### 5. Renaming `trait_name` to `phenotype_name` would break named callers

The plan says to rename the public argument (`plans/fix_add_ebv.md:89-92`, `138-150`) and later says simple trait call sites continue to work (`plans/fix_add_ebv.md:247-250`). Positional calls like `add_ebv("ADG")` would work, but named calls like `add_ebv(trait_name = "ADG", ...)` would break.

This repo already documents named `trait_name` usage in `plans/run_blupf90_manual.md`, and users may have scripts doing the same.

Corrective plan change:

- Prefer a backward-compatible signature, for example:
  - keep `trait_name` as a deprecated alias, or
  - add `phenotype_name = trait_name` while still accepting `trait_name`.
- Update docs to say the argument now refers to observed phenotype names, even if the compatibility alias remains.

## Medium Severity Issues

### 6. `formula_tbv` resolution silently drops unsupported contributor types

The plan extracts only `self` and `dam` refs from `.walk_formula_tbv_ast()` (`plans/fix_add_ebv.md:117-125`). The walker can return `sire`, `group_sum`, and `group_mean` too (`R/formula_helpers.R:140-225`).

For a social genetic effect like:

```r
ADG_direct + group_sum(ADG_social, pen_id)
```

the plan would not include `ADG_social` in `component_traits`. That is wrong if parent-average mode is meant to provide "general composite support", and dangerous if BLUPF90 support is narrower than the phenotype DSL.

Corrective plan change:

- Explicitly validate component roles.
- For BLUPF90, support only `self + dam(...)` initially, and error on `sire`, `group_sum`, `group_mean`, scalar weights, nonlinear arithmetic, repeated refs, and anything else not intentionally modeled.
- For parent average, decide whether all refs should become independent component EBVs or whether some roles should error.

### 7. Derived formula support is overclaimed

The plan says derived phenotypes such as FCR should resolve formula symbols with `.extract_all_symbols()` and return parent averages for those components (`plans/fix_add_ebv.md:131-134`, `196-220`). This is plausible for parent-average mode when all formula symbols are simple phenotypes. It is not generally correct.

Problems:

- `formula` symbols are phenotype names, not guaranteed trait names (`R/define_phenotype.R:87-98`).
- A formula symbol can refer to a composite phenotype, which then needs recursive resolution.
- `derived_formula` phenotypes have no independent residual variance by design (`R/define_phenotype.R:261-289`), so BLUPF90 mode cannot just use the derived phenotype without a defined residual/R strategy.
- Math functions appear in formulas; `.extract_all_symbols()` returns names that must be filtered against `.FORMULA_MATH_WHITELIST` (`R/formula_helpers.R:104-135`). The plan's snippet does not mention this filter.

Corrective plan change:

- Scope derived-formula support to `parent_avg = TRUE` unless BLUPF90 semantics are explicitly designed.
- Resolve formula symbols recursively through phenotype metadata.
- Validate that final resolved components correspond to genetic `trait_meta` rows before using them for `ind_ebv`.

### 8. Multi-phenotype and mixed model calls need a clear support boundary

Current `add_ebv()` accepts a vector. The plan does not define behavior for calls such as:

```r
add_ebv(c("WW", "ADG"), software = "blupf90")
```

where `"WW"` expands to `"WWD"` and `"WWM"` but `"ADG"` is simple. R would be keyed by observed phenotypes (`WW`, `ADG`), G by genetic components (`WWD`, `WWM`, `ADG`), and solutions parsing would need a nontrivial mapping across observed trait number and effect number.

Corrective plan change:

- Either support only one maternal/composite phenotype per BLUPF90 call initially, or fully specify the multi-phenotype mapping.
- Keep parent-average mode vectorized; it is simpler because it only reads/writes `ind_ebv`.

### 9. Empty phenotype filters currently fall back to all records

`add_ebv()` warns when a phenotype filter matches no records (`R/add_ebv.R:183-189`), but `build_data_file()` only applies `pheno_ids` when `length(pheno_ids) > 0` (`R/blupf90_helpers.R:112-123`). If the filter matches zero rows, `pheno_clause` becomes empty and all phenotype records are included.

This is not mentioned in the plan, but it is adjacent to the BLUPF90 data-file fix and should be fixed while touching `build_data_file()`.

Corrective plan change:

- If `phenotype` is supplied and matches zero rows, either stop before BLUPF90 mode or pass a sentinel that guarantees zero rows.

### 10. SQL construction should use existing escaping helpers

The plan continues the current style of building SQL with `paste0("'", x, "'")`. Some values are validated identifiers, but IDs and database-sourced values can still contain characters that make ad hoc SQL brittle.

Corrective plan change:

- Use existing helpers such as `sql_in_list()` / `format_sql_value()` where appropriate.
- At minimum, escape all string literals with `gsub("'", "''", x)`.

## Lower Severity / Documentation Issues

### 11. Function/file references are stale

The plan says `ebv_parent_avg()` is in `R/blupf90_helpers.R` (`plans/fix_add_ebv.md:233`), but it is currently in `R/add_ebv.R:295-394`.

### 12. Regenerated docs list is incomplete

The plan only lists `man/add_ebv.Rd` for regeneration (`plans/fix_add_ebv.md:234`). Roxygen docs exist for internal helpers too, and this repo already has man pages for them:

- `man/build_data_file.Rd`
- `man/write_renum_par.Rd`
- `man/write_meta_file.Rd`
- `man/parse_blupf90_solutions.Rd`
- `man/update_covars_from_blupf90.Rd`

Any signature changes in `R/blupf90_helpers.R` or `R/add_ebv.R` should regenerate all affected `.Rd` files.

### 13. Test coverage is missing from the plan

The plan should add tests that do not require BLUPF90 binaries:

- `build_data_file()` uses `phenotype_name` and no longer crashes on `ind_phenotype`.
- Fixed effects are loaded from `trait_effects.phenotype_name`.
- A zero-row phenotype filter does not include all records.
- Component resolver handles simple, formula_tbv maternal, formula_tbv unsupported roles, components-based composites, and derived formulas.
- Parent-average mode returns component rows for maternal/composite phenotypes.
- `write_renum_par()` emits the exact expected maternal `renum.par` for a small synthetic case.
- `parse_blupf90_solutions()` maps synthetic direct and maternal solution rows to the intended component trait names.
- Named `trait_name = ...` remains accepted or gives a deliberate deprecation warning.

## Suggested Revised Scope

For a safe first implementation:

1. Fix schema crashes for `ind_phenotype.phenotype_name` and `trait_effects.phenotype_name`.
2. Rename concepts internally, but keep `trait_name` as a compatibility alias for the public API.
3. Add a resolver that returns a structured mapping:

```r
list(
  phenotype_name = "WW",
  observed_names = "WW",
  components = tibble::tibble(
    trait_name = c("WWD", "WWM"),
    role = c("self", "dam")
  ),
  model_class = "maternal_self_dam"
)
```

4. Implement parent-average composite support first.
5. For BLUPF90 maternal support, restrict to one phenotype with exactly `self` plus `dam` components and no weights/nonlinear expressions until the generated `renum.par` and `solutions.orig` mapping are validated.
6. Add explicit errors for unsupported composite forms instead of silently treating them as simple traits.


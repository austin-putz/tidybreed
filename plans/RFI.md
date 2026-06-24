# RFI (Residual Feed Intake) Simulation — Design Plan

## What RFI Is

RFI is a statistical adjustment, not a genetic component. Its genetic architecture
emerges automatically from the ADFI genetics once you partial out covariate phenotypes.
It is the residual from regressing ADFI on growth and body composition covariates:

```
ADFI_ij = μ + b1·ADG_ij + b2·MBW_ij + b3·backfat_ij + sex_i + pen_j + e_ij

RFI_ij = ê_ij  (the residual)
```

Key properties:
- All covariates are **phenotypic**, not genetic
- No `residual_var` specified by the user — residual variance emerges from the regression
- RFI mean is ~0 by construction
- No `define_trait("RFI")` call — there are no QTL for RFI, only for its components
- By construction, RFI is genetically uncorrelated with the covariate traits when
  the regression is done at the phenotypic level

---

## Workflow: How Users Simulate RFI

1. `define_trait()` + `define_additive_effects()` for ADFI, ADG, MBW, backfat
2. `define_effect_cov_matrix("gen_add", G)` — with realistic genetic covariance
   matrix among ADFI and all covariates
3. `define_phenotype()` for ADFI, ADG, MBW, backfat (all standard continuous traits)
4. `define_phenotype("RFI", trait_type = "regression_residual", ...)` — RFI definition
5. `add_phenotype()` — simulate ADFI and covariate phenotypes first, then RFI

---

## Proposed API: `define_phenotype()` Extensions

```r
define_phenotype(
  pop,
  phenotype_name        = "RFI",
  trait_type            = "regression_residual",  # new type
  regression_response   = "ADFI",                 # response phenotype (must be in phenotype_meta)
  regression_covariates = list(                   # named list: name → effect type
    ADG     = "fixed_cov",
    MBW     = "fixed_cov",
    backfat = "fixed_cov",
    sex     = "fixed_class"   # pulled from ind_meta, not a phenotype
  ),
  regression_random     = c("pen"),   # character vector → triggers lmer; NULL → lm()
  store_betas           = FALSE,      # TRUE: fit once, reuse across generations
  expressed_sex         = "both"
  # no mean, no residual_var — both emerge from the regression
)
```

Covariate source resolution (automatic):
- Name found in `phenotype_meta` → pull from `ind_phenotype`
- Name found as column in `ind_meta` → pull from `ind_meta`
- Otherwise → error

---

## What `add_phenotype("RFI")` Does Internally

1. Verify all response and covariate phenotypes have records for the individual subset
2. Collect response + covariates into a flat data frame (one row per individual)
3. Build formula dynamically from stored `phenotype_regression_spec` metadata
4. If `regression_random` is NULL → `lm()`; if non-empty → `lme4::lmer()`
5. Extract residuals → write to `ind_phenotype` as `pheno_value`
6. If `store_betas = TRUE` → write coefficients to `phenotype_regression_betas`;
   subsequent calls use `predict()` with stored betas instead of refitting

---

## The `store_betas` Question (Important for Multi-Generation Simulations)

Without stored betas, each generation's RFI is computed from a regression fit on
that generation's phenotypes. Betas shift across generations as the population
changes genetically — **RFI is not comparable across generations**.

With `store_betas = TRUE`, the model is fit once (typically on the base population),
coefficients are stored, and applied to all future generations. RFI is then
comparable across generations — which is what most research protocols do.

- Default `store_betas = FALSE`: fine for within-generation selection
- `store_betas = TRUE`: required for multi-generation breeding program simulations

---

## New Storage Tables

### `phenotype_regression_spec`

Stores covariate specification for each regression_residual phenotype.

| Column             | Type    | Notes                                              |
|--------------------|---------|----------------------------------------------------|
| id_reg_spec        | INTEGER | PK (auto-incrementing)                             |
| phenotype_name     | VARCHAR | FK to `phenotype_meta`                             |
| response_phenotype | VARCHAR | e.g. "ADFI"                                        |
| covariate_name     | VARCHAR | e.g. "ADG", "sex", "pen"                           |
| covariate_source   | VARCHAR | "phenotype" or "ind_meta"                          |
| covariate_type     | VARCHAR | "fixed_cov", "fixed_class", "random"               |

### `phenotype_regression_betas`

Stores fitted regression coefficients when `store_betas = TRUE`.

| Column         | Type      | Notes                                              |
|----------------|-----------|----------------------------------------------------|
| id_reg_beta    | INTEGER   | PK (auto-incrementing)                             |
| phenotype_name | VARCHAR   | FK to `phenotype_meta` ("RFI")                     |
| term           | VARCHAR   | e.g. "ADG", "MBW", "(Intercept)"                   |
| estimate       | DOUBLE    | Regression coefficient                             |
| std_error      | DOUBLE    | SE from model                                      |
| fit_n          | INTEGER   | Sample size used for fitting                       |
| fit_timestamp  | TIMESTAMP | When betas were stored                             |

---

## Package Dependency: `lme4`

`lme4` belongs in `Suggests`, not `Imports` — it is heavy and most users will not
need it. Use a runtime check:

```r
if (!requireNamespace("lme4", quietly = TRUE)) {
  stop(
    "Package 'lme4' is required for regression_random effects in RFI. ",
    "Install it with: install.packages('lme4')",
    call. = FALSE
  )
}
```

Only triggered when `regression_random` is non-empty.

---

## Summary of Code Changes Needed

| What                                               | Where                       |
|----------------------------------------------------|-----------------------------|
| New `trait_type = "regression_residual"` value     | `phenotype_meta` schema      |
| New args: `regression_response`, `regression_covariates`, `regression_random`, `store_betas` | `define_phenotype()` |
| New table `phenotype_regression_spec`              | `initialize_genome()`        |
| New table `phenotype_regression_betas`             | `initialize_genome()`        |
| New dispatch branch for regression_residual type  | `R/add_phenotype.R`          |
| Runtime `lme4` check                              | Inside the new dispatch branch |

Complexity: roughly equivalent to the SGE group contributor implementation.

---

## Open Questions

### 1. Within-group vs. whole-subset regression

In real data analysis, the regression is often fit within contemporary group
(pen, batch, test period). Two options:

- **Option A (default)**: fit across the whole individual subset passed to `add_phenotype()`
- **Option B**: fit separately within levels of a `within_group` column (e.g. pen)

Suggestion: default to Option A (whole subset) with an optional `within_group`
argument for Option B. Keeps the simple case simple.

### 2. Order of `add_phenotype()` calls

When the user calls `add_phenotype()` without a `phenotype_name` (i.e. simulate
all phenotypes), `add_phenotype()` must detect that "RFI" depends on "ADFI" and
covariate phenotypes being simulated first and order the dispatch accordingly.
Currently dispatch order is not defined. This needs a dependency-resolution step.

### 3. Comparable across generations: how to expose stored betas to the user

When `store_betas = TRUE`, should the user be able to inspect, reset, or override
the stored betas? For example:
- `get_table("phenotype_regression_betas")` — already works generically
- An explicit `reset_regression_betas(pop, "RFI")` helper?
- Allow passing `betas = <data frame>` directly to `add_phenotype()` to override?

### 4. What happens when a covariate record is missing for an individual?

If an individual has ADFI but no ADG record, should they:
- Be excluded from the regression and receive `NA` for RFI (`"skip"`)
- Cause an error (`"error"`)

Suggestion: reuse the existing `missing_component_action` field already in
`phenotype_meta` — the same skip/error logic applies here.

### 5. lmer random effect grouping column source

For a random effect like `pen`, where does the group column come from? Currently
`ind_meta` is the natural place, but users may have pen stored elsewhere.
Mirror the `source_table` pattern from `trait_effects` — add a
`covariate_table` column to `phenotype_regression_spec`.

### 6. Should RFI TBV be computable?

It could be useful to compute a "true RFI" — the expected residual under the
true (known) regression coefficients from the simulation. This would require
storing the population-level phenotypic regression slope as a known quantity.
Probably out of scope for the first implementation but worth noting.

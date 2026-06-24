# Formula-Based Phenotype Specification — Design Plan

Covers a unified formula string API for two distinct phenotype assembly problems:
(1) how additive TBV components compose into a composite genetic value, and
(2) how phenotype records combine arithmetically into a derived value. These are
different operations that reference different tables and require different parsing
machinery.

See `RATIO_TRAITS.md` for the motivating FCR case (type 2 only).

---

## The Three-Layer Phenotype Model

Every phenotype observation decomposes into three layers:

```
P  =  TBV composition  +  fixed/random effects  +  residual
P  =  (WWD + 0.5·WWM)  +  (sex + litter + pe)   +  e
```

Formula strings are appropriate for **layer 1** (TBV composition) and for
**post-hoc arithmetic derivation** over phenotype records such as ratio
traits (e.g. feed conversion ratio). They are not appropriate for layers 2 and 3:

- Random effects and the residual require **sampling from distributions** — they
  are not arithmetic over existing stored values. `pe_dam` is a draw from
  `N(0, σ²_pe)` keyed to a dam ID, not a column you compute row-wise.
- Multi-trait residual and random-effect correlations are **matrix quantities
  spanning phenotypes** — they cannot be expressed inside a single-phenotype
  formula. The existing `define_effect_cov_matrix()` / `phenotype_residual_cov`
  infrastructure handles this correctly and is not replaced.

Layers 2 and 3 remain in the existing `define_effect_*` / `define_residual_cov()`
infrastructure.

---

## Two Formula Types

### `formula_tbv` — TBV composition layer

Specifies how additive TBV components combine to produce the genetic portion of
the phenotype. Symbols are **trait names from `trait_meta`**. Values come from
`ind_tbv`. Contributor type (self / dam / sire / group) is expressed through a
set of reserved DSL functions.

The name `formula_tbv` is intentionally specific: it reads from `ind_tbv` (additive
true breeding values). Dominance deviations and epistatic interactions are separate
quantities not currently stored in `ind_tbv` and are not covered here. No `tbv()`
wrapper is needed inside the expression — bare symbols already imply a TBV lookup;
the contributor functions specify *whose* TBV to fetch.

Used for: maternal models (WW), social genetic effect models (ADG), any
composite phenotype where the genetic architecture spans multiple `trait_meta`
entries.

`formula_tbv` and the `components` data frame are mutually exclusive in a single
`define_phenotype()` call. When `formula_tbv` is set, `phenotype_components`
has **no rows** for that phenotype — the formula is the sole specification and
`add_phenotype()` evaluates it directly without touching `phenotype_components`.
The `components` data frame remains available for advanced cases that the formula
syntax cannot express (per-individual covariate weights, Legendre polynomial
weights — see "Relationship to the `components` Data Frame" below).

### `formula` — derived arithmetic layer

Specifies post-hoc arithmetic over **existing phenotype records**. Symbols are
**phenotype names from `phenotype_meta`**. Values come from `ind_phenotype`.
Evaluated row-wise using R's `eval(parse())` after collecting component values
into a wide data frame.

Used for: ratio traits (FCR), scaled traits, power traits (Kleiber), RFI with
user-supplied regression coefficients, any derived quantity computed from
already-simulated phenotypes.

See `RATIO_TRAITS.md` for full detail on this type.

---

## `formula_tbv` — DSL Grammar

The formula is a standard R arithmetic expression. Four reserved DSL functions
extend it with contributor-aware TBV lookup. All other symbols are treated as
self-TBV trait names.

| Syntax                              | Meaning                                              |
|-------------------------------------|------------------------------------------------------|
| `ADG_direct`                        | Self TBV for trait `"ADG_direct"`                    |
| `self(ADG_direct)`                  | Explicit self (identical to bare symbol)             |
| `dam(WWM)`                          | Dam's TBV for trait `"WWM"` (joined via `id_parent_2`) |
| `sire(WWS)`                         | Sire's TBV for trait `"WWS"` (joined via `id_parent_1`) |
| `group_sum(ADG_SGE, pen_id)`        | Sum of group-mates' TBVs for `"ADG_SGE"`, grouped by `pen_id` in `ind_meta` (self excluded) |
| `group_mean(ADG_SGE, pen_id)`       | Mean of group-mates' TBVs (self excluded)            |
| `0.5`, `100`, `0.036`               | Scalar numeric literals — no lookup, used as weights |

Standard R arithmetic (`+`, `-`, `*`, `/`, `^`, parentheses) composes all
terms. All DSL functions are pre-evaluated to per-individual scalar vectors
before `eval()` is called on the remaining expression.

### Examples

```r
# Maternal model — weaning weight
define_phenotype(pop, "WW", trait_type = "composite",
                 formula_tbv = "WWD + 0.5 * dam(WWM)")

# SGE / social genetic effects — average daily gain
define_phenotype(pop, "ADG_obs", trait_type = "composite",
                 formula_tbv = "ADG_direct + group_sum(ADG_SGE, pen_id)")

# SGE with mean aggregation
define_phenotype(pop, "ADG_obs", trait_type = "composite",
                 formula_tbv = "ADG_direct + group_mean(ADG_SGE, pen_id)")

# Sire model — a trait expressed in the progeny but partly from sire
define_phenotype(pop, "SC", trait_type = "composite",
                 formula_tbv = "SC_direct + 0.25 * sire(SC_sire)")

# Two separate group effects — not yet valid (see Open Question #5)
# formula_tbv = "ADG_direct + group_sum(SGE_pen, pen_id) + group_sum(SGE_barn, barn_id)"
```

---

## How `formula_tbv` Is Parsed and Evaluated

### At `define_phenotype()` time

1. Parse the formula string: `parse(text = formula_tbv)[[1]]`
2. Walk the AST. Detect:
   - Bare symbols → treat as self-TBV trait names; validate against `trait_meta`
   - `self(x)` calls → same as bare symbol
   - `dam(x)` / `sire(x)` calls → validate `x` against `trait_meta`; note contributor
   - `group_sum(x, col)` / `group_mean(x, col)` calls → validate `x` against `trait_meta`;
     note `col` as a group column (validated against `ind_meta` columns at call time or
     deferred — see Open Question #2)
3. Validate no unknown symbols remain (after stripping DSL function names and numeric
   literals)
4. Store the raw formula string in `phenotype_meta.formula_tbv`

### At `add_phenotype()` time

1. Re-parse formula → extract all DSL function calls and bare symbols
2. For each term, pre-fetch the required TBV vector aligned to the individual subset:
   - Bare symbol / `self(x)` → fetch `ind_tbv` rows for this individual subset
   - `dam(x)` → join individual subset to `ind_meta` on `id_parent_2`, fetch dam TBVs
   - `sire(x)` → join on `id_parent_1`, fetch sire TBVs
   - `group_sum(x, col)` → join individual subset to `ind_meta` to get group memberships;
     fetch all group-member TBVs; sum per individual (excluding self); return aligned vector
   - `group_mean(x, col)` → same but mean
3. Construct an environment mapping each term's placeholder name to its pre-fetched vector
4. Replace each DSL call in the AST with its placeholder symbol
5. `eval(modified_expr, envir = pre_fetched_env)` — produces the composite TBV vector
6. This vector replaces what `.assemble_composite_tbv()` currently computes from the
   `phenotype_components` rows
7. Downstream assembly (fixed effects, residual, liability conversion) proceeds unchanged

---

## `formula` — Derived Arithmetic (FCR Type)

See `RATIO_TRAITS.md` for full specification. Summary:

- Symbols are phenotype names from `phenotype_meta`; values from `ind_phenotype`
- Standard R arithmetic + whitelisted math functions (`sqrt`, `log`, `exp`, `abs`, etc.)
- No contributor functions (`dam()`, `group_sum()`, etc.) — those are `formula_tbv` only
- Evaluated after all component phenotypes have been simulated
- No `ind_tbv` rows written; no genetic architecture implied

```r
define_phenotype(pop, "FCR",     trait_type = "derived_formula", formula = "ADFI / ADG")
define_phenotype(pop, "Kleiber", trait_type = "derived_formula", formula = "ADG / MBW^0.75")
define_phenotype(pop, "RFI",     trait_type = "derived_formula", formula = "ADFI - 0.036*ADG - 0.0072*MBW")
```

---

## Schema Changes

Two new columns added to `phenotype_meta`:

| Column             | Type    | Notes                                                                |
|--------------------|---------|----------------------------------------------------------------------|
| `formula_tbv`      | VARCHAR | NULL unless `trait_type = "composite"` with formula shorthand        |
| `formula`          | VARCHAR | NULL unless `trait_type = "derived_formula"`                         |

`phenotype_arithmetic_spec` (previously proposed in `RATIO_TRAITS.md`) is **not
created** — both formula types are stored directly in `phenotype_meta`.

When `formula_tbv` is set, `phenotype_components` has **no rows** for that
phenotype — the formula is the sole specification. `phenotype_components` is
only populated when the `components` data frame is passed explicitly for
advanced weighting cases.

**On non-additive genetics**: `formula_tbv` is deliberately scoped to additive
TBVs from `ind_tbv`. If dominance deviations or epistatic interactions are ever
added to the package, they would live in a separate table (e.g. `ind_dom`) and
require a separate formula column (e.g. `formula_dom`). The naming makes this
extension path explicit rather than hiding it behind a generic `formula_genetic`.

---

## Relationship to the `components` Data Frame

`formula_tbv` is the **primary** way to define composite phenotypes. It is not
syntactic sugar that compiles down to `phenotype_components` — it is a separate
evaluation path that `add_phenotype()` executes directly from the stored formula
string. `phenotype_components` has no rows for formula-based phenotypes.

`formula_tbv` cannot express:

- Per-individual covariate weights (`weight_type = "covariate"`)
- Legendre polynomial weights (`weight_type = "legendre"`)
- Group columns from tables other than `ind_meta` (see Open Question #4)

For those rare advanced cases the `components` data frame remains available.
`formula_tbv` and `components` in the same `define_phenotype()` call → error.

---

## Dispatch Ordering for `add_phenotype()`

When `add_phenotype()` is called with no `phenotype_name`, all phenotypes in
`phenotype_meta` are dispatched. The required ordering is:

1. Composite phenotypes with `formula_tbv` (or `components`) — after their
   source traits have TBVs in `ind_tbv`
2. Derived arithmetic phenotypes with `formula` — after all referenced phenotype
   records are in `ind_phenotype`

A dependency graph is built by parsing all `formula_tbv` and `formula`
strings, then topologically sorted. Cycles are detected and raise an error.

---

## Open Questions

### 1. `dam()` / `sire()` vs. `_dam` / `_sire` suffix for contributor syntax

Two syntaxes are possible for single-parent contributors:

- **Function syntax**: `dam(WWM)`, `sire(WWS)` — unambiguous, R-parseable,
  requires AST inspection
- **Suffix syntax**: `WWM_dam`, `WWS_sire` — more compact, but conflicts with
  trait names that genuinely end in `_dam` or `_sire`, and requires a naming
  convention rule across all of `trait_meta`

Recommendation is function syntax (`dam()`, `sire()`) but this needs a decision.

**DECISION:** Go with recommended of `dam()` and `sire()` etc to be explicit, much better

### 2. When to validate group column names (`pen_id`, `barn_id`)

`group_sum(ADG_SGE, pen_id)` specifies `pen_id` as a column in `ind_meta`. Should
this column's existence be validated at `define_phenotype()` time (requires querying
`ind_meta` schema immediately) or deferred to `add_phenotype()` time (fails later
but avoids a schema dependency at definition time)?

Option: validate at `add_phenotype()` time only, with a clear error message
referencing the formula.

**DECISION:** I think validate only at `add_phenotype()` to ERROR, but provide a
clear warning to the user at `define_phenotype()` as well if possible please. 

### 3. Does `phenotype_components` still need to exist?

**Decision**: `phenotype_components` is NOT used for `formula_tbv` phenotypes.
`add_phenotype()` reads `formula_tbv` from `phenotype_meta` directly and
evaluates it — no rows are written to or read from `phenotype_components` for
those phenotypes.

`phenotype_components` is retained solely for the advanced cases that
`formula_tbv` cannot express: per-individual covariate weights
(`weight_type = "covariate"`) and Legendre polynomial weights
(`weight_type = "legendre"`). For most users these cases will never arise and
`phenotype_components` becomes an internal detail of the advanced path only.

Open sub-question: should `phenotype_components` eventually be deprecated in
favour of extending `formula_tbv` syntax (e.g. `group_sum(SGE, pen_id, weight = days_on_test)`)?
Or kept as the permanent escape hatch for things the formula DSL cannot cover?

**DECISION:** I think we leave it available for now but most not used unless those
edge cases need it later, we've not implemented random regression models yet
so let's wait to remove until then please. 


### 4. Group table other than `ind_meta`

Current syntax `group_sum(ADG_SGE, pen_id)` implicitly assumes `pen_id` is in
`ind_meta`. If the group column is in another table (e.g. a housing table), the
syntax needs extension. Options:

- `group_sum(ADG_SGE, pen_id, table = "housing_meta")` — named argument
- Require group columns to always be in `ind_meta` (users must join them there first)
- `group_sum(ADG_SGE, housing_meta.pen_id)` — dot-qualified table reference

Requiring the column to be in `ind_meta` (with the user pre-joining if needed) is
the simplest rule and consistent with how all other effect lookups work.

**DECISION:** Let's use the table argument, but default is "ind_meta" and not
required. Provide ERROR if not found and give user to state the argument name `table=`
for users to learn. 


### 5. Multiple group aggregation terms in one formula

Is `"ADG_direct + group_sum(SGE_pen, pen_id) + group_sum(SGE_barn, barn_id)"` valid
(two different group effects, two different group columns)? This requires two
separate group aggregation passes at assembly time. Should it be supported from
the start, or restricted to one group term per formula?

**DECISION:** Yes we need to allow the user to specify what they want so yes 
allow please. 


### 6. `formula_tbv` with scalar weights: interaction with genetic variance

When the user writes `"WWD + 0.5 * dam(WWM)"`, the `0.5` halves the contribution
of the maternal TBV to the composite genetic value. This is biologically meaningful
(the dam contributes half its additive genetic value). Should `define_phenotype()`
document or store this weight for reference (e.g. in a derived-variance table)?
Or is it purely a runtime detail with no stored metadata?

**DECISION:** If it doesn't impact computational time, I don't think we need to store it
outside the formula, I suspect the evaluation of the formula to extract the 0.5
is not meaningful in terms of added time. 

### 7. Mixing `formula_tbv` and scalar fixed-effect adjustments

Could a user write `"WWD + dam(WWM) - 10"` — subtracting a scalar 10 from the
composite TBV? Syntactically valid and arithmetically unambiguous, but semantically
odd (why adjust TBVs by a constant?). Should scalar arithmetic on the whole
expression be blocked, allowed silently, or allowed with a warning?

**DECISION:** Allowed with a warning, but allowed for sure, we want `tidybreed`
to be maximally flexible even if it doesn't make sense to us and assume
the user knows what he/she is doing. 

### 8. Dependency ordering across both formula types

A phenotype could have `formula_tbv` (composite TBV) and then another
phenotype uses it in `formula` (derived arithmetic). The topological sort must
handle both formula types in the same graph. Is there any case where a
`formula_tbv` phenotype could reference a `formula` phenotype's output?
(Probably not, but worth confirming the ordering rule is one-directional:
`formula_tbv` → phenotype records → `formula`.)

**DECISION:** Good question but I can't think of an example where this
could happen yet. So don't worry about it but do consider order for
implementation just in case it comes up in testing. 

### 9. Error messaging for formula parse failures

If the user writes `"WWD + dam(WM)"` and `"WM"` is not in `trait_meta`, the
error should name the unrecognised symbol and the formula string. Should
validation also suggest close matches (e.g. "Did you mean `WWM`?")?

**DECISION:** Yes, please throw a meaningful error message for users and alert them it needs
to be found in `trait_meta` and defined above in the script. You can try to find
a similar name but not sure how reliable this is, but you can try, certainly 
small typos from capitalization or spelling is important, but don't need to 
go overboard (not certain how easy that is to do). 


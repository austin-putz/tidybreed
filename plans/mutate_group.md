# Design Plan: `mutate_group_*()` Functions

## Overview

Three generic functions for assigning group membership to individuals in any
tidybreed table. All functions accept a `tidybreed_table` from `get_table()`
(optionally filtered) as their first argument and return `tidybreed_pop`.
Group assignment is always random within the filtered set unless otherwise noted.
Unfiltered rows receive `NULL` for new columns or retain existing values.

These functions only write table columns — they do not register model effects
or group semantics. Connecting a group column to `define_phenotype()` components
or `define_effect_random()` is the user's responsibility after assignment.

## Shared Conventions

- **First argument**: `tidybreed_table` (from `get_table()` + optional `filter()`)
- **`col_name`**: name of the new (or existing) column being created or appended to.
  First call creates the column; subsequent calls on new filtered rows add values
  to existing `NULL` rows in that column.
- **`overwrite = FALSE`** (all three functions): when `FALSE` (default), any
  **filtered row** that already has a non-`NULL` value in `col_name` raises an
  error before any sampling occurs. This check is on the rows being written only —
  the same label values appearing in unfiltered rows elsewhere in the column are
  not a collision and do not require `overwrite = TRUE`. Set `TRUE` to replace
  existing values in the filtered rows.
- **Row ordering**: filtered primary keys are collected and sorted by the target
  table's registered primary key before random shuffling. This makes assignment
  stable regardless of row insertion order. No upstream `arrange()` is applied —
  always PK order in v1. An explicit `order_by` argument is deferred to a future
  release.
- **Empty filter**: when the filter matches zero rows the function warns and
  returns the population unchanged, matching `mutate_table()`. A new column is
  NOT created from a zero-row filter — a typo in a filter should not silently
  alter schema.
- **Column type enforcement**: if `col_name` already exists, the function checks
  that its DuckDB type is compatible before writing. `mutate_group_seq()` errors
  on an existing non-INTEGER column; `mutate_group_named()` and
  `mutate_group_concatenate()` error on an existing non-VARCHAR column.
- **Validate before sampling**: all argument validation, overwrite checks, and
  column-type checks happen before any call to `sample()` or `set.seed()`, so
  that a failed call does not advance the caller's RNG state.
- **Seed handling**: when `seed` is supplied the implementation uses
  `withr::with_seed(seed, ...)` for local RNG-state handling so the caller's
  global RNG stream is not permanently reset.
- **Diagnostics**: all three functions message the number of rows assigned, rows
  left `NULL`, and group size summary (min/mean/max). Pass `quiet = TRUE` to
  suppress all messages; default is `quiet = FALSE`.
- **Primary key resolution**: the keyed tibble path uses the target table's
  registered primary key via `TABLE_PRIMARY_KEYS[[table_name]]`, not a
  hard-coded `id_ind`. Functions error early if the table has no registered
  primary key.
- All functions return `tidybreed_pop` invisibly.
- All functions use `mutate_table()` internally to write results; internal
  `mutate_table()` messages are suppressed so only group-level diagnostics print.

### Modeling note

These functions work on any tidybreed table, but group columns have different
meanings depending on where they live:

- `ind_meta`: individual-level group (pen, farm, barn, CG) — the right table
  for `define_effect_random(source_column = "pen")` and phenotype components.
- `ind_phenotype`: observation-level group — not equivalent to an individual
  group; use with care if the column will feed model paths.
- `genome_meta`: locus-level group — valid for locus partitions; not useful
  for animal-level social or environmental effects.

---

## 1. `mutate_group_seq()`

Assigns individuals to sequentially numbered integer groups. Each new group
gets the next integer after the current `MAX(col_name)` in the table by default,
enabling multiple calls on different filtered subsets to build up a complete
column without overlap. The starting number can be overridden with `start`.

### seq — Use cases

- Pen assignment (add pigs to pens 1, 2, 3 … across generations)
- Test group numbering
- Restarting pen numbers per generation (`start = 1L` with a new filter)
- Any unnamed, auto-incrementing group structure

### seq — Signature

```r
mutate_group_seq(
  tbl,
  col_name,
  n_per_group    = NULL,  # scalar or length-2 vector c(min, max) for variable sizes
  n_groups       = NULL,  # alternative: split filtered set into N equal groups
  start          = NULL,  # starting group number; NULL = MAX(col_name) + 1
  include_leftover = FALSE, # TRUE = remainder form a smaller final group
                            # FALSE = remainder animals get NULL (DEFAULT)
  seed           = NULL,  # integer RNG seed; NULL = caller's RNG state
  overwrite      = FALSE, # FALSE (default) = error if filtered rows already have non-NULL values
  quiet          = FALSE  # TRUE = suppress diagnostic messages
)
```

### seq — Arguments

| Argument | Type | Description |
|---|---|---|
| `tbl` | `tidybreed_table` | From `get_table()` + optional `filter()` |
| `col_name` | character | Column name to create or append to; always INTEGER |
| `n_per_group` | integer scalar or `c(min, max)` | Target group size. Scalar = fixed; vector = sample uniformly from `min:max` for each group |
| `n_groups` | integer scalar | Alternative to `n_per_group`; splits filtered set into N groups |
| `start` | integer or NULL | Starting group number. `NULL` = `MAX(col_name) + 1` or `1` if column is new/all NULL. Label reuse in unaffected rows does not require `overwrite = TRUE`. |
| `include_leftover` | logical | `FALSE` (default): remainder animals receive `NULL` with a warning. `TRUE`: remainder form a smaller final group. |
| `seed` | integer or NULL | RNG seed; `NULL` uses caller's RNG state. |
| `overwrite` | logical | `FALSE` (default): error if any filtered row already has non-`NULL` value. `TRUE`: overwrite. |
| `quiet` | logical | `FALSE` (default): print assignment diagnostics. `TRUE`: suppress. |

### seq — Constraints

- `n_per_group` and `n_groups` are mutually exclusive; exactly one must be provided.
- Column type is always `INTEGER`. No prefix/VARCHAR option.
- Edge cases that error: `n_groups > n` (more groups than animals), `n_per_group <= 0`,
  `n_groups <= 0`, variable bounds where `min > max` or either is non-positive.

### seq — Group size precision

**`n_groups` mode:**
- `include_leftover = FALSE`: group size = `floor(n / n_groups)`; the remaining
  `n %% n_groups` animals (chosen after shuffling) receive `NULL` and a warning
  prints the count.
- `include_leftover = TRUE`: the first `n %% n_groups` groups (after shuffling)
  receive one extra animal so total allocation = n. Groups are as balanced as
  possible — no single "dump" group.

**`n_per_group = c(min, max)` mode:** group sizes are sampled independently from
`sample(min:max, 1)` until animals run out. If the cumulative size of sampled
groups exceeds n before all groups are filled, the last group is truncated to
the remaining animals. `include_leftover` applies only to the final partial
group once the total exceeds n.

### seq — Examples

```r
# Fixed pen size of 10; continues from MAX pen number on second call
pop |> get_table("ind_meta") |> filter(gen == 1L, sex == "M") |>
  mutate_group_seq(col_name = "pen", n_per_group = 10L)
pop |> get_table("ind_meta") |> filter(gen == 1L, sex == "F") |>
  mutate_group_seq(col_name = "pen", n_per_group = 10L)

# Variable pen size 10-15, reproducible
pop |> get_table("ind_meta") |> filter(gen == 2L, sex == "M") |>
  mutate_group_seq(col_name = "pen", n_per_group = c(10L, 15L), seed = 42L)

# Split into 5 equal groups; remainder excluded
pop |> get_table("ind_meta") |> filter(gen == 2L) |>
  mutate_group_seq(col_name = "test_group", n_groups = 5L, include_leftover = FALSE)

# Restart pen numbering for each generation (pen 1-N per gen, not globally unique)
pop |> get_table("ind_meta") |> filter(gen == 2L) |>
  mutate_group_seq(col_name = "pen", n_per_group = 10L, start = 1L)
# This is valid even though pen 1-N already exist for gen == 1L rows —
# start = 1L is label reuse, not row overwrite; overwrite = TRUE not required.
```

---

## 2. `mutate_group_named()`

Assigns individuals to a fixed set of named groups. Groups are always VARCHAR.
Allocation is controlled by `counts` (exact animals per group) or `proportions`
(probability-based). Typical use: farm splits, barn assignments, sex labeling,
gate cuts.

### named — Use cases

- Farm or barn assignment ("farm_A", "farm_B")
- Sexed semen allocation ("intact", "castrate")
- Sale group splits ("keep", "sell")
- Any scenario where the group names are known in advance

### named — Signature

```r
mutate_group_named(
  tbl,
  col_name,
  group_names,                  # character vector of group labels
  counts          = NULL,       # integer vector, same length as group_names
  proportions     = NULL,       # numeric vector summing to ~1, same length as group_names
  prop_method     = "exact",    # proportions sampling (only when proportions used):
                                #   "exact"  = sample without replacement from pre-built vector (DEFAULT)
                                #   "random" = sample with replacement, each animal independently
  include_leftover = FALSE,     # TRUE = remaining animals added to last group
                                # FALSE = remaining animals get NULL (DEFAULT)
  underfill_action = "error",   # what to do when n < sum(counts):
                                #   "error"        = stop with clear message (DEFAULT)
                                #   "proportional" = rescale counts to available n
                                #   "sequential"   = fill groups in order until run out
  seed            = NULL,       # integer RNG seed; NULL = caller's RNG state
  overwrite       = FALSE,      # FALSE = error if filtered rows already have non-NULL values
  quiet           = FALSE       # TRUE = suppress diagnostic messages
)
```

`include_leftover` (n > sum(counts)): when `FALSE` (default) extras receive
`NULL` and a warning is printed stating the count and suggesting
`include_leftover = TRUE`. When `TRUE` extras are added to the last group.

`underfill_action` (n < sum(counts)): when fewer animals are available than
counts require (e.g., 15 animals but `counts = c(10, 10)` requires 20):

- `"error"` (default): stop before sampling with a message stating available vs.
  required count and suggesting `underfill_action = "proportional"`.
- `"proportional"`: rescale counts preserving the ratio. Uses largest-remainder
  rounding; tie-breaker is stable left-to-right order in `group_names`.
- `"sequential"`: fill groups in order until animals run out. Later groups get
  fewer or zero animals — biases later groups, use with caution.

`prop_method` (applies to `proportions` only):

- `"exact"` (default): builds a pre-sized label vector from scaled counts,
  shuffles, assigns. Group counts are deterministic to the proportion. Rounding
  uses largest-remainder; tie-breaker is stable left-to-right order in
  `group_names`. This means for `n = 5`, `proportions = c(1/3, 1/3, 1/3)` the
  first two groups in `group_names` get 2 animals and the third gets 1.
- `"random"`: each animal drawn independently via `sample(..., replace = TRUE,
  prob = proportions)`. Converges to target proportion over many replicates but
  any single run varies. Use for stochastic processes (natural sex ratios, etc.).

`proportions` tolerance: accepted within `sqrt(.Machine$double.eps)` of 1.0 and
normalized automatically. Negative, `NA`, `Inf`, or all-zero values error.

### named — Arguments

| Argument | Type | Description |
|---|---|---|
| `tbl` | `tidybreed_table` | From `get_table()` + optional `filter()` |
| `col_name` | character | Column name to create or append to; always VARCHAR |
| `group_names` | character vector | Names of the groups in assignment order |
| `counts` | integer vector | Exact number of animals per group; length must equal `length(group_names)` |
| `proportions` | numeric vector | Probability per group; must sum to ~1; length must equal `length(group_names)` |
| `prop_method` | character | Proportions strategy: `"exact"` (default) or `"random"`. Ignored when `counts` is used. |
| `include_leftover` | logical | When n > sum(counts): `FALSE` (default) = extras get `NULL` with warning. `TRUE` = extras added to last group. |
| `underfill_action` | character | When n < sum(counts): `"error"` (default), `"proportional"`, or `"sequential"`. |
| `seed` | integer or NULL | RNG seed; `NULL` uses caller's RNG state. |
| `overwrite` | logical | `FALSE` (default): error if any filtered row already has non-`NULL` value. `TRUE`: overwrite. |
| `quiet` | logical | `FALSE` (default): print assignment diagnostics. `TRUE`: suppress. |

### named — Constraints

- `counts` and `proportions` are mutually exclusive; exactly one must be provided.
- `include_leftover` and `underfill_action` only apply with `counts`; with
  `proportions` all filtered animals are always assigned.
- `prop_method` only applies with `proportions`; ignored when `counts` is used.
- No auto-continuation; names are static.

### named — Examples

```r
# Equal proportions, exact allocation; c(1/3, 1/3, 1/3) accepted via tolerance
pop |> get_table("ind_meta") |> filter(gen == 1L) |>
  mutate_group_named(
    col_name    = "replicate",
    group_names = c("A", "B", "C"),
    proportions = c(1/3, 1/3, 1/3)
  )

# Stochastic proportions, reproducible
pop |> get_table("ind_meta") |> filter(gen == 1L) |>
  mutate_group_named(
    col_name    = "farm",
    group_names = c("farm_A", "farm_B"),
    proportions = c(0.5, 0.5),
    prop_method = "random",
    seed        = 99L
  )

# Fixed capacity: 600 sows to each farm; leftover excluded, error if underfilled
pop |> get_table("ind_meta") |> filter(sex == "F", gen == 2L) |>
  mutate_group_named(
    col_name         = "farm",
    group_names      = c("farm_A", "farm_B"),
    counts           = c(600L, 600L),
    include_leftover = FALSE,
    underfill_action = "error"
  )

# Sexed semen: 75% intact, 25% castrate
pop |> get_table("ind_meta") |> filter(sex == "M") |>
  mutate_group_named(
    col_name    = "castration_group",
    group_names = c("intact", "castrate"),
    proportions = c(0.75, 0.25)
  )
```

---

## 3. `mutate_group_concatenate()`

Creates a new VARCHAR column by concatenating values from one or more existing
columns in the **same table**. Deterministic — no random assignment. Every row
in the filtered set receives a value defined by its own covariate combination.

Cross-table concatenation is not supported. To combine columns from different
tables, use `collect()` + join in R, update the target table with
`mutate_table()` (keyed tibble path), then call `mutate_group_concatenate()`.

### concat — Use cases

- Contemporary group IDs (`sex + barn + generation`)
- Full-sib group labels (`id_parent_1 + id_parent_2`)
- Paternal half-sib labels (`id_parent_1` only)
- Maternal half-sib labels (`id_parent_2` only)
- Any composite identifier from existing columns in the same table

### concat — Signature

```r
mutate_group_concatenate(
  tbl,
  col_name,
  source_columns,             # character vector of column names to concatenate
  sep         = "_",          # separator between values
  null_action = "propagate",  # "propagate" = any NULL source → NULL result (DEFAULT)
                              # "skip"      = NULL sources omitted (CONCAT_WS behavior)
                              # "literal"   = NULL sources replaced with "NA"
  overwrite   = FALSE,        # FALSE = error if filtered rows already have non-NULL values
  quiet       = FALSE         # TRUE = suppress diagnostic messages
)
```

### concat — Arguments

| Argument | Type | Description |
|---|---|---|
| `tbl` | `tidybreed_table` | From `get_table()` + optional `filter()` |
| `col_name` | character | Column name to create or append to; always VARCHAR |
| `source_columns` | character vector | Column names from the same table to concatenate, in order |
| `sep` | character | Length-one non-missing separator; default `"_"` |
| `null_action` | character | `"propagate"` (default): any NULL source → NULL. `"skip"`: NULLs omitted. `"literal"`: NULLs written as `"NA"`. |
| `overwrite` | logical | `FALSE` (default): error if any filtered row already has non-`NULL` value. `TRUE`: overwrite. |
| `quiet` | logical | `FALSE` (default): print assignment diagnostics. `TRUE`: suppress. |

### concat — Constraints

- No `include_leftover` — every filtered row gets a value by definition.
- All `source_columns` must exist in the table; error if any are missing.
- Output is always VARCHAR regardless of source column types.
- Same-table only. For cross-table CGs: collect source table, join in R, write
  the joined column into the target table with `mutate_table()`, then concatenate.
- `null_action = "propagate"` (default) is safest for sib and CG labels where
  a partial key must not silently match another partial key.

### concat — Examples

```r
# Contemporary group from sex, barn, and generation
pop |> get_table("ind_meta") |> filter(gen >= 1L) |>
  mutate_group_concatenate(
    col_name       = "cg",
    source_columns = c("sex", "barn", "gen")
  )

# Full-sib group; founders get NULL (no parents) via null_action = "propagate"
pop |> get_table("ind_meta") |>
  mutate_group_concatenate(
    col_name       = "full_sib_group",
    source_columns = c("id_parent_1", "id_parent_2")
  )

# Paternal half-sib group
pop |> get_table("ind_meta") |> filter(!is.na(id_parent_1)) |>
  mutate_group_concatenate(
    col_name       = "pat_hs_group",
    source_columns = c("id_parent_1")
  )

# Custom separator
pop |> get_table("ind_meta") |>
  mutate_group_concatenate(
    col_name       = "cg",
    source_columns = c("sex", "gen"),
    sep            = ":"
  )
```

---

## Function Comparison

| | `mutate_group_seq()` | `mutate_group_named()` | `mutate_group_concatenate()` |
|---|---|---|---|
| Output type | INTEGER | VARCHAR | VARCHAR |
| Assignment | Random | Random | Deterministic |
| Group names | Auto (1, 2, 3…) | User-defined | Derived from data |
| Auto-continuation | Yes (MAX + 1 or `start`) | No | No |
| `include_leftover` | Yes | Yes (counts only) | No |
| `underfill_action` | No | Yes (counts only) | No |
| `prop_method` | No | Yes (proportions only) | No |
| `null_action` | No | No | Yes |
| `seed` | Yes | Yes | No |
| `overwrite` | Yes | Yes | Yes |
| `quiet` | Yes | Yes | Yes |
| Use cases | Pens, test groups | Farms, barns, gate cuts | CGs, sib groups |

---

## Implementation Notes

- All three functions validate arguments, check overwrite status, and verify
  column types **before** any call to `sample()` or `set.seed()`, so a failed
  call does not advance caller RNG state.
- Use `withr::with_seed(seed, { ... })` for scoped seed handling so the caller's
  global RNG is restored after the call.
- Primary keys are resolved via `TABLE_PRIMARY_KEYS[[table_name]]`, not
  hard-coded `id_ind`. Error early if the table has no registered primary key.
- The keyed tibble path in `mutate_table()` is used for all per-row writes:
  `tibble(pk_col = pks, col_name = values)`. This makes the ID-to-label join
  explicit and safe regardless of row order.
- Internal `mutate_table()` messages are suppressed; only group-level diagnostics
  are emitted by the group function itself.
- **`mutate_group_concatenate()` NULL handling**: DuckDB's `CONCAT_WS` skips
  `NULL` arguments rather than propagating them. For `null_action = "propagate"`,
  use an explicit SQL pattern:
  `CASE WHEN (col1 IS NULL OR col2 IS NULL) THEN NULL ELSE CONCAT_WS(sep, CAST(col1 AS VARCHAR), ...) END`.
  All source columns cast to `VARCHAR` before concatenation.
- **Shared internal helpers** to avoid duplicating fragile table logic:
  - `resolve_group_target(tbl, col_name, output_type, overwrite)` — checks type, schema, and overwrite
  - `collect_group_pks(tbl, pk_col)` — collects and PK-sorts filtered IDs
  - `check_existing_group_values(conn, table_name, col_name, pk_col, pks)` — overwrite pre-check
  - `write_group_values(tbl, col_name, pk_col, values_by_pk)` — keyed write via `mutate_table()`
  - `summarize_group_assignment(values)` — counts and group-size summary for diagnostics

---

## Critique Responses

### Second-pass critique — accepted and incorporated above

**Row ordering contradiction resolved**: removed the "arrange() is respected"
promise. v1 always sorts filtered PKs by primary key before shuffling.
Deterministic, unambiguous, testable. An `order_by` argument is deferred.

**Primary key hardcoding fixed**: implementation notes now use
`TABLE_PRIMARY_KEYS[[table_name]]` throughout. Error raised if no PK is
registered.

**`start` collision semantics clarified**: `overwrite = FALSE` only checks
whether the **filtered rows** already have non-`NULL` values. The same label
integers existing in unfiltered rows (e.g., gen == 1 having pen 1-5 while
gen == 2 also starts at pen 1) is label reuse, not row overwrite, and does not
require `overwrite = TRUE`. This is the most important behavior fix from the
second critique.

**`n_groups` + `include_leftover` precision added**: exact floor/remainder
rules documented. Balanced distribution for `include_leftover = TRUE`; strict
floor + NULL remainder for `include_leftover = FALSE`.

**Proportions tie-breaking documented**: largest-remainder with stable
left-to-right `group_names` order as tie-breaker. `n = 5`,
`proportions = c(1/3, 1/3, 1/3)` → groups 1 and 2 get 2 animals, group 3
gets 1.

**Misleading cross-table example replaced**: removed the non-working
`mutate_table(barn = barn)` code. Replaced with a plain-language description
of the collect/join/`mutate_table()` pattern.

**`quiet = FALSE` added in v1**: added to all three signatures rather than
deferring. Tests can suppress messages with `quiet = TRUE`.

**Empty-filter behavior documented**: warn and return unchanged; do not create
new columns from a zero-row filter.

**Modeling note added**: short section in Shared Conventions explaining which
tables produce model-valid group columns vs. observation or locus groups.

**Validate before sampling**: documented in Implementation Notes and Shared
Conventions.

**`withr::with_seed()` for RNG isolation**: documented in Implementation Notes.

**Column type enforcement**: added to Shared Conventions.

**Diagnostic duplication**: internal `mutate_table()` messages suppressed; only
group-level diagnostics print.

### First-pass critique — still rejected

**Restricting random assignment to `id_ind` tables only**: still rejected.
Generic function design is a core principle. The modeling note in Shared
Conventions addresses the user expectation gap without restricting the API.

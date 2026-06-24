# Second-Pass Critique of `mutate_group_*()` Plan

## Source Note

This is a review of `plans/mutate_group.md` as updated after the earlier
`define_group` critique. The plan now directly addresses many of the original
concerns, especially reproducibility, overwrite safety, sequential restart,
NULL handling, proportion tolerance, diagnostics, and the same-table boundary
for concatenation.

## Overall Assessment

The revised plan is much stronger and is close to implementable. The
three-function split still makes sense:

- `mutate_group_seq()` for random assignment to generated integer groups.
- `mutate_group_named()` for random assignment to fixed character labels.
- `mutate_group_concatenate()` for deterministic character groups from existing
  columns.

The most important design improvement is that the functions are now clearly
data mutation helpers, not model-definition helpers. That fits tidybreed's
existing `get_table()` / `mutate_table()` style and avoids overloading
`define_*()` semantics.

The remaining issues are no longer about the high-level API shape. They are
mostly consistency and implementation-contract problems that should be resolved
before coding, because they affect reproducibility, keyed updates, and user
expectations.

## Resolved From Prior Critique

The updated plan satisfactorily resolves these earlier points:

- **Naming and scope**: `mutate_group_*()` is the right name for functions that
  write table columns. The plan now explicitly says they do not register model
  effects or group semantics.
- **Reproducibility**: `seed = NULL` is present for the random functions, with
  the expected R convention that `NULL` uses the caller's RNG state.
- **Overwrite safety**: `overwrite = FALSE` is now the default on all three
  functions, with an error for filtered rows that already have non-`NULL`
  values.
- **Sequential restart**: `start = NULL` keeps `MAX(col_name) + 1`, while
  `start = 1L` supports reuse of local group IDs by generation or another
  caller-controlled context.
- **Concatenation NULL behavior**: `null_action` is included and
  `"propagate"` is the default. The plan correctly warns not to rely on bare
  DuckDB `CONCAT_WS` for this case.
- **Proportion tolerance**: the plan accepts sums within
  `sqrt(.Machine$double.eps)` and normalizes before assignment.
- **Diagnostics**: the plan now requires row counts, leftovers, and group-size
  summaries.
- **Cross-table concatenation boundary**: same-table-only behavior is now
  explicit, with a materialization workflow deferred to `mutate_table()`.
- **Deferred features**: `strata`, `balance_by`, `by`, group metadata,
  observation-level grouping, and sequence label templates are appropriately
  out of v1 scope.

## High-Priority Remaining Issues

### 1. Row ordering contract is internally contradictory

The plan says:

- filtered primary keys are "collected and sorted by primary key before random
  shuffling"; and
- upstream `arrange()` is respected.

Both cannot be the canonical pre-shuffle order unless the implementation
explicitly distinguishes arranged and unarranged tables. If you always sort by
primary key, `arrange()` is ignored. If you always collect the lazy table's
current order, unrelated database row order may leak in when the user did not
arrange.

Recommendation:

- Define the rule precisely:
  - default: collect the primary key ordered by primary key;
  - if the `tidybreed_table` has an explicit arrange/order operation, collect in
    that order.
- Verify whether the current `tidybreed_table` object retains enough metadata
  to know that `arrange()` was explicitly called. `arrange.tidybreed_table()`
  currently modifies the lazy table, but it does not appear to mark the object
  as user-arranged.
- If explicit arrange detection is not easy, choose one simpler v1 behavior:
  always order by primary key and remove the "arrange is respected" promise.

This should be settled before implementation because reproducibility tests will
depend on it.

### 2. The keyed update examples need to use the target table's primary key

The implementation notes say to use the `mutate_table()` tibble path and show an
`id_ind`-keyed tibble. That is safe for `ind_meta`, but not for "any tidybreed
table." Local `mutate_table()` requires the tibble to contain the target table's
registered primary key from `TABLE_PRIMARY_KEYS`.

Examples:

- `ind_meta` uses `id_ind`.
- `genome_meta` uses `locus_id`.
- `ind_phenotype` uses `id_phenotype`.
- user-defined tables may have their own optional primary key.

Recommendation:

- Replace "keyed tibble path (`tibble(id_ind = ids, col = values)`)" with
  "keyed tibble path using the target table's registered primary key."
- Add an implementation helper that resolves `pk_col <- TABLE_PRIMARY_KEYS[[table_name]]`.
- Error early if the table has no registered primary key. This matters for
  user-defined tables created without `primary_key`.

Without this, the generic "any tidybreed table" promise is easy to implement
incorrectly.

### 3. `start` collision semantics conflict with intended label reuse

The plan says `start = 1L` lets users restart pen numbering per generation, but
also says a supplied `start` that collides with existing values errors unless
`overwrite = TRUE`.

Those rules conflict for normal reuse:

```r
pop |> get_table("ind_meta") |> filter(gen == 1L) |>
  mutate_group_seq(col_name = "pen", n_per_group = 10L, start = 1L)

pop |> get_table("ind_meta") |> filter(gen == 2L) |>
  mutate_group_seq(col_name = "pen", n_per_group = 10L, start = 1L)
```

The second call should be valid if it only writes previously `NULL` rows, even
though values `1`, `2`, ... already exist in the same column for generation 1.
`overwrite = TRUE` is about replacing existing values in affected rows; it
should not be required merely because group labels are reused in unaffected
rows.

Recommendation:

- Distinguish **row overwrite collisions** from **label reuse collisions**.
- Keep the overwrite error when filtered rows already have non-`NULL`
  `col_name` values.
- Do not error on label reuse elsewhere in the table by default.
- If global uniqueness is desired, add a separate future argument such as
  `unique_labels = TRUE`; do not overload `overwrite`.

This is the most important behavior fix in the current plan.

### 4. `n_groups` and `leftover` semantics need a sharper definition

The plan says `n_groups` splits the filtered set into `N` equal groups and that
`leftover = FALSE` leaves remainder animals `NULL`. That is understandable but
still ambiguous in edge cases.

For `n = 11`, `n_groups = 5`, `leftover = FALSE`, possible outcomes include:

- five groups of 2 and one unassigned animal;
- four groups of 2, one group of 3, and no unassigned animals;
- two groups of 3 and three groups of 2, and no unassigned animals.

The plan implies the first, but it should say so explicitly.

Recommendation:

- For `n_groups` with `leftover = FALSE`: group size is `floor(n / n_groups)`;
  `n %% n_groups` animals are left `NULL`.
- For `n_groups` with `leftover = TRUE`: use balanced sizes where the first
  remainder groups get one extra animal after shuffling, or explicitly say the
  final group gets all leftovers.
- Add tests for `n = 11`, `n_groups = 5`, both leftover settings.

The same clarity is needed for variable `n_per_group = c(min, max)`: specify
whether the sampled sequence of group capacities can overshoot `n`, and exactly
how the last partially filled sampled group is handled when `leftover = FALSE`.

### 5. `mutate_group_named()` should define exact rounding tie-breaking

Largest-remainder rounding is a good choice for `method = "exact"`, but the
plan does not define tie-breaking. This affects reproducibility for inputs like
`n = 5`, `proportions = c(1/3, 1/3, 1/3)`.

Recommendation:

- Use stable group-name order as the tie-breaker for equal remainders.
- Document that the last group only absorbs rounding in the `underfill =
  "proportional"` path for `counts`, not in the proportions largest-remainder
  path.
- Add tests for tiny `n`: 1 animal / 3 groups, 2 animals / 3 groups, and equal
  fractional remainders.

### 6. The cross-table example is currently misleading

The plan's cross-table example says to bring `barn` from `ind_phenotype` into
`ind_meta`, but the code calls:

```r
pop <- pop |>
  get_table("ind_phenotype") |>
  filter(phenotype_name == "ADG") |>
  mutate_table(barn = barn)
```

That mutates `ind_phenotype`, not `ind_meta`. It also assumes `barn` already
exists in `ind_phenotype`, so it does not demonstrate materializing a column
across tables.

Recommendation:

- Replace this with a conceptual example or a real supported helper path.
- If no existing helper supports a keyed cross-table update from
  `ind_phenotype` to `ind_meta`, state the limitation plainly and avoid showing
  non-working code.
- Consider pointing users to a collect/join/keyed-`mutate_table()` pattern only
  if that pattern is actually supported for the target table's primary key.

## Medium-Priority Issues

### Diagnostics should avoid duplicate noise from `mutate_table()`

`mutate_table()` already messages when it adds a column, fills `NULL`s, or
replaces existing values. The group functions also plan to message rows
assigned, rows left `NULL`, and group sizes.

Recommendation:

- Decide whether group functions will let both message layers print, or whether
  they need an internal quiet path in `mutate_table()`.
- If both print, keep group-level diagnostics compact and avoid repeating the
  same row-count sentence twice.
- If a future `quiet = FALSE` is likely, consider adding it in v1 rather than
  promising a "future argument" that tests and examples cannot use.

### Empty-filter behavior should be specified

`mutate_table()` warns and returns unchanged when a filter matches zero rows.
The group functions should specify whether they mirror that behavior or error.

Recommendation:

- Warn and return the population unchanged for zero filtered rows, matching
  `mutate_table()`.
- Do not create a new group column on an empty filtered subset unless the user
  supplies an explicit typed schema path. Otherwise a typo in a filter can
  silently alter schema.

### Table support should be generic but honest

The revised plan rejects restricting random assignment to `id_ind` tables, which
is reasonable for a generic mutation helper. But the documentation should still
make the model-facing limitation prominent:

- group columns used for individual-level phenotypes should usually live in
  `ind_meta` or another table that can be joined by `id_ind`;
- group columns on `ind_phenotype` are observation-record groups, not
  individual groups;
- group columns on `genome_meta` are locus groups, not animal groups.

Recommendation:

- Keep the functions generic.
- Add a short "Modeling note" section with examples of valid source-table
  choices for `define_effect_random(source_column = ...)` and phenotype formulas.

### Existing-value checks should happen before random sampling

With `overwrite = FALSE`, the function should check for non-`NULL` values in the
affected rows before consuming RNG state. Otherwise a failed call changes the
caller's random stream.

Recommendation:

- Validate arguments, target column state, existing values, and table keys
  before calling `set.seed()` or `sample()`.
- If `seed` is supplied, use `withr::with_seed()` or equivalent local RNG-state
  handling so the caller's RNG stream is not accidentally reset.
- Add tests that failed calls do not advance RNG state if that guarantee is
  intended.

### Column type behavior for existing columns needs enforcement

The plan says sequential output is always `INTEGER`, named and concatenated
outputs are always `VARCHAR`. That is clear for newly created columns, but not
for existing columns.

Recommendation:

- If `col_name` already exists, check that its DuckDB type is compatible before
  writing.
- Error if `mutate_group_seq()` targets an existing `VARCHAR` column, or if
  named/concatenate target an incompatible numeric column.
- Add tests for attempting to append to an existing wrong-type column.

## Implementation Recommendations

Use small shared helpers so the three functions do not duplicate fragile table
logic:

```r
resolve_group_target(tbl, col_name, output_type, overwrite = FALSE)
collect_group_pks(tbl, pk_col, ordering = c("pk", "arranged"))
check_existing_group_values(conn, table_name, col_name, pk_col, pks)
write_group_values(tbl, col_name, pk_col, values_by_pk)
summarize_group_assignment(values)
```

Implementation details to keep explicit:

- Use the target table's registered primary key, not hard-coded `id_ind`.
- Use the keyed tibble path in `mutate_table()` for per-row values.
- Preserve R RNG state when `seed` is supplied.
- Check overwrite/type/key errors before sampling.
- For `null_action = "propagate"`, implement an explicit `CASE WHEN source IS
  NULL THEN NULL ELSE ... END` expression or do the concatenation after
  collecting keyed source values.
- Validate all user-supplied SQL identifiers with existing SQL helper rules.

## Test Recommendations

Add focused unit tests for:

- same seed + same filtered rows gives identical assignment;
- different upstream insertion/order does not affect default assignment;
- explicit arrange behavior, if the plan keeps that promise;
- `overwrite = FALSE` errors before sampling and before writing;
- appending to only `NULL` rows succeeds;
- reusing `start = 1L` on a different filtered subset succeeds without
  requiring `overwrite = TRUE`;
- `n_groups` edge cases with leftovers;
- variable `n_per_group` edge cases;
- named `counts` underfill modes;
- named `proportions` exact rounding, tolerance, invalid values, and tiny `n`;
- wrong existing column type errors;
- zero-row filters;
- concatenation `null_action = "propagate"`, `"skip"`, and `"literal"`;
- use on `ind_meta`, `genome_meta`, `ind_phenotype`, and a user-defined table
  with and without a primary key.

Add at least one integration test where a generated `pen_id` or contemporary
group column is consumed by a model-facing path, such as a phenotype formula or
`define_effect_random(source_column = ...)`. That catches type and table
assumption mismatches that simple count tests will miss.

## Bottom Line

The updated plan has incorporated the major critique and is directionally ready
for implementation. Before coding, fix the remaining contract details:

- make row ordering unambiguous;
- use the target table's primary key everywhere;
- allow label reuse across unaffected rows when `start` is supplied;
- precisely define `n_groups` leftover behavior;
- document proportions tie-breaking;
- replace the misleading cross-table example.

Once those are addressed, the v1 scope is appropriately narrow and should fit
tidybreed's existing table-mutation architecture.

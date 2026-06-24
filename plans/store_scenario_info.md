# Store Scenario Info — Design Plan

**Status:** Initial draft for discussion  
**Related:** `plans/YAML_PLAN.md` (covers the YAML file format; this covers DB storage)

---

## Why Store Scenario Params in a Table?

The YAML file is the *authoring* format — easy to edit and diff. But once
a simulation runs, the parameters should live inside the `.duckdb` file alongside
the results so that:

1. **Reproducibility** — open a `.duckdb` file six months later and know
   exactly what was run.  No hunting for the original YAML.
2. **`restore_pop()` completeness** — design principle: restarts must not require
   the user to re-supply options. If scenario parameters are in the DB, they
   are available to any restore logic.
3. **Cross-scenario queries** — if multiple scenarios are stored or summarized
   together, you can ask "which scenarios used `n_sires > 5`?" in plain SQL.

---

## The Shape of the Problem

Consider the `age_at_puberty_sexed_semen.yaml` file already in use:

```yaml
general:
  n_reps: 5
  start_date: "2026-01-01"
  cycle_len: 7
  mean_puberty_age: 220

selection:
  n_sires: 10
  n_dams_per_breeding: 25
  genomic_selection: false

traits:
  - name: "AP"
  - name: "ADG"
```

The parameters are:

- **Heterogeneous types** — integers, doubles, strings, booleans, dates, and
  lists (the `traits:` section)
- **Grouped** — natural sections like `general`, `selection`, `testing`,
  `culling`
- **Some lists / nested** — `traits:` is a list of objects, `chr_len_Mb` can be
  a vector of per-chromosome values

Your proposed columns — `scenario_name`, `group`, `argument_name`,
`argument_value` — map onto this naturally for scalar values. The main question
is how to handle types and lists.

---

## Three Design Options

### Option A: EAV + `argument_type` column (recommended starting point)

```
scenario_params

  id_scenario_param  INTEGER   PK (auto-incrementing)
  scenario_name      VARCHAR   e.g. "age_at_puberty"
  group_name         VARCHAR   e.g. "general", "selection", "mating"
  argument_name      VARCHAR   e.g. "n_sires"
  argument_type      VARCHAR   "integer" | "double" | "character" |
                               "boolean" | "date" | "json"
  argument_value     VARCHAR   always stored as text; cast on retrieval
```

Every value is stored as a VARCHAR string. The `argument_type` column carries
enough information to safely cast back in R:

```r
# retrieval helper sketch
parse_scenario_value <- function(val, type) {
  switch(type,
    integer   = as.integer(val),
    double    = as.double(val),
    character = val,
    boolean   = as.logical(val),
    date      = as.Date(val),
    json      = jsonlite::fromJSON(val)
  )
}
```

For lists and nested objects (the traits section, per-chromosome lengths), 
serialize to JSON and set `argument_type = "json"`:

```
scenario_name  group_name  argument_name  argument_type  argument_value
age_at_puberty general     n_reps         integer        5
age_at_puberty general     start_date     date           2026-01-01
age_at_puberty general     cycle_len      integer        7
age_at_puberty selection   n_sires        integer        10
age_at_puberty selection   genomic_sel    boolean        FALSE
age_at_puberty traits      trait_list     json           ["AP","ADG","BF","NW"]
age_at_puberty genome      n_loci         integer        10000
age_at_puberty genome      chr_len_Mb     json           [50,50,50,...] (18 values)
```

**Pros:** single table, simple schema, survives arbitrary new groups/keys,
easy to add new scenarios without schema changes.  
**Cons:** all SQL queries need to `CAST(argument_value AS INTEGER)` etc.; no
column-level type enforcement at the DB layer.

---

### Option B: Typed columns (wide format)

```
scenario_params

  id_scenario_param  INTEGER   PK
  scenario_name      VARCHAR
  group_name         VARCHAR
  argument_name      VARCHAR
  int_value          INTEGER
  dbl_value          DOUBLE
  chr_value          VARCHAR
  bool_value         BOOLEAN
  date_value         DATE
  json_value         VARCHAR   (for lists / nested objects)
```

One column per type; only the matching column is populated for each row.

**Pros:** SQL aggregations work without casting (`WHERE int_value > 5`).  
**Cons:** many NULLs per row; awkward to read; `infer_duckdb_type()` would 
need a five-column insert path instead of one.

---

### Option C: JSON blob per group

```
scenario_groups

  id_scenario_group  INTEGER   PK
  scenario_name      VARCHAR
  group_name         VARCHAR
  params_json        VARCHAR   full JSON object for this group
```

Each group becomes one row containing all its key-value pairs as a JSON blob.

**Pros:** fewest rows; trivially handles nested/list structures; easy to
round-trip a whole group.  
**Cons:** no SQL-level filtering on individual params; opaque to anything that
reads the DB without JSON parsing.

---

## Recommendation for the Initial Plan

**Start with Option A** (EAV + `argument_type`). It is the most SQL-friendly
for per-parameter queries while remaining schema-stable as the YAML evolves.
JSON serialization handles the two edge cases (list values, nested objects)
without a separate table.

If cross-scenario parameter queries turn out to be rare in practice, Option C
(JSON blob per group) would also be workable and far simpler to write.

---

## The Traits Section: Special Case

The `traits:` section in the YAML has its own sub-structure (trait names,
heritabilities, QTL counts, etc. per the full YAML_PLAN schema). Two options:

1. **Serialize as JSON** into a single row in `scenario_params`:
   `group = "traits"`, `argument_name = "trait_list"`, `argument_type = "json"`,
   `argument_value = '[{"name":"ADG","h2":0.25,"n_qtl":500},...]'`

2. **Separate table** `scenario_traits` with one row per (scenario × trait):

   ```
   scenario_traits

     id_scenario_trait  INTEGER   PK
     scenario_name      VARCHAR
     trait_name         VARCHAR
     argument_name      VARCHAR   e.g. "h2", "n_qtl", "index_weight"
     argument_type      VARCHAR
     argument_value     VARCHAR
   ```

Option 2 makes it easy to join `scenario_traits` against `trait_meta` to cross-
reference what parameters were used for each trait. Worth doing if traits are
queried per-trait often; otherwise JSON in Option A is enough.

---

## Integration: When Does the Table Get Written?

The logical place is inside `load_scenario()` (per YAML_PLAN.md), but only
**after** `initialize_genome()` creates the DuckDB file.  A possible split:

1. `load_scenario(path)` → validates YAML, returns `tidybreed_scenario` object
   (no DB writes yet; DB may not exist)
2. `initialize_from_scenario(scenario, ...)` → creates the DB, calls 
   `initialize_genome()`, then immediately calls an internal
   `.persist_scenario_params(pop, scenario)` that flattens the scenario list
   into `scenario_params` rows

This keeps the table write co-located with DB creation and avoids writing to
a DB that doesn't exist yet.

---

## Open Questions

1. **Multiple runs / overwrite?** If the user runs the same scenario YAML twice
   (re-run with a different seed), should the old rows be overwritten, or should
   rows be appended with a `run_id` / `eval_timestamp`? Appending makes the table
   a full audit trail but complicates "what were the params for this pop?"

2. **Scenario name vs. file name** — the YAML already has a `scenario_name:` key.
   Should the DB use that string as the FK, or derive it from the file path?
   Currently the vignette YAML has `scenario_name: "age_at_puberty"` which is
   clean. File-path derivation would work as a fallback.

3. **Free-form `simulation:` section** — per YAML_PLAN, this section is not
   validated by the package. Should it still be persisted verbatim (as a JSON
   blob row with `group_name = "simulation"`) for reproducibility?  Almost
   certainly yes.

4. **Covariance matrices** — these are 2-D arrays. The `covariances:` section in
   the full YAML schema has named matrices. These don't fit a single row well.
   Options: JSON serialize the full matrix, or store as a long-format sub-table
   `scenario_cov_matrix (scenario_name, effect_name, trait_name_1, trait_name_2, value)`.
   The latter mirrors the existing `trait_var_comp` design and may be redundant
   with it.

5. **Schema version** — should a `schema_version` row be written so that future
   `open_pop()` / `restore_pop()` can detect and handle format changes?

---

## Minimal Viable Table (first implementation target)

If the goal is just reproducibility and restore completeness, the smallest useful
starting point is:

```
scenario_params: id_scenario_param, scenario_name, group_name,
                 argument_name, argument_type, argument_value
```

Populated once at `initialize_from_scenario()` time for all scalar and string
values in the YAML. Lists and matrices serialized to JSON with
`argument_type = "json"`. The traits list gets one JSON row rather than a
separate table.

That covers reproducibility. The separate `scenario_traits` table and typed
column Option B can be added later if cross-scenario per-trait queries become
common.

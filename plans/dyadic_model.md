# Dyadic Phenotype Model — Design Proposal

**Context**: Implementing simulation support for dyadic social interaction traits
(e.g., pig attack duration) as requested by Iowa State University collaborator.

---

## What a Dyadic Phenotype Is

A dyadic phenotype is an observation that belongs to a **pair of individuals**, not a
single individual. The canonical example is "time attacking" (minutes per day) observed
for every directed pair within a pen on day 1 post-mixing.

The dataset structure is all ordered pairwise interactions within pen:
- A attacks B → one record (A as actor, B as partner)
- B attacks A → separate record (B as actor, A as partner)
- N animals in a pen → N×(N-1) directed records

**The phenotype model:**

```
y_ij = mean + TBV_actor(i) + TBV_partner(j) + PE_actor(i) + PE_partner(j) + ε_ij
```

Where:
- `i` = actor (attacker)
- `j` = partner (receiver)
- `TBV_actor(i)` = additive genetic value of animal i for the "attacker" trait
- `TBV_partner(j)` = additive genetic value of animal j for the "receiver" trait
- `PE_actor(i)` = permanent environment effect of animal i in the actor role
- `PE_partner(j)` = permanent environment effect of animal j in the partner role
- `ε_ij` = observation-specific residual (pair × time point)

Two covariance structures:
- **2×2 G matrix** between `(attack_direct, attack_received)` genetic traits
- **2×2 PE matrix** between `(actor PE, partner PE)` per-individual role effects

---

## What Already Works — Don't Change

- The 2×2 genetic covariance matrix for `(attack_direct, attack_received)`:
  `define_effect_cov_matrix("gen_add", G)` already handles this.
- Correlated QTL effect sampling: `define_additive_effects(c("attack_direct", "attack_received"), G = G)`
  already works for the multi-trait case.
- The `phenotype_components` mechanism with `contributor_type` is the right abstraction —
  we just need two new valid type values (see below).

---

## Schema Changes

### New table: `dyadic_phenotype`

Do **not** add a nullable `id_partner` to `ind_phenotype`. A nullable column that
fundamentally changes the row's semantics forces every downstream query to filter on
`WHERE id_partner IS NULL` vs. `IS NOT NULL`, and silently breaks naive queries.
A separate table gives explicit, unambiguous semantics.

```sql
CREATE TABLE dyadic_phenotype (
  id_dyadic      INTEGER PRIMARY KEY,
  id_ind_1       VARCHAR NOT NULL,    -- actor / initiator
  id_ind_2       VARCHAR NOT NULL,    -- partner / receiver
  phenotype_name VARCHAR NOT NULL,
  pheno_value    DOUBLE,
  pheno_number   INTEGER DEFAULT 1,   -- increments per (id_ind_1, id_ind_2, phenotype_name)
  group_name     VARCHAR              -- resolved group value (e.g. "pen_42") shared by both individuals
  -- user columns added via ... in add_phenotype() or mutate_table() afterwards
)
```

`add_phenotype()` dispatches to this table when `phenotype_class = "dyadic"`.

**`group_name` rationale**: storing the resolved group value (e.g. `"pen_42"`, not the
column name) in every row lets users:
- Filter or join by pen without re-joining `ind_meta`
- Add a random pen-level variance component via `define_effect_random()` referencing
  this column as the grouping variable
- Compute pair-level fixed covariates (e.g. weight difference between id_ind_1 and
  id_ind_2) via `mutate_table()` after phenotype generation, stored directly in this table

---

### New table: `ind_pe`

PE effects are individual-level — the same animal always has the same tendency to be an
aggressive attacker or reactive receiver. They must be sampled once and reused (not
resampled each call to `add_phenotype()`). They don't belong in `ind_meta` (identity
metadata). A dedicated table mirrors the `ind_tbv` design.

```sql
CREATE TABLE ind_pe (
  id_pe          INTEGER PRIMARY KEY,
  id_ind         VARCHAR NOT NULL,
  phenotype_name VARCHAR NOT NULL,   -- which dyadic phenotype these PE effects belong to
  role_name      VARCHAR NOT NULL,   -- "dyadic_actor" or "dyadic_partner"
  pe_value       DOUBLE NOT NULL
)
```

Logical key: `(id_ind, phenotype_name, role_name)`.

Sampled jointly from `MVN(0, PE_cov)` the first time `add_phenotype()` runs for a
given phenotype. On subsequent calls, already-sampled individuals are skipped
(`overwrite_pe = FALSE` default). Set `overwrite_pe = TRUE` to resample (e.g. after
updating the PE covariance).

---

### Changes to `phenotype_meta` — 4 new columns

```sql
phenotype_class VARCHAR DEFAULT 'individual',  -- 'individual' or 'dyadic'
group_column    VARCHAR,    -- column in group_table that assigns individuals to groups
group_table     VARCHAR DEFAULT 'ind_meta',    -- table containing group_column
directed        BOOLEAN DEFAULT TRUE           -- TRUE = A→B and B→A are separate records
```

**Why `group_column` and `group_table` belong in `phenotype_meta`** (not in
`phenotype_components`): they govern pair *enumeration* — a phenotype-level concern
that happens before any component is evaluated. `add_phenotype()` reads these first
to know where to look for group membership and generate all pairs. This mirrors how
the existing SGE group logic uses `group_column`/`group_table` in `phenotype_components`,
but here it applies to the whole dyadic phenotype.

`group_column` references the column written by the `define_group_*()` functions
(e.g. `pen_id` in `ind_meta`). No coupling between `define_phenotype()` and the group
machinery — just a name reference.

---

### Changes to `phenotype_var_comp` — 2 new nullable columns

```sql
role_1  VARCHAR,   -- NULL for residual and other non-PE effects
role_2  VARCHAR    -- NULL for residual and other non-PE effects
```

PE covariance is stored with `effect_name = "pe"`:

| effect_name | phenotype_name_1   | phenotype_name_2   | role_1         | role_2         | cov_value |
|-------------|--------------------|--------------------|----------------|----------------|-----------|
| pe          | time_attacking     | time_attacking     | dyadic_actor   | dyadic_actor   | 0.15      |
| pe          | time_attacking     | time_attacking     | dyadic_actor   | dyadic_partner | 0.04      |
| pe          | time_attacking     | time_attacking     | dyadic_partner | dyadic_actor   | 0.04      |
| pe          | time_attacking     | time_attacking     | dyadic_partner | dyadic_partner | 0.12      |

Both (i,j) and (j,i) pairs are stored (consistent with the existing convention).
`role_1` and `role_2` are `NULL` for all existing residual and named random effect rows —
fully backward compatible.

The existing `define_effect_cov_matrix()` routing is extended: when `effect_name = "pe"`,
it routes to `phenotype_var_comp` and writes `role_1`/`role_2` from the matrix
`dimnames`. The `phenotype_name` argument identifies which dyadic phenotype the PE
structure belongs to.

---

### Changes to `phenotype_components` — no schema change

Two new valid values for the existing `contributor_type` column:
- `"dyadic_actor"` — pull TBV of `id_ind_1` for `source_trait_name`
- `"dyadic_partner"` — pull TBV of `id_ind_2` for `source_trait_name`

These are additive to the existing vocabulary (`"self"`, `"dam"`, `"sire"`, `"group"`).

**How TBV sourcing differs from other contributor types:**
- `"self"`, `"dam"`, `"sire"` — pull one TBV per focal individual from `ind_tbv`
- `"group"` — aggregate TBVs of all group-mates (excluding self) from `ind_tbv`
- `"dyadic_actor"` — pull TBV of `id_ind_1` (the attacker in a specific pair)
- `"dyadic_partner"` — pull TBV of `id_ind_2` (the receiver in a specific pair)

The last two are fundamentally pair-level lookups, not individual-level. For a given
directed pair (A→B), `add_phenotype()` looks up A's TBV for the actor trait and B's TBV
for the partner trait separately — there is no aggregation. This makes dyadic TBV
sourcing unique even compared to the maternal weaning weight model (dam/sire TBV is
still one lookup per focal individual; dyadic requires two different lookups per pair).

---

## API

### Defining a dyadic phenotype

```r
# ── Genetic layer (already works today) ───────────────────────────────────────

G <- matrix(c(0.25, 0.06, 0.06, 0.20), 2, 2,
            dimnames = list(c("attack_direct", "attack_received"),
                            c("attack_direct", "attack_received")))

pop <- pop |>
  define_trait("attack_direct",   target_add_var = 0.25) |>
  define_trait("attack_received", target_add_var = 0.20) |>
  define_effect_cov_matrix("gen_add", G)

pop <- pop |>
  get_table("genome_meta") |>
  define_additive_effects(c("attack_direct", "attack_received"), G = G)

# ── Observation layer (new dyadic support) ────────────────────────────────────

pop <- pop |>
  define_phenotype(
    phenotype_name = "time_attacking",
    class          = "dyadic",         # new: "individual" (default) or "dyadic"
    type           = "continuous",
    mean           = 30,               # population mean in minutes
    group_column   = "pen_id",         # column in group_table holding group membership
    group_table    = "ind_meta",       # default; where group_column lives
    directed       = TRUE,             # default; A→B and B→A are separate records
    components     = data.frame(
      source_trait_name = c("attack_direct",  "attack_received"),
      contributor_type  = c("dyadic_actor",   "dyadic_partner"),
      weight            = c(1, 1)
    ),
    residual_var = 5
  )

# ── PE covariance (new routing in define_effect_cov_matrix) ───────────────────

PE <- matrix(c(0.15, 0.04, 0.04, 0.12), 2, 2,
             dimnames = list(c("dyadic_actor", "dyadic_partner"),
                             c("dyadic_actor", "dyadic_partner")))

pop <- pop |>
  define_effect_cov_matrix("pe", PE, phenotype_name = "time_attacking")

# ── Simulate dyadic phenotypes ─────────────────────────────────────────────────

pop <- pop |>
  get_table("ind_meta") |>
  filter(gen == 1L) |>
  add_phenotype("time_attacking")

# ── Add fixed/random effects using stored group_name ──────────────────────────

# Random pen effect (iid variance not explained by genetics/PE)
pop <- pop |>
  define_effect_random("time_attacking", "pen",
    source_column = "group_name",
    source_table  = "dyadic_phenotype",   # group_name is stored here
    variance      = 2.0)

# Fixed covariate: weight difference added to dyadic_phenotype via mutate_table()
pop |>
  get_table("dyadic_phenotype") |>
  filter(phenotype_name == "time_attacking") |>
  mutate_table(wt_diff = ...)   # user computes weight of id_ind_1 minus id_ind_2
```

---

## What `add_phenotype()` Does Internally for Dyadic

When `phenotype_class = "dyadic"`, dispatch to internal helper `.add_dyadic_phenotype()`.
All existing individual-class code paths are unchanged.

### Step 1 — Validate metadata

Read `phenotype_meta` for the requested phenotype. Confirm `phenotype_class = "dyadic"`,
`group_column` is set, and `group_table` is set.

### Step 2 — Enumerate pairs

Resolve group membership for all candidate individuals in the filter set:

```r
# Join filtered individuals to their group values from group_table
focal_groups <- DBI::dbGetQuery(pop$db_conn,
  'SELECT id_ind, "{group_column}" AS group_name
   FROM {group_table}
   WHERE id_ind IN (...)'
)

# Cross-join within the same group_name, exclude self-pairs
pairs <- inner_join(focal_groups, focal_groups, by = "group_name") |>
  filter(id_ind.x != id_ind.y) |>
  rename(id_ind_1 = id_ind.x, id_ind_2 = id_ind.y)

# If directed = FALSE, keep lexicographically first of each unordered pair
if (!directed) pairs <- pairs |> filter(id_ind_1 < id_ind_2)
```

Individuals with a NULL group value are excluded with a warning (count + up to 5
example IDs), same pattern as `missing_component_action = "skip"`.

### Step 3 — Compute TBVs

Call `add_tbv()` (internal) for all required source traits for the union of all
`id_ind_1` and `id_ind_2` values across all pairs. Uses the same pre-computation
pattern as the existing SGE group logic — all TBVs are fetched once before assembly.

### Step 4 — Sample PE effects

1. Read PE covariance from `phenotype_var_comp` where `effect_name = "pe"` and
   `phenotype_name_1 = phenotype_name`. Reconstruct the 2×2 matrix from the 4 rows.
2. Find individuals not yet in `ind_pe` for this phenotype.
3. Sample jointly from `MVN(0, PE_cov)` for new individuals only.
4. Write new rows to `ind_pe`; skip existing individuals unless `overwrite_pe = TRUE`.
5. If no PE covariance is defined, skip this step entirely (PE contribution = 0).

### Step 5 — Assemble phenotype values

For each pair `(i, j)`:

```
y_ij = mean
     + TBV_actor(i)   [from ind_tbv: id_ind=i, trait_name="attack_direct"]
     + TBV_partner(j) [from ind_tbv: id_ind=j, trait_name="attack_received"]
     + PE_actor(i)    [from ind_pe:  id_ind=i, role_name="dyadic_actor"]
     + PE_partner(j)  [from ind_pe:  id_ind=j, role_name="dyadic_partner"]
     + ε_ij           [drawn from N(0, residual_var)]
```

Multiple `"dyadic_actor"` or `"dyadic_partner"` rows in `phenotype_components` are
summed with their respective weights, same as composite traits.

### Step 6 — Write to `dyadic_phenotype`

Append rows with columns: `id_ind_1`, `id_ind_2`, `phenotype_name`, `pheno_value`,
`pheno_number` (auto-incremented per pair), `group_name`, plus any scalar `...` from
the `add_phenotype()` call.

---

## The `filter` Semantics for Dyadic

`get_table("ind_meta") |> filter(gen == 1L) |> add_phenotype("time_attacking")`

The filter specifies **which individuals participate in pair enumeration**. All ordered
pairs within `(filtered_set × filtered_set)` sharing the same group value are generated.
Individuals outside the filter are not included as actor or partner.

---

## Fixed Effects for Dyadic Phenotypes

Dyadic phenotypes introduce three distinct kinds of fixed effects that don't all fit
the existing `define_effect_fixed_class()` / `define_effect_fixed_cov()` pattern,
which assumes the source column is looked up by a single `id_ind`.

### Type 1 — Individual-level, applied to actor or partner

These are properties of `id_ind_1` or `id_ind_2` individually, looked up from
`ind_meta` (or any `id_ind`-keyed table). Examples: sex of the attacker, age of
the receiver, generation of either animal.

Proposed: extend `define_effect_fixed_class()` and `define_effect_fixed_cov()` with
a `contributor_role` argument (`"dyadic_actor"` or `"dyadic_partner"`) that tells
`add_phenotype()` which individual's column to fetch.

```r
# Sex of the attacker as a fixed class effect
pop <- pop |>
  define_effect_fixed_class(
    "time_attacking",
    effect_name      = "sex_actor",
    source_column    = "sex",
    source_table     = "ind_meta",
    contributor_role = "dyadic_actor",    # look up by id_ind_1
    levels_json      = '{"M": 5, "F": 0}'
  )
```

### Type 2 — Pair-level: a property of the (i, j) combination

The strongest and most important example: **were these two animals previous pen
mates?** This is not a property of i alone or j alone — it is a property of the
specific pair `(id_ind_1, id_ind_2)` and must be looked up by both IDs together.

This is a genuinely new concept relative to the existing individual-level effect
machinery.

**Proposed schema — `ind_pair_meta` table (new, user-populated):**

```sql
CREATE TABLE ind_pair_meta (
  id_ind_1  VARCHAR NOT NULL,
  id_ind_2  VARCHAR NOT NULL,
  -- user columns added via mutate_table() or direct insert
  -- e.g.:
  -- prev_pen_mates  BOOLEAN
  -- familiarity_days INTEGER
  PRIMARY KEY (id_ind_1, id_ind_2)
)
```

This mirrors the `ind_meta` pattern — a wide table keyed by the pair, with user
columns added freely. The user populates it before calling `add_phenotype()`.

When `source_table = "ind_pair_meta"` is supplied to `define_effect_fixed_class()`,
`add_phenotype()` joins on `(id_ind_1, id_ind_2)` instead of a single `id_ind`.

**Full previous pen mates example:**

```r
# 1. Create the pair metadata table and populate it
#    (user computes which pairs were previously housed together from simulation history)
pop <- pop |>
  get_table("ind_pair_meta") |>
  mutate_table(prev_pen_mates = NA)   # pre-declare schema

# ... user fills in TRUE/FALSE for each (id_ind_1, id_ind_2) pair based on
#     their group history from a prior generation's group assignments

# 2. Define the fixed class effect
pop <- pop |>
  define_effect_fixed_class(
    "time_attacking",
    effect_name   = "prev_pen_mates",
    source_column = "prev_pen_mates",
    source_table  = "ind_pair_meta",   # signals pair-level lookup
    levels_json   = '{"true": -15, "false": 0}'   # familiar pairs fight less
  )

# 3. add_phenotype() automatically joins ind_pair_meta on (id_ind_1, id_ind_2)
#    and applies the shift before sampling residuals
pop <- pop |>
  get_table("ind_meta") |>
  filter(gen == 1L) |>
  add_phenotype("time_attacking")
```

The model then becomes:

```text
y_ij = mean
     + TBV_actor(i) + TBV_partner(j)
     + PE_actor(i)  + PE_partner(j)
     + prev_pen_mates_shift(i,j)     ← −15 if familiar, 0 if strangers
     + ε_ij
```

### Type 3 — Interaction effects between actor and partner properties

Not yet designed. Example: "do littermates from the same dam fight less?" —
this is a function of both animals' `id_parent_2`. Could be handled by a
derived column in `ind_pair_meta` computed by the user before `add_phenotype()`.

---

## Open Questions for Professor Review

1. **PE persistence across generations** — PE is sampled once per individual per
   phenotype definition. If gen-2 animals are added to the same pen as gen-1, gen-1 PE
   values are reused. Is this biologically correct, or should PE be re-expressed each
   generation?

2. **Singleton handling** — an individual with no group-mates (pen of size 1) produces
   zero dyadic records. Is that a warning, an error, or silently skipped?

3. **Mixed-generation pens** — if a pen contains gen-1 and gen-2 animals and you filter
   to `gen == 1L`, do you want only gen-1 vs. gen-1 pairs, or gen-1 animals paired with
   all pen-mates (including gen-2)? The current design gives gen-1 × gen-1 only.

4. **Undirected dyadic traits** — `directed = FALSE` stores one record per unordered
   pair. Which individual becomes `id_ind_1` vs. `id_ind_2`? Lexicographic order of
   `id_ind`? Does this matter for the PE model if roles are symmetric?

5. **Export format** — downstream BLUP software (e.g. BLUPF90) expects a specific flat
   file structure for dyadic records. Does this table structure (separate
   `dyadic_phenotype`) make the export function straightforward, or is a combined
   view more useful?

6. **Threshold/categorical dyadic traits** — e.g. "did A attack B? (yes/no)" would be
   a binary dyadic phenotype. Does the liability threshold model apply here the same way
   it does for individual phenotypes?

7. **Random pen effect** — the `group_name` column in `dyadic_phenotype` enables a
   pen-level random effect (variance shared within a barn that is not genetics or PE).
   Should we provide a convenience path for this, or leave it as user-driven via
   `define_effect_random()` targeting `dyadic_phenotype.group_name`?

---

## What Does NOT Change

- `define_group_*()` functions write group assignment columns to `ind_meta`; `define_phenotype()`
  just references the column by name. No coupling.
- `add_tbv()`, `add_ebv()`, index machinery — untouched.
- All individual (`class = "individual"`) phenotype paths — fully backward compatible.
- The `define_phenotype()` entry point is retained; `class = "dyadic"` is a new argument,
  not a new function.
- `define_trait()` and `define_additive_effects()` — unchanged (the genetic layer already
  handles correlated `attack_direct` / `attack_received` effects via the G matrix).

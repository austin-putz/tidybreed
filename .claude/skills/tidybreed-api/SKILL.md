---
name: tidybreed-api
description: Reference for every implemented tidybreed function — arguments, semantics and internals of open_pop/define_genome, add_founders, define_chromosome, get_table subsets, schema(), define_trait/define_additive_effects/define_genome_effects, define_genome_effect_terms, define_phenotype, add_phenotype stages, add_tgv and the evaluator, extract_genetic_variance, index functions. Load before changing or explaining any exported function.
---

# tidybreed Implemented Functions

### `open_pop()` + `define_genome()` (genome setup)

`R/open_pop.R`, `R/define_genome.R`

The current surface for creating a population and its genome is
`open_pop() |> define_genome(...)`. `define_genome()` populates the genome tables:

- Genome: `genome_meta` (physical `pos_bp`, `locus_id` `PRIMARY KEY`), `genome_map` (default map), `ind_haplotype` (empty), `ind_genotype` (empty), `ind_crossover` (empty), `chr_inheritance` + `chr_recombination` (default autosome rows)
- Effects: `genome_effects`, `genome_effect_members`, `genome_effect_member_origins` (all empty), plus the `genome_effect_terms` and `genome_effect_loci` views. **These are created here, not in `open_pop()`** — `genome_effect_members` declares a foreign key to `genome_meta.locus_id`, and DuckDB refuses a foreign key to a column that does not exist yet

`define_genome()` key params: `pop`, `n_loci`, `n_chr`, `chr_len_Mb` (finite,
strictly positive), `cM_per_Mb` (genetic-map rate, cM per Mb; scalar or
length-`n_chr`, finite, strictly positive, default `1.0` →
`pos_cM = pos_bp/1e6 * cM_per_Mb`), `locus_names`, `chr_names`,
`recombines_M`/`recombines_F` (genome-wide per-parent-sex recombination defaults,
both `TRUE`; set one `FALSE` for a whole-genome achiasmatic sex, seeded into
`chr_recombination`). Calling
`define_genome()` on a population where **any** of those ten tables or two views
already exists is a hard error (no partial re-definition).

### `add_founders()`

`R/add_founders.R`

Samples haplotypes for each founder individual from the `founder_haplotypes`
pool. Appends rows to `ind_meta` (core 6 cols, including `ploidy`) and
`ind_haplotype` (long: one row per (individual x haplotype x locus), row count
per chromosome driven by the resolved `chr_inheritance`
`from_parent_1`/`from_parent_2` for the founder's sex — 2 rows/locus for a plain
autosome (`1, 1`), 1 for a hemizygous sex chromosome (e.g. `0, 1`), 0 for an
absent chromosome (`0, 0`); `line_origin` = the founder's line, `strand = 1`).
Does **not** write
`ind_genotype` (on-demand via `add_dosage()`). ID format: `{line_name}_{n}`
(e.g. `Libra_1`).

Accepts `...` for custom `ind_meta` columns written atomically with the new
rows (see **Custom field forwarding** below).

Key params: `n_males`, `n_females`, `line_name`, `ploidy` (must be `2` in this
version), then `...` for custom fields.

### Custom field forwarding in `add_*` functions

`add_founders()`, `add_phenotype()`, `add_tgv()`, and `add_ebv()` all accept
`...` for optional custom columns written to the target table at the same time
the new rows are inserted. This avoids a redundant second `mutate_table()` step.

```r
# Single call — gen and farm written with the founders
pop <- pop |>
  add_founders(n_males = 10, n_females = 100, line_name = "A",
               gen = 0L, farm = "Iowa")

# add_phenotype: scalar only (broadcast to all phenotype records)
pop <- pop |>
  get_table("ind_meta") |>
  add_phenotype("ADG", test_env = "barn_A")
```

**Argument disambiguation**: R's standard matching routes explicit formal params
(`n_males`, `n_females`, `line_name`, etc.) to their positions; anything else
falls into `...` and is treated as a custom column. Reserved column names
(`id_ind`, `sex`, `line_name`, etc.) are blocked with an error.

**Type safety** — column types are inferred from the R value via
`infer_duckdb_type()`. Use R's type suffixes to get the right DuckDB type:

| R value | DuckDB type |
|---------|-------------|
| `0L`, `NA_integer_` | `INTEGER` |
| `0`, `0.0`, `NA_real_` | `DOUBLE` |
| `TRUE`/`FALSE`, `NA` (bare) | `BOOLEAN` |
| `"text"`, `NA_character_` | `VARCHAR` |
| `as.Date(...)` | `DATE` |
| `Sys.time()`, `as.POSIXct(...)` | `TIMESTAMP` |

Common pitfall: `gen = 0` gives DOUBLE, not INTEGER. Use `gen = 0L`.

**Pre-declaring a column schema** before data exists (typed-NA workflow):

```r
# 1. After open_pop() |> define_genome(), declare column types on the empty table
pop <- pop |>
  get_table("ind_meta") |>
  mutate_table(gen = NA_integer_, farm = NA_character_)

# 2. add_founders fills values in; types already match
pop <- pop |>
  add_founders(n_males = 10, n_females = 100, line_name = "A",
               gen = 0L, farm = "Iowa")
```

The shared internal helper `prepare_extra_cols()` in `R/sql_utils.R` handles
validation, type inference, ALTER TABLE, and scalar broadcast for all `add_*`
functions.

**Scalar vs. vector**:
- `add_founders()` / `add_offspring()`: scalars broadcast; vectors must have
  length `n_males + n_females` (or `n_offspring`).
- `add_phenotype()` / `add_tgv()` / `add_ebv()`: scalar only (row count per
  trait varies). Use `mutate_table()` afterwards for per-record vectors.

### `define_chromosome()`

`R/define_chromosome.R`

Sets a non-default rule for one chromosome — sex chromosomes (X/Y, Z/W, X0/Z0)
and organelles (MT, plastids). Call before `add_founders()` for the chromosomes
it configures. **Each call sets exactly one concern:** supply `from_parent_1` +
`from_parent_2` to write a `chr_inheritance` row (keyed by `offspring_sex`), or
supply `recombines` to write a `chr_recombination` row (keyed by `parent_sex`) —
never both in one call (mixing them would make "which sex" ambiguous). Writes are
transactional (delete-then-insert with `IS NOT DISTINCT FROM`, then validate, then
commit; rollback on failure). `overwrite = TRUE` (default) upserts; `overwrite =
FALSE` errors if the exact NULL-safe key already exists.

```r
pop <- pop |>
  # Mammal X/Y — inheritance (override only the deviating sexes)
  define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1) |>
  define_chromosome("Y", offspring_sex = "M", from_parent_1 = 1, from_parent_2 = 0) |>
  define_chromosome("Y", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 0) |>
  define_chromosome("Y", recombines = FALSE)   # recombination (both parent sexes)
```

`add_founders()` and `add_offspring()` resolve `chr_inheritance`/
`chr_recombination` per chromosome: a chromosome whose resolved inheritance is
`1,1` for both offspring sexes **and** `recombines` for both parent sexes goes
through the original, unchanged diploid path; any other chromosome routes through
a separate branch that writes only the applicable `(sex, parent_origin)` rows and
— for non-recombining or single-copy inheritance (Y, W, MT) — passes the parent's
stored copy straight through instead of simulating a crossover. Copy resolution
uses the **offspring's** sex/line; recombination resolution uses the **producing
parent's** sex/line.

Real polyploidy (ploidy > 2, uneven-ploidy crosses) is not yet supported;
`ind_meta.ploidy` must be `2` for every individual.

### `mutate_table()`

`R/mutate_table.R`

Generic function. Adds or updates columns in **any** database table. Chain after
`get_table()` (and optionally `filter()`) then call `mutate_table(col = value)`.
Scalar values are broadcast to all (or filtered) rows; vectors must match the
effective row count. Type is inferred via `infer_duckdb_type()`. Reserved
columns are blocked via `TABLE_RESERVED_COLS` in `R/sql_utils.R`. Returns `pop`
invisibly.

When called on an **empty table**, `mutate_table()` still creates the column
schema via `ALTER TABLE ADD COLUMN` (no rows are updated). This is the mechanism
for pre-declaring typed column schemas before data arrives.

```r
# All rows
pop <- pop |> get_table("ind_meta") |> mutate_table(gen = 1L)

# Filtered rows (unmatched rows get NULL for new columns)
pop <- pop |>
  get_table("ind_meta") |>
  filter(sex == "M") |>
  mutate_table(gen = 2L)

# Pre-declare schema on empty table
pop <- pop |>
  get_table("ind_ebv") |>
  mutate_table(model_version = NA_character_)
```

### `define_chip()`

`R/define_chip.R`

Marks the loci in a filtered `genome_meta` table as members of a named chip,
writing a `BOOLEAN` column `is_{chip_name}` to `genome_meta`. All other loci
receive `FALSE`. Accepts a `tidybreed_table` from `get_table("genome_meta")`
(optionally filtered) as its first argument; returns `tidybreed_pop`.

```r
# Filter first, then define chip
pop |> get_table("genome_meta") |> filter(chr %in% 1:5) |> define_chip("chr1to5")

# Random chip: sample locus names, then filter
sel <- pop |> get_table("genome_meta") |> collect() |>
  slice_sample(n = 500) |> pull(locus_name)
pop |> get_table("genome_meta") |> filter(locus_name %in% sel) |> define_chip("50K")

# Complement of an existing chip
pop |> get_table("genome_meta") |> filter(is_50K == FALSE) |> define_chip("non50K")
```

### `get_table()` / `close_pop()` / `print.tidybreed_pop()`

`R/tidybreed_pop.R`

`get_table(pop, "table_name")` returns a `tidybreed_table` S3 object that
carries the pop reference, table name, lazy dplyr tbl, and pending filter.
Supports `filter()`, `collect()`, `select()`, `arrange()`, `pull()`, `count()`,
and `mutate_table()`. `close_pop()` safely closes the DuckDB connection.

**Subset selection for action functions** (`add_phenotype`,
`add_tgv`, `add_ebv`, `add_dosage`, `add_genotypes`, `extract_genotypes`)
requires `get_table()` as the first step. `filter()` is called on the
`tidybreed_table`, not on the pop directly. **The individuals acted on are
the distinct `id_ind` values present in the (filtered) table**, whatever the
table: `ind_meta` for dates/generation, `ind_genotype`/`ind_haplotype` for
marker-assisted pre-selection, `ind_ebv`/`ind_index` for EBV- or index-based
selection, `ind_phenotype` for prior records. An unfiltered `ind_meta` is
everyone; an unfiltered `ind_ebv` is *the animals that have an EBV*, not
everyone. A table without `id_ind` is an error, filtered or not.

All seven go through one internal helper, `resolve_subset_ids(tbl, what,
all_if_null)` in `R/sql_utils.R`: it re-applies the stashed filter quosures
to a fresh un-projected lazy tbl (so `select()` cannot hide `id_ind`), renders
it with `dbplyr::sql_render()`, and runs one `SELECT DISTINCT id_ind ... JOIN
ind_meta` in DuckDB — only the id vector is collected, never the table. It
returns `NULL` for an unfiltered `ind_meta` so callers can keep an unrestricted
SQL fast path (`add_genotypes()`'s `UPDATE` without `WHERE`); callers that need
a concrete vector pass `all_if_null = TRUE`. `add_index()` and `remove_rows()`
deliberately do **not** use it — they act on the rows of the passed table
itself, not on a derived animal set.

```r
# All individuals
pop |> get_table("ind_meta") |> add_phenotype("ADG")

# Filtered by ind_meta
pop |>
  get_table("ind_meta") |>
  dplyr::filter(sex == "F", gen == 1L) |>
  add_phenotype("ADG")

# Marker-assisted pre-selection (add_dosage() first; ind_genotype is a cache)
pop |>
  get_table("ind_genotype") |>
  dplyr::filter(locus_name == "Locus_10", dosage_value == 2L) |>
  add_phenotype("ADG")

# EBV-based: the animals above a threshold in one evaluation
pop |>
  get_table("ind_ebv") |>
  dplyr::filter(trait_name == "ADG", eval_number == 3L, ebv_value > 0.5) |>
  add_phenotype("ADG")

# Pre-select top performers from a prior phenotype
pop |>
  get_table("ind_phenotype") |>
  dplyr::filter(pheno_value > 500) |>
  add_phenotype("ADG2")
```

### `schema()` / `describe_table()`

`R/schema.R`

`schema(pop)` returns a tibble of every table with its display group, row count,
column count and description, and prints it grouped by pipeline stage under
section headings. `describe_table(pop, "name")` drills into one table's columns.
Descriptions live in the `_schema_meta` table and travel with the `.duckdb` file.

```r
schema(pop)                                   # grouped, empty tables collapsed
schema(pop, show_empty = TRUE)                # one row per table
schema(pop, include_system = TRUE)            # also list _schema_meta
schema(pop, order = "rows")                   # flat, biggest first
schema(pop, order = "size", sizes = TRUE)     # on-disk bytes (issues CHECKPOINT)
subset(schema(pop), table_group == "Genome")  # the grouping is data, not text
```

The header reports whole-database size from `PRAGMA database_size` — the file
size plus the WAL when it is uncheckpointed, or memory usage for an in-memory
population. `sizes = TRUE` is opt-in because per-table sizes require a
`CHECKPOINT`, which is a write; the resulting column always prints its caveat
footnote (256 KiB block quantization, and per-table sizes not summing to the
file total).

**Maintenance obligation when adding a table.** Two hard-coded lists in
`R/schema.R` must name every table:

1. `.schema_table_order()` — the display group and in-group workflow position.
2. The matching `.<group>_descriptions()` helper, aggregated by
   `.all_schema_descriptions()` and registered once by `open_pop()`.

A table missing from the first prints under **User tables**; a table missing from
the second prints `(no description)`. Both are deliberately visible failures
rather than silent misfiling, and `tests/testthat/test-schema-print.R` asserts
that the two lists and `SYSTEM_TABLES` name the same tables.

### `define_trait()` / `define_additive_effects()`

`R/define_trait.R`, `R/define_additive_effects.R`

- `define_trait()` — **genetic layer only**. Writes one row to `trait_meta` and
  a global `(index_name = NULL, trait_name, economic_weight = 0)` row to
  `index_meta`. Accepted arguments: `trait_name`, `description`, `units`,
  `overwrite`. It writes **no target** (0.73.0): targets enter only through
  `define_effect_cov_matrix()` or a generator's `G =`. `define_trait_simple()`
  was removed; chain `define_trait()` → `define_additive_effects(G = )` →
  `define_phenotype()`. **Never** pass observation-layer arguments here
  (`type`, `mean`, `expressed_sex`, `residual_var`, etc.) — those belong
  in `define_phenotype()`. `overwrite = FALSE` (default) errors if the trait
  already exists; `overwrite = TRUE` replaces both the `trait_meta` row and its
  `index_meta` entry.
- `define_additive_effects()` — accepts a `tidybreed_table` from
  `get_table("genome_meta")` (optionally filtered) as its **first argument**.
  `trait_name` accepts a scalar **or vector** of trait names. One flow for
  k = 1 and k >= 2 (0.73.0). **Generated means calibrated** (0.74.1, Q21):
  there is no `effects =` or `scale_to_target =`; every call samples and
  calibrates, so a `"generated"` term always delivers its stored target.
  Known coefficients (GWAS estimates, a QTL map, test oracles) go through
  `define_genome_effect_terms()` + `ad_terms()` under a user owner.
  1. **Validate** everything; no RNG use before the draw (`seed` is applied
     after every input check and the anchor-rank check; only the
     architecture-rank and verification failures come after it).
     **Owner content (0.76.0)**: a trait whose `generated` model has anything
     but order-one `additive` terms (a `define_genome_effects()` model) is
     refused after argument validation and before the target
     (`.dae_refuse_nonadditive_model()`); the error gives the additive-only
     `define_genome_effects(..., trait_var_comp_tbl = <additive rows>)`
     re-run, and says later calls need the same filter while the
     non-additive targets stay stored.
     `parent_origin` must be exactly `1`/`2` before coercion.
     `anchor = "realised"` needs an individuals `base_tbl` and the common
     scope. Sex-linked / organelle QTL are refused (`assert_qtl_autosomal()`;
     the error points at the writer).
  2. **Target** (`.dae_resolve_target()`, plan §6C): a passed `G` (matrix, or a
     number for one trait; dimnames checked, never relabelled; PSD) is written
     with the terms and refused over any stored block (whole-table check, even
     identical; the error gives a working `remove_rows()` call and says to
     re-run the same call). A new `G` is also refused when generated terms
     of a line with no target of its own fell back to it (they would keep
     the old one; route: that line's own `G` first), and it warns, naming
     them, for other-`parent_origin` terms at the same target scope (no
     per-parent target exists, so refusing would deadlock;
     `.dae_target_dependents()`). An explicit `trait_var_comp_tbl` picks
     the effect and the traits, never the scope: its block must be the one
     the call's scope reads (`.tvc_resolve_line()`: the line's own block,
     else population-wide), or the call is refused (0.74.5). `G = NULL`
     reads the stored rows (line block, else population-wide), or
     `trait_var_comp_tbl` rows. Refused: a stored block pairing a call trait
     with an outside trait; a stored `dominance` / `additive_by_additive` block
     for the traits (**with or without** a passed `G`; with `G` and generated
     additive terms already present, the error routes through removing the
     non-additive block, since `define_effect_cov_matrix()` would be
     refused; the error also names `define_genome_effects()`); a partial
     block; two candidate sets. The block itself comes from the shared
     per-block resolver `.tvc_block_from_rows()` (both generators). Targets
     are validated by `.qtl_target_std()`: rank and
     PSD on the **correlation scale**, never relative to `G`'s largest
     eigenvalue (that depends on the traits' units). A singular target gets a
     `message()` naming the reason (`.qtl_rank_note()`, 0.76.0; also from
     `define_effect_cov_matrix()` and `define_genome_effects()`).
  3. **Draw** the architecture with `.draw_additive_architecture()` (today's
     draws: rnorm / signed gamma for k = 1, `MASS::mvrnorm(G)` rows for k >= 2,
     masked per trait under `"union"`). The same seed gives the same `B0`.
  4. **Calibrate** with `.dae_calibrate_shared()` → `.qtl_calibrate()` →
     `.qtl_congruence()` (the anchor from `.dae_anchor_at()`, the frames from
     `.dae_build_traits()`: the three helpers the additive-only route of
     `define_genome_effects()` shares, gate C4 (b))
     (`R/qtl_congruence.R`, ported from the source method): `B = B0 A` with
     `B' M B = G`, run on the correlation scale with `B0` columns normalised
     to unit anchor variance, then **verified** against the stored `G` at
     `QTL_CALIBRATION_TOL = 1e-8` (correlation scale); a miss errors before
     any write. "exact" in the message is this check, never assumed. `M` is the
     anchor: `"genic"` = `diag(n_eligible p q)` at the base frequencies;
     `"realised"` = `Cov(X)` of the base individuals (collected in R, size
     guard `QTL_REALISED_MAX_CELLS`). k = 1 reduces to the scalar rescale.
     Two distinct rank errors: rank(G) > rank(M) (checked before the seed),
     rank(B0'MB0) < rank(G). `"union"` (k >= 2) scales each trait alone and
     warns "approximate" whenever the delivered covariance misses `G`
     (including a zero target covariance on overlapping sets); a trait with
     positive target variance and no QTL is an error; one with target 0 gets
     no terms, and its old variant at the scope is deleted like every
     other trait's (0.74.5).
  5. **Diagnostics computed** (common scope only, nothing stored), **before**
     the commit so they cannot fail after it: the relative spectrum of what
     another population sees vs the target — pool expectation `2 Cov(H)` with
     divisor `n_h` (founder base; `add_founders()` samples with replacement;
     identity columns read from the un-projected filtered pool), observed
     `Cov(X)` with `n - 1` (individuals base), or the genic limit (realised
     anchor).
  6. **Commit** terms and any new target in one transaction:
     `.ge_commit(..., before_commit = function(conn) .tvc_write_block(...))`,
     then report the diagnostics. The founder-pool
     comparison is always a `message()` (sampling LD, Q22), adding the
     realised-anchor hint outside `warn_bounds`; the other two warn outside
     `warn_bounds` (default `c(0.8, 1.25)`, `NULL` = off).
  7. **Message** says "exact" / "approximate" with the
     delivered covariance; a line-scoped call adds the line-mean message.

  Writes one order-one `additive` term per locus under the reserved effect
  owner `generated`, with `center_value` = the base allele
  frequency. `line_name` and `parent_origin` compose into the single origin row
  an additive member may carry:

  | `line_name` | `parent_origin` | Stored scope |
  |---|---|---|
  | `NULL` | `NULL` | no origin rows (the common scope) |
  | `"A"` | `NULL` | `('exact', 'A', parent NULL, copy_count = 1)` |
  | `NULL` | `1` / `2` | `('any', NULL, parent, copy_count = 1)` |
  | `"A"` | `1` / `2` | `('exact', 'A', parent, copy_count = 1)` |

  Re-calling replaces **only the variant at the same scope**, so successive
  common / line-A / line-B calls each keep the others — the per-copy fallback
  that makes crossbred breeding values correct needs all of them standing.
  Changing `parent_origin` on a re-run therefore **adds** a variant rather than
  replacing one; that is a legal containment pair and rarely intended, so the
  function warns on exactly that case. `parent_origin` is **per trait** (scalar
  recycled, positional vector, or named by trait); a call mixing origins across
  traits is rejected because one call has one anchor, defined for one set of
  inherited copies (no mixed-scope anchor yet). The covariance is not always
  zero: paternal-only vs maternal-only is zero under random mating, but
  both-parents vs paternal-only is `pq` per locus.

  The calibration is origin-aware:
  `V_A = Σ_j n_eligible,j · p_j q_j a_j²`, `n_eligible` = 2 unparented, 1
  parent-qualified.

  **The base population is `base_tbl`, a filtered `tidybreed_table`** — the
  same two-table shape as `add_ebv(tbl, phenotype = )`: `tbl` says which loci,
  `base_tbl` says which allele copies define `p`. Three shapes, dispatched on
  `table_name` and checked for projected columns: `founder_haplotypes` (the
  pool), `ind_haplotype` (these copies — `filter(line_origin == "Duroc")` is
  Duroc copies at any cross depth), or any table with `id_ind` (these
  individuals, semi-joined on `DISTINCT id_ind`). `p` always comes from
  `extract_allele_freq()`, so a base selection means the same population in
  every writer. `base_tbl = NULL` is the population the effect applies to,
  resolved with the `line → NULL` precedence of `resolve_genome_map()`: the
  line's own founder pool, else the shared (`NULL`) pool, else an error listing
  the pools that exist. Only a population-wide effect on a multi-pool founder
  table warns (Wahlund); an explicit `base_tbl` — including the whole founder
  table, which is how the common fallback variant of a crossbreeding model is
  defined — never warns. A selected QTL with no copies in the base is an error
  naming the loci, never silently centred at `p = 0`.

  ```r
  # Single trait
  pop |> get_table("genome_meta") |> filter(chr %in% 1:5) |> define_additive_effects("ADG")

  # Multiple correlated traits (shared QTL set)
  G <- matrix(c(0.25, 0.10, 0.10, 0.30), 2, 2,
              dimnames = list(c("ADG", "BW"), c("ADG", "BW")))
  pop |> get_table("genome_meta") |> filter(chr %in% 1:5) |>
    define_additive_effects(c("ADG", "BW"), G = G)

  # Generation-0 animals define base allele frequencies
  pop |> get_table("genome_meta") |> filter(...) |>
    define_additive_effects("ADG",
      base_tbl = get_table(pop, "ind_meta") |> filter(gen == 0L))

  # Crossbreeding: common fallback pooled on purpose (no warning), then each
  # line centred on its own founder pool by default
  gm <- pop |> get_table("genome_meta") |> filter(chr %in% 1:5)
  gm |> define_additive_effects("ADG", base_tbl = get_table(pop, "founder_haplotypes"))
  gm |> define_additive_effects("ADG", line_name = "Duroc")
  gm |> define_additive_effects("ADG", line_name = "Landrace")

  # Imprinting: paternal expression, per line rather than trait-wide
  pop |> get_table("genome_meta") |>
    define_additive_effects("IMP", line_name = "Duroc", parent_origin = 1)
  ```

### `define_genome_effects()`

`R/define_genome_effects.R`, calibration in `R/genome_effects_calibration.R`
(0.76.0; plan §9, phase-5 plan 5b)

`define_genome_effects(tbl, trait_name, G_A = NULL, G_D = NULL, G_AA = NULL,
trait_var_comp_tbl = NULL, pairs = NULL, n_pairs = NULL, anchor =
c("genic", "realised"), dominance_degree_mean = 0.19, dominance_degree_sd =
0.097, inbreeding_depression = NULL, base_tbl = NULL, warn_bounds = c(0.8,
1.25))` — samples additive, dominance and A×A effects and **calibrates** them so
`G_A`, `G_D`, `G_AA` are delivered exactly under the anchor; writes under
`"generated"` with `mode = "replace_owner"` (the trait's whole generated model,
line-scoped Part A variants included, counted in a message). Common scope
only; no `seed` (callers use `set.seed()`), no `effects`. Same pipe subject,
base resolution (`.dae_resolve_base()`, Wahlund warning) and size limit as
`define_additive_effects()`.

1. **Validate, no RNG, no write** (C20; the order of the plan's 5b.2):
   arguments → **targets per block** (`.dge_resolve_targets()`, D8: each of
   `additive` / `dominance` / `additive_by_additive` from exactly one source —
   a passed matrix, refused if that block is stored for any call trait
   anywhere at `line_name IS NULL` or also in the explicit table; else the
   explicit `trait_var_comp_tbl` rows or the stored `NULL`-line rows, through
   the shared `.tvc_block_from_rows()` with its scope check) → route (D5:
   **additive-only** when no D / A×A block or only zero ones; otherwise `G_A`
   must be positive definite on its correlation scale, an
   implementation-limit error) → zero-degree refusal (non-zero `G_D` with
   both degree parameters 0) → owner counts → loci, base, autosomal check,
   `n` by `COUNT` → pairs (supplied: unknown / outside / self / repeated
   refused by name; random: `m >= 2`, `n_pairs <= floor(m/2)`;
   `rank(G_AA) <= r`; realised `rank(G_c) <= n - 1`) → realised dosage guard
   counting kept designs (`.dge_dosage_guard()`: `n (m + m_D + r)`, `m_D`
   and `r` only for non-zero blocks — only those get a design
   (`.na_anchors(dominance =)`, A×A anchor only when live), so the count is
   exact; Part A's `n m` on the additive-only route) → anchor ranks per
   block.
2. **Draw** (`.dge_draw()`, D2 order): `B_a` by
   `.draw_additive_architecture(mask, "normal", G_A)` (C4); with a dominance
   block `z` (m×k standard normal; degrees `mean + sd z`); with an A×A block
   and no supplied pairs one `sample()` matching, canonicalised
   (`.dge_canonical_pairs()`: C-locale order within a pair, `locus_id` across
   pairs); then `B_aa`. Zero blocks still draw.
3. **Calibrate.** Additive-only route: `.dae_calibrate_shared()` +
   `.dae_build_traits()` — Part A's path, row-identical (C4 (b)). Otherwise
   `.na_calibrate()` (pure): A×A by `.qtl_calibrate()`; dominance
   `(mean + sd z)|B_a|`, with each trait named in `inbreeding_depression`
   given its solved mean (`.na_solve_dd_mean()`, robust: solved in
   `x = mu/sd`, degenerate / linear / quadratic branches, vertex at a
   repeated root, positive `V_D` required, verified), then `.qtl_calibrate()`;
   additive by `.na_additive_stage()` on the correlation scale of `G_A`: the
   floor `C_res' M C_res` from the residual coupling (PSD by construction),
   rounding budget `1e-10 max(1, ||floor_s||)`, `T` maximising `tr(T)`. The
   coupling `C = b B_d + E_c` is one `rowsum()` (hubs exact). Anchors are
   objects (`.na_anchors()`, `.na_aa_anchor()`; `cross()` added to both anchor
   kinds, design rank cached). A zero D or A×A block gets zero coefficients
   without a calibration or anchor. Every present block is verified at
   `QTL_CALIBRATION_TOL`. No eligibility zeroing (D3 (a)). Then
   `.na_store_alpha()`: the stored (HWE-at-`p`) additive coefficient is
   `B_alpha` itself under genic and `B_alpha + Delta` (coupling of the
   HWE-minus-observed `b`, `c`) under realised — never `(B_alpha - C) + C`,
   which cancels when `G_A` ≪ `G_D` / `G_AA` — re-verified as stored; a miss
   is refused. `.dge_build_nonadditive()` writes these in place of
   `.noia_to_stored()`'s recovered `alpha`.
4. **Diagnostics** (D6, before the commit): `.dae_diagnostics()` on the
   additive-only route; otherwise `.dge_diagnostics()` — realised anchor: each
   block's genic limit; genic + individuals: each block realised on them (A×A
   via `.egv_aa_values()` chunks); genic + founder pool: additive pool
   expectation only, a message.
5. **Store**: additive-only route Part A's frames; otherwise per trait
   `.noia_to_stored(a, d, pairs(e), p_base)` → `.noia_terms()` (Cockerham
   `ad_terms()` + `aa_terms()` at the base `p`, exact zeros dropped). Under
   `"realised"` the stored split is HWE-referenced; the total is exact and the
   realised extractor gives the targets back. One `.ge_commit()` with
   `before_commit` writing each **passed** block (`.tvc_write_block()`).
6. **Messages** (`.dge_messages()`): replaced counts, the random-pair default,
   delivered blocks, the floor (conditional on the draw), inbreeding
   depression (requested vs delivered; approximate for k >= 2, D4; or
   implied), the functional summary, diagnostics, rank notes.

Gates: `tests/testthat/test-define_genome_effects.R` (C1–C18, C20, G1–G8) and
`tests/testthat/test-genome-effects-calibration.R` (the source suite ported,
G4, G6), with the source generator copied in
`tests/testthat/helper-nonadd-generator-oracle.R` (isolated environment).
End to end (5c): `tests/testthat/test-define_genome_effects-integration.R`
(C13 phenotypes, prevalence on a generated A + D + A×A model, removal, C19,
C20, with an exact HWE + LE factorial fixture). Benchmark
`dev/benchmarks/benchmark_define_genome_effects.R`.

**Vignette** `vignettes/genetic-models.Rmd` ("Genetic models") walks the four
paths — `define_genome_effect_terms()` with the builders,
`define_additive_effects()`, `define_genome_effects()` (with the floor), and
`extract_genetic_variance()` joined to `trait_var_comp` — and states the scope
promises. Keep it in step when any of the four changes.

### `remove_generated_effects()`

`R/remove_generated_effects.R` (0.74.5)

`remove_generated_effects(pop, trait_name, line_name = NULL, parent_origin =
NULL)` — the **only** route that deletes generated terms other than re-running
their scope: `define_genome_effect_terms()` and `remove_rows()` refuse the
reserved owner. Deletes, per trait, the generated terms of **every kind** (additive,
dominance, interaction) at exactly the scope `.dae_scope(line_name,
parent_origin)` (via `.ge_resolve_deletes(..., "replace_scope")`, which
matches origin rows, not contrasts), in one transaction
validated by `validate_genome_effects()`. Errors, deleting nothing (all or
nothing across traits), when a trait has none at the scope. Targets and
`ind_tgv` are untouched; values are stale until `add_tgv()` re-evaluates. The
typical use is a variant added by a re-run with a new `parent_origin`, which
`.dae_warn_parent_only()` now names. A generated model is removed whole, never
per component (decided 2026-10-05): a `define_genome_effects()` model is
common-scope, so `remove_generated_effects(pop, trait)` removes all of it, after
which `define_additive_effects()` accepts the trait again.

### `extract_allele_freq()`

`R/extract_allele_freq.R`

`extract_allele_freq(tbl)` — the single place a population selection becomes
per-locus allele frequency. Takes the three `base_tbl` shapes above; returns
one row per `genome_meta` locus in `locus_id` order (`locus_id`, `locus_name`,
`allele_freq`), `NA` at a locus the selection has no copies for (never `0`),
an error if no locus is covered at all. One SQL statement with the filter
rendered as a subquery via `dbplyr::sql_render()`; nothing else is collected.
Never warns, never writes. Users call it to obtain `p` for `ad_terms()`. Also
holds `.validate_base_tbl()`, shared by both genome-effect writers.

### `extract_genetic_variance()`

`R/extract_genetic_variance.R` (0.75.1). `extract_genetic_variance(tbl,
trait_name = NULL, base_tbl = NULL, anchor = c("realised", "genic"))` — the
measuring instrument. **Read-only**: values come from `.gev_evaluate()`, which
writes nothing; `ind_tgv` is not touched.

- `tbl` selects the cohort through `resolve_subset_ids(all_if_null = TRUE)`;
  fewer than 2 individuals is an error.
- `trait_name = NULL` is every trait **with stored terms**, in `id_trait` order
  (`.egv_traits()`; deliberately not `.gev_resolve_traits(NULL)`, which returns
  every trait and which `add_tgv()` relies on). A named trait without terms
  errors.
- `base_tbl` is genic only (non-`NULL` with `"realised"` errors). `NULL` is the
  cohort's **whole-genotype** frequencies via a registered id view, whatever
  table selected the cohort (an `ind_haplotype` filter picks individuals, not
  copies). An explicit `base_tbl` uses `extract_allele_freq()` semantics. A
  decomposed locus with no copies in the base errors, naming it.
- **Output**: tibble `effect_name, trait_name_1, trait_name_2, cov_value,
  n_ind, decomposition, anchor`; full square per block, effect order
  `additive, dominance, additive_by_additive, unpartitioned,
  between_components, total`. The A×A name is `additive_by_additive`
  (`trait_var_comp` vocabulary), not the `ind_tgv` component `interaction`.
  A `message()` names the population (decision D1); the tibble carries only
  `anchor`.
- **Classification** (`.egv_classify()`): a term is covered when it is
  order-one `additive` / `dominance` / diploid `indicator` (any state) or
  two-member `additive × additive`, has no origin rows, and every member locus
  is `1, 1` for both sexes in every line (`.egv_diploid_loci()`). A family
  (`family_key`) is covered only as a whole (**family rule**: scope variants
  compete). Cases per trait: `full` (all covered), `additive_only` (every term
  order-one additive, not full: the evaluated `additive` component),
  `partial` (covered part projected, the rest evaluated alone into
  `unpartitioned`, absent individuals 0). Off-diagonals carry the less
  complete label. `"genic"` on a non-full trait errors naming `"realised"`.
- **Computation**: covered terms per trait → `.stored_to_functional()` (owners
  summed after classification) → coefficient matrices on the covered-locus
  index (`.egv_coefficients()`). Realised: the source's
  `nonadd_covariates()` algebra without its dense anchor matrices — cohort
  `p`, observed regression `b` (0 at monomorphic loci), `alpha = a + b d +
  sum e c` (the induced `e·c` of a pair with a fixed member is **kept**), value
  matrices `Z_A alpha`, `Z_D d`, and A×A accumulated in deterministic pair
  chunks (`.egv_aa_values()`, chunk `QTL_REALISED_MAX_CELLS / n`, no `n × r`
  matrix). Covariances divisor `n − 1`. `between_components` = sum of
  `Cov(b_t1, b'_t2)` over every ordered pair of different blocks
  (`unpartitioned` included), computed directly. An internal assertion
  requires the centred blocks to sum to the centred evaluated total (mixed
  tolerance) or errors. Genic: closed forms at the base `p`, `b = q − p`; no
  `between_components`.
- **Block availability** follows the canonical model, never stored contrast
  names: `additive` for every trait with a covered term; `dominance` iff some
  canonical `d ≠ 0`; `additive_by_additive` iff some canonical `e ≠ 0`; an
  off-diagonal block row only when both traits have the block. A supported
  block of variance 0 is reported as 0. `.stored_to_functional()` sets a `d`
  or `e` to exactly 0 when it is cancellation residue, `|x| ≤ n·eps·Σ|contrib|`
  (`.cancelled()`), so a linear surface in any row order has no dominance row;
  a single small contribution is kept.
- **Meaning under LD** (Codex implementation review finding 1): the realised
  blocks are covariances of the NOIA *contrast components* (each locus'
  heterozygosity regressed on its own dosage only), not the cohort's joint
  least-squares additive projection. Under LD `additive` can differ from
  `var(fitted(lm(g ~ dosages)))` even with `between_components = 0`; `full`
  means supported term shapes, not recovered breeding-value variance. Any
  breeding-value export must state which of the two it is.
- **Resources**: the `n × m` dosage guard (`.dosage_guard()`) runs before any
  evaluation; dosages come from the shared `.collect_dosages()` in
  `R/genome_effects_helpers.R` (also behind `define_additive_effects(anchor =
  "realised")`, each caller with its own message pieces).
- Gates B1–B18 in `tests/testthat/test-extract_genetic_variance.R`; the source
  oracle is copied verbatim in `tests/testthat/helper-nonadd-oracle.R`.

**How the two writers relate.** `define_genome_effect_terms()` writes any effect
you supply; `define_*_effects()` functions sample effects of one shape and
write them through the same engine (`.ge_build → .ge_read_model →
.ge_resolve_deletes → .ge_commit`); `define_genome_effects()` builds its
non-additive rows with `.noia_terms()`. `define_additive_effects()` is provably
sugar over the writer — `tests/testthat/test-genome-effects-writer.R`
("generator == writer") reproduces its output exactly through
`define_genome_effect_terms()` with the reserved owner, `replace_scope`, and the
same `base_tbl`. `base_tbl = NULL` deliberately differs: the generator has a
domain default; the writer fills nothing, so a missing centre is an error.

### `define_genome_effect_terms()` / `ad_terms()` / `genotype_terms()`

`R/define_genome_effect_terms.R`, `R/genome_effect_terms_builders.R`

`define_genome_effect_terms(pop, trait_name, terms, effect_owner = "custom", mode =
c("append", "replace_scope", "replace_owner", "replace_trait"), origin = NULL,
base_tbl = NULL, require_complete = FALSE)` — the
general writer
for arbitrary genome effects. `terms` is a **long data frame, one row per
(term × locus)**: `term_id` (user-facing only; never stored), `genome_value`,
`effect_name`, `locus_name`, `contrast_name`, `center_value`,
`copy_count_value`, `dosage_value`. A single-term call may omit `term_id`.
Scope lives in a separate `origin` argument — `NULL` (the common scope), a
named scalar list applied to every member, or a data frame keyed by
`locus_name` for the exact multisets a genotype member takes. The reserved owner `generated` is refused
with **no exported override** (Q23): only the internal engine
`.ge_write_terms(..., allow_reserved_owner = TRUE)` may write it, so "every
active term is `generated`" proves calibration. `replace_trait` refuses to
delete generated terms; they are replaced only by re-running the generator.

The writer resolves `locus_name` → `locus_id`, canonicalizes members by
ascending `locus_id` and origin rows by the sorted tuple, infers
`copy_count_value` for indicator input at diploid-autosomal loci, assigns ids
via `next_int_id()`, and writes in **one transaction** that validates the whole
table set before `COMMIT`. Every message about malformed input names the
`term_id` the user typed, never an `id_genome_effect` they have not seen.
`require_complete = TRUE` demands every reachable `(copy_count, dosage)` state
on every member of an indicator surface — including `copy_count_value = 0`
where a chromosome can be absent. `base_tbl` (a filtered `tidybreed_table`;
see `extract_allele_freq()`) fills Cockerham `center_value` on any `additive`
or `dominance` member that has none — column omitted or `NA`. An explicit
centre always wins, `indicator` members are never touched, the base is
validated whenever supplied but queried only if some centre is missing, and
without `base_tbl` a missing centre is an error. The fill happens inside
`.ge_build()` before member validation, whose per-row message says when the
base had no copies at that locus. One `base_tbl` gives one `p` per locus per
call, so line-scoped surfaces are written one line at a time with
`mode = "append"`.

Three builders produce `terms`, because a surface is rows, not a second
representation. **All three return one fixed column set** (0.75.0, gate C16):
`term_id, locus_name, contrast_name, center_value, copy_count_value,
dosage_value, genome_value, effect_name` — character, character, character,
double, integer, integer, double, character — with a typed `NA` where a column
does not apply (built by the internal `.terms_frame()`), so their outputs
`rbind()` in any combination. **`term_id`s are collision-free** (`.term_id()`):
a builder prefix plus length-prefixed locus names plus a suffix, e.g.
`"ad:2:L1#a"`, `"aa:1:A|3:BxC"`, `"geno:1:A|1:B#4"`. Joining names with a
delimiter is not enough — `c("A","B")` and the single locus `"AxB"` would
share an id and the writer would merge two surfaces into one product term.
Bound outputs share a `term_id` only when they describe the same term on the
same loci; overlapping definitions are then refused by the writer's family
rule, never merged.

- `ad_terms(locus_name, a, d, p, coding = c("functional", "cockerham"))` —
  expands an (a, d) pair. Functional coding is `additive`@`0.5` plus
  `indicator`@`(2, 1)`; Cockerham is `additive`@`p` plus `dominance`@`p`. It
  **reports** the implied genetic mean `μ = a(p − q) + 2pq·d` and writes it
  nowhere — putting it in `phenotype_meta.mean` would double-count once
  non-additive genetic values reach the phenotype layer.
- `aa_terms(locus_name_1, locus_name_2, e, p_1, p_2, coding, effect_name,
  report)` — two `additive` members per pair, one coefficient; centres 0.5
  (functional) or `p_1`, `p_2` (Cockerham). Pair order is canonicalised
  (C-locale radix; each `p` moves with its locus); a self-pair or a repeated
  pair is refused; inputs are validated **before** `e = 0` pairs are dropped;
  all-zero `e` is an error. Functional coding reports each pair's share of μ,
  `e(2p_1 − 1)(2p_2 − 1)`. Cockerham `ad_terms()` takes α, and functional `a`
  is not α once pairs exist — feeding it changes genotypic values, not just the
  component split.
- `genotype_terms(genotypes, value, copy_count = NULL, drop_zero = TRUE)` —
  turns a genotype-by-value table into `indicator` terms, one term per row and
  one member per locus column. Dosages and `copy_count` must be non-negative
  whole numbers, checked before `as.integer()` and before zero rows drop
  (a `2.9` would otherwise truncate into a different, valid state).

**The NOIA conversion pair** (internal, same file; plan Q13):
`.stored_to_functional(terms, members)` converts covered common-scope terms
(order-one `additive` / `dominance` / diploid `indicator` in any of the three
states, two-member `additive × additive`) to functional `(a, d, e)` keyed by
`locus_id`, plus `kappa` (**stored = functional + kappa**). `.noia_to_stored(a,
d, pairs, p)` is the forward map (`α = a + (q − p)d + Σ e(2p_l − 1)`, `μ` the
functional HWE/LE mean; **functional = statistical + μ**, so a round trip gives
`kappa = −μ`); it returns coefficient data, and `.noia_terms()` turns it into a
writer frame. The input is sparse: a locus missing from `a` or `d` has that
coefficient 0, and the result spans every locus the model names, pair-only
loci included (a pair induces α at both loci). `.noia_terms()` refuses
`alpha`/`d`/`p` that are not named alike or miss a pair locus. `ad_terms()` is not rewired through them; its no-pair μ is
asserted equal. Gates N1–N3 in `tests/testthat/test-genome-effect-terms-builders.R`
check every conversion row against the real evaluator.

```r
# One dominance term, Cockerham coding at p = 0.3
pop |> define_genome_effect_terms("ADG", data.frame(
  locus_name = "Locus_10", contrast_name = "dominance",
  center_value = 0.3, genome_value = 0.8))

# A 3x3 A x A surface: nine cells, nine terms, two members each
cells <- expand.grid(Locus_10 = 0:2, Locus_44 = 0:2)
pop |> define_genome_effect_terms(
  "ADG", genotype_terms(cells, c(0, 0, 0, 0, 1.4, 2.1, 0, 2.1, 3.6)),
  effect_owner = "epistasis_AxA")

# Reciprocal dominance: the F1 value depends on which parent gave which line
pop |> define_genome_effect_terms(
  "ADG",
  terms  = data.frame(term_id = 1L, locus_name = "Locus_10",
                      contrast_name = "dominance", center_value = 0.3,
                      genome_value = 1.2),
  origin = data.frame(term_id = 1L, locus_name = "Locus_10",
                      line_match_type = "exact",
                      line_name     = c("Duroc", "Landrace"),
                      parent_origin = c(1L, 2L), copy_count = c(1L, 1L)),
  effect_owner = "reciprocal")
```

### `define_effect_cov_matrix()` / `define_effect_random()` / `define_effect_fixed_class()` / `define_effect_fixed_cov()` / `define_effect_intercept()`

`R/define_effect_cov_matrix.R`, `R/define_effect_random.R`, `R/define_effect_fixed_class.R`, `R/define_effect_fixed_cov.R`, `R/define_effect_intercept.R`

- `define_effect_cov_matrix(pop, effect_name, cov_matrix, trait_name = NULL, line_name = NULL)` — **single entry
  point for all variance/covariance data**. `cov_matrix` may be a number when
  one name is given; a named matrix must match `trait_name` in order (never
  relabelled). Routes by `effect_name`:
  genetic effects (`GENETIC_EFFECT_NAMES`: `"additive"`, `"dominance"`,
  `"additive_by_additive"`) → `trait_var_comp` through `.tvc_write_block()`
  (PSD, `%.17g` full precision, one transaction, **never overwrites** a stored
  block for the same `effect_name` × any trait × `line_name`; `line_name` is
  genetic-only). The exported function also **refuses a genetic block when a
  trait already has `"generated"` terms of that kind at the block's scope**
  (`.tvc_refuse_under_generated()`, 0.74.1, Q21): line-C terms for `"C"`;
  for `line_name = NULL`, population-wide terms **and** the terms of any line
  with no block of its own (they fell back to the population-wide target); still refused after the old
  block is removed. The route is `remove_rows()` then
  `define_additive_effects(G = )`, which writes target and terms together. A
  line's target before that line's effects is accepted; a generator's `G =`
  goes through `.tvc_write_block()` and is never refused here;
  `"residual"` → `define_residual_cov()` → `phenotype_var_comp`;
  any other name → `phenotype_var_comp` with that `effect_name`.
  `GENETIC_EFFECT_NAMES_FUTURE` (`additive_by_dominance`,
  `dominance_by_dominance`) → "not yet supported"; `DERIVED_EFFECT_NAMES`
  (`total`, `unpartitioned`, `between_components`) → refused. The phenotype
  layer (`define_effect_random()`, `define_effect_fixed_*()`,
  `write_phenotype_cov_block()`) refuses all three sets as effect names
  (`.check_effect_name_input()`).
  Readers `get_trait_var()` / `load_trait_cov()` take `line_name = NULL`
  (population-wide rows; a named line uses its own rows when it has any,
  else falls back; lines are never mixed).
  Can be called before `define_trait()` or `define_effect_random()`.
- `define_effect_random()` — `variance = NULL` (default) requires a value
  already in `phenotype_var_comp`; a number writes a 1 × 1 block and is an error
  when the phenotype is already in a multi-phenotype block for that effect
  (redeclare it with `define_effect_cov_matrix()`). Once in a block of two or
  more, the row must be `distribution = "normal"` and share the block's
  `(source_column, source_table)`. `overwrite = TRUE` discards that phenotype's
  stored draws for the effect. One transaction. A level's draw is persistent
  (see `phenotype_random_effects`): an effect that should be re-realized per
  batch needs the batch in the level (`pen_batch`), not a new feature.
- `define_effect_fixed_class()` — discrete level → shift mapping.
- `define_effect_fixed_cov()` — linear regression term (`slope * (x - center)`).
- `define_effect_intercept()` — sets the phenotype intercept (`phenotype_meta.mean`).

### `define_phenotype()` / `define_residual_cov()`

`R/define_phenotype.R`, `R/define_residual_cov.R`

- `define_phenotype(pop, phenotype_name, type, mean, expressed_sex, repeatable, ...)` —
  registers an observed phenotype in `phenotype_meta`. For simple traits
  `phenotype_name` matches the `trait_name` already in `trait_meta`. For
  composite phenotypes (e.g. weaning weight, SGE ADG) the name is new and no
  prior `define_trait()` call is needed for it.

  Key arguments:
  - `residual_var` — scalar; writes a 1 × 1 unconditional residual block to
    `phenotype_var_comp` (with `effect_name = 'residual'`). Error if the
    phenotype is already in a multi-phenotype residual block (redeclare it with
    `define_residual_cov()`) or if its residual has realized draws.
    `overwrite = TRUE` without `residual_var` leaves `phenotype_var_comp`
    untouched. Every defined member of the phenotype's residual block must
    share its `condition_change_action`; checked before anything is written.
  - `components` — data frame with columns `source_trait_name` and
    `contributor_type` (`"self"`, `"dam"`, `"sire"`, `"group"`). Optional
    columns: `weight`, `weight_type`, `aggregation`, `group_column`,
    `group_table`, `covariate_name`, etc. Writes to `phenotype_components`.
    `NULL` (default) = simple single-self trait. `component_names` (default
    `"total"`) is which genetic value a row reads: `"total"` or a
    comma-separated list of `TGV_COMPONENT_NAMES`, validated (with every other
    component check) before anything is written. A `group` contributor's
    mate sum accumulates exactly (`GEV_ACC_TYPE`, `.group_mate_tgv()`), so
    `group_sum()` / `group_mean()` are bit-identical across thread counts.
  - `formula_tgv` — DSL shorthand for a composite genetic value (Q18). A bare trait symbol is `self(trait)`;
    `self(trait)` / `dam(trait)` / `sire(trait)` take one positional trait,
    `group_sum(trait, col)` / `group_mean(trait, col)` two. Every call takes a
    named-only `component =` (`TGV_COMPONENT_NAMES` or `"total"`, default
    `"total"`); the group calls a named-only `table =` (default
    `"ind_meta"`). `.walk_formula_tgv_ast()` validates and substitutes in one
    pass (each reference becomes a `.tgv_<n>` placeholder, so two references
    differing in any argument stay distinct); `.validate_formula_tgv()` also
    checks the traits and that each group table exists with `id_ind` and the
    column. Any other call (only `+ - * / ^`, parentheses, numbers and the
    math whitelist are allowed — the expression is `eval()`ed), an unknown
    or extra positional argument, or a non-identifier `col` / `table` is an
    error in `define_phenotype()`, before anything is written. Numbers
    (weights, offsets) are accepted silently. Table and column names match
    exactly; a case-only mismatch gets a "did you mean" hint (`.case_hint()`).
  - `prevalence` (categorical, two categories) — the fraction strictly
    above the threshold `mean + qnorm(1 - prevalence) * sqrt(Vg + Vr + Ve)`,
    a Gaussian approximation at an HWE/LE reference with orthogonal
    components (fixed effects are not in it). `Vg` is the active-block sum
    (`.ap_prevalence_genetic_var()`): the trait's stored population-wide
    `additive`, `dominance` and `additive_by_additive` diagonals, each counted
    only if the model has terms of that kind, and **every term must be owned
    by `"generated"`** (the owner rule, 0.74.1, Q21: only generator terms are
    known to deliver the stored target). `Vr` is every named random effect's
    stored variance (`.ap_prevalence_env_var()`, 0.74.5; `normal` and
    `uniform`); `Ve` the unconditional residual. Refused with `components` /
    `formula_tgv` (no stored variance describes a composite liability: use
    `thresholds`). `add_phenotype()` errors in PLAN (`.ap_check_prevalence()`,
    before any write or draw) when any term is not `"generated"`, a kind of
    term has no stored target, the model has terms outside the three kinds,
    one kind at one line has generated variants for two parent scopes
    (`.gev_term_parent()`; each was calibrated alone), a random effect is
    `gamma`, or the residual has conditional strata; and in RESOLVE when the
    total is 0. There is no silent `Vg = 0`.
    Skipped for `user_values` calls, which place no threshold.
  - `thresholds` — finite, strictly ascending; validated before any write. A
    liability exactly on a cutpoint stays in the lower category
    (`liability_to_categorical()`, `findInterval(left.open = TRUE)`).
  - `formula_tgv` evaluation (`.eval_formula_tgv()`): a constant is
    broadcast; an `Inf` / `NaN` the arithmetic produced (not an `NA` from a
    missing contributor) is an error in PLAN naming the individuals.
  - **Atomic.** `mean` (one finite number) and `thresholds` are checked first;
    the delete of an overwritten definition, the `phenotype_meta` insert, the
    residual (`.pvc_write_block()`, transaction-free) and the components are
    one transaction.
  - `missing_component_action` — `"skip"` (default) or `"error"`. Stored in
    `phenotype_meta` and applied uniformly by `add_phenotype()` for **any**
    missing composite piece (missing group assignment, missing dam/sire genetic value,
    etc.). `"skip"` excludes the individual and warns with a count + up to 5
    example IDs. `"error"` stops immediately.

- `define_residual_cov(pop, phenotype_names, cov_matrix, condition_column = NULL, ...)` —
  writes one stratum of a residual covariance block to `phenotype_var_comp`
  (always with `effect_name = 'residual'`). Supply a named matrix for
  multi-phenotype correlated residuals, or call once per sex/group level with
  `condition_column = "sex"` and `condition_level = "M"` / `"F"` (both together)
  for heterogeneous residuals. The block rules under `phenotype_var_comp` apply:
  whole block per call, one condition column, same phenotypes in every stratum,
  locked once realized. Rejected calls change nothing.

- `define_condition_change_action(pop, phenotype_name, condition_change_action)` —
  sets `phenotype_meta.condition_change_action` on **every member** of the
  named phenotype's residual covariance block, in one transaction. D6 requires
  the members to agree, so once a block has two or more there is no ordering of
  `define_phenotype()` calls that changes it — each single-member flip is the
  disagreeing state D6 refuses. This writer is the way (D6 mutability). It
  touches one column and nothing else, so unlike
  `define_phenotype(overwrite = TRUE)` it cannot reset the rest of the row. It
  is **not** locked by realized draws: D3 locks the covariance *matrix*, while
  the action only governs how future records condition on stored residuals.

### `add_phenotype()` / `add_tgv()`

`R/add_phenotype.R`, `R/add_tgv.R`, `R/genome_effects_eval.R`

Both functions accept a `tidybreed_table` (from `get_table()` + optional
`filter()`) as their first argument and return `tidybreed_pop`.

- `add_phenotype()` — the workhorse. `phenotype_name` (formerly `trait_name`)
  defaults to all phenotypes in `phenotype_meta` when omitted. Runs in
  **three stages** (`R/add_phenotype_stages.R`, `?add_phenotype_stages`):
  1. **PLAN** (`.ap_plan()`, no RNG, no writes except the `add_tgv()`
     prerequisite): sorted subset, metadata, topological sort of derived
     formulas, sex expression, repeatable guard, fixed-effect terms with
     `null_class_action`, the genetic value (simple: the trait's total,
     `ind_tgv_total`; composite via `.assemble_composite_tgv()`, each row
     reading its `component_names`, `"total"` by default; `formula_tgv` via
     the DSL, each reference's `component =`, `"total"` by default) with
     `missing_component_action`, `pheno_number`, the residual condition value
     and the random-effect level of every planned record.
  2. **RESOLVE** (`.ap_resolve()`, RNG, no writes): every draw in a fixed
     order, through two adapters over `find_covariance_blocks()` and
     `resolve_correlated_draws()`. First the **named-effect adapter**
     (`.ap_resolve_named_effects()` → `.ap_named_effect_block()`): effects
     in byte-sorted `effect_name` order, blocks in loader order, entity =
     the level, coordinates = the block's phenotypes; a level draws its
     planned coordinates conditional on the draws it already has stored in
     `phenotype_random_effects` for the block's other members, one resolver
     call per sample-set group; `validate_named_effect_block()` re-runs per
     block as the §5.6 backstop; a 1 × 1 `gamma`/`uniform` block keeps its
     marginal sampler. Then the **residual adapter**
     (`.ap_resolve_residuals()` → `.ap_residual_block()`): one residual
     covariance block at a time, entity = `(id_ind, pheno_number)`; each
     entity draws its planned coordinates from the stratum its condition
     value selects (unconditional `R` as fallback, stored as
     `residual_condition_level = NULL`; error if there is none), conditional
     on the residuals it has already realized — stored on disk for any
     block member at the same `pheno_number`, or fixed by `user_residual` —
     one call per `(stratum, sample set)` group. D6 (agreement) and D2
     (stratum change: error, or drop with a warning under `'independent'`)
     run here on the stored coordinates. Then liability and type
     conversion, all in memory.
  3. **COMMIT** (`.ap_commit()`, writes, no RNG): one transaction, register +
     `INSERT` into `phenotype_random_effects` and `ind_phenotype`; rollback on
     failure.

  Rules that follow: planned ids never appear in SQL text (`.ap_read_by_id()`
  registers them); never `dbWriteTable()` in this path (it advances the RNG);
  sort before every RNG-consuming step; an individual without a record must
  not consume RNG or leave stochastic state. Records are ordered by `id_ind`
  within a phenotype — that is the positional order for `user_values` /
  `user_residual` (a plain vector when one phenotype is model-generated,
  else a named list that may name a subset; the rest are drawn conditional
  on it). `residual_value` / `residual_condition_level` are written for
  every model-path record. See `plans/sample_correlated_effects.md` §5.5.

  **Failure contract (D7)**: the database is atomic, the RNG is not. Any
  error — Stage-1 rejection, Stage-2 error after some draws, failed Stage-3
  write — leaves `ind_phenotype` and `phenotype_random_effects` untouched
  (the RNG-independent `add_tgv()` write is the one write that remains), and
  `.Random.seed` advanced by exactly the draws made before the error.
  Nothing in `R/` touches `.Random.seed`; never add seed restoration to one
  function — if the package ever adopts it, it is a package-wide policy.
  `tests/testthat/test-add_phenotype_failure_contract.R` asserts both halves.
  **Phenotypes see the total genetic value** (0.74.0, plan §6A). Simple
  phenotypes need a trait with at least one term (any kind, any owner). A
  `prevalence` threshold uses the active-block rule
  (`.ap_prevalence_genetic_var()`): the sum of the population-wide stored
  diagonals of `additive`, `dominance` and `additive_by_additive`, each only if
  the model has terms of that kind, plus the named random effects' variances
  and the residual; a term not owned by `"generated"`, an `indicator` surface
  or another interaction, a kind with no stored target, two parent scopes of
  one kind at one line, a `gamma` effect or conditional residual strata is an
  error naming `thresholds =`.
- `add_tgv(tbl, trait_name = NULL, index_names = NULL, weight_type =
  c("index", "economic", "both"), component_name = "additive",
  overwrite_index = FALSE, ...)` — the one table of true genetic values.
  Evaluates **every** term of a trait (every owner) and writes `ind_tgv`, one
  row per (individual × trait × `component_name`): `additive`, `dominance`,
  `indicator` (a one-locus term's contrast) or `interaction` (two or more
  loci) — `TGV_COMPONENT_NAMES`. **The breeding value is `additive`** for
  generated effects (statistical coding at one base `p`); for hand-written
  functional terms it is the functional additive effect (`α = a + d(q − p)`).
  The raw sum of the stored terms; **no mean is added**. Each allele copy
  takes the most specific variant whose origin predicate matches its
  `(line_origin, parent_origin)` label, falling back per copy to the common
  variant (crossbreeding, imprinting). Re-evaluating an (individual, trait)
  deletes components the model no longer produces and upserts the rest on
  `(id_ind, trait_name, component_name)`, so custom columns survive. Total via
  the `ind_tgv_total` view (fixed-order sum, bit-identical at any thread
  count). `...` writes scalar custom columns. True index:
  - `index_names` — named indices; multiplies per-trait genetic values by the
    index weights and writes `ind_true_index`.
  - `weight_type` — `"index"` (default, `index_weight`), `"economic"`, or
    `"both"`.
  - `component_name` — which value is weighted: `"additive"` (default, the
    breeding value), another component, or `"total"`; stored in
    `ind_true_index.component_name`, so additive and total indices coexist.
    `"additive"` is structural, so it **warns** when an index trait has
    `indicator` terms or hand-written interactions, whose additive value it
    misses (`.tgv_warn_structural_additive()`, 0.74.5).
  - `overwrite_index = FALSE` — skips individuals that already have a row for
    `(index_name, weight_type, component_name)`; `TRUE` recomputes.
  Consumers read values through `.tgv_read()` / `.tgv_by_id()` /
  `.group_mate_tgv()` (`components = "total"` or a listed set).

#### The evaluator (`R/genome_effects_eval.R`)

Both functions run one evaluator, built on the fact that an origin predicate
reads a copy's `(line_origin, parent_origin)` **label** and nothing else. The
winning variant is therefore a function of the label, evaluation tuples group by
label-vector, and the inner sum factors inside each group. Three artifacts:

| Artifact | Grain | Built |
|---|---|---|
| Label alphabet | one row per distinct label | one `DISTINCT` per member kind |
| Resolved variant map | `(family, label-vector) → id_genome_effect` | in R, from stored rows only — never per individual |
| Member reduction | one row per `(id_ind, id_genome_effect, member_slot, label)` | SQL |

Containment search runs **only** while building the map; it never runs during
evaluation, and no tie can reach it because overlapping-but-incomparable scopes
are refused at write time. Evaluation is a fixed number of statements whatever
the population size (five for an additive model), and individual identifiers
never appear in the SQL text.

Two shortcuts keep the map small, and both are the plan's fast path rather than
special cases: a family no variant scopes resolves to itself for every
label-vector, and a **member** no variant scopes carries the sentinel label
`"*"`, reducing over every unit at once. Without the second, a 50-locus
unscoped dominance term would enumerate `|labels|^50` label-vectors. The
resolution of a family is cached on a signature that excludes `locus_id`, so a
500-QTL model with common/line-A/line-B variants solves one problem, not 500.

`options(tidybreed.label_vector_warn)` (default `1e4`) and
`options(tidybreed.label_vector_max)` (default `1e6`) bound the map.

Reference implementations live in `tests/testthat/helper-genome-effects.R`
(Phase A: two independent evaluators, hand-computed fixtures, no database);
`tests/testthat/test-genome-effects-eval.R` asserts the SQL evaluator agrees
with them for every fixture and every individual.

### `define_index()` / `add_index()`

`R/define_index.R`, `R/add_index.R`

- `define_index(pop, index_name, trait_names, index_wts, economic_wts = NULL, overwrite = FALSE, ...)` —
  registers a named selection index in `index_meta`. `overwrite = FALSE` (default)
  is a no-op when `(index_name, trait_name)` already exists; `overwrite = TRUE`
  updates weights and economic weights in place. `economic_wts` is an optional
  numeric vector (same length as `trait_names`; some values may be 0) written to
  `index_meta.economic_weight`. Extra `...` columns are broadcast or per-trait.
- `add_index(tbl, index_name, value_col = NULL, overwrite_index = FALSE, delete_all = FALSE, ...)` —
  accepts a `tidybreed_table` from `get_table()` (optionally filtered). Any table
  with `id_ind`, `trait_name`, and a numeric value column is accepted: `ind_ebv`,
  `ind_phenotype`, `ind_tgv` (filter to one `component_name`), or a user-defined
  table. `value_col` is auto-detected from the table name (`ind_ebv` →
  `"ebv_value"`, `ind_phenotype` → `"pheno_value"`, `ind_tgv` → `"tgv_value"`);
  supply it explicitly for unknown tables.
  Multiplies each individual's values by the index weights in `index_meta` and
  appends to `ind_index`. Every individual must have exactly one value per index
  trait — an error is thrown if duplicates are found (filter to a single model /
  `eval_number` / `pheno_number` / `component_name` first); values are never
  summed across rows. Issues a warning when no filter is applied.
  `overwrite_index = TRUE` clears prior runs for the named index; `delete_all = TRUE`
  clears all of `ind_index`.

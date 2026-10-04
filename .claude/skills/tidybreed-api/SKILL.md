---
name: tidybreed-api
description: Reference for every implemented tidybreed function — arguments, semantics and internals of open_pop/define_genome, add_founders, define_chromosome, get_table subsets, schema(), define_trait/define_additive_effects, define_genome_effect_terms, define_phenotype, add_phenotype stages, add_tbv/add_tgv and the evaluator, index functions. Load before changing or explaining any exported function.
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

`add_founders()`, `add_phenotype()`, `add_tbv()`, and `add_ebv()` all accept
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
- `add_phenotype()` / `add_tbv()` / `add_ebv()`: scalar only (row count per
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

**Subset selection for action functions** (`add_phenotype`, `add_tbv`,
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
  k = 1 and k >= 2 (0.73.0):
  1. **Validate** everything; no RNG use before the draw (`seed` is applied
     after every check). `G` with manual `effects` or `scale_to_target = FALSE`
     is refused (nothing would be calibrated to it). `anchor = "realised"`
     needs an individuals `base_tbl` and the common scope.
  2. **Target** (`.dae_resolve_target()`, plan §6C): a passed `G` (matrix, or a
     number for one trait; dimnames checked, never relabelled; PSD) is written
     with the terms and refused over any stored block (whole-table check, even
     identical; the error gives a working `remove_rows()` call). `G = NULL`
     reads the stored rows (line block, else population-wide), or
     `trait_var_comp_tbl` rows. Refused: a stored block pairing a call trait
     with an outside trait; a stored `dominance` / `additive_by_additive` block
     for the traits; a partial block; two candidate sets. Manual effects take
     no target and skip these checks; unscaled k >= 2 draws read it only as
     the draw's Sigma.
  3. **Draw** the architecture with `.draw_additive_architecture()` (today's
     draws: rnorm / signed gamma for k = 1, `MASS::mvrnorm(G)` rows for k >= 2,
     masked per trait under `"union"`). The same seed gives the same `B0`.
  4. **Calibrate** with `.qtl_congruence()` (`R/qtl_congruence.R`, ported from
     the source method): `B = B0 A` with `B' M B = G` exactly. `M` is the
     anchor: `"genic"` = `diag(n_eligible p q)` at the base frequencies;
     `"realised"` = `Cov(X)` of the base individuals (collected in R, size
     guard `QTL_REALISED_MAX_CELLS`). k = 1 reduces to the scalar rescale.
     Two distinct rank errors: rank(G) > rank(M), rank(B0'MB0) < rank(G).
     `"union"` (k >= 2) scales each trait alone and warns "approximate" for a
     non-zero covariance.
  5. **Commit** terms and any new target in one transaction:
     `.ge_commit(..., before_commit = function(conn) .tvc_write_block(...))`.
  6. **Diagnostics** (common scope only, nothing stored): the relative
     spectrum of what another population sees vs the target — pool
     expectation `2 Cov(H)` (founder base), observed `Cov(X)` (individuals
     base), or the genic limit (realised anchor). The founder-pool
     comparison is always a `message()` (sampling LD, Q22), adding the
     realised-anchor hint outside `warn_bounds`; the other two warn outside
     `warn_bounds` (default `c(0.8, 1.25)`, `NULL` = off).
  7. **Message** says "exact" / "approximate" / "not calibrated" with the
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
  traits while supplying `G` is rejected, because the genetic covariance
  between a paternal-only and a maternal-only trait is zero under random mating
  and the requested off-diagonal is unobtainable, not merely approximate.

  `scale_to_target` is origin-aware:
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

**How the two writers relate.** `define_genome_effect_terms()` writes any effect
you supply; `define_*_effects()` functions sample effects of one shape and
write them through the same engine (`.ge_build → .ge_read_model →
.ge_resolve_deletes → .ge_commit`). `define_additive_effects()` is provably
sugar over the writer — `tests/testthat/test-genome-effects-writer.R`
("generator == writer") reproduces its output exactly through
`define_genome_effect_terms()` with the reserved owner, `replace_scope`, and the
same `base_tbl`. `base_tbl = NULL` deliberately differs: the generator has a
domain default; the writer fills nothing, so a missing centre is an error.

### `define_genome_effect_terms()` / `ad_terms()` / `genotype_terms()`

`R/define_genome_effect_terms.R`, `R/genome_effect_terms_builders.R`

`define_genome_effect_terms(pop, trait_name, terms, effect_owner = "custom", mode =
c("append", "replace_scope", "replace_owner", "replace_trait"), origin = NULL,
base_tbl = NULL, require_complete = FALSE, allow_reserved_owner = FALSE)` — the
general writer
for arbitrary genome effects. `terms` is a **long data frame, one row per
(term × locus)**: `term_id` (user-facing only; never stored), `genome_value`,
`effect_name`, `locus_name`, `contrast_name`, `center_value`,
`copy_count_value`, `dosage_value`. A single-term call may omit `term_id`.
Scope lives in a separate `origin` argument — `NULL` (the common scope), a
named scalar list applied to every member, or a data frame keyed by
`locus_name` for the exact multisets a genotype member takes.

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

Two builders produce `terms`, because a surface is rows, not a second
representation:

- `ad_terms(locus_name, a, d, p, coding = c("functional", "cockerham"))` —
  expands an (a, d) pair. Functional coding is `additive`@`0.5` plus
  `indicator`@`(2, 1)`; Cockerham is `additive`@`p` plus `dominance`@`p`. It
  **reports** the implied genetic mean `μ = a(p − q) + 2pq·d` and writes it
  nowhere — putting it in `phenotype_meta.mean` would double-count once
  non-additive genetic values reach the phenotype layer.
- `genotype_terms(genotypes, value, copy_count = NULL, drop_zero = TRUE)` —
  turns a genotype-by-value table into `indicator` terms, one term per row and
  one member per locus column.

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
  genetic-only);
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
    `NULL` (default) = simple single-self trait. A `group` contributor's
    mate sum accumulates exactly (`GEV_ACC_TYPE`, `.group_mate_tbv()`), so
    `group_sum()` / `group_mean()` are bit-identical across thread counts.
  - `prevalence` (categorical, two categories) — the threshold is
    `mean + qnorm(1 - prevalence) * sqrt(Va + Ve)`, with `Va` the trait's
    stored `additive` diagonal. Refused with `components` / `formula_tbv`
    (no stored variance describes a composite liability: use `thresholds`).
    `add_phenotype()` errors in PLAN (`.ap_check_prevalence()`, before any
    write or draw) when the trait has no stored `additive` row; there is no
    silent `Va = 0`. Skipped for `user_values` calls, which place no threshold.
  - `missing_component_action` — `"skip"` (default) or `"error"`. Stored in
    `phenotype_meta` and applied uniformly by `add_phenotype()` for **any**
    missing composite piece (missing group assignment, missing dam/sire TBV,
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

### `add_phenotype()` / `add_tbv()` / `add_tgv()`

`R/add_phenotype.R`, `R/add_tbv.R`, `R/add_tgv.R`, `R/genome_effects_eval.R`

Both functions accept a `tidybreed_table` (from `get_table()` + optional
`filter()`) as their first argument and return `tidybreed_pop`.

- `add_phenotype()` — the workhorse. `phenotype_name` (formerly `trait_name`)
  defaults to all phenotypes in `phenotype_meta` when omitted. Runs in
  **three stages** (`R/add_phenotype_stages.R`, `?add_phenotype_stages`):
  1. **PLAN** (`.ap_plan()`, no RNG, no writes except the `add_tbv()`
     prerequisite): sorted subset, metadata, topological sort of derived
     formulas, sex expression, repeatable guard, fixed-effect terms with
     `null_class_action`, TBV (simple from `ind_tbv`; composite via
     `.assemble_composite_tbv()`; `formula_tbv` via the DSL) with
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
  (the RNG-independent `add_tbv()` upsert is the one write that remains), and
  `.Random.seed` advanced by exactly the draws made before the error.
  Nothing in `R/` touches `.Random.seed`; never add seed restoration to one
  function — if the package ever adopts it, it is a package-wide policy.
  `tests/testthat/test-add_phenotype_failure_contract.R` asserts both halves.
- `add_tbv()` — TBV-only; no phenotype records. **One filtered call into the
  same evaluator `add_tgv()` uses** — reserved owner, order-one, contrast
  `additive` — never a second implementation of the effect math. Computes
  centered TBV from the
  order-one `additive` terms owned by `generated`: each allele copy
  takes the most specific variant whose origin predicate matches its
  `(line_origin, parent_origin)` label, falling back per copy to the common
  variant. This is what makes crossbreeding TBV correct (e.g. a Duroc × Landrace
  F1 centered against each parent line's own QTL effects and base allele
  frequency), and it is also how **imprinting** works now: a term scoped to one
  `parent_origin` reads only that parent's copies, per locus and per line
  rather than per trait. `trait_name` also defaults to all traits in
  `trait_meta` when omitted. Optional arguments for true index computation:
  - `index_names` — character vector of named indices; when supplied, multiplies
    per-trait TBVs by the index weights and writes results to `ind_true_index`.
    `NULL` (default) skips index computation.
  - `type` — `"index"` (default, uses `index_weight`), `"economic"` (uses
    `economic_weight`), or `"both"` (writes two rows per individual distinguished
    by `weight_type`).
  - `overwrite_index = FALSE` — when `FALSE`, skips individuals that already have
    a value in `ind_true_index` for the given `(index_name, weight_type)`. Set
    `TRUE` to recompute (e.g. after updating index weights).

  The filter is not conservatism. Under functional `(a, d)` input the stored
  coefficient is `a` while the breeding-value coefficient in a diploid HWE base
  is `α = a + d(q − p)`; under epistasis, average effects depend on other loci
  and on LD. So terms written through `define_genome_effect_terms()` move `ind_tgv`
  and never silently redefine `ind_tbv`, and additive members sitting inside an
  interaction are ignored.

  **`add_tbv()` warns once per trait when the coefficients it reads have
  stopped being average effects** — a non-reserved order-one `additive` term (it
  is part of A and is skipped), an `indicator` surface, or an interaction. An
  order-one `dominance` term centred where the additive term is centred is the
  exception and stays **silent**: Cockerham coding is HWE-orthogonal, so it
  contributes nothing to A and `tbv_value` remains exact. The warning is about
  *which terms were read*, never a claim that the arithmetic is wrong.
- `add_tgv()` — evaluates **every** term of a trait and writes `ind_tgv`, one
  row per (individual × trait × `component_name`). The raw sum of the stored
  terms; **no mean is added**. Idempotent per (individual, trait) — the delete
  is by trait, not by component, so a component that leaves the model leaves
  `ind_tgv` with it. Total via the `ind_tgv_total` view.

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
  `ind_phenotype`, `ind_tbv`, or a user-defined table. `value_col` is auto-detected
  from the table name (`ind_ebv` → `"ebv_value"`, `ind_phenotype` → `"pheno_value"`,
  `ind_tbv` → `"tbv_value"`); supply it explicitly for unknown tables.
  Multiplies each individual's values by the index weights in `index_meta` and
  appends to `ind_index`. Every individual must have exactly one value per index
  trait — an error is thrown if duplicates are found (filter to a single model /
  `eval_number` / `pheno_number` first). Issues a warning when no filter is applied.
  `overwrite_index = TRUE` clears prior runs for the named index; `delete_all = TRUE`
  clears all of `ind_index`.

# tidybreed — Developer & AI Context

## What This Package Does

`tidybreed` simulates breeding programs. It provides a
pipe-friendly API where all genomic and individual data is stored in a
file-based DuckDB database (not in R memory). Users build a `tidybreed_pop`
object step by step, query tables with `dplyr`, and eventually run selection
and mating cycles.

## Development Status: Pre-1.0.0 — Break Freely

**The package has no users yet. Until `1.0.0`, there is NO backward-compatibility
obligation of any kind.** Make whatever breaking changes produce the best design:
rename or remove columns, functions, and arguments; change schemas; alter default
behavior. **Do not** write compatibility shims, deprecation aliases, or "legacy"
code paths, and **do not** preserve old behavior for its own sake. If a redesign is
cleaner, take it — a breaking change is not a cost to weigh, it is the expected mode
of work in this phase.

No one should be treated as a downstream user before `1.0.0`; the few people who
have seen the package have been told not to build against it. Before `1.0.0`,
breaking changes are preferred whenever they improve the long-term schema, API,
correctness, or implementation.

When a function, argument, table, or column is renamed before `1.0.0`, remove the
old name completely from the codebase. Do not keep deprecated wrappers, aliases,
manual-page examples, roxygen examples, tests, comments, or generated docs that
teach the old name unless the explicit task is to document a migration that still
exists in current code.

Leftover compatibility files are technical debt, not harmless history. If an old
entry point has been replaced, delete the old R file, remove it from `NAMESPACE`,
remove its manual page, update examples/tests/vignettes, and make the new path
the only documented path. Do not leave "deprecated for later" functions around
before `1.0.0`; there is no downstream compatibility contract to protect.

**The only reproducibility contract is forward-looking, not cross-version:**

1. **Same seed reproduces within the current code.** A given `set.seed()` (or base
   seed) must produce identical output on repeated runs of the *current*
   implementation. We do **not** care whether it matches any previous version's
   output — never write a test that compares against pre-change/"golden-from-old"
   output, and never contort a formula to stay "byte-identical to today."
2. **R ↔ Rcpp parity.** Where the same algorithm exists in both R and C++, a given
   seed must produce identical output in both (this is a *within-current-code*
   guarantee, and the reason RNG choices like `dqrng` matter).

Before `1.0.0`, algorithms may change and seeded output may change across
versions. Within a single implementation, seeded runs should be reproducible.
After `1.0.0`, users must be able to reproduce simulations exactly for the same
package version, seed, inputs, and supported platform. Any algorithm implemented
in both R and Rcpp/C++ must have explicit parity tests showing identical results
for the same seed.

When `1.0.0` ships, this section is replaced by a normal semantic-versioning
compatibility policy. Until then: design first, break as needed.

## Design Principles

1. **Database-first** — all data lives in DuckDB, not R objects
2. **File-based by default** — enables simulations larger than RAM; supports
   replicates, resumable runs, and sharing easily
3. **Lazy evaluation** — use `dplyr::tbl()` / `get_table()` and filter before
   `collect()`-ing into R
4. **Pipe-friendly** — most exported functions accept a `tidybreed_pop` and
   return a `tidybreed_pop`; action functions (`add_phenotype`, `add_tbv`,
   `add_genotypes`, `extract_genotypes`, `define_chip`, `define_additive_effects`)
   accept a `tidybreed_table` from `get_table()` and return `tidybreed_pop`
5. **Type-safe** — all table columns have explicit DuckDB types; user-added
   columns are inferred via `infer_duckdb_type()`
6. **Disdain and intolerance for storing metadata** - storing data such as
   n_loci for loci count is silly when the function that needs that data
   can run a basic SQL query to pull out info and makes the program truly modular
7. **No implicit ordering** — never accept positional vectors (logical TRUE/FALSE
   or integer indices) to select rows from a database table. Row order is not
   guaranteed and changes silently when rows are added or reordered. All
   selection flows through named identifiers (`id_ind`, `locus_name`, `locus_id`
   as PK) or explicit SQL filter predicates via `get_table() |> filter()`.

If metadata is required, it needs to be considered carefully and likely 
should be stored in a table and never in the R object itself, for restarts
and for `restore_pop()` to be complete without having to recreate it or
given as options again. 

## Schema Design Bias

Schema choices should optimize for future breeding-program designs, especially
crossbreeding. Before adding or reshaping a table, ask whether the schema can
later support multiple lines, line-specific and population-wide effects,
reciprocal crosses, sex-specific recombination maps, non-additive effects, and
contributor-specific phenotype components without another fundamental rewrite.

Do not implement future biology before it is needed. It is enough to reserve
clean dimensions now when the schema would be painful to alter later. Prefer
long tables with explicit dimensions such as `line_name`, `sex`, `map_name`,
`genome_effect_type`, and `effect_name`. Use `NULL` deliberately for
shared/default behavior, such as population-wide genome effects or maps applying
to all lines.

## Naming Consistency Rules

1. **Function argument names must match the database column they populate exactly.**
   Never use a different name for the same thing (e.g. argument `line_name` → column
   `line_name`, not `line`).

2. **All primary numeric value columns follow the `{prefix}_value` pattern.**
   Examples: `pheno_value` in `ind_phenotype`, `tbv_value` in `ind_tbv`,
   `ebv_value` in `ind_ebv`, `cov_value` in `trait_var_comp`,
   `index_value` in `ind_index`.

3. **All name/label columns end in `_name`.**
   Examples: `trait_name`, `locus_name`, `chr_name`, `line_name`, `effect_name`,
   `index_name`. Never abbreviate to just `trait`, `line`, `chr`, etc. when
   used as a column or function parameter that maps to one of these columns.

4. **All ID foreign-key columns start with `id_`.**
   Examples: `id_ind`, `id_trait`, `id_ebv`, `id_tbv`. No `phenotype_id`-style names.

5. **No abbreviations in column names when the full word is unambiguous.**
   `index_weight` not `index_wt`; `trait_name_1`/`trait_name_2` not `trait_1`/`trait_2`.

## Two-Layer Phenotype Design (v0.31.0+)

The model is split into two distinct layers with a strict boundary between them:

**Genetic component layer** — managed by `define_trait()`:
- One row in `trait_meta` per underlying genetic quantity (e.g. `ADG_direct`, `ADG_social`, `WWD`, `WWM`)
- Has QTL effects in `genome_effects`, TBVs in `ind_tbv`, additive variance in `trait_var_comp`
- Arguments: `target_add_var`, `target_add_mean`, `expressed_parent`, `description`, `units`
- No phenotype-level information at all — no mean, no residual, no type, no expressed_sex

**Observation layer** — managed by `define_phenotype()`:
- One row in `phenotype_meta` per observed phenotype individuals receive records for (e.g. `ADG`, `WW`, `mortality`)
- For simple traits, `phenotype_name` equals the `trait_name` of its single genetic component
- For composite traits (maternal, SGE), `phenotype_name` is new and one or more `trait_meta` rows feed into it via `phenotype_components`
- Arguments: `type`, `mean`, `expressed_sex`, `repeatable`, `min_value`, `max_value`, `prevalence`, `thresholds`, `cat_values`, `cat_names`, `store_liability`, `residual_var`, `components`, `formula_tbv`, `formula`, `missing_component_action`

**The rule**: if an argument describes the genetics (variance, QTL structure, parent-of-origin), it belongs in `define_trait()`. If it describes what observers record (mean, distribution, sex expression, residual noise, how to assemble from components), it belongs in `define_phenotype()`. Never put observation-layer arguments on `define_trait()` or genetic-layer arguments on `define_phenotype()`.

## Function Naming Convention

| Prefix        | Meaning                                        |
|---------------|------------------------------------------------|
| `open_`       | Opens or creates a population database/session |
| `restore_`    | Restores an existing population database       |
| `add_`        | Inserts simulation output rows                 |
| `define_`     | Writes model configuration / metadata          |
| `mutate_`     | Adds or updates columns in an existing table   |
| `extract_`    | Returns analysis/export data without changing simulation state |
| `remove_`     | Deletes selected rows                          |
| `archive_`    | Moves/stamps completed replicate data          |

**`add_*` vs `define_*` rule**: if the function writes rows that represent
simulation *output* (data produced by running the model), use `add_`. If the
function writes rows that configure *how* the model runs (parameters, weights,
effect definitions), use `define_`.

Examples: `add_founders()`, `add_phenotype()`, `add_tbv()`, `add_ebv()`,
`add_index()` — all write simulation output.  
`define_trait()`, `define_additive_effects()`, `define_effect_cov_matrix()`,
`define_chip()`, `define_index()` — all write model configuration.

## Database Schema

### `genome_meta`

Locus-level metadata. One row per locus. Holds the **physical** coordinate only
(`pos_bp`); the **genetic** map (`pos_cM`) lives in the separate `genome_map` table
(see below), mirroring how QTL effects live in `genome_effects`.

| Column     | Type    | Notes                          |
|------------|---------|--------------------------------|
| locus_id   | INTEGER | Primary key (1, 2, …, n_loci)  |
| locus_name | VARCHAR | e.g. "Locus_1", "rs12345"; validated unique |
| chr        | INTEGER | Chromosome number              |
| chr_name   | VARCHAR | Chromosome name string         |
| pos_bp     | BIGINT  | Physical position, base pairs, **1-based** (VCF/PLINK convention). Created via explicit typed `CREATE TABLE` so the type is enforced; user row-adds of an integer widen to `BIGINT` automatically |
| founder_allele_freq | DOUBLE | Base allele frequency; added by `define_founder_haplotypes()` (via `ALTER TABLE`, never a table rewrite — preserves `pos_bp` `BIGINT`) |
| *user cols*| any     | Added via `mutate_table()` or `define_chip()` |

**Reserved** (cannot be modified): `locus_id`, `locus_name`, `chr`, `chr_name`, `pos_bp`

Example user columns: `is_50K BOOLEAN`, `is_HD BOOLEAN`

**Note**: QTL effects are **not** stored as columns in `genome_meta`. They live in
the `genome_effects` table (see below). There are no `add_{trait}`, `is_QTL_{trait}`,
or `base_allele_freq_{trait}` columns. **Genetic-map positions (`pos_cM`) are not in
`genome_meta`** either — use `genome_map`, joined via `locus_id`. (`pos_Mb` and
`introduced_gen` were removed in v0.50.0.)

### `genome_map`

The genetic map, in **long** format. One row per (locus × sex × line × map) with a
defined genetic position. Populated by `define_genome()` (a single default map) and,
later, by a `define_genetic_map()`-style writer for sex/line/version-specific maps.
Adding a map dimension is **rows, never a schema change** — the same precedent as
`genome_effects`.

| Column        | Type    | Notes                                                        |
|---------------|---------|--------------------------------------------------------------|
| id_genome_map | INTEGER | Surrogate PK (assigned via `next_int_id()`, not DB auto-increment) |
| locus_id      | INTEGER | FK to `genome_meta.locus_id`; internal join/order key        |
| locus_name    | VARCHAR | FK to `genome_meta.locus_name`; denormalized                 |
| sex           | VARCHAR | `NULL` = both sexes; `'M'`/`'F'` = sex-specific map           |
| line_name     | VARCHAR | `NULL` = all lines; set for line-specific maps               |
| map_name      | VARCHAR | Map version/identity; default `"default"`                    |
| pos_cM        | DOUBLE  | Genetic-map position, centiMorgans                           |

**Logical key** `(locus_id, sex, line_name, map_name)` (nullable `sex`/`line_name`
enforced in R by `validate_genome_map()`). **Reserved**: all columns.

Two internal helpers are the single source of genetic positions for all
distance-driven code (founder LD, recombination):
- `resolve_genome_map(conn, sex, line_name, map_name)` — returns exactly one row per
  locus, `locus_id`-ordered, applying per-locus precedence `(sex=S,line=L)` →
  `(sex=S,NULL)` → `(NULL,line=L)` → `(NULL,NULL)`; errors on a missing locus or a
  `pos_cM` that is non-monotonic within a chromosome after fallback.
- `validate_genome_map(conn)` — logical-key uniqueness (NULL-normalized),
  agreement with `genome_meta`, valid `sex`/`map_name`. Run after every map write.

### `genome_effects`

QTL effect data. One row per (locus × trait × effect type × line). Populated by
`define_additive_effects()` (single or multi-trait).

| Column             | Type    | Notes                                                      |
|--------------------|---------|------------------------------------------------------------|
| id_genome_effect   | INTEGER | Primary key assigned by tidybreed via `next_int_id()`       |
| locus_name         | VARCHAR | FK to `genome_meta.locus_name`                             |
| line_name          | VARCHAR | NULL = population-wide; set for line-specific effects      |
| trait_name         | VARCHAR | FK to `trait_meta.trait_name`                              |
| genome_effect_type | VARCHAR | `"additive"` now; `"dominance"` and others later           |
| genome_value       | DOUBLE  | Effect size                                                |
| base_allele_freq   | DOUBLE  | Base allele frequency used for TBV centering (Falconer)    |

**Reserved**: all columns (the table is managed exclusively by `define_additive_effects()`).

QTL membership is **implicit**: a locus is a QTL for a trait if it has a row in
`genome_effects` for that `(trait_name, genome_effect_type)`. No separate boolean
flag is stored.

For crossbreeding, `line_name = NULL` means a population-wide/common effect shared
across lines; non-NULL `line_name` rows are line-specific effects. Current code
implements additive effects only, but the long `genome_effect_type` dimension is
reserved so dominance, epistasis, and other non-additive effects can be added later
without a table rewrite.

### `ind_haplotype`

Phased haplotypes in **long** format. One row per (individual × haplotype ×
locus). Populated by `add_founders()` and `add_offspring()`. Row count per
individual per chromosome follows the resolved `chr_inheritance`
`from_parent_1`/`from_parent_2` for that individual's sex: 2 rows/locus for a
plain autosome (`1, 1`, the default), 1 for a hemizygous sex chromosome (e.g.
`0, 1`), 0 for an absent chromosome (`0, 0`, e.g. Y in females). See
`chr_inheritance` below and `define_chromosome()`.

| Column        | Type     | Notes                                                       |
|---------------|----------|-------------------------------------------------------------|
| id_ind        | VARCHAR  | FK to `ind_meta.id_ind`; part of composite PK               |
| parent_origin | UTINYINT | 1 = from parent_1 (sire), 2 = from parent_2 (dam); PK part  |
| strand        | UTINYINT | Copy within a parent's contribution; always 1 for diploids; PK part |
| line_origin   | VARCHAR  | Founding line this allele traces to; used by `add_tbv()` for line-specific crossbreeding TBV |
| locus_id      | INTEGER  | FK to `genome_meta.locus_id`; physical sort/PK key          |
| locus_name    | VARCHAR  | FK to `genome_meta.locus_name`; denormalized for direct `genome_effects` joins |
| allele        | UTINYINT | 0 or 1 (phased)                                             |

**Primary key**: `(id_ind, parent_origin, strand, locus_id)`.

### `ind_genotype`

Genotype dosages in **long** format, 0/1/2 encoding. One row per (individual ×
locus). **On-demand cache** — starts empty and is populated only by
`add_dosage()` (never by `add_founders()`/`add_offspring()`). May be empty or
partial.

| Column       | Type     | Notes                                        |
|--------------|----------|----------------------------------------------|
| id_ind       | VARCHAR  | FK to `ind_meta.id_ind`; part of composite PK |
| locus_id     | INTEGER  | FK to `genome_meta.locus_id`; part of composite PK |
| locus_name   | VARCHAR  | FK to `genome_meta.locus_name`               |
| dosage_value | UTINYINT | Sum of alleles across strands (0/1/2 diploid) |

**Primary key**: `(id_ind, locus_id)`. Populated via `INSERT OR REPLACE`
(idempotent).

Haplotypes are the source of truth; dosage is derived on demand (SUM of alleles)
rather than auto-stored, because dosage is cheap to recompute and only needed for
specific downstream analyses (MAS, GBLUP export, allele frequencies).

### `chr_inheritance`

Per-chromosome copy counts, in **long** format, keyed by **offspring** sex. One
row per `(chr_name, offspring_sex, line_name)`. Created by `define_genome()` with
one seeded default row per chromosome (`from_parent_1 = from_parent_2 = 1`, a
plain diploid autosome). `define_chromosome()` sets non-default rules (sex
chromosomes, organelles). Answers: *"an offspring of sex S inherits N copies of
this chromosome from each parent."* Real polyploidy (ploidy > 2) is not yet
supported — see `ind_meta.ploidy`.

| Column        | Type     | Notes                                                     |
|---------------|----------|-----------------------------------------------------------|
| chr_name      | VARCHAR  | FK to `genome_meta.chr_name` (R-enforced)                 |
| offspring_sex | VARCHAR  | `NULL` = all/default; `'M'`/`'F'` — the **offspring** (carrier) sex |
| line_name     | VARCHAR  | `NULL` = all lines; reserved for line-specific rules (crossbreeding) |
| from_parent_1 | UTINYINT | Absolute copies inherited from parent_1 (sire), at ploidy 2 |
| from_parent_2 | UTINYINT | Absolute copies inherited from parent_2 (dam), at ploidy 2  |

**Logical key** `(chr_name, offspring_sex, line_name)`, NULL-normalized in R.
Counts are **absolute** (correct at ploidy 2, enforced): autosome `1,1`; male's X
`0,1`; male's Y `1,0`; female's Y `0,0`; maternal mito `0,1`. Row-local
`CHECK`/`NOT NULL` constraints enforce `offspring_sex IN ('M','F')` (NULL passes),
non-negative counts, and `from_parent_1 + from_parent_2 <= 2` (a diploid-release
constraint). `from_parent_1`/`from_parent_2` map directly onto
`ind_haplotype.parent_origin` (1 = sire, 2 = dam) and `strand`.

### `chr_recombination`

Per-chromosome recombination, in **long** format, keyed by **producing-parent**
sex. One row per `(chr_name, parent_sex, line_name)`. Seeded by `define_genome()`
from its genome-wide `recombines_M`/`recombines_F` defaults (one `parent_sex =
NULL` row per chromosome when both agree, else a `'M'` and an `'F'` row).
`define_chromosome()` sets non-default rules (Y, W, achiasmy). Answers: *"when a
parent of sex S makes gametes, does this chromosome recombine?"*

| Column     | Type    | Notes                                                        |
|------------|---------|--------------------------------------------------------------|
| chr_name   | VARCHAR | FK to `genome_meta.chr_name` (R-enforced)                    |
| parent_sex | VARCHAR | `NULL` = both parents; `'M'`/`'F'` — the **producing-parent** sex |
| line_name  | VARCHAR | `NULL` = all lines; reserved for line-specific rules         |
| recombines | BOOLEAN | TRUE if the chromosome recombines in that parent sex's meiosis |

**Logical key** `(chr_name, parent_sex, line_name)`, NULL-normalized in R.

**Why two tables:** "which sex" means the **offspring** for copy count but the
**producing parent** for recombination — one table would force one `sex` column to
mean both. Splitting them lets each column say what it means (`offspring_sex` vs
`parent_sex`), and the two concerns resolve **independently** (a copy rule can
never shadow a recombination rule).

**Resolution.** Two internal resolvers mirror `resolve_genome_map()`'s
priority-window fallback `(sex=S,line=L) → (sex=S,NULL) → (NULL,line=L) →
(NULL,NULL)`:
- `resolve_chr_inheritance(conn, offspring_sex, line_name)` — called with the
  **offspring's** sex and line; returns `(from_parent_1, from_parent_2)` per chr.
- `resolve_chr_recombination(conn, parent_sex, line_name)` — called with the
  **producing parent's** sex and line; returns `recombines` per chr.

`validate_chr_inheritance()`/`validate_chr_recombination()` run inside every write
transaction: NULL-normalized key uniqueness, orphan-`chr_name` (R-enforced FK),
valid sex, `sum ≤ 2`, resolvability for both `M` and `F`, and a deterministic
sex-vs-line shadowing check. A chromosome takes the fast autosome path only when
its resolved inheritance is `1,1` for both offspring sexes **and** `recombines`
for both parent sexes. `2,0` (uniparental disomy) is storage-expressible but
errors at the `add_offspring()`/`add_founders()` kernel boundary (unimplemented
transmission mechanism).

### `founder_haplotypes`

Founder haplotype pool in **long** format. Created by
`define_founder_haplotypes()`, sampled by `add_founders()`.

| Column       | Type     | Notes                                                    |
|--------------|----------|----------------------------------------------------------|
| line_name    | VARCHAR  | NULL = shared pool; set for line-specific pools; logical key part |
| haplotype_id | INTEGER  | Sequential within the pool (unique per `line_name`); logical key part |
| locus_name   | VARCHAR  | FK to `genome_meta.locus_name`; logical key part         |
| allele       | INTEGER  | 0 or 1                                                   |

Logical key `(line_name, haplotype_id, locus_name)` enforced in R (nullable
`line_name`), matching the `index_meta` convention.

**Reserved**: all columns (the table is a sampling pool written exclusively by
`define_founder_haplotypes()` and read by `add_founders()`).

### `ind_meta`

Individual-level metadata. Created empty by `open_pop()`; rows
populated by `add_founders()` and `add_offspring()`.

| Column      | Type    | Notes                               |
|-------------|---------|-------------------------------------|
| id_ind      | VARCHAR  | Primary key, format `{line_name}_{n}` (e.g. `Libra_1020`) |
| id_parent_1 | VARCHAR  | NA for founders                     |
| id_parent_2 | VARCHAR  | NA for founders                     |
| line_name   | VARCHAR  | Genetic line name                   |
| sex         | VARCHAR  | "M" or "F"                          |
| ploidy      | UTINYINT | Genome ploidy; declared at `add_founders()` time (must be `2` in this version), computed at `add_offspring()` time as the sum of each parent's gamete contribution (`own_ploidy / 2` per parent). Default `2`. |
| *user cols* | any      | Added via `mutate_table()` or `...` in `add_founders()` |

**Reserved**: `id_ind`, `id_parent_1`, `id_parent_2`, `line_name`, `sex`, `ploidy`

### `trait_meta`

One row per **genetic component trait**. Populated by `define_trait()`.
Contains only genetic-layer information — no phenotype-level metadata.
Observation-layer metadata lives in `phenotype_meta`.

| Column          | Type    | Notes                                                              |
|-----------------|---------|--------------------------------------------------------------------|
| id_trait        | INTEGER | Primary key assigned by tidybreed via `next_int_id()`               |
| trait_name      | VARCHAR | Unique identifier; equals `phenotype_name` for simple traits       |
| description     | VARCHAR | Free text                                                          |
| units           | VARCHAR | e.g. `"kg"`, `"g/day"`                                             |
| expressed_parent| VARCHAR | `"both"` (default), `"parent_1"` (paternal), `"parent_2"` (maternal) — imprinting |
| target_add_mean | DOUBLE  | TBV centering mean for the base population; default `0`            |

**What does NOT belong here** (all moved to `phenotype_meta` in v0.31.0):
`type`, `expressed_sex`, `repeatable`, `mean`, `min_value`, `max_value`,
`prevalence`, `thresholds`, `cat_values`, `cat_names`, `residual_var`,
`index_weight`, `economic_value`.

### `trait_var_comp`

Genetic-layer variance components. One row per (effect_name, trait_name_1, trait_name_2).
Both `(i,j)` and `(j,i)` pairs stored. Populated by `define_effect_cov_matrix()` and
`define_trait()`. Stores **only** genetic effects — no phenotype-level variances.

Valid `effect_name` values: `"gen_add"` (additive genetic G matrix);
future: `"dominance"`, `"epistasis"`. Named random effects (HYS, litter, pen)
go to `phenotype_var_comp`, not here.

| Column           | Type    | Notes                                              |
|------------------|---------|----------------------------------------------------|
| id_trait_var_comp| INTEGER | Primary key assigned by tidybreed via `next_int_id()` |
| effect_name      | VARCHAR | `"gen_add"`; future: `"dominance"`, `"epistasis"`  |
| trait_name_1     | VARCHAR |                                                    |
| trait_name_2     | VARCHAR |                                                    |
| cov_value        | DOUBLE  | Variance (diagonal) or covariance (off-diagonal)   |

### `phenotype_meta`

Observed phenotype definitions. One row per phenotype name. Populated by
`define_phenotype()`. Analogous to `trait_meta` but for the observation layer —
simple traits have the same name in both tables; composite phenotypes (e.g. WW,
SGE ADG) appear only here.

| Column                   | Type    | Notes                                                         |
|--------------------------|---------|---------------------------------------------------------------|
| id_phenotype_meta        | INTEGER | Primary key assigned by tidybreed via `next_int_id()`          |
| phenotype_name           | VARCHAR | Unique. Equals `trait_name` for simple traits.                |
| type                     | VARCHAR | `"continuous"`, `"count"`, `"categorical"`, `"derived_formula"` |
| mean                     | DOUBLE  | Phenotypic population mean / liability intercept              |
| expressed_sex            | VARCHAR | `"both"`, `"M"`, or `"F"`                                     |
| repeatable               | BOOLEAN | Repeated records allowed?                                     |
| min_value / max_value    | DOUBLE  | Clipping bounds for count traits                              |
| prevalence               | DOUBLE  | For 2-category categorical traits                             |
| thresholds               | VARCHAR | Comma-separated liability cutpoints for K-category traits     |
| cat_values               | VARCHAR | Comma-separated numeric phenotype values per category         |
| cat_names                | VARCHAR | Comma-separated labels per category                           |
| store_liability          | BOOLEAN | Write raw liability to `ind_phenotype.liability_value`        |
| missing_component_action | VARCHAR | `"skip"` (default) or `"error"` — what to do when any component of a composite phenotype cannot be resolved for an individual |

**Reserved**: all columns (managed by `define_phenotype()`).

### `phenotype_components`

Component definitions for composite phenotypes. One row per (phenotype ×
component). Populated by `define_phenotype(..., components = ...)`. Simple
(non-composite) phenotypes have no rows here.

| Column             | Type    | Notes                                                              |
|--------------------|---------|--------------------------------------------------------------------|
| id_phenotype_comp  | INTEGER | Primary key assigned by tidybreed via `next_int_id()`               |
| phenotype_name     | VARCHAR | FK to `phenotype_meta.phenotype_name`                              |
| source_trait_name  | VARCHAR | FK to `trait_meta.trait_name` — the genetic component trait        |
| contributor_type   | VARCHAR | `"self"`, `"dam"`, `"sire"`, or `"group"`                          |
| group_column       | VARCHAR | Column in `group_table` that holds group membership (required for `"group"`) |
| group_table        | VARCHAR | Table containing `group_column`; default `"ind_meta"`              |
| aggregation        | VARCHAR | `"sum"` (default) or `"mean"` — how group-mates' TBVs are combined |
| weight             | DOUBLE  | Scalar multiplier; default `1.0`                                   |
| weight_type        | VARCHAR | `"fixed"` (default), `"covariate"`, `"legendre"`, `"raw_poly"`     |
| covariate_name     | VARCHAR | Covariate column when `weight_type = "covariate"`                  |
| covariate_table    | VARCHAR | Table containing `covariate_name`                                  |
| poly_order         | INTEGER | Polynomial basis order                                             |
| poly_scale_min/max | DOUBLE  | Legendre scaling bounds                                            |
| genome_effect_types| VARCHAR | Default `"additive"`                                               |
| missing_action     | VARCHAR | Per-component fallback (currently unused; use `phenotype_meta.missing_component_action`) |
| contributor_filter | VARCHAR | Reserved for spatial/neighborhood lookup — not yet implemented     |

**Note on SGE (Social Genetic Effects / Bijma model)**: for `contributor_type = "group"`,
`add_phenotype()` aggregates group-mates' TBVs (excluding self). A singleton (no
group-mates) receives a social contribution of 0 and is not excluded. An individual
with no group assignment receives `NA` and is handled by `missing_component_action`.

**Reserved**: all columns (managed exclusively by
`define_phenotype(..., components = ...)`).

### `phenotype_effects`

Non-additive-genetic, non-residual terms in the phenotype model (fixed and
random effects). One row per (phenotype × effect). An **observation-layer**
table, keyed by `phenotype_name` (FK to `phenotype_meta`), not by `trait_name`.

| Column            | Type    | Notes                                                    |
|-------------------|---------|----------------------------------------------------------|
| phenotype_name    | VARCHAR | FK to `phenotype_meta.phenotype_name`; PK part           |
| effect_name       | VARCHAR | e.g. "sex", "gen", "litter"; PK part                     |
| effect_class      | VARCHAR | `"fixed_class"`, `"fixed_cov"`, or `"random"`            |
| source_column     | VARCHAR | Column in source table used as grouping variable         |
| source_table      | VARCHAR | Table containing `source_column` (default `"ind_meta"`)  |
| distribution      | VARCHAR | For random effects: `"normal"`, `"gamma"`, `"uniform"`   |
| levels_json       | VARCHAR | For fixed_class effects: JSON `{"M":30,"F":0}`           |
| slope             | DOUBLE  | For fixed_cov effects: regression coefficient            |
| center            | DOUBLE  | For fixed_cov effects: centering value                   |
| value             | DOUBLE  | Rarely used scalar                                       |
| poly_order        | INTEGER | Polynomial order for covariate effects; default 1        |
| null_class_action | VARCHAR | Behavior when the grouping column is NULL; default `"skip"` |

**Primary key**: `(phenotype_name, effect_name)`.

**Reserved**: all columns (managed by `define_effect_fixed_class()`,
`define_effect_fixed_cov()`, and `define_effect_random()`).

### `phenotype_random_effects`

Sampled draws for the random effects declared in `phenotype_effects`. One row per
(phenotype × effect × level), written by `define_effect_random()` and read by
`add_phenotype()`. Pure observation-layer noise — no genetic content.

| Column         | Type    | Notes                                              |
|----------------|---------|-----------------------------------------------------|
| phenotype_name | VARCHAR | FK to `phenotype_meta.phenotype_name`; PK part     |
| effect_name    | VARCHAR | FK to `phenotype_effects.effect_name`; PK part     |
| level          | VARCHAR | Level of the grouping column; PK part              |
| draw_value     | DOUBLE  | Sampled deviation for that level                   |
| date_sampled   | DATE    | When the draw was taken                            |

**Primary key**: `(phenotype_name, effect_name, level)`.

**Reserved**: all columns (managed exclusively by `define_effect_random()`).

### `phenotype_var_comp`

Phenotype-layer variance components. One row per (effect_name, phenotype₁, phenotype₂,
condition). Both (i,j) and (j,i) pairs stored for off-diagonal entries. Populated by
`define_phenotype(..., residual_var = ...)`, `define_residual_cov()`, and
`define_effect_random()`.

`effect_name = 'residual'` is reserved for residual noise. All named random effects
(HYS, litter, pen, etc.) use their own string (e.g. `'hys'`, `'litter'`).
The `condition_column` / `condition_level` columns are used only for `'residual'`
to model heterogeneous residual variance by sex, group, etc.

| Column               | Type    | Notes                                                              |
|----------------------|---------|--------------------------------------------------------------------|
| id_phenotype_var_comp| INTEGER | Primary key assigned by tidybreed via `next_int_id()`              |
| effect_name          | VARCHAR | `'residual'` or any named random effect (e.g. `'hys'`, `'litter'`)|
| phenotype_name_1     | VARCHAR |                                                                    |
| phenotype_name_2     | VARCHAR |                                                                    |
| cov_value            | DOUBLE  |                                                                    |
| condition_column     | VARCHAR | NULL = unconditional; used only for `effect_name = 'residual'`     |
| condition_table      | VARCHAR | Default `"ind_meta"`                                               |
| condition_level      | VARCHAR | Value of `condition_column` for this row                           |
| weight_type          | VARCHAR | Default `"fixed"`                                                  |
| poly_order           | INTEGER | Polynomial order for `"legendre"` weight type                      |

### `ind_phenotype`

Phenotype records in long format. Populated by `add_phenotype()`.

| Column         | Type    | Notes                                             |
|----------------|---------|---------------------------------------------------|
| id_phenotype   | INTEGER | Primary key assigned by tidybreed via `next_int_id()` |
| id_ind         | VARCHAR |                                                   |
| phenotype_name | VARCHAR | FK to `phenotype_meta.phenotype_name`             |
| pheno_value    | DOUBLE  | Phenotype value                                   |
| pheno_number   | INTEGER | 1 = first record for this individual × trait, etc.|
| *user cols*    | any     | Added via `mutate_table()` or scalar `...` in `add_phenotype()` |

### `ind_tbv`

True breeding values (simulation ground truth). Populated by
`add_phenotype()` and `add_tbv()`. Logical key `(id_ind, trait_name)` unique.

| Column     | Type    | Notes                                  |
|------------|---------|----------------------------------------|
| id_tbv     | INTEGER | Primary key assigned by tidybreed via `next_int_id()` |
| id_ind     | VARCHAR |                                        |
| trait_name | VARCHAR |                                        |
| tbv_value  | DOUBLE  |                                        |

### `ind_ebv`

Estimated breeding values from external BLUP / GBLUP runs. Logical key
`(id_ind, trait_name, model, eval_number)` unique. Populated by `add_ebv()`.

| Column      | Type    | Notes                                                   |
|-------------|---------|--------------------------------------------------------------|
| id_ebv      | INTEGER | Primary key assigned by tidybreed via `next_int_id()`         |
| id_ind      | VARCHAR |                                                              |
| trait_name  | VARCHAR |                                                              |
| model       | VARCHAR | User label, e.g. "ssGBLUP_v1"                               |
| ebv_value   | DOUBLE  |                                                              |
| acc         | DOUBLE  | Optional accuracy                                            |
| se          | DOUBLE  | Optional standard error                                      |
| eval_number | INTEGER | Auto-incrementing counter per trait (global across models); 1 = first evaluation |

### `index_meta`

Selection index definitions. One row per (index × trait). Populated by
`define_index()`. A special row with `index_name = NULL` is written by
`define_trait()` to record the global economic weight for each trait.

| Column          | Type    | Notes                                                         |
|-----------------|---------|---------------------------------------------------------------|
| id_index_name   | INTEGER | Primary key assigned by tidybreed via `next_int_id()`          |
| index_name      | VARCHAR | NULL = global/default entry written by `define_trait()`; named index written by `define_index()` |
| trait_name      | VARCHAR | FK to `trait_meta.trait_name`                                 |
| index_weight    | DOUBLE  | Selection index weight (NULL for global rows)                 |
| economic_weight | DOUBLE  | Economic value per unit of the trait                          |
| *user cols*     | any     | Added via `...` in `define_index()`                           |

**Reserved**: `id_index_name`, `index_name`, `trait_name`, `index_weight`, `economic_weight`

**Unique constraint**: `(index_name, trait_name)` (NULL `index_name` uniqueness enforced in R code).

### `ind_index`

Computed selection index values. One row per (individual × index × run).
Populated by `add_index()`.

| Column       | Type    | Notes                                          |
|--------------|---------|------------------------------------------------|
| id_index     | INTEGER | Primary key assigned by tidybreed via `next_int_id()` |
| id_ind       | VARCHAR |                                                |
| index_name   | VARCHAR | FK to `index_meta.index_name`                  |
| index_number | INTEGER | Auto-incrementing run counter per individual   |
| index_value  | DOUBLE  | Computed index value                           |
| *user cols*  | any     | Added via `...` in `add_index()`               |

### `ind_true_index`

True selection index values computed from TBVs (simulation ground truth).
Populated by `add_tbv()` when `index_names` is supplied.

| Column           | Type    | Notes                                                        |
|------------------|---------|--------------------------------------------------------------|
| id_true_index    | INTEGER | Primary key assigned by tidybreed via `next_int_id()`         |
| id_ind           | VARCHAR |                                                              |
| index_name       | VARCHAR | FK to `index_meta.index_name`                                |
| weight_type      | VARCHAR | `"index"` (uses `index_weight`) or `"economic"` (uses `economic_weight`) |
| true_index_value | DOUBLE  | Weighted sum: `sum(weight_i * tbv_i)` across index traits    |

**Reserved**: all columns (managed exclusively by `add_tbv()` when `index_names` is supplied).

Logical row key `(id_ind, index_name, weight_type)` — one true index value per
individual × index × weight type. No SQL `UNIQUE` constraint; uniqueness enforced
in R via DELETE + INSERT when `overwrite_index = TRUE`.

## Implemented Functions (Phase 1 Complete)

### `open_pop()` + `define_genome()` (genome setup)

`R/open_pop.R`, `R/define_genome.R`

The current surface for creating a population and its genome is
`open_pop() |> define_genome(...)`. `define_genome()` populates the genome tables:

- Genome: `genome_meta` (physical `pos_bp`), `genome_map` (default map), `ind_haplotype` (empty), `ind_genotype` (empty), `chr_inheritance` + `chr_recombination` (default autosome rows)

`define_genome()` key params: `pop`, `n_loci`, `n_chr`, `chr_len_Mb` (finite,
strictly positive), `cM_per_Mb` (genetic-map rate, cM per Mb; scalar or
length-`n_chr`, finite, strictly positive, default `1.0` →
`pos_cM = pos_bp/1e6 * cM_per_Mb`), `locus_names`, `chr_names`,
`recombines_M`/`recombines_F` (genome-wide per-parent-sex recombination defaults,
both `TRUE`; set one `FALSE` for a whole-genome achiasmatic sex, seeded into
`chr_recombination`). Calling
`define_genome()` on a population that already has a non-empty `genome_meta` is a
hard error (no partial re-definition).

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
`add_genotypes`, `extract_genotypes`) requires `get_table()` as the first step.
`filter()` is called on the `tidybreed_table`, not on the pop directly.
The unique `id_ind` values from the collected filtered table are used as the
candidate set. Any table that has an `id_ind` column can be used (e.g.
`ind_meta`, `ind_phenotype`, `ind_haplotype`, `ind_genotype`).

```r
# All individuals
pop |> get_table("ind_meta") |> add_phenotype("ADG")

# Filtered by ind_meta
pop |>
  get_table("ind_meta") |>
  dplyr::filter(sex == "F", gen == 1L) |>
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
  `index_meta`. Accepted arguments: `trait_name`, `target_add_var` (writes to
  `trait_var_comp`), `target_add_mean`, `expressed_parent`, `description`,
  `units`, `overwrite`. **Never** pass observation-layer arguments here
  (`type`, `mean`, `expressed_sex`, `residual_var`, etc.) — those belong
  in `define_phenotype()`. `overwrite = FALSE` (default) errors if the trait
  already exists; `overwrite = TRUE` replaces both the `trait_meta` row and its
  `index_meta` entry.
- `define_additive_effects()` — accepts a `tidybreed_table` from
  `get_table("genome_meta")` (optionally filtered) as its **first argument**.
  `trait_name` accepts a scalar **or vector** of trait names:
  - **Single trait** — manual (`effects` vector) or sampled (`distribution =
    "normal"/"gamma"`) with optional Falconer rescale. Reads `target_add_var`
    from `trait_var_comp`.
  - **Multiple traits** — draws correlated effects from `MVN(0, G)` via
    `MASS::mvrnorm`. `G = NULL` reads from `trait_var_comp`. `method =
    "shared"` (all traits use the filtered loci) or `"union"` (per-trait QTL
    sets from existing `genome_effects` rows, restricted to the filtered loci).

  Re-calling for the same trait replaces existing rows in `genome_effects`.

  ```r
  # Single trait
  pop |> get_table("genome_meta") |> filter(chr %in% 1:5) |> define_additive_effects("ADG")

  # Multiple correlated traits (shared QTL set)
  G <- matrix(c(0.25, 0.10, 0.10, 0.30), 2, 2,
              dimnames = list(c("ADG", "BW"), c("ADG", "BW")))
  pop |> get_table("genome_meta") |> filter(chr %in% 1:5) |>
    define_additive_effects(c("ADG", "BW"), G = G)

  # Use generation-0 animals to define base allele frequencies
  gen0 <- get_table(pop, "ind_meta") |> filter(gen == 0L)
  pop |> get_table("genome_meta") |> filter(...) |>
    define_additive_effects("ADG", base = "current_pop", base_tbl = gen0)
  ```

### `define_effect_cov_matrix()` / `define_effect_random()` / `define_effect_fixed_class()` / `define_effect_fixed_cov()` / `define_effect_intercept()`

`R/define_effect_cov_matrix.R`, `R/define_effect_random.R`, `R/define_effect_fixed_class.R`, `R/define_effect_fixed_cov.R`, `R/define_effect_intercept.R`

- `define_effect_cov_matrix(pop, effect_name, cov_matrix)` — **single entry
  point for all variance/covariance data**. Routes by `effect_name`:
  genetic effects (`"gen_add"`, `"dominance"`, `"epistasis"`) → `trait_var_comp`;
  `"residual"` → `define_residual_cov()` → `phenotype_var_comp`;
  any other name → `phenotype_var_comp` with that `effect_name`.
  Can be called before `define_trait()` or `define_effect_random()`.
- `define_effect_random()` — `variance` optional if already in `phenotype_var_comp`.
- `define_effect_fixed_class()` — discrete level → shift mapping.
- `define_effect_fixed_cov()` — linear regression term (`slope * (x - center)`).
- `define_effect_intercept()` — sets intercept (`target_add_mean`) for a trait.

### `define_phenotype()` / `define_residual_cov()`

`R/define_phenotype.R`, `R/define_residual_cov.R`

- `define_phenotype(pop, phenotype_name, type, mean, expressed_sex, repeatable, ...)` —
  registers an observed phenotype in `phenotype_meta`. For simple traits
  `phenotype_name` matches the `trait_name` already in `trait_meta`. For
  composite phenotypes (e.g. weaning weight, SGE ADG) the name is new and no
  prior `define_trait()` call is needed for it.

  Key arguments:
  - `residual_var` — scalar; writes one unconditional diagonal entry to
    `phenotype_var_comp` (with `effect_name = 'residual'`). For correlated or
    heterogeneous residuals use `define_residual_cov()` afterwards.
  - `components` — data frame with columns `source_trait_name` and
    `contributor_type` (`"self"`, `"dam"`, `"sire"`, `"group"`). Optional
    columns: `weight`, `weight_type`, `aggregation`, `group_column`,
    `group_table`, `covariate_name`, etc. Writes to `phenotype_components`.
    `NULL` (default) = simple single-self trait.
  - `missing_component_action` — `"skip"` (default) or `"error"`. Stored in
    `phenotype_meta` and applied uniformly by `add_phenotype()` for **any**
    missing composite piece (missing group assignment, missing dam/sire TBV,
    etc.). `"skip"` excludes the individual and warns with a count + up to 5
    example IDs. `"error"` stops immediately.

- `define_residual_cov(pop, phenotype_names, cov_matrix, condition_column = NULL, ...)` —
  writes conditional or unconditional residual (co)variance entries to
  `phenotype_var_comp` (always with `effect_name = 'residual'`). Supply a named
  matrix for multi-phenotype correlated residuals, or call once per sex/group
  level with `condition_column = "sex"` and `condition_level = "M"` / `"F"` for
  heterogeneous residuals.

### `add_phenotype()` / `add_tbv()`

`R/add_phenotype.R`, `R/add_tbv.R`

Both functions accept a `tidybreed_table` (from `get_table()` + optional
`filter()`) as their first argument and return `tidybreed_pop`.

- `add_phenotype()` — the workhorse. `phenotype_name` (formerly `trait_name`)
  defaults to all phenotypes in `phenotype_meta` when omitted. Internally calls
  `add_tbv()` first for all required source traits (including all composite
  components and group-member IDs), then assembles the composite TBV via
  `.assemble_composite_tbv()`. For group contributors (SGE model), all
  group-member TBVs are pre-fetched so group aggregation is a single in-memory
  pass. Adds fixed/random covariate contributions, samples residuals (joint
  `MVN(0, R)` when multiple phenotypes share the subset and `R` is stored;
  otherwise independent). Converts liability to phenotype per `type`.
  Writes `ind_phenotype` rows and updates `ind_tbv`.
- `add_tbv()` — TBV-only; no phenotype records. Computes centered TBV by joining
  `ind_haplotype` to `genome_effects` (`genome_effect_type = "additive"`) on
  `(locus_name, line_origin)`: a line-specific effect row (`genome_effects.line_name
  = ind_haplotype.line_origin`) is preferred, falling back per-locus to the
  population-wide row (`line_name IS NULL`) only when no line-specific row exists
  for that locus/line. This is what makes crossbreeding TBV correct (e.g. a Duroc ×
  Landrace F1 centered against each parent line's own QTL effects and base allele
  frequency). Imprinted traits (`expressed_parent = "parent_1"`/`"parent_2"`)
  restrict the join to that parent's `parent_origin` before the same line-matching
  logic applies. `trait_name` also defaults to all traits in `trait_meta` when
  omitted. Optional arguments for true index computation:
  - `index_names` — character vector of named indices; when supplied, multiplies
    per-trait TBVs by the index weights and writes results to `ind_true_index`.
    `NULL` (default) skips index computation.
  - `type` — `"index"` (default, uses `index_weight`), `"economic"` (uses
    `economic_weight`), or `"both"` (writes two rows per individual distinguished
    by `weight_type`).
  - `overwrite_index = FALSE` — when `FALSE`, skips individuals that already have
    a value in `ind_true_index` for the given `(index_name, weight_type)`. Set
    `TRUE` to recompute (e.g. after updating index weights).

### `define_trait_simple()`

`R/define_trait_simple.R`

Convenience wrapper that chains `define_trait()` and `define_additive_effects()` for
a single uncorrelated trait. QTL are always placed randomly (n = `n_qtl`).
For non-random QTL placement or correlated multi-trait effects, use the
functions individually with `get_table("genome_meta") |> filter(...)  |>
define_additive_effects()`.

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

## Roadmap

### Longer-Term

- `select_parents()` — selection index or truncation selection
- Export: PLINK `.bed/.bim/.fam`, VCF
- Visualization helpers
- Dominance and epistasis effects (currently only additive)

## Future Compiled Code Policy

When adding C++ code, prefer standard Rcpp/Rcpp Attributes and CRAN-compatible
dependencies that install cleanly on Linux, Windows, and macOS. Avoid required
system libraries, non-portable compiler flags, and architecture-specific code in
the main path.

Architecture-specific acceleration such as CUDA, Metal, or platform-specific
SIMD should be isolated behind optional backends with feature detection and an
R/Rcpp fallback. Do not make CUDA, Metal, or another accelerator required for
installation unless the package design deliberately splits those backends later.

## License

`tidybreed` is licensed under the **MIT License**. Use `License: MIT + file LICENSE`
in `DESCRIPTION`. The `LICENSE` file must contain the standard MIT license text with
the year and copyright holder.

## Versioning Policy

Use **three-part semantic versioning**: `MAJOR.MINOR.PATCH` (e.g. `0.0.1`).
Do **not** use the four-part devtools convention (`0.0.0.9000`).

**Pre-1.0.0 (current phase): version bumps are bookkeeping, not compatibility
promises.** Per "Development Status: Pre-1.0.0 — Break Freely" above, breaking
schema/API changes are normal development work. Document meaningful changes in
`NEWS.md`, but do not add compatibility layers. Reserve the `0 → 1` MAJOR bump
for the deliberate `1.0.0` stabilization. The table below is the policy that
takes effect **at and after 1.0.0**:

| Part  | Bump when… (≥ 1.0.0)                                 |
|-------|------------------------------------------------------|
| PATCH | Bug fixes, doc updates, minor internal changes       |
| MINOR | New exported functions or non-breaking feature additions |
| MAJOR | Breaking API changes                                 |

**Before every commit + push, update:**
1. `DESCRIPTION` — `Version:` field
2. `NEWS.md` — add an entry under the new version heading

## Key File Locations

| Path                       | Contents                          |
|----------------------------|-----------------------------------|
| `R/`                       | All exported functions            |
| `tests/testthat/`          | Formal testthat test suite        |
| `tests/test_*.R`           | Manual/dev test scripts           |
| `man/`                     | Roxygen-generated documentation   |

## Design Rationale

**Why DuckDB?** Columnar, embedded (no server), SQL via dbplyr, handles
datasets larger than RAM, excellent R integration.

**Why 0/1/2 encoding?** Standard in genomics (PLINK, VCF). Easy to interpret
(count of alternate allele). Efficient integer storage.

**Why store both haplotypes and genotypes?** Haplotypes are required for
recombination and phased exports. Genotypes are required for GWAS and genomic
prediction. Computing genotypes on the fly during every query would be wasteful.

## Development Workflow

For focused changes, prefer `pkgload::load_all(".", quiet = TRUE)` plus targeted
`testthat::test_file()` calls over running the full suite first. Run broader
tests before declaring shared schema, RNG, or cross-module changes done.

Performance work must start from a small reproducible benchmark or profiling
script under `dev/benchmarks/`. Keep benchmarks deterministic, small enough to
run during development, and scalable enough to expose the intended bottleneck.

### Test Coverage (`covr`)

`covr` is in `Suggests` and is run **locally only** — there is no
`test-coverage.yaml` workflow and no Codecov integration or badge. Coverage is a
diagnostic for finding untested and dead code, not a metric to publish or chase.

Coverage measures which lines *execute* during the suite, not whether anything is
asserted about them. Treat a low number as reliable bad news and a high number as
weak good news.

**The invocation that works** (both deviations from the obvious call are load-bearing):

```r
cov <- covr::package_coverage(
  path = ".",
  type = "none",
  code = 'testthat::test_dir("tests/testthat", package = "tidybreed",
                             load_package = "installed", reporter = "summary")'
)
covr::report(cov)   # interactive HTML, uncovered lines in red
```

- **Keep `type = "none"` and scope explicitly to `tests/testthat`.** `covr`
  installs from source (`R CMD INSTALL` ignores `.Rbuildignore`), so anything
  parked in `tests/` gets executed even when `R CMD check` would skip it. The
  `.Rbuildignore` entry `^tests/test_[^/]*\.R$` still guards that slot — do not
  add unassertive demo scripts there again. The formal suite in `tests/testthat/`
  is the only test surface.
- **Do not use `testthat::test_local()`.** It calls `pkgload::load_all()`
  internally, which replaces `covr`'s instrumented package and silently reports
  **0% for every R file** while the C++ still reads ~99% (gcov instrumentation
  lives in the `.so` and survives). This fails silently as a plausible-looking
  result, not as an error. `test_dir(load_package = "installed")` uses the
  instrumented install.

A full instrumented run takes **~12 minutes** (the suite is DuckDB file-backed
and instrumentation adds 2–5×). Do not put it in a fast edit-test loop.

**Baseline at v0.63.0** — 80.7% across R code; `src/make_gametes.cpp` 98.9%.
Remaining gaps:

| File | Coverage | Why |
|------|----------|-----|
| `blupf90_helpers.R` | 22% | Needs the external BLUPF90 binary |
| `add_ebv.R` | 29% | Same — external solver dependency |
| `define_effect_cov_matrix.R` | 54% | Routing branches per `effect_name` |
| `tidybreed-package.R` | 67% | Mostly `.onLoad`/startup paths |
| `restore_pop.R` | 68% | Reopen-and-resume paths |
| `add_phenotype.R` | 71% | Composite/SGE and distribution branches |

The BLUPF90 paths are expected to stay low without a CI solver; do not chase
them. `add_tbv.R` was the standing concern at 37% and is now 99.6% — see
`tests/testthat/test-add_tbv_index.R`, which covers the `index_names` /
`weight_type` block, and `test-add_tbv.R`, which covers the line-precedence
crossbreeding join.

## Development Environment

### Running R Commands

The R executable path is platform-specific. When running R or Rscript via the
Bash tool, use the appropriate path:

**Windows (Hendrix Genetics AVD):**
- Check if working directory contains "Hendrix"
- R: `"/c/Program Files/R/R-4.5.1/bin/x64/R.exe"`
- Rscript: `"/c/Program Files/R/R-4.5.1/bin/x64/Rscript.exe"`

**Mac/Linux:**
- Use standard shell commands: `R` or `r` and `Rscript`

Example test command:
```bash
"/c/Program Files/R/R-4.5.1/bin/x64/Rscript.exe" -e "print('Hello from R')"
```

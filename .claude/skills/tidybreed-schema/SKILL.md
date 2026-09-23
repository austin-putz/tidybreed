---
name: tidybreed-schema
description: Full tidybreed DuckDB schema reference — every table's columns, types, keys, reserved columns and invariants (genome_meta, genome_map, genome_effects*, ind_haplotype, chr_inheritance, phenotype_*, ind_tbv/tgv/ebv, index tables). Load before reading, querying, or changing any table, column, view, or DDL.
---

# tidybreed Database Schema

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

### `genome_effects` / `genome_effect_members` / `genome_effect_member_origins`

Genome effects are stored as **terms**, not one row per locus. A term is one
coefficient over one or more loci, each locus contributing a named basis function
(`contrast_name`), each optionally scoped to allele copies of a given line and/or
parent of origin. Created by `define_genome()` (not `open_pop()` — the `locus_id`
foreign key needs `genome_meta` to exist first). Written by
`define_genome_effects()` and, for the reserved `generated_additive_tbv` owner,
by `define_additive_effects()`.

**`genome_effects`** — one row per term.

| Column           | Type    | Notes                                                     |
|------------------|---------|-----------------------------------------------------------|
| id_genome_effect | INTEGER | Primary key assigned via `next_int_id()`                  |
| trait_name       | VARCHAR | R-enforced FK to `trait_meta.trait_name`                  |
| effect_owner     | VARCHAR | Which writer owns these rows, **for replacement only**. Owners always sum and are never selected between: `"generated_additive_tbv"` is reserved for `define_additive_effects()`, `"custom"` is the `define_genome_effects()` default |
| effect_name      | VARCHAR | Optional per-term label; no mathematical meaning          |
| genome_value     | DOUBLE  | The term's coefficient                                    |

**`genome_effect_members`** — one row per (term × locus), canonicalized by
ascending `locus_id` with `member_slot` running `1..n`.

| Column           | Type     | Notes                                                    |
|------------------|----------|----------------------------------------------------------|
| id_genome_effect | INTEGER  | FK to `genome_effects`; PK part                          |
| member_slot      | INTEGER  | Position in the term; PK part                            |
| locus_id         | INTEGER  | FK to `genome_meta.locus_id`                             |
| contrast_name    | VARCHAR  | `"additive"` (per allele copy), `"dominance"` (Cockerham, diploid), `"indicator"` (one genotype state) |
| copy_count_value | UTINYINT | Indicator state: realized copy count. Required with `dosage_value`, because dosage alone conflates "no copy", "one allele-0 copy" and "two allele-0 copies" |
| dosage_value     | UTINYINT | Indicator state: dosage of allele 1                      |
| center_value     | DOUBLE   | Per-copy centring constant: `p` under Cockerham coding, `0.5` under functional coding. Required for non-indicator contrasts |

**`genome_effect_member_origins`** — the scope of a member, as a predicate over
allele copies. **No rows = the common scope**, which matches every copy.

| Column           | Type     | Notes                                                    |
|------------------|----------|----------------------------------------------------------|
| id_genome_effect | INTEGER  | PK part; composite FK to `genome_effect_members`         |
| member_slot      | INTEGER  | PK part; composite FK to `genome_effect_members`         |
| origin_slot      | INTEGER  | PK part; canonicalized by sorting the origin tuple       |
| line_match_type  | VARCHAR  | `"exact"` (named line), `"unknown"` (copies with no line), `"any"` (additive members only, and must carry a `parent_origin`) |
| line_name        | VARCHAR  | Set only for `"exact"`                                   |
| parent_origin    | UTINYINT | 1 = sire, 2 = dam, NULL = either. Non-NULL on an additive member is how **imprinting** is expressed |
| copy_count       | INTEGER  | Copies demanded; always 1 on an additive member          |

**Reserved**: all columns of all three. Row deletion is refused — effect
definitions are configuration and are replaced through
`define_genome_effects(mode = ...)`, not row-deleted.

**No foreign keys *inside* the set** (members → effects, origins → members),
deliberately. DuckDB 1.5.5 refuses to delete a parent row inside an explicit
transaction whose children were deleted earlier in that same transaction — for
single-column and composite keys alike, in either delete order — which makes
every replace mode unwritable as one transaction. Since a half-replaced effect
model is a *different* model rather than a weaker one, the FKs go and
`validate_genome_effects()` reports orphans in both directions before every
`COMMIT` instead. The `locus_id` → `genome_meta` key **stays**: `genome_meta`
rows are never deleted, so it never sits in the failing position. See
`tests/testthat/test-genome-effects-schema.R` for the pinned DuckDB behaviour.

**Rules.** An additive member takes at most one origin row (`"A or B"` is
expanded into separate variants); a genotype member takes an exact multiset whose
`copy_count`s sum to the state's copy count (2 for `dominance`, the declared
`copy_count_value` for `indicator`). Terms sharing a **family signature** —
trait, owner, and the ordered member states, exposed as `family_key` on the
`genome_effect_terms` view — are scope variants of one term and **compete**: the
most specific matching scope wins, and overlapping-but-incomparable scopes are
refused at write time. Terms with different keys **sum**.

Causal-locus membership is **implicit**: a locus is causal for a trait if it
appears as a member of any term for it. No boolean flag is stored.

### Genome-effect views

| View                   | Grain                      | What it is for                                 |
|------------------------|----------------------------|-------------------------------------------------|
| `genome_effect_terms`  | one row per term           | `effect_order`, `contrast_signature`, `family_key`, `scope_description` — all derived, never stored. `family_key` is how you see which terms compete and which sum |
| `genome_effect_loci`   | one row per (term × locus) | `locus_name` joined from `genome_meta` and the term's `genome_value` repeated on every member row; the place to ask which loci are causal, and the only relation at locus grain that can be filtered by effect size. `locus_name` lives only here, so there is no id/name agreement invariant in the base tables. Never `SUM(genome_value)` — an interaction term would be counted once per member |

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
| locus_name    | VARCHAR  | FK to `genome_meta.locus_name`; denormalized so exports and user queries read without joining `genome_meta`. The effect tables key on `locus_id`, not on this column |
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
| condition_change_action  | VARCHAR | `"error"` (default) or `"independent"` — what to do when a correlated phenotype's stored residual was drawn under a different residual `condition_level` than the current record resolves to (see `plans/sample_correlated_effects.md` D2/D6). Must agree across every phenotype in one residual covariance block, so it is **block-scoped**: `define_phenotype()` sets it while the phenotype is still a block of one, and `define_condition_change_action()` changes it afterwards, on every member in one transaction |

**Reserved**: all columns (managed by `define_phenotype()`, except
`condition_change_action`, which `define_condition_change_action()` also
writes — block-scoped, one column, one transaction).

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
| weight_type        | VARCHAR | `"fixed"` (default) or `"covariate"`. Those are the only two implemented, and `define_phenotype()` rejects anything else |
| covariate_name     | VARCHAR | Covariate column when `weight_type = "covariate"`                  |
| covariate_table    | VARCHAR | Table containing `covariate_name`                                  |
| poly_order         | INTEGER | Polynomial basis order                                             |
| poly_scale_min/max | DOUBLE  | Legendre scaling bounds                                            |
| component_names    | VARCHAR | Comma-separated `ind_tgv.component_name` values this component draws from; default `"order1_additive"`. **Reserved** — `add_phenotype()` reads only the additive breeding value today. This is the one reserved column here, kept because its counterpart (`ind_tgv.component_name`) is already written by `add_tgv()` |

**Note on SGE (Social Genetic Effects / Bijma model)**: for `contributor_type = "group"`,
`add_phenotype()` aggregates group-mates' TBVs (excluding self). A singleton (no
group-mates) receives a social contribution of 0 and is not excluded. An individual
with no group assignment receives `NA` and is handled by `missing_component_action`.
`group_table` must have exactly one row per focal individual (error otherwise).
All contributor lookups — self, dam, sire, group, and `formula_tbv`'s
`dam()`/`sire()`/`group_sum()`/`group_mean()` — go through `R/contributor_tbv.R`
(`.tbv_by_id()`, `.group_mate_tbv()`, `.group_members()`), one registered-view
SQL each; ids never enter SQL text.

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
(phenotype × effect × level), written by `add_phenotype()` (Stage 3) and
deleted by `define_effect_random(overwrite = TRUE)` or `remove_rows()`. Pure
observation-layer noise — no genetic content. A level is a **persistent
entity**: its draw is realized the first time a planned record touches it
and reused by every later record with that level, in every later call. In a
covariance block (`define_effect_cov_matrix(effect_name, …)`) a level's draw
for one phenotype is drawn conditional on the draws it already has stored
for the block's other phenotypes; a member with no draw stays latent (no
row) until a record needs it.

| Column         | Type    | Notes                                              |
|----------------|---------|-----------------------------------------------------|
| phenotype_name | VARCHAR | FK to `phenotype_meta.phenotype_name`; PK part     |
| effect_name    | VARCHAR | FK to `phenotype_effects.effect_name`; PK part     |
| level          | VARCHAR | Level of the grouping column; PK part              |
| draw_value     | DOUBLE  | Sampled deviation for that level                   |
| date_sampled   | DATE    | When the draw was taken                            |

**Primary key**: `(phenotype_name, effect_name, level)`.

**Reserved**: all columns (written only by `add_phenotype()`).

### `phenotype_var_comp`

Phenotype-layer variance components. One row per (effect_name, phenotype₁, phenotype₂,
condition). Both (i,j) and (j,i) pairs stored for off-diagonal entries. Populated by
`define_phenotype(..., residual_var = ...)`, `define_residual_cov()`, and
`define_effect_random()`.

`effect_name = 'residual'` is reserved for residual noise. All named random effects
(HYS, litter, pen, etc.) use their own string (e.g. `'hys'`, `'litter'`).
The `condition_column` / `condition_level` columns are used only for `'residual'`
to model heterogeneous residual variance by sex, group, etc.

**Covariance blocks.** For one `effect_name`, the phenotypes joined by any stored
pair row (an explicit `0` counts) form a *block*, and every writer goes through
`write_phenotype_cov_block()` / `validate_phenotype_cov_block()` in
`R/phenotype_cov_block.R` inside one transaction:

- A block is **declared in one call, as a complete matrix** (D1). A call naming
  a fragment or a strict subset of an existing block is an error naming the
  omitted phenotypes; the matrix must be symmetric, finite and PSD.
- Conditional `'residual'` rows form **strata** `(condition_table,
  condition_column, condition_level)` of the same block: every stratum names the
  same phenotypes and a block has one condition column. Unconditional rows have
  `condition_column`, `condition_table` and `condition_level` all `NULL`.
- A block is **locked once realized** (D3): any `ind_phenotype` row of a member
  with `residual_value IS NOT NULL`, or any `phenotype_random_effects` row for
  `(effect_name, member)`. The error gives the `remove_rows()` call that clears
  the realizations. There is no `force`.
- Every defined member of a residual block carries the same
  `phenotype_meta.condition_change_action` (D6); in a named-effect block of two
  or more, every `phenotype_effects` row is `random`, `normal`, and reads the
  same `(source_column, source_table)`.

Stored matrices are exactly symmetric (the writer stores `(M + t(M)) / 2`).
`find_covariance_blocks(conn, effect_name, phenotype_names)` in
`R/correlated_draws.R` is the one reader that turns stored rows back into
matrices — one entry per block touching the targets, every stratum assembled,
invariants re-checked — and `resolve_correlated_draws()` beside it is the pure
conditional-MVN sampler the phenotype layer will draw through.

See `plans/sample_correlated_effects.md` §5.2, §5.4, §5.9 and D1/D3/D5/D6.

| Column               | Type    | Notes                                                              |
|----------------------|---------|--------------------------------------------------------------------|
| id_phenotype_var_comp| INTEGER | Primary key assigned by tidybreed via `next_int_id()`              |
| effect_name          | VARCHAR | `'residual'` or any named random effect (e.g. `'hys'`, `'litter'`)|
| phenotype_name_1     | VARCHAR |                                                                    |
| phenotype_name_2     | VARCHAR |                                                                    |
| cov_value            | DOUBLE  |                                                                    |
| condition_column     | VARCHAR | NULL = unconditional stratum; used only for `effect_name = 'residual'` |
| condition_table      | VARCHAR | Table holding `condition_column`; `NULL` on unconditional rows      |
| condition_level      | VARCHAR | Value of `condition_column` for this stratum; `NULL` on unconditional rows |
| weight_type          | VARCHAR | Default `"fixed"`                                                  |
| poly_order           | INTEGER | Polynomial order for `"legendre"` weight type                      |

### `ind_phenotype`

Phenotype records in long format. Populated by `add_phenotype()`.

| Column                   | Type    | Notes                                             |
|--------------------------|---------|---------------------------------------------------|
| id_phenotype             | INTEGER | Primary key assigned by tidybreed via `next_int_id()` |
| id_ind                   | VARCHAR |                                                   |
| phenotype_name           | VARCHAR | FK to `phenotype_meta.phenotype_name`             |
| pheno_value              | DOUBLE  | Phenotype value                                   |
| pheno_number             | INTEGER | 1 = first record for this individual × trait, etc. Ordinal identity, **not** simulated time |
| liability_value          | DOUBLE  | Raw liability for categorical phenotypes with `store_liability = TRUE`; NULL otherwise |
| cat_name                 | VARCHAR | Category label for categorical phenotypes defined with `cat_names`; NULL otherwise |
| residual_value           | DOUBLE  | Realized **liability-scale** residual for model-generated and `user_residual` records; NULL for `user_values` / `derived_formula` records. Conditions later draws of correlated phenotypes (`plans/sample_correlated_effects.md`) |
| residual_condition_level | VARCHAR | `condition_level` of the residual (co)variance stratum the residual was drawn under; NULL when the unconditional `R` was used (the *selected* stratum, not the raw column value) |
| *user cols*              | any     | Added via `mutate_table()` or scalar `...` in `add_phenotype()` |

All nine columns are in the base `CREATE TABLE` (in `ensure_trait_tables()`);
nothing is added by on-demand `ALTER TABLE`. **Reserved**: all nine.

### `ind_tbv`

True breeding values (simulation ground truth). Populated by
`add_phenotype()` and `add_tbv()`. Logical key `(id_ind, trait_name)` unique.

| Column     | Type    | Notes                                  |
|------------|---------|----------------------------------------|
| id_tbv     | INTEGER | Primary key assigned by tidybreed via `next_int_id()` |
| id_ind     | VARCHAR |                                        |
| trait_name | VARCHAR |                                        |
| tbv_value  | DOUBLE  |                                        |

### `ind_tgv`

True **genetic** values (simulation ground truth): the total genotypic value,
split by declared model structure. One row per (individual × trait × component).
Populated by `add_tgv()`. Created in `define_trait()`'s lazy DDL block beside
`ind_tbv`.

| Column         | Type    | Notes                                              |
|----------------|---------|-----------------------------------------------------|
| id_tgv         | INTEGER | Primary key assigned via `next_int_id()`            |
| id_ind         | VARCHAR |                                                     |
| trait_name     | VARCHAR |                                                     |
| component_name | VARCHAR | `"order1_additive"`, `"order1_dominance"`, `"order1_other"` (a hand-entered order-1 indicator surface), or `"interaction"` (any term with ≥ 2 members) |
| tgv_value      | DOUBLE  | Raw sum of the contributing terms; **no mean is added** |

**Unique**: `(id_ind, trait_name, component_name)`. **Reserved**: all columns.

`component_name` records **how a term was declared, not a variance component**.
A functional A×A term contributes to A, D *and* I in the statistical sense; the
names carry the order precisely so they cannot be misread as `V_A` / `V_D` /
`V_I`. There is deliberately **no** `replicate` column: like `ind_tbv`, that
column exists only in the archive copy, and `archive_replicate()` refuses to
stamp a table that already has one.

The total is the derived view **`ind_tgv_total`** (`id_ind`, `trait_name`,
`tgv_total`), never a stored `'total'` row — a stored total would make every
`SUM(tgv_value)` double-count.

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

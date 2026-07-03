# Critique: `plans/refactor_haplotype.md`

## Executive summary

The proposed move from wide haplotype/genotype tables to long tables is directionally correct. It directly addresses several real limitations in the current codebase: loci are currently hard-coded as `locus_1 ... locus_n` columns, `add_offspring()` assumes exactly two haplotype rows per individual, `extract_genotypes()` and BLUPF90 export depend on wide dosage rows, and `add_tbv()` still pulls whole genotype matrices into R.

However, the current plan is too broad to implement safely as one refactor. It combines at least six design changes:

1. Wide-to-long storage.
2. Table renaming from `genome_*` to `ind_*`.
3. On-demand dosage caching.
4. Allele line-of-origin tracking for crossbreeding TBV.
5. Sex chromosome and polyploid inheritance rules.
6. Future sparse mutation storage and future Rcpp acceleration.

Those changes are related, but they do not all need to ship together. The most important recommendation is to split the refactor into a minimal, correct, diploid long-format implementation first. Add sex chromosomes, polyploid meiosis, mutation sparsity, and Rcpp only after the core long schema has tests proving parity with the current behavior.

## Highest-risk issues

### 1. `parent_origin` and `strand` semantics are under-specified for meiosis

This is the largest correctness risk in the plan.

In the current diploid model, a parent has two haplotype rows:

- `parent_origin = 1`: the haplotype that individual inherited from its sire.
- `parent_origin = 2`: the haplotype that individual inherited from its dam.

When producing a gamete from that parent, recombination must operate across both of that parent's homologs. For example, a sire's gamete is made by recombining across the sire's own `parent_origin = 1` and `parent_origin = 2` rows. The resulting offspring row is then written with `parent_origin = 1` because it came from the sire.

The plan sometimes reads as if `parent_origin` can be used to fetch "the sire contribution" or "the dam contribution" during gamete generation. That is not correct. Before the offspring exists, `parent_origin` in the parent's rows is grandparental origin, not mating-side origin.

The plan should explicitly separate these concepts:

- Within an existing individual, `parent_origin` describes the origin of that individual's chromosome copy.
- In a newly written offspring row, `parent_origin` describes which mating parent contributed the gamete.
- Gamete generation for a parent must load all chromosome copies carried by that parent for the chromosome, not only rows with one `parent_origin`.

This becomes more important for polyploids. A tetraploid parent has four homologs. Labeling them as `parent_origin in {1,2}` x `strand in {1,2}` is a reasonable storage scheme, but the meiosis algorithm must know which homologs can pair and how many gamete copies to emit. It cannot simply recombine within one `parent_origin` group.

Recommendation: for v1, keep the schema columns `(parent_origin, strand)` if desired, but implement and test only the diploid case where each individual has exactly two rows per locus. Document polyploid support as schema-ready but algorithm-pending.

### 2. The plan mixes dense and sparse haplotype semantics

The proposed schema says founding loci are stored densely, while novel mutations are sparse: absence of a row means allele `0`. Later the plan rejects sparse haplotype storage for v1 and says dense storage is the adopted design.

This contradiction matters because several proposed SQL snippets assume dense rows:

- `add_dosage()` aggregates from `ind_haplotype` and will not emit zero-dosage rows for requested individual x locus combinations that are absent.
- `extract_genotypes()` via `PIVOT` will produce missing values, missing columns, or incomplete rows for sparse mutation loci unless it cross joins target individuals and target loci first.
- Allele frequency denominators cannot be derived from row counts once any locus is sparse.
- Parent haplotype loading for recombination must materialize implicit zeros for loci absent from one or both parents.

Recommendation: choose one of these for v1:

1. Fully dense `ind_haplotype` for all stored loci, including future mutations with backfilled zeros. This is simple but expensive for mutations.
2. Fully sparse `ind_haplotype` with explicit cross-join/coalesce logic everywhere dosage or export needs a dense matrix.
3. Hybrid dense founding loci plus sparse introduced loci, but with an `is_sparse` or `introduced_gen`-driven query convention and tests proving every downstream query handles both cases.

The current plan describes option 3, but most implementation details are written as if option 1 is true. Do not start coding until this is resolved.

### 3. `line_origin` should not be a foreign key to `genome_effects.line_name`

The proposed `line_origin` column is biologically useful, but the FK target is wrong. `genome_effects.line_name` is not a line dimension table. It may contain zero rows for a line until effects are defined, may contain multiple rows per line per trait/locus, and may legitimately have `NULL` for population-wide effects.

`line_origin` should mean "the founding genetic line this allele traces to." That meaning is independent of whether any line-specific effect has been defined.

Recommendation:

- Do not FK `ind_haplotype.line_origin` to `genome_effects.line_name`.
- Either leave it unconstrained, FK it to a future `line_meta(line_name)`, or use `ind_meta.line_name` values by convention.
- Always set `line_origin` for founders when a line is known.
- Reserve `NULL` for genuinely unknown/untracked origin, not for "population-wide effects mode."

### 4. Crossbreeding TBV SQL omits centering and trait rules

The proposed long-format TBV SQL computes `SUM(effect * allele)`, but current `add_tbv()` centers additive TBV using `base_allele_freq`:

- Non-imprinted current formula: `genotype %*% a - 2 * sum(p_base * a)`.
- Imprinted current formula: `haplotype %*% a - sum(p_base * a)`.

The critique is not that the new SQL cannot do this. It can. The issue is that the plan does not specify the centering query, especially for line-specific effects. If `genome_effects` has per-line effects, each allele's centering constant probably needs to use the matching line-specific `base_allele_freq`, not a population-wide fallback.

The new TBV implementation also needs to preserve these current behaviors:

- filtered `tidybreed_table` subsets;
- `trait_meta.expressed_parent` for parent-specific expression;
- `trait_meta.expressed_sex` indirectly through `add_phenotype()` workflows;
- existing `ind_tbv` upsert semantics;
- index calculations in `add_tbv(index_names = ...)`;
- formula/composite phenotype code paths that depend on TBVs.

Recommendation: add an explicit TBV design section with centered SQL for:

- population-wide additive effects;
- line-specific additive effects with fallback;
- parent-specific expression;
- dense and sparse haplotype cases, if sparse remains in scope.

### 5. Founder haplotypes lose required line-pool information

The proposed long `founder_haplotypes` table has:

```text
haplotype_id, locus_name, allele
```

The current implementation supports line-specific founder haplotype pools using `founder_haplotypes.line_name`, and `add_founders()` is intentionally piped from a filtered `get_table("founder_haplotypes")`. Removing `line_name` would break multi-line setup and crossbreeding simulations.

Recommendation: long founder haplotypes should include `line_name` or a pool identifier:

```text
line_name      VARCHAR NULL
haplotype_id   VARCHAR or INTEGER
locus_name     VARCHAR or INTEGER
allele         UTINYINT
```

Use a key like `(line_name, haplotype_id, locus_name)` if haplotype IDs are only unique within line, or a globally unique `haplotype_id` plus a separate `line_name` column.

Also decide what replaces the current `genome_meta.founder_allele_freq` behavior. Today it is "last founder pool wins", which is already weak for multi-line simulations. The refactor is a good chance to replace it with a long table such as:

```text
founder_allele_freq(line_name, locus_name, allele_freq)
```

or to compute base frequencies directly from the selected founder/current population in SQL.

### 6. Table renaming needs a compatibility decision

The plan proposes `ind_haplotype` and `ind_genotype`, while the current package and tests use `genome_haplotype` and `genome_genotype` in many places:

- `define_genome()`;
- `add_founders()`;
- `add_offspring()`;
- `add_tbv()`;
- `extract_genotypes()`;
- `blupf90_helpers.R`;
- `remove_rows()`;
- `summary_pop()`;
- `schema.R`;
- `sql_utils.R`;
- user-facing examples and tests.

The new names are arguably clearer because the rows are individual-level, but a rename is a breaking change. It also increases the size of the refactor.

Recommendation: choose one of these before implementation:

1. Keep table names `genome_haplotype` and `genome_genotype`, but change their shape to long. This minimizes API breakage.
2. Rename to `ind_haplotype` and `ind_genotype`, and provide compatibility views or clear migration errors for old names.

If you choose option 2, include views only if they are cheap and semantically honest. A view that pretends long data is wide will be expensive and fragile for large genomes.

### 7. Sex chromosome support is not yet complete enough to code safely

The `chr_meta` idea is good, but the proposed columns are not enough to implement all listed cases without hidden assumptions:

```text
chr_name, copies_M, copies_F, hemi_parent, recombines
```

Problems to resolve:

- The text says `copies_M` and `copies_F` are `0, 1, or 2`, then later uses `4, 6, 8` for polyploids.
- `copies / 2` only works for even-copy biparental inheritance. It does not apply to `copies = 1`, and it may not apply to odd ploidies.
- Organellar inheritance with `copies_M = 1, copies_F = 1, hemi_parent = parent_2` means both male and female offspring get one copy from the dam, but the code path must not try to recombine diploid-like homologs.
- Pseudoautosomal regions, no-recombination regions, and chromosome segments with different rules are not representable if rules are only chromosome-level.
- Standard GBLUP dosage coding and TBV centering become sex- and chromosome-copy dependent.

Recommendation: do not include sex chromosomes in the first implementation unless the immediate goal requires them. For the first refactor, create `chr_meta` with default diploid autosome rows only, but keep `define_chr()` and non-autosomal inheritance as a follow-up.

### 8. Performance estimates need benchmarking before the schema is locked

The storage and performance section makes plausible claims, but they are not yet reliable enough to drive the design.

For 5,000 individuals x 50,000 loci x 2 haplotypes, `ind_haplotype` is 500 million rows. Even if DuckDB compresses the columns well, this changes:

- insert volume;
- checkpoint/write-ahead-log behavior;
- query planning;
- R memory pressure when collecting rows;
- `PIVOT` cost for genotype export;
- `remove_rows()` cost;
- archive/reset workflows.

The statement that `locus_name VARCHAR` has little penalty due to dictionary encoding may be true on disk, but joins and pivots still pay for string keys. A storage key of `locus_id INTEGER` is worth benchmarking. The convenience of joining directly to `genome_effects.locus_name` should not be the only reason to choose string keys for the hottest table.

Recommendation: before implementing the full package refactor, write a standalone DuckDB benchmark script comparing:

- wide current format;
- long with `locus_name VARCHAR`;
- long with `locus_id INTEGER`;
- dense long vs sparse or hybrid;
- `GROUP BY` TBV on QTL-only loci;
- `PIVOT` export for chip loci.

Use at least one realistic scale, such as 5,000 individuals x 50,000 loci, and one CI-friendly scale for regression tests.

## Detailed design comments

### Schema

The proposed primary key `(id_ind, parent_origin, strand, locus_name)` is reasonable for dense haplotypes. It is not sufficient to describe sparse mutation absence unless the rest of the system consistently treats missing rows as allele zero.

If `locus_name` remains the key, ensure `genome_meta.locus_name` is unique and validated. Current `define_genome()` accepts custom names but does not appear to enforce uniqueness. A duplicated `locus_name` would corrupt joins in the long design.

Consider using `locus_id` as the internal key in hot tables and keeping `locus_name` for user-facing filters. That would require either adding `locus_id` to `genome_effects` or joining through `genome_meta`, but it gives stable integer keys, easier sort order, and safer pivot/export behavior.

If `ind_genotype` is a cache, its partial nature needs to be obvious in the table description. A partially populated cache is easy for users and internal code to misuse.

`INSERT OR REPLACE` on `ind_genotype` depends on a declared primary key or unique constraint. The DDL must create that constraint; otherwise idempotence is only aspirational.

### API surface

`add_dosage()` is a good idea, but its relationship with existing `add_genotypes()` must be clear:

- `add_genotypes()` currently means "mark animals as physically genotyped on a chip" by writing `has_<chip>` to `ind_meta`.
- `add_dosage()` would mean "materialize genotype dosage rows" in `ind_genotype`.

Those names are close enough that users may confuse them. The documentation should explicitly contrast them. Alternatively, consider `materialize_dosage()` if this is intended as an internal/cache operation.

The plan says `add_dosage()` accepts any table with `id_ind`. That is consistent with `add_tbv()` and `add_phenotype()`, but the implementation must use distinct IDs and must define what happens when the input table has multiple rows per individual, such as `ind_phenotype`.

The plan recommends no `overwrite` argument because dosage is idempotent. That is true for fixed haplotypes, but cache scope still matters. If a user materializes all loci and later materializes one chip, should old non-chip rows remain? Usually yes, but the docs need to say `ind_genotype` accumulates rows. A `delete_prior` option may still be useful for cache management.

### TBV and effects

Line-specific fallback should be implemented carefully. The proposed `NOT EXISTS` fallback is conceptually right, but there are edge cases:

- What if a locus has line-specific effects for one line but not another?
- What if `h.line_origin` is `NULL`?
- What if both a line-specific and population-wide effect row exist for the same line/locus/trait?
- What if `base_allele_freq` differs by line?

These should become explicit tests.

Current `define_additive_effects()` supports `base = "founder_haplotypes"` and `base = "current_pop"`. Both helper paths read wide tables today. The plan should list `R/define_additive_effects.R` helpers as a major rewrite target, not just `add_tbv.R`.

### Genotype export

`extract_genotypes()` and BLUPF90 helpers should not depend on `ind_genotype` if the cache is optional. That is a good design.

But the `PIVOT` query must guarantee:

- one row for every requested individual;
- one column for every requested locus;
- zeros for missing sparse haplotype rows;
- stable column order by `genome_meta.locus_id`;
- safe handling of arbitrary `locus_name` values;
- ploidy-aware dosage for non-diploid or sex-linked loci.

The simple `PIVOT` in the plan does not yet guarantee all of these.

### Migration

The proposed migration uses `tidyr::pivot_longer()` on entire wide tables. That will not be viable for large databases.

Migration should be transactional and SQL-first:

1. Create new long tables with temporary names.
2. Use DuckDB `UNPIVOT` or chunked SQL to populate them.
3. Validate row counts and checksums against the old tables.
4. Rename old tables to backup names.
5. Rename new tables into place.
6. Update `_schema_meta`, `pop$tables`, `TABLE_ROW_KEYS` assumptions, and any archive metadata.
7. Drop backups only after explicit confirmation or a successful second validation.

Migration also needs to handle:

- `founder_haplotypes`;
- chip columns in `genome_meta`;
- `has_<chip>` columns in `ind_meta`;
- existing `genome_effects`;
- archive/reset workflows;
- old databases that do not have line-aware founder haplotypes.

### Rcpp section

The Rcpp plan is useful as future direction, but it should not be part of the first refactor acceptance criteria. Linking directly against DuckDB from Rcpp and sharing connections safely is not a small detail. The R implementation should establish correctness first, then the C++ backend can target the same behavior.

Recommendation: move Rcpp details to a separate plan or appendix marked "future optimization."

## Recommended staged implementation

### Stage 1: Long diploid core, no new biology

Goal: preserve current behavior exactly while changing storage shape.

Scope:

- Keep diploid autosomes only.
- Keep current founder generation methods.
- Keep `add_founders()` and `add_offspring()` behavior equivalent to current tests.
- Store haplotypes long.
- Either keep table names or add compatibility handling.
- Rewrite `get_genotype_matrix()`, `get_haplotype_matrix()`, `add_tbv()`, `extract_genotypes()`, and BLUPF90 export against the long source of truth.
- Do not implement sparse mutation semantics yet.
- Do not implement polyploid or sex chromosome behavior yet.

Exit criteria:

- Current tests pass after being updated for expected table shape.
- A small deterministic simulation produces identical haplotypes, genotypes, TBVs, and exports before and after the refactor when using the same seed.
- A benchmark exists for insert, TBV, and export.

### Stage 2: Line-origin TBV

Goal: enable crossbreeding additive TBV by allele line of origin.

Scope:

- Add and populate `line_origin`.
- Update recombination to return both allele and `line_origin` vectors.
- Add centered SQL for population-wide and line-specific additive effects.
- Add tests for F1, F2, and backcross animals.

Exit criteria:

- A hand-calculated F1 TBV test passes.
- An F2 recombination test proves `line_origin` follows the inherited segment.
- Fallback from line-specific to population-wide effects is deterministic and documented.

### Stage 3: Dosage cache and export hardening

Goal: make marker-assisted selection and external software export efficient without making cache population a hidden requirement.

Scope:

- Add `add_dosage()` or `materialize_dosage()`.
- Treat `ind_genotype` as a documented partial cache.
- Ensure `extract_genotypes()` and BLUPF90 export read from haplotypes directly or have a clear cache fallback.

Exit criteria:

- Cache and direct extraction produce identical matrices.
- Repeated cache materialization is idempotent.
- Partial cache contents cannot silently change TBV or phenotype results.

### Stage 4: Non-diploid inheritance

Goal: add `chr_meta`, `define_chr()`, sex chromosomes, organelles, and later polyploid rules.

Scope:

- Start with sex chromosome and organellar cases.
- Defer autopolyploid/allopolyploid pairing until a clear meiosis model exists.
- Make TBV centering and dosage export ploidy-aware.

Exit criteria:

- X/Y or Z/W hand tests pass.
- Absent chromosomes produce no haplotype rows and no accidental dosage.
- Export either handles sex-linked loci correctly or warns/errors when used for standard GBLUP.

### Stage 5: Mutation sparsity

Goal: support introduced loci without backfilling every existing individual.

Scope:

- Decide fully sparse vs hybrid dense/sparse.
- Add query helpers that always coalesce missing alleles to zero.
- Update dosage, TBV, allele frequency, export, and recombination parent loading.

Exit criteria:

- A mutation in one individual can be inherited by descendants.
- Non-carriers export as dosage zero.
- TBV centering remains correct for carriers and non-carriers.

## Concrete edits suggested for `refactor_haplotype.md`

1. Add a "Non-goals for v1" section excluding polyploid meiosis, sex chromosomes, sparse mutation storage, and Rcpp unless you intentionally want them in the first implementation.
2. Clarify whether the table names will stay `genome_haplotype`/`genome_genotype` or become `ind_haplotype`/`ind_genotype`.
3. Fix the `parent_origin` explanation so gamete generation loads all homologs carried by a parent.
4. Add `line_name` or `pool_name` back to `founder_haplotypes`.
5. Remove the FK from `line_origin` to `genome_effects.line_name`.
6. Add explicit centered TBV SQL.
7. Decide dense, sparse, or hybrid haplotype semantics before describing `add_mutation()`.
8. Replace the R `pivot_longer()` migration sketch with a SQL/chunked migration plan.
9. Add `define_additive_effects.R`, `phenotype_helpers.R`, and `founder_haplotype_helpers.R` to the file change table.
10. Add a benchmark plan before accepting storage and performance assumptions.

## Questions to answer before coding

These are design decisions rather than blockers for this critique:

1. Is preserving the existing table names more important than the cleaner `ind_*` naming?
2. Should v1 be diploid-only, with `chr_meta` created but non-default rules disabled?
3. Should `line_origin` always be populated from founder line, even when no line-specific effects exist? My recommendation is yes.
4. Should `founder_allele_freq` move out of `genome_meta` into a line-aware long table?
5. Should hot haplotype tables use `locus_id` internally and expose `locus_name` through joins?
6. Is mutation support actually part of this refactor, or should the schema avoid blocking it while implementation waits?

## Bottom line

The core idea is worth doing, but the current plan is not yet safe to implement as written. The most likely ways to create a mess are:

- coding polyploid/sex chromosome behavior before `parent_origin` semantics are unambiguous;
- implementing hybrid sparse/dense rows without every query using the same missing-as-zero convention;
- losing line-specific founder pool information;
- rewriting TBV without preserving centering and current trait rules;
- renaming major tables without an explicit compatibility/migration strategy.

Narrow the first refactor to "long diploid haplotypes with current behavior preserved." Once that passes parity tests and benchmarks, add line-origin TBV, dosage caching, and non-diploid inheritance as separate, testable layers.

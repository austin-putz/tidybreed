# Post-Implementation Critique: `plans/position_coordinates.md`

Review date: 2026-07-06

Scope: current implementation review only. The earlier pre-implementation critique
has been removed. I did not change implementation code or `CLAUDE.md`.

## Overall Assessment

The core coordinate migration is mostly implemented correctly:

- `define_genome()` now writes `genome_meta.pos_bp` as an explicitly typed DuckDB
  `BIGINT`, with no `pos_Mb`, no `pos_cM`, and no `introduced_gen`.
- `define_genome()` writes a separate long `genome_map` table with default rows
  `(sex = NULL, line_name = NULL, map_name = "default")`, and derives `pos_cM`
  from sampled `pos_bp`.
- `cM_per_Mb` is implemented and validated as finite/positive.
- `chr_meta.recombines` was split into `recombines_M` and `recombines_F`.
- Founder-LD methods now read genetic distance from `resolve_genome_map()`.
- The `founder_allele_freq` write is now `ALTER TABLE` + `UPDATE`, avoiding a
  full `genome_meta` rewrite and preserving `pos_bp BIGINT`.
- Schema registries include `genome_map`, reserved columns, and system-table
  protection.

Focused verification run:

```sh
Rscript -e 'pkgload::load_all(".", quiet = TRUE); testthat::test_file("tests/testthat/test-genome_map.R"); testthat::test_file("tests/testthat/test-define_genome.R"); testthat::test_file("tests/testthat/test-define_chr.R"); testthat::test_file("tests/testthat/test-sex_chromosomes.R")'
```

Result: pass. I did not run the full package test suite or `R CMD check`.

## Findings

### 1. High: `add_offspring()` does not actually resolve maps per gamete

The plan says recombination readers should fetch parent `ind_meta.line_name` and
resolve the genetic map for the gamete-producing parent:

```text
resolve_genome_map(conn, sex = parent sex, line_name = parent line)
```

The implementation fetches `parent_line`, but then resolves the map once with no
sex or line:

- `R/add_offspring.R:205-209` fetches and stores `parent_line`.
- `R/add_offspring.R:220-225` calls `resolve_genome_map(pop$db_conn)` once.
- `R/add_offspring.R:413` and `R/add_offspring.R:417` pass the same
  `autosome_chr_info` to both sire and dam gametes.
- `R/add_offspring.R:456` uses the same default-map `special_chr_local_info` for
  recombining special chromosomes.

That means manually adding valid sex-specific or line-specific `genome_map` rows
will not change `add_offspring()` recombination. The resolver itself supports
sex fallback, but `add_offspring()` does not use it.

This matters because `plans/optimize_add_offspring.md` explicitly assumes the
caller resolves the map per gamete before passing Morgans-native positions to the
kernel. If this first-stage migration is supposed to provide that seam, it is not
complete yet.

Suggested fix direction: cache resolved maps/`chr_info` by unique
`(parent_sex, parent_line, map_name)` among parents, then choose the sire/dam
`chr_info` inside the offspring loop. For special chromosomes, use the
contributing parent's resolved local chromosome info when `k_p == 2L` and
`recombines_*` is true.

Add a regression test that inserts a monotonic male-specific or line-specific map
override, runs `add_offspring()`, and proves the override changes recombination
behavior or at least proves `build_chr_info()`/`make_gamete()` receives the
parent-specific resolved map. The current tests only exercise
`resolve_genome_map()` directly.

### 2. Medium: `CLAUDE.md` mostly matches the schema, but the `initialize_genome()` section is inaccurate

`CLAUDE.md` correctly documents the new `genome_meta`, `genome_map`, and
`chr_meta` shapes. The mismatch is the implemented-functions section for
`initialize_genome()`:

- `CLAUDE.md:551-553` lists `cM_per_Mb` as an `initialize_genome()` key param.
- `R/initialize_genome.R:32-39` has no `cM_per_Mb` argument.
- `R/initialize_genome.R:63-70` calls `define_genome()` without forwarding
  `cM_per_Mb`.

Since `initialize_genome()` is deprecated, the cleanest fix is probably to make
`CLAUDE.md` describe `define_genome()` as the current implementation surface and
keep `initialize_genome()` as a deprecated wrapper. If the wrapper is expected to
remain functional for this migration, it should accept and pass `cM_per_Mb`.

### 3. Medium: user-facing recombination docs still describe the old Mb model

The code now uses `pos_cM` in `make_gamete()`, but `add_offspring()`'s roxygen
still says:

```text
Poisson(chr_len_Mb[i] / 100), assuming approximately 1 Morgan per 100 Mb
```

See `R/add_offspring.R:40-44` and the generated `man/add_offspring.Rd`. This is
now wrong for non-default `cM_per_Mb` and contradicts the reason for this
migration. It should describe `Poisson(chr_len_cM / 100)`, where `chr_len_cM`
comes from the resolved genetic map.

Minor related stale wording:

- `R/recombination_helpers.R:87-90` still says `chr_meta.recombines = FALSE`
  instead of `recombines_M`/`recombines_F`.
- `R/add_offspring.R:389-393` still claims the autosome path is
  "BYTE-IDENTICAL to the pre-Stage-4 code", which conflicts with the
  pre-1.0.0 stance and is no longer useful design guidance.

### 4. Medium: `genome_map` validation is good, but not quite strict enough for the table contract

`validate_genome_map()` correctly catches duplicate nullable logical keys, orphan
`locus_id`s, locus-name mismatches, invalid sex, and `NULL map_name`.

Remaining gaps:

- It does not validate `pos_cM` is non-missing and finite for every stored row.
  `resolve_genome_map()` later treats missing resolved `pos_cM` as "no row", but
  a map write should fail earlier because `genome_map` is defined as rows with a
  genetic position.
- It does not validate `id_genome_map` uniqueness, even though the plan and
  `CLAUDE.md` call it a surrogate PK. Either enforce uniqueness in the validator
  or soften the docs to say the id is R-managed but not DB-constrained.
- The tests do not cover invalid sex, `NULL map_name`, missing resolved loci, or
  line-specific precedence. The code appears to handle invalid sex and
  `NULL map_name`, but the contract is important enough to test.

### 5. Medium: `genome_map` has no DB-level primary key/unique constraint despite "PK" wording

`R/define_genome.R:191-195` creates `genome_map` with typed columns but no
`PRIMARY KEY`, and `validate_genome_map()` does not check `id_genome_map`.
This may be consistent with the package's current R-side enforcement style, but
the plan/CLAUDE wording says "Surrogate PK".

This is not blocking the current default map, because internal writes assign
`id_genome_map` with `next_int_id()`. It is still worth resolving before a future
map writer exists, because `remove_rows()`/`mutate_table()` treat
`id_genome_map` as the row key.

### 6. Low: `define_genome()` validates `cM_per_Mb` carefully, but not `chr_len_Mb`

`R/define_genome.R:77-81` validates `cM_per_Mb` as finite and positive, but
`chr_len_Mb` only gets a type/length check. Negative, `NA`, `NaN`, or `Inf`
physical lengths produce downstream errors or confusing "too short" messages.

This is not a regression from the coordinate plan, but it is part of the same
user-facing argument surface and easy to harden while touching `define_genome()`:
require finite, non-missing, strictly positive `chr_len_Mb`.

### 7. Low: `define_genome()` can partially overwrite an existing genome if called twice

`R/define_genome.R:152-170` uses `CREATE OR REPLACE` for `genome_meta`, then
`R/define_genome.R:172-185` uses plain `CREATE TABLE` for `ind_haplotype` and
`ind_genotype`. If a user accidentally calls `define_genome()` on a population
that already has genome tables, `genome_meta`/`genome_map` can be replaced before
the later `CREATE TABLE` fails.

This is outside the plan's main coordinate goal, but it is a sharp edge around
the newly explicit table creation. A preflight "genome already defined" error
before any writes would avoid a partial replacement.

### 8. Low: crossover length semantics should be kept explicit for the next plan

`build_chr_info()` sets chromosome length to `max(pos_cM)` for the chromosome.
That matches the current observed-marker inheritance model: crossovers after the
last marker do not change stored haplotypes. It is not the same as full
chromosome genetic length unless the last locus is at the chromosome end.

This is fine for current haplotype sampling, but `plans/optimize_add_offspring.md`
adds optional crossover-event storage. If that feature is intended to store all
biological crossovers, the next plan needs an explicit chromosome-end map length
or terminal map point. If it only stores crossovers that affect represented loci,
the current `max(pos_cM)` convention is acceptable and should be stated.

## Test Gaps To Add Before Treating This Stage As Fully Closed

- `add_offspring()` uses a parent-specific resolved map when sex/line-specific
  `genome_map` rows exist.
- `resolve_genome_map()` line precedence:
  `(sex=S,line=L)` > `(sex=S,NULL)` > `(NULL,line=L)` > `(NULL,NULL)`.
- `resolve_genome_map()` errors when a resolved locus is missing.
- `validate_genome_map()` catches invalid sex and `NULL map_name`.
- `validate_genome_map()` catches missing/non-finite `pos_cM` if that stricter
  contract is adopted.
- `genome_map` is present in schema metadata tests, not just `pop$tables`.

## Bottom Line

The storage migration itself is in good shape: `pos_bp`/`genome_map` are real,
`BIGINT` is preserved, founder LD moved to cM, and old `genome_meta` columns are
gone from active code.

The important missing piece is the recombination reader seam: current
`add_offspring()` still uses one default resolved map for all gametes. That should
be fixed before starting the optimized Morgans-native kernel, otherwise the next
stage will be built on an API contract the current R implementation does not
actually exercise.

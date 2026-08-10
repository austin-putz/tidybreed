# define_genome() — Review Findings & Bug-Fix Plan

**Source:** Code review of [`R/define_genome.R`](/Users/austinputz/Claude/tidybreed/R/define_genome.R) (2026-08-09) plus `R/genome_map_helpers.R`, `R/chr_meta_helpers.R`, `R/open_pop.R`, `R/schema.R`, `R/sql_utils.R`, and tests in `tests/testthat/test-define_genome.R` / `test-genome_map.R`. No code was changed in this pass — findings only.

**Status:** Draft — awaiting review. Pre-1.0.0, so breaking fixes are preferred.

---

## Summary

`define_genome()` is correct on the happy path but has **2 functional bugs** (non-atomic multi-table write, zero-loci chromosome crash) and **5 high-priority validation gaps** that leave the DB in a partially-created state and allow silently corrupt inputs. The fix is small and contained to `R/define_genome.R`.

---

## 1. Critical — must fix

### 1A. Not atomic: partial writes persist on error

*Location:* [`R/define_genome.R:202-351`](/Users/austinputz/Claude/tidybreed/R/define_genome.R:202) — seven `CREATE`/`INSERT` statements plus three `validate_*` calls run with DuckDB autocommit; no `BEGIN`/`COMMIT`/`ROLLBACK`.

*Observed behavior:* `genome_meta` and `genome_map` use `CREATE OR REPLACE` (lines 202, 262) while `ind_haplotype`/`ind_genotype`/`ind_crossover`/`chr_inheritance`/`chr_recombination` (lines 230, 236, 251, 292, 301) use plain `CREATE`. If any `INSERT` or validator (`validate_genome_map` at 284, `validate_chr_*` at 350) fails, the first tables remain. The next call is blocked by the preflight at line 112 (`genome_meta COUNT>0`) and would also hit `table already exists` on the haplotype tables — no clean retry. Test at [`tests/testthat/test-genome_map.R:50`](/Users/austinputz/Claude/tidybreed/tests/testthat/test-genome_map.R:50) expects `genome_meta` absent after a `too short` error; it passes today only because the `pos_bp` guard (line 187) fires before the `CREATE`, but any later failure would violate the expectation.

*Would-be fix:* Wrap the entire body (after argument validation / preflight) in one transaction — `DBI::dbExecute(conn,"BEGIN TRANSACTION")`, `on.exit` rollback guard, `COMMIT` on success. Pattern already used by `define_chromosome()` at [`R/define_chromosome.R:222`](/Users/austinputz/Claude/tidybreed/R/define_chromosome.R:222). On rollback the DB is as if `define_genome()` was never called; on success nothing changes externally.

*Risk of not fixing:* Unrecoverable DB after any mid-flight error; requires manual `DROP TABLE` or file delete.

### 1B. Zero-loci chromosome: wrong `pos_bp` derivation

*Location:* [`R/define_genome.R:168-182`](/Users/austinputz/Claude/tidybreed/R/define_genome.R:168)

```r
loci_per_chr   <- diff(round(seq(0, n_loci, length.out = n_chr + 1)))
...
for (i in seq_len(n_chr)) {
  chr_loci <- which(chr_assignment == i)
  n_chr_loci <- length(chr_loci)
  chr_len_bp_i <- round(chr_len_Mb[i] * 1e6)
  pos_bp_i <- round(seq(0, chr_len_bp_i, length.out = n_chr_loci + 2)[2:(n_chr_loci + 1)])
```

When `n_loci < n_chr`, some `n_chr_loci == 0`. Then `length.out = 2` and `[2:1]` in R yields `c(2,1)` — a length-2 vector assigned to a 0-length slot, and the `any(pos_bp_i <=0)` guard (line 187) fires incorrectly.

*Would-be fix:* Early `if (n_chr_loci == 0L) next` (or validate `n_loci >= n_chr` with a clear error — decision needed). No behavior change when every chromosome has ≥1 locus.

---

## 2. High — validation gaps (silent corruption)

### 2A. Preflight checks only `genome_meta`

*Location:* [`R/define_genome.R:112-117`](/Users/austinputz/Claude/tidybreed/R/define_genome.R:112)

*Gap:* After a partial failure that leaves `ind_haplotype` etc created but `genome_meta` empty, `COUNT==0` passes but subsequent plain `CREATE` fails. *Would-be fix:* Check `DBI::dbListTables()` for any of the seven genome tables (`genome_meta`, `genome_map`, `ind_haplotype`, `ind_genotype`, `ind_crossover`, `chr_inheritance`, `chr_recombination`) before proceeding.

### 2B. `n_loci` / `n_chr` accept non-integers / NA / Inf

*Location:* [`R/define_genome.R:93-94`](/Users/austinputz/Claude/tidybreed/R/define_genome.R:93)

`stopifnot(is.numeric(n_loci), length(n_loci)==1, n_loci>0)` allows `100.5`, `NA`, `NaN`, `Inf`. `NA>0` is `NA`, `stopifnot` emits `missing value where TRUE/FALSE needed`. Downstream `diff(round(seq(...)))` and `locus_id = seq_len(n_loci)` then misbehave. *Would-be fix:* `is.finite`, `!is.na`, `== as.integer`, `>0L` with explicit message, matching the roxygen `Integer scalar` type.

### 2C. `locus_names` allows NA / empty strings

*Location:* [`R/define_genome.R:126-141`](/Users/austinputz/Claude/tidybreed/R/define_genome.R:126)

Only `length` and `anyDuplicated` are checked. `NA` or `""` locus names would corrupt denormalized joins to `genome_effects`/`ind_haplotype`/`genome_map`. `chr_names` has the correct guard at lines 157-165 (`anyNA`, `nchar==0`, duplicates). *Would-be fix:* Apply the same `anyNA`/`nchar==0` check to `locus_names`.

### 2D. `pop` not validated via `validate_tidybreed_pop()`

*Location:* [`R/define_genome.R:90`](/Users/austinputz/Claude/tidybreed/R/define_genome.R:90)

Uses `stopifnot(inherits(pop,"tidybreed_pop"))`; the codebase standard (e.g. `define_chromosome`, `restore_pop`) is `validate_tidybreed_pop(pop)` which also checks `dbIsValid`. Without it, a closed connection errors later as a cryptic DBI failure.

### 2E. `pop$tables` updated by blind `c()`

*Location:* [`R/define_genome.R:357`](/Users/austinputz/Claude/tidybreed/R/define_genome.R:357)

`pop$tables <- c(pop$tables, ...)` can duplicate on retry and drift from `DBI::dbListTables()`. `restore_pop()` at [`R/restore_pop.R:80`](/Users/austinputz/Claude/tidybreed/R/restore_pop.R:80) refreshes from the DB. *Would-be fix:* `union()` or `pop$tables <- DBI::dbListTables(db_conn)` after commit; add final `validate_tidybreed_pop(pop)` before return.

---

## 3. Medium — consistency nits (low risk)

* **Mixed DDL modes:** `CREATE OR REPLACE` vs plain `CREATE` — intentional `CREATE OR REPLACE` for idempotency on `genome_meta`/`genome_map` but plain `CREATE` for others. Under a transaction the distinction matters less; without it, retry leaves orphans. Worth making uniform (all `CREATE OR REPLACE` or all guarded by transaction).
* **Error messages via `stopifnot`:** `locus_names`/`chr_names` length check (lines 129, 147) uses `stopifnot` → unhelpful message. Elsewhere explicit `stop(..., call=FALSE)` is used.
* **Missing final validation:** No `validate_tidybreed_pop(pop)` / schema check before return.
* **`chr_len_Mb`/`cM_per_Mb` scalar detection:** After expansion at lines 120-123, `length(unique(chr_len_Mb))==1` at line 361 cannot distinguish scalar input from `c(100,100,100)` — cosmetic only.

---

## 4. What is NOT a bug

* `pos_bp` `BIGINT` via explicit typed `CREATE` + `duckdb_register`/`INSERT` (lines 202-220) — correct and RNG-neutral (avoids `dbWriteTable`'s random temp name). Validated by test at [`tests/testthat/test-define_genome.R:30`](/Users/austinputz/Claude/tidybreed/tests/testthat/test-define_genome.R:30).
* `pos_cM = pos_bp/1e6 * cM_per_Mb` derivation (line 194) and per-chromosome `cM_per_Mb` — correct.
* `genome_map` long format with `id_genome_map` via `next_int_id` + `validate_genome_map` at `R/genome_map_helpers.R:108` — correct.
* `chr_inheritance`/`chr_recombination` split (lines 292-337) — correct per CLAUDE.md (offspring sex vs producing-parent sex).
* `founder_allele_freq` deferred to `define_founder_haplotypes()` via `ALTER TABLE` — correct.

---

## 5. Proposed plan (no code changed in this review)

1. `R/define_genome.R:90` — replace `stopifnot(inherits(...))` with `validate_tidybreed_pop(pop)` + `stopifnot(DBI::dbIsValid(...))` guard at entry.
2. `R/define_genome.R:93-104` — tighten `n_loci`/`n_chr` to finite integer scalar `>0`; explicit errors.
3. `R/define_genome.R:126-141` — add `anyNA`/`nchar==0` check for `locus_names`; improve length-mismatch message.
4. `R/define_genome.R:112-117` — broaden preflight to all seven genome tables.
5. `R/define_genome.R:168-195` — handle `n_chr_loci==0` (skip) and/or validate `n_loci >= n_chr`; keep existing `pos_bp_i` positivity/collision guard.
6. `R/define_genome.R:202-351` — wrap table creation + inserts + validators in `BEGIN`/`COMMIT`/`ROLLBACK` (mirror `define_chromosome()`); unify `CREATE` modes.
7. `R/define_genome.R:357` — refresh `pop$tables` from DB (`union` or `dbListTables`) and add `validate_tidybreed_pop(pop)` before return.

Estimated scope: ~30 lines in one file; no schema migration; existing tests green. Atomicity fix should add a new test: `define_genome()` with `chr_len_Mb=0.001` leaves no genome tables behind (already partially asserted at `test-genome_map.R:58`).

---

## Appendix: review prompt

> review the `define_genome` function and let me know if you spot any errors, bugs, or critical updates it needs to make

Artifacts inspected: `R/define_genome.R`, `R/genome_map_helpers.R`, `R/chr_meta_helpers.R`, `R/define_chromosome.R`, `R/open_pop.R`, `R/schema.R`, `R/sql_utils.R`, `R/tidybreed_pop.R`, `tests/testthat/test-define_genome.R`, `tests/testthat/test-genome_map.R`.

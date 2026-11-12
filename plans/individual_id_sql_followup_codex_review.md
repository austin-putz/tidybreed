# Review of `individual_id_sql_followup.md`

Updated 2026-11-04 after reviewing Claude's revised plan. The initial review
was dated 2026-11-03; the working tree still reports version 0.76.14.

**Recommendation: proceed with the revised plan. No blocking findings remain
at the design stage.** The revision addresses all four original findings and
explicitly includes the adjacent deterministic-mean fix with its own tests and
NEWS entry. This approves the plan, not an implementation: `.with_id_view()`
and the query changes are still unimplemented.

## Disposition of the original findings

| Original finding | Status in the revised plan |
|---|---|
| Set semantics and empty selections | Resolved. The helper registers `unique(ids)`, accepts character vectors without `NA`, rejects `NULL`, and gives `character(0)` an empty `VARCHAR` column. Public guards stay in place; the internal empty-ploidy behavior change is explicit. |
| Capture tests need branch-specific results | Resolved. The result matrix covers overwrite dosage, selection from the cache being overwritten, chip intersections, parent evaluations, pedigree depth, fixed-effect reads, and genotype files. Capture must establish that each intended statement actually ran. |
| Registration lifetime and actual nesting | Resolved. The plan distinguishes each changing ID set, preserves the dosage snapshot, uses fresh names without overwrite, installs cleanup after successful registration, and directly tests nested success/error lifetimes and RNG neutrality. Eager callback execution and return values are documented. |
| Bound the no-ID-in-SQL claim | Resolved. Scope is manual interpolation of resolved individual-ID vectors in the six modules. Rendered user predicates are explicitly excluded, capture fixtures use non-ID predicates, and the hard-rule documentation will state the exception. |

The six-module inventory still matches the direct ID-list construction sites
in `R/`. Keep the founder/offspring writers and existing phenotype helpers
outside this change, as the revised plan does. No schema migration is needed.

## Decisions and pushback I agree with

- **This is rule enforcement, not a demonstrated injection bug in current
  public paths.** These resolved IDs come from package-written individuals;
  the current writers restrict line names. The unescaped sites still need
  removal to satisfy the stated transport rule, but quoting alone is not the
  reason to characterize today's paths as broken.
- **Fresh names protect a shared helper's lifetime contract.** The current
  BLUPF90 helper sequence is sequential. DuckDB R 1.5.5 refuses duplicate
  registration names by default; it does not silently replace an outer view.
  The revised explanation correctly avoids claiming an existing nesting bug.
- **User-filter translation needs its own design.** Narrowing the hard rule
  explicitly is appropriate here. It does not weaken the requirement to
  register the resolved vectors this patch owns.
- **Keep `sql_in_list()` for non-individual values.** Its remaining trait,
  locus, chromosome, component, index, and owner callers justify retaining it.
  Remove its individual-ID callers and examples as planned.
- **A solver and full coverage run are unnecessary for these helper tests.**
  Parent averaging and the file builders are directly exercisable without
  BLUPF90. The specified focused regressions followed by the full suite are
  appropriate for the eventual cross-module implementation.

## Reproducibility changes are now adequately scoped

The plan now makes output ordering deliberate for genotype extraction,
pedigree output, and genotype-file rows. Testing the expected order directly
is stronger than hoping a small unordered query changes order between runs.
The remaining matrix aggregates are placed by ID/locus keys, and
`build_data_file()` already orders its per-animal result through `merge()`.

I also agree with including the repeated-record mean as a separate fix within
this work. `list_avg(list(pheno_value ORDER BY id_phenotype))` gives each
animal/phenotype group a total accumulation order using the primary key. It
follows the existing ordered-list pattern in `ind_tgv_total` without imposing
the decimal conversion and finite accumulator range of `GEV_ACC_TYPE`.
The earlier review did not require that accumulator; fixed-order evaluation
is an appropriate solution.

I checked the proposed expression with the installed DuckDB R 1.5.5. With
100,000 animals, five records each, and seeded normal-valued phenotypes:

| Physical record layout | Plain `AVG`: differing means, 1 vs 4 threads | Plain `AVG`: differing means, 1 vs 8 threads |
|---|---:|---:|
| Records contiguous by animal | 2 | 2 |
| Same records interleaved by record position | 19,457 | 27,793 |

The ordered-list mean was `identical()` across both layouts and all three
thread counts. Separate probes confirmed exact preservation of a singleton
`0.30000000000000004`, ignoring a `NULL` in a mixed group, and an all-`NULL`
result matching `AVG()`.

These are aggregate probes, not an end-to-end solver test or an independent
reproduction of every count and file-rounding claim in Claude's plan.

## Small test clarifications to carry into implementation

These do not block the revised design:

1. **Treat the row count as a fixture starting point, not a failure
   threshold.** The plan's “a smaller table does not diverge” is too absolute.
   Values, physical layout, and parallel execution matter as well as size.
   In additional probes, R-generated `runif(..., 0, 100)` fixtures with
   500,000 and 1,000,000 rows showed no thread drift, while the normal-valued
   fixture above exposed it. Use seeded values and distribute repeated
   records across scan chunks so the regression actually exercises a changing
   reduction order. Retain the plan's explicit requirement to demonstrate
   failure against plain `AVG()` before relying on the test; merely reaching
   100,000 rows is insufficient. That mutation check concerns the current
   algorithm's determinism, not matching historical output.
2. **Pin the mean's edge semantics in the builder test.** Alongside the large
   determinism fixture, include singleton, mixed-`NULL`, and all-`NULL` groups.
   Check exact singleton preservation and the builder's existing missing-value
   conversion to zero. Compare returned `data` with `expect_identical()` and
   compare `data.txt` contents across thread counts, so the test covers the
   solver input as well as the in-memory result. Do not compare to a stored
   pre-change output fixture.

The API-skill update should describe both the shared helper contract and the
builder's deliberately fixed-order mean. The existing schema skill needs no
change because the tables, keys, and schema invariants stay the same.

## Validation and limits of this update

- Compared the revised plan and its working-tree diff with the initial review
  and `CLAUDE.md`; consulted the project API and relevant schema references.
- Rechecked the direct ID-list inventory, dosage snapshot behavior, subset
  resolution, parent-average queries, genotype row construction, BLUPF90 file
  builders, and the existing determinism-test pattern against current code.
- Ran the DuckDB aggregate and mean edge-case probes described above.

Only this review document was changed. The main plan and production code were
not edited. The follow-up helper does not yet exist, so its cleanup, nesting,
SQL capture, and branch results remain implementation acceptance criteria.
No solver was launched and the package suite was not run for this
documentation-only update.

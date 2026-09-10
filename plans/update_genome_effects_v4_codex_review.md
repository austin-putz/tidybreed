# Final review of `update_genome_effects_v4.md` after v4.6

**Review date:** 2026-09-06  
**Reviewed proposal:** `plans/update_genome_effects_v4.md`, revision v4.6  
**Status:** the Q1 decision is approved, but the replacement contract needs two
correctness edits before implementation.

## Executive verdict

Deleting `trait_meta.expressed_parent` is the right decision. Imprinting belongs on
the affected additive effect, not on the whole trait, and retaining both mechanisms
would create two sources of truth. There is no migration or compatibility requirement
before 1.0.0; the old column and argument should simply disappear.

Adding `line_match_type = 'any'` is also necessary. Without a stored `ANY` line value,
the schema cannot represent `(ANY line, parent 1)` because the only previous spelling
of `ANY` was zero origin rows, which has nowhere to store `parent_origin`.

The containment lattice remains coherent after that addition. In particular:

- common `(ANY, ANY)` contains parent-only `(ANY, 1)`;
- parent-only `(ANY, 1)` contains line-and-parent `(exact A, 1)`;
- line-only `(exact A, ANY)` and parent-only `(ANY, 1)` overlap without containment
  and are correctly rejected inside one fallback family;
- reciprocal genotype predicates remain exact multisets and do not need `'any'`.

The v4.6 edit did, however, expose two contracts that are not yet safe to implement:

1. `define_additive_effects()` must compose `line_name` and `parent_origin` in one
   origin row; it must not stamp an `'any'` row unconditionally.
2. `scale_to_target = TRUE` must use the number of eligible parental copies. The
   current `2pq` variance factor is wrong when only one parent is expressed.

The shape of `parent_origin` for multi-trait calls also needs to be declared. These
corrections require no new table or column.

## What v4.6 gets right

- `expressed_parent` is deleted rather than deprecated or aliased.
- Parent-of-origin becomes effect-local and can vary by locus and effect set.
- Parent-only additive matching has a normalized stored representation:
  `('any', NULL, parent_origin, copy_count = 1)`.
- `'any'` is restricted to additive members. Genotype contrasts still require an
  exact origin multiset, so dominance is never computed from one allele copy.
- `('any', NULL)` is rejected as a redundant spelling of the common scope.
- The old trait-wide numerical behavior can be represented without changing the
  three-table design.
- Phase C removes the legacy flag before Phase D rewires `add_tbv()`.

## Blocking issue 1: the wrapper currently erases line specificity

The Writer API section says that `define_additive_effects(parent_origin = ...)`
“stamps an `('any', parent_origin, copy_count = 1)` origin row on every term.” That is
correct only when `line_name = NULL`.

When `line_name = "Duroc"` and `parent_origin = 1`, the one permitted additive origin
row must be:

```text
(line_match_type = 'exact', line_name = 'Duroc', parent_origin = 1,
 copy_count = 1)
```

Writing an `'any'` row would lose the Duroc condition. A second call for Landrace
would then produce the same family and scope, either failing as a duplicate or
replacing the Duroc variant. This breaks the existing line-specific imprinted case in
`tests/testthat/test-add_tbv.R` and defeats per-line centering.

The wrapper needs this explicit mapping:

| `line_name` | `parent_origin` | Stored additive scope |
|---|---:|---|
| `NULL` | `NULL` | zero origin rows (common) |
| `A` | `NULL` | `(exact A, ANY parent)` |
| `NULL` | `1` or `2` | `(ANY line, exact parent)` using `'any'` |
| `A` | `1` or `2` | `(exact A, exact parent)` |

This is one row in every scoped case, consistent with the additive-member invariant.
Update the Writer API, the imprinting decision, and gate 40 to use this mapping.

Gate 40 should not test only a population-wide `('any', 1)` term. It should preserve
the stronger existing fixture: Duroc and Landrace line-specific variants, different
coefficients and centering constants, an F1, and expression from one parent. That
fixture detects both a lost line dimension and an incorrect parent filter.

## Blocking issue 2: variance targeting must become origin-aware

`define_additive_effects(scale_to_target = TRUE)` currently rescales with:

```text
V_A = sum(2 p_j q_j a_j^2)
```

The factor 2 is the number of independent autosomal copies contributing to an
ordinary additive value. A parent-qualified additive term contributes one copy, so
under the same HWE assumptions its variance is:

```text
V_A,parent = sum(p_j q_j a_j^2)
```

If v4.6 preserves `2pq` after moving imprinting into the effect rows, an imprinted
model requested at `target_add_var = V` will have expected variance `V / 2`. This was
latent in the old trait-wide implementation, but the new API now owns the scope and
must not carry that bug into the rewritten writer.

For the reserved order-one additive set, state the scaling contract as:

```text
V_A = sum(n_eligible,j p_j q_j a_j^2)
```

where `n_eligible,j = 2` for an unparented autosomal scope and `1` for an exact
`parent_origin`. The line-specific base frequency remains the `p_j` used by that
variant. If a more general origin mixture cannot be reduced safely by the helper, the
helper should reject `scale_to_target = TRUE` rather than silently use `2pq`.

Add an empirical or enumerated acceptance gate showing that ordinary and
parent-qualified generated effects each hit the requested target variance under the
plan's base-population assumptions.

## Blocking issue 3: multi-trait `parent_origin` is unspecified

`define_additive_effects()` accepts a vector of trait names, but v4.6 merely says it
“gains `parent_origin = NULL`.” The plan must state whether this argument is:

- one scalar applied to every trait; or
- a named/per-trait vector, with `NULL`/`1`/`2` resolved separately for each trait.

The second form preserves the expressiveness of the removed per-trait metadata. If
the first form is intentional, mixed parent expression in one correlated multi-trait
call must fail clearly rather than recycle ambiguously.

There is also a covariance consequence. Under random mating, a trait using only the
paternal homolog and another using only the maternal homolog do not obtain their
requested off-diagonal genetic covariance merely because their sampled coefficients
are correlated. The plan should either implement an origin-aware covariance
calculation or reject incompatible mixed-origin `G`/`scale_to_target` calls. A scalar
parent applied to every trait avoids this particular case but should be documented as
a deliberate limitation.

## Important issue 4: put the new row-local invariant in the DDL

The plan says every row-local invariant is a SQL constraint, but the displayed DDL
allows:

```text
line_match_type = 'any', line_name = NULL, parent_origin = NULL
```

and relies on the R validator to reject it. Because “`'any'` requires a non-NULL
parent” is row-local, incorporate it into the `CHECK`, for example by giving
`'exact'`, `'unknown'`, and `'any'` separate branches. The rule that `'any'` is valid
only for an additive member remains cross-table validation in R.

Gate 41 should test both the public writer and a direct invalid SQL insert, matching
the plan's declared division of responsibility.

## Important issue 5: decide the no-eligible-copy behavior

The current `add_tbv()` errors when an individual has no haplotype row matching any
additive effect, including a parent-qualified effect on a chromosome not inherited
from that parent. The new general evaluator says “no match contributes 0.” Those are
different contracts.

For example, a paternally qualified X-linked effect in a male has no eligible copy.
Returning zero is defensible and arguably more natural for an effect-local predicate,
but it is not numerical-and-diagnostic parity with the current engine. The plan must
choose explicitly between:

- zero contribution for a valid but absent eligible copy; or
- retaining an error when every reserved additive term is ineligible for an
  evaluated individual.

Whichever behavior is chosen needs a sex-chromosome gate. Gate 40's autosomal parity
case cannot detect it.

## Important issue 6: the deletion inventory is incomplete

The Phase-C table covers the functional R code, registries, generated man pages, and
three directly affected tests. Repository-wide search finds additional live
documentation/examples that would be left stale:

- `CLAUDE.md` has four current-contract references, including the detailed
  `define_trait()` and `add_tbv()` sections, not only the two sites named in the plan.
- `vignettes/tidybreed-introduction.Rmd` lists `expressed_parent` as a
  `define_trait()` argument.
- `vignettes/swine/swine-time-based-age-at-puberty-sex-semen.R` passes
  `expressed_parent = "both"` in two executable examples; those calls will error once
  the argument is removed.

Historical plans and the historical NEWS entry may remain unchanged because they
describe earlier revisions. The live vignette sources and current package contract
must be included in Phase C.

## Pre-existing DDL issue found during the re-audit

The `genome_effect_members` `CHECK` does not actually require `center_value` for
`additive` or `dominance`. In SQL, `CHECK` accepts `UNKNOWN`, so
`center_value BETWEEN 0 AND 1` permits NULL unless the branch also says
`center_value IS NOT NULL`. DuckDB 1.5.5 was verified to accept a NULL under a bare
`CHECK(center BETWEEN 0 AND 1)`.

Add `center_value IS NOT NULL` to the non-indicator branch and a direct constraint
test. This was not caused by Q1, but it contradicts the plan's implementable-readiness
claim and should be fixed in the same final semantic edit.

## Minor corrections

1. The title says v4.6 but the metadata line still says `Revision: v4.5`; change it
   to v4.6.
2. The acceptance subsection titled “Origin containment and states (26–39)” now
   contains gates 26–41.
3. The blast-radius count should include the live vignette sources identified above.

## Acceptance-gate revisions

Keep gates 1–39. Revise and add the Q1 gates as follows:

40. **Full imprinting parity.** Preserve the current line-specific F1 fixture:
    Duroc and Landrace variants with distinct values and centers, exact line plus
    parent in the same origin row, and equality to an independent per-copy oracle.
41. **`'any'` constraints.** Reject `('any', NULL parent)` through both the writer and
    direct SQL; reject `'any'` on a genotype member through the writer.
42. **Wrapper scope matrix.** Round-trip all four combinations of `line_name` and
    `parent_origin` in the mapping above, including `replace_scope` isolation.
43. **Origin-aware variance target.** Ordinary two-copy and parent-qualified one-copy
    generated models each hit `target_add_var` under the stated assumptions.
44. **Multi-trait origin contract.** Validate the chosen scalar/per-trait argument
    shape and reject any unsupported mixed-origin covariance request.
45. **No eligible copy.** Lock the chosen zero-or-error behavior on a sex-chromosome
    fixture.
46. **Required center.** Direct insertion of an additive or dominance member with
    NULL `center_value` fails.

## Final decision summary

| Decision | Verdict |
|---|---|
| Delete `trait_meta.expressed_parent` | **Approved** |
| Migration or compatibility shim | **None; pre-1.0 removal is intentional** |
| Add `line_match_type = 'any'` | **Approved and necessary** |
| Restrict `'any'` to additive members with a parent | **Approved** |
| Predicate-containment lattice | **Still approved** |
| Wrapper's unconditional `'any'` stamping | **Incorrect; compose line and parent in one row** |
| Imprinted variance scaling | **Must use one eligible copy, not `2pq`** |
| Multi-trait `parent_origin` | **Specify before implementation** |
| New schema objects | **None required** |

## Final recommendation

Keep the v4.6 Q1 decision and the new `'any'` lattice value. Before Phase A is treated
as closed, revise the wrapper mapping, make generated-effect scaling origin-aware,
declare the multi-trait argument contract, and add gates 40–46 above. Then the plan is
implementable without reopening its architecture or adding migration work.

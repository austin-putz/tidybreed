# Codex critique: `update_chr_meta.md` v5

**Reviewer**: Codex  
**Status**: Revised review after the two-table redesign  
**Recommendation**: **Approve the architecture after one remaining semantic
decision and several plan clarifications.** The v5 split is a substantial
improvement. It resolves the central defects in the earlier proposal rather than
documenting around them.

---

## Executive summary

The new design is the right shape:

```text
chr_inheritance(chr_name, offspring_sex, line_name,
                from_parent_1, from_parent_2)

chr_recombination(chr_name, parent_sex, line_name, recombines)
```

Separating inheritance from recombination gives each sex dimension one clear
meaning, prevents copy rules from shadowing recombination rules, and makes parent
contribution explicit. Removing `copy_mode` and `hemi_parent` also eliminates the
old validation contradiction for absent Y/W chromosomes.

The original blocking issues now dispose as follows:

| Earlier issue | v5 disposition |
|---|---|
| `copy_mode="none"` contradicted `hemi_parent` validation | **Resolved structurally**: both concepts are gone; absence is `0,0`. |
| Parent contribution was implicit and too coarse | **Resolved for explicit copy counts** by `from_parent_1` / `from_parent_2`. |
| One `sex` meant offspring sex and producing-parent sex | **Resolved for those two contexts** by separate tables and names. |
| `line_name` was claimed to model crossbred structural variants | **Resolved by narrowing scope** to pure-line/founder defaults and documenting per-copy SVs as future work. |
| Complete-row fallback caused cross-concern shadowing | **Resolved structurally** by independent tables/resolvers. |
| “Never update later” conflicted with known future columns | **Resolved by reframing** the goal as avoiding another fundamental rewrite. |

I no longer object to implementing two long tables without the speculative
polyploid, compartment, and structural-state columns. Under this design those are
ordinary additive extensions rather than another normalization rewrite.

---

## Remaining decision before implementation: parent sex is not gamete role

The split resolves **offspring sex vs producing-parent sex**, but it does not
resolve the third context from the earlier review: **gamete role**.

`chr_recombination.parent_sex` works for the package's current dioecious M/F
model, where a male parent occupies `parent_1` and a female parent occupies
`parent_2`. It does not express different recombination behavior for pollen vs
ovule production by the same monoecious or hermaphroditic individual. In that
case the individual's `ind_meta.sex` cannot select both maps.

This conflicts with two v5 claims:

- “Hermaphrodite / monoecious ... ✓ trivial.”
- “the two explicit tables express the storage for every system above.”

Choose and state one of these scopes:

1. **Recommended narrow scope for this refactor:** keep `parent_sex`, explicitly
   limit recombination resolution to the current `{M,F}` dioecious parent model,
   and move sex-role-decoupled selfing/monoecy to Boundaries. Autosome inheritance
   remains representable; role-specific recombination is not yet supported.
2. Key recombination by `parent_role` (`parent_1` / `parent_2`) now. This naturally
   matches mating position and selfing, but changes the meaning of the existing
   M/F API and map lookup.
3. Add both producer sex and gamete role as independent dimensions. This is most
   general but is unnecessary if the current package cannot create such crosses.

This is not an argument to reunify the tables. It is a request to make the
remaining dimension choice honest. I would approve option 1 for the current
diploid refactor.

---

## Major clarification: copy counts do not define the transmission mechanism

The storage can represent `2,0`, but the implementation plan does not yet explain
how `add_offspring()` obtains two child haplotypes from one diploid parent's
gamete. The file-by-file section says only:

> parent-origin slot 1 gets `from_parent_1` strands

and:

> `from_parent_1` from the sire's gamete

A normal diploid gamete contains one homolog. Producing two child copies from one
parent requires a rule: duplicate one transmitted haplotype, transmit two
distinct parental homologs, bypass meiosis, or invoke another mechanism. Those
choices have different genetic consequences.

Therefore:

- `2,0` is **storage-expressible**, but not biologically implemented merely by
  changing row counts.
- The plan should not present diploid uniparental inheritance as supported until
  the strand-selection behavior is specified and tested.
- For this refactor, validate only the count patterns the kernel implements, or
  allow `2,0` in storage but fail clearly at the kernel boundary.

The same distinction should be applied to the polyploid coverage language:
counts describe desired output cardinality; they do not by themselves implement
pairing, gamete formation, or selection of homologs.

---

## Major clarification: the ploidy validation rule has no single database value

The proposed validator says:

> `from_parent_1 + from_parent_2 ≤ ind_meta.ploidy`

But `ind_meta.ploidy` is per individual, may contain multiple values, and may be
empty when chromosome rules are defined. A chromosome configuration row cannot
be validated against “the” value of that column without specifying a population
or context.

For this explicitly diploid release, use a schema/version capability check such
as:

```text
from_parent_1 + from_parent_2 <= 2
```

and retain the existing `assert_ploidy_2()` boundary for consumers. When
non-diploid support lands, validation must occur against the particular offspring
and parents being processed, or against a separate genome-level ploidy contract;
it cannot be a timeless invariant of these config rows.

Also add a kernel-boundary assertion that the resolved total equals the number of
haplotype rows actually created. A `<= 2` check alone permits legitimate absence
but does not catch an implementation that silently emits the wrong number of
strands.

---

## Resolver and line semantics need two precise statements

### 1. Define which line selects each table

The plan says both maps are keyed by `(sex, line_name)`, but not always whose
`line_name` is passed:

- inheritance should be explicit about whether it uses the **offspring line**;
- recombination should be explicit about whether it uses the **producing
  parent's line**.

Those are normally different in a cross. This should be pinned in
`test-resolver-context.R` alongside the sex-context tests. Passing the offspring
line to recombination would silently defeat founder/pure-line recombination
defaults in crossbred offspring.

### 2. Specify the within-table shadowing validator algorithm

“A line-specific row whose value would be out-prioritized” is not yet an
implementable invariant. State exactly which contexts are enumerated and what
constitutes an error. A defensible rule is:

1. enumerate every configured non-NULL line and both supported sexes;
2. if both `(sex,NULL)` and `(NULL,line)` match and differ, require an explicit
   `(sex,line)` row;
3. apply this independently to each table;
4. include chromosome, sex, and line in the error.

This makes the intended anti-shadowing behavior deterministic and testable.

---

## API and transaction details to finish

### `overwrite` is still unspecified

The plan now clearly says every call is a one-row upsert, which fixes the earlier
ambiguity about clearing sibling rows. It still retains `overwrite = TRUE`
without defining `FALSE`.

State that `overwrite = FALSE` errors if the exact NULL-safe logical key already
exists and otherwise inserts. Test all four combinations: NULL/non-NULL sex and
NULL/non-NULL line. If that is not useful behavior, remove the argument.

### Validation must be inside the write transaction

The plan calls the upsert transaction-backed and then says to run the matching
validator. Make explicit that delete, insert, and validation occur in the **same
transaction**, with rollback on validation failure. Otherwise a bad replacement
can remain stored after the validator errors.

### Add database constraints where cheap

R validation is necessary for fallback resolution, but the tables should still
use straightforward database constraints for row-local invariants:

- non-null `chr_name`, counts, and `recombines`;
- counts restricted to the supported diploid range;
- sex restricted to `NULL`, `M`, or `F`;
- foreign key behavior consistent with the rest of the package.

DuckDB NULL semantics mean a normal composite `UNIQUE` constraint may not enforce
the desired NULL-normalized logical key, so the explicit validator remains
necessary unless an expression index/normalized stored key is chosen.

---

## Coverage wording to tighten

Most of v5's boundary language is now appropriately careful. I recommend four
small edits:

1. Change “express the storage for every system above” to “express the current
   diploid copy-count and whole-chromosome recombination configuration for the
   supported subset described above.”
2. Mark monoecious/hermaphrodite **role-specific recombination** unsupported if
   retaining `parent_sex`.
3. Mark `2,0` and autopolyploid counts as schema-expressible but kernel-future.
4. Keep the current qualifications for organelle proportions, PARs,
   B chromosomes, homeologous exchange, and crossbred SVs; those corrections are
   good and should remain.

The allopolyploid wording is much better than in v4. Calling standard disomic
allopolyploids an approximation is accurate enough for this plan.

---

## Tests I would add or sharpen

The proposed test list is strong. Add these cases:

1. Recombination uses the **producing parent's** sex and line, while inheritance
   uses the **offspring's** sex and line, in an actual cross where all differ.
2. Conflicting `(sex,NULL)` and `(NULL,line)` rows require `(sex,line)`; matching
   values do not spuriously fail.
3. Failed validation rolls back delete-then-insert and preserves the old row.
4. `overwrite = FALSE` respects NULL-safe composite keys.
5. Resolved copy totals exactly match emitted `ind_haplotype` row counts for
   autosome, X/Y or Z/W, absent chromosome, and organelle cases.
6. A currently unsupported `2,0` mechanism fails explicitly rather than
   duplicating or dropping strands accidentally, unless this refactor implements
   and tests its exact semantics.
7. If `parent_sex` is retained, document/test rejection of unsupported
   gamete-role-only recombination instead of implying selfing makes it work.

---

## Non-blocking observations

- No surrogate primary key is reasonable for these small managed configuration
  tables. The NULL-safe logical-key tooling and tests matter more.
- Four resolver queries in the common batch case are reasonable. The plan should
  cache by the actual distinct contexts needed by that batch, including parent
  line for recombination.
- Keeping resolvers and validators in `chr_meta_helpers.R` is coherent. Rename
  that helper during the later schema-wide naming sweep.
- The separate `chr` → `chromosome` follow-up commit is acceptable now that Q7
  is a decided action rather than an indefinite possibility.
- One concern per `define_chromosome()` call is slightly verbose but protects the
  key conceptual distinction and is worth it.

---

## Bottom line

The two-table redesign is the architectural correction I wanted from the first
review. It makes the common XY/ZW/X0, whole-chromosome achiasmy, and simple
uniparental-organelle cases substantially clearer, and it removes the most likely
sources of implementation error.

Before implementation, settle or narrow the `parent_sex` versus gamete-role
scope. Then clarify ploidy validation, `2,0` kernel semantics, line context,
shadowing validation, `overwrite = FALSE`, and transactional rollback. With
those changes, I recommend approval; speculative future biology columns do not
need to be pre-created.

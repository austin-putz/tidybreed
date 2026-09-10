# Review of `update_genome_effects_v3_parsimony.md`

**Review date:** 2026-09-05  
**Reviewed proposal:** `plans/update_genome_effects_v3_parsimony.md`  
**Review scope:** whether v3 provides a parsimonious, understandable foundation for
additive, dominance, and epistatic genome effects across lines, including effects
supplied by users.

## Executive verdict

v3 is a substantial improvement over the 16-table v2 design. Its central
representation—one coefficient term in `genome_effects` with one or more locus
members in `genome_effect_members`—is the right foundation. Indicator terms are a
complete representation for arbitrary functions of finite, unphased genotype-state
sets, so separate contrast-registry and genotype-surface table families are not
needed for the package's present biallelic scope.

However, v3 is **not implementation-ready as written**. It currently overclaims
three outcomes:

1. Non-additive terms are representable in storage but are explicitly not evaluated.
2. A scalar term-level `line_name` cannot express the line-origin semantics of a
   dominance or multi-locus coefficient.
3. No public general writer is planned, so the claimed user-injection path would
   require manually maintaining normalized IDs and invariants.

The recommendation is **not to return to v2**. Keep the two-table term/member core,
correct the contrast and model semantics, add a general writer and evaluator, and
use at most one additional origin-composition table if line- or cross-specific
non-additive effects are required now. Three coherent tables are preferable to two
tables with an overloaded `line_name` whose meaning cannot be evaluated.

## What v3 gets right

### 1. Term plus members is the correct normalization boundary

Epistasis forces one coefficient to refer to multiple loci. Separating a term's
coefficient from its locus members solves exactly that problem without constructing
a general state-registry framework.

The following representation is fundamentally sound:

```text
term contribution = coefficient x product(member contrast values)
model value        = sum(term contributions)
```

It naturally supports order-one effects, pairwise interactions, and interactions of
arbitrary order without another schema change.

### 2. Indicator completeness removes most of v2 legitimately

For a finite product of unphased genotype-state sets, indicator cells can represent
any lookup function exactly. Consequently, separate surface, cell, state, and
contrast-value table families are not necessary merely to let a user enter a 3 x 3
genotypic table.

The scope of that statement must remain precise. It covers functions of the modeled
per-locus genotype states. It does not cover phased haplotype combinations,
multiallelic identity, origin composition, sex, environment, or other dimensions
unless those dimensions are represented separately.

### 3. `base_allele_freq` belongs with the contrast member

A frequency-dependent coding is defined at a locus and can vary by line or base
context. Keeping the copied frequency on the member makes the numerical definition
restart-safe and fixes the current permissive `COALESCE(base_allele_freq, 0)` path in
`R/add_tbv.R:210`.

### 4. Member count should be derived

Not storing `effect_order` is consistent with the package's aversion to redundant
metadata. A view can expose member count and causal-locus membership.

### 5. Canonical member ordering is necessary

Canonicalizing by `locus_id` is the right way to make term matching independent of
input order. The tests should compare canonical query results rather than claiming
physical or byte-level database identity, because surrogate IDs and physical row
order are not semantic contracts.

## Blocking findings

### 1. The implementation plan does not deliver non-additive genome effects

The acceptance fixtures show that dominance and epistasis can be *stored*, but
Phases C and D write and evaluate only the existing order-one additive path. The
proposal then explicitly places non-additive evaluation, `ind_genetic_value`, and
phenotype integration out of scope (`plans/update_genome_effects_v3_parsimony.md:405-406`).

That means a user could successfully insert a dominance or epistatic term that has
no effect on any simulated genetic value or phenotype. Storage success is not an
adequate acceptance gate for an effect model.

The plan must choose one honest boundary:

- **Schema-foundation release:** state explicitly that only common and line-specific
  order-one additive effects are executable, and do not advertise day-one custom
  gene action.
- **General-effects release:** include a public general writer, a term evaluator,
  an output contract for total genetic values versus breeding values, and phenotype
  integration.

The stated goal favors the second option.

Non-additive totals should not silently be written to `ind_tbv`. A TBV is not the
same thing as an individual's total genotypic value once dominance and epistasis are
present. The implementation needs either `ind_genetic_value` or an equally explicit
result contract.

### 2. A scalar term-level `line_name` cannot mean allele-copy origin for general terms

The current additive evaluator joins each `ind_haplotype` row to an effect whose
`line_name` matches that copy's `line_origin`, with common fallback separately for
each locus and copy (`R/add_tbv.R:210-222`). That works because the current
coefficient multiplies one allele copy at one locus.

It does not extend unchanged to either of the new cases:

- A dominance or indicator contrast evaluates a collapsed local genotype. An A/B
  genotype contains two origins and has no single `line_name`.
- A multi-locus coefficient may combine genotype states with different origin
  compositions at different loci, especially after recombination. It likewise has
  no single allele-copy origin.

Therefore the claims at `plans/update_genome_effects_v3_parsimony.md:234-237` that
the scalar column expresses common, line-A, and line-B effects, and at `:215` that
the existing fallback is unchanged, are true only for order-one copy-additive
terms.

There are two valid parsimonious choices:

#### Restricted two-table release

Keep `line_name`, but enforce that non-NULL values are allowed only on an order-one
copy-additive term. Dominance, indicator, and epistatic terms are common-only until
origin composition is implemented.

#### Three-table general cross-line release

Add one origin-composition child table and design its matching and precedence
resolver in the same plan. The resolver must define common, within-line, A/B,
reciprocal A/B versus B/A, unknown origin, haploid, diploid, and supported polyploid
behavior.

The proposed assertion that this future table changes "zero existing readers" is
incorrect. Existing rows need not change, but every general evaluator must be
extended to resolve the new scope rows. In addition, the interaction between the
old scalar `line_name` and the new origin rows must be defined or the scalar column
must be removed.

Because the requested outcome includes effects "across lines," the three-table
option is the more faithful design unless that phrase is intentionally narrowed to
copy-additive effects only.

### 3. The contrast definitions are not yet executable specifications

#### Hard-coded diploidy

The proposed genotype-level additive coding is `dosage - 2p`. Tidybreed already
supports chromosomes with one or zero inherited copies, and `ind_haplotype` is the
authoritative copy-level representation (`R/define_genome.R:292-303`). For a
hemizygous locus the centered value is `dosage - p`, not `dosage - 2p`.

The general definition must use the individual's realized copy count:

```text
additive contrast = dosage - copy_count x base_allele_freq
```

Equivalently, it is the sum of `(allele - p)` over eligible allele copies.

#### Undefined aggregation for epistatic copy-additive terms

For an additive-by-additive term, the intended value should be stated explicitly as:

```text
[sum copies at locus 1 (allele - p1)] x
[sum copies at locus 2 (allele - p2)] x coefficient
```

The evaluator must aggregate each member to one individual-locus contrast value
before multiplying members. Multiplying raw haplotype rows directly would produce
ambiguous pairing and duplication behavior.

#### Underspecified dominance coding

"HWE Cockerham `x_D`" is not a sufficient stored contract. The exact value at each
diploid dosage state and the allele-frequency orientation must be documented and
tested. The evaluator must also define what happens at haploid loci and at copy
counts greater than two. A reasonable initial rule is:

- named `dominance` is valid only for realized diploid genotypes;
- arbitrary haploid or polyploid gene action uses indicators until a named coding is
  deliberately defined.

The proposed functional-to-Cockerham equivalence test also needs the population-mean
intercept. With allele-1 frequency `p`, allele substitution effect
`alpha = a + d(q - p)` and dominance coefficient `d` reproduce the *centered*
genotypic values. Reproducing raw `(-a, d, a)` values requires adding the appropriate
mean/intercept term.

#### Redundant additive contrast names

The current distinction between `'allele'` and `'additive'` is an implementation-unit
distinction, not necessarily two user concepts. A smaller closed vocabulary can use
one ploidy-aware `additive` contrast defined as a sum across eligible copies, plus
`dominance` and `indicator`. If both names remain, their different behavior must be
observable, necessary, and fully specified.

### 4. `model_name` is stored but has no selection semantics

The schema permits multiple models for the same trait. The planned replacement key
also treats each model independently. But the current evaluator selects by trait and
effect type, while `ind_tbv` has uniqueness only on `(id_ind, trait_name)`
(`R/define_trait.R:221-228`). No phase explains which model is active, whether models
are summed, or how outputs from multiple models coexist.

Choose one of the following:

- Wire `model_name` through writer arguments, evaluator selection, result identity,
  phenotype components, and replacement rules.
- Remove `model_name` for now and add it when multiple executable models are needed.

The second option follows v3's own rule that an easily added dimension should be
deferred. The first is justified only if alternate named models are an actual
day-one feature.

### 5. The plan lacks the user-facing injection path required by its goal

Manually inserting into two normalized tables requires users to coordinate:

- surrogate effect IDs;
- canonical member slots;
- term and member rows;
- logical uniqueness across a parent and its child rows;
- conditional `genotype` and frequency requirements;
- model replacement and transaction boundaries.

That is not a simple custom-effect interface, even though it is a small database
schema.

Add an exported `define_genome_effects()` writer. It should accept named loci and a
tidy or list-based term specification, resolve locus IDs, canonicalize members,
validate the entire model, and replace or append in one transaction. Convenience
helpers may convert common inputs such as a 3 x 3 genotype surface into indicator
terms, but those helpers do not require additional database tables.

Direct SQL writes should not be the documented extension point. If direct writes
remain supported, the DDL needs substantially more `CHECK`, foreign-key, and
logical-uniqueness protection than the sketch currently contains.

## Important non-blocking corrections

### 1. `effect_class` is underspecified and is not always biologically derivable

An indicator for one genotype state is a basis cell, not intrinsically an additive
or dominance effect. Its decomposition depends on the rest of the model and the
chosen basis. Consequently, "validated against the ordered member contrasts" is not
enough to define `effect_class`.

Derive unambiguous properties such as member count and contrast signature in views.
If a user-facing label is useful, store an optional `effect_name` or `term_name` and
do not treat it as mathematical identity. If `effect_class` remains, define a small
closed vocabulary including a neutral custom/cell class and specify its validation
rules.

### 2. `locus_name` need not be duplicated in `genome_effect_members`

`locus_id` is the stable join key and `ind_haplotype` already carries it. The public
writer can accept `locus_name`, resolve it once, and the display view can join
`genome_meta` to expose it. Removing the duplicated name makes the schema smaller
and eliminates a locus ID/name agreement invariant.

### 3. `genotype` should describe what it stores

For the current biallelic dosage-state model, `dosage_value` or `genotype_state` is
clearer than `genotype`. If the column is intended to hold copy count without an
artificial ceiling, use `INTEGER`; otherwise acknowledge the `UTINYINT` limit rather
than claiming unrestricted ploidy storage.

### 4. Indicator omission is zero only under explicit component semantics

Sparse indicator storage can unambiguously mean that an omitted cell contributes
zero. That is a reasonable internal rule. The writer still needs to distinguish the
common user intentions:

- these cells are additional model components;
- these values define the complete genotypic surface and replace overlapping terms.

This distinction may be a writer operation rather than a persisted
`value_semantics` column, but it cannot simply disappear from the user API. Optional
surface-completeness validation is also appropriate at write time.

### 5. The DDL comments are not foreign keys

The proposed schema describes `trait_name` and `locus_id` as foreign keys but does
not declare them as constraints. Today `genome_effects` is created in `open_pop()`
before `genome_meta` exists (`R/open_pop.R:285-295`). The plan should state whether
these are SQL constraints or transactionally enforced R invariants and adjust table
creation order if actual foreign keys are desired.

### 6. The upstream multiallelic audit should be reworded

`UTINYINT` itself does not prevent storing allele codes greater than one. The actual
block is that founder generation, allele frequencies, dosage calculation, and effect
semantics are biallelic. For example, `add_dosage()` currently computes dosage as
`SUM(h.allele)` (`R/add_dosage.R:126-130`), which is not a multiallelic genotype
representation. The conclusion that multiallelic support is a larger upstream
project is correct; the datatype-specific rationale is not.

## Recommended minimal v4 shape

The following is a direction, not final DDL. Exact constraints should be written
after the evaluator semantics are fixed.

### Required configuration tables

```sql
CREATE TABLE genome_effects (
  id_genome_effect INTEGER PRIMARY KEY,
  trait_name       VARCHAR NOT NULL,
  effect_name      VARCHAR,
  genome_value     DOUBLE NOT NULL
);

CREATE TABLE genome_effect_members (
  id_genome_effect INTEGER NOT NULL,
  member_slot      INTEGER NOT NULL,
  locus_id         INTEGER NOT NULL,
  contrast_name    VARCHAR NOT NULL,
  dosage_value     INTEGER,
  base_allele_freq DOUBLE,
  PRIMARY KEY (id_genome_effect, member_slot)
);
```

Suggested initial contrast vocabulary:

- `additive`: sum `(allele - p)` across eligible copies; naturally ploidy-aware;
- `dominance`: exact documented diploid coding;
- `indicator`: one for an exact local dosage state, otherwise zero.

`model_name` and `effect_class` should be included only after their executable and
validation semantics are specified.

### Conditional third table

If line- and cross-specific dominance or epistasis are required, add one
origin-composition table. Do not finalize its columns independently of its resolver.
The schema and resolver truth table are one design decision, because stored origin
rows without precedence and matching semantics are not usable model definitions.

If those effects are deferred, retain the scalar `line_name` only as a strictly
validated order-one copy-additive shortcut and document the limitation plainly.

### Required API and execution pieces

1. `define_genome_effects()` as the only documented general writer.
2. A member evaluator that reduces haplotype rows to one value per
   individual/effect/member before multiplying members.
3. A model evaluator that sums terms into a clearly named total genetic value.
4. A separate additive breeding-value path where TBV semantics are required.
5. An explicit model-selection rule if `model_name` is retained.
6. Views for readable terms and causal loci.

## Required acceptance gates

Storage-only fixtures are insufficient. Each fixture must be evaluated against a
hand-computed expected value.

At minimum, test:

1. Current order-one copy-additive values, including per-copy line fallback.
2. A hemizygous locus proving centering uses realized copy count rather than `2p`.
3. Exact diploid dominance values at dosages 0, 1, and 2.
4. A pairwise additive-by-dominance term.
5. A three-way term.
6. A complete and a sparse 3 x 3 indicator surface.
7. Functional versus Cockerham coding including the required intercept.
8. `{A, NULL}` origin handling.
9. A/B and reciprocal origin cases if the origin table is included.
10. Canonical equality for reversed input member order.
11. Duplicate logical-term rejection.
12. Transaction rollback after any invalid member.
13. Model selection and replacement if `model_name` is retained.
14. Proof that an injected dominance or epistatic term changes the computed genetic
    value and, when requested, the resulting phenotype.

The existing independent additive oracle in `tests/testthat/test-add_tbv.R:16-40`
should remain a committed mathematical regression test. It is not prohibited legacy
golden output; it independently checks the current formula.

## Decision summary

| v3 decision | Review disposition |
|---|---|
| Term plus member tables | **Keep** |
| Indicator basis for arbitrary genotype-state functions | **Keep, with a precise scope claim** |
| Separate surface/state/contrast registry families | **Do not restore** |
| `base_allele_freq` copied onto frequency-dependent members | **Keep** |
| Derived effect order | **Keep** |
| Scalar `line_name` for every contrast and interaction | **Reject** |
| Defer all origin scopes while claiming effects across lines | **Reject or narrow the claim** |
| `'additive' = dosage - 2p` | **Replace with realized-copy-count coding** |
| Unspecified Cockerham dominance label | **Specify exact values and domain** |
| Both `'allele'` and `'additive'` | **Reconsider; one ploidy-aware additive contrast may suffice** |
| `effect_class` as a validated materialized field | **Redesign or derive** |
| `model_name` without model selection | **Remove or wire through completely** |
| Duplicated `locus_name` in members | **Remove; expose through a view** |
| Raw table insertion as custom-effect workflow | **Reject; add `define_genome_effects()`** |
| Storage representability as Phase A gate | **Replace with numerical evaluation fixtures** |
| Two tables under a restricted line scope | **Acceptable** |
| Three tables with a complete origin resolver | **Preferred for the stated across-lines goal** |

## Final recommendation

Adopt the v3 term/member architecture, but revise the plan before implementation.
The minimal durable design is two effect-definition tables plus, only if required by
the intended first release, one origin-composition table. The difficult part is no
longer table normalization; it is defining evaluation, line-origin precedence, and
the user-facing writer precisely enough that a stored effect always has an observable
and reproducible meaning.


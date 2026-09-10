# Review of `plans/update_genome_effects.md`

**Review version:** v1

**Date:** 2026-09-04
**Scope:** effect-definition storage only. Simulation, TBV calculation, phenotype
assembly, variance targeting, and individual-level output storage are discussed
only where they expose a schema requirement.

## Executive verdict

The plan identifies the correct central abstraction: an effect is not always a
property of one locus, so the package needs a term header plus a variable-length
member list. That is the right replacement for the current single
`locus_name` column, and wide `locus_name_2`, `locus_name_3`, ... columns should
be rejected.

I would **not implement the plan's locked schema as written**, however. It
generalizes the number of loci in a term, but it still hard-codes the state
space at each locus to one biallelic dosage and treats free-text contrast names
as if they defined contrast functions. It also leaves the multi-line problem
on a single nullable `line_name`, even though a dominance effect can belong to
an A/B line-origin pair and an epistatic effect can involve a different
line-origin composition at each locus. Those limitations are exactly the kind
that would force another structural rewrite for polyploids, multiallelic loci,
or cross-specific effects.

The durable design should have two supported representations:

1. **Contrast terms:** coefficient times a product of explicitly defined
   single-locus contrasts. This is compact for additive, dominance, and named
   epistatic decompositions of any order.
2. **Genotype-value surfaces:** an arbitrary value for each normalized
   multi-locus genotype state. This is the model-free escape hatch for any gene
   action that is awkward or impossible to name as additive/dominance
   contrasts.

Both representations need normalized genotype states, immutable reference
population identities, and a line-origin scope that can represent tuples or
multisets of lines rather than one line. Labels such as `add_dom` should be
descriptive classifications, not the definition or identity of an effect.

The recommendation below preserves the plan's good hyperedge core but adds the
missing semantic layers. It also removes `ind_genetic_value`, phenotype
assembly, and variance-component renaming from this schema change; those are
downstream projects, not effect-definition storage.

## 1. What the package does today

### 1.1 Physical schema

`R/open_pop.R:285-295` creates:

```sql
CREATE TABLE genome_effects (
  id_genome_effect   INTEGER PRIMARY KEY,
  locus_name         VARCHAR NOT NULL,
  line_name          VARCHAR,
  trait_name         VARCHAR NOT NULL,
  genome_effect_type VARCHAR NOT NULL,
  genome_value       DOUBLE  NOT NULL,
  base_allele_freq   DOUBLE
)
```

The current row grain is intended to be:

```text
one locus x one trait x one effect type x one applicable line
```

In practice, only `genome_effect_type = 'additive'` is written or consumed.
There are no database foreign keys, checks on frequency/value domains, or a
logical unique constraint. The package relies on writer behavior and R-side
validation.

The table is created before `genome_meta`, so a SQL foreign key to a locus is
not currently convenient. That does not prevent an R-enforced relationship,
but it is one reason the current DDL has none.

### 1.2 Writer behavior

`define_additive_effects()` is the only writer.

- It selects loci from a collected `genome_meta` query and orders them by
  `locus_id` when available (`R/define_additive_effects.R:229-250`).
- It computes one alternate-allele frequency per selected locus from founder
  haplotypes or current individual haplotypes
  (`R/define_additive_effects.R:553-644`).
- It optionally rescales sampled coefficients using
  `sum(2 * p * (1 - p) * a^2)`
  (`R/define_additive_effects.R:534-542`).
- Before writing, it deletes every existing additive row for the same trait and
  `line_name`, then appends one row per selected locus
  (`R/define_additive_effects.R:297-322` and `:497-517`).
- The multi-trait `method = 'union'` path reads QTL membership from existing
  additive rows (`R/define_additive_effects.R:404-420`).

This means replacement scope is broader than a logical row key: redefining the
additive architecture for a trait/line replaces the whole locus set for that
trait/line. That behavior should become an explicit **model/version replacement
rule** in a redesigned schema, not remain an accidental consequence of a
DELETE predicate.

### 1.3 Reader behavior

`add_tbv()` is the substantive reader. For each haplotype copy it finds an
additive row at the same `locus_name`, preferring a row whose `line_name`
matches that copy's `line_origin`; it uses the common row only when no matching
line-specific row exists. It then sums:

```text
(allele - base_allele_freq) * genome_value
```

See `R/add_tbv.R:174-226`.

Important consequences:

- Effects are currently **per allele copy**, not per collapsed individual
  genotype.
- A crossbred individual can use a different coefficient and centering
  frequency for each inherited copy at the same locus.
- The meaning of `line_name` is not simply “individual belongs to this line.”
  It is “this effect applies to an allele copy with this founding line origin.”
- `NULL base_allele_freq` is silently treated as zero. That permissive behavior
  is not a good invariant for a future statistical coding.

Other readers are shallow but form part of the migration surface:

- `add_phenotype()` requires at least one common additive row for a simple
  phenotype (`R/add_phenotype.R:247-264`).
- `extract_genotypes()` accepts only a `genome_effects` table, expects it to
  contain `locus_name`, and uses those names as a QTL locus set
  (`R/extract_genotypes.R:119-126`, `:203-218`).
- The population print method counts distinct additive `locus_name` values
  (`R/tidybreed_pop.R:158-159`).
- Schema registries and descriptions enumerate all seven current columns.

### 1.4 Current biological domain

The effect schema cannot be reviewed in isolation from the genotype domain it
references:

- `ind_haplotype.allele` is documented and generated as 0/1.
- `ind_genotype.dosage_value` is the count of allele 1.
- Founders and offspring are currently constrained to organismal ploidy 2,
  although chromosome copy count can be 0, 1, or 2 and the haplotype row shape
  has a `strand` dimension.
- `line_origin` is carried per inherited copy through recombination.

Thus the present implementation is specifically a biallelic, additive,
copy-based model. `genome_effect_type` being a VARCHAR did not by itself reserve
enough structure for dominance or epistasis.

### 1.5 Dominance is not “end-to-end” today

The plan says dominance variance-component plumbing works end to end. That is
too strong. `define_effect_cov_matrix()` recognizes the literal string
`'dominance'` and will write a covariance matrix to `trait_var_comp`
(`R/define_effect_cov_matrix.R:118-148`), and generic lookup helpers can retrieve
an explicitly named matrix. No effect writer, effect reader, simulation path,
phenotype path, or result table uses that matrix. The accurate statement is:

> The covariance table can already store a matrix labeled dominance, but no
> functional dominance pipeline exists.

That distinction matters because it argues against coupling the effect-term
schema to the current variance-component naming scheme.

## 2. What the plan gets right

### 2.1 The term/member split is the correct foundation

A variable-length member table is the clean relational representation of a
set/product of locus factors. It supports order 1, order 2, and higher order
without adding columns. The plan is right to reject wide slots and right that
an order-1 effect is structurally the same family of object as a higher-order
effect.

### 2.2 A direct genotype-value representation is necessary

Not every useful genotype-to-value map should be forced through a named
additive/dominance decomposition. A direct table is valuable for:

- a hand-entered 3 x 3 two-locus diploid table;
- threshold, duplicate-factor, complementary, and sign-epistasis models;
- validating coefficient-based representations;
- future models whose preferred contrasts have not been implemented.

The plan is right that this is an addition to, not a replacement for, compact
contrast terms.

### 2.3 Frequencies do not belong on each coefficient row

`base_allele_freq` is a property of a reference population, locus, allele, and
possibly line—not of a trait coefficient. Duplicating it across every trait
allows inconsistent snapshots and cannot naturally represent an interaction's
multiple loci. Moving it out is correct.

### 2.4 Functional and statistical parameters must not be confused

The plan is right that the stored number's meaning depends on its coding and
reference. A coefficient without its basis is incomplete data. This is one of
the most important observations in the proposal.

### 2.5 It correctly finds the multi-line non-additive ambiguity

Collapsing copies to dosage discards the line-origin composition. A diploid
dosage of one could be A/A-origin, A/B-origin, B/B-origin, or include an unknown
origin. There is no principled way to select one scalar line-specific dominance
coefficient after that information has been erased. The plan is right to make
this explicit rather than silently using `ind_meta.line_name`.

### 2.6 It protects the additive hot path conceptually

Isolating the additive migration and requiring numerical invariance is sound.
Even though this review is scoped to storage, compatibility with the current
per-copy fallback semantics is a real acceptance criterion for the new data
model.

## 3. Blocking issues in the proposed locked design

The following are schema blockers, not optional refinements.

### 3.1 `locus_contrast` names do not define a contrast

The proposed member row stores only:

```text
locus_name + locus_contrast = 'additive'/'dominance'/...
```

That says what a contrast is called, but not what value it takes for each
genotype state. The plan then puts formulas in prose and expects future code to
know the formula associated with each free-text label. This breaks for:

- different functional additive reference alleles;
- Cockerham/F2 versus NOIA versus user-defined statistical bases;
- the multiple additive and dominance degrees of freedom at a multiallelic
  locus;
- the multiple non-additive contrasts at a polyploid locus;
- non-Hardy-Weinberg populations;
- line-origin- or parent-origin-aware contrasts;
- custom contrasts entered before the package knows how to simulate them.

If arbitrary future contrast labels are allowed, their actual state-to-number
mapping must be data. Otherwise adding a contrast still requires new code and
the stored effects are not self-describing after restore.

**Required change:** add a contrast definition table and a contrast-value table
that records the numeric contrast value for normalized genotype states. A
`contrast_class` such as `additive` or `dominance` remains useful for filtering,
but is not the mathematical definition.

### 3.2 `dosage_key` is biallelic, lossy, and unnecessarily serialized

The proposed direct payload uses strings such as `'2_0'`. This has several
problems:

1. One dosage per locus identifies a genotype only when there are exactly two
   alleles and “dosage” means copies of the designated alternate allele.
2. A scalar dosage does not identify a multiallelic polyploid genotype. In a
   tetraploid, allele counts A=2/B=1/C=1 and A=2/B=2/C=0 can share the same
   dosage of one chosen allele while being different genotypes.
3. The same dosage can mean different states at different local copy numbers.
   Dosage 1 is a haploid alternate genotype, a diploid heterozygote, or one
   alternate copy in a tetraploid.
4. It cannot retain line-origin composition or reciprocal parent origin.
5. Compound keys encoded into VARCHAR require parsing and duplicate database
   structure in a string.

`(ploidy + 1)^k` is therefore the cell count only for **biallelic loci of the
same known ploidy**. It is not a general polyploid genotype representation.

**Required change:** normalize genotype states. A state has a locus, local copy
number, and rows of `(allele_code, copy_count)`. A multi-locus cell then has one
state ID per member slot. Optional line/parent-origin scope must be represented
separately or as explicit state dimensions.

### 3.3 A single `line_name` cannot represent general non-additive line effects

The proposal recognizes the ambiguity but leaves exactly one nullable
`line_name` on the term. That can express:

- common;
- line A;
- line B.

It cannot express:

- dominance specific to the A/B origin pair;
- different A/B and B/A reciprocal effects;
- an A/A genotype at locus 1 interacting with an A/B genotype at locus 2;
- a tetraploid A/A/B/B line-origin composition;
- one line condition per member of a higher-order interaction.

The suggested compute restriction (`line_name IS NULL` for order >= 2) avoids a
wrong calculation, but it does not make the schema future-proof. Lifting that
restriction later would require new columns/tables—the exact outcome this
project is meant to prevent.

The proposed order-1 heuristic is also not a sufficient definition. In its SQL
sketch, `COUNT(DISTINCT line_origin)` ignores NULL. An origin set `{A, NULL}`
would report one distinct non-NULL line and `MIN(line_origin) = A`, potentially
masquerading as a within-A genotype.

**Required change:** replace the effect's scalar line with a non-null scope ID.
The scope must be able to hold a normalized line-origin multiset per locus
member, with optional parent-origin roles for reciprocal/imprinting models.
Provide a distinguished common scope instead of using a nullable key.

### 3.4 The coding distinction is too coarse and stored at the wrong grain

`effect_coding = 'functional' | 'statistical'` does not uniquely identify a
basis. “Statistical” could mean an F2 coding, a Hardy-Weinberg Cockerham coding,
a NOIA coding based on observed genotype frequencies, or a custom orthogonal
basis. Polyploid and multiallelic bases add more choices.

The formulas in the plan are a valid biallelic coding under HWE, but calling
them “NOIA-Cockerham” is misleading. General NOIA constructs design matrices
from genotype frequencies and supports departure from Hardy-Weinberg
proportions; allele frequency `p` alone is insufficient. See Álvarez-Castro and
Carlborg's unified model and the later multiallelic extension:

- [A unified model for functional and statistical epistasis](https://pmc.ncbi.nlm.nih.gov/articles/PMC1894581/)
- [Multiallelic models of genetic effects and variance decomposition](https://pmc.ncbi.nlm.nih.gov/articles/PMC3247674/)

Coding is normally shared by a coherent model/basis or by a contrast
definition, not independently repeated on every coefficient term. Repetition
permits one interaction term to claim a different basis from its main effects
without an intentional model boundary.

**Required change:** use an open `coding_name`/version as provenance and store
the actual contrast-state values. Attach the immutable reference population to
the contrast/model definition. Do not make two broad strings carry the
semantics.

### 3.5 `genome_effect_type` is not safely derivable as proposed

The plan says to derive the label by joining member labels in slot order. This
is sound only if slots are canonical. Without a canonical locus order, these
two rows describe the same product but receive different labels:

```text
slot 1: locus j, additive    slot 2: locus k, dominance
slot 1: locus k, dominance  slot 2: locus j, additive
```

Multiplication is commutative; user input order cannot define distinct biology.
A x D and D x A are distinct only **relative to a canonical ordering of the
loci**. The lower-`locus_id` member must always be slot 1, for example.

There are further problems:

- A free label such as `dom_2` contains the proposed delimiter, so concatenated
  type strings are not reliably parseable.
- Multiple contrasts can have the same broad class in multiallelic/polyploid
  systems. A class sequence does not identify the actual bases.
- The plan says all `'genotype'` member contrasts derive
  `'genotype_table'`; ordinary concatenation would instead produce labels such
  as `genotype_genotype`. That is already a special case.
- `effect_order` and `slot` are `UTINYINT`, which technically contradicts “all
  orders” by imposing a limit of 255. There is no benefit to that limit in
  small configuration tables.

**Required change:** canonicalize members by stable `locus_id`. Keep a
materialized `effect_class` only as a convenience validated from broad contrast
classes. Never use it as the term's identity or as the source of its formula.
Derive order from member count or store it only in a validated view/cache.

### 3.6 `base_name` is not an immutable reference population

A user string such as `'current_pop'` does not identify *which* current
population, at what generation/replicate, with which filter, or at which point
in time. If rows under the same name are overwritten, old statistical effects
change meaning. If they are appended, the logical key collides.

The proposed frequency table also remains biallelic: it has no allele code.
For a multiallelic locus it cannot store a frequency vector. For general
statistical coding away from HWE, allele frequencies alone do not determine the
genotype-state distribution. For exact multilocus orthogonality in the
presence of LD, even marginal single-locus distributions may not be enough.

**Required change:** create an immutable reference header with an integer ID.
Store allele frequencies by allele and, where needed, genotype-state
frequencies. The already materialized contrast values are the authoritative
evaluation semantics; reference frequencies provide provenance and permit
validation/transformation.

### 3.7 The proposed direct table is an invalid tagged union

In the term header, `genome_value` is non-NULL for contrast terms and NULL for
genotype-table terms, while another table holds values for the latter. This
makes a `genome_effects` row mean two different things and permits empty or
partially populated objects unless several cross-table rules are perfect.

The design also does not state whether a 3 x 3 table contains total genotypic
values, deviations from a mean, or an epistatic residual to add on top of main
effects. Those interpretations produce different totals.

**Required change:** keep contrast terms and genotype surfaces as distinct,
strongly typed payloads under a shared model header. Every numeric payload
column should remain NOT NULL. Give a surface explicit additive semantics such
as `component` (one matched cell contributes to the model sum) and a declared
missing-cell action.

### 3.8 The plan omits the zero-order term/intercept

A complete genotype-to-value model has a reference/mean. An arbitrary genotype
table may contain absolute values, while contrast terms typically contain
deviations. The hyperedge definition starts at one member and offers no
zero-order coefficient.

**Required change:** allow an explicitly typed, scope-aware intercept term with
zero members. A model-header scalar would not be enough for line- or
cross-specific intercepts. Do not infer it from `trait_meta.target_add_mean`;
that field has different existing semantics.

### 3.9 `locus_name` should not remain the internal relationship key

The package already treats `genome_meta.locus_id` as its internal order/join
key and `locus_name` as the user-facing identifier (`R/define_genome.R:270-289`).
The proposed member table keeps only `locus_name`, even though member
canonicalization requires a stable order and the new tables will already need
joins.

**Required change:** store `locus_id` in effect member/state tables. Join to
`genome_meta` to display names. If profiling later proves that a denormalized
name materially helps a hot query, it can be added and validated, but it should
not be the sole relationship key.

### 3.10 The breeding-value discussion is incomplete when epistasis exists

For a functional one-locus diploid model, the average substitution effect
`alpha = a + d(q - p)` is appropriate under the stated assumptions. Once
functional epistasis exists, average effects can also depend on interaction
coefficients, allele/genotype frequencies at other loci, and LD. Switching
`add_tbv()` to only the one-locus alpha formula does not generally turn a
functional epistatic genotypic model into a true breeding value.

This is primarily a future computation issue, but it reinforces the storage
requirement: preserve the full model basis and immutable reference context.
Do not encode the promise that an `additive`-classified coefficient is always
the breeding value.

## 4. Important scope problems in the plan

The source plan mixes four projects:

1. storage of effect definitions;
2. simulation/evaluation of those definitions;
3. storage of individual genetic-value outputs;
4. variance-component and phenotype-component APIs.

That breadth makes it harder to tell which decisions must be locked now.

### Move out of this schema project

- `ind_genetic_value` and whether it includes an additive component;
- the future non-additive `add_tbv()`/genotypic-value dispatcher;
- `phenotype_components.genome_effect_types` activation;
- the precise dominance line-fallback calculation;
- variance targeting and samplers;
- renaming `trait_var_comp.effect_name` values;
- a prefix rule that routes all `gen_*` names.

These topics deserve follow-up plans after the definition schema is stable.
The prefix routing proposal is particularly risky: `effect_name` is currently
an open user-facing label, and lexical prefix is not a relational type. A
separate `effect_class`/destination argument or normalized variance-component
type is safer.

### Keep as schema requirements, not implementations

- order-1 additive data must migrate without loss;
- copy-specific line-origin applicability must remain representable;
- arbitrary multi-locus genotype surfaces must be storable;
- local copy number and future polyploid states must be representable;
- effect definitions must remain interpretable after database restore;
- a complete model/version must be atomically replaceable.

## 5. Recommended data model

The following is an illustrative logical schema. Exact names can change, but
the separations and invariants should not.

### 5.1 Model header

```sql
CREATE TABLE genome_effect_models (
  id_genome_effect_model INTEGER PRIMARY KEY,
  model_name             VARCHAR NOT NULL,
  trait_name             VARCHAR NOT NULL,
  model_version          INTEGER NOT NULL DEFAULT 1,
  notes                  VARCHAR
);
```

Logical key: `(trait_name, model_name, model_version)`.

This groups all terms/surfaces that form one coherent architecture. It solves
the current ambiguous DELETE scope and allows a future caller to compare or
retain multiple architectures without adding columns. `model_name = 'default'`
is enough for the existing API.

Intercepts are zero-member rows in `genome_effects`, not a model-header value,
so they can use the same common, line-specific, or cross-specific scopes as
other terms.

Whether multiple models can be active at once is an API decision. Do not add an
`is_active` flag unless runtime selection actually needs persistent state.

### 5.2 Immutable reference populations

```sql
CREATE TABLE genome_effect_references (
  id_genome_effect_reference INTEGER PRIMARY KEY,
  reference_name             VARCHAR NOT NULL,
  source_type                VARCHAR NOT NULL,
  source_description         VARCHAR,
  created_generation         INTEGER,
  created_replicate          INTEGER,
  notes                      VARCHAR
);

CREATE TABLE genome_effect_reference_alleles (
  id_genome_effect_reference INTEGER NOT NULL,
  locus_id                   INTEGER NOT NULL,
  allele_code                VARCHAR NOT NULL,
  line_name                  VARCHAR,
  allele_freq                DOUBLE NOT NULL
);
```

R-side logical key:
`(id_genome_effect_reference, locus_id, allele_code, line_name)` with NULL
normalized. Validate `0 <= allele_freq <= 1` and per-locus/context sums.

Optionally add reference genotype-state frequencies when a statistical basis
requires them:

```sql
CREATE TABLE genome_effect_reference_states (
  id_genome_effect_reference INTEGER NOT NULL,
  id_genome_state            INTEGER NOT NULL,
  line_name                  VARCHAR,
  state_freq                 DOUBLE NOT NULL
);
```

References are immutable once a contrast/model uses them. Creating another
snapshot creates another ID even if the user-facing name is similar.

### 5.3 Normalized local genotype states

```sql
CREATE TABLE genome_states (
  id_genome_state INTEGER PRIMARY KEY,
  locus_id        INTEGER NOT NULL,
  copy_number     USMALLINT NOT NULL,
  state_domain    VARCHAR NOT NULL
);

CREATE TABLE genome_state_alleles (
  id_genome_state_allele INTEGER PRIMARY KEY,
  id_genome_state        INTEGER NOT NULL,
  allele_code            VARCHAR NOT NULL,
  line_name              VARCHAR,
  parent_origin          UTINYINT,
  copy_count             USMALLINT NOT NULL
);
```

Suggested `state_domain` values are `allele`, `allele_line`, and
`allele_line_parent`. They state whether the nullable origin columns
participate in state identity; NULL must not be used as an undocumented
wildcard. The R-side logical key for count rows is
`(id_genome_state, allele_code, line_name, parent_origin)`, NULL-normalized.
The surrogate row ID allows the same allele to occur in copies from multiple
lines within one state.

Validation:

- `copy_number >= 0`;
- `copy_count > 0`;
- sum of `copy_count` equals `copy_number`;
- no duplicate normalized state at one locus;
- allele codes must be valid for the locus once the genome schema gains an
  allele dimension.

For today's origin-agnostic biallelic diploid locus, the states are `{0:2}`,
`{0:1,1:1}`, and `{1:2}`. A hemizygous state has `copy_number = 1`. A
tetraploid AABC state is `A:2, B:1, C:1`. An origin-aware state can additionally
distinguish “allele 1 from line A, allele 0 from line B” from the reverse. No
serialized dosage key is needed.

This table does not require tidybreed to simulate multiallelic or polyploid
genomes now. It merely gives effect definitions a stable state vocabulary.

### 5.4 Explicit contrast definitions

```sql
CREATE TABLE genome_contrasts (
  id_genome_contrast         INTEGER PRIMARY KEY,
  locus_id                   INTEGER NOT NULL,
  contrast_name              VARCHAR NOT NULL,
  contrast_class             VARCHAR NOT NULL,
  coding_name                VARCHAR NOT NULL,
  evaluation_unit            VARCHAR NOT NULL,
  id_genome_effect_reference INTEGER,
  notes                      VARCHAR
);

CREATE TABLE genome_contrast_values (
  id_genome_contrast INTEGER NOT NULL,
  id_genome_state    INTEGER NOT NULL,
  contrast_value     DOUBLE NOT NULL,
  PRIMARY KEY (id_genome_contrast, id_genome_state)
);
```

Examples of `contrast_class` are `additive`, `dominance`, and
`non_additive`. Examples of `coding_name` are `functional_diploid_v1`,
`cockerham_hwe_v1`, `noia_v1`, or `custom`. These names aid interpretation;
the rows in `genome_contrast_values` are authoritative.

`evaluation_unit` distinguishes a contrast evaluated once on a whole local
genotype from one evaluated per allele copy. That difference is required to
preserve today's line-specific per-copy additive semantics. A conventional
genotype-level additive contrast used inside A x D epistasis can be a separate
contrast ID of the same broad `additive` class.

A tetraploid locus can have A, D1, D2, and D3 contrasts as separate IDs without
renaming a schema-wide enum. A multiallelic locus can have several additive and
dominance contrast IDs of the same broad class. Future code can evaluate a
stored contrast without hard-coding its label.

If storage size becomes a concern, built-in basis families can be represented
by a versioned family plus parameters, but custom contrasts still need explicit
state values. Given that these are QTL/configuration tables, explicit values are
the safer first design.

### 5.5 General applicability scopes

```sql
CREATE TABLE genome_effect_scopes (
  id_genome_effect_scope INTEGER PRIMARY KEY,
  scope_name             VARCHAR NOT NULL,
  match_unit             VARCHAR NOT NULL,
  notes                  VARCHAR
);

CREATE TABLE genome_effect_scope_origins (
  id_genome_effect_scope_origin INTEGER PRIMARY KEY,
  id_genome_effect_scope INTEGER NOT NULL,
  member_slot            INTEGER,
  origin_slot            INTEGER NOT NULL,
  line_match_type        VARCHAR NOT NULL,
  line_name              VARCHAR,
  parent_origin          UTINYINT,
  copy_count             USMALLINT NOT NULL
);
```

The R-side logical key is
`(id_genome_effect_scope, member_slot, origin_slot)`, NULL-normalized. The
surrogate primary key avoids pretending that nullable/model-wide member scopes
can participate in a SQL primary key.

`line_match_type` makes `exact`, `any`, and `unknown` different stored
conditions. `exact` requires `line_name`; `any` and `unknown` require it to be
NULL. This avoids assigning wildcard semantics to SQL NULL by accident.

Suggested `match_unit` values:

- `common`: no origin rows; applies regardless of origin;
- `allele_copy`: current additive line behavior—apply to matching copies;
- `locus_genotype`: match the line-origin multiset at each member locus;
- `individual`: explicitly match a whole-individual/cross context if ever
  required.

The exact matching rules need a dedicated small design review before this table
is implemented. The non-negotiable point is the cardinality: scope must support
zero, one, or many line-origin conditions per locus member. One `line_name` on
the effect cannot do that.

A distinguished common scope gives every effect a non-null scope ID and avoids
NULL-sensitive logical keys. An A/B dominance scope can carry two origin rows at
member slot 1. Parent-origin values can distinguish reciprocal A-sire/B-dam
from B-sire/A-dam when desired. Tetraploid scopes use `copy_count` rather than
four repeated rows.

Scopes handle reusable applicability and fallback rules. Origin-aware
`genome_states` are available when the value depends on the *binding* between
allele and origin (for example, allele 1 from A versus allele 1 from B), which a
separate allele-only state plus an unordered A/B scope cannot express.

### 5.6 Contrast coefficient terms (`genome_effects`)

```sql
CREATE TABLE genome_effects (
  id_genome_effect       INTEGER PRIMARY KEY,
  id_genome_effect_model INTEGER NOT NULL,
  id_genome_effect_scope INTEGER NOT NULL,
  effect_class           VARCHAR NOT NULL,
  genome_value           DOUBLE NOT NULL
);

CREATE TABLE genome_effect_members (
  id_genome_effect   INTEGER NOT NULL,
  member_slot        INTEGER NOT NULL,
  locus_id           INTEGER NOT NULL,
  id_genome_contrast INTEGER NOT NULL,
  PRIMARY KEY (id_genome_effect, member_slot)
);
```

Semantics:

```text
term contribution = genome_value * product(member contrast values)
```

Core invariants:

- A coefficient term has at least one member. An `effect_class = 'intercept'`
  term has exactly zero; mathematically its empty contrast product is one.
- A locus occurs at most once in a term.
- Members are canonicalized by ascending `locus_id`; `member_slot` is assigned
  after canonicalization, never from input order.
- Each contrast's `locus_id` equals the member's `locus_id`.
- `effect_class` is materialized only for convenient filtering/display and is
  validated from the ordered broad classes (`additive`, `dominance`, etc.).
- Term identity uses the model, scope, and ordered contrast IDs—not the class
  label.
- No `effect_order` is needed in base storage. A view can expose
  `COUNT(members)` and a derived class label. If later profiling proves a cached
  order valuable, add and validate it then.

Canonical ordering makes A(lower locus) x D(higher locus) distinct from
D(lower locus) x A(higher locus), without treating arbitrary input order as
biology.

### 5.7 Direct genotype-value surfaces

```sql
CREATE TABLE genome_effect_surfaces (
  id_genome_effect_surface INTEGER PRIMARY KEY,
  id_genome_effect_model   INTEGER NOT NULL,
  id_genome_effect_scope   INTEGER NOT NULL,
  surface_name             VARCHAR NOT NULL,
  missing_state_action     VARCHAR NOT NULL,
  default_value            DOUBLE,
  notes                    VARCHAR
);

CREATE TABLE genome_effect_surface_members (
  id_genome_effect_surface INTEGER NOT NULL,
  member_slot              INTEGER NOT NULL,
  locus_id                 INTEGER NOT NULL,
  PRIMARY KEY (id_genome_effect_surface, member_slot)
);

CREATE TABLE genome_effect_surface_cells (
  id_genome_effect_cell    INTEGER PRIMARY KEY,
  id_genome_effect_surface INTEGER NOT NULL,
  genome_value             DOUBLE NOT NULL
);

CREATE TABLE genome_effect_surface_cell_states (
  id_genome_effect_cell INTEGER NOT NULL,
  member_slot           INTEGER NOT NULL,
  id_genome_state       INTEGER NOT NULL,
  PRIMARY KEY (id_genome_effect_cell, member_slot)
);
```

Semantics: for each surface, exactly one matching cell contributes its
`genome_value` to the model. Multiple surfaces can be added as components.
Common or scope-specific intercepts are zero-member coefficient terms.

Validation:

- Surface members are canonicalized by `locus_id`.
- Every cell has exactly one state for every surface member.
- Each state belongs to the member's locus.
- Cell tuples are unique.
- `missing_state_action` is explicit (`error`, `zero`, or `default`); `default`
  requires a non-NULL `default_value`.
- A writer can require a complete Cartesian table or intentionally allow a
  sparse table. Missing cells must never be silently interpreted.

For the current diploid biallelic two-locus case this stores the same nine cells
as the plan's `dosage_key`, but it also supports unequal local copy number,
polyploidy, and multiallelic states without a new column or string grammar.

## 6. Why two representations are better than one nullable union

Contrast terms and genotype surfaces answer different questions:

| Property | Contrast term | Genotype surface |
|---|---|---|
| Compact additive architecture | Excellent | Wasteful |
| Named A/D/epistatic decomposition | Explicit | Not intrinsic |
| Arbitrary gene action | Only with a full basis | Exact |
| User-defined polyploid basis | Via contrast-state rows | Directly |
| Variance-component classification | Available as metadata | Requires decomposition |
| Missing genotype states | Contrast definition validation | Surface missing policy |

Putting both behind a shared model header allows them to coexist without
making `genome_value` nullable or overloading one row type. For example, a model
can contain compact additive main effects plus one hand-entered two-locus
surface representing a special interaction residual.

## 7. Examples proving the schema's range

These examples describe stored objects, not promised compute functions.

### 7.1 Existing additive effect

- Model: trait ADG, name `default`; an optional common zero-member term stores
  its intercept.
- Contrast at locus 10: current centered allele-copy coding over copy-number-1
  states, with the relevant immutable reference ID and
  `evaluation_unit = 'allele_copy'`.
- Scope: `allele_copy`, line A.
- Effect: coefficient 0.3.
- Member: slot 1, locus 10, that additive contrast ID.

This can reproduce current per-copy line behavior. A common fallback is a
separate effect with common scope. Precedence belongs in a future resolver.

### 7.2 Diploid dominance

- Contrast values for locus 10 over states 0/0, 0/1, 1/1 are 0, 1, 0 for a
  functional heterozygote indicator.
- One coefficient term refers to that contrast.
- Common dominance uses common scope.
- A/B-specific dominance uses a `locus_genotype` scope with one A-origin copy
  and one B-origin copy at member slot 1.

No meaning is inferred from the word `dominance`; the three stored state values
define it.

### 7.3 Pairwise A x D epistasis

- Members are canonicalized by locus ID.
- Lower-ID locus uses an additive contrast; higher-ID locus uses a dominance
  contrast.
- `effect_class = 'add_dom'` is a validated display/filter value.
- The two contrast IDs, not the label, identify the term.

Reversing input locus order produces the same stored term. Applying dominance
to the lower-ID locus and additive to the higher-ID locus produces the distinct
`dom_add` term.

### 7.4 Autotetraploid locus

For a biallelic tetraploid locus, define its five normalized states and as many
independent contrast IDs as the chosen basis requires. Store the numeric value
of every contrast at every state. The effect tables need no new columns and do
not need to agree in advance on labels such as `dom_2`, `dom_3`, and `dom_4`.

This matters because polyploid contrast systems are basis-dependent. A free
VARCHAR reserves labels; explicit state values reserve semantics.

### 7.5 Arbitrary two-locus table

- Create one surface with two canonical member loci.
- Insert one cell per pair of local genotype-state IDs.
- Store each requested component value.
- Set `missing_state_action = 'error'` for a complete manual table.

The same four tables hold 3 x 3 diploid, 5 x 5 tetraploid, mixed-copy-number,
or multiallelic tables.

### 7.6 Reciprocal multi-line effect

Create two scopes:

- locus member 1: parent 1 / line A and parent 2 / line B;
- locus member 1: parent 1 / line B and parent 2 / line A.

Attach different coefficient terms or surfaces to those scopes. This stores a
reciprocal effect without adding `line_name_1`, `line_name_2`, or a future
reciprocal-cross column.

## 8. Treatment of base/reference frequencies

The plan is correct to normalize frequencies but should revise the ownership
model:

- A reference is a snapshot, not a mutable name.
- Biallelic alternate-allele frequency is the current special case of an
  allele-frequency vector.
- Statistical contrast definitions point to the reference ID used to construct
  them.
- The numeric contrast-state values are stored, so evaluation does not depend
  on re-running a future formula or reconstructing an old population.
- Genotype-state frequencies are stored when the chosen basis needs departure
  from HWE represented.
- Joint multilocus frequencies should be added only for models/operations that
  require them; do not claim full epistatic orthogonality under LD from marginal
  allele frequencies alone.

For current additive migration, create reference snapshots from the exact
`base_allele_freq` values already attached to rows. If identical current rows
disagree for the same `(base context, locus)`, migration should stop and report
the conflict instead of choosing one.

## 9. Multi-line semantics: what must be decided now

The final evaluator can be deferred, but storage must distinguish these
concepts now:

1. **Allele-copy origin scope.** This is what current additive `line_name`
   actually means.
2. **Within-locus origin composition.** Needed for dominance in A/A, A/B, and
   B/B origin pairs, and higher-ploidy multisets.
3. **Per-member origin composition.** Needed when each locus in an epistatic
   term has a different composition.
4. **Parent-origin role.** Needed for reciprocal crosses and imprinting.
5. **Common fallback.** Must be a deliberate scope and precedence rule, not a
   side effect of SQL NULL logic.

The plan's proposed “both copies same line, otherwise common” rule is a
reasonable first evaluator policy, but it should not define the storage limit.
It also should not be described as a complete heterosis model: it permits
heterosis caused by directional dominance and divergent allele frequencies,
but it cannot represent an explicit A/B-specific dominance deviation while the
schema has only one line.

Before implementation, write a truth table for common, A-specific, B-specific,
A/B-specific, reciprocal A/B, unknown origin, haploid, diploid, and tetraploid
cases. The scope tables should be accepted only if all rows have an unambiguous
stored representation.

## 10. Views and user-facing ergonomics

The normalized base tables should not force every user to write six joins.
Provide read-only views/helpers after the storage design is stable:

- `genome_effect_terms`: one row per coefficient term with trait/model, order,
  effect class, coefficient, and a human-readable member description;
- `genome_effect_loci`: one row per model/effect/locus, suitable for QTL-set
  extraction;
- `genome_effect_surfaces_long`: surface/cell/member states in tidy long form;
- optionally an order-1 compatibility view during migration.

`extract_genotypes()` should eventually accept the locus view or any filtered
table containing `locus_id`; it should not require that the source table be
named exactly `genome_effects`. That existing nominal type check makes the
schema harder to evolve.

QTL terminology should also be precise:

- a locus is causal for a model if it appears in any coefficient term or
  surface;
- an additive QTL appears in an order-1 additive-class term;
- an epistatic-only locus can be causal without an order-1 row.

## 11. Required validation and transaction rules

Because several constraints span parent/child tables, they will be R-enforced.
Every public write must be one DuckDB transaction and run a single comprehensive
validator before commit.

Minimum validator checks:

### Model/reference

- logical keys unique;
- referenced traits/loci/scopes/references exist;
- references used by contrasts are immutable;
- numeric values are finite unless missingness is explicitly meaningful.

### States/contrasts

- state allele counts sum to local copy number;
- no duplicate normalized state;
- state origin columns agree with `state_domain`;
- contrast states belong to the contrast locus;
- required state coverage is complete for built-in contrasts;
- statistical contrasts have a valid immutable reference;
- reference allele frequencies sum correctly within locus/context.

### Terms

- every non-intercept term has members and every intercept has none;
- member loci are unique and canonical;
- contrast locus agrees with member locus;
- no duplicate term identity within model/scope;
- materialized effect class agrees with ordered contrast classes.

### Surfaces

- members are unique and canonical;
- cells contain exactly the member slots of their surface;
- cell states agree with member loci;
- no duplicate state tuple;
- completeness or sparsity agrees with `missing_state_action`.

### Scopes

- common scopes have no origin rows;
- copy counts are positive;
- member slots exist on the attached term/surface when slot-specific;
- parent-origin values are valid;
- normalized origin tuples are unique;
- `line_match_type` and `line_name` agree;
- unknown/wildcard semantics are explicit and never inferred from a missing
  SQL comparison.

## 12. Migration from the current table

The migration should be lossless and deliberately narrower than the source
plan's six-step implementation roadmap.

### Phase A: finalize semantics with fixtures

Before DDL, build small table fixtures for:

- common additive effect;
- line-A additive plus common fallback;
- line-A and line-B additive effects with different base frequencies;
- an A/B dominance scope;
- an A x D pair;
- a tetraploid contrast;
- a direct 3 x 3 surface;
- a reciprocal line effect.

The proposed schema must represent all fixtures without serialized compound
keys or special columns.

### Phase B: add definition/reference tables

Create the new tables and registries. Keep downstream readers on the old table
temporarily if needed. Populate built-in biallelic states and explicit additive
contrast definitions.

### Phase C: migrate additive rows

For every existing trait/line architecture:

1. create a default model/version;
2. create immutable reference snapshot rows from stored base frequencies;
3. create common or allele-copy line scopes;
4. create order-1 coefficient terms with one canonical member each;
5. preserve exact coefficients;
6. validate counts and every old-to-new mapping.

Do not infer that rows with the same frequency necessarily came from the same
base selection. When provenance is unknowable, label the migrated reference as
legacy and preserve the values rather than inventing provenance.

### Phase D: switch readers through views

First expose an order-1 flattened view and verify the current independent TBV
oracle against it. Then update QTL extraction and summaries to use locus views.
This contains the migration blast radius while keeping the normalized storage
authoritative.

### Phase E: remove old columns/table shape

Only after parity checks should the old row layout disappear. This is pre-1.0,
so a clean break is reasonable, but a staged development migration makes errors
easier to localize.

Non-additive computation, output tables, phenotype inclusion, and samplers are
separate phases after this project.

## 13. Acceptance criteria for the storage design

The schema is ready only if all of the following are true:

- Any number of distinct loci can belong to one coefficient term.
- Member order is canonical and cannot create duplicate biology.
- A contrast's numeric meaning is recoverable entirely from stored data.
- Multiple contrasts of the same broad class can exist at one locus.
- Local genotype state records copy number, an arbitrary allele-count vector,
  and optionally the allele-to-line/parent-origin binding.
- Direct genotype surfaces contain no parsed compound key.
- A term/surface can be common, copy-line-specific, cross-line-specific,
  reciprocal, or have per-locus line-origin composition without schema change.
- Statistical effects reference an immutable base snapshot.
- Current biallelic additive rows migrate without coefficient or frequency loss.
- The schema permits hemizygous and polyploid local states.
- Missing surface cells have explicit behavior.
- Contrast coefficients and direct surface values cannot be confused through
  nullable payload columns.
- The schema supports common and scope-specific intercept/reference values.
- QTL loci are discoverable through a normalized view regardless of effect
  order or representation.
- All writes and model replacements are atomic.

## 14. Section-by-section disposition of the source plan

| Source-plan proposal | Disposition | Reason |
|---|---|---|
| Term header + member list | **Keep** | Correct hyperedge abstraction |
| `effect_order` stored | **Prefer derived view** | Member count is authoritative; cache only if measured |
| `slot UTINYINT` | **Change to INTEGER** | Avoid artificial all-order ceiling |
| `locus_name` in members | **Change to `locus_id`** | Stable canonical relationship key |
| Free `locus_contrast` label | **Replace/augment with contrast ID** | A label does not define state values |
| Derived/materialized `genome_effect_type` | **Keep only as validated class** | Useful filter, unsafe identity/semantics |
| `effect_coding` functional/statistical | **Replace with versioned coding provenance + stored values** | Too coarse and repeated at wrong grain |
| `base_name` + `base_line_name` | **Replace with immutable reference ID** | Names do not identify snapshots |
| `genome_base_freq` | **Generalize by allele and optional state frequencies** | Current proposal remains biallelic/HWE-limited |
| `dosage_key` | **Reject** | Serialized, biallelic, loses copy number/origin |
| Direct genotype table | **Keep as normalized surface tables** | Essential general escape hatch |
| Single header `line_name` | **Reject for final schema** | Cannot represent line tuples/multisets/reciprocals |
| Common fallback | **Keep as explicit scope policy** | Valid behavior, should not rely on NULL logic |
| `ind_genetic_value` | **Move to later plan** | Simulation output, not effect-definition storage |
| `trait_var_comp` prefix rewrite | **Move/reconsider** | Lexical prefix is not a relational type |
| Compute query sketches | **Move to evaluator plan** | Useful requirements, premature implementation |
| Functional default | **Do not lock here** | Storage supports multiple explicit bases |
| `define_genome_effects()` | **Likely keep** | Appropriate general writer after input shape is designed |

## 15. Recommended decision

Adopt the source plan's **unified term/member direction**, but reopen the
“locked design” before implementation.

The minimum changes needed to make it genuinely durable are:

1. group effects into named/versioned models;
2. use stable `locus_id` and canonical member order;
3. define contrasts by IDs plus stored genotype-state values;
4. normalize genotype states instead of using `dosage_key`;
5. replace scalar `line_name` with normalized origin-aware scopes;
6. replace mutable base names with immutable reference IDs;
7. store direct genotype surfaces in strongly typed, separate tables;
8. keep individual outputs and evaluator behavior out of this schema project.

If the team wants a smaller first implementation, the safe subset is model
header + coefficient terms + members + immutable references + biallelic state
and contrast tables. It is acceptable to postpone direct surfaces and complex
scope writers **only if their final table shapes are agreed and reserved now**.
It is not safe to ship the proposed `dosage_key` or single-line schema as the
general endpoint.

The term/member split prevents a rewrite for higher interaction order. The
additional state, contrast, reference, and scope normalization prevents the
next rewrites—for polyploid bases, multiallelic genotypes, and multi-line gene
action. That is the level at which this schema should be considered general.

---

# Response to this review

**Responder:** tidybreed maintainer (with Claude)
**Date:** 2026-09-04
**Result:** `plans/update_genome_effects.md` rewritten as **v2**. The review is
accepted on almost every substantive point. Two items rejected, three accepted
with amendment, four escalated to maintainer decision.

## Verification of factual claims

Every code claim in the review was checked against the source and holds:

| Claim | Verified |
|---|---|
| `extract_genotypes()` enforces a nominal `table_name == "genome_effects"` check | ✅ `R/extract_genotypes.R:124-127` |
| `genome_effects` is created before `genome_meta`, so a SQL FK to a locus is not currently possible | ✅ `open_pop.R:286` vs. `define_genome.R:272` |
| `genome_meta.locus_id` is the internal order/join key | ✅ `define_genome.R:270-289` |
| `define_additive_effects()` DELETE scope is a whole-architecture replacement | ✅ `:297-322`, `:497-517` |
| `add_tbv()` resolves effects **per allele copy**, not per genotype | ✅ `add_tbv.R:206-227` |
| `NULL base_allele_freq` silently treated as zero | ✅ `add_tbv.R:212` `COALESCE(..., 0)` |
| Only `define_effect_cov_matrix.R:118-148` recognizes `'dominance'`; nothing consumes it | ✅ |

## Disposition

### Accepted in full (18)

| § | Item |
|---|---|
| 1.5 | "The covariance table can already store a matrix labeled dominance, but no functional dominance pipeline exists" — adopted verbatim as the accurate wording |
| 3.1 | Contrast definition + contrast-value tables. A label does not define a contrast; state→value rows must be data |
| 3.2 | Normalized genotype states (`genome_states` + `genome_state_alleles`) replacing `dosage_key` |
| 3.3 | Origin-aware scope tables replacing scalar `line_name` — **including the `COUNT(DISTINCT line_origin)` NULL bug**, which was a genuine defect in the v1 SQL sketch: `{A, NULL}` reported one distinct line with `MIN = A` and masqueraded as a within-A genotype |
| 3.4 | Coding as versioned provenance at contrast grain, not two broad strings repeated per term; and the correction that the v1 formulas are the **HWE Cockerham** special case, not general NOIA (which builds design matrices from observed genotype frequencies and tolerates non-HWE) |
| 3.5 | Canonical member ordering by `locus_id`; `effect_class` materialized but never identity; no stored `effect_order`; `INTEGER` not `UTINYINT` |
| 3.6 | Immutable integer reference IDs; allele-keyed frequencies; reserved genotype-state frequencies for non-HWE bases |
| 3.7 | Separate strongly-typed surface tables; every payload column `NOT NULL`; explicit declared value semantics |
| 3.8 | Zero-member intercept terms, scope-aware, not inferred from `trait_meta.target_add_mean` |
| 3.9 | `locus_id` as the member relationship key |
| 3.10 | Breeding value under epistasis: the one-locus α is insufficient once functional epistasis exists. Storage implication adopted — do not encode the promise that an `additive`-classified coefficient is the breeding value |
| 5.1 | Model/version header, turning today's accidental DELETE scope into explicit atomic model replacement |
| 5.4 | `evaluation_unit` (`allele_copy` vs `locus_genotype`) — this is what preserves today's per-copy line semantics |
| 5.5 | `line_match_type` as explicit stored conditions (`exact`/`any`/`unknown`), never wildcard-by-NULL |
| 10 | Ergonomic views; removing the `extract_genotypes()` nominal type check; precise QTL terminology (causal locus vs. additive QTL vs. epistatic-only locus) |
| 11 | Comprehensive validator, one DuckDB transaction per public write |
| 12 (A) | **Fixtures before DDL.** Judged the single most valuable contribution of this review — the cheapest possible way to discover the schema is wrong |
| 13 | Acceptance criteria adopted as the plan's design invariants |

### Accepted with amendment (3)

**A1. `genome_states.locus_id` → nullable.** The review specifies `NOT NULL`.
That forces the three structurally identical biallelic diploid states
(`{0:2}`, `{0:1,1:1}`, `{1:2}`) to be written once *per locus*: 3,000 rows for a
1,000-QTL model, 150,000 for a 50K-locus polygenic model. Made nullable, where
NULL means a locus-independent template. Locus-specific states are then needed
only when allele codes are locus-specific (multiallelic). No loss of generality;
the review's own §5.3 tetraploid and origin-aware examples still work.

**A2. `genome_contrasts.locus_id` → nullable.** Same problem, sharper. A
*functional* coding is the same function at every locus (`x_A(g) = g − 1`), so
`NOT NULL` forces one contrast row per locus to store a locus-invariant
function. With nullable `locus_id`, the entire functional diploid additive basis
is **1 contrast row + 3 value rows for the whole genome**. Frequency-dependent
codings still get one contrast per locus, because `p_j` genuinely differs — and
they still reuse the three shared template states, so the state table never
amplifies in the biallelic case.

**A3. Origin double-specification rule added.** The review introduces line
origin in **two** places — inside a state (`state_domain = 'allele_line'`) and
inside a scope origin row — and acknowledges the overlap only in passing at the
end of §5.5. Leaving which is canonical unstated would reintroduce precisely the
ambiguity the redesign exists to remove ("do I express A/B dominance as a state
or as a scope?"). The plan adds an explicit rule — **state = genotype identity,
and carries origin only when the contrast value depends on the allele↔origin
binding; scope = applicability, for everything else** — plus a validator check
rejecting any member that uses both at once.

### Rejected (2)

**R1. Removing Part 2 (`§4 Move out of this schema project`) from the plan.**

The review is right that v1 mixed four projects, and right that the mixing made
it unclear which decisions were being locked. That diagnosis is **accepted** and
implemented as an explicit Part 1 (lock) / Part 2 (sketch, non-binding) split,
with every named item — `ind_genetic_value`, the evaluator, the
`phenotype_components.genome_effect_types` activation, the dominance line
fallback, samplers, variance-component renaming — labelled "own plan".

What is rejected is **deleting** the downstream content. Two reasons:

1. A storage design's acceptance criteria are unfalsifiable in isolation. "Is
   this schema sufficient?" can only be answered against how it will be read.
   The review's own §7 examples exist for exactly this reason, and its §13
   criteria ("QTL loci are discoverable", "missing surface cells have explicit
   behavior") are statements about evaluation.
2. The maintainer's request was explicitly for a way to *include* dominance and
   epistasis — the whole picture — not for effect-definition DDL alone.

The demotion, not the deletion, is what fixes the stated problem.

**R2. Free-text discriminator columns.**

§3.1 rejects a free-text `locus_contrast` on the principle that "a label does not
define state values" — correct, and accepted. But the recommended schema then
introduces seven free-text discriminators of its own: `state_domain`,
`match_unit`, `line_match_type`, `evaluation_unit`, `missing_state_action`,
`source_type`, `contrast_class`. Each is a string that code must branch on, and
each is presented with "suggested values".

The principle is right; it was applied unevenly. Resolution adopted in v2: those
columns are **closed sets validated in R**, because they name *structure the
package must branch on* and their valid set ships with the package.
`contrast_name` and `coding_name` remain open, because they are *provenance*
whose semantics live in `genome_contrast_values`. That distinction is now stated
explicitly in the schema preamble rather than left implicit.

### Escalated to maintainer decision (4)

These are genuine judgement calls, not review defects. All four are in the v2
plan's Open Questions with options and a recommendation.

| ID | Question | Recommendation |
|---|---|---|
| **D2** | Build-now scope. The design goes from 1 table to 16, each needing six registry entries. Tier 2 (15 tables, everything but reference states) / Tier 1 (11, the review's "safe subset", coefficient terms only) / Tier 0 (v1's simpler schema) | **Tier 2**, plus a registration helper built first. Surfaces are the *reason* the project exists; postponing them postpones the payoff. Tier 0 is unsafe and the review is right to say so |
| **D3** | The scope **resolver** — which scope wins when several match. §5.5 states plainly that these rules "need a dedicated small design review" | **Defer, gated on a truth table** written during Phase A (common / A / B / A-B / reciprocal / unknown / haploid / diploid / tetraploid). Accept the scope tables only if every row has an unambiguous stored representation |
| **D4** | Migrate existing `.duckdb` populations, or regenerate? `CLAUDE.md` owes no user migration pre-1.0, but the maintainer may have their own runs | **Regenerate if cheap.** If migrating: §12's rule that conflicting `base_allele_freq` values for the same `(base context, locus)` must **stop and report**, never be silently resolved, is adopted; migrated references get `source_type = 'legacy'` with no invented provenance |
| **D5** | Contrast values always materialized, or computed at read time for built-in codings? | **Always materialized**, with built-in codings as write-time generators. Computing at read time from `coding_name` reintroduces the exact failure mode §3.1 identifies: a restored database interpretable only by a version that still knows the label |

## Section 14 disposition table — line-by-line agreement

| Review's row | Review's disposition | v2 outcome |
|---|---|---|
| Term header + member list | Keep | ✅ kept |
| `effect_order` stored | Prefer derived view | ✅ removed; view |
| `slot UTINYINT` | Change to INTEGER | ✅ done |
| `locus_name` in members | Change to `locus_id` | ✅ done |
| Free `locus_contrast` label | Replace with contrast ID | ✅ done |
| Derived `genome_effect_type` | Keep only as validated class | ✅ `effect_class`, non-identity |
| `effect_coding` functional/statistical | Replace with versioned provenance + stored values | ✅ done |
| `base_name` + `base_line_name` | Replace with immutable reference ID | ✅ done |
| `genome_base_freq` | Generalize by allele + optional state freqs | ✅ done; state freqs reserved not built |
| `dosage_key` | Reject | ✅ rejected |
| Direct genotype table | Keep as normalized surfaces | ✅ done, + `value_semantics` |
| Single header `line_name` | Reject | ✅ rejected; scopes |
| Common fallback | Keep as explicit scope policy | ✅ distinguished common scope |
| `ind_genetic_value` | Move to later plan | ⚠️ demoted to Part 2, not deleted (**R1**) |
| `trait_var_comp` prefix rewrite | Move/reconsider | ✅ prefix rule rejected outright — a user naming a random effect `gen_flock` would silently route into `trait_var_comp`; lexical prefix is not a relational type. Renaming moved to its own plan |
| Compute query sketches | Move to evaluator plan | ⚠️ demoted to Part 2, not deleted (**R1**) |
| Functional default | Do not lock here | ⚠️ partially accepted — storage is basis-agnostic, but `define_additive_effects()` still needs a default when the user says nothing. Raised as **D8** rather than silently omitted |
| `define_genome_effects()` | Likely keep | ⚠️ deferred to **D6** — with normalized states/contrasts/scopes its argument list is no longer obvious; design after Phase A fixtures |

## What the review changed most

The single most consequential correction is **§3.1 + §3.2 together**: v1
generalized the *number of loci* in a term but left the *state space at each
locus* hard-coded to one biallelic dosage and the *contrast semantics* in prose.
That combination would have reserved labels rather than semantics, and would
have forced a rewrite for the first polyploid, multiallelic, or non-HWE basis —
the exact outcome the project exists to prevent. Normalized states plus explicit
contrast values fix both at once, and the amendments A1/A2 make them cheap.

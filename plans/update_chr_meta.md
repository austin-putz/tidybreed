# Refactor: `chr_meta` → explicit per-parent inheritance + recombination tables, and `define_chr()` → `define_chromosome()`

**Status**: Draft — for review (**v6**). v5 was the redesign: the single wide
`chr_meta` table is **replaced by two explicit long tables** — `chr_inheritance`
(copy count / origin, keyed by *offspring* sex) and `chr_recombination`
(recombines, keyed by *parent* sex). `copy_mode` and `hemi_parent` are removed in
favour of explicit per-parent copy counts (`from_parent_1`, `from_parent_2`),
resolving Codex blocking issues #2/#3/#5 structurally. **v6 folds in the v5 review
round** (Codex + Gemini, both approving): DB `CHECK` constraints, `≤ 2` ploidy
invariant, offspring-vs-parent line context, deterministic shadowing algorithm,
`overwrite=FALSE`, transactional rollback, `2,0`/gamete-role scope limits — see
"Response to the v5 review round." (v4 kept one table with `copy_mode`/
`hemi_parent`; v3 folded the first Gemini review; v2 added cross-species coverage.)

**Author intent**: make chromosome inheritance **explicit and learnable**. The
schema should say, in plain terms, *"an offspring of sex S receives N copies from
each parent, and a parent of sex S does/doesn't recombine when making gametes."*
Eliminate the abbreviation `chr`, eliminate the sex dimension baked into column
names, and — critically — **stop overloading one `sex` column to mean two
different individuals**. Reserve `line_name` so crossbreeding-specific rules are
*rows, not a rewrite*, and confirm the schema expresses the diploid inheritance
systems of every livestock, plant, and insect species we intend to support.

**Goal wording (Codex review):** the aim is **"no future fundamental rewrite,"
not "never add a column."** Per CLAUDE.md ("reserve clean dimensions now *only
when* altering later would be painful"), adding a nullable column to a *long*
table later is deliberately cheap. Known future columns are documented
(Finding E), not pre-built.

---

## The core problem this redesign fixes: `sex` meant two different individuals

`chr_meta` is the **only** table in the schema that encodes a dimension in its
column *names* instead of its *rows* (`copy_mode_M`/`copy_mode_F`,
`recombines_M`/`recombines_F`). That alone violates the Schema Design Bias. But
the deeper, subtler defect — the one that makes sex chromosomes genuinely hard to
reason about — is that **"which sex" means two different things depending on which
column you read:**

- **Copy count** is a property of the **offspring** (a *male* is XY → his X is
  single-copy). `copy_mode_M` keys on the **carrier/offspring** sex.
- **Recombination** is a property of the **parent producing the gamete** (a
  *Drosophila male parent* is achiasmatic). `recombines_M` keys on the **parent**
  sex.

So in the current wide table, the `_M` in `copy_mode_M` and the `_M` in
`recombines_M` **refer to different individuals**. A naive long refactor that
collapses both into one `sex` column per row makes this *worse* — one `sex` value
would have to mean "offspring" for `copy_mode` and "parent" for `recombines`
simultaneously. That incoherence is exactly what makes building a sex chromosome
confusing.

**The fix is to separate the two questions into two tables, each with a single,
unambiguous meaning of sex, and to state copy origin explicitly instead of
encoding it in `copy_mode` + `hemi_parent`.**

**Pre-1.0 break** — the old `chr_meta` table, `copy_mode`/`hemi_parent`, and the
old function name are removed completely, no shims.

---

## Target schema — two explicit long tables

### `chr_inheritance` — copy count & origin, keyed by *offspring* sex

> "If an offspring is sex S, how many copies of this chromosome does it inherit
> from each parent?"

```text
chr_inheritance
  chr_name        VARCHAR   -- FK genome_meta.chr_name
  offspring_sex   VARCHAR   -- NULL = all/default; 'M' / 'F'
  line_name       VARCHAR   -- NULL = all lines (reserved, mirrors genome_map)
  from_parent_1   UTINYINT  -- copies inherited from parent_1 (sire / male-role)
  from_parent_2   UTINYINT  -- copies inherited from parent_2 (dam / female-role)
```

- **Logical key**: `(chr_name, offspring_sex, line_name)`, NULL-normalized in R.
- `from_parent_1` / `from_parent_2` map **directly** onto
  `ind_haplotype.parent_origin` (1 = sire, 2 = dam) and `strand`: a value of `k`
  means *k* haplotype copies in that parent-origin slot. Total copies =
  `from_parent_1 + from_parent_2`.
- **No `copy_mode`, no `hemi_parent`.** Origin is explicit in which count is
  non-zero; count is the number itself.

### `chr_recombination` — recombines, keyed by *parent* (gamete-producer) sex

> "When a parent of sex S produces gametes, does this chromosome recombine?"

```text
chr_recombination
  chr_name        VARCHAR   -- FK genome_meta.chr_name
  parent_sex      VARCHAR   -- NULL = both; 'M' / 'F' (the gamete-producing parent)
  line_name       VARCHAR   -- NULL = all lines (reserved)
  recombines      BOOLEAN
```

- **Logical key**: `(chr_name, parent_sex, line_name)`, NULL-normalized in R.
- **No surrogate PK** on either table (small, fully-managed config tables).
  *(Open question 1.)*
- **Reserved**: all columns of both tables.

**Row-local DB constraints (v5 review #8).** Both tables carry cheap DuckDB
`CHECK`/`NOT NULL` constraints at `CREATE TABLE` for row-local invariants, on top
of the R validators (which own the NULL-normalized logical key and fallback
resolution that a plain `UNIQUE` can't express under DuckDB NULL semantics):

```sql
CREATE TABLE chr_inheritance (
  chr_name       VARCHAR  NOT NULL,
  offspring_sex  VARCHAR  CHECK (offspring_sex IN ('M','F')),  -- NULL passes (default row)
  line_name      VARCHAR,
  from_parent_1  UTINYINT NOT NULL CHECK (from_parent_1 >= 0),
  from_parent_2  UTINYINT NOT NULL CHECK (from_parent_2 >= 0),
  CHECK (from_parent_1 + from_parent_2 <= 2)   -- diploid-release constraint; relaxed at ploidy > 2
);
CREATE TABLE chr_recombination (
  chr_name    VARCHAR NOT NULL,
  parent_sex  VARCHAR CHECK (parent_sex IN ('M','F')),          -- NULL passes (both parents)
  line_name   VARCHAR,
  recombines  BOOLEAN NOT NULL
);
```

The `chr_name → genome_meta.chr_name` **FK is enforced in R** (validator orphan
check), **not** as a DB constraint — matching tidybreed's convention (cf.
`validate_genome_map()`); a DB-level FK here would diverge from the rest of the
package (v5 review #8, refined).

**Why two tables, in one sentence:** "which sex" means the **child** for copying
and the **parent** for recombination — one table forces one `sex` column to mean
both (the thing that confused users); two tables let each column *say what it
means* (`offspring_sex` vs `parent_sex`).

**Bonus:** the two concerns now resolve **independently**, which structurally
eliminates the cross-concern shadowing bug (Codex #5 / Q10) — a sex-specific copy
rule and a line-specific recombination rule live in different tables and can never
out-prioritize each other.

### Reading the counts (why this is clearer than `copy_mode`/`hemi_parent`)

| Case | old (`copy_mode`,`hemi_parent`) | new (`from_parent_1`,`from_parent_2`) |
|---|---|---|
| Autosome | `full`, — | `1, 1` |
| Male's X (XY) | `half`, `parent_2` | `0, 1` |
| Male's Y | `half`, `parent_1` | `1, 0` |
| Female's Y (absent) | `none`, — | `0, 0` |
| Maternal mitochondria | `half`, `parent_2` | `0, 1` |
| Paternal chloroplast (conifer) | `half`, `parent_1` | `1, 0` |
| Biparental/heteroplasmic organelle | `full`, — | `1, 1`, `recombines=FALSE` |
| **Diploid uniparental** (both copies one parent) | **inexpressible** | `2, 0` — *storage only; kernel-future (v5 review #2)* |

`copy_mode` is *relative to ploidy*; the new counts are **absolute** and
correct at the current ploidy (`ind_meta.ploidy = 2`, enforced). When ploidy > 2
lands, a reserved `copy_count_basis` column distinguishes
`relative_to_ploidy` from `absolute` (Finding E) — documented, not built. Absolute
integer counts were chosen over the relative model because clarity for the task
users actually do today (sex chromosomes, organelles) outweighs auto-scaling for
unbuilt polyploid/haplodiploid biology, which needed special handling either way
(Finding C).

---

## Response to Gemini critique (pre-redesign round)

*(This section addresses the earlier review that prompted the v5 redesign; the
current v5 reviews are answered in "Response to the v5 review round" below.)*

"Strongly Approve with Enhancements." Disposition under the v5 two-table design:

| Point | Decision |
|---|---|
| **Biparental/heteroplasmic organelle** | **Accepted** — now `from_parent_1=1, from_parent_2=1, recombines=FALSE` (one copy from each parent, non-recombining). Even cleaner than the old `copy_mode=full` phrasing. Caveat unchanged: models haplotype presence, not copy-number proportions. |
| **Homeologous exchange** (Ph1-mutant / synthetic wheat) | **Accepted** — real gap; still a future `ploidy=6` + `homeology_group` + `inheritance_mode` item (Finding E). Additive nullable columns. |
| **Q2 — NULL-safe upsert** with `IS NOT DISTINCT FROM` | **Accepted** — applies to both tables' logical keys. Q2 decided. |
| **Q3 — pre-resolve maps, index `[[sex]][[chr]]`** | **Accepted, generalized** — now **two** pre-resolved maps (inheritance by offspring_sex, recombination by parent_sex), each keyed by `(sex, line_name)`. Q3 decided. |
| **Haplodiploidy** `own_ploidy %/% 2` gives `0` for drones | **Accepted as a note** — an `add_offspring()`/`ind_meta.ploidy` boundary, not a chromosome-table concern (Finding C). |
| Q1 no PK · Q5 hard error · Q6 accepted · Q7 follow-up commit | **Concur** — Q5 hard error now applies per table (resolve for both offspring sexes and both parent sexes). |

---

## Response to Codex critique (pre-redesign round)

*(This section addresses the earlier review that prompted the v5 redesign; the
current v5 reviews are answered in "Response to the v5 review round" below.)*

Codex approved the long-table direction and pushed for a much more explicit schema.
**v5 adopts the substance of its three architectural objections** (#2/#3/#5),
which is why the schema changed from v4's single table:

| Codex point | v5 disposition |
|---|---|
| **#1 `hemi_parent` validation bug** for `copy_mode='none'` | **Moot — `hemi_parent` removed entirely.** Origin is now explicit per parent (`from_parent_1`/`from_parent_2`); an absent chromosome is simply `0, 0`. The whole class of bug disappears. |
| **#2 (Q8) `copy_mode`+`hemi_parent` too coarse / not future-proof** | **Accepted.** Replaced by explicit `from_parent_1`/`from_parent_2`, which map directly to `ind_haplotype.parent_origin`/`strand` (exactly Codex's recommended direction) and express diploid uniparental (`2,0`) and future polyploid per-parent counts that the old model could not. |
| **#3 (Q9) one `sex` conflates carrier vs gamete-producer** | **Accepted.** Split into `chr_inheritance.offspring_sex` and `chr_recombination.parent_sex` — two tables, two column names, impossible to conflate. |
| **#5 (Q10) resolver shadowing** (sex seed shadows a line recombines override) | **Accepted (cross-concern case resolved structurally).** Copy and recombination are separate tables resolved independently, so a copy rule can never shadow a recombination rule. The remaining *within-table* sex-vs-line precedence is handled by a validator shadowing check (dormant until `line_name` is populated). This is Codex's own preferred "option 1 (split tables)." |
| **#4 `line_name` can't model crossbred heterozygous SVs** | **Accepted — claim narrowed.** `line_name` is a pure-line/founder default; per-copy inversion-heterozygote suppression depends on the *pair* of homolog structures (`ind_haplotype.line_origin`) and is a documented future extension, not built here (Q4). |
| **Organelle count tied to nuclear ploidy** | **Accepted as future item.** With absolute counts at ploidy 2 this is correct today; `copy_count_basis`/`genome_compartment` reserved for ploidy > 2 (Finding E). |
| **Allopolyploid framing / PAR / B-chromosomes / NULL-safe row tooling / `assert_qtl_autosomal` & `is_plain_autosome` on resolved rules / Smaller #1–#6** | **Accepted** — see Findings, Boundaries, and file-by-file. Row tooling uses `IS NOT DISTINCT FROM`; the autosome checks read *resolved* inheritance+recombination, not raw columns; resolvers/helpers live in one file. |
| **"never update later" goal wording** | **Accepted (reframed)** — "no fundamental rewrite," not "never add a nullable column." |

**Still declined — do not pre-build (CLAUDE.md).** Codex's "add all reserved
columns now" (`inheritance_mode`, `pairing_group`, `pairing_param`,
`centromere_pos_bp`, `structural_state_name`, `copy_count_basis`,
`genome_compartment`) conflicts with *"Do not implement future biology before it
is needed."* On **long** tables a later nullable column is cheap. We adopted the
three changes that improve *today's* clarity (`from_parent_1`, `from_parent_2`,
and the `offspring_sex`/`parent_sex` split); we **document** the rest (Finding E)
and add each when its feature lands. The one class we chose to build now is the
per-parent contribution (#2) because it *simplifies* the current model rather than
anticipating unbuilt biology.

---

## Response to the v5 review round (Codex + Gemini)

Both reviewers assessed the two-table v5 design. **Codex now approves the
architecture** ("the architectural correction I wanted from the first review");
**Gemini strongly approves**. Neither asks to change the schema shape — the two
tables, explicit per-parent counts, and `offspring_sex`/`parent_sex` split all
stand. What remains is a set of **implementation-precision refinements** the two
reviewers independently converged on. Disposition:

| # | Point (both reviewers agree unless noted) | Disposition |
|---|---|---|
| 1 | **`parent_sex` is not gamete role** — one hermaphrodite doing both pollen and ovule meiosis with *different* recombination can't be selected by `ind_meta.sex`. Both recommend **narrow scope (Option 1)**: keep `parent_sex` = producing parent's genetic sex `{M,F}`; role-specific recombination *within one* hermaphrodite → Boundaries/unsupported. | **Accepted (narrow scope).** Conservative, reversible (a `gamete_role` dimension is additive later), and CLAUDE.md-aligned (don't build unbuilt biology). Coverage claim for hermaphrodites narrowed; new Boundaries row. **Flagged for you** — see summary; the only genuine scope call in this round. |
| 2 | **`2,0` (uniparental disomy) is storage-expressible but not kernel-implemented** — one diploid parent yielding two child copies needs an unspecified mechanism (duplicate a homolog? transmit both? bypass meiosis?). | **Accepted.** `2,0`/`0,2` are allowed in storage but `add_offspring()` throws an explicit "uniparental/non-standard diploid transmission not implemented" error at the kernel boundary. Coverage marks it schema-expressible / kernel-future. My v5 wording implied it "works" — corrected. |
| 3 | **Ploidy validation has no single DB value** — `from_parent_1+from_parent_2 ≤ ind_meta.ploidy` can't be checked against a per-individual, possibly-empty column at config time. | **Accepted — real spec bug.** For this diploid release the invariant is the constant `≤ 2` (R validator + DB `CHECK`), plus the existing `assert_ploidy_2()` consumer guard, plus a **kernel-boundary assertion that resolved total == emitted `ind_haplotype` rows**. When ploidy > 2 lands, validation moves to the actual offspring/parents being processed. |
| 4 | **Resolver line context** — whose `line_name` selects each table wasn't pinned. | **Accepted.** **Inheritance resolves by the *offspring's* line** (its expected karyotype); **recombination resolves by the *producing parent's* line** (meiosis is a property of the parent's germline). Pinned in `test-resolver-context.R`. Passing the offspring line to recombination would silently defeat pure-line recombination defaults in crossbreds. |
| 5 | **Shadowing validator must be a concrete algorithm**, not "a row that would be out-prioritized." | **Accepted.** Deterministic rule per table: enumerate every configured non-NULL `line` × supported sex `{M,F}`; if `(sex,NULL)` and `(NULL,line)` both match and **differ** and no explicit `(sex,line)` row exists → error naming chr, sex, line. |
| 6 | **`overwrite = FALSE` undefined.** | **Accepted.** `FALSE` = error if the exact NULL-safe logical key already exists, else insert. Tested across all four NULL/non-NULL `(sex, line)` combinations. |
| 7 | **Validation must be inside the write transaction.** | **Accepted.** `dbBegin` → delete-then-insert → run validator → `dbCommit`; any failure → `dbRollback`, preserving the prior valid row. |
| 8 | **Add cheap DB-level constraints** (NOT NULL, `sex IN ('M','F')` [NULL passes], counts ≥ 0, `sum ≤ 2`, `recombines` NOT NULL). | **Accepted, with one refinement.** Add the row-local `CHECK`/`NOT NULL` constraints (Gemini's DDL). **Declined the DB-level foreign key** on `chr_name → genome_meta`: tidybreed enforces FKs *in R validators* (e.g. `validate_genome_map()`'s orphan check), not as DB constraints — a real FK here would diverge from the package convention. The orphan check stays in `validate_chr_*()`. The `sum ≤ 2` CHECK is a **diploid-release constraint**, dropped/relaxed when ploidy > 2 lands. |
| — | **Coverage wording** (Codex's 4 edits) — soften "express the storage for every system," qualify hermaphrodite recombination, mark `2,0`/autopolyploid kernel-future, keep organelle/PAR/B-chr/homeologous/SV caveats. | **Accepted** — applied in the coverage verdict and Boundaries. |
| — | **New tests** — cross where offspring vs parent sex+line all differ; shadowing require/allow; rollback preserves old row; `overwrite=FALSE`; resolved totals == emitted haplotype rows; `2,0` fails at kernel; gamete-role rejection. | **Accepted** — added to the test table. |
| — | Non-blocking: no surrogate PK fine; cache by actual distinct contexts incl. parent line; rename helper file in the Q7 sweep; one-concern-per-call worth the verbosity. | **Concur** — no change beyond wording. |

**Net:** nothing here reopens the architecture. Every point is an accepted
precision fix except **8's FK** (refined to match the package's R-enforced-FK
convention) and **1's scope** (accepted narrow, but surfaced for your call). These
are folded into the schema, validator, resolver, API, Boundaries, and test
sections below.

---

## The key design property (why sex chromosomes "just work")

**Each sex chromosome is its own `chr_name` with per-offspring-sex copy counts.**
Because Z and W (or X and Y) are *separate chromosomes*, the presence/absence of a
sex chromosome in an individual *falls out of* `from_parent_1`/`from_parent_2` for
that offspring sex. Combined with tidybreed's **sex-first** assignment (offspring
sex is drawn by ratio, then chromosomes follow), inheritance is automatically
consistent with no mechanistic sex-determination model:

- A ZW **son** (M): Z `from_parent_1=1, from_parent_2=1` (one from each parent) →
  **ZZ**; W `0, 0` → no W. Correct.
- A ZW **daughter** (F): Z `from_parent_1=1, from_parent_2=0` → Z from sire only;
  W `from_parent_1=0, from_parent_2=1` → W from dam. Correct **ZW**.

No row needs to "know" the system is ZW vs XY — the counts encode it. The one
thing sex-first assignment cannot model is sex *determined by* genotype (csd in
bees, meiotic drive, polygenic SD); that is an `add_offspring()` sex-assignment
concern, not the chromosome tables (see Boundaries).

---

## Cross-species coverage review

Legend: **✓ schema fits** = the two tables + `define_chromosome()` can store the
rule; the *mechanism* may be future work (noted). Representations give
`(from_parent_1, from_parent_2)` per offspring sex, and recombination per parent
sex where non-default. Recall **parent_1 = sire (male-role), parent_2 = dam
(female-role)**.

### Sex-determination systems

| System | Example species | Representation | Status |
|---|---|---|---|
| **XY, ♂ heterogametic** | cattle, pig, sheep, goat, horse, dog, mouse, human; Cannabis, asparagus, papaya | X: M `0,1`, F default `1,1`; Y: M `1,0`, F `0,0` (recombines FALSE) | ✓ |
| **ZW, ♀ heterogametic** | chicken, turkey, duck, quail; silkworm; some fish/reptiles; strawberry, poplar | Z: M default `1,1`, F `1,0`; W: M `0,0`, F `0,1` (recombines FALSE) | ✓ heterogamety = which sex carries the single copy |
| **X0 / XO** | grasshoppers, crickets, locusts; C. elegans (♂ X0, hermaphrodite XX); some Hemiptera | X: M `0,1`, F `1,1`; **no Y defined** | ✓ absent chromosome = simply not defined |
| **Z0** | some Lepidoptera moths | Z: M `1,1`, F `1,0`; no W | ✓ |
| **Haplodiploidy** | honeybee, ants, parasitoid wasps (Nasonia), thrips, some mites | No sex chromosomes; ♂ haploid via `ind_meta.ploidy=1`, ♀ diploid; autosomes default `1,1` for ♀, single-parent copy for ♂ drones | ✓ *schema fits with zero sex-chr rows*; per-sex **ploidy** + male-offspring copy default is future `ind_meta`/`add_offspring` work (Finding C) — nothing structural in these tables |
| **Multiple / neo-sex chromosomes** | platypus (X1Y1…X5Y5), some rodents/insects (X1X2Y) | each sex chromosome = its own `chr_name` with per-offspring-sex counts | ✓ storage fits; **chain segregation** is an offspring-kernel concern (Boundaries) |
| **TSD / ESD** | many reptiles, some fish | no sex-chr rows; sex by ratio/rule | ✓ sex-first; environment→sex is `add_offspring` policy |
| **Polygenic / homomorphic SD** | zebrafish, tilapia, Aedes (M-locus) | no heteromorphic rows; autosomal + sex by ratio, or a single SD locus | ✓ |
| **Hermaphrodite / monoecious** | Arabidopsis, rice, wheat, tomato, soybean; maize | all autosomes default `1,1`; no sex chromosomes; an individual acts as both sire & dam (mating logic) | ✓ inheritance trivial. **Role-specific recombination** (different pollen vs ovule maps *within one* individual) is **unsupported** — `parent_sex` keys the producer's genetic sex, not gamete role (v5 review #1; Boundaries) |

### Recombination systems

| Feature | Example | Representation | Status |
|---|---|---|---|
| **♂ achiasmy** (genome-wide) | *Drosophila* ♂; many Diptera (poss. BSF) | `recombines_M = FALSE` at `define_genome()` → seeds `parent_sex='M', recombines=FALSE` genome-wide | ✓ one line (Finding B) |
| **♀ achiasmy** (genome-wide) | silkworm *Bombyx* ♀; other Lepidoptera | `recombines_F = FALSE` at `define_genome()` | ✓ (Finding B) |
| **Non-recombining single chromosome** | *Drosophila* dot chr4 (both sexes) | `chr_recombination`: `parent_sex=NULL, recombines=FALSE` — one row | ✓ |
| **Non-recombining sex-chr body** | Y, W | that chr: `recombines=FALSE` | ✓ |
| **Heterochiasmy (different *rates*)** | human ♀ ≈1.6× ♂ | **not these tables** — per-sex `pos_cM` in `genome_map` | ✓ correct home is `genome_map` (Boundaries) |
| **Pseudoautosomal region (PAR)** | mammalian X-Y PAR recombines, rest doesn't | sub-chromosomal — not a whole-chr boolean | ⚠ boundary — `genome_map` / PAR-as-tiny-chromosome is an *approximation* (Boundaries) |

### Organellar & uniparental genomes

| Organelle / mode | Example | Representation | Status |
|---|---|---|---|
| **Mitochondria, maternal** | ~all animals; most plants | `from_parent_1=0, from_parent_2=1, recombines=FALSE` (one row, both sexes) | ✓ |
| **Chloroplast, maternal** | most angiosperms | `0, 1, recombines=FALSE` | ✓ |
| **Chloroplast, paternal** | conifers / Pinaceae | CP: `1, 0`; MT: `0, 1` (separate chromosomes → separate rows) | ✓ paternal organelles are just `from_parent_1=1` |
| **Biparental / heteroplasmic** | *Pelargonium*, alfalfa; bivalve DUI | `1, 1, recombines=FALSE` — one copy from each parent | ✓ (Gemini). Caveat: models haplotype *presence*, not copy-number/proportion (out of scope) |

### Ploidy systems

| System | Example | Where it lives | Status |
|---|---|---|---|
| **Diploid** | most livestock, Arabidopsis, mouse | counts `1,1` etc. at `ind_meta.ploidy=2` | ✓ |
| **Haploid sex (haplodiploidy)** | honeybee/ant/wasp ♂ | per-sex `ind_meta.ploidy` + male-offspring single-parent copy default | ✓ schema fits; **ploidy work is future** (Finding C) |
| **Autopolyploidy** | potato (4N), autotetraploid alfalfa, banana (3N) | `ind_meta.ploidy=4`; counts `2,2`; **pairing mode** = future column (Finding E) | ✓ counts fit once `copy_count_basis`/absolute counts generalize; pairing/segregation is the future extension |
| **Allopolyploidy (subgenomes)** | bread wheat (AABBDD), canola, cotton | diploid-with-more-chromosomes under disomic inheritance — 2 copies of each chr (`1,1`); subgenome tagged by `line_name`/user col | ✓ *approximation*, no ploidy machinery (Finding E) |

### Model-species spot check

Arabidopsis (5 autosomes, selfing) ✓ · Mouse (XY) ✓ · *Drosophila* (XY, ♂
achiasmy, non-recombining chr4) ✓ · C. elegans (XX/X0) ✓ · Chicken (ZW) ✓ ·
Maize (monoecious) ✓ · Silkworm (ZW, ♀ achiasmy) ✓ · Honeybee (haplodiploid) ✓
(via ploidy) · Nasonia ✓ · Zebrafish ✓ · Tribolium (XY) ✓.

**Verdict (scoped per v5 review): the two tables express the current diploid
copy-count and whole-chromosome recombination configuration for the supported
subset above.** No additional column is required for diploid support. Two honest
limits, both *storage-expressible but kernel-future*: uniparental disomy (`2,0`)
and autopolyploid counts (`2,2`) store fine but the meiosis kernel does not yet
implement them (it errors at the boundary — v5 review #2); role-specific
recombination within one hermaphrodite is unsupported (v5 review #1). Finding B
(genome-wide per-parent-sex recombination default) is accepted; Finding E's
`copy_count_basis`/`inheritance_mode`/etc. remain cheap additive future columns.

---

## Findings from the review

### Finding A — explicit per-parent origin subsumes `hemi_parent` (and the paternal-organelle win)

The wide table's single shared `hemi_parent` could not cleanly express "MT from
dam **and** CP from sire" in one genome. v5 removes `hemi_parent` altogether:
origin is explicit per parent, so conifer-style paternal-chloroplast /
maternal-mitochondria is two ordinary rows (`CP: 1,0` and `MT: 0,1`). This is the
same free win as before, now with fewer concepts.

### Finding B — genome-wide per-parent-sex achiasmy needs a genome-wide default  **[ACCEPTED]**

Silkworm females and *Drosophila* males are achiasmatic across the **entire
genome**. Without a genome-wide default you would write one recombination row per
chromosome for the achiasmatic sex — a `define_chromosome(chr, parent_sex = "F",
recombines = FALSE)` call for each of the silkworm's ~28 chromosomes — which is
ergonomically bad. **`define_genome()` gains `recombines_M` / `recombines_F`
(both default `TRUE`)** — genome-wide recombination defaults for **male-parent**
and **female-parent** meiosis. The seed step materializes them into
`chr_recombination`:

- if `recombines_M == recombines_F` → one `parent_sex=NULL` row per chromosome;
- if they differ → a `parent_sex='M'` and a `parent_sex='F'` row per chromosome.

Bombyx/Drosophila become a one-liner at genome setup; `define_chromosome()` is
reserved for genuinely per-chromosome overrides (Y/W/chr4). The `_M`/`_F` on
`define_genome()` here unambiguously mean **parent sex** (recombination is a
parental-meiosis property), consistent with `chr_recombination.parent_sex`.

### Finding C — haplodiploidy & polyploidy are `ind_meta.ploidy`/`add_offspring`, not these tables

Honeybee support needs **zero** sex-chromosome rows — drones are haploid by
whole-genome ploidy. The chromosome tables do not block it. Concrete gotcha
(Gemini): `add_offspring()`'s gamete-ploidy computation (`own_ploidy %/% 2L`)
gives `1 %/% 2 = 0` for a haploid drone, so he would contribute nothing — but a
drone passes his entire single genome to every daughter. Under **absolute**
per-parent counts this is even more clearly an `add_offspring`/`ind_meta.ploidy`
concern: a drone-son's autosome is `from_parent_1=0, from_parent_2=1` (all from
the queen), best expressed via a genome-wide male-offspring default when
per-individual ploidy lands — not a per-chromosome rule. Action: **none in the
chromosome tables**; recorded as the first `add_offspring` fix for haplodiploidy.

### Finding D — seed-row lifecycle & the `NULL` fallback (per table)

Each table is seeded with a `sex=NULL`/`line=NULL` default row per chromosome
(inheritance `1,1`; recombination `recombines_*`). The resolver's precedence means
sex-specific overrides compose cleanly, and you often override **only the
deviating sex**:

- **X** (mammal): the inheritance seed `(X, NULL, 1,1)` already gives correct
  *female* behavior. Override **only** the male:
  `define_chromosome("X", offspring_sex="M", from_parent_1=0, from_parent_2=1)`.
  → **one override call.**
- **Y**: the `1,1` seed is wrong for both sexes, so both offspring-sex inheritance
  rows are overridden (`M: 1,0`; `F: 0,0`), plus one recombination row
  (`parent_sex=NULL, recombines=FALSE`).

**Each call is a plain per-row upsert of one logical key** — no "reset the whole
chromosome" magic. `offspring_sex='M'` upserts the `(chr, 'M', line)` inheritance
row and touches nothing else; `offspring_sex=NULL` upserts the `(chr, NULL, line)`
default row and leaves any `'M'`/`'F'` rows in place (the resolver's precedence
still lets those override the default, exactly as intended). This is more
predictable than a reset and matches how `genome_map`'s writer treats each row as
independent. To change a sex-specific rule you overwrite that row; there is no
implicit clearing of sibling rows.

`validate_chr_inheritance()`/`validate_chr_recombination()` guarantee every
chromosome resolves for both `M` and `F` in each table (the seeds make this
automatic).

### Finding E — polyploidy (4N+): scope of the absolute-count model

1. **Allopolyploids are diploids-with-more-chromosomes.** Bread wheat (AABBDD)
   carries exactly 2 copies of each chromosome under disomic inheritance —
   mechanically a diploid with 21 chromosome names, subgenome-tagged by
   `line_name`/user column. No ploidy machinery. This is an **approximation, not
   full allopolyploid coverage** (Codex): it discards that the organism is 6N and
   cannot model homeologous exchange without reserved columns.
2. **Only autopolyploidy genuinely needs ploidy > 2.** Potato (4N), banana (3N):
   truly homologous copies. `ind_meta.ploidy=4` + per-parent counts `2,2` is
   correct once counts generalize past ploidy 2. A mixed genome (some chr 4N,
   others 2N) breaks the single scalar `ploidy` — exotic, out of scope.
3. **The missing rule is pairing/segregation mode — future, additive.** Absolute
   counts give *how many* homologs but not *how they pair/assort* (disomic vs
   polysomic vs preferential, double reduction). This is an inheritance rule that
   belongs beside these tables as a future nullable `inheritance_mode`
   (and perhaps `pairing_param`) column — cheap under long format, not built now.
4. **`copy_count_basis` for counts.** Today counts are **absolute** and correct at
   ploidy 2. When ploidy > 2 lands, a `copy_count_basis`
   (`relative_to_ploidy` vs `absolute`) column lets autosomes scale with ploidy
   while organelles stay absolute (a maternal mito must remain 1 copy regardless
   of nuclear ploidy — Codex). Reserved, not built.
5. **`homeology_group`** (Gemini): the disomic-as-diploid trick cannot model
   homeologous exchange (1A↔1B↔1D in *Ph1*-mutant/synthetic wheat) — a
   21-separate-diploid layout never lets homeologs pair. Needs `ploidy=6` +
   `homeology_group` + `inheritance_mode`. Additive nullable columns; not needed
   for standard (Ph1-intact) allopolyploids.
6. **`genome_compartment`** (`nuclear`/`mitochondrial`/`plastid`) — reserved to
   keep organelle handling from being tied forever to nuclear rules.
7. **Genuinely not these tables (kernel / other tables):** multivalent formation &
   double reduction (needs centromere → `genome_meta`/`genome_map`); unreduced (2n)
   gametes and polyploid *formation* (→ `add_offspring()` gamete/ploidy);
   triploid/odd-ploidy meiosis; aneuploidy.

**Net for 4N+:** the schema is sound; the allopolyploid clarification *reduces*
what needs new machinery; the future additions (`inheritance_mode`,
`pairing_param`, `homeology_group`, `copy_count_basis`, `genome_compartment`) are
nullable additive columns, documented not pre-built. Nothing in this diploid-only
refactor blocks any of it.

---

## Boundaries — what the chromosome tables deliberately do **not** model

| Concern | Correct home |
|---|---|
| **Recombination *rate* differences by sex (heterochiasmy)** | `genome_map` per-sex `pos_cM` maps |
| **Pseudoautosomal region (sub-chromosomal recombination on/off)** | interval-level metadata (`genome_map`) or a compound-chromosome kernel. **PAR-as-tiny-chromosome is an *approximation*** — it breaks physical linkage to the X/Y body (Codex) |
| **Per-individual variable chromosome presence** (B chromosomes, aneuploidy) | individual-level mechanism (stochastic transmission rule or an `ind`-level chromosome-copy table) — **not** a per-`(sex,line)` rule (Codex) |
| **Crossbred heterozygous structural variants** (inversion het suppression) | offspring kernel inspecting both homolog states (per-copy `ind_haplotype.line_origin`) + a future per-copy structural-state concept — **not** a single `line_name` (Codex #4) |
| **Role-specific recombination within one individual** (different pollen vs ovule maps in a hermaphrodite/monoecious plant) | future `gamete_role` dimension. `chr_recombination.parent_sex` keys the producer's genetic sex, not gamete role — role-decoupled selfing recombination is **unsupported** (v5 review #1) |
| **Uniparental disomy / non-standard diploid transmission** (`2,0`, `0,2`) | offspring meiosis kernel. **Storage-expressible but kernel-future** — `add_offspring()` errors explicitly rather than guessing how one diploid parent yields two child copies (v5 review #2) |
| **Sex *determined by* genotype** (csd, meiotic drive, polygenic SD) | `add_offspring()` sex-assignment policy (tidybreed is sex-first) |
| **Chain segregation of multi-sex-chromosome systems** (platypus) | offspring meiosis kernel |
| **Polysomic / multivalent pairing — the *mechanism*** | offspring meiosis kernel. The *config* ("polysomic, double-reduction rate α") is a future `inheritance_mode` column (Finding E) |
| **Double reduction / centromere position** | `genome_meta`/`genome_map` + offspring kernel |
| **Unreduced (2n) gametes / polyploid formation** | `add_offspring()` gamete + ploidy computation |
| **Per-generation / life-cycle mating-mode switches** (aphid cyclical parthenogenesis, apomixis, selfing vs outcrossing) | mating / `add_offspring()` logic |
| **Whole-genome ploidy** (haplodiploidy, polyploidy) | `ind_meta.ploidy` |

Rule of thumb: **`chr_inheritance` answers "for an offspring of sex S (line L),
how many copies of this chromosome come from each parent?"; `chr_recombination`
answers "when a parent of sex S makes gametes, does this chromosome recombine?" —
nothing about *rates*, *sub-chromosomal structure*, *whole-genome ploidy*, or
*how sex itself is decided*.**

---

## Naming change: `define_chr()` → `define_chromosome()`

Per tidybreed's explicit-naming ethos, `chr` is an unacceptable abbreviation for a
core entry point. Rename the function and its file now:

- `R/define_chr.R` → `R/define_chromosome.R`
- `define_chr()` → `define_chromosome()`
- `man/define_chr.Rd` → `man/define_chromosome.Rd`
- `NAMESPACE`, all `@seealso`/examples/vignettes/tests, CLAUDE.md.

The **table/column** rename (`chr_inheritance`/`chr_recombination` →
`chromosome_inheritance`/`chromosome_recombination`, `chr_name` →
`chromosome_name`, `genome_meta.chr` → `chromosome_number`) stays a **separate
follow-up commit** (Q7) to keep diffs reviewable — the two new tables are created
with the current `chr_*` convention now and swept with the rest of the genome
schema then.

---

## Resolvers + validators (mirror `genome_map_helpers.R`)

Added to `R/chr_meta_helpers.R` (one home for resolvers + helpers, Codex Smaller #1;
file keeps its name until the Q7 sweep).

**Whose sex, whose line? (v5 review #4)** The two resolvers are called in
different contexts, and passing the wrong entity's `line_name` silently corrupts
crossbred offspring:

- **`resolve_chr_inheritance`** is called with the **offspring's** sex **and the
  offspring's line** — the karyotype an individual of that line/sex should carry.
- **`resolve_chr_recombination`** is called with the **producing parent's** sex
  **and the producing parent's line** — meiotic recombination is a property of the
  parent's germline, not the offspring's. In an A×B cross the sire's gamete
  recombines under Line A's rule, the dam's under Line B's.

`test-resolver-context.R` pins this in a cross where offspring and both parents
differ in sex *and* line.

### `resolve_chr_inheritance(conn, offspring_sex, line_name = NULL)`

One row per `chr_name`: `(from_parent_1, from_parent_2)`. Precedence
`(offspring_sex=S,line=L) → (offspring_sex=S,NULL) → (NULL,line=L) → (NULL,NULL)`
via the same priority-window SQL as `resolve_genome_map()`:

```sql
WITH cand AS (
  SELECT chr_name, from_parent_1, from_parent_2,
         (CASE WHEN offspring_sex IS NOT NULL THEN 2 ELSE 0 END) +
         (CASE WHEN line_name     IS NOT NULL THEN 1 ELSE 0 END) AS prio
  FROM chr_inheritance
  WHERE (offspring_sex IS NULL OR offspring_sex = $sex)
    AND (line_name     IS NULL OR line_name     = $line)
),
ranked AS (
  SELECT *, ROW_NUMBER() OVER (PARTITION BY chr_name ORDER BY prio DESC) AS rn
  FROM cand
)
SELECT chr_name, from_parent_1, from_parent_2 FROM ranked WHERE rn = 1
```

Errors if any `genome_meta.chr_name` has no resolved row.

### `resolve_chr_recombination(conn, parent_sex, line_name = NULL)`

Identical shape on `chr_recombination`, keyed by `parent_sex`, returning
`recombines`. Errors on any unresolved `chr_name`.

### `validate_chr_inheritance(conn)`

Run after every write (inside the write transaction — see API). Checks:
(1) logical-key uniqueness `(chr_name, offspring_sex, line_name)` NULL-normalized;
(2) no orphan `chr_name` (the R-enforced FK to `genome_meta`, v5 review #8);
(3) valid `offspring_sex` ∈ {NULL,M,F}; (4) `from_parent_1`, `from_parent_2` are
non-negative integers with `from_parent_1 + from_parent_2 ≤ 2` — the **constant**
diploid-release invariant, not a per-individual `ind_meta.ploidy` lookup (v5
review #3); (5) every chromosome resolves for **both** `offspring_sex` `'M'` and
`'F'` (hard error, Q5); (6) **deterministic shadowing check** (v5 review #5):

> For each configured non-NULL `line` value and each supported sex `S ∈ {M,F}`:
> if a `(sex=S, line=NULL)` row **and** a `(sex=NULL, line=L)` row both exist and
> their `(from_parent_1, from_parent_2)` **differ**, and no explicit
> `(sex=S, line=L)` row is present → **error** naming the chromosome, `S`, and `L`,
> instructing the user to materialize the explicit `(S, L)` row.

This is dormant until `line_name` is populated (no non-NULL lines → no check).

### `validate_chr_recombination(conn)`

Same shape on `chr_recombination`: (1) key uniqueness `(chr_name, parent_sex,
line_name)`; (2) no orphan `chr_name`; (3) valid `parent_sex` ∈ {NULL,M,F};
(4) `recombines` non-null boolean; (5) resolves for both `parent_sex` `'M'` and
`'F'` (hard error); (6) the same deterministic shadowing check on `parent_sex`
vs `line_name`, comparing the `recombines` value.

**Cross-concern shadowing (Codex #5) cannot occur** — copy and recombination are
separate tables resolved by separate calls, so a copy rule never shadows a
recombination rule. Only the within-table sex-vs-line case remains, covered by
each validator's check (6).

---

## `define_chromosome()` API

```r
define_chromosome(pop, chr_name,
                  offspring_sex = NULL,  # inheritance row key (child sex); NULL = all
                  parent_sex    = NULL,  # recombination row key (parent sex); NULL = both
                  line_name     = NULL,  # reserved; NULL = all lines
                  from_parent_1 = NULL,  # copies from parent_1 (sire / male-role)
                  from_parent_2 = NULL,  # copies from parent_2 (dam / female-role)
                  recombines    = NULL,  # per parent_sex
                  overwrite     = TRUE)
```

**Routing — one function, but one concern per call:**

- Supplying `from_parent_1` + `from_parent_2` writes/updates a
  **`chr_inheritance`** row keyed by `(chr_name, offspring_sex, line_name)`. This
  is an *inheritance* call — `parent_sex` and `recombines` must **not** be given.
- Supplying `recombines` writes/updates a **`chr_recombination`** row keyed by
  `(chr_name, parent_sex, line_name)`. This is a *recombination* call —
  `offspring_sex`, `from_parent_1`, and `from_parent_2` must **not** be given.
- **A single call sets exactly one concern.** This is deliberate: mixing them in
  one call would pair an `offspring_sex` with a `recombines` value, and since
  `recombines` is keyed by *parent* sex, a reader could not tell whether that
  sex meant the child or the parent — the exact ambiguity this whole redesign
  removes. Setting both concerns is two clear calls, not one confusing one.

**Validation (errors, not silent coercion):**

- Exactly one concern per call: either both `from_parent_*` **or** `recombines`
  is supplied — never both, never neither.
- `from_parent_1`/`from_parent_2` must be supplied **together**; each is a
  non-negative integer and their sum is `≤ 2` — the **constant** diploid-release
  invariant, **not** a lookup against `ind_meta.ploidy` (which is per-individual
  and may be empty at config time — v5 review #3).
- `offspring_sex` (inheritance) / `parent_sex` (recombination) ∈ {`NULL`,`'M'`,`'F'`};
  the *other* concern's sex argument must be left `NULL`.
- Each call is a **plain per-row upsert** of its logical key (Finding D) — no
  clearing of sibling sex rows.
- `overwrite = FALSE` (v5 review #6): **error** if the exact NULL-safe logical key
  already exists; otherwise insert. (`overwrite = TRUE` upserts.) **Note:**
  `define_genome()` pre-seeds a `(chr, NULL, NULL)` row for every chromosome, so an
  `overwrite = FALSE` call whose key resolves to that default (`offspring_sex`/
  `parent_sex` and `line_name` all `NULL`) collides with the seed **by design** —
  `overwrite = FALSE` is for guarding *new* sex/line-specific rows, not for setting
  a chromosome's default. The `overwrite = FALSE` test must seed first, then assert
  the collision, across all four NULL/non-NULL `(sex, line)` combinations.

**Transaction & rollback (v5 review #7).** The whole write is one transaction:
`dbBegin` → delete-then-insert (RNG-neutral `duckdb_register`, key matched with
`IS NOT DISTINCT FROM`, Q2) → run the matching validator → `dbCommit`. Any write
**or** validator failure triggers `dbRollback`, so a rejected replacement never
leaves the config table in a corrupt state.

### Worked examples (proof the API covers the review)

One concern per call: `from_parent_*` lines set inheritance (keyed by
`offspring_sex`); `recombines` lines set recombination (keyed by `parent_sex`).

```r
# Mammal (mouse/cattle) — inheritance: override only the deviating sexes
pop <- pop |>
  define_chromosome("X", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 1) |>
  define_chromosome("Y", offspring_sex = "M", from_parent_1 = 1, from_parent_2 = 0) |>
  define_chromosome("Y", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 0) |>
  # recombination: Y non-recombining (parent_sex = NULL → both parent sexes)
  define_chromosome("Y", recombines = FALSE) |>
  # mitochondria: maternal (inheritance) + non-recombining (recombination)
  define_chromosome("MT", from_parent_1 = 0, from_parent_2 = 1) |>
  define_chromosome("MT", recombines = FALSE)

# Bird / silkworm (ZW)
pop <- pop |>
  define_chromosome("Z", offspring_sex = "F", from_parent_1 = 1, from_parent_2 = 0) |>
  define_chromosome("W", offspring_sex = "M", from_parent_1 = 0, from_parent_2 = 0) |>
  define_chromosome("W", offspring_sex = "F", from_parent_1 = 0, from_parent_2 = 1) |>
  define_chromosome("W", recombines = FALSE)

# Silkworm genome-wide female achiasmy — one line at genome setup (Finding B)
pop <- open_pop(...) |> define_genome(..., recombines_F = FALSE)

# Conifer organelles — paternal chloroplast, maternal mitochondria (both non-recombining)
pop <- pop |>
  define_chromosome("CP", from_parent_1 = 1, from_parent_2 = 0) |>
  define_chromosome("CP", recombines = FALSE) |>
  define_chromosome("MT", from_parent_1 = 0, from_parent_2 = 1) |>
  define_chromosome("MT", recombines = FALSE)

# Honeybee — NO chromosome-table sex rows; haplodiploidy lives in ind_meta.ploidy (future)
# Arabidopsis / maize — NO sex-chromosome rows; autosomes use the seed default (1,1)
```

---

## File-by-file changes

| File | Change |
|---|---|
| `R/define_genome.R` (`~L251-278`) | Replace the single `CREATE TABLE chr_meta` with `CREATE TABLE chr_inheritance` + `CREATE TABLE chr_recombination` (typed). Seed one `offspring_sex=NULL` inheritance row per chr (`from_parent_1=1, from_parent_2=1`). Add genome-wide `recombines_M`/`recombines_F` params (Finding B); seed `chr_recombination` (one `parent_sex=NULL` row when equal, else per-parent-sex rows). Run both validators after seeding. |
| `R/define_chr.R` → `R/define_chromosome.R` | Rename file + function; rewrite to the per-parent `from_parent_1`/`from_parent_2` + `offspring_sex`/`parent_sex` API. **One concern per call** (inheritance *or* recombination — validate and error on mixed/empty calls); plain per-row upsert (Finding D). |
| `R/chr_meta_helpers.R` | **One home** (Codex Smaller #1). Remove `resolve_chr_copy_count()` (counts are now absolute in `chr_inheritance`). Add `resolve_chr_inheritance()`, `resolve_chr_recombination()`, `validate_chr_inheritance()`, `validate_chr_recombination()`. `get_chr_meta_map()` → returns **two** pre-resolved maps: inheritance by `(offspring_sex[,line])` and recombination by `(parent_sex[,line])`. `is_plain_autosome()` operates on **resolved** rules (`from_parent_1==1 && from_parent_2==1` for both offspring sexes **and** `recombines` for both parent sexes). `assert_qtl_autosomal()` checks **resolved** inheritance (full diploid) for the relevant carrier sexes — not a raw column scan (Codex Smaller #4/#5). |
| `R/add_founders.R` (`~L268-306`) | Replace `keep_po_for()`'s `copy_mode`/`hemi_parent` switch with resolved inheritance: parent-origin slot 1 gets `from_parent_1` strands, slot 2 gets `from_parent_2` — a founder's single-copy sex chromosome sits in the non-zero slot automatically. |
| `R/add_offspring.R` (`~L275,397,682-711`) | Copy counts from `resolve_chr_inheritance(offspring_sex=off_sex, line_name=offspring_line)` (`from_parent_1` from the sire's gamete, `from_parent_2` from the dam's). Recombination from `resolve_chr_recombination(parent_sex=parent_sex[[pid]], line_name=parent_line[[pid]])` — **producing parent's sex *and* line** (v5 review #4). **Kernel-boundary guards (v5 review #2/#3):** if any producing diploid parent's resolved contribution is not `0` or `1` (e.g. `2,0`), throw an explicit "uniparental/non-standard diploid transmission not implemented" error; assert the resolved total copy count equals the number of `ind_haplotype` rows actually written. |
| `R/schema.R` (`~L128-142`) | Replace the `chr_meta` `.sm_col`/`.sm_tbl` entries with `chr_inheritance` (`chr_name`,`offspring_sex`,`line_name`,`from_parent_1`,`from_parent_2`) and `chr_recombination` (`chr_name`,`parent_sex`,`line_name`,`recombines`). |
| `R/sql_utils.R` (`L84,117,153`) | Remove `chr_meta`; add `TABLE_RESERVED_COLS`/`TABLE_PRIMARY_KEYS`/`TABLE_ROW_KEYS` entries for both new tables' logical keys. **NULL-safe logical-key predicate helper**; route all row matching/deletion on both tables through `IS NOT DISTINCT FROM` (Codex Medium). |
| `R/open_pop.R` (`~L245`) | Roxygen comment lists `chr_meta` among the tables `define_genome()` creates → replace with `chr_inheritance`/`chr_recombination`. (No functional change; `.create_core_tables()` does not touch these.) |
| `R/ploidy_helpers.R` (`~L4-8,34`) | `assert_ploidy_2()` docstring + error message name `chr_meta (copy_mode_M/copy_mode_F)` → reword to `chr_inheritance`/`from_parent_*`. Function logic unchanged (guards organism ploidy, not chr shape). |
| `NAMESPACE` | `define_chr` → `define_chromosome`. |
| `CLAUDE.md` | Replace the `### chr_meta` section with `### chr_inheritance` + `### chr_recombination`; rewrite the `define_chromosome()` section; add the Boundaries table. |
| `NEWS.md` / `DESCRIPTION` | Version bump + breaking-change entry. |

**Access-pattern decision (Q3):** `get_chr_meta_map()` runs each priority-window
query once per distinct context in the batch — **inheritance by distinct
`(offspring_sex, offspring_line)`, recombination by distinct `(parent_sex,
parent_line)`** (v5 review #4) — and returns two nested lookups, so consumers index
`inheritance[[off_sex]][[chr]]` and `recombination[[parent_sex]][[chr]]` with
**zero** SQL in the per-individual loop. For the common no-line-specific case that
is four queries total (M/F × two tables), independent of individual count; it grows
only with the number of lines that actually carry line-specific rows.

---

## Old-schema cleanup — delete completely, no shims

Per CLAUDE.md ("remove the old name completely … delete the old R file, remove it
from NAMESPACE, remove its manual page … Leftover compatibility files are technical
debt"). **`chr_meta` is gone** — split into `chr_inheritance` + `chr_recombination`.
`copy_mode_M`/`copy_mode_F`, `hemi_parent`, `recombines_M`/`recombines_F` (as
`chr_meta` columns), `define_chr()`, and `resolve_chr_copy_count()` **cease to
exist**. Nothing is aliased, wrapped, or kept "for later." This inventory is
grounded in a full-tree grep (below); it is exhaustive by construction.

**Delete outright:**
- `R/define_chr.R` (superseded by `R/define_chromosome.R`).
- `man/define_chr.Rd` and `man/resolve_chr_copy_count.Rd` (their functions no
  longer exist — do **not** let `document()` leave stale `.Rd` behind).
- The `CREATE TABLE chr_meta` statement and `resolve_chr_copy_count()`.
- Every read of `copy_mode_M/_F`, `hemi_parent`, `recombines_M/_F` **as `chr_meta`
  columns** (they are replaced by resolved `from_parent_*` / `recombines`).

**Scrub old names from prose/comments/errors** (covered in file-by-file):
`R/open_pop.R` (table-list comment), `R/ploidy_helpers.R` (docstring + error text).

**Tests** touching old names: the eight in Test changes, plus rename
`test-define_chr.R` → `test-define_chromosome.R`.

**Vignette** `vignettes/tidybreed-introduction.Rmd` — prose (L156-157, 461, 597),
the `?define_chr` cross-reference, and the **live `{r chr_meta}` chunk
(`get_table("chr_meta") |> collect()`, L160-161) which will error post-refactor** —
rewrite to `define_chromosome()` and query `chr_inheritance`/`chr_recombination`.

**Regenerate (roxygen, never hand-edit):** `man/` for `add_offspring`,
`add_founders`, `define_genome`, `assert_ploidy_2`, `is_plain_autosome`,
`get_chr_meta_map`, `assert_qtl_autosomal`, `dot-create_core_tables`.

**Deliberately retained** (the substring `chr_meta` survives here *by design* until
the Q7 rename sweep, so they are the *only* allowed matches): the file
`R/chr_meta_helpers.R` and the function `get_chr_meta_map()`. These do not use the
bare table name.

**Verification gate — must return nothing before the commit lands:**

```bash
grep -rnE '\bchr_meta\b|copy_mode_[MF]|recombines_[MF]|\bhemi_parent\b|\bdefine_chr\b|\bresolve_chr_copy_count\b' \
  R/ tests/ man/ vignettes/ NAMESPACE CLAUDE.md
```

`\bchr_meta\b` matches the bare table name but **not** `chr_meta_helpers` or
`get_chr_meta_map` (underscore is a word char, so no boundary), and `\bdefine_chr\b`
does **not** match `define_chromosome`. Any hit is leftover old-schema debt to fix
before merging.

---

## Test changes

| Test | Change |
|---|---|
| `test-define_chr.R` → `test-define_chromosome.R` | Rename; per-parent API; routing to the correct table; upsert/fallback (Finding D); `from_parent_*` sum `≤ 2` validation; **`overwrite = FALSE`** errors on an existing NULL-safe key across all four NULL/non-NULL `(sex,line)` combos (v5 #6); **transaction rollback** — an invalid write leaves the prior valid row intact (v5 #7). |
| `test-sex_chromosomes.R` | New call style; behavior unchanged (row counts per sex, inheritance, achiasmy). |
| `test-define_genome.R` | Two seeded tables; genome-wide `recombines_*` seeding both branches (Finding B); inheritance seed `1,1`. |
| `test-make_gametes_parity.R`, `test-long-schema.R`, `test-ploidy_helpers.R`, `test-add_tbv.R`, `test-define_additive_effects.R`, `test-extract_genotypes.R`, `test-add_dosage.R`, `test-haplotype_write_helpers.R` | Audit for direct `chr_meta`/`copy_mode_*`/`recombines_*`/`hemi_parent`/`resolve_chr_copy_count`/`define_chr()` reads; update to the new tables/API. Seeded output stays reproducible (within-current-code). |
| `test-chr_meta_resolve.R` (new) | Both resolvers' precedence (all four tiers); missing-chr error; every validator check incl. the deterministic shadowing algorithm (conflicting `(sex,NULL)`+`(NULL,line)` require `(sex,line)`; **matching values do not spuriously fail** — v5 #5) and NULL-key dup; targeted `(chr,NULL,line)` delete. |
| `test-resolver-context.R` (new) | Pin both contexts **in a cross where offspring and both parents differ in sex *and* line** (v5 #4): inheritance resolved by the **offspring's** sex+line, recombination by the **producing parent's** sex+line; prove a sex-specific copy rule and a line-specific recombination rule do **not** interact. |
| `test-kernel-boundary.R` (new) | Resolved copy total exactly equals emitted `ind_haplotype` row count for autosome / X-Y / Z-W / absent-chr / organelle (v5 #2/#3); a stored `2,0` fails with the explicit "not implemented" error rather than duplicating/dropping strands; if `parent_sex` is retained, an attempt to get gamete-role-specific recombination from one hermaphrodite is rejected, not silently accepted (v5 #1). |
| `test-species_coverage.R` (new, recommended) | Build ZW, X0, conifer-organelle, and *Drosophila* (♂ achiasmy + chr4) setups from this review; assert resolved per-parent copy counts and recombination. A living record that the coverage claims hold. |

Reproducibility: storage/API reshape, not an algorithm change. Same seed
reproduces; R↔Rcpp parity holds; **no** comparison against pre-change output.

---

## Rollout

1. Resolvers + validators (+ tests) — pure addition to `chr_meta_helpers.R`.
2. `define_genome()`: two `CREATE TABLE`s (with the row-local `CHECK`/`NOT NULL`
   constraints, v5 #8) + seeds + Finding-B genome-wide `recombines_*`.
3. Rename + rewrite `define_chromosome()` (routing to both tables; transaction +
   rollback; `overwrite=FALSE` semantics).
4. Update consumers (`add_founders`, `add_offspring` incl. kernel-boundary guards,
   `assert_qtl_autosomal`, `get_chr_meta_map`, `is_plain_autosome`; remove
   `resolve_chr_copy_count`).
5. `schema.R`, `sql_utils.R`.
6. Tests — gate on `test-sex_chromosomes.R`, `test-make_gametes_parity.R`,
   `test-long-schema.R`, `test-resolver-context.R`, `test-kernel-boundary.R`,
   `test-species_coverage.R`.
7. **Old-schema cleanup** (see the dedicated section): delete `R/define_chr.R`,
   `man/define_chr.Rd`, `man/resolve_chr_copy_count.Rd`; scrub `open_pop.R`,
   `ploidy_helpers.R`, the three extra tests, and the vignette; then run the
   **verification grep gate** — it must return nothing but the deliberately-retained
   `chr_meta_helpers.R`/`get_chr_meta_map()` (which the gate's `\b` anchors already
   exclude). This gate is a required check before the commit lands.
8. Docs: `CLAUDE.md`, `NEWS.md`, `DESCRIPTION`, `man/` regen.
9. *(Separate follow-up commit — Q7)* `chr` → `chromosome` column/table rename:
   `chr_inheritance`/`chr_recombination` → `chromosome_inheritance`/
   `chromosome_recombination`, `chr_name` → `chromosome_name`,
   `genome_meta.chr` → `chromosome_number`, joins, helpers, and tests.

---

## Open questions (for review)

1. **Surrogate PK?** Proposed none on either table (logical keys only).
2. **Upsert mechanism**: **DECIDED — transaction-backed delete-then-insert using
   `IS NOT DISTINCT FROM`** for NULL-safe matching on each logical key.
3. **Resolve placement**: **DECIDED — two pre-resolved maps** in
   `get_chr_meta_map()` (Q3), zero SQL in per-individual loops.
4. **`line_name` now or later?** **DECIDED — reserve now, scoped.** Nullable
   `line_name` on both tables + resolver support + `define_chromosome(line_name=)`.
   It is a **pure-line / founder default**, validated for shadowing, **not** a
   crossbred structural-variant model (per-copy SV support is future — Codex #4).
5. **Validator both-sexes coverage**: **DECIDED — hard error**, applied per table
   (inheritance resolves for both offspring sexes; recombination for both parent
   sexes). Tracks the sex domain `ind_meta.sex` admits (today `{M,F}`).
6. **Finding B**: **DECIDED — accepted.** Genome-wide `recombines_M`/`recombines_F`
   on `define_genome()`, seeded into `chr_recombination`.
7. **`chr` → `chromosome` schema-wide rename**: **DECIDED — yes, separate follow-up
   commit** (now includes the two new tables).
8. **(Codex #2) Parent contribution explicit?** **DECIDED — yes.** `copy_mode`/
   `hemi_parent` removed; replaced by explicit `from_parent_1`/`from_parent_2`.
9. **(Codex #3) Split the `sex` meanings?** **DECIDED — yes, two tables.**
   `chr_inheritance.offspring_sex` vs `chr_recombination.parent_sex`.
10. **(Codex #5) Resolver shadowing?** **DECIDED — two tables** give independent
    resolution, eliminating cross-concern shadowing; within-table sex-vs-line
    shadowing is caught by each validator's deterministic check (dormant until
    `line_name` is used).
11. **(v5 #1) `parent_sex` vs gamete role?** **RECOMMENDED — narrow scope** (keep
    `parent_sex` = producer's genetic sex; role-specific recombination in one
    hermaphrodite → Boundaries). Both reviewers concur. **Surfaced for your call**
    — the only scope decision this round; a `gamete_role` dimension is additive
    later if wanted.
12. **(v5 #2) `2,0` uniparental disomy?** **DECIDED — storage-expressible,
    kernel-future.** Allowed in the table; `add_offspring()` errors explicitly.
13. **(v5 #3) Ploidy validation value?** **DECIDED — constant `≤ 2`** (R + DB
    `CHECK`) for this diploid release, plus `assert_ploidy_2()` and a kernel
    total-copies == emitted-rows assertion. Moves to per-offspring/parent context
    when ploidy > 2 lands.
14. **(v5 #4) Which line selects each table?** **DECIDED — inheritance by the
    offspring's line; recombination by the producing parent's line.** Pinned in
    `test-resolver-context.R`.
15. **(v5 #6) `overwrite = FALSE`?** **DECIDED — error if the exact NULL-safe key
    exists, else insert;** tested across all four `(sex,line)` NULL combinations.
16. **(v5 #7) Transaction safety?** **DECIDED — write + validator in one
    transaction; rollback on any failure.**
17. **(v5 #8) DB constraints / FK?** **DECIDED — add row-local `CHECK`/`NOT NULL`
    constraints; keep the `chr_name` FK as an R validator orphan check** (matches
    the package convention), not a DB-level foreign key.

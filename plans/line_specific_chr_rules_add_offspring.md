# Plan: honor line-specific chromosome rules in `add_offspring()`

Status: proposed, for review (user + Codex). Nothing implemented yet.
Target version for the implementation: 0.76.12 (this plan is committed as 0.76.11).

## 1. The bug

`add_offspring()` ignores a line-specific `chr_inheritance` or
`chr_recombination` rule on any chromosome that is a plain autosome under the
line-agnostic rules.

How it happens (`R/add_offspring.R`, §8, around lines 266–294):

- `chr_rules_map <- get_chr_rules_map(pop$db_conn)` resolves only the
  `line_name IS NULL` rows.
- `is_autosome_id` classifies every chromosome from that line-agnostic map with
  `is_plain_autosome()`.
- A chromosome classed as a plain autosome goes to the gamete kernel
  (`make_gametes_batch()`, C++/R parity pair). The kernel always gives every
  offspring one recombined copy from each parent. It never looks at a rule.
- Only the special-chromosome path (`special_chr_names`) resolves rules per line
  through `chr_rules_for()`: inheritance at the offspring's line, recombination
  at the producing parent's line.

The code comment says so directly: "In v1 no line-specific chromosome rules
exist ... the seam is reserved, not yet exercised." The rules can be written now
(`define_chromosome(line_name = )`), so the seam is reachable.

### Observed

- Chromosome `"1"` is a default autosome.
- Before any F1 animal exists, `define_chromosome(pop, "1", offspring_sex = "M",
  line_name = "F1", from_parent_1 = 0, from_parent_2 = 1)`.
- `add_offspring()` of F1 males: each male gets two recombined copies of chr 1
  (one per parent). There is no error.

### Second symptom (same cause)

`add_founders()` *does* resolve rules at the founders' line
(`get_chr_rules_map(conn, line_name)`, `R/add_founders.R:272`). So a founder
line with a line rule that makes chr 1 hemizygous gets founders with one copy.
When one of them becomes a parent, `add_offspring()` treats chr 1 as an
autosome and stops with "Parent '...' has N autosomal ind_haplotype rows;
expected 2N". The same rule therefore works for founders and fails for their
offspring.

### What is *not* affected

- Chromosomes that are already special under the line-agnostic rules (X, Y, Z,
  W, MT set up for all lines). They take the special path, which resolves line
  rules correctly today.
- Every population with no line-specific rows. That is all current tests,
  vignettes and examples.

## 2. The fix

Classify a chromosome as a plain autosome only if it is plain for **every line
that has rules**, not just under the line-agnostic default.

```r
# R/chr_meta_helpers.R (new internal helper)
#' Chromosomes that are plain autosomes for every line with rules
.plain_autosome_chrs <- function(conn) {
  lines <- DBI::dbGetQuery(conn, paste0(
    "SELECT line_name FROM chr_inheritance WHERE line_name IS NOT NULL ",
    "UNION ",
    "SELECT line_name FROM chr_recombination WHERE line_name IS NOT NULL ",
    "ORDER BY line_name"))$line_name
  maps <- lapply(c(list(NULL), as.list(lines)),
                 function(ln) get_chr_rules_map(conn, ln))
  chrs <- names(maps[[1L]])
  chrs[vapply(chrs, function(cn) {
    all(vapply(maps, function(m) is_plain_autosome(m[[cn]]), logical(1)))
  }, logical(1))]
}
```

In `add_offspring()` §8:

```r
plain_chrs     <- .plain_autosome_chrs(pop$db_conn)
is_autosome_id <- unname(chr_name_by_id[chr_id_keys]) %in% plain_chrs
```

The line-agnostic `chr_rules_map` local is then unused and goes. `chr_rules_for()`
stays as it is. Rewrite the stale comment to say that classification covers
every line with rules and that per-individual decisions resolve through
`chr_rules_for()`.

Why this works: the special path already does the right thing for any rule.

- It reads the offspring's copies from `chr_rules_for(offspring line)`.
- It reads each parent's own copies from `chr_rules_for(parent line)`.
- It reads recombination from `chr_rules_for(parent line)`.
- A chromosome that is `1,1` and recombining for one line and special for
  another simply takes the special path for everyone. A `k_p = 2`, recombining
  parent goes through `.draw_chr_recombination()`, which is a correct
  recombined gamete.

So nothing new has to be written for the per-individual logic. Only the
classification changes.

### Design choices for review

**D1. Classify over all lines with rules, not over the lines in this call.**

Classifying from only the lines present in the current mating would make one
offspring's chromosomes depend on which other matings are in the same call.
That breaks the rule that seeded output is a function of stored inputs. The
autosome kernel draws all autosomes from one kind = 1 stream per gamete, so
moving a chromosome in or out of the autosome set changes the draws for the
others.

With D1, classification is a function of the stored rule tables alone. The
consequence is that adding a rule for line B changes line A's seeded offspring
from then on (the chromosome moves to the special path for every line, so the
draws change). That is allowed by the reproducibility contract (forward-looking,
within the current stored inputs) and is the price of a stable classification.

**D2. Resolve rules, not raw rows.**

The alternative is "plain iff every stored row for the chromosome is `1,1` /
`TRUE`". It is simpler SQL but stricter: a row fully shadowed by a more specific
row would still mark the chromosome special. Resolving per line over
`{NULL} ∪ lines with rules` is exact and matches `.egv_diploid_loci()` in
`R/extract_genetic_variance.R`, which already classifies this way. Either is
correct (the special path handles a plain chromosome correctly), so this is
only about how often the slower path is used.

**D3. No per-line kernel.**

Another option would be to give the kernel a per-gamete set of autosomes, so a
chromosome stays on the fast path for lines where it is plain. That changes the
C++/R parity pair and the stream layout. It isn't worth it until someone has a
real use case with large special sets. The cost of D1 is speed only: the
special path is an R loop per offspring, and it only runs for chromosomes that
some line has made special.

## 3. Reproducibility

- **No line-specific rows:** `.plain_autosome_chrs()` returns exactly the
  chromosomes `is_plain_autosome(get_chr_rules_map(conn)[[cn]])` did. The
  classification, the kernel inputs and the RNG streams are the same, so output
  is bit-identical to today. `test-parity.R` and
  `test-genome-effects-determinism.R` keep covering this. Nothing in them
  changes.
- **With line-specific rows:** output changes, which is the point of the fix.
  No golden-from-old comparison.
- **RNG:** classification is pure SQL reads (`dbGetQuery`), which does not
  touch R's RNG. The stream keys (`.gamete_stream_id(o, r, kind)`) are
  unchanged.
- **R ↔ C++ parity:** the kernel's inputs come after classification, and the
  kernel itself is untouched.

## 4. Tests (`tests/testthat/test-add_offspring.R`, new block)

Every test builds the population with `make_test_pop()` (or the local
equivalent the file already uses), and sets line rules before the line has
individuals, because the 0.76.10 guard refuses anything else.

1. **Line inheritance rule on a default autosome is honored.**
   - Rule: chr 1, `offspring_sex = "M"`, `line_name = "F1"`, `0,1`.
   - Founders of line A (males and females). `add_offspring()` with both F1
     and A matings in one call.
   - F1 males: chr 1 rows only for `parent_origin = 2`.
   - F1 females and all A offspring: chr 1 rows for both origins.
   - Other chromosomes: two copies for everyone.
2. **A parent made hemizygous by a founder line rule can breed.**
   - Rule set before founders of line H (`offspring_sex = "M"`, `0,1`).
   - H founder males carry one copy (already true today). Crossing an H male
     with a female works (today it errors with the "autosomal ind_haplotype
     rows" message). The offspring get the expected copies.
3. **Line recombination rule on a default autosome is honored.**
   - Rule: chr 1, `parent_sex = "M"`, `line_name = "A"`, `recombines = FALSE`.
     It may be set after founders (the guard only covers inheritance).
   - `add_offspring(store_crossovers = TRUE)`: no `ind_crossover` rows for chr 1
     with `parent_origin = 1` from line-A sires. Each offspring's paternal chr 1
     equals one of its sire's two haplotypes exactly.
   - Dams, and sires of another line, still recombine on chr 1 (at least one
     crossover across the offspring, with a seed and enough offspring that this
     is not flaky; or check against a line without the rule).
4. **Same seed, same output** with line rules present: two identical runs,
   `expect_identical()` on `ind_haplotype` and `ind_crossover`.
5. **Classification is a function of stored rules only (pins D1).**
   - Build two identical populations, both with a line-B rule on chr 1.
   - In one, call `add_offspring(seed = s)` with line-A matings only. In the
     other, use the same A matings followed by line-B matings.
   - The A matings are rows `o = 1..n` in both calls, and streams are keyed on
     `o`, so the A offspring's haplotypes must be identical.
   - A per-call classification would fail this: chr 1 would be on the kernel in
     the first call and on the special path in the second.
6. **No line rules: unchanged classification.** A direct test of the helper
   that `.plain_autosome_chrs()` equals the line-agnostic classification on a
   population with X/Y/MT rules and no line rows.

Then re-run the suites that touch chromosomes and offspring:
`test-define_chromosome.R`, `test-kernel-boundary.R`, `test-add_founders.R`,
`test-add_offspring*.R`, `test-parity.R`, `test-extract_genotypes.R`,
`test-add_dosage.R`, and the full suite before the commit.

## 5. Docs and bookkeeping

- `?define_chromosome`: the 0.76.9 page already says `line_name` applies in
  `add_founders()` and `add_offspring()`. After the fix that is true; re-read
  it once. Add one sentence: a chromosome that any line makes special is drawn
  on the slower per-offspring path for every line.
- `.claude/skills/tidybreed-schema/SKILL.md` (chr rule resolution) and
  `.claude/skills/tidybreed-api/SKILL.md` (`add_offspring()`): describe the
  classification over all lines with rules.
- `R/add_offspring.R` comments in §8 and §9 that say "v1 has only the default"
  or "reserved, not yet exercised": rewrite them to match the code.
- `DESCRIPTION` 0.76.12, `NEWS.md` entry under bug fixes.

## 6. Out of scope (flagged, not fixed here)

Line-agnostic checks elsewhere also read only the `line_name IS NULL` rules:

- `assert_qtl_autosomal()` (`R/chr_meta_helpers.R`): the generators' "QTL must
  be diploid" check.
- `.ge_validate_dominance_ploidy()` (`R/genome_effects_helpers.R`): dominance
  needs a diploid locus.
- `.ge_reachable_copy_counts()` (`R/define_genome_effect_terms.R`): which copy
  counts a dosage state can take.

Each is documented as population-level on purpose ("an effect definition is
population-level"). But with a line rule that makes chr 1 hemizygous in line
H, these checks would accept a population-wide dominance term or a generated
QTL on chr 1, and H animals would then carry one copy. `.egv_diploid_loci()`
already classifies over all lines with rules, so the code base is not
consistent here.

Question for review: should these use the same "every line with rules" rule
(probably via a shared helper next to `.plain_autosome_chrs()`)? I'd do that as
a separate change after this one, since it changes what the genome-effect
writers accept rather than how offspring are made.

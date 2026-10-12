# tidybreed documentation proposal

Drafted 2026-10-12 against the v0.76.4 working tree. This is a proposed structure
and writing sequence, not an implemented book or an additional genetics roadmap.

## Recommendation

Create **one Quarto book, "Simulating Breeding Programs with tidybreed"**, with
three parts: a step-by-step cookbook, methodology, and validation/case studies.
Keep pkgdown as the function reference. Readers should be able to follow the
cookbook without reading every derivation, and jump from a recipe to the exact
method and assumptions behind it.

The methodology part should be complete enough to support a methods paper.
A later paper can draw on its derivations, citations and validation figures;
it need not become a second independently maintained technical manual.
Start with HTML, with a downloadable PDF once the content settles. Quarto books
support multiple output formats and chapter cross-references:
[official book documentation](https://quarto.org/docs/books/).

## What already exists

| Material | Use in the new documentation |
|---|---|
| `vignettes/tidybreed-introduction.Rmd` | Starting material for the first complete workflow; it executes code, covers the two object types, traits versus observations, offspring, genotypes, indexes and close/restore |
| `vignettes/genetic-models.Rmd` | Starting material for known coefficients, the two generators, the additive floor and target comparisons |
| Roxygen / `man/` / `_pkgdown.yml` | Canonical argument and return-value reference; link to it instead of duplicating every signature in the book |
| `tests/testthat/` | Evidence for numerical contracts, scope resolution, failure behaviour and determinism; turn selected cases into readable explanations |
| `dev/benchmarks/` | Starting scripts for explicitly versioned performance results |
| `plans/import_qtl_effect_methods*.md` | Design rationale and derivation leads; verify against implemented code and original sources before turning them into methods text |
| `vignettes/swine/swine-time-based-age-at-puberty-sex-semen.R` | An advanced case-study candidate after auditing it against the current API and extracting a small runnable example |
| `tools/quarto/legacy/` | Historical material only; do not use its old API examples as current documentation |

The current introduction provides one generation of offspring, but does not teach
a full repeated selection-and-mating cycle. The main missing tutorial is that
complete loop, including how to inspect and interpret its results.

## Part I: Cookbook

Use a small livestock example throughout the core chapters: one quantitative
trait, two sexes, a modest genome, and a few generations. Extend that same
example to correlated traits and non-additive genetics. Give independent recipes
a visible setup so readers can start there without executing earlier chapters.

| Chapter | What the reader finishes with |
|---|---|
| 1. Your first breeding program | A complete small simulation: genome, founder pool, animals, trait, additive effects, phenotypes, parent selection, matings, offspring, and a genetic-trend plot |
| 2. Working with a population database | Confidence with `tidybreed_pop`, `tidybreed_table`, lazy filters, identifiers, custom columns, schema inspection, file-backed runs and close/restore |
| 3. Genomes and founder populations | Chosen locus/map settings, founder frequencies and LD; understanding pool haplotypes versus founder individuals and line fallback |
| 4. Choosing a genetic model | Known coefficients versus sampled/calibrated effects; additive and multi-trait models; dominance and A x A; pair designs; genic versus realised anchors |
| 5. Measuring the genetics you generated | Targets joined to measured covariance blocks, missing-block checks, cohort versus generation base, `between_components`, and what "exact" means |
| 6. Recording phenotypes | Genetics versus observation model; continuous/count/categorical records, means/residuals, fixed/random effects, repeated records, conditional residuals, liability and derived phenotypes |
| 7. Selection and mating over generations | Candidate filters, phenotype/additive-value/index selection, explicit parent IDs, mating plans and offspring; genetic gain, diversity and replicate uncertainty |
| 8. Chips, genotypes and genetic evaluation | Chip design, recording/extracting genotypes, and an optional BLUPF90 workflow with its own setup and model assumptions |
| 9. More complex programs | Maternal/social components, line-specific additive variants, crossbreeding and common non-additive effects; separate recipes rather than one enormous script |
| 10. Scenarios, replicates and analysis | Ordinary R configuration, separate run directories, archiving, restarting, comparing designs, and practical memory/runtime limits |

Do not imply that user-managed YAML is a package execution engine, that mutation
is implemented, or that generated line-specific dominance/A x A effects already
exist. State the supported inheritance/ploidy combinations from the current
validation rules. Treat pedigree/inbreeding analysis as explicit analysis code
or an external-tool recipe wherever the package has no corresponding function.

Every recipe uses the same short structure:

1. **Goal and biological assumptions.** What question the simulation answers.
2. **Inputs and setup.** Libraries, seed, population state and dependencies.
3. **Runnable steps.** Explain why each public function is called and what it writes.
4. **Check the result.** A small table/plot and a check that could detect a mistake.
5. **Interpretation and variations.** What the result means, what to change next,
   and a link to the corresponding methods section.

Keep code visible by default for the core steps. Offer each complete recipe as
a downloadable R script generated from the same source, including setup and
cleanup. Avoid hidden helper functions that make a copied recipe impossible to
run independently.

In the persistence/replicate recipes, distinguish reopening a working database
from analyzing an archive. Default archiving discards haplotype/genotype rows;
use `(replicate, id_ind)` in archive joins because individual IDs repeat between
replicates. Fix the missing archive schema descriptions/misleading restore
warning recorded in TODO before using that round trip as the polished example.

## Part II: Methodology

Organize this by scientific concepts, rather than R filenames or development
phases. Define notation once and use it consistently with the tables/output.

| Chapter | Required content |
|---|---|
| Simulation state and inheritance | State transitions; locus/map representation; allele-copy origins; founder sampling; meiosis, recombination and supported ploidy rules |
| Genetic-value model | Functional versus statistical coding; additive, dominance, indicator and interaction terms; centres, scope precedence and total genetic values |
| Additive calibration | Drawn architecture versus calibrated coefficients; `B' M B = G`; genic and realised anchors; rank/PSD feasibility; exact shared versus approximate union models |
| Dominance and A x A calibration | NOIA re-expression, induced additive effects, highest-order-first calibration, the architecture-specific additive floor, pairs/hubs and inbreeding-depression limitations |
| Measurement and population change | Cohort re-projection; HWE/LE assumptions; LD and cross-component covariance; changing frequencies; sample versus population variance; partial/scoped decomposition |
| Observation models | Trait/phenotype separation; residual and random-effect covariance; repeated/conditional draws; contributor formulas; count transformations, categorical thresholds and prevalence approximation |
| Evaluation and selection | Generation truth versus evaluation assumptions; BLUPF90 inputs/outputs; selection indexes and the meaning of reported accuracy; mating designs supplied by the user |
| Computation and reproducibility | DuckDB storage and lazy queries, in-memory calibration limits, transactions/failure contracts, RNG streams, exact accumulation, version/platform conditions and external binaries |

For each method: model/equations, assumptions, algorithm, a worked numerical
example, an independent check, limitations and references. Distinguish delivered
generation targets from quantities merely reported or approximated.

Important explanations to get right after the effect-function update:

- Generation targets are not automatically the covariance of finite founders,
  later generations, selected animals or a different line.
- A realised-anchor calibration and stored HWE-referenced components are
  different descriptions; explain the extractor's re-projection.
- Under LD, the NOIA additive contrast component is not the joint least-squares
  breeding-value projection. Neither the `additive` name nor a zero
  `between_components` row proves otherwise.
- Manual functional models can have functional `ind_tgv` components rather
  than statistical breeding values; their total remains the evaluated model.
- A generated owner's presence alone does not certify a target comparison;
  anchor, reference population, active variants and measured quantity must match.
- Explicit zero blocks and absent blocks differ. Generated model replacement,
  target persistence and removal must be explained together.
- Liability prevalence uses a distributional approximation; a target variance
  does not guarantee the requested prevalence for a few major QTL or a selected
  population.

Build `references.bib` from verified original papers and software citations.
Trace calibration/NOIA claims to their actual sources, separate mathematical
results from tidybreed implementation choices, and acknowledge imported methods.
The import plan's section 14 records unreported upstream comparison findings:
resolve their evidence/reporting status before publishing those comparisons.
Do not turn development notes into unsupported claims about another package.

## Part III: Validation and case studies

This is the research-paper evidence, presented with runnable code and modest
default sizes. Larger experiments use explicit separate benchmark scripts.

- Hand-calculated one-/two-locus models versus evaluated values.
- Multi-trait target covariance versus delivered covariance under each anchor.
- Dominance/A x A calibration, additive-floor failures and numerical limits.
- LD, drift and selection: how genic and realised reports differ over time.
- Common versus line/parent-scoped models, including reciprocal crosses where
  supported; no claim that scoped non-additive generation is implemented.
- Phenotype moments, correlated residuals, maternal/social components and
  liability frequencies with uncertainty.
- A multi-generation selection study with replicated genetic trends.
- A separately executed BLUPF90 example, recording binary versions and inputs.
- Runtime, RAM and database size against animals, loci, traits and pairs. Record
  hardware, thread settings, versions and seeds; avoid blanket "larger than RAM"
  promises for operations that materialize design matrices.

Use assertions appropriate to each claim: algebraic checks for exact contracts,
sampling uncertainty for stochastic results. Cross-package comparisons are
optional and need matched biological models/coding/anchors, verified software
versions and shared inputs; seeded draws alone are not a valid comparison.

## Source layout and execution

Keep `book/` beside `R/`, `tests/` and `vignettes/`, with `_quarto.yml`,
`index.qmd`, `cookbook/`, `methods/`, `validation/`, `references.bib` and derived
downloadable scripts. Build to `book/_book/`; exclude book sources/artifacts from
the R source package with `.Rbuildignore` and generated output from Git.

Install the current checkout for CI and render against that installation. Each
chapter runs in a fresh session and creates its own small population, with no
dependence on a live connection or a `.duckdb` left by an earlier render. Use
temporary output locations and close connections explicitly. Include one visible
file-backed persistence recipe; other examples can be in memory for speed.

Allow freezing/caching for local writing, but execute all core chapters from
scratch when validating package changes. Quarto's `freeze: auto` detects document
changes, not every change to package code or external inputs, so a frozen build
alone cannot check the examples:
[official execution documentation](https://quarto.org/docs/projects/code-execution.html).
Keep expected error demonstrations local and explicit; ordinary execution errors
must fail the build.

Make external BLUPF90 setup optional for the core book. Publish its executed
example through a separate integration run once the evaluation design is settled;
show dependency requirements and recorded outputs rather than silently skipping
an advertised check. Do not download proprietary/external binaries during normal
book rendering.

Publish the book at the existing site's `/book/` path and add a pkgdown navbar
link. Stage pkgdown and book outputs together so one deployment cannot overwrite
the other. PRs validate builds; main/releases publish. Identify the package
version/commit used for each published book; retain a versioned release artifact
for citations and reproducibility. Update `inst/CITATION` when there is an actual
book URL/release to cite.

## Writing sequence and completion criteria

1. Clean the current README/citation/build issues and decide scope/RNG wording.
2. Deliver a small vertical slice: first complete breeding cycle, genetic-model
   recipe, variance-check recipe and matching calibration/measurement methods.
   Render them from a clean session, review outputs and make scripts runnable.
3. Expand the core workflow through repeated selection, file persistence and
   replicates. Adapt the existing vignettes; keep short package entry points and
   derive shared examples rather than maintaining copied implementations.
4. Add phenotype, crossbreeding and advanced contributor recipes alongside their
   methods. Complete Q20 before finalizing configurable evaluation documentation.
5. Add validation figures, citations, benchmark provenance and the optional solver
   integration example. Then produce the downloadable PDF and paper outline.

The first documentation milestone is usable when a reader can install the
package, run a complete breeding program, change its genetic model, verify the
results, and find the assumptions without reading internal plans. The full
methods milestone is complete when each supported biological/computational
method has equations, references, runnable evidence and stated limits.

Mutation, additional architecture samplers, new extractors, higher-order
epistasis and other deferred biology are not requirements for these milestones.

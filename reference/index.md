# Package index

## Population lifecycle

Create, reopen, and close the DuckDB-backed population that holds every
table.

- [`open_pop()`](https://austin-putz.github.io/tidybreed/reference/open_pop.md)
  : Open a new breeding population
- [`restore_pop()`](https://austin-putz.github.io/tidybreed/reference/restore_pop.md)
  : Reconnect to an existing tidybreed database
- [`close_pop()`](https://austin-putz.github.io/tidybreed/reference/close_pop.md)
  : Close tidybreed population database connection
- [`archive_replicate()`](https://austin-putz.github.io/tidybreed/reference/archive_replicate.md)
  : Archive a completed replicate and reset the working database

## Genome and genetic map

Define loci, chromosomes, per-chromosome inheritance rules, and the
founder haplotype pools that individuals are sampled from.

- [`define_genome()`](https://austin-putz.github.io/tidybreed/reference/define_genome.md)
  : Define the genome structure of a breeding population
- [`define_chr()`](https://austin-putz.github.io/tidybreed/reference/define_chr.md)
  : Define a chromosome's inheritance rule
- [`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)
  : Define founder haplotypes for a tidybreed population

## Individuals

Create founders and simulate meiosis to produce offspring from a mating
plan.

- [`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md)
  : Add founder individuals to population
- [`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md)
  : Add offspring via recombination

## Traits — genetic layer

Define the underlying genetic quantities: additive variance, QTL
effects, and genetic covariances between traits. Nothing here describes
an observation.

- [`define_trait()`](https://austin-putz.github.io/tidybreed/reference/define_trait.md)
  : Define a genetic component trait
- [`define_trait_simple()`](https://austin-putz.github.io/tidybreed/reference/define_trait_simple.md)
  : Define a trait with QTL and sampled effects in one call
- [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
  : Define additive QTL effects for one or more traits
- [`define_effect_cov_matrix()`](https://austin-putz.github.io/tidybreed/reference/define_effect_cov_matrix.md)
  : Define a variance-covariance matrix for any named effect

## Phenotypes — observation layer

Define what is actually recorded on an animal: distribution, mean,
residual variance, sex expression, and any fixed or random effects.

- [`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
  : Define an observed phenotype
- [`define_residual_cov()`](https://austin-putz.github.io/tidybreed/reference/define_residual_cov.md)
  : Define residual covariance entries for observed phenotypes
- [`define_effect_intercept()`](https://austin-putz.github.io/tidybreed/reference/define_effect_intercept.md)
  : Set the intercept (population mean) for a phenotype
- [`define_effect_fixed_class()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_class.md)
  : Define a discrete fixed-class effect in a phenotype model
- [`define_effect_fixed_cov()`](https://austin-putz.github.io/tidybreed/reference/define_effect_fixed_cov.md)
  : Define a continuous covariate (regression) effect in a phenotype
  model
- [`define_effect_random()`](https://austin-putz.github.io/tidybreed/reference/define_effect_random.md)
  : Define a random group effect in a phenotype model

## Simulation output

Compute and store true breeding values, phenotype records, and estimated
breeding values from an external evaluation.

- [`add_tbv()`](https://austin-putz.github.io/tidybreed/reference/add_tbv.md)
  : Compute and store true breeding values without writing phenotypes
- [`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
  : Generate phenotype records for a subset of individuals
- [`add_ebv()`](https://austin-putz.github.io/tidybreed/reference/add_ebv.md)
  : Add estimated breeding values to a tidybreed population

## Genotypes and SNP chips

Define chips, record which animals were genotyped, and pull dosage
matrices out for analysis or export.

- [`define_chip()`](https://austin-putz.github.io/tidybreed/reference/define_chip.md)
  : Define a SNP chip

- [`add_genotypes()`](https://austin-putz.github.io/tidybreed/reference/add_genotypes.md)
  : Mark animals as genotyped on a SNP chip

- [`add_dosage()`](https://austin-putz.github.io/tidybreed/reference/add_dosage.md)
  :

  Materialize genotype dosage values into `ind_genotype`

- [`extract_genotypes()`](https://austin-putz.github.io/tidybreed/reference/extract_genotypes.md)
  : Extract genotype data for individuals, by chip and/or QTL loci

## Selection index

Register index weights and apply them to EBVs, TBVs, or any table of
values.

- [`define_index()`](https://austin-putz.github.io/tidybreed/reference/define_index.md)
  : Define a selection index
- [`add_index()`](https://austin-putz.github.io/tidybreed/reference/add_index.md)
  : Compute a selection index from a tidybreed table

## Tables, columns, and queries

Read any table lazily, add or update columns, build group-level
variables, and delete rows. This is how simulation state stays in the
database rather than in parallel R objects.

- [`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md)
  : Get a table reference from a tidybreed population
- [`mutate_table()`](https://austin-putz.github.io/tidybreed/reference/mutate_table.md)
  : Add or modify columns in any population database table
- [`mutate_derived()`](https://austin-putz.github.io/tidybreed/reference/mutate_derived.md)
  : Compute and write a derived column from a cross-table join
- [`mutate_group_seq()`](https://austin-putz.github.io/tidybreed/reference/mutate_group_seq.md)
  : Assign sequential integer group IDs to filtered rows
- [`mutate_group_named()`](https://austin-putz.github.io/tidybreed/reference/mutate_group_named.md)
  : Assign named group labels to filtered rows by count or proportion
- [`mutate_group_concatenate()`](https://austin-putz.github.io/tidybreed/reference/mutate_group_concatenate.md)
  : Create a composite group column by concatenating existing columns
- [`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md)
  : Remove rows from one or more population tables
- [`define_table()`](https://austin-putz.github.io/tidybreed/reference/define_table.md)
  : Define a custom table in the population database

## Schema documentation

Inspect what every table and column means, and attach your own
descriptions that travel with the `.duckdb` file.

- [`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md)
  : View table-level descriptions for a tidybreed population
- [`describe_table()`](https://austin-putz.github.io/tidybreed/reference/describe_table.md)
  : View column-level descriptions for a tidybreed table
- [`define_schema_description()`](https://austin-putz.github.io/tidybreed/reference/define_schema_description.md)
  : Define or update a description for a table or column

## dplyr methods for tidybreed_table

Standard dplyr verbs work on the lazy table reference returned by
[`get_table()`](https://austin-putz.github.io/tidybreed/reference/get_table.md).
Filtering happens inside DuckDB; `collect()` pulls into R.

- [`filter(`*`<tidybreed_table>`*`)`](https://austin-putz.github.io/tidybreed/reference/filter.tidybreed_table.md)
  : Filter method for tidybreed_table
- [`select(`*`<tidybreed_table>`*`)`](https://austin-putz.github.io/tidybreed/reference/select.tidybreed_table.md)
  : Select method for tidybreed_table
- [`arrange(`*`<tidybreed_table>`*`)`](https://austin-putz.github.io/tidybreed/reference/arrange.tidybreed_table.md)
  : Arrange method for tidybreed_table
- [`collect(`*`<tidybreed_table>`*`)`](https://austin-putz.github.io/tidybreed/reference/collect.tidybreed_table.md)
  : Collect method for tidybreed_table
- [`pull(`*`<tidybreed_table>`*`)`](https://austin-putz.github.io/tidybreed/reference/pull.tidybreed_table.md)
  : Pull method for tidybreed_table
- [`count(`*`<tidybreed_table>`*`)`](https://austin-putz.github.io/tidybreed/reference/count.tidybreed_table.md)
  : Count method for tidybreed_table
- [`slice_head(`*`<tidybreed_table>`*`)`](https://austin-putz.github.io/tidybreed/reference/slice_head.tidybreed_table.md)
  : slice_head method for tidybreed_table
- [`slice_tail(`*`<tidybreed_table>`*`)`](https://austin-putz.github.io/tidybreed/reference/slice_tail.tidybreed_table.md)
  : slice_tail method for tidybreed_table
- [`slice_min(`*`<tidybreed_table>`*`)`](https://austin-putz.github.io/tidybreed/reference/slice_min.tidybreed_table.md)
  : slice_min method for tidybreed_table
- [`slice_max(`*`<tidybreed_table>`*`)`](https://austin-putz.github.io/tidybreed/reference/slice_max.tidybreed_table.md)
  : slice_max method for tidybreed_table
- [`slice_sample(`*`<tidybreed_table>`*`)`](https://austin-putz.github.io/tidybreed/reference/slice_sample.tidybreed_table.md)
  : slice_sample method for tidybreed_table

## Print and summary methods

Display methods for the package’s S3 classes.

- [`print(`*`<tidybreed_pop>`*`)`](https://austin-putz.github.io/tidybreed/reference/print.tidybreed_pop.md)
  : Print method for tidybreed_pop
- [`summary(`*`<tidybreed_pop>`*`)`](https://austin-putz.github.io/tidybreed/reference/summary.tidybreed_pop.md)
  : Summarize a tidybreed population
- [`print(`*`<tidybreed_table>`*`)`](https://austin-putz.github.io/tidybreed/reference/print.tidybreed_table.md)
  : Print method for tidybreed_table
- [`print(`*`<tidybreed_schema>`*`)`](https://austin-putz.github.io/tidybreed/reference/print.tidybreed_schema.md)
  : Print method for tidybreed_schema
- [`print(`*`<tidybreed_summary>`*`)`](https://austin-putz.github.io/tidybreed/reference/print.tidybreed_summary.md)
  : Print a tidybreed_summary object
- [`print(`*`<tidybreed_table_desc>`*`)`](https://austin-putz.github.io/tidybreed/reference/print.tidybreed_table_desc.md)
  : Print method for tidybreed_table_desc

# Pivot ind_phenotype rows wide for formula evaluation.

Uses pheno_number = 1 by default. If, across the entire `rows` input,
not a single row has `pheno_number == 1` (e.g. first records were
deleted), falls back to the minimum `pheno_number` per (id_ind,
phenotype_name) pair instead. This fallback is all-or-nothing over the
whole input, not decided per phenotype or per individual.

## Usage

``` r
.pivot_pheno_wide(rows, ids, pheno_names)
```

## Arguments

- rows:

  Data frame from ind_phenotype query (id_ind, phenotype_name,
  pheno_value, pheno_number).

- ids:

  Character vector. Ordered individual IDs.

- pheno_names:

  Character vector. Phenotype column names to create.

## Value

Data frame with one row per id in ids, columns = pheno_names.

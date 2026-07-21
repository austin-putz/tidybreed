# Compute the next pheno_number for each individual for a given phenotype

Compute the next pheno_number for each individual for a given phenotype

## Usage

``` r
next_pheno_numbers(pop, phenotype_name, ids)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- phenotype_name:

  Character. Phenotype name.

- ids:

  Character vector of individual IDs.

## Value

Integer vector, same length and order as `ids`.

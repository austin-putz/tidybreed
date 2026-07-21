# Print method for tidybreed_pop

Prints a grouped, human-readable summary of the population: a header
rule with the population name, the database location and connection
status, then one section each for the **genome** (loci, chromosomes,
physical/genetic length, founder haplotype pool), the **model**
(genetic-component traits, observed phenotypes, selection indices, QTL),
the **individuals** (total, broken down by sex and by line), and the
**records** written so far (phenotypes, TBVs, EBVs, index values).

Sections are shown only when their underlying data exists, so a freshly
opened population collapses to just the header, the database line, and a
pointer to
[`schema()`](https://austin-putz.github.io/tidybreed/reference/schema.md)
/
[`describe_table()`](https://austin-putz.github.io/tidybreed/reference/describe_table.md).
Counts are queried live from the DuckDB connection; a closed connection
prints as `[disconnected]` without running any queries.

## Usage

``` r
# S3 method for class 'tidybreed_pop'
print(x, ...)
```

## Arguments

- x:

  A tidybreed_pop object

- ...:

  Additional arguments (not used)

## Value

`x`, invisibly.

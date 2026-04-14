# tidybreed

<!-- badges: start -->
<!-- badges: end -->

## Overview

tidybreed is a next-generation R package for breeding program simulation that embraces tidyverse principles. It provides a flexible, user-friendly framework for simulating animal and plant breeding programs with:

- **Tidyverse-style design** using pipe operators (`%>%`)
- **Memory-efficient storage** using DuckDB backends
- **Flexible trait architectures** supporting multiple QTL and SNP panels
- **Custom metadata** allowing users to define their own fields and workflows
- **Compatible with dplyr** for easy data manipulation and analysis

## Installation

You can install the development version of tidybreed from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("YOUR-USERNAME/tidybreed")
```

## Example

This is a basic example of simulating a breeding program:

``` r
library(tidybreed)

# Initialize a genome
pop_A <- initialize_genome(
  pop_name = "A",
  n_loci = 100,
  n_chr = 2,
  chr_len_Mb = 100
) %>%
  add_founders(
    n_males = 10,
    n_females = 100,
    line_name = "A"
  ) %>%
  add_trait(
    name = "ADG",
    n_qtl = 10,
    target_add_var = 1,
    res_var = 1
  )

# Phenotype and select
pop_A %>%
  filter(Sex == "Males") %>%
  add_phenotype(name = "ADG")
```

## Documentation

See the package vignettes for detailed examples and workflows.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

MIT License

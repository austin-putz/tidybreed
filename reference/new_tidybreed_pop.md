# Create a new tidybreed population object

Constructor for the `tidybreed_pop` S3 class. This class wraps a DuckDB
connection and provides a tidy interface for breeding program
simulation.

## Usage

``` r
new_tidybreed_pop(
  db_conn,
  pop_name,
  db_path,
  tables = character(),
  run_dirs = character(0)
)
```

## Arguments

- db_conn:

  A DuckDB connection object

- pop_name:

  Character string naming the population

- db_path:

  Path to the DuckDB database file

- tables:

  Character vector of table names in the database

- run_dirs:

  Named character vector mapping tool names to their pre-created layer-4
  output directories. `character(0)` for in-memory databases or when no
  tools were declared. Always includes a `"base"` entry pointing to the
  layer-3 scenario directory when non-empty.

## Value

A `tidybreed_pop` object

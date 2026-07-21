# Assert that `name` is a valid, non-reserved SQL identifier.

Used by every function that interpolates user input into SQL to protect
against SQL injection.

## Usage

``` r
validate_sql_identifier(name, what = "identifier", reserved = character())
```

## Arguments

- name:

  Character scalar to validate.

- what:

  Human-readable role label used in error messages (e.g. `"trait name"`,
  `"effect name"`).

- reserved:

  Optional character vector of additional reserved names for the calling
  context (e.g. reserved column names in a table).

## Value

Invisible `NULL` on success; errors otherwise.

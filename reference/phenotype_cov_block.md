# Covariance blocks in `phenotype_var_comp`

A *covariance block* is a connected component of the graph whose
vertices are phenotypes and whose edges are stored pair rows in
`phenotype_var_comp` for one `effect_name` — a pair row with
`cov_value = 0` is still an edge. Every writer of that table goes
through
[`write_phenotype_cov_block()`](https://austin-putz.github.io/tidybreed/reference/write_phenotype_cov_block.md),
which runs
[`validate_phenotype_cov_block()`](https://austin-putz.github.io/tidybreed/reference/validate_phenotype_cov_block.md)
inside its transaction before any row is deleted, so a rejected
declaration never leaves a half-written block.

The rules, from `plans/sample_correlated_effects.md`:

- **D1 — a block is declared in one call, as a complete matrix.** For a
  call over phenotypes `N`, let `U` be `N` plus every member of every
  existing block that touches `N`. The call is accepted only when
  `N == U`; a fragment or a strict subset of an existing block is an
  error naming the omitted phenotypes. The matrix must be symmetric,
  finite, and positive semi-definite.

- **Strata (residual only).** Conditional rows partition into strata
  `(condition_table, condition_column, condition_level)`. Every stratum
  of a block names the same phenotypes and a block has at most one
  condition column, so which stratum an individual resolves to never
  changes block membership.

- **D3 — realization lock.** A block cannot be redefined once a draw
  exists under it: for the residual effect, any `ind_phenotype` row of a
  member with `residual_value IS NOT NULL`; for a named effect, any
  `phenotype_random_effects` row for `(effect_name, member)`. The error
  gives the
  [`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md)
  call that clears the realizations.

- **D6 — agreement.** Every phenotype in a residual block that has a
  `phenotype_meta` row carries the same `condition_change_action`.

- **Named effects (§5.6).** In a block of two or more phenotypes every
  `phenotype_effects` row for the effect is `random`, uses
  `distribution = "normal"`, and reads the same
  `(source_column, source_table)`. A 1 × 1 block is exempt, so a gamma
  or uniform effect on a single phenotype stays legal until something
  tries to join it to a second phenotype.

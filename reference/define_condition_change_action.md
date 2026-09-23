# Set `condition_change_action` for a whole residual covariance block

`condition_change_action` decides what
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)
does when a correlated phenotype's stored residual was drawn under a
different residual `condition_level` than the current record resolves to
— an animal that moved farms between two records, say. It is a property
of the **residual covariance block**, not of one phenotype: D6 requires
every member of a block to carry the same value, so this function is the
only way to change it once a block of two or more exists.

## Usage

``` r
define_condition_change_action(
  pop,
  phenotype_name,
  condition_change_action = c("error", "independent")
)
```

## Arguments

- pop:

  A `tidybreed_pop`.

- phenotype_name:

  Character scalar. Any member of the block; the value is written to
  every member.

- condition_change_action:

  `"error"` (stop when a stored residual comes from another stratum) or
  `"independent"` (drop it from the conditioning set and warn).

## Value

The `tidybreed_pop`, invisibly.

## Details

[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md)
sets the value at registration time, when the phenotype is still a block
of one. Afterwards the two are locked together — flipping either one
alone through `define_phenotype(overwrite = TRUE)` would leave the block
disagreeing, which D6 refuses, and restating the phenotype would reset
every other column of its `phenotype_meta` row along the way. This
function changes exactly one column, on every member, in one
transaction.

Unlike the covariance matrix itself, the action is **not** locked by
realized draws (D3). It governs how *future* records condition on stored
residuals and says nothing about the ones already drawn, so changing it
mid-simulation is legitimate and leaves `ind_phenotype` and
`phenotype_random_effects` untouched.

## See also

[`define_phenotype()`](https://austin-putz.github.io/tidybreed/reference/define_phenotype.md),
[`define_residual_cov()`](https://austin-putz.github.io/tidybreed/reference/define_residual_cov.md),
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# A and B share a farm-conditioned residual block, both at the default
# "error". An animal that changes farms between records should now be
# drawn independently of its stale residual rather than stopping the run.
pop <- define_condition_change_action(pop, "A", "independent")
# -> also set on B (same residual covariance block)
} # }
```

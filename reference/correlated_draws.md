# Correlated draws: block loading and the conditional resolver

The two database-independent halves of sequential correlated sampling
(`plans/sample_correlated_effects.md` §5.2–5.4):

- [`find_covariance_blocks()`](https://austin-putz.github.io/tidybreed/reference/find_covariance_blocks.md)
  reads `phenotype_var_comp` once and returns the complete covariance
  block(s) touching a set of phenotypes, one validated matrix per
  stratum. It is the only place a stored block becomes a matrix.

- [`resolve_correlated_draws()`](https://austin-putz.github.io/tidybreed/reference/resolve_correlated_draws.md)
  takes one such matrix and, for a set of entities with possibly
  different observed coordinates, returns draws of the requested
  coordinates from the Gaussian conditional on what each entity has
  already realized. It has no notion of individuals, strata or tables —
  the two adapters in `R/add_phenotype_stages.R` own those: the residual
  adapter
  ([`.ap_resolve_residuals()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_resolve_residuals.md),
  entity `(id_ind, pheno_number)`) and the named-effect adapter
  ([`.ap_resolve_named_effects()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_resolve_named_effects.md),
  entity `(effect_name, level)`).

Neither function writes. The resolver is the only one that consumes RNG,
and it does so only after every check has passed and only when there is
something to draw.

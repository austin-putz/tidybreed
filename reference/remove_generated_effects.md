# Remove one scope of generated effects

Deletes every term a generator wrote for `trait_name` at exactly one
scope, `(line_name, parent_origin)` — the scope a re-run of
[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)
would replace — **whatever its kind**: additive, dominance and
interaction terms at that scope go together. A generated model is
calibrated as a whole, so it is removed as a whole and never one
component at a time (for "A without D", re-run the generator with an
additive-only target, which re-calibrates). Every other scope stands.
This is the only way to remove a generated variant:
[`define_genome_effect_terms()`](https://austin-putz.github.io/tidybreed/reference/define_genome_effect_terms.md)
and
[`remove_rows()`](https://austin-putz.github.io/tidybreed/reference/remove_rows.md)
refuse the reserved owner `"generated"`, so that `"generated"` keeps
meaning "calibrated to the stored target".

The usual reason is a variant added by mistake, such as a re-run with a
new `parent_origin`, which **adds** a variant beside the old one rather
than replacing it. Two parent scopes at one line are a legal model, but
no stored target describes their combined variance, so
`define_phenotype(prevalence = )` refuses such a trait.

Nothing else changes. The stored targets in `trait_var_comp` stay (a
target with no terms of its kind is ignored by the `prevalence`
threshold and blocks nothing), and so do the values already in
`ind_tgv`: they describe the old model until
[`add_tgv()`](https://austin-putz.github.io/tidybreed/reference/add_tgv.md)
(or
[`add_phenotype()`](https://austin-putz.github.io/tidybreed/reference/add_phenotype.md),
which calls it) re-evaluates the individuals. Removing a line's variant
makes that line's copies fall back to the common variant, where there is
one.

## Usage

``` r
remove_generated_effects(
  pop,
  trait_name,
  line_name = NULL,
  parent_origin = NULL
)
```

## Arguments

- pop:

  A `tidybreed_pop` object.

- trait_name:

  Character vector. The trait(s) whose variant is removed.

- line_name:

  `NULL` (the population-wide scope) or one line name, as passed to
  [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md).

- parent_origin:

  `NULL` (both parents' copies), `1` (sire) or `2` (dam), as passed to
  [`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md).
  A generator with no scope arguments writes the common scope,
  `line_name = NULL, parent_origin = NULL`.

## Value

The `tidybreed_pop`, invisibly. An error, with nothing deleted, when a
trait has no generated terms at that scope.

## See also

[`define_additive_effects()`](https://austin-putz.github.io/tidybreed/reference/define_additive_effects.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# A paternal-only variant added by mistake next to the common one
pop <- remove_generated_effects(pop, "ADG", parent_origin = 1)
} # }
```

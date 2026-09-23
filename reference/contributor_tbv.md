# Contributor lookups for composite and `formula_tbv` phenotypes

The one place a contributor of a composite phenotype — the individual
itself, its dam or sire, or its group-mates — becomes a per-individual
value. Both
[`.assemble_composite_tbv()`](https://austin-putz.github.io/tidybreed/reference/dot-assemble_composite_tbv.md)
(`phenotype_components`) and
[`.build_tbv_env()`](https://austin-putz.github.io/tidybreed/reference/dot-build_tbv_env.md)
(`formula_tbv`) read through these helpers, and
[`.ap_materialize_tbvs()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_materialize_tbvs.md)
uses the same group lookup to decide whose TBVs to compute first.
Individual ids never appear in SQL text: every lookup joins a registered
view.

Group semantics (SGE / Bijma): a focal's group-mates are the *other*
individuals with the same value of `group_column` in `group_table`; the
aggregate is over the mates that have a TBV; a focal with no mates gets
`0`; a focal whose group value is `NULL` gets `NA` (a missing
component). `group_table` must have exactly one row per focal
individual.

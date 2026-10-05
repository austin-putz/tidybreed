# Contributor lookups for composite and `formula_tgv` phenotypes

The one place a contributor of a composite phenotype — the individual
itself, its dam or sire, or its group-mates — becomes a per-individual
value. Both
[`.assemble_composite_tgv()`](https://austin-putz.github.io/tidybreed/reference/dot-assemble_composite_tgv.md)
(`phenotype_components`) and
[`.build_tgv_env()`](https://austin-putz.github.io/tidybreed/reference/dot-build_tgv_env.md)
(`formula_tgv`) read through these helpers, and
[`.ap_materialize_tgvs()`](https://austin-putz.github.io/tidybreed/reference/dot-ap_materialize_tgvs.md)
uses the same group lookup to decide whose genetic values to compute
first. Every lookup reads `ind_tgv`: the total (`ind_tgv_total`) by
default, or the listed components. Individual ids never appear in SQL
text: every lookup joins a registered view.

Group semantics (SGE / Bijma): a focal's group-mates are the *other*
individuals with the same value of `group_column` in `group_table`; the
aggregate is over the mates that have a genetic value; a focal with no
mates gets `0`; a focal whose group value is `NULL` gets `NA` (a missing
component). `group_table` must have exactly one row per focal
individual.

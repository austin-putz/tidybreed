# Haplotype long-write + memory-batching helpers

Internal helpers shared by the three haplotype write paths
([`add_offspring()`](https://austin-putz.github.io/tidybreed/reference/add_offspring.md),
[`add_founders()`](https://austin-putz.github.io/tidybreed/reference/add_founders.md),
and
[`define_founder_haplotypes()`](https://austin-putz.github.io/tidybreed/reference/define_founder_haplotypes.md)'s
pool write). They convert the write path to direct **long** inserts and
bound peak memory with RAM-aware batching. No values change — only the
write shape.

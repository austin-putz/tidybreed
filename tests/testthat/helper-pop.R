# Shared test population builders — automatically sourced by testthat before
# any test file. Use these for minimal, fast in-memory populations.
#
# For domain-specific helpers (specific traits, phenotypes, chips, etc.) define
# local helpers inside the relevant test file.

#' Minimal in-memory population (genome defined, no founders)
make_pop_base <- function(pop_name = "t", n_loci = 100, n_chr = 2,
                          chr_len_Mb = 100) {
  open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = n_chr, chr_len_Mb = chr_len_Mb)
}

#' In-memory population with founder haplotypes and founders added
make_test_pop <- function(pop_name     = "t",
                          n_loci       = 100,
                          n_chr        = 2,
                          chr_len_Mb   = 100,
                          n_males      = 5,
                          n_females    = 5,
                          n_haplotypes = 50) {
  open_pop(pop_name = pop_name, db_name = ":memory:") |>
    define_genome(n_loci = n_loci, n_chr = n_chr, chr_len_Mb = chr_len_Mb) |>
    define_founder_haplotypes(n_haplotypes = n_haplotypes) |>
    add_founders(n_males = n_males, n_females = n_females, line_name = "A")
}

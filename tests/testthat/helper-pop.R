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
    get_table("founder_haplotypes") |>
    add_founders(n_males = n_males, n_females = n_females, line_name = "A")
}

#' Define a trait and store its population-wide additive target
#'
#' Test-only shorthand for `define_trait()` followed by
#' `define_effect_cov_matrix("additive", var, trait_name =)`: since 0.73.0
#' targets enter only through `trait_var_comp`'s writers. Neither call draws
#' from the RNG. Not a package function: `define_trait_simple()` was removed
#' on purpose (plans/import_qtl_effect_methods.md §6C).
with_additive_target <- function(pop, trait_name, var, ...) {
  pop <- define_trait(pop, trait_name, ...)
  suppressMessages(
    define_effect_cov_matrix(pop, "additive", var, trait_name = trait_name))
}

#' The breeding values: `ind_tgv` rows of the `additive` component
#'
#' Test-only shorthand for
#' `get_table(pop, "ind_tgv") |> filter(component_name == "additive")`,
#' collected, with columns `id_ind`, `trait_name`, `tgv_value`. Since 0.74.0
#' `ind_tgv` is the one table of genetic values (the old breeding-value table is gone);
#' for generated effects its `additive` component is the breeding value.
tgv_additive <- function(pop, trait_name = NULL) {
  out <- dplyr::collect(get_table(pop, "ind_tgv"))
  out <- out[out$component_name == "additive", , drop = FALSE]
  if (!is.null(trait_name)) out <- out[out$trait_name %in% trait_name, , drop = FALSE]
  out <- out[order(out$trait_name, out$id_ind), c("id_ind", "trait_name", "tgv_value")]
  tibble::as_tibble(out)
}

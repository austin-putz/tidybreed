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

#' Write additive terms with known values, centred like the generator
#'
#' Test-only stand-in for the manual `effects =` that `define_additive_effects()`
#' had before 0.74.1 (Q21: the generator now always samples and calibrates).
#' Writes one order-one `additive` term per selected locus through
#' `define_genome_effect_terms()`, with the generator's conventions: values in
#' ascending `locus_id` order, `center_value` = the base allele frequency
#' resolved exactly as the generator resolves it (`base_tbl = NULL` is the
#' founder pool of the line the effect applies to, with the Wahlund warning),
#' the generator's scope for (`line_name`, `parent_origin`), and
#' `mode = "replace_scope"`, so a re-run replaces only its own scope.
#'
#' The owner is `"custom"`: these values were not calibrated to any target,
#' which is exactly what the owner records. `add_tgv()` evaluates every owner,
#' so breeding-value oracles are unchanged.
#'
#' @param tbl `get_table(pop, "genome_meta")`, optionally filtered.
#' @param effects Numeric, one per selected locus, in `locus_id` order.
#' @return The `tidybreed_pop`, invisibly.
with_additive_terms <- function(tbl, trait_name, effects, line_name = NULL,
                                parent_origin = NULL, base_tbl = NULL,
                                effect_owner = "custom") {
  pop  <- tbl$pop
  conn <- pop$db_conn
  sel  <- dplyr::collect(tbl)$locus_name
  go   <- DBI::dbGetQuery(conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id")
  rows <- go$locus_name %in% sel
  if (length(effects) != sum(rows)) {
    stop("`effects` length (", length(effects), ") must equal the number of ",
         "selected loci (", sum(rows), ").", call. = FALSE)
  }
  base <- .dae_resolve_base(pop, base_tbl, line_name)
  .dae_require_base_at(base$p_base, rows, go$locus_name)
  terms <- data.frame(term_id       = go$locus_name[rows],
                      locus_name    = go$locus_name[rows],
                      contrast_name = "additive",
                      center_value  = base$p_base[rows],
                      genome_value  = as.numeric(effects),
                      stringsAsFactors = FALSE)
  .ge_write_terms(pop, trait_name, terms, effect_owner, "replace_scope",
                  .dae_scope(line_name, parent_origin), NULL, FALSE,
                  allow_reserved_owner = identical(effect_owner, "generated"))
  invisible(pop)
}

#' Plant uncalibrated terms under the reserved `"generated"` owner
#'
#' Test-only, through the internal writer (`allow_reserved_owner = TRUE`), as
#' the generator == writer test in `test-genome-effects-writer.R` does. Used
#' where a test needs generator-owned terms with chosen loci or values that no
#' calibrated call produces: per-trait QTL sets for `method = "union"`, or
#' fixtures for the owner rule. These terms are **not** calibrated, so they
#' break the package's "generated means calibrated" invariant on purpose.
plant_generated_additive <- function(tbl, trait_name, effects, ...) {
  with_additive_terms(tbl, trait_name, effects, ..., effect_owner = "generated")
}

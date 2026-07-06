#' Initialize a breeding population with genome (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `initialize_genome()` is deprecated as of tidybreed 0.40.0. Replace with
#' the new two-step API:
#'
#' ```r
#' pop <- open_pop(pop_name = "A", db_name = "A_tidybreed.duckdb") |>
#'   define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100)
#' ```
#'
#' @param pop_name Character scalar. Name of the population.
#' @param n_loci Integer scalar. Total number of loci to simulate.
#' @param n_chr Integer scalar. Number of chromosomes.
#' @param chr_len_Mb Numeric scalar or vector. Chromosome length(s) in Mb.
#' @param db_path Character scalar. Path to the DuckDB file, or `":memory:"`.
#'   Passed through to [open_pop()]'s `db_name` argument, so a relative,
#'   non-`":memory:"` value is placed inside [open_pop()]'s auto-created
#'   layer-2/3 folder structure (see `?open_pop`), not directly in the working
#'   directory. Default (when `NULL`) is `{pop_name}_tidybreed.duckdb`.
#' @param locus_names Character vector of length `n_loci` or `NULL`.
#' @param chr_names Character vector of length `n_chr` or `NULL`.
#' @param overwrite Logical. If `TRUE`, delete and recreate an existing database.
#'
#' @return A `tidybreed_pop` object.
#'
#' @seealso [open_pop()], [define_genome()]
#'
#' @export
initialize_genome <- function(pop_name,
                               n_loci,
                               n_chr,
                               chr_len_Mb,
                               db_path     = NULL,
                               locus_names = NULL,
                               chr_names   = NULL,
                               overwrite   = FALSE) {

  .Deprecated(
    new = "open_pop() |> define_genome()",
    msg = paste0(
      "'initialize_genome()' is deprecated as of tidybreed 0.40.0.\n",
      "Replace with:\n",
      "  open_pop(pop_name = \"", pop_name, "\") |>\n",
      "    define_genome(n_loci = ", n_loci,
      ", n_chr = ", n_chr,
      ", chr_len_Mb = ", chr_len_Mb[1], ")"
    )
  )

  if (is.null(db_path))
    db_path <- paste0(pop_name, "_tidybreed.duckdb")

  pop <- open_pop(
    pop_name = pop_name,
    db_name  = db_path,  # open_pop() handles absolute paths and :memory: directly
    tools    = NULL,
    clean    = overwrite
  )

  define_genome(
    pop,
    n_loci      = n_loci,
    n_chr       = n_chr,
    chr_len_Mb  = chr_len_Mb,
    locus_names = locus_names,
    chr_names   = chr_names
  )
}

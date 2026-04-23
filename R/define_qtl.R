#' Define QTL loci for a trait
#'
#' @description
#' Marks a subset of loci as QTL for a named trait by writing a BOOLEAN
#' column `is_QTL_{trait}` to `genome_meta`. Mirrors [define_chip()]: the same
#' six selection methods are supported (by count + `method`, by logical
#' vector, by locus IDs, by locus names).
#'
#' The trait must already exist in `trait_meta` (see [add_trait()]).
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character. Name of an existing trait in `trait_meta`.
#' @param n Integer. Number of QTL to select (with `method`). Mutually
#'   exclusive with `locus_tf`, `locus_ids`, `locus_names`.
#' @param method Character. Selection method when using `n`: `"random"`
#'   (default), `"even"`, or `"chromosome_even"`.
#' @param locus_tf Logical. TRUE/FALSE vector of length `n_loci` indicating
#'   QTL membership.
#' @param locus_ids Integer vector of locus IDs to mark as QTL.
#' @param locus_names Character vector of locus names to mark as QTL.
#' @param col_name Character. Column name in `genome_meta`. Default:
#'   `paste0("is_QTL_", trait_name)`.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_chip()], [set_qtl_effects()], [add_trait()]
#'
#' @examples
#' \dontrun{
#' pop <- pop |>
#'   add_trait("ADG", target_add_var = 0.25, residual_var = 0.75) |>
#'   define_qtl("ADG", n = 100, method = "chromosome_even")
#' }
#' @export
define_qtl <- function(pop,
                       trait_name,
                       n           = NULL,
                       method      = "random",
                       locus_tf    = NULL,
                       locus_ids   = NULL,
                       locus_names = NULL,
                       col_name    = paste0("is_QTL_", trait_name)) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(trait_name, what = "trait name")

  if (!"trait_meta" %in% pop$tables) {
    stop("No traits defined yet. Call add_trait() first.", call. = FALSE)
  }
  trait_exists <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM trait_meta WHERE trait_name = '",
           trait_name, "'")
  )$n
  if (trait_exists == 0) {
    stop("Trait '", trait_name, "' not found in trait_meta.", call. = FALSE)
  }

  methods_provided <- sum(
    !is.null(n),
    !is.null(locus_tf),
    !is.null(locus_ids),
    !is.null(locus_names)
  )
  if (methods_provided == 0) {
    stop(
      "Must specify one selection method: n, locus_tf, locus_ids, or locus_names",
      call. = FALSE
    )
  }
  if (methods_provided > 1) {
    stop(
      "Cannot specify multiple selection methods. ",
      "Choose one of: n, locus_tf, locus_ids, locus_names",
      call. = FALSE
    )
  }

  if (!is.null(n)) {
    valid_methods <- c("random", "even", "chromosome_even")
    if (!method %in% valid_methods) {
      stop(
        "Invalid method: '", method, "'. ",
        "Must be one of: ", paste(valid_methods, collapse = ", "),
        call. = FALSE
      )
    }
  }

  stopifnot(is.character(col_name), length(col_name) == 1, nchar(col_name) > 0)

  if (!"genome_meta" %in% pop$tables) {
    stop("genome_meta table does not exist. Call initialize_genome() first.",
         call. = FALSE)
  }

  genome <- dplyr::collect(get_table(pop, "genome_meta"))
  n_loci <- nrow(genome)
  if (n_loci == 0) {
    stop("genome_meta table is empty. Cannot define QTL.", call. = FALSE)
  }

  qtl_indicator <- rep(FALSE, n_loci)

  if (!is.null(n)) {
    qtl_indicator <- select_by_n(genome, n, method)
    message("Defined QTL for '", trait_name, "' with ", sum(qtl_indicator),
            " loci (method: ", method, ")")
  } else if (!is.null(locus_tf)) {
    if (!is.logical(locus_tf)) {
      stop("locus_tf must be a logical vector", call. = FALSE)
    }
    if (length(locus_tf) != n_loci) {
      stop("locus_tf length (", length(locus_tf),
           ") must equal number of loci (", n_loci, ")",
           call. = FALSE)
    }
    qtl_indicator <- locus_tf
    message("Defined QTL for '", trait_name, "' with ", sum(qtl_indicator),
            " loci (by locus_tf)")
  } else if (!is.null(locus_ids)) {
    qtl_indicator <- select_by_locus_ids(genome, locus_ids)
    message("Defined QTL for '", trait_name, "' with ", sum(qtl_indicator),
            " loci (by locus IDs)")
  } else if (!is.null(locus_names)) {
    qtl_indicator <- select_by_locus_names(genome, locus_names)
    message("Defined QTL for '", trait_name, "' with ", sum(qtl_indicator),
            " loci (by locus names)")
  }

  args   <- setNames(list(qtl_indicator), col_name)
  result <- do.call(mutate_table, c(list(tbl_obj = get_table(pop, "genome_meta")), args))

  invisible(result)
}

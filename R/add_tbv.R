#' Compute and store true breeding values without writing phenotypes
#'
#' @description
#' Computes the true breeding value (TBV) for each individual in the current
#' subset and each requested trait, and writes them to `ind_tbv`. Uses the
#' same math as the TBV step inside [add_phenotype()]:
#'
#' \preformatted{
#'   TBV_i = sum over QTL of add_{trait} * dose_i
#' }
#'
#' where `dose_i` is the 0/1/2 genotype for non-imprinted traits, or the 0/1
#' haplotype dose from the relevant parent for imprinted traits.
#'
#' The subset is taken from [filter.tidybreed_pop()] stashed predicates or
#' defaults to all individuals in `ind_meta`. The `expressed_sex` rule from
#' `trait_meta` is applied on top.
#'
#' Useful for tracking genetic trend across generations without collecting
#' phenotypes.
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait Character vector of trait name(s).
#' @param date_calc Date stored in `ind_tbv.date_calc`.
#'
#' @return The modified `tidybreed_pop` (invisibly) with `pending_filter`
#'   cleared.
#'
#' @seealso [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' pop |>
#'   dplyr::filter(gen == 2L) |>
#'   add_tbv(c("ADG", "BW"))
#' }
#' @export
add_tbv <- function(pop, trait, date_calc = Sys.Date()) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  stopifnot(is.character(trait), length(trait) >= 1)
  lapply(trait, validate_sql_identifier, what = "trait name")

  resolved <- resolve_pending_filter(pop)
  pop <- resolved$pop
  subset_ids <- resolved$ids

  ind_meta_subset <- if (is.null(subset_ids)) {
    dplyr::collect(get_table(pop, "ind_meta"))
  } else {
    get_table(pop, "ind_meta") |>
      dplyr::filter(.data$id_ind %in% !!subset_ids) |>
      dplyr::collect()
  }
  if (nrow(ind_meta_subset) == 0) {
    warning("No individuals matched; no TBVs computed.", call. = FALSE)
    return(invisible(pop))
  }

  meta_rows <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT trait_name, expressed_sex, expressed_parent ",
           "FROM trait_meta WHERE trait_name IN (",
           paste0("'", trait, "'", collapse = ", "), ")")
  )
  missing_t <- setdiff(trait, meta_rows$trait_name)
  if (length(missing_t) > 0) {
    stop("Traits not found: ", paste(missing_t, collapse = ", "),
         call. = FALSE)
  }
  meta_rows <- meta_rows[match(trait, meta_rows$trait_name), , drop = FALSE]

  genome_cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  for (t in trait) {
    if (!paste0("add_", t) %in% genome_cols) {
      stop("Additive-effect column 'add_", t, "' not found. Call ",
           "set_qtl_effects('", t, "', ...) first.", call. = FALSE)
    }
  }

  genome <- dplyr::collect(get_table(pop, "genome_meta"))
  geno_mat_full <- get_genotype_matrix(pop, subset_ids = subset_ids)

  for (t in trait) {
    m <- meta_rows[meta_rows$trait_name == t, ]
    if (m$expressed_sex != "both") {
      ids_t <- ind_meta_subset$id_ind[ind_meta_subset$sex == m$expressed_sex]
    } else {
      ids_t <- ind_meta_subset$id_ind
    }
    if (length(ids_t) == 0) next

    a <- genome[[paste0("add_", t)]]
    a[is.na(a)] <- 0

    if (m$expressed_parent == "both") {
      rows_idx <- match(ids_t, rownames(geno_mat_full))
      tbv <- as.numeric(geno_mat_full[rows_idx, , drop = FALSE] %*% a)
    } else {
      parent_origin <- if (m$expressed_parent == "parent_1") 1L else 2L
      hap_mat <- get_haplotype_matrix(pop, parent_origin, ids_t)
      tbv <- as.numeric(hap_mat %*% a)
    }

    tbv_df <- tibble::tibble(
      id_ind     = ids_t,
      trait_name = t,
      tbv        = tbv,
      date_calc  = as.Date(date_calc)
    )
    upsert_ind_tbv(pop, tbv_df)
    message("Computed TBV for ", length(ids_t), " individuals on trait '",
            t, "'.")
  }

  invisible(pop)
}

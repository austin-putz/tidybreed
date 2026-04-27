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
#' Pipe a `tidybreed_table` (from [get_table()] and optionally [filter()]) as
#' the first argument to select individuals. The `expressed_sex` rule from
#' `trait_meta` is applied on top.
#'
#' Useful for tracking genetic trend across generations without collecting
#' phenotypes.
#'
#' @param tbl A `tidybreed_table` object from [get_table()] (optionally piped
#'   through [filter()]). The table must contain an `id_ind` column.
#' @param trait Character vector of trait name(s).
#' @param date_calc Date stored in `ind_tbv.date_calc`.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 2L) |>
#'   add_tbv(c("ADG", "BW"))
#' }
#' @export
add_tbv <- function(tbl, trait, date_calc = Sys.Date(), ...) {

  stopifnot(inherits(tbl, "tidybreed_table"))
  pop <- tbl$pop
  validate_tidybreed_pop(pop)
  stopifnot(is.character(trait), length(trait) >= 1)
  lapply(trait, validate_sql_identifier, what = "trait name")

  extra_cols <- list(...)
  if (length(extra_cols) > 0) {
    for (nm in names(extra_cols)) {
      if (length(extra_cols[[nm]]) != 1L) {
        stop("Custom field '", nm, "' in add_tbv() must be a scalar ",
             "(broadcast to all records). Supply per-record vectors with ",
             "mutate_table() after the call.", call. = FALSE)
      }
    }
  }

  if (length(tbl$pending_filter) == 0) {
    subset_ids <- NULL
  } else {
    collected <- dplyr::collect(tbl)
    if (!"id_ind" %in% names(collected)) {
      stop("Filtered table '", tbl$table_name,
           "' must contain 'id_ind' to subset individuals for TBV computation.",
           call. = FALSE)
    }
    subset_ids <- unique(collected[["id_ind"]])
  }

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

    existing_ids <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT id_ind FROM ind_tbv WHERE trait_name = '", t,
             "' AND id_ind IN (", paste0("'", ids_t, "'", collapse = ", "), ")")
    )$id_ind
    if (length(existing_ids) > 0) {
      message(length(existing_ids), " individual(s) already have a TBV for '",
              t, "' and will be skipped.")
      ids_t <- setdiff(ids_t, existing_ids)
    }
    if (length(ids_t) == 0) next

    a <- genome[[paste0("add_", t)]]
    a[is.na(a)] <- 0

    p_base_col <- paste0("base_allele_freq_", t)
    p_base <- if (p_base_col %in% names(genome)) {
      pv <- as.numeric(genome[[p_base_col]])
      pv[is.na(pv)] <- 0
      pv
    } else {
      rep(0, length(a))
    }

    if (m$expressed_parent == "both") {
      rows_idx <- match(ids_t, rownames(geno_mat_full))
      tbv <- as.numeric(geno_mat_full[rows_idx, , drop = FALSE] %*% a) -
             2 * sum(p_base * a)
    } else {
      parent_origin <- if (m$expressed_parent == "parent_1") 1L else 2L
      hap_mat <- get_haplotype_matrix(pop, parent_origin, ids_t)
      tbv <- as.numeric(hap_mat %*% a) - sum(p_base * a)
    }

    tbv_df <- tibble::tibble(
      id_ind     = ids_t,
      trait_name = t,
      tbv        = tbv,
      date_calc  = as.Date(date_calc)
    )
    if (length(extra_cols) > 0) {
      prepped <- prepare_extra_cols(extra_cols, nrow(tbv_df), "ind_tbv",
                                   pop$db_conn)
      for (nm in names(prepped)) tbv_df[[nm]] <- prepped[[nm]]
    }
    upsert_ind_tbv(pop, tbv_df)
    message("Computed TBV for ", length(ids_t), " individuals on trait '",
            t, "'.")
  }

  invisible(pop)
}

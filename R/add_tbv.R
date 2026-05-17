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
#' Optionally computes true selection index values by multiplying per-trait TBVs
#' by weights from named indices defined with [define_index()], and writes them
#' to `ind_true_index`.
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
#' @param trait_name Character vector of trait name(s). When `NULL` (default),
#'   all traits currently in `trait_meta` are used (in `id_trait` order).
#' @param index_names Character vector of named index(es) from `index_meta` for
#'   which true index values should be computed from TBVs and written to
#'   `ind_true_index`. When `NULL` (default), no true index computation is
#'   performed. All index traits must be included in `trait_name` (or all traits
#'   when `trait_name = NULL`).
#' @param type Which weight column from `index_meta` to use: `"index"` uses
#'   `index_weight`, `"economic"` uses `economic_weight`, `"both"` computes and
#'   stores both (distinguished by the `weight_type` column in `ind_true_index`).
#'   Defaults to `"index"`.
#' @param overwrite_index Logical. When `FALSE` (default), individuals that
#'   already have a true index value in `ind_true_index` for the given
#'   `(index_name, weight_type)` combination are skipped — avoids redundant
#'   recomputation across generations. When `TRUE`, existing rows are deleted
#'   and recomputed (use when index weights have changed).
#' @param ... Optional extra columns written to `ind_tbv` (scalars only;
#'   broadcast to all records).
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [add_phenotype()], [define_index()], [add_index()]
#'
#' @examples
#' \dontrun{
#' # TBVs only
#' pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 2L) |>
#'   add_tbv(c("ADG", "BW"))
#'
#' # TBVs + true index values (index weights)
#' pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 2L) |>
#'   add_tbv(c("ADG", "BW"), index_names = "terminal", type = "both")
#' }
#' @export
add_tbv <- function(tbl, trait_name = NULL,
                    index_names    = NULL,
                    type           = c("index", "economic", "both"),
                    overwrite_index = FALSE,
                    ...) {

  stopifnot(inherits(tbl, "tidybreed_table"))
  pop <- tbl$pop
  validate_tidybreed_pop(pop)

  if (is.null(trait_name)) {
    trait_name <- DBI::dbGetQuery(
      pop$db_conn,
      "SELECT trait_name FROM trait_meta ORDER BY id_trait"
    )$trait_name
    if (length(trait_name) == 0L)
      stop("No traits found in trait_meta. Define traits with add_trait() first.",
           call. = FALSE)
  }

  stopifnot(is.character(trait_name), length(trait_name) >= 1)
  lapply(trait_name, validate_sql_identifier, what = "trait name")
  trait <- trait_name

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

  # Locus order needed to map genome_effects locus_name to genotype matrix columns
  locus_order <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id"
  )
  n_loci <- nrow(locus_order)

  geno_mat_full <- get_genotype_matrix(pop, subset_ids = subset_ids)

  for (t in trait) {
    m <- meta_rows[meta_rows$trait_name == t, ]
    if (m$expressed_sex != "both") {
      ids_t <- ind_meta_subset$id_ind[ind_meta_subset$sex == m$expressed_sex]
    } else {
      ids_t <- ind_meta_subset$id_ind
    }
    if (length(ids_t) == 0) next

    # Load additive effects from genome_effects (population-wide: line_name IS NULL)
    effects_df <- DBI::dbGetQuery(
      pop$db_conn,
      paste0(
        "SELECT locus_name, genome_value, base_allele_freq ",
        "FROM genome_effects ",
        "WHERE trait_name = '", t, "' ",
        "AND genome_effect_type = 'additive' ",
        "AND line_name IS NULL"
      )
    )
    if (nrow(effects_df) == 0L) {
      stop(
        "No additive effects found for trait '", t, "' in genome_effects. ",
        "Call add_additive_effects() first.",
        call. = FALSE
      )
    }

    # Map locus_name → position in the genotype matrix (locus_id order)
    a      <- rep(0, n_loci)
    p_base <- rep(0, n_loci)
    idx    <- match(effects_df$locus_name, locus_order$locus_name)
    a[idx]      <- effects_df$genome_value
    baf <- effects_df$base_allele_freq
    p_base[idx] <- ifelse(is.na(baf), 0, baf)

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
      tbv_value  = tbv
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

  # --- True index computation from TBVs ---
  if (!is.null(index_names)) {
    stopifnot(is.character(index_names), length(index_names) >= 1)
    lapply(index_names, validate_sql_identifier, what = "index name")
    type        <- match.arg(type)
    weight_types <- switch(type,
      "index"    = "index",
      "economic" = "economic",
      "both"     = c("index", "economic")
    )
    all_subset_ids <- ind_meta_subset$id_ind

    for (idx_name in index_names) {
      idx_check <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT COUNT(*) AS n FROM index_meta WHERE index_name = '",
               idx_name, "'")
      )$n
      if (idx_check == 0L)
        stop("Index '", idx_name, "' not found in index_meta. ",
             "Define it with define_index() first.", call. = FALSE)

      idx_meta <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT trait_name, index_weight, economic_weight ",
               "FROM index_meta WHERE index_name = '", idx_name, "' ",
               "ORDER BY trait_name")
      )
      idx_traits <- idx_meta$trait_name

      for (wt in weight_types) {
        wt_col  <- if (wt == "index") "index_weight" else "economic_weight"
        wt_vals <- idx_meta[[wt_col]]

        if (any(is.na(wt_vals)))
          stop("Some ", wt_col, " values are NA for index '", idx_name, "' ",
               "trait(s): ", paste(idx_traits[is.na(wt_vals)], collapse = ", "),
               ". Supply them via define_index().", call. = FALSE)

        wt_vec <- setNames(as.numeric(wt_vals), idx_traits)

        # Skip individuals that already have a value when overwrite_index = FALSE
        target_ids <- all_subset_ids
        if (!overwrite_index) {
          id_list <- paste0("'", target_ids, "'", collapse = ", ")
          existing_ids <- DBI::dbGetQuery(
            pop$db_conn,
            paste0("SELECT DISTINCT id_ind FROM ind_true_index ",
                   "WHERE id_ind IN (", id_list, ") ",
                   "AND index_name = '", idx_name, "' ",
                   "AND weight_type = '", wt, "'")
          )$id_ind
          target_ids <- setdiff(target_ids, existing_ids)
        }

        if (length(target_ids) == 0L) {
          message("True index '", idx_name, "' (", wt, ") already exists for ",
                  "all individuals; skipping (overwrite_index = FALSE).")
          next
        }

        # Read TBVs from ind_tbv for target individuals and index traits
        id_list    <- paste0("'", target_ids, "'", collapse = ", ")
        trait_list <- paste0("'", idx_traits, "'", collapse = ", ")
        tbv_data <- DBI::dbGetQuery(
          pop$db_conn,
          paste0("SELECT id_ind, trait_name, tbv_value FROM ind_tbv ",
                 "WHERE id_ind IN (", id_list, ") ",
                 "AND trait_name IN (", trait_list, ")")
        )

        # Build n_ind × n_traits matrix in consistent column order
        ind_order <- sort(unique(tbv_data$id_ind))
        tbv_mat <- matrix(
          unlist(lapply(idx_traits, function(t) {
            sub <- tbv_data[tbv_data$trait_name == t, , drop = FALSE]
            sub$tbv_value[match(ind_order, sub$id_ind)]
          })),
          nrow = length(ind_order),
          ncol = length(idx_traits),
          dimnames = list(ind_order, idx_traits)
        )

        if (any(is.na(tbv_mat))) {
          miss_idx    <- which(is.na(tbv_mat), arr.ind = TRUE)
          miss_traits <- unique(idx_traits[miss_idx[, 2L]])
          stop("TBVs missing for index '", idx_name, "' trait(s): ",
               paste(miss_traits, collapse = ", "),
               ". Ensure those traits are included in the trait_name argument ",
               "of this add_tbv() call.", call. = FALSE)
        }

        true_index_values <- as.numeric(tbv_mat %*% wt_vec)

        # Delete existing rows before reinsert when overwrite_index = TRUE
        if (overwrite_index) {
          id_list2 <- paste0("'", target_ids, "'", collapse = ", ")
          DBI::dbExecute(
            pop$db_conn,
            paste0("DELETE FROM ind_true_index ",
                   "WHERE id_ind IN (", id_list2, ") ",
                   "AND index_name = '", idx_name, "' ",
                   "AND weight_type = '", wt, "'")
          )
        }

        ti_df <- tibble::tibble(
          id_ind           = ind_order,
          index_name       = idx_name,
          weight_type      = wt,
          true_index_value = true_index_values
        )
        upsert_ind_true_index(pop, ti_df)
        message("Computed true index '", idx_name, "' (", wt, ") for ",
                nrow(ti_df), " individuals.")
      }
    }
  }

  invisible(pop)
}


upsert_ind_true_index <- function(pop, df) {
  if (nrow(df) == 0L) return(invisible(NULL))
  start <- next_int_id(pop$db_conn, "ind_true_index", "id_true_index")
  df <- tibble::add_column(
    df,
    id_true_index = seq.int(start, start + nrow(df) - 1L),
    .before = 1
  )
  tmp <- paste0("_ti_tmp_", as.character(round(as.numeric(Sys.time()) * 1000)))
  duckdb::duckdb_register(pop$db_conn, tmp, as.data.frame(df))
  on.exit(duckdb::duckdb_unregister(pop$db_conn, tmp), add = TRUE)
  DBI::dbExecute(pop$db_conn,
    paste0("INSERT INTO ind_true_index SELECT * FROM ", tmp))
  invisible(NULL)
}

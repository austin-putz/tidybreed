#' Compute and store true breeding values without writing phenotypes
#'
#' @description
#' Computes the true breeding value (TBV) for each individual in the current
#' subset and each requested trait, and writes them to `ind_tbv`. This is the
#' exact function [add_phenotype()] calls internally (once for every source
#' trait it needs) before assembling phenotype records — there is no separate
#' "TBV math" duplicated elsewhere.
#'
#' TBV is the Falconer-centered sum, across every `ind_haplotype` row (one per
#' allele copy, not genotype dosage) for the individual, of:
#'
#' \preformatted{
#'   TBV_i = sum over haplotype rows of (allele - base_allele_freq) * genome_value
#' }
#'
#' `genome_value` and `base_allele_freq` are read from `genome_effects`
#' (`genome_effect_type = "additive"`). For each haplotype row, a
#' **line-specific** effect (`genome_effects.line_name` matching that row's
#' `line_origin`) is preferred; the **population-wide** effect
#' (`genome_effects.line_name IS NULL`) is used only when no line-specific row
#' exists for that locus/line. This per-locus fallback is what makes
#' crossbreeding TBV correct — e.g. a Duroc x Landrace F1 is centered against
#' each parent line's own QTL effects and base allele frequency (see the
#' "Crossbreeding TBV" example below). For **imprinted** traits
#' (`trait_meta.expressed_parent` = `"parent_1"` or `"parent_2"`), only
#' haplotype rows from that parent's `parent_origin` are summed before the
#' same line-matching logic applies.
#'
#' Optionally computes true selection index values by multiplying per-trait TBVs
#' by weights from named indices defined with [define_index()], and writes them
#' to `ind_true_index`.
#'
#' Pipe a `tidybreed_table` (from [get_table()] and optionally [dplyr::filter()]) as
#' the first argument to select individuals. The `expressed_sex` rule from
#' `trait_meta` is applied on top.
#'
#' Useful for tracking genetic trend across generations without collecting
#' phenotypes.
#'
#' @param tbl A `tidybreed_table` object from [get_table()] (optionally piped
#'   through [dplyr::filter()]). The table must contain an `id_ind` column.
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
#' # Crossbreeding TBV: line-specific additive effects for two pure lines, then
#' # a Duroc x Landrace F1 centered against each parent line's own effects and
#' # base allele frequency (see the Description above for the matching rule)
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   define_additive_effects("ADG", effects = duroc_effects, line_name = "Duroc")
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   define_additive_effects("ADG", effects = landrace_effects, line_name = "Landrace")
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(line_name == "F1") |>
#'   add_tbv("ADG")
#'
#' # TBVs only, for a generation subset
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 2L) |>
#'   add_tbv(c("ADG", "BW"))
#'
#' # TBVs + true index values (both index and economic weights) written to
#' # ind_true_index, distinguished by weight_type
#' pop <- pop |>
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
      stop("No traits found in trait_meta. ",
           "Define traits with define_trait() first.", call. = FALSE)
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
    paste0("SELECT trait_name, expressed_parent ",
           "FROM trait_meta WHERE trait_name IN (",
           paste0("'", trait, "'", collapse = ", "), ")")
  )
  missing_t <- setdiff(trait, meta_rows$trait_name)
  if (length(missing_t) > 0) {
    stop("Traits not found: ", paste(missing_t, collapse = ", "),
         call. = FALSE)
  }
  meta_rows <- meta_rows[match(trait, meta_rows$trait_name), , drop = FALSE]

  for (t in trait) {
    m     <- meta_rows[meta_rows$trait_name == t, ]
    ids_t <- ind_meta_subset$id_ind
    if (length(ids_t) == 0) next

    effect_count <- DBI::dbGetQuery(
      pop$db_conn,
      paste0(
        "SELECT COUNT(*) AS n FROM genome_effects ",
        "WHERE trait_name = '", t, "' ",
        "AND genome_effect_type = 'additive'"
      )
    )$n
    if (effect_count == 0L) {
      stop(
        "No additive effects found for trait '", t, "' in genome_effects. ",
        "Call define_additive_effects() first.",
        call. = FALSE
      )
    }

    id_list <- paste0("'", ids_t, "'", collapse = ", ")
    parent_filter <- if (m$expressed_parent == "both") {
      ""
    } else {
      parent_origin <- if (m$expressed_parent == "parent_1") 1L else 2L
      paste0("AND h.parent_origin = ", parent_origin, " ")
    }

    # Centered TBV, folding (allele - base_allele_freq) into the summed term:
    # this generalizes Falconer centering to per-line base_allele_freq without
    # a separate centering constant. Line-specific genome_effects rows
    # (matched on line_origin) take precedence; a population-wide row
    # (line_name IS NULL) is only used when no line-specific row exists for
    # that locus/line (per-locus fallback, via NOT EXISTS). COALESCE(...,0)
    # on base_allele_freq mirrors the prior R-side `ifelse(is.na(baf), 0, baf)`
    # -- a NULL base_allele_freq must not null out the whole summed term.
    tbv_sql <- DBI::dbGetQuery(
      pop$db_conn,
      paste0(
        "SELECT h.id_ind, ",
        "SUM((h.allele - COALESCE(e.base_allele_freq, 0)) * e.genome_value) AS tbv_value ",
        "FROM ind_haplotype h ",
        "JOIN genome_effects e ",
        "  ON h.locus_name = e.locus_name ",
        " AND e.trait_name = '", t, "' ",
        " AND e.genome_effect_type = 'additive' ",
        " AND ( e.line_name = h.line_origin ",
        "       OR (e.line_name IS NULL AND NOT EXISTS ( ",
        "             SELECT 1 FROM genome_effects e2 ",
        "             WHERE e2.locus_name = h.locus_name ",
        "               AND e2.trait_name = '", t, "' ",
        "               AND e2.genome_effect_type = 'additive' ",
        "               AND e2.line_name = h.line_origin)) ) ",
        "WHERE h.id_ind IN (", id_list, ") ",
        parent_filter,
        "GROUP BY h.id_ind"
      )
    )

    tbv <- tbv_sql$tbv_value[match(ids_t, tbv_sql$id_ind)]
    if (anyNA(tbv)) {
      missing_ids <- ids_t[is.na(tbv)]
      stop(
        "No haplotype rows matched additive effects for trait '", t,
        "' for individual(s): ",
        paste(utils::head(missing_ids, 5), collapse = ", "),
        if (length(missing_ids) > 5) ", ..." else "",
        ". This indicates missing ind_haplotype rows (v1 storage is fully ",
        "dense, so every individual should have a row at every QTL locus).",
        call. = FALSE
      )
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

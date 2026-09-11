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
#'   TBV_i = sum over allele copies of (allele - center_value) * genome_value
#' }
#'
#' `genome_value` and `center_value` come from the **order-one `additive`
#' terms** under the reserved effect owner `generated_additive_tbv`, the terms
#' [define_additive_effects()] writes. `center_value` is that variant's base
#' allele frequency.
#'
#' **This is one filtered call into the same evaluator [add_tgv()] uses**, not a
#' second implementation: `add_tbv()` is `add_tgv()` restricted to the reserved
#' owner and to single-member additive terms.
#'
#' The filter is deliberate and not merely conservative. Under functional
#' \eqn{(a, d)} input the stored coefficient is \eqn{a}, while the
#' breeding-value coefficient in a diploid HWE base is
#' \eqn{\alpha = a + d(q - p)}; under epistasis, average effects depend on other
#' loci and on LD. So arbitrary terms written through [define_genome_effects()]
#' contribute to `ind_tgv` but **never silently redefine the breeding value**,
#' additive members appearing inside interactions are ignored, and `ind_tbv`
#' keeps its exact meaning. Deriving average effects from a general
#' non-additive model is a separate calculation.
#'
#' Each allele copy takes the **most specific** variant whose origin predicate
#' matches its `(line_origin, parent_origin)` label, falling back per copy to
#' the common variant. This per-copy fallback is what makes crossbreeding TBV
#' correct — e.g. a Duroc x Landrace F1 is centered against each parent line's
#' own effects and base allele frequency (see the "Crossbreeding TBV" example
#' below). **Imprinting** is a property of the effect, not of the trait: a term
#' scoped to one `parent_origin` (see [define_additive_effects()]) reads only
#' that parent's allele copies, per locus and per line.
#'
#' Optionally computes true selection index values by multiplying per-trait TBVs
#' by weights from named indices defined with [define_index()], and writes them
#' to `ind_true_index`.
#'
#' Pipe a `tidybreed_table` (from [get_table()] and optionally [dplyr::filter()]) as
#' the first argument to select individuals. Every individual in that subset
#' receives a TBV for every requested trait — unlike [add_phenotype()], no
#' sex-expression rule is applied here (`expressed_sex` is an observation-layer
#' property of `phenotype_meta`, not of a genetic component trait).
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
#' @seealso [add_tgv()] for every component of the genetic value,
#'   [add_phenotype()], [define_index()], [add_index()]
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
  conn <- pop$db_conn

  trait <- .gev_resolve_traits(conn, trait_name)

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

  ids_t <- .gev_subset_ids(tbl, "TBV")
  if (length(ids_t) == 0) {
    warning("No individuals matched; no TBVs computed.", call. = FALSE)
    return(invisible(pop))
  }

  model <- .gev_read_model(conn, trait, GE_ADDITIVE_OWNER,
                           order1_additive_only = TRUE)
  for (t in trait) .gev_require_terms(model, t, GE_ADDITIVE_OWNER, TRUE)

  res <- .gev_evaluate(conn, ids_t, trait, GE_ADDITIVE_OWNER,
                       order1_additive_only = TRUE, model = model)

  for (t in trait) {
    sub <- res[res$trait_name == t, , drop = FALSE]
    .gev_require_contribution(ids_t, sub$id_ind, t)

    tbv_df <- tibble::tibble(
      id_ind     = ids_t,
      trait_name = t,
      tbv_value  = sub$tgv_value[match(ids_t, sub$id_ind)]
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
    all_subset_ids <- ids_t

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

        # An individual with no ind_tbv row for *any* index trait never reaches
        # tbv_mat at all, so the anyNA() check below cannot see it -- it would
        # be dropped silently (zero rows written, no error). Catch that here.
        dropped <- setdiff(target_ids, ind_order)
        if (length(dropped) > 0) {
          stop("No TBVs found for index '", idx_name, "' trait(s) (",
               paste(idx_traits, collapse = ", "), ") for individual(s): ",
               paste(utils::head(sort(dropped), 5), collapse = ", "),
               if (length(dropped) > 5) ", ..." else "",
               ". Ensure those traits are included in the trait_name argument ",
               "of this add_tbv() call.", call. = FALSE)
        }

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

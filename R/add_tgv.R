#' Compute and store true genetic values
#'
#' @description
#' Evaluates the stored effect model for each individual in the current subset
#' and writes the result to `ind_tgv`, one row per
#' (individual x trait x **component**). The total is the derived view
#' `ind_tgv_total`, never a stored row — a stored total would make every
#' `SUM(tgv_value)` double-count. This is the one table of true genetic values:
#' [add_phenotype()] calls this function for every trait a phenotype reads, and
#' reads the total (or the components `phenotype_components.component_names`
#' lists) from it.
#'
#' Every term of the trait is evaluated: additive, dominance, hand-entered
#' genotype surfaces and multi-locus interactions, under every effect owner.
#' `component_name` records how each term was declared:
#'
#' | `component_name` | Terms it collects |
#' |---|---|
#' | `"additive"` | one-locus terms whose contrast is `additive` |
#' | `"dominance"` | one-locus `dominance` terms |
#' | `"indicator"` | one-locus `indicator` terms (a hand-entered genotype surface) |
#' | `"interaction"` | any term over two or more loci |
#'
#' **The breeding value is `component_name = "additive"`** for effects written
#' in statistical coding at a single base allele frequency per locus, which is
#' what the generators ([define_additive_effects()]) write. For hand-written
#' functional terms ([ad_terms()] with `coding = "functional"`, or an
#' `indicator` surface) `additive` is the functional additive effect, not the
#' breeding value: under functional coding the average effect is
#' \eqn{\alpha = a + d(q - p)}, and under epistasis it depends on other loci
#' and on LD.
#'
#' `ind_tgv` stores the **raw sum of the stored terms — no mean is added.** A
#' pure Cockerham model yields centered deviations; a raw `indicator` surface
#' yields absolute genotypic values with a non-zero mean by construction, and is
#' not silently re-centered. See [ad_terms()], which reports the implied genetic
#' mean and writes it nowhere.
#'
#' Each allele copy takes the **most specific** variant whose origin predicate
#' matches its `(line_origin, parent_origin)` label, falling back per copy to
#' the common variant. This per-copy fallback is what makes crossbred genetic
#' values correct — e.g. a Duroc x Landrace F1 is centered against each parent
#' line's own effects and base allele frequency. **Imprinting** is a property of
#' the effect, not of the trait: a term scoped to one `parent_origin` reads
#' only that parent's allele copies.
#'
#' Re-evaluating replaces an individual's rows for a trait: a component the
#' model no longer has is deleted, so no stale component survives, and a
#' component that is still there is updated in place, keeping any custom
#' columns added to it.
#'
#' Optionally computes true selection index values by multiplying per-trait
#' genetic values (`component_name`, default the breeding value) by weights
#' from named indices defined with [define_index()], and writes them to
#' `ind_true_index`.
#'
#' Every individual in the subset receives a value for every requested trait —
#' unlike [add_phenotype()], no sex-expression rule is applied here
#' (`expressed_sex` is a property of `phenotype_meta`, not of a genetic trait).
#'
#' @param tbl A `tidybreed_table` from [get_table()], optionally piped through
#'   [dplyr::filter()]. Any table with an `id_ind` column is accepted; the
#'   individuals acted on are the distinct `id_ind` values present in the
#'   (filtered) table. An unfiltered `ind_meta` selects every individual; an
#'   unfiltered `ind_ebv`, `ind_index`, `ind_genotype`, ... selects only the
#'   individuals that have rows there. A table without `id_ind` is an error.
#' @param trait_name Character vector of trait name(s). When `NULL` (default),
#'   all traits in `trait_meta` are used (in `id_trait` order).
#' @param index_names Character vector of named index(es) from `index_meta` for
#'   which true index values are computed and written to `ind_true_index`. When
#'   `NULL` (default), no true index is computed. Every index trait must have an
#'   `ind_tgv` row for every individual of the subset (computed by this call or
#'   an earlier one).
#' @param weight_type Which weight column from `index_meta` to use: `"index"`
#'   uses `index_weight`, `"economic"` uses `economic_weight`, `"both"` computes
#'   and stores both (distinguished by `ind_true_index.weight_type`). Defaults to
#'   `"index"`.
#' @param component_name Which genetic value the true index weights: one of
#'   `"additive"` (default — selection indices are on breeding values),
#'   `"dominance"`, `"indicator"`, `"interaction"`, or `"total"` (the
#'   `ind_tgv_total` view). A component the trait's model has no terms for
#'   contributes 0. Stored in `ind_true_index.component_name`.
#' @param overwrite_index Logical. When `FALSE` (default), individuals that
#'   already have a true index value for the given
#'   `(index_name, weight_type, component_name)` are skipped. When `TRUE`,
#'   existing rows are deleted and recomputed (use when index weights have
#'   changed).
#' @param ... Optional extra columns written to `ind_tgv` (scalars only;
#'   broadcast to every row written).
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @section Resource guard:
#' Evaluation enumerates **label-vectors** — one distinct
#' `(line_origin, parent_origin)` label per member — not allele copies per
#' individual, so its cost is a function of the model rather than of population
#' size. The count is estimated per fallback family before any work runs; it
#' warns above `getOption("tidybreed.label_vector_warn", 1e4)` and stops above
#' `getOption("tidybreed.label_vector_max", 1e6)`, so an accidental high-order
#' scoped term fails loudly instead of appearing to hang.
#'
#' @seealso [define_genome_effect_terms()] and [define_additive_effects()] for
#'   writing the terms this evaluates, [add_phenotype()], [define_index()],
#'   [add_index()].
#'
#' @examples
#' \dontrun{
#' # Every component of the model, for one generation
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 2L) |>
#'   add_tgv("ADG")
#'
#' # The breeding values, and the derived total
#' pop |> get_table("ind_tgv") |>
#'   dplyr::filter(component_name == "additive") |> dplyr::collect()
#' pop |> get_table("ind_tgv_total") |> dplyr::collect()
#'
#' # Genetic values + true index values (index and economic weights) on the
#' # breeding values, written to ind_true_index
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 2L) |>
#'   add_tgv(c("ADG", "BW"), index_names = "terminal", weight_type = "both")
#' }
#' @export
add_tgv <- function(tbl, trait_name = NULL,
                    index_names     = NULL,
                    weight_type     = c("index", "economic", "both"),
                    component_name  = "additive",
                    overwrite_index = FALSE,
                    ...) {

  stopifnot(inherits(tbl, "tidybreed_table"))
  pop <- tbl$pop
  validate_tidybreed_pop(pop)
  conn <- pop$db_conn

  traits <- .gev_resolve_traits(conn, trait_name)

  extra_cols <- list(...)
  if (length(extra_cols) > 0L) {
    if (is.null(names(extra_cols)) || any(names(extra_cols) == "")) {
      stop("Custom fields in add_tgv() must be named.", call. = FALSE)
    }
    for (nm in names(extra_cols)) {
      if (length(extra_cols[[nm]]) != 1L) {
        stop("Custom field '", nm, "' in add_tgv() must be a scalar ",
             "(broadcast to all records). Supply per-record vectors with ",
             "mutate_table() after the call.", call. = FALSE)
      }
    }
  }

  # Validate the index request before anything is written.
  if (!is.null(index_names)) {
    weight_type <- match.arg(weight_type)
    .tgv_check_component(component_name, "component_name")
    if (!is.logical(overwrite_index) || length(overwrite_index) != 1L ||
        is.na(overwrite_index)) {
      stop("`overwrite_index` must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.character(index_names) || length(index_names) == 0L ||
        anyNA(index_names)) {
      stop("`index_names` must be a character vector of index names.",
           call. = FALSE)
    }
    for (idx in index_names) validate_sql_identifier(idx, what = "index name")
    known <- DBI::dbGetQuery(conn, paste0(
      "SELECT DISTINCT index_name FROM index_meta WHERE index_name IN (",
      sql_in_list(index_names, what = "index name"), ")"))$index_name
    unknown <- setdiff(index_names, known)
    if (length(unknown) > 0L) {
      stop("Index '", unknown[1], "' not found in index_meta. ",
           "Define it with define_index() first.", call. = FALSE)
    }
  }

  ids <- resolve_subset_ids(tbl, "TGV computation", all_if_null = TRUE)
  if (length(ids) == 0L) {
    warning("No individuals matched; no TGVs computed.", call. = FALSE)
    return(invisible(pop))
  }

  model <- .gev_read_model(conn, traits)
  for (t in traits) .gev_require_terms(model, t)

  res <- .gev_evaluate(conn, ids, traits, model = model)
  for (t in traits) {
    .gev_require_contribution(ids, res$id_ind[res$trait_name == t], t)
  }

  extras <- prepare_extra_cols(extra_cols, nrow(res), "ind_tgv", conn)
  .gev_write_tgv(conn, res, traits, ids, extras)
  for (t in traits) {
    message("Computed TGV for ", length(ids), " individuals on trait '", t,
            "' (", paste(sort(unique(res$component_name[res$trait_name == t])),
                         collapse = ", "), ").")
  }

  if (!is.null(index_names)) {
    weight_types <- switch(weight_type,
                           index    = "index",
                           economic = "economic",
                           both     = c("index", "economic"))
    for (idx in index_names) {
      .tgv_true_index(conn, ids, idx, weight_types, component_name,
                      overwrite_index)
    }
  }

  invisible(pop)
}


#' Validate a genetic-value component selector
#'
#' One of [TGV_COMPONENT_NAMES], or `"total"` (the `ind_tgv_total` view).
#' @keywords internal
#' @noRd
.tgv_check_component <- function(x, arg) {
  ok <- c(TGV_COMPONENT_NAMES, "total")
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% ok) {
    stop("`", arg, "` must be one of ", paste0("\"", ok, "\"", collapse = ", "),
         ".", call. = FALSE)
  }
  invisible(x)
}


#' Genetic values of individuals for a set of traits, by component
#'
#' The one reader of `ind_tgv` for consumers (phenotypes, true indices).
#' `components = "total"` reads the `ind_tgv_total` view; otherwise the listed
#' components are summed in `component_name` order (deterministic, as the
#' view), and a listed component the individual has no row for contributes 0. An individual with **no** `ind_tgv`
#' row for a trait is absent from the result — the caller decides whether that
#' is missing or an error. Ids are registered as a view, never written into SQL.
#'
#' @param conn A DBI connection.
#' @param ids Character vector of individuals.
#' @param trait_names Character vector of traits.
#' @param components `"total"`, or a character vector of
#'   [TGV_COMPONENT_NAMES].
#' @return Data frame `id_ind`, `trait_name`, `value`.
#' @keywords internal
#' @noRd
.tgv_read <- function(conn, ids, trait_names, components = "total") {
  empty <- data.frame(id_ind = character(0), trait_name = character(0),
                      value = numeric(0), stringsAsFactors = FALSE)
  ids <- unique(ids[!is.na(ids)])
  if (length(ids) == 0L) return(empty)

  tmp <- paste0("__tgv_ids_", as.character(round(as.numeric(Sys.time()) * 1000)))
  duckdb::duckdb_register(conn, tmp, data.frame(id_ind = ids,
                                                stringsAsFactors = FALSE))
  on.exit(try(duckdb::duckdb_unregister(conn, tmp), silent = TRUE), add = TRUE)

  tr <- sql_in_list(trait_names, what = "trait name")
  # The selected individuals are joined first and aggregated after, so a call
  # never aggregates the whole of ind_tgv. "total" is the same ordered sum the
  # ind_tgv_total view computes, so the two agree bit for bit.
  value_sql <- if (identical(components, "total")) {
    "list_sum(list(t.tgv_value ORDER BY t.component_name))"
  } else {
    bad <- setdiff(components, TGV_COMPONENT_NAMES)
    if (length(bad) > 0L) {
      stop("Internal error: unknown ind_tgv component(s) ",
           paste0("'", bad, "'", collapse = ", "), ".", call. = FALSE)
    }
    .tgv_component_sum_sql(components, "t.")
  }
  sql <- paste0("SELECT t.id_ind, t.trait_name, ", value_sql, " AS value ",
                "FROM ind_tgv AS t JOIN ", tmp, " AS f USING (id_ind) ",
                "WHERE t.trait_name IN (", tr, ") ",
                "GROUP BY t.id_ind, t.trait_name")
  out <- DBI::dbGetQuery(conn, sql)
  out[order(out$trait_name, out$id_ind), , drop = FALSE]
}


#' SQL for the sum of listed `ind_tgv` components of one group
#'
#' Added in `component_name` order, as the `ind_tgv_total` view does, so the
#' result does not depend on DuckDB's thread count; an unlisted component adds
#' an exact 0, so one listed component is returned bit for bit.
#' @keywords internal
#' @noRd
.tgv_component_sum_sql <- function(components, prefix = "") {
  paste0("list_sum(list(CASE WHEN ", prefix, "component_name IN (",
         sql_in_list(components, what = "component name"), ") THEN ", prefix,
         "tgv_value ELSE 0 END ORDER BY ", prefix, "component_name))")
}


#' Compute and write one index's true index values
#'
#' @keywords internal
#' @noRd
.tgv_true_index <- function(conn, ids, idx_name, weight_types, component_name,
                            overwrite_index) {
  idx_lit  <- DBI::dbQuoteLiteral(conn, idx_name)
  comp_lit <- DBI::dbQuoteLiteral(conn, component_name)

  idx_meta <- DBI::dbGetQuery(conn, paste0(
    "SELECT trait_name, index_weight, economic_weight FROM index_meta ",
    "WHERE index_name = ", idx_lit, " ORDER BY trait_name"))
  if (nrow(idx_meta) == 0L) {
    stop("Index '", idx_name, "' not found in index_meta. ",
         "Define it with define_index() first.", call. = FALSE)
  }
  idx_traits <- idx_meta$trait_name

  tmp <- paste0("__ti_ids_", as.character(round(as.numeric(Sys.time()) * 1000)))
  on.exit(try(duckdb::duckdb_unregister(conn, tmp), silent = TRUE), add = TRUE)

  for (wt in weight_types) {
    wt_col  <- if (wt == "index") "index_weight" else "economic_weight"
    wt_vals <- idx_meta[[wt_col]]
    if (anyNA(wt_vals)) {
      stop("Some ", wt_col, " values are NA for index '", idx_name, "' ",
           "trait(s): ", paste(idx_traits[is.na(wt_vals)], collapse = ", "),
           ". Supply them via define_index().", call. = FALSE)
    }
    wt_vec  <- stats::setNames(as.numeric(wt_vals), idx_traits)
    key_sql <- paste0("index_name = ", idx_lit, " AND weight_type = ",
                      DBI::dbQuoteLiteral(conn, wt), " AND component_name = ",
                      comp_lit)

    target_ids <- ids
    if (!overwrite_index) {
      duckdb::duckdb_register(conn, tmp, data.frame(id_ind = target_ids,
                                                    stringsAsFactors = FALSE))
      have <- DBI::dbGetQuery(conn, paste0(
        "SELECT DISTINCT t.id_ind FROM ind_true_index AS t JOIN ", tmp,
        " AS f USING (id_ind) WHERE ", key_sql))$id_ind
      duckdb::duckdb_unregister(conn, tmp)
      target_ids <- setdiff(target_ids, have)
    }
    if (length(target_ids) == 0L) {
      message("True index '", idx_name, "' (", wt, ", ", component_name,
              ") already exists for all individuals; skipping ",
              "(overwrite_index = FALSE).")
      next
    }

    vals      <- .tgv_read(conn, target_ids, idx_traits, component_name)
    ind_order <- sort(unique(target_ids))
    val_mat <- matrix(
      unlist(lapply(idx_traits, function(t) {
        sub <- vals[vals$trait_name == t, , drop = FALSE]
        sub$value[match(ind_order, sub$id_ind)]
      })),
      nrow = length(ind_order), ncol = length(idx_traits),
      dimnames = list(ind_order, idx_traits))

    if (anyNA(val_mat)) {
      miss_idx    <- which(is.na(val_mat), arr.ind = TRUE)
      miss_traits <- unique(idx_traits[miss_idx[, 2L]])
      miss_ids    <- unique(ind_order[miss_idx[, 1L]])
      stop("No genetic values in ind_tgv for index '", idx_name,
           "' trait(s) ", paste(miss_traits, collapse = ", "),
           " for individual(s): ",
           paste(utils::head(miss_ids, 5), collapse = ", "),
           if (length(miss_ids) > 5) ", ..." else "",
           ". Include those traits in the trait_name argument of this ",
           "add_tgv() call.", call. = FALSE)
    }

    ti_df <- data.frame(
      id_true_index    = seq.int(next_int_id(conn, "ind_true_index",
                                             "id_true_index"),
                                 length.out = length(ind_order)),
      id_ind           = ind_order,
      index_name       = idx_name,
      weight_type      = wt,
      component_name   = component_name,
      true_index_value = as.numeric(val_mat %*% wt_vec),
      stringsAsFactors = FALSE)

    duckdb::duckdb_register(conn, tmp, ti_df)
    tryCatch({
      DBI::dbExecute(conn, "BEGIN TRANSACTION")
      if (overwrite_index) {
        DBI::dbExecute(conn, paste0(
          "DELETE FROM ind_true_index WHERE ", key_sql,
          " AND id_ind IN (SELECT id_ind FROM ", tmp, ")"))
      }
      DBI::dbExecute(conn, paste0(
        "INSERT INTO ind_true_index (id_true_index, id_ind, index_name, ",
        "weight_type, component_name, true_index_value) ",
        "SELECT id_true_index, id_ind, index_name, weight_type, ",
        "component_name, true_index_value FROM ", tmp))
      DBI::dbExecute(conn, "COMMIT")
    }, error = function(e) {
      try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE)
      stop(e)
    })
    duckdb::duckdb_unregister(conn, tmp)
    message("Computed true index '", idx_name, "' (", wt, ", ", component_name,
            ") for ", nrow(ti_df), " individuals.")
  }
  invisible(NULL)
}


#' Replace an individual's `ind_tgv` rows for a trait, in one transaction
#'
#' A plain upsert would leave a stale component behind when the model changes
#' shape — drop the dominance terms and the old `'dominance'` row would survive
#' and keep being summed by `ind_tgv_total`. So every row of a re-evaluated
#' (individual, trait) whose component is no longer produced is deleted, and
#' the rest are upserted on (individual, trait, component), which keeps any
#' custom columns on a surviving row.
#'
#' @param extras Named list of prepared custom-column vectors (one value per
#'   row of `res`), from [prepare_extra_cols()].
#' @keywords internal
#' @noRd
.gev_write_tgv <- function(conn, res, traits, ids, extras = list()) {
  stamp   <- as.character(round(as.numeric(Sys.time()) * 1000))
  ind_tmp <- paste0("_tgv_ind_", stamp)
  row_tmp <- paste0("_tgv_row_", stamp)

  out <- data.frame(
    id_tgv         = seq.int(next_int_id(conn, "ind_tgv", "id_tgv"),
                             length.out = nrow(res)),
    id_ind         = res$id_ind,
    trait_name     = res$trait_name,
    component_name = res$component_name,
    tgv_value      = as.numeric(res$tgv_value),
    stringsAsFactors = FALSE)
  for (nm in names(extras)) out[[nm]] <- extras[[nm]]

  duckdb::duckdb_register(conn, ind_tmp,
                          data.frame(id_ind = ids, stringsAsFactors = FALSE))
  duckdb::duckdb_register(conn, row_tmp, out)
  on.exit({
    try(duckdb::duckdb_unregister(conn, ind_tmp), silent = TRUE)
    try(duckdb::duckdb_unregister(conn, row_tmp), silent = TRUE)
  }, add = TRUE)

  cols   <- paste0("\"", names(out), "\"", collapse = ", ")
  update <- paste0("\"", c("tgv_value", names(extras)), "\" = EXCLUDED.\"",
                   c("tgv_value", names(extras)), "\"", collapse = ", ")

  tryCatch({
    DBI::dbExecute(conn, "BEGIN TRANSACTION")
    DBI::dbExecute(conn, paste0(
      "DELETE FROM ind_tgv AS t WHERE t.trait_name IN (",
      sql_in_list(traits, what = "trait name"), ") ",
      "AND t.id_ind IN (SELECT id_ind FROM ", ind_tmp, ") ",
      "AND NOT EXISTS (SELECT 1 FROM ", row_tmp, " AS r ",
      "WHERE r.id_ind = t.id_ind AND r.trait_name = t.trait_name ",
      "AND r.component_name = t.component_name)"))
    if (nrow(out) > 0L) {
      DBI::dbExecute(conn, paste0(
        "INSERT INTO ind_tgv (", cols, ") SELECT ", cols, " FROM ", row_tmp,
        " ON CONFLICT (id_ind, trait_name, component_name) DO UPDATE SET ",
        update))
    }
    DBI::dbExecute(conn, "COMMIT")
  }, error = function(e) {
    try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE)
    stop(e)
  })
  invisible(NULL)
}

#' Compute and store true genetic values
#'
#' @description
#' Evaluates the stored effect model for each individual in the current subset
#' and writes the result to `ind_tgv`, one row per
#' (individual x trait x **component**). The total is the derived view
#' `ind_tgv_total`, never a stored row — a stored total would make every
#' `SUM(tgv_value)` double-count.
#'
#' Unlike [add_tbv()], this evaluates **every** term of the trait: additive,
#' dominance, hand-entered genotype surfaces and multi-locus interactions, under
#' every effect owner. `component_name` records how each term was *declared*:
#'
#' | `component_name` | Terms it collects |
#' |---|---|
#' | `"order1_additive"` | single-member terms whose contrast is `additive` |
#' | `"order1_dominance"` | single-member `dominance` terms |
#' | `"order1_other"` | single-member `indicator` terms (a hand-entered surface) |
#' | `"interaction"` | any term with two or more members |
#'
#' These are **model-structure components, not variance components.** A
#' functional A x A term contributes to \eqn{V_A}, \eqn{V_D} *and* \eqn{V_I} in
#' the statistical sense; the names carry the declared order precisely so they
#' cannot be read as an orthogonal decomposition.
#'
#' `ind_tgv` stores the **raw sum of the stored terms — no mean is added.** A
#' pure Cockerham model yields centered deviations; a raw `indicator` surface
#' yields absolute genotypic values with a non-zero mean by construction, and is
#' not silently re-centered. See [ad_terms()], which reports the implied genetic
#' mean and writes it nowhere.
#'
#' Writes are idempotent: an individual's rows for a trait are replaced, so
#' re-evaluating after changing the effect model never leaves stale components
#' behind.
#'
#' @param tbl A `tidybreed_table` from [get_table()] (optionally piped through
#'   [dplyr::filter()]). The table must contain an `id_ind` column.
#' @param trait_name Character vector of trait name(s). When `NULL` (default),
#'   all traits in `trait_meta` are used (in `id_trait` order).
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
#' @seealso [add_tbv()] for the breeding value, [define_genome_effects()] and
#'   [define_additive_effects()] for writing the terms this evaluates.
#'
#' @examples
#' \dontrun{
#' # Every component of the model, for one generation
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 2L) |>
#'   add_tgv("ADG")
#'
#' # The components, and the derived total
#' pop |> get_table("ind_tgv") |> dplyr::collect()
#' pop |> get_table("ind_tgv_total") |> dplyr::collect()
#' }
#' @export
add_tgv <- function(tbl, trait_name = NULL) {

  stopifnot(inherits(tbl, "tidybreed_table"))
  pop <- tbl$pop
  validate_tidybreed_pop(pop)
  conn <- pop$db_conn

  traits <- .gev_resolve_traits(conn, trait_name)
  ids    <- .gev_subset_ids(tbl, "TGV")
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

  .gev_write_tgv(conn, res, traits, ids)
  for (t in traits) {
    message("Computed TGV for ", length(ids), " individuals on trait '", t,
            "' (", paste(sort(unique(res$component_name[res$trait_name == t])),
                         collapse = ", "), ").")
  }
  invisible(pop)
}


#' Replace an individual's `ind_tgv` rows for a trait, in one transaction
#'
#' A plain upsert would leave a stale component behind when the model changes
#' shape — drop the dominance terms and the old `'order1_dominance'` row would
#' survive and keep being summed by `ind_tgv_total`. So the delete is by
#' (individual, trait), not by (individual, trait, component).
#'
#' @keywords internal
#' @noRd
.gev_write_tgv <- function(conn, res, traits, ids) {
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

  duckdb::duckdb_register(conn, ind_tmp,
                          data.frame(id_ind = ids, stringsAsFactors = FALSE))
  duckdb::duckdb_register(conn, row_tmp, out)
  on.exit({
    try(duckdb::duckdb_unregister(conn, ind_tmp), silent = TRUE)
    try(duckdb::duckdb_unregister(conn, row_tmp), silent = TRUE)
  }, add = TRUE)

  tryCatch({
    DBI::dbExecute(conn, "BEGIN TRANSACTION")
    DBI::dbExecute(conn, paste0(
      "DELETE FROM ind_tgv WHERE trait_name IN (",
      sql_in_list(traits, what = "trait name"), ") ",
      "AND id_ind IN (SELECT id_ind FROM ", ind_tmp, ")"))
    if (nrow(out) > 0L) {
      DBI::dbExecute(conn, paste0(
        "INSERT INTO ind_tgv (id_tgv, id_ind, trait_name, component_name, ",
        "tgv_value) SELECT id_tgv, id_ind, trait_name, component_name, ",
        "tgv_value FROM ", row_tmp))
    }
    DBI::dbExecute(conn, "COMMIT")
  }, error = function(e) {
    try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE)
    stop(e)
  })
  invisible(NULL)
}

#' Contributor lookups for composite and `formula_tbv` phenotypes
#'
#' @description
#' The one place a contributor of a composite phenotype — the individual
#' itself, its dam or sire, or its group-mates — becomes a per-individual
#' value. Both `.assemble_composite_tbv()` (`phenotype_components`) and
#' `.build_tbv_env()` (`formula_tbv`) read through these helpers, and
#' `.ap_materialize_tbvs()` uses the same group lookup to decide whose TBVs
#' to compute first. Individual ids never appear in SQL text: every lookup
#' joins a registered view.
#'
#' Group semantics (SGE / Bijma): a focal's group-mates are the *other*
#' individuals with the same value of `group_column` in `group_table`; the
#' aggregate is over the mates that have a TBV; a focal with no mates gets
#' `0`; a focal whose group value is `NULL` gets `NA` (a missing component).
#' `group_table` must have exactly one row per focal individual.
#'
#' @name contributor_tbv
#' @keywords internal
NULL


#' One value of `column` from `table` per id, in `ids` order
#'
#' Requires exactly one row per id; zero or several rows is an error naming
#' the table, the column and up to five example ids. The column is returned
#' in its native type (`NA` for `NULL`).
#'
#' @param what Prefix for error messages, e.g. `"Residual condition lookup"`.
#' @keywords internal
.read_one_per_id <- function(conn, table, column, ids, what) {
  if (!table %in% DBI::dbListTables(conn)) {
    stop(what, ": table '", table, "' does not exist.", call. = FALSE)
  }
  if (!column %in% DBI::dbListFields(conn, table)) {
    stop(what, ": column '", column, "' not found in table '", table, "'.",
         call. = FALSE)
  }
  rows    <- .ap_read_by_id(conn, table, ids, column)
  counts  <- table(rows$id_ind)
  missing <- setdiff(ids, rows$id_ind)
  several <- names(counts)[counts > 1L]
  if (length(missing) > 0L || length(several) > 0L) {
    bad <- c(missing, several)
    stop(
      what, ": '", table, "' must have exactly one row per individual to ",
      "read '", column, "', but ",
      if (length(missing) > 0L) paste0(length(missing), " planned individual(s) have no row"),
      if (length(missing) > 0L && length(several) > 0L) " and ",
      if (length(several) > 0L) paste0(length(several), " have several rows"),
      " (e.g. ", paste(utils::head(bad, 5L), collapse = ", "),
      if (length(bad) > 5L) paste0(" ... +", length(bad) - 5L, " more") else "",
      ").", call. = FALSE)
  }
  rows[[column]][match(ids, rows$id_ind)]
}


#' `ind_tbv.tbv_value` of a trait per id (`NA` for a `NA` id or no TBV)
#' @keywords internal
.tbv_by_id <- function(conn, trait_name, ids) {
  ids  <- as.character(ids)
  have <- ids[!is.na(ids)]
  out  <- rep(NA_real_, length(ids))
  if (length(have) == 0L) return(out)
  rows <- .ap_read_by_id(conn, "ind_tbv", have, "tbv_value",
                         where = paste0("t.trait_name = ",
                                        DBI::dbQuoteLiteral(conn, trait_name)))
  out[!is.na(ids)] <- rows$tbv_value[match(have, rows$id_ind)]
  out
}


#' The group value of each focal individual, registered as `__ap_focal_groups`
#'
#' Focals with a `NULL` group value are left out of the view. The caller
#' unregisters the view.
#' @return Logical: which focals have a group value.
#' @keywords internal
.register_focal_groups <- function(conn, focal_ids, group_column, group_table,
                                   what) {
  gv  <- .read_one_per_id(conn, group_table, group_column, focal_ids, what)
  has <- !is.na(gv)
  f <- data.frame(id_ind = focal_ids[has], stringsAsFactors = FALSE)
  f$group_val <- gv[has]
  duckdb::duckdb_register(conn, "__ap_focal_groups", f)
  has
}

# The members of a group table, one row per (individual, group value).
.group_members_sql <- function(group_column, group_table) {
  paste0("(SELECT DISTINCT g.id_ind, g.\"", group_column, "\" AS group_val ",
         "FROM ", group_table, " AS g WHERE g.\"", group_column,
         "\" IS NOT NULL)")
}


#' Every individual sharing a group with a focal individual (focals included)
#' @keywords internal
.group_members <- function(conn, focal_ids, group_column, group_table, what) {
  if (length(focal_ids) == 0L) return(character(0))
  has <- .register_focal_groups(conn, focal_ids, group_column, group_table, what)
  on.exit(try(duckdb::duckdb_unregister(conn, "__ap_focal_groups"),
              silent = TRUE), add = TRUE)
  if (!any(has)) return(character(0))
  DBI::dbGetQuery(conn, paste0(
    "SELECT DISTINCT m.id_ind FROM ", .group_members_sql(group_column, group_table),
    " AS m JOIN (SELECT DISTINCT group_val FROM __ap_focal_groups) AS f ",
    "ON m.group_val = f.group_val"))$id_ind
}


#' Aggregated group-mate TBV per focal individual
#'
#' @param aggregation `"sum"` or `"mean"` over the mates that have a TBV.
#' @return Numeric per focal: the aggregate, `0` with no such mates, `NA`
#'   with no group value.
#' @keywords internal
.group_mate_tbv <- function(conn, trait_name, focal_ids, group_column,
                            group_table, aggregation, what) {
  if (!aggregation %in% c("sum", "mean")) {
    stop(what, ": aggregation must be 'sum' or 'mean', not '", aggregation,
         "'.", call. = FALSE)
  }
  out <- rep(NA_real_, length(focal_ids))
  if (length(focal_ids) == 0L) return(out)
  has <- .register_focal_groups(conn, focal_ids, group_column, group_table, what)
  on.exit(try(duckdb::duckdb_unregister(conn, "__ap_focal_groups"),
              silent = TRUE), add = TRUE)
  if (!any(has)) return(out)
  agg <- DBI::dbGetQuery(conn, paste0(
    "SELECT f.id_ind, SUM(t.tbv_value) AS total, COUNT(t.tbv_value) AS n_mates ",
    "FROM __ap_focal_groups AS f ",
    "LEFT JOIN ", .group_members_sql(group_column, group_table), " AS m ",
    "ON m.group_val = f.group_val AND m.id_ind <> f.id_ind ",
    "LEFT JOIN ind_tbv AS t ON t.id_ind = m.id_ind AND t.trait_name = ",
    DBI::dbQuoteLiteral(conn, trait_name), " ",
    "GROUP BY f.id_ind"))
  k <- match(focal_ids[has], agg$id_ind)
  total <- agg$total[k]
  n_mates <- agg$n_mates[k]
  value <- ifelse(n_mates == 0L, 0,
                  if (aggregation == "mean") total / n_mates else total)
  out[has] <- value
  out
}

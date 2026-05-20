#' Add residual covariance entries for observed phenotypes
#'
#' @description
#' Writes rows to `phenotype_residual_cov` representing the residual
#' (co)variance matrix for one or more phenotypes, optionally conditioned on a
#' group variable (e.g. farm, sex). Both `(i,j)` and `(j,i)` pairs are stored.
#'
#' Three typical call patterns:
#'
#' 1. **Called by [define_phenotype()] internally** when `residual_var` is
#'    supplied (scalar diagonal, unconditional).
#' 2. **Called by [define_effect_cov_matrix()]** when `effect_name = "residual"`
#'    to store a full multi-trait unconditional R matrix.
#' 3. **Called directly** to add group-specific (heterogeneous) residual rows
#'    after [define_phenotype()] has written the unconditional default.
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_names Character vector of phenotype names (must match
#'   `rownames(cov_matrix)` and `colnames(cov_matrix)`). For a single
#'   phenotype, a scalar is accepted.
#' @param cov_matrix Numeric matrix. Must be symmetric with
#'   `dimnames(cov_matrix)` matching `phenotype_names`. For a single phenotype
#'   the matrix is `1×1`.
#' @param condition_column Character or `NULL`. Column in `condition_table` used
#'   to look up group membership at phenotype time. `NULL` (default) stores an
#'   unconditional default row (the fallback when no group-specific row exists).
#' @param condition_table Character. Table containing `condition_column`.
#'   Default `"ind_meta"`.
#' @param condition_level Character or `NULL`. Specific level of `condition_column`
#'   this row applies to. `NULL` (default) = unconditional fallback row.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_phenotype()], [define_effect_cov_matrix()]
#'
#' @examples
#' \dontrun{
#' # Heterogeneous residual by sex for turkey BW
#' pop <- pop |>
#'   define_phenotype("BW_turkey",
#'     trait_type   = "continuous",
#'     mean         = 8000,
#'     residual_var = 600) |>   # unconditional default
#'   add_residual_cov("BW_turkey",
#'     cov_matrix       = matrix(900, 1, 1, dimnames = list("BW_turkey","BW_turkey")),
#'     condition_column = "sex",
#'     condition_level  = "M") |>
#'   add_residual_cov("BW_turkey",
#'     cov_matrix       = matrix(400, 1, 1, dimnames = list("BW_turkey","BW_turkey")),
#'     condition_column = "sex",
#'     condition_level  = "F")
#' }
#' @export
add_residual_cov <- function(pop,
                             phenotype_names,
                             cov_matrix,
                             condition_column = NULL,
                             condition_table  = "ind_meta",
                             condition_level  = NULL) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)

  phenotype_names <- as.character(phenotype_names)
  stopifnot(length(phenotype_names) >= 1)

  # Validate cov_matrix
  if (!is.matrix(cov_matrix) || !is.numeric(cov_matrix)) {
    stop("`cov_matrix` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(cov_matrix) != ncol(cov_matrix)) {
    stop("`cov_matrix` must be square.", call. = FALSE)
  }
  rn <- rownames(cov_matrix)
  cn <- colnames(cov_matrix)
  if (is.null(rn) || is.null(cn) || !identical(sort(rn), sort(phenotype_names)) ||
      !identical(sort(cn), sort(phenotype_names))) {
    stop(
      "`cov_matrix` dimnames must match `phenotype_names`: ",
      paste(phenotype_names, collapse = ", "),
      call. = FALSE
    )
  }

  # Symmetry check
  tol <- 1e-9
  if (length(phenotype_names) > 1 &&
      max(abs(cov_matrix - t(cov_matrix))) > tol) {
    stop("`cov_matrix` must be symmetric.", call. = FALSE)
  }

  pop <- ensure_trait_tables(pop)

  # Delete existing rows for this exact (phenotype pair, condition) combination
  pn1_in <- paste0("'", gsub("'", "''", phenotype_names), "'", collapse = ", ")
  cond_col_sql <- if (is.null(condition_column)) "IS NULL"
                  else paste0("= '", gsub("'", "''", condition_column), "'")
  cond_lvl_sql <- if (is.null(condition_level)) "IS NULL"
                  else paste0("= '", gsub("'", "''", condition_level), "'")

  DBI::dbExecute(pop$db_conn, paste0(
    "DELETE FROM phenotype_residual_cov ",
    "WHERE phenotype_name_1 IN (", pn1_in, ")",
    "  AND phenotype_name_2 IN (", pn1_in, ")",
    "  AND condition_column ", cond_col_sql,
    "  AND condition_level  ", cond_lvl_sql
  ))

  # Build all n² rows (both directions for off-diagonal)
  n <- length(phenotype_names)
  pn_order <- phenotype_names

  rows <- list()
  first_id <- next_int_id(pop$db_conn, "phenotype_residual_cov", "id_residual_cov")
  k <- 0L

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      k <- k + 1L
      rows[[k]] <- list(
        id_residual_cov  = first_id + k - 1L,
        phenotype_name_1 = pn_order[i],
        phenotype_name_2 = pn_order[j],
        cov_value        = cov_matrix[pn_order[i], pn_order[j]],
        condition_column = if (is.null(condition_column)) NA_character_
                           else condition_column,
        condition_table  = condition_table,
        condition_level  = if (is.null(condition_level)) NA_character_
                           else as.character(condition_level)
      )
    }
  }

  cov_df <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  cov_df$id_residual_cov  <- as.integer(cov_df$id_residual_cov)
  cov_df$cov_value        <- as.numeric(cov_df$cov_value)
  cov_df$condition_column <- as.character(cov_df$condition_column)
  cov_df$condition_table  <- as.character(cov_df$condition_table)
  cov_df$condition_level  <- as.character(cov_df$condition_level)

  DBI::dbWriteTable(pop$db_conn, "phenotype_residual_cov", cov_df, append = TRUE)

  invisible(pop)
}

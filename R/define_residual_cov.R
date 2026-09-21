#' Define residual covariance entries for observed phenotypes
#'
#' @description
#' Writes rows to `phenotype_var_comp` (with `effect_name = "residual"`)
#' representing the residual (co)variance matrix for one or more phenotypes,
#' optionally conditioned on a group variable (e.g. farm, sex). Both `(i,j)`
#' and `(j,i)` pairs are stored.
#'
#' Three typical call patterns:
#'
#' 1. **Called by [define_phenotype()] internally** when `residual_var` is
#'    supplied (scalar diagonal, unconditional).
#' 2. **Called by [define_effect_cov_matrix()]** when `effect_name = "residual"`
#'    to store a full multi-phenotype unconditional R matrix.
#' 3. **Called directly** to declare group-specific (heterogeneous) residual
#'    strata, one call per `condition_level`.
#'
#' @section A block is declared in one call, as a complete matrix:
#' The phenotypes that share any stored residual pair row form a *covariance
#' block*; a pair written as an explicit `0` still joins its two phenotypes.
#' A call must name a whole block: declaring `{A, B}` and then `{B, C}` is an
#' error (the block would be `{A, B, C}` with `Cov(A, C)` undeclared), and so is
#' redeclaring `{A, B}` or `A` alone once `{A, B, C}` exists. Redeclare the
#' complete block instead, writing `0` for uncorrelated pairs. The matrix must
#' be symmetric, finite and positive semi-definite; a rejected call changes
#' nothing.
#'
#' **Strata.** Conditional calls (`condition_column` + `condition_level`) add
#' strata to the block. Every stratum names the same phenotypes and a block has
#' one `condition_column`; to grow a block that already has several strata,
#' clear its rows from `phenotype_var_comp` with [remove_rows()] and redeclare
#' each stratum.
#'
#' **How the block is sampled.** [add_phenotype()] draws each record's
#' residual from this block conditional on the residuals the same individual
#' has already realized for the block's other phenotypes at the same
#' `pheno_number` — in the same call or any earlier one — so the declared
#' covariance holds whether the phenotypes are recorded together or a
#' hundred simulated days apart with culling in between. With strata, each
#' record draws from the stratum its `condition_column` value selects,
#' falling back to the unconditional stratum when the value is `NULL` or
#' matches none (an error if there is no unconditional stratum).
#'
#' **Realized draws lock the block.** Once any `ind_phenotype` row of a member
#' has a non-`NULL` `residual_value`, the block cannot be redefined; the error
#' gives the [remove_rows()] call that clears those rows. Every phenotype in
#' the block that is already defined must carry the same
#' `condition_change_action` (see [define_phenotype()]).
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_names Character vector of phenotype names (must match
#'   `rownames(cov_matrix)` and `colnames(cov_matrix)`). For a single
#'   phenotype, a scalar is accepted.
#' @param cov_matrix Numeric matrix with `dimnames(cov_matrix)` matching
#'   `phenotype_names`. For a single phenotype the matrix is `1×1`.
#' @param condition_column Character or `NULL`. Column in `condition_table` used
#'   to look up group membership at phenotype time. `NULL` (default) declares
#'   the unconditional stratum. Must be supplied together with
#'   `condition_level`.
#' @param condition_table Character. Table containing `condition_column`.
#'   Default `"ind_meta"`. Ignored (stored as `NULL`) for the unconditional
#'   stratum.
#' @param condition_level Character or `NULL`. Level of `condition_column` this
#'   stratum applies to. `NULL` (default) = unconditional stratum.
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
#'     type         = "continuous",
#'     mean         = 8000,
#'     residual_var = 600) |>   # unconditional default
#'   define_residual_cov("BW_turkey",
#'     cov_matrix       = matrix(900, 1, 1, dimnames = list("BW_turkey","BW_turkey")),
#'     condition_column = "sex",
#'     condition_level  = "M") |>
#'   define_residual_cov("BW_turkey",
#'     cov_matrix       = matrix(400, 1, 1, dimnames = list("BW_turkey","BW_turkey")),
#'     condition_column = "sex",
#'     condition_level  = "F")
#' }
#' @export
define_residual_cov <- function(pop,
                                phenotype_names,
                                cov_matrix,
                                condition_column = NULL,
                                condition_table  = "ind_meta",
                                condition_level  = NULL) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)

  phenotype_names <- as.character(phenotype_names)
  if (length(phenotype_names) < 1L) {
    stop("`phenotype_names` must name at least one phenotype.", call. = FALSE)
  }
  lapply(phenotype_names, validate_sql_identifier, what = "phenotype name")

  if (is.null(condition_column)) {
    condition_table <- NULL
  } else {
    validate_sql_identifier(condition_column, what = "condition_column")
    validate_sql_identifier(condition_table,  what = "condition_table")
  }
  if (!is.null(condition_level) &&
      (length(condition_level) != 1L || is.na(condition_level))) {
    stop("`condition_level` must be a single non-missing value.", call. = FALSE)
  }

  pop <- ensure_trait_tables(pop)

  write_phenotype_cov_block(
    pop$db_conn, "residual", phenotype_names, cov_matrix,
    condition_column = condition_column,
    condition_table  = condition_table,
    condition_level  = condition_level,
    caller           = "define_residual_cov()"
  )

  invisible(pop)
}

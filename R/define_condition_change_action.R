#' Set `condition_change_action` for a whole residual covariance block
#'
#' `condition_change_action` decides what [add_phenotype()] does when a
#' correlated phenotype's stored residual was drawn under a different residual
#' `condition_level` than the current record resolves to — an animal that moved
#' farms between two records, say. It is a property of the **residual
#' covariance block**, not of one phenotype: D6 requires every member of a
#' block to carry the same value, so this function is the only way to change it
#' once a block of two or more exists.
#'
#' [define_phenotype()] sets the value at registration time, when the phenotype
#' is still a block of one. Afterwards the two are locked together — flipping
#' either one alone through `define_phenotype(overwrite = TRUE)` would leave
#' the block disagreeing, which D6 refuses, and restating the phenotype would
#' reset every other column of its `phenotype_meta` row along the way. This
#' function changes exactly one column, on every member, in one transaction.
#'
#' Unlike the covariance matrix itself, the action is **not** locked by
#' realized draws (D3). It governs how *future* records condition on stored
#' residuals and says nothing about the ones already drawn, so changing it
#' mid-simulation is legitimate and leaves `ind_phenotype` and
#' `phenotype_random_effects` untouched.
#'
#' @param pop A `tidybreed_pop`.
#' @param phenotype_name Character scalar. Any member of the block; the value
#'   is written to every member.
#' @param condition_change_action `"error"` (stop when a stored residual comes
#'   from another stratum) or `"independent"` (drop it from the conditioning
#'   set and warn).
#'
#' @return The `tidybreed_pop`, invisibly.
#'
#' @seealso [define_phenotype()], [define_residual_cov()], [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' # A and B share a farm-conditioned residual block, both at the default
#' # "error". An animal that changes farms between records should now be
#' # drawn independently of its stale residual rather than stopping the run.
#' pop <- define_condition_change_action(pop, "A", "independent")
#' # -> also set on B (same residual covariance block)
#' }
#' @export
define_condition_change_action <- function(pop, phenotype_name,
                                           condition_change_action =
                                             c("error", "independent")) {
  stopifnot(inherits(pop, "tidybreed_pop"))
  if (!is.character(phenotype_name) || length(phenotype_name) != 1L ||
      is.na(phenotype_name)) {
    stop("'phenotype_name' must be a single non-missing phenotype name.",
         call. = FALSE)
  }
  condition_change_action <- match.arg(condition_change_action)
  conn <- pop$db_conn

  known <- DBI::dbGetQuery(conn, sprintf(
    "SELECT phenotype_name FROM phenotype_meta WHERE phenotype_name = %s",
    DBI::dbQuoteLiteral(conn, phenotype_name)))$phenotype_name
  if (length(known) == 0L) {
    stop("Phenotype '", phenotype_name, "' not found in phenotype_meta. ",
         "Register it with define_phenotype() first.", call. = FALSE)
  }

  # The whole block, not just the named member: the value is block-scoped, and
  # writing one member would leave the block in the state D6 refuses. A
  # phenotype in no block is a block of one and this reduces to a single row.
  members <- .pvc_block_members(conn, "residual", phenotype_name)
  in_list <- .pvc_in_list(conn, members)

  before <- DBI::dbGetQuery(conn, sprintf(
    "SELECT phenotype_name, condition_change_action FROM phenotype_meta
     WHERE phenotype_name IN (%s) ORDER BY phenotype_name", in_list))
  before$condition_change_action[is.na(before$condition_change_action)] <- "error"

  # A block member that has no phenotype_meta row yet cannot be written here;
  # it will take the value from its own define_phenotype() call, and the D6
  # check there compares it against the members written below.
  undefined <- setdiff(members, before$phenotype_name)

  changed <- before$phenotype_name[
    before$condition_change_action != condition_change_action]
  if (length(changed) == 0L && length(undefined) == 0L) {
    message("condition_change_action is already '", condition_change_action,
            "' for ", .pvc_set(members), "; nothing to do.")
    return(invisible(pop))
  }

  DBI::dbExecute(conn, "BEGIN TRANSACTION")
  ok <- FALSE
  on.exit(if (!ok) try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE),
          add = TRUE)
  DBI::dbExecute(conn, sprintf(
    "UPDATE phenotype_meta SET condition_change_action = %s
     WHERE phenotype_name IN (%s)",
    DBI::dbQuoteLiteral(conn, condition_change_action), in_list))
  DBI::dbExecute(conn, "COMMIT")
  ok <- TRUE

  message("condition_change_action = '", condition_change_action, "' set on ",
          .pvc_set(before$phenotype_name),
          if (nrow(before) > 1L) " \u2014 one residual covariance block" else "",
          ".")
  if (length(undefined) > 0L) {
    message("Block member(s) ", .pvc_set(sort(undefined)),
            " have no phenotype_meta row yet; set the same value in their ",
            "define_phenotype() call.")
  }
  invisible(pop)
}

#' Remove one scope of generated effects
#'
#' @description
#' Deletes every term a generator wrote for `trait_name` at exactly one
#' scope, `(line_name, parent_origin)` — the scope a re-run of
#' [define_additive_effects()] would replace — **whatever its kind**:
#' additive, dominance and interaction terms at that scope go together. A
#' generated model is calibrated as a whole, so it is removed as a whole and
#' never one component at a time (for "A without D", re-run the generator with
#' an additive-only target, which re-calibrates). Every other scope stands.
#' This is the only way to remove a generated variant:
#' [define_genome_effect_terms()] and
#' [remove_rows()] refuse the reserved owner `"generated"`, so that
#' `"generated"` keeps meaning "calibrated to the stored target".
#'
#' The usual reason is a variant added by mistake, such as a re-run with a new
#' `parent_origin`, which **adds** a variant beside the old one rather than
#' replacing it. Two parent scopes at one line are a legal model, but no
#' stored target describes their combined variance, so
#' `define_phenotype(prevalence = )` refuses such a trait.
#'
#' A [define_genome_effects()] model is common-scope (no `line_name`, no
#' `parent_origin`), so `remove_generated_effects(pop, trait_name)` removes it
#' whole: its additive, dominance and additive-by-additive terms together.
#'
#' Nothing else changes:
#' * **The stored targets in `trait_var_comp` stay.** A target with no terms
#'   of its kind is ignored by the `prevalence` threshold, but it still
#'   counts for generation: [define_additive_effects()] refuses a trait with
#'   a stored `dominance` or `additive_by_additive` target unless
#'   `trait_var_comp_tbl` selects the additive rows alone (or those targets
#'   are removed with [remove_rows()]); [define_genome_effects()] calibrates
#'   to all of them again.
#' * **The values already in `ind_tgv` stay.** They describe the old model
#'   until [add_tgv()] (or [add_phenotype()], which calls it) re-evaluates the
#'   individuals. If the removal leaves the trait with **no terms at all**,
#'   that re-evaluation is an error ("No genome effects found"), so the old
#'   values cannot be refreshed and are still what [get_table()] shows: write
#'   a new model before computing genetic values or recording phenotypes
#'   again. No new record is ever made from them.
#'
#' Removing a line's variant makes that line's copies fall back to the common
#' variant, where there is one.
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character vector. The trait(s) whose variant is removed.
#' @param line_name `NULL` (the population-wide scope) or one line name, as
#'   passed to [define_additive_effects()].
#' @param parent_origin `NULL` (both parents' copies), `1` (sire) or `2`
#'   (dam), as passed to [define_additive_effects()]. A generator with no
#'   scope arguments writes the common scope, `line_name = NULL,
#'   parent_origin = NULL`.
#'
#' @return The `tidybreed_pop`, invisibly. An error, with nothing deleted,
#'   when a trait has no generated terms at that scope.
#'
#' @seealso [define_additive_effects()], [define_genome_effects()]
#'
#' @examples
#' \dontrun{
#' # A paternal-only variant added by mistake next to the common one
#' pop <- remove_generated_effects(pop, "ADG", parent_origin = 1)
#' }
#' @export
remove_generated_effects <- function(pop, trait_name, line_name = NULL,
                                     parent_origin = NULL) {
  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  conn <- pop$db_conn

  if (!is.character(trait_name) || length(trait_name) == 0L ||
      anyNA(trait_name) || anyDuplicated(trait_name)) {
    stop("`trait_name` must be a character vector of distinct trait names.",
         call. = FALSE)
  }
  if (!is.null(line_name) &&
      (!is.character(line_name) || length(line_name) != 1L ||
       is.na(line_name) || !nzchar(line_name))) {
    stop("`line_name` must be NULL or one line name.", call. = FALSE)
  }
  if (!is.null(parent_origin) &&
      (length(parent_origin) != 1L || !parent_origin %in% c(1, 2))) {
    stop("`parent_origin` must be NULL, 1 (sire) or 2 (dam).", call. = FALSE)
  }
  if (!is.null(parent_origin)) parent_origin <- as.integer(parent_origin)

  scope <- .ge_scope_from_origin(.dae_scope(line_name, parent_origin),
                                 "replace_scope")
  model <- .ge_read_model(conn)
  drop  <- integer(0)
  for (t in trait_name) {
    ids <- .ge_resolve_deletes(model, t, GE_GENERATED_OWNER, "replace_scope",
                               scope, TRUE)
    if (length(ids) == 0L) {
      stop("Trait '", t, "' has no generated effects at scope ",
           .dae_scope_label(line_name, parent_origin), ". Nothing was ",
           "removed.", call. = FALSE)
    }
    drop <- c(drop, ids)
  }

  DBI::dbExecute(conn, "BEGIN TRANSACTION")
  ok <- FALSE
  on.exit(if (!ok) try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE),
          add = TRUE)
  .ge_delete_terms(conn, unique(drop))
  validate_genome_effects(conn)
  DBI::dbExecute(conn, "COMMIT")
  ok <- TRUE

  message("Removed ", length(unique(drop)), " generated term(s) for trait(s) ",
          paste(trait_name, collapse = ", "), " at scope ",
          .dae_scope_label(line_name, parent_origin), ".")
  invisible(pop)
}

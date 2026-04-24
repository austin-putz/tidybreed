#' Check that a trait exists in trait_meta
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character scalar.
#' @keywords internal
.check_trait_exists <- function(pop, trait_name) {
  n <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM trait_meta WHERE trait_name = '",
           trait_name, "'")
  )$n
  if (n == 0) {
    stop("Trait '", trait_name, "' not found. Call add_trait() first.",
         call. = FALSE)
  }
  invisible(NULL)
}


#' Handle duplicate effect_name: stop or delete existing row
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character scalar.
#' @param effect_name Character scalar.
#' @param overwrite Logical.
#' @keywords internal
.handle_effect_overwrite <- function(pop, trait_name, effect_name, overwrite) {
  existing <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM trait_effects WHERE trait_name = '",
           trait_name, "' AND effect_name = '", effect_name, "'")
  )$n
  if (existing > 0 && !overwrite) {
    stop("Effect '", effect_name, "' already exists for trait '", trait_name,
         "'. Use overwrite = TRUE.", call. = FALSE)
  }
  if (existing > 0 && overwrite) {
    DBI::dbExecute(
      pop$db_conn,
      paste0("DELETE FROM trait_effects WHERE trait_name = '", trait_name,
             "' AND effect_name = '", effect_name, "'")
    )
    # Also remove any stored random draws for this effect so they're resampled
    if ("trait_random_effects" %in% DBI::dbListTables(pop$db_conn)) {
      DBI::dbExecute(
        pop$db_conn,
        paste0("DELETE FROM trait_random_effects WHERE trait_name = '",
               trait_name, "' AND effect_name = '", effect_name, "'")
      )
    }
  }
  invisible(NULL)
}

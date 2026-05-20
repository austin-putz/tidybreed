#' Check that a phenotype exists in phenotype_meta
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_name Character scalar.
#' @keywords internal
.check_phenotype_exists <- function(pop, phenotype_name) {
  n <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM phenotype_meta WHERE phenotype_name = '",
           gsub("'", "''", phenotype_name), "'")
  )$n
  if (n == 0) {
    stop("Phenotype '", phenotype_name, "' not found in phenotype_meta. ",
         "Call define_phenotype() first.",
         call. = FALSE)
  }
  invisible(NULL)
}


#' Handle duplicate effect_name: stop or delete existing row
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_name Character scalar.
#' @param effect_name Character scalar.
#' @param overwrite Logical.
#' @keywords internal
.handle_effect_overwrite <- function(pop, phenotype_name, effect_name, overwrite) {
  pn_safe <- gsub("'", "''", phenotype_name)
  en_safe <- gsub("'", "''", effect_name)
  existing <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM trait_effects ",
           "WHERE phenotype_name = '", pn_safe, "'",
           " AND effect_name = '", en_safe, "'")
  )$n
  if (existing > 0 && !overwrite) {
    stop("Effect '", effect_name, "' already exists for phenotype '", phenotype_name,
         "'. Use overwrite = TRUE.", call. = FALSE)
  }
  if (existing > 0 && overwrite) {
    DBI::dbExecute(
      pop$db_conn,
      paste0("DELETE FROM trait_effects WHERE phenotype_name = '", pn_safe,
             "' AND effect_name = '", en_safe, "'")
    )
    if ("trait_random_effects" %in% DBI::dbListTables(pop$db_conn)) {
      DBI::dbExecute(
        pop$db_conn,
        paste0("DELETE FROM trait_random_effects WHERE phenotype_name = '",
               pn_safe, "' AND effect_name = '", en_safe, "'")
      )
    }
  }
  invisible(NULL)
}

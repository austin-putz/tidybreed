#' Store externally computed estimated breeding values (EBVs)
#'
#' @description
#' Inserts per-individual EBVs — typically the output of a BLUP / GBLUP /
#' ssGBLUP run — into the `ind_ebv` table. Each insertion is tagged with a
#' user-supplied `model` label so multiple evaluation runs can coexist for
#' the same individual and trait.
#'
#' Existing rows with the same `(id_ind, trait_name, model)` are replaced.
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait Character. Trait name (must exist in `trait_meta`).
#' @param ebv_df Data frame with columns `id_ind` and `ebv`. Optional
#'   columns: `acc` (accuracy, 0-1), `se` (standard error).
#' @param model Character. Label identifying the evaluation run, e.g.
#'   `"BLUP_2026Q1"`.
#' @param date_calc Date stored in `ind_ebv.date_calc`.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @examples
#' \dontrun{
#' ebv_df <- tibble::tibble(
#'   id_ind = c("A-1", "A-2", "A-3"),
#'   ebv    = c(0.3, -0.1, 0.5),
#'   acc    = c(0.6, 0.55, 0.7)
#' )
#' pop |> add_ebv(trait = "ADG", ebv_df = ebv_df, model = "ssGBLUP_v1")
#' }
#' @export
add_ebv <- function(pop, trait, ebv_df, model, date_calc = Sys.Date()) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(trait, what = "trait name")
  validate_sql_identifier(model, what = "model label")
  stopifnot(is.data.frame(ebv_df))

  if (!"ind_ebv" %in% pop$tables) {
    pop <- ensure_trait_tables(pop)
  }

  trait_exists <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM trait_meta WHERE trait_name = '",
           trait, "'")
  )$n
  if (trait_exists == 0) {
    stop("Trait '", trait, "' not found in trait_meta.", call. = FALSE)
  }

  required <- c("id_ind", "ebv")
  missing_cols <- setdiff(required, names(ebv_df))
  if (length(missing_cols) > 0) {
    stop("ebv_df is missing required columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  row <- tibble::tibble(
    id_ind     = as.character(ebv_df$id_ind),
    trait_name = trait,
    model      = model,
    ebv        = as.numeric(ebv_df$ebv),
    acc        = if ("acc" %in% names(ebv_df)) as.numeric(ebv_df$acc)
                 else NA_real_,
    se         = if ("se" %in% names(ebv_df)) as.numeric(ebv_df$se)
                 else NA_real_,
    date_calc  = as.Date(date_calc)
  )

  # Remove prior rows for the same (id_ind, trait, model)
  ids_sql <- paste0("'", row$id_ind, "'", collapse = ", ")
  DBI::dbExecute(
    pop$db_conn,
    paste0("DELETE FROM ind_ebv WHERE trait_name = '", trait,
           "' AND model = '", model,
           "' AND id_ind IN (", ids_sql, ")")
  )
  DBI::dbWriteTable(pop$db_conn, "ind_ebv", row, append = TRUE)

  message("Wrote ", nrow(row), " EBV rows for trait '", trait,
          "' under model '", model, "'.")
  invisible(pop)
}

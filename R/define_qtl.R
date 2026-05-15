#' Define QTL loci for one or more traits
#'
#' @description
#' Marks the loci in a filtered `genome_meta` table as QTL for the specified
#' trait(s), writing a `BOOLEAN` column `is_QTL_{trait_name}` per trait to
#' `genome_meta`. All other loci receive `FALSE`.
#'
#' Pipe `get_table("genome_meta")` — optionally through `dplyr::filter()` —
#' into this function. The unique `locus_id` values in the collected result
#' determine QTL membership for every trait named in `trait_name`.
#'
#' When `trait_name` is omitted, all traits currently in `trait_meta` are used.
#'
#' Returns the modified `tidybreed_pop` so the result can be piped directly
#' into [set_qtl_effects()].
#'
#' @param tbl A `tidybreed_table` from `get_table("genome_meta")` (optionally
#'   piped through `dplyr::filter()`). Must contain a `locus_id` column.
#' @param trait_name Character vector of trait name(s) that exist in
#'   `trait_meta`. Default: all traits in `trait_meta` (in id_trait order).
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_chip()], [set_qtl_effects()], [set_qtl_effects_multi()],
#'   [add_trait()]
#'
#' @examples
#' \dontrun{
#' # Define QTL for ADG on chromosome 1 loci, then set effects
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(chr == 1) |>
#'   define_qtl("ADG") |>
#'   set_qtl_effects("ADG", distribution = "normal")
#'
#' # Define the same loci as QTL for multiple traits in one call
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(locus_name %in% qtl_loci) |>
#'   define_qtl(c("ADG", "BW"))
#'
#' # Shared QTL pleiotropy: give BW the same QTL set as ADG
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(is_QTL_ADG == TRUE) |>
#'   define_qtl("BW")
#'
#' # Omit trait_name to define QTL for every trait in trait_meta
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(locus_name %in% qtl_loci) |>
#'   define_qtl()
#' }
#' @export
define_qtl <- function(tbl, trait_name = NULL) {

  stopifnot(inherits(tbl, "tidybreed_table"))
  pop <- tbl$pop
  validate_tidybreed_pop(pop)

  if (!"trait_meta" %in% pop$tables) {
    stop("No traits defined yet. Call add_trait() first.", call. = FALSE)
  }

  # Resolve trait_name: default to all traits in trait_meta
  if (is.null(trait_name)) {
    trait_name <- DBI::dbGetQuery(
      pop$db_conn,
      "SELECT trait_name FROM trait_meta ORDER BY id_trait"
    )$trait_name
    if (length(trait_name) == 0L) {
      stop("No traits found in trait_meta. Call add_trait() first.", call. = FALSE)
    }
  } else {
    stopifnot(is.character(trait_name), length(trait_name) >= 1L)
    lapply(trait_name, validate_sql_identifier, what = "trait name")

    existing_traits <- DBI::dbGetQuery(
      pop$db_conn,
      "SELECT trait_name FROM trait_meta"
    )$trait_name
    bad_traits <- setdiff(trait_name, existing_traits)
    if (length(bad_traits) > 0L) {
      stop(
        "Trait(s) not found in trait_meta: ",
        paste(bad_traits, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (!"genome_meta" %in% pop$tables) {
    stop("genome_meta table does not exist. Call initialize_genome() first.",
         call. = FALSE)
  }

  # Collect filtered locus IDs
  collected <- dplyr::collect(tbl)
  if (!"locus_id" %in% names(collected)) {
    stop(
      "The filtered table must contain 'locus_id'. ",
      "Pipe get_table(\"genome_meta\") into define_qtl().",
      call. = FALSE
    )
  }
  if (nrow(collected) == 0L) {
    stop(
      "No loci selected — your filter returned zero rows. ",
      "Check the filter before piping into define_qtl().",
      call. = FALSE
    )
  }

  selected_ids <- collected$locus_id
  ids_sql      <- paste(selected_ids, collapse = ", ")

  # For each trait: add BOOLEAN column (DEFAULT FALSE) and flip QTL loci to TRUE
  existing_cols <- DBI::dbListFields(pop$db_conn, "genome_meta")
  for (tn in trait_name) {
    qtl_col <- paste0("is_QTL_", tn)
    if (qtl_col %in% existing_cols) {
      stop(
        "Column '", qtl_col, "' already exists in genome_meta. ",
        "QTL for trait '", tn, "' have already been defined.",
        call. = FALSE
      )
    }
    DBI::dbExecute(
      pop$db_conn,
      paste0("ALTER TABLE genome_meta ADD COLUMN ", qtl_col, " BOOLEAN DEFAULT FALSE")
    )
    DBI::dbExecute(
      pop$db_conn,
      paste0("UPDATE genome_meta SET ", qtl_col, " = TRUE WHERE locus_id IN (", ids_sql, ")")
    )
    message("Defined QTL for '", tn, "' with ", length(selected_ids), " loci")
  }

  invisible(pop)
}

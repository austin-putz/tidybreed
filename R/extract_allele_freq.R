#' Extract per-locus allele frequencies from a filtered table
#'
#' @description
#' Computes the frequency of allele 1 at every locus in `genome_meta` from the
#' allele copies a filtered `tidybreed_table` selects. The table's identity
#' says *what kind of thing* is being selected; the user's [dplyr::filter()]
#' says *which ones*:
#'
#' | Piped table | Meaning | Copies counted |
#' |---|---|---|
#' | `founder_haplotypes` | the founder pool | the filtered pool rows |
#' | `ind_haplotype` | these allele copies | the filtered haplotype rows |
#' | any other table with `id_ind` | these individuals | every `ind_haplotype` row of the **distinct** selected `id_ind` |
#'
#' The frequency is computed in one SQL statement with the filter rendered as
#' a subquery; nothing but the per-locus result is collected into R. For an
#' `id_ind` table the frequency depends only on which individuals are selected,
#' never on how many rows each has — `ind_phenotype` with five records per
#' animal gives the same answer as `ind_meta` for the same animals.
#'
#' This is the single place tidybreed turns a population selection into `p`.
#' [define_additive_effects()] and [define_genome_effects()] both call it for
#' their `base_tbl`, so a base selection means the same population in either.
#' Users call it to obtain `p` for [ad_terms()].
#'
#' @section Two notions of line:
#' `ind_meta.line_name` is a pedigree label — an F1 is whatever it was called;
#' `ind_haplotype.line_origin` is the founding line each allele copy traces to,
#' carried through every [add_offspring()] call. They coincide for purebreds and
#' diverge for crosses. `get_table(pop, "ind_meta") |> filter(line_name ==
#' "F1")` counts every copy those animals carry, whichever line it came from;
#' `get_table(pop, "ind_haplotype") |> filter(line_origin == "Duroc")` counts
#' Duroc copies wherever they sit, at any cross depth.
#'
#' @param tbl A `tidybreed_table` from [get_table()], optionally filtered. One
#'   of the three shapes above. Columns removed with [dplyr::select()] are
#'   checked before any SQL runs.
#'
#' @return A tibble with exactly one row per `genome_meta` locus in ascending
#'   `locus_id` order: `locus_id` (integer), `locus_name` (character),
#'   `allele_freq` (double). `allele_freq` is `NA` at a locus the selection has
#'   no copies for; it is never `0` in that case, and an observed frequency of
#'   exactly `0` or `1` is a real value and is kept. Errors if no locus has a
#'   copy at all. Never warns and never writes.
#'
#' @seealso [define_additive_effects()], [define_genome_effects()], [ad_terms()].
#'
#' @examples
#' \dontrun{
#' # Founder pool of one line
#' p <- pop |> get_table("founder_haplotypes") |>
#'   dplyr::filter(line_name == "Duroc") |> extract_allele_freq()
#'
#' # Generation-0 animals
#' p <- pop |> get_table("ind_meta") |> dplyr::filter(gen == 0L) |>
#'   extract_allele_freq()
#'
#' # Duroc copies wherever they sit, including inside crossbreds
#' p <- pop |> get_table("ind_haplotype") |>
#'   dplyr::filter(line_origin == "Duroc") |> extract_allele_freq()
#'
#' # Feed p to ad_terms()
#' loci <- c("Locus_10", "Locus_44")
#' tt <- ad_terms(loci, a = c(0.4, 0.1), d = c(0.2, 0.05),
#'                p = p$allele_freq[match(loci, p$locus_name)],
#'                coding = "cockerham")
#' }
#' @export
extract_allele_freq <- function(tbl) {
  .validate_base_tbl(tbl, pop = NULL, arg = "tbl")   # no second pop to agree with
  conn <- tbl$pop$db_conn
  sub  <- as.character(dbplyr::sql_render(tbl$tbl))

  freq_sql <- switch(tbl$table_name,
    founder_haplotypes = paste0(
      "SELECT gm.locus_id, AVG(CAST(b.allele AS DOUBLE)) AS allele_freq ",
      "FROM (", sub, ") b JOIN genome_meta gm ON b.locus_name = gm.locus_name ",
      "GROUP BY gm.locus_id"),
    ind_haplotype = paste0(
      "SELECT b.locus_id, AVG(CAST(b.allele AS DOUBLE)) AS allele_freq ",
      "FROM (", sub, ") b GROUP BY b.locus_id"),
    # Individuals: a semi-join on the distinct selected ids, so row
    # multiplicity in the source table cannot weight anyone.
    paste0(
      "SELECT h.locus_id, AVG(CAST(h.allele AS DOUBLE)) AS allele_freq ",
      "FROM ind_haplotype h ",
      "JOIN (SELECT DISTINCT id_ind FROM (", sub, ") b) ids USING (id_ind) ",
      "GROUP BY h.locus_id"))

  # LEFT JOIN to genome_meta is what guarantees one row per locus and lets
  # "no copies" (NA) be told apart from an observed frequency of 0.
  out <- DBI::dbGetQuery(conn, paste0(
    "SELECT gm.locus_id, gm.locus_name, f.allele_freq ",
    "FROM genome_meta gm LEFT JOIN (", freq_sql, ") f USING (locus_id) ",
    "ORDER BY gm.locus_id"))

  if (nrow(out) == 0L) {
    stop("genome_meta has no loci; call define_genome() first.", call. = FALSE)
  }
  if (all(is.na(out$allele_freq))) {
    # Proves no usable locus matched -- not that the source had zero rows (a
    # founder pool whose locus_name keys miss genome_meta looks the same).
    if (tbl$table_name == "founder_haplotypes") {
      .founder_base_empty_error(conn, "this selection")
    }
    stop("The filtered base (", tbl$table_name, ") contains no allele copies ",
         "matching genome_meta loci.", call. = FALSE)
  }
  f <- out$allele_freq[!is.na(out$allele_freq)]
  if (any(!is.finite(f) | f < 0 | f > 1)) {
    stop("Allele frequencies outside [0, 1] were computed -- ",
         "are there allele values other than 0/1 in the haplotype rows?",
         call. = FALSE)
  }

  out$locus_id    <- as.integer(out$locus_id)
  out$allele_freq <- as.numeric(out$allele_freq)
  tibble::as_tibble(out)
}


#' Validate a `base_tbl` argument
#'
#' One validator for every place a population selection enters a genome-effect
#' writer, and for [extract_allele_freq()] itself. Checks the class, the pop,
#' that it is on the **same connection** as the writer's pop, and that the
#' columns the generated SQL touches are still projected. The last check exists
#' because `select.tidybreed_table()` replaces the lazy query and returns the
#' same wrapper, so a table whose `table_name` is `"ind_haplotype"` may no
#' longer expose `allele`; without it the failure would surface as a DuckDB
#' binder error inside generated SQL.
#'
#' `line_name` is deliberately not required for `founder_haplotypes`: nothing
#' in the frequency query reads it, and the two line-aware diagnostics
#' (`.dae_default_base()`'s Wahlund check, `.founder_base_empty_error()`) query
#' the physical table.
#'
#' @param base_tbl The object to validate.
#' @param pop The writer's `tidybreed_pop`, or `NULL` to use `base_tbl$pop`
#'   (the public-helper case, where there is no second pop to agree with).
#' @param arg Name used in messages.
#' @return `base_tbl`, invisibly.
#' @keywords internal
#' @noRd
.validate_base_tbl <- function(base_tbl, pop = NULL, arg = "base_tbl") {
  if (!inherits(base_tbl, "tidybreed_table")) {
    stop("'", arg, "' must be a tidybreed_table from get_table() |> filter(...). ",
         "Accepted shapes: founder_haplotypes, ind_haplotype, or any table ",
         "with an id_ind column.", call. = FALSE)
  }
  validate_tidybreed_pop(base_tbl$pop)
  if (!is.null(pop) && !identical(base_tbl$pop$db_conn, pop$db_conn)) {
    stop("'", arg, "' must be piped from get_table() on the same pop as 'tbl'.",
         call. = FALSE)
  }

  cols <- colnames(base_tbl$tbl)            # the projected columns, after select()
  need <- switch(base_tbl$table_name,
    founder_haplotypes = c("locus_name", "allele"),
    ind_haplotype      = c("locus_id", "allele"),
    "id_ind")
  miss <- setdiff(need, cols)
  if (length(miss)) {
    stop("'", arg, "' (", base_tbl$table_name, ") is missing column(s) ",
         toString(miss), " needed to compute allele frequencies. ",
         "Accepted shapes: founder_haplotypes, ind_haplotype, or any table ",
         "with an id_ind column.", call. = FALSE)
  }
  invisible(base_tbl)
}


#' Error for an empty founder-pool selection, listing the pools that exist
#'
#' Must be loud: a typo'd line name that silently yielded no rows would centre
#' every allele at 0 and contribute nothing to the Falconer `V_A`, producing
#' plausible-looking output. Queries the physical `founder_haplotypes` table.
#'
#' @param conn DBI connection.
#' @param what Description of the failed selection for the message:
#'   `"this selection"` (from the helper) or `"line 'A'"` (from the default
#'   resolver in `define_additive_effects()`).
#' @keywords internal
#' @noRd
.founder_base_empty_error <- function(conn, what) {
  avail <- DBI::dbGetQuery(conn,
    "SELECT DISTINCT line_name FROM founder_haplotypes ORDER BY line_name")$line_name
  if (length(avail) == 0L) {
    stop("founder_haplotypes table is empty.", call. = FALSE)
  }
  named   <- stats::na.omit(avail)
  parts   <- if (length(named)) paste0("'", named, "'", collapse = ", ") else NULL
  if (any(is.na(avail))) parts <- c(parts, "an unnamed (line_name = NULL) pool")
  stop("No founder_haplotypes rows for ", what, ". Available: ",
       paste(parts, collapse = " and "), ".", call. = FALSE)
}

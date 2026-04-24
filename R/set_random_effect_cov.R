#' Store a covariance matrix for correlated random effects across traits
#'
#' @description
#' Saves a covariance matrix `R` into `trait_random_effect_cov`, enabling
#' [add_phenotype()] to draw the named random effect jointly across traits via
#' `MVN(0, R)`. The diagonal of `R` must match the `variance` stored in
#' `trait_effects` for each `(trait, effect_name)` pair (within `tol`).
#'
#' Analogous to [set_residual_cov()] but for named group-level random effects.
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character. The effect name shared across traits (must
#'   exist in `trait_effects` for each trait listed in `traits`).
#' @param traits Character vector of at least two trait names.
#' @param R Numeric matrix of size `length(traits) x length(traits)`. Must be
#'   symmetric. Rows and columns are matched to `traits` in order.
#' @param tol Numeric. Tolerance for symmetry and diagonal-variance checks.
#' @param overwrite Logical. Replace existing entries for this
#'   `(effect_name, traits)` combination.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [add_effect_random()], [set_residual_cov()], [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' pop <- pop |>
#'   add_effect_random("ADG", "herd", source_column = "herd_id", variance = 150) |>
#'   add_effect_random("BW",  "herd", source_column = "herd_id", variance = 100) |>
#'   set_random_effect_cov("herd",
#'     traits = c("ADG", "BW"),
#'     R = matrix(c(150, 60, 60, 100), 2, 2))
#' }
#' @export
set_random_effect_cov <- function(pop,
                                  effect_name,
                                  traits,
                                  R,
                                  tol       = 1e-6,
                                  overwrite = FALSE) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(effect_name, what = "effect name")
  stopifnot(is.character(traits), length(traits) >= 2)
  lapply(traits, validate_sql_identifier, what = "trait name")

  n <- length(traits)
  if (!is.matrix(R) || !all(dim(R) == n)) {
    stop("`R` must be a ", n, " x ", n, " matrix.", call. = FALSE)
  }
  if (!isSymmetric(R, tol = tol)) {
    stop("`R` must be symmetric.", call. = FALSE)
  }
  if (any(diag(R) < 0)) {
    stop("All diagonal elements of `R` must be non-negative.", call. = FALSE)
  }

  if (!"trait_effects" %in% pop$tables) {
    stop("No trait tables exist. Call add_trait() first.", call. = FALSE)
  }

  # Validate diagonal matches variance in trait_effects
  traits_sql <- paste0("'", traits, "'", collapse = ", ")
  for (i in seq_along(traits)) {
    t <- traits[i]
    eff_row <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT variance FROM trait_effects WHERE trait_name = '", t,
             "' AND effect_name = '", effect_name, "'")
    )
    if (nrow(eff_row) == 0) {
      stop("Effect '", effect_name, "' not found for trait '", t,
           "'. Add it with add_effect_random() first.", call. = FALSE)
    }
    stored_var <- eff_row$variance[1]
    if (!is.na(stored_var) && abs(R[i, i] - stored_var) > tol) {
      stop("R[", i, ",", i, "] = ", R[i, i],
           " but trait_effects variance for '", t, "' / '", effect_name,
           "' = ", stored_var, ". Diagonal must match.", call. = FALSE)
    }
  }

  pop <- ensure_trait_tables(pop)

  if (overwrite) {
    DBI::dbExecute(
      pop$db_conn,
      paste0("DELETE FROM trait_random_effect_cov WHERE effect_name = '",
             effect_name, "' AND trait_1 IN (", traits_sql,
             ") AND trait_2 IN (", traits_sql, ")")
    )
  } else {
    existing <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT COUNT(*) AS n FROM trait_random_effect_cov ",
             "WHERE effect_name = '", effect_name, "' AND trait_1 IN (",
             traits_sql, ") AND trait_2 IN (", traits_sql, ")")
    )$n
    if (existing > 0) {
      stop("Covariance entries already exist for effect '", effect_name,
           "'. Use overwrite = TRUE.", call. = FALSE)
    }
  }

  rows <- vector("list", n * n)
  k <- 1L
  for (i in seq_along(traits)) {
    for (j in seq_along(traits)) {
      rows[[k]] <- tibble::tibble(
        effect_name = effect_name,
        trait_1     = traits[i],
        trait_2     = traits[j],
        cov         = R[i, j]
      )
      k <- k + 1L
    }
  }
  DBI::dbWriteTable(pop$db_conn, "trait_random_effect_cov",
                    do.call(rbind, rows), append = TRUE)

  message("Stored random-effect covariance for '", effect_name,
          "' across traits: ", paste(traits, collapse = ", "), ".")
  invisible(pop)
}

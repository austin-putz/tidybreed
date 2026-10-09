#' Define a variance-covariance matrix for any named effect
#'
#' @description
#' Single entry point for storing all variance and covariance data in tidybreed.
#' Routes to `trait_var_comp` for genetic effects and to `phenotype_var_comp` for
#' phenotype-level effects.
#'
#' Common `effect_name` values:
#'
#' * `"additive"` — additive genetic (co)variances (G matrix). Written to
#'   `trait_var_comp`. Used by [define_additive_effects()] when rescaling to
#'   target variance and as the sampling distribution for multi-trait draws.
#' * `"dominance"`, `"additive_by_additive"` — reserved genetic effects; no
#'   generator calibrates them yet. Written to
#'   `trait_var_comp`. Row/column names are trait names.
#' * `"residual"` — residual (co)variances (R matrix). Routed to
#'   `phenotype_var_comp` with `effect_name = "residual"`. Row/column names are
#'   phenotype names. Equivalent to calling [define_residual_cov()] with
#'   `condition_column = NULL`. Use this for a multi-phenotype correlated
#'   residual matrix; for a single scalar residual use `residual_var` in
#'   [define_phenotype()] instead.
#' * Any named random effect (`"hys"`, `"litter"`, `"pen"`, …) — written to
#'   `phenotype_var_comp`. Must match the `effect_name` used in
#'   [define_effect_random()]. Row/column names are phenotype names. Each
#'   level of the effect (each pen) then carries one draw per phenotype with
#'   this covariance, realized sequentially: whichever phenotype
#'   [add_phenotype()] generates first for a level draws marginally, and the
#'   others are later drawn conditional on what the level has stored —
#'   however many calls apart (see [define_effect_random()]).
#'
#' `define_effect_cov_matrix()` can be called **before** [define_trait()] or
#' [define_effect_random()] — no prior setup is required.
#'
#' All n² pairs are stored. For `phenotype_var_comp` effects the names form a
#' *covariance block* that is declared in one call, as a complete matrix: a
#' call that names a fragment or a strict subset of an existing block is an
#' error, the matrix must be positive semi-definite, and a block cannot be
#' redefined once draws exist under it (the error gives the [remove_rows()]
#' call that clears them). In a block of two or more phenotypes every
#' [define_effect_random()] row for the effect must use
#' `distribution = "normal"` and read the same `(source_column, source_table)`.
#' See [define_residual_cov()] for the full rules. A rejected call changes
#' nothing.
#'
#' **Genetic blocks are written once.** A genetic block (`"additive"`,
#' `"dominance"`, `"additive_by_additive"`) is validated as positive
#' semi-definite, stored at full double precision, and never overwritten: if
#' any row already exists for that `effect_name`, any of the named traits and
#' the same `line_name`, the call is an error, even when the matrix is
#' identical. The error gives the [remove_rows()] call that clears the stored
#' block. `trait_var_comp` is the single source of generation targets; the
#' effect generators read it and never overwrite it either.
#'
#' **A target is written before its effects are generated, or with them.** A
#' genetic block is refused when any of its traits already has terms of that
#' kind written by a generator (owner `"generated"`) at the block's scope:
#' line-`"C"` terms for `line_name = "C"`; for `line_name = NULL`, the
#' population-wide terms and the terms of every line with no block of its own
#' (the generator calibrated those to the population-wide target). Those terms were calibrated to the target they were
#' generated with, and a new target would not describe them; the prevalence
#' threshold of [define_phenotype()] trusts the stored target for exactly that
#' reason. The refusal holds after the old block is removed too. To change
#' the target, remove the old block and call [define_additive_effects()] with
#' `G =`, which writes the new target and re-draws the terms in one
#' transaction. A line's target written before that line's effects are
#' generated is accepted.
#'
#' `"additive_by_dominance"` and `"dominance_by_dominance"` are reserved for
#' future generators and refused. `"total"`, `"unpartitioned"` and
#' `"between_components"` are output names of the variance extractor and
#' refused as input.
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character. Label for the variance component, e.g.
#'   `"additive"`, `"residual"`, `"hys"`.
#' @param cov_matrix A numeric square matrix, or a single number when one
#'   trait/phenotype is named. Must be symmetric within `tol`. A named matrix
#'   must carry the same names as `trait_name`, in the same order (it is never
#'   relabelled); an unnamed one is taken in `trait_name` order.
#' @param trait_name Character vector of trait/phenotype names (length
#'   == `nrow(cov_matrix)`). Optional when the matrix has names.
#' @param line_name Character or `NULL` (default). Genetic effects only: the
#'   line whose generation target this is. `NULL` is the population-wide
#'   target, which a line without its own block falls back to.
#' @param tol Numeric. Tolerance for symmetry check (default `1e-9`).
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_trait()], [define_effect_random()], [define_additive_effects()],
#'   [add_phenotype()]
#'
#' @examples
#' \dontrun{
#' # Additive genetic covariance matrix → trait_var_comp
#' G <- matrix(c(100, -20, -20, 50), 2, 2,
#'             dimnames = list(c("ADG", "BF"), c("ADG", "BF")))
#' pop <- pop |>
#'   define_effect_cov_matrix("additive", G)
#'
#' # One trait: a number is a 1 x 1 matrix
#' pop <- pop |>
#'   define_effect_cov_matrix("additive", 0.25, trait_name = "WW")
#'
#' # Residual → phenotype_var_comp (effect_name = "residual")
#' R <- matrix(c(30, 5, 5, 10), 2, 2,
#'             dimnames = list(c("ADG", "BF"), c("ADG", "BF")))
#' pop <- pop |>
#'   define_effect_cov_matrix("residual", R)
#'
#' # Multi-phenotype HYS covariance → phenotype_var_comp (effect_name = "hys")
#' R_hys <- matrix(c(0.2, 0.05, 0.05, 0.3), 2, 2,
#'                 dimnames = list(c("ADG", "BF"), c("ADG", "BF")))
#' pop <- pop |>
#'   define_effect_cov_matrix("hys", R_hys)
#' }
#' @export
define_effect_cov_matrix <- function(pop,
                                     effect_name,
                                     cov_matrix,
                                     trait_name  = NULL,
                                     line_name   = NULL,
                                     tol         = 1e-9) {
  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)
  validate_sql_identifier(effect_name, what = "effect name")
  .check_effect_name_input(effect_name, genetic_ok = TRUE)

  cov_matrix <- .check_cov_dimnames(cov_matrix, trait_name, "cov_matrix")
  trait_name <- rownames(cov_matrix)

  if (!isSymmetric(unname(cov_matrix), tol = tol)) {
    stop("`cov_matrix` must be symmetric (max discrepancy: ",
         max(abs(cov_matrix - t(cov_matrix))), ").", call. = FALSE)
  }
  if (any(diag(cov_matrix) < 0)) {
    stop("Diagonal entries (variances) must be non-negative.", call. = FALSE)
  }

  is_genetic <- effect_name %in% GENETIC_EFFECT_NAMES
  if (!is.null(line_name)) {
    if (!is_genetic) {
      stop("`line_name` applies to genetic effects (",
           paste0("'", GENETIC_EFFECT_NAMES, "'", collapse = ", "),
           ") only, not to '", effect_name, "'.", call. = FALSE)
    }
    if (!is.character(line_name) || length(line_name) != 1L || is.na(line_name)) {
      stop("`line_name` must be a single character string or NULL.", call. = FALSE)
    }
    validate_sql_identifier(line_name, what = "line name")
  }

  if (is_genetic) {
    conn <- pop$db_conn
    DBI::dbExecute(conn, "BEGIN TRANSACTION")
    committed <- FALSE
    on.exit(if (!committed) try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE),
            add = TRUE)
    .tvc_refuse_under_generated(conn, effect_name, trait_name, line_name)
    .tvc_write_block(conn, effect_name, cov_matrix, line_name)
    DBI::dbExecute(conn, "COMMIT")
    committed <- TRUE
  } else if (identical(effect_name, "residual")) {
    define_residual_cov(pop,
      phenotype_names  = trait_name,
      cov_matrix       = cov_matrix,
      condition_column = NULL)
  } else {
    # Named random effect: one complete, unconditional block (D1/D3/section 5.6)
    write_phenotype_cov_block(
      pop$db_conn, effect_name, trait_name, cov_matrix,
      caller = "define_effect_cov_matrix()", tol = tol)
  }

  message("Stored '", effect_name, "' covariance matrix for: ",
          paste(trait_name, collapse = ", "),
          if (!is.null(line_name)) paste0(" (line '", line_name, "')"), ".")
  if (is_genetic) .qtl_rank_note(cov_matrix, effect_name)
  invisible(pop)
}


# ---------------------------------------------------------------------------
# Effect-name vocabulary (plans/import_qtl_effect_methods.md §6B rule 5)
# ---------------------------------------------------------------------------

#' Genetic variance components a target can be stored for
#'
#' Routed to `trait_var_comp` by [define_effect_cov_matrix()]. A finite list,
#' not a pattern: a name moves here when a generator can calibrate it.
#' @keywords internal
GENETIC_EFFECT_NAMES <- c("additive", "dominance", "additive_by_additive")

#' Genetic variance components reserved for future generators
#' @keywords internal
GENETIC_EFFECT_NAMES_FUTURE <- c("additive_by_dominance", "dominance_by_dominance")

#' Output-only names of the genetic variance extractor
#' @keywords internal
DERIVED_EFFECT_NAMES <- c("total", "unpartitioned", "between_components")

#' Refuse a reserved name where a user-named effect is expected
#'
#' `genetic_ok = TRUE` (only [define_effect_cov_matrix()]) lets the genetic
#' names through, because that function routes them to `trait_var_comp`. Every
#' phenotype-layer definer passes `FALSE`: a random or fixed effect named
#' `"additive"` would collide with the genetic vocabulary.
#' @keywords internal
.check_effect_name_input <- function(effect_name, genetic_ok = FALSE) {
  if (effect_name %in% GENETIC_EFFECT_NAMES_FUTURE) {
    stop("effect_name '", effect_name, "' is reserved but not yet supported: ",
         "no generator calibrates it. Write such effects by hand with ",
         "define_genome_effect_terms().", call. = FALSE)
  }
  if (effect_name %in% DERIVED_EFFECT_NAMES) {
    stop("effect_name '", effect_name, "' is reserved for derived output ",
         "(variance extraction) and cannot be defined.", call. = FALSE)
  }
  if (!genetic_ok && effect_name %in% GENETIC_EFFECT_NAMES) {
    stop("effect_name '", effect_name, "' is a reserved genetic variance ",
         "component (", paste0("'", GENETIC_EFFECT_NAMES, "'", collapse = ", "),
         ") and cannot name a phenotype-level effect. Choose another name.",
         call. = FALSE)
  }
  invisible(effect_name)
}


#' Check a covariance matrix's names against `trait_name`, never relabelling
#'
#' A single number is a 1 x 1 matrix when one name is given. A named matrix
#' (row names, and column names when present) must equal `trait_name` in
#' order; an unnamed one is taken in `trait_name` order. With no
#' `trait_name`, the row names are the names. Before 0.73.0 the names were
#' assigned over whatever the matrix carried, so a matrix named
#' `c("BF", "ADG")` passed with `trait_name = c("ADG", "BF")` was stored with
#' its rows the wrong way round.
#'
#' @return The matrix with `dimnames = list(names, names)`.
#' @keywords internal
.check_cov_dimnames <- function(x, trait_name, arg = "cov_matrix") {
  if (is.numeric(x) && is.null(dim(x)) && length(x) == 1L) {
    if (is.null(trait_name) || length(trait_name) != 1L) {
      stop("A single number for `", arg, "` needs exactly one `trait_name`.",
           call. = FALSE)
    }
    x <- matrix(x, 1L, 1L)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`", arg, "` must be a numeric matrix.", call. = FALSE)
  }
  n <- nrow(x)
  if (ncol(x) != n) stop("`", arg, "` must be square.", call. = FALSE)
  rn <- rownames(x); cn <- colnames(x)
  if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) {
    stop("`", arg, "` row names (", paste(rn, collapse = ", "),
         ") and column names (", paste(cn, collapse = ", "), ") differ.",
         call. = FALSE)
  }
  nm <- if (!is.null(rn)) rn else cn
  if (!is.null(trait_name)) {
    if (length(trait_name) != n) {
      stop("`trait_name` length (", length(trait_name),
           ") must equal matrix dimension (", n, ").", call. = FALSE)
    }
    if (!is.null(nm) && !identical(as.character(nm), as.character(trait_name))) {
      stop("`", arg, "` is named (", paste(nm, collapse = ", "),
           ") but `trait_name` is (", paste(trait_name, collapse = ", "),
           "). Names must match in order; the matrix is never relabelled.",
           call. = FALSE)
    }
    nm <- trait_name
  }
  if (is.null(nm) || any(is.na(nm)) || any(!nzchar(nm))) {
    stop("`", arg, "` must have row names, or supply `trait_name`.",
         call. = FALSE)
  }
  nm <- as.character(nm)
  lapply(nm, validate_sql_identifier, what = "trait name")
  if (anyDuplicated(nm)) {
    stop("`trait_name` must not contain duplicates.", call. = FALSE)
  }
  dimnames(x) <- list(nm, nm)
  x
}


# ---------------------------------------------------------------------------
# trait_var_comp: the one writer and the readers
# ---------------------------------------------------------------------------

#' SQL predicate for a `line_name` value, NULL-safe
#' @keywords internal
.tvc_line_sql <- function(conn, line_name) {
  if (is.null(line_name)) "line_name IS NULL"
  else paste0("line_name = ", DBI::dbQuoteLiteral(conn, line_name))
}

#' The traits of the stored block(s) touching `traits`
#'
#' A block is found from the rows themselves: the traits linked by
#' off-diagonal rows within one `effect_name` x `line_name`. No block id is
#' stored. Returns the connected closure of `traits`, restricted to traits
#' that have at least one stored row.
#' @keywords internal
.tvc_block_traits <- function(conn, effect_name, line_name, traits) {
  rows <- DBI::dbGetQuery(conn, paste0(
    "SELECT trait_name_1, trait_name_2 FROM trait_var_comp WHERE effect_name = ",
    DBI::dbQuoteLiteral(conn, effect_name), " AND ",
    .tvc_line_sql(conn, line_name)))
  stored <- unique(c(rows$trait_name_1, rows$trait_name_2))
  block  <- intersect(traits, stored)
  repeat {
    hit <- rows$trait_name_1 %in% block | rows$trait_name_2 %in% block
    grown <- union(block, c(rows$trait_name_1[hit], rows$trait_name_2[hit]))
    if (length(grown) == length(block)) break
    block <- grown
  }
  sort(block)
}

#' The `remove_rows()` call that clears one stored block
#' @keywords internal
.tvc_removal_call <- function(effect_name, line_name, block) {
  line_pred <- if (is.null(line_name)) "is.na(line_name)"
               else paste0('line_name == "', line_name, '"')
  paste0(
    '  get_table(pop, "trait_var_comp") |>\n',
    '    dplyr::filter(effect_name == "', effect_name, '", ', line_pred, ',\n',
    '                  trait_name_1 %in% c(',
    paste0('"', block, '"', collapse = ", "), ')) |>\n',
    '    remove_rows()')
}

#' Write one genetic covariance block to `trait_var_comp`
#'
#' The single write path for generation targets, shared by
#' [define_effect_cov_matrix()] and the generators' `G =`. It
#' * validates the block as finite, symmetric and positive semidefinite, on
#'   the correlation scale so the check does not depend on the traits' units;
#' * refuses when **any** row exists for `effect_name`, any of the block's
#'   traits and the same `line_name` (NULL-safe), checked on the whole table,
#'   even for an identical matrix;
#' * inserts all n^2 rows with `%.17g` literals, which round-trip a double
#'   exactly, through `dbExecute()` (never `dbWriteTable()`, which advances
#'   the RNG).
#'
#' It opens **no** transaction: the caller owns one, so a generator commits
#' the target together with its terms.
#'
#' @param G Named square matrix (`dimnames` = trait names).
#' @keywords internal
.tvc_write_block <- function(conn, effect_name, G, line_name = NULL) {
  traits <- rownames(G)
  what <- paste0("'", effect_name, "' covariance for ",
                 paste(traits, collapse = ", "))
  if (anyNA(G) || any(!is.finite(G))) {
    stop("The ", what, " must contain only finite values.", call. = FALSE)
  }
  .qtl_target_std(unname(G), name = paste0("the ", what))
  G <- (G + t(G)) / 2

  block <- .tvc_block_traits(conn, effect_name, line_name, traits)
  if (length(block)) {
    stop("A '", effect_name, "' block is already stored for ",
         paste(block, collapse = ", "),
         if (is.null(line_name)) " (population-wide)"
         else paste0(" (line '", line_name, "')"),
         ". Stored targets are never overwritten, even by an identical matrix. ",
         "To replace it, remove it first:\n",
         .tvc_removal_call(effect_name, line_name, block),
         "\nthen write the new one",
         if (effect_name == "additive") paste0(
           " with define_additive_effects(..., G = ), which writes the target ",
           "and re-draws the generated terms in one transaction, or with ",
           "define_effect_cov_matrix() if no generated terms exist at this ",
           "scope yet") else " with define_effect_cov_matrix()",
         ".", call. = FALSE)
  }

  n     <- length(traits)
  start <- next_int_id(conn, "trait_var_comp", "id_trait_var_comp")
  q     <- function(x) DBI::dbQuoteLiteral(conn, x)
  ln    <- if (is.null(line_name)) "NULL" else q(line_name)
  ij    <- expand.grid(j = seq_len(n), i = seq_len(n))
  vals  <- sprintf("(%d, %s, %s, %s, %s, %s)",
                   start + seq_len(n * n) - 1L, q(effect_name), ln,
                   vapply(traits[ij$i], q, character(1)),
                   vapply(traits[ij$j], q, character(1)),
                   sprintf("%.17g", G[cbind(ij$i, ij$j)]))
  DBI::dbExecute(conn, paste0(
    "INSERT INTO trait_var_comp (id_trait_var_comp, effect_name, line_name, ",
    "trait_name_1, trait_name_2, cov_value) VALUES ",
    paste(vals, collapse = ", ")))
  invisible(traits)
}

#' The traits with generated terms of one kind calibrated to one target scope
#'
#' A `line_name = "C"` target covers line-C terms. A `line_name = NULL` target
#' covers the population-wide terms (common or parent-only) **and** the
#' line-scoped terms of every line that has no stored block of its own for
#' that kind and trait: the generator resolves a line's target with the
#' `line -> NULL` fallback (`.tvc_resolve_line()`), so those terms were
#' calibrated to the population-wide target. A line block cannot be added
#' under existing line terms (the refusal below), so a line block present now
#' was present when its terms were generated.
#' @return Sorted character vector of trait names (empty when none).
#' @keywords internal
.tvc_generated_traits <- function(conn, effect_name, traits, line_name) {
  g <- .tvc_generated_terms(conn, effect_name, traits, line_name)
  sort(unique(g$model$terms$trait_name[g$hit]))
}

#' The generated terms calibrated to one target scope, term by term
#'
#' The work behind `.tvc_generated_traits()`, which see for the scope rule.
#' @return list(model = the `.gev_read_model()` result for the generated
#'   owner, hit = logical per term, line = each term's line or `NA`).
#' @keywords internal
.tvc_generated_terms <- function(conn, effect_name, traits, line_name) {
  model <- .gev_read_model(conn, traits, effect_owner = GE_GENERATED_OWNER)
  if (nrow(model$terms) == 0L) {
    return(list(model = model, hit = logical(0), line = character(0)))
  }
  kind <- .gev_target_kind(model)
  tl   <- .gev_term_line(model)
  tr   <- model$terms$trait_name
  at <- if (!is.null(line_name)) {
    !is.na(tl) & tl == line_name
  } else {
    is.na(tl) | vapply(seq_along(tl), function(i) {
      !is.na(tl[i]) && is.null(.tvc_resolve_line(conn, effect_name, tr[i], tl[i]))
    }, logical(1))
  }
  list(model = model, hit = !is.na(kind) & kind == effect_name & at, line = tl)
}

#' Refuse a genetic target under terms a generator calibrated (Q21)
#'
#' A `"generated"` term is calibrated to the target it was generated with, and
#' the prevalence threshold trusts the stored target for that reason. Writing a
#' target under such terms (even after removing the old one) would break that.
#' Only the exported [define_effect_cov_matrix()] calls this; a generator's
#' `G =` writes through `.tvc_write_block()` together with the terms it
#' calibrates, so it is never refused here.
#'
#' Scope as in `.tvc_generated_traits()`.
#' @keywords internal
.tvc_refuse_under_generated <- function(conn, effect_name, traits, line_name) {
  hit_traits <- .tvc_generated_traits(conn, effect_name, traits, line_name)
  if (length(hit_traits) == 0L) return(invisible(NULL))
  scope <- if (is.null(line_name)) {
    "population-wide, or scoped to a line with no target of its own"
  } else paste0("line '", line_name, "'")
  block <- .tvc_block_traits(conn, effect_name, line_name, traits)
  remove <- if (length(block)) paste0(
    "remove the stored block:\n",
    .tvc_removal_call(effect_name, line_name, block), "\nthen ")
  regenerate <- if (effect_name == "additive") paste0(
    "call define_additive_effects(..., G = ) with the new matrix",
    if (!is.null(line_name)) paste0(' and line_name = "', line_name, '"'),
    ", which writes the target and re-draws the terms in one transaction",
    if (is.null(line_name)) paste0(
      " (it re-draws one scope per call: re-run the others, each line and ",
      "parent_origin, afterwards so they read the new target)"),
    ".")
  else paste0("re-run define_genome_effects() with the new target (it ",
              "re-draws the trait's whole generated model).")
  stop("Trait(s) ", paste(hit_traits, collapse = ", "), " already have ",
       "generated '", effect_name, "' terms (", scope, "), calibrated to ",
       "the target they were generated with. A target written now would not ",
       "describe them, and define_phenotype(prevalence = ) trusts the stored ",
       "target. To change it, ", remove, regenerate, call. = FALSE)
}

#' The `line_name` whose rows a reader should use
#'
#' `NULL` reads the population-wide rows. A named line reads its own rows when
#' it has any for `effect_name` and these traits, and otherwise falls back to
#' the population-wide rows. The fallback is decided per `effect_name`.
#' @keywords internal
.tvc_resolve_line <- function(conn, effect_name, traits, line_name) {
  if (is.null(line_name)) return(NULL)
  n <- DBI::dbGetQuery(conn, paste0(
    "SELECT COUNT(*) AS n FROM trait_var_comp WHERE effect_name = ",
    DBI::dbQuoteLiteral(conn, effect_name), " AND ",
    .tvc_line_sql(conn, line_name), " AND trait_name_1 IN (",
    paste(DBI::dbQuoteLiteral(conn, traits), collapse = ", "), ")"))$n
  if (n > 0) line_name else NULL
}


#' Get the variance (diagonal) for one trait from trait_var_comp
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character.
#' @param trait_name Character.
#' @param line_name `NULL` (population-wide rows) or a line, which falls back
#'   to the population-wide rows when it has none of its own.
#' @return Numeric scalar, or `NA_real_` if not found.
#' @keywords internal
get_trait_var <- function(pop, effect_name, trait_name, line_name = NULL) {
  conn <- pop$db_conn
  line_name <- .tvc_resolve_line(conn, effect_name, trait_name, line_name)
  tn  <- DBI::dbQuoteLiteral(conn, trait_name)
  row <- DBI::dbGetQuery(conn, paste0(
    "SELECT cov_value FROM trait_var_comp WHERE effect_name = ",
    DBI::dbQuoteLiteral(conn, effect_name), " AND ",
    .tvc_line_sql(conn, line_name), " AND trait_name_1 = ", tn,
    " AND trait_name_2 = ", tn))
  if (nrow(row) == 0L) NA_real_ else row$cov_value[[1L]]
}


#' Get the variance (diagonal) for one phenotype from phenotype_var_comp
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character.
#' @param phenotype_name Character.
#' @return Numeric scalar, or `NA_real_` if not found.
#' @keywords internal
get_phenotype_var <- function(pop, effect_name, phenotype_name) {
  conn <- pop$db_conn
  pn   <- DBI::dbQuoteLiteral(conn, phenotype_name)
  row  <- DBI::dbGetQuery(conn, paste0(
    "SELECT cov_value FROM phenotype_var_comp ",
    "WHERE effect_name = ", DBI::dbQuoteLiteral(conn, effect_name), " ",
    "AND phenotype_name_1 = ", pn, " AND phenotype_name_2 = ", pn, " ",
    "AND condition_column IS NULL"))
  if (nrow(row) == 0L) NA_real_ else row$cov_value[[1L]]
}


#' Load a full covariance matrix from trait_var_comp
#'
#' @param pop A `tidybreed_pop` object.
#' @param effect_name Character.
#' @param trait_names Character vector of trait names.
#' @param line_name `NULL` (population-wide rows) or a line, which falls back
#'   to the population-wide rows when it has none of its own. Lines are never
#'   mixed.
#' @return Named numeric matrix, or `NULL` if any entry is missing.
#' @keywords internal
load_trait_cov <- function(pop, effect_name, trait_names, line_name = NULL) {
  conn <- pop$db_conn
  line_name <- .tvc_resolve_line(conn, effect_name, trait_names, line_name)
  n <- length(trait_names)
  R <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(trait_names, trait_names))
  tn <- paste(DBI::dbQuoteLiteral(conn, trait_names), collapse = ", ")
  rows <- DBI::dbGetQuery(conn, paste0(
    "SELECT trait_name_1, trait_name_2, cov_value FROM trait_var_comp ",
    "WHERE effect_name = ", DBI::dbQuoteLiteral(conn, effect_name), " AND ",
    .tvc_line_sql(conn, line_name), " AND trait_name_1 IN (", tn, ") ",
    "AND trait_name_2 IN (", tn, ")"))
  if (nrow(rows) == 0L) return(NULL)
  for (i in seq_len(nrow(rows))) {
    R[rows$trait_name_1[i], rows$trait_name_2[i]] <- rows$cov_value[i]
  }
  if (any(is.na(R))) return(NULL)
  R
}

#' Covariance blocks in `phenotype_var_comp`
#'
#' @description
#' A *covariance block* is a connected component of the graph whose vertices
#' are phenotypes and whose edges are stored pair rows in `phenotype_var_comp`
#' for one `effect_name` — a pair row with `cov_value = 0` is still an edge.
#' Every writer of that table goes through [write_phenotype_cov_block()], which
#' runs [validate_phenotype_cov_block()] inside its transaction before any row
#' is deleted, so a rejected declaration never leaves a half-written block.
#'
#' The rules, from `plans/sample_correlated_effects.md`:
#'
#' * **D1 — a block is declared in one call, as a complete matrix.** For a call
#'   over phenotypes `N`, let `U` be `N` plus every member of every existing
#'   block that touches `N`. The call is accepted only when `N == U`; a fragment
#'   or a strict subset of an existing block is an error naming the omitted
#'   phenotypes. The matrix must be symmetric, finite, and positive
#'   semi-definite.
#' * **Strata (residual only).** Conditional rows partition into strata
#'   `(condition_table, condition_column, condition_level)`. Every stratum of a
#'   block names the same phenotypes and a block has at most one condition
#'   column, so which stratum an individual resolves to never changes block
#'   membership.
#' * **D3 — realization lock.** A block cannot be redefined once a draw exists
#'   under it: for the residual effect, any `ind_phenotype` row of a member
#'   with `residual_value IS NOT NULL`; for a named effect, any
#'   `phenotype_random_effects` row for `(effect_name, member)`. The error
#'   gives the `remove_rows()` call that clears the realizations.
#' * **D6 — agreement.** Every phenotype in a residual block that has a
#'   `phenotype_meta` row carries the same `condition_change_action`.
#' * **Named effects (§5.6).** In a block of two or more phenotypes every
#'   `phenotype_effects` row for the effect is `random`, uses
#'   `distribution = "normal"`, and reads the same
#'   `(source_column, source_table)`. A 1 × 1 block is exempt, so a gamma or
#'   uniform effect on a single phenotype stays legal until something tries to
#'   join it to a second phenotype.
#'
#' @name phenotype_cov_block
#' @keywords internal
NULL


# ── Block discovery ─────────────────────────────────────────────────────────

.pvc_in_list <- function(conn, x) {
  paste(vapply(as.character(x), function(v) DBI::dbQuoteLiteral(conn, v),
               character(1)), collapse = ", ")
}

.pvc_set <- function(x) paste0("{", paste(x, collapse = ", "), "}")

#' Members of the block(s) touching a set of phenotypes
#'
#' The transitive closure over pair-row existence for `effect_name`, across
#' every stratum. Always contains `phenotype_names` itself, so a phenotype with
#' no stored rows is its own (empty) block.
#'
#' @param conn A DBI connection.
#' @param effect_name Character scalar.
#' @param phenotype_names Character vector.
#' @return Character vector of block members, sorted.
#' @keywords internal
.pvc_block_members <- function(conn, effect_name, phenotype_names) {
  eff_lit <- DBI::dbQuoteLiteral(conn, effect_name)
  members <- unique(as.character(phenotype_names))
  repeat {
    in_list <- .pvc_in_list(conn, members)
    rows <- DBI::dbGetQuery(conn, sprintf(
      "SELECT DISTINCT phenotype_name_1, phenotype_name_2 FROM phenotype_var_comp
       WHERE effect_name = %s
         AND (phenotype_name_1 IN (%s) OR phenotype_name_2 IN (%s))",
      eff_lit, in_list, in_list))
    found <- unique(c(rows$phenotype_name_1, rows$phenotype_name_2))
    added <- setdiff(found, members)
    if (length(added) == 0L) break
    members <- c(members, added)
  }
  sort(members)
}

.pvc_block_rows <- function(conn, effect_name, members) {
  in_list <- .pvc_in_list(conn, members)
  DBI::dbGetQuery(conn, sprintf(
    "SELECT phenotype_name_1, phenotype_name_2, cov_value,
            condition_column, condition_table, condition_level
     FROM phenotype_var_comp
     WHERE effect_name = %s
       AND phenotype_name_1 IN (%s) AND phenotype_name_2 IN (%s)",
    DBI::dbQuoteLiteral(conn, effect_name), in_list, in_list))
}

# One string per row identifying its stratum; "" is the unconditional stratum.
.pvc_stratum_key <- function(condition_table, condition_column, condition_level) {
  ifelse(is.na(condition_column) | !nzchar(condition_column), "",
         paste(condition_table, condition_column, condition_level, sep = "\r"))
}

.pvc_writer_call <- function(effect_name, members) {
  names_r <- paste0('c(', paste0('"', members, '"', collapse = ", "), ')')
  if (identical(effect_name, "residual")) {
    sprintf("define_residual_cov(pop, %s, R)", names_r)
  } else {
    sprintf('define_effect_cov_matrix(pop, "%s", R)   # dimnames(R) = %s',
            effect_name, names_r)
  }
}


# ── Checks ──────────────────────────────────────────────────────────────────

.pvc_check_matrix <- function(cov_matrix, phenotype_names, caller, tol = 1e-8) {
  n <- length(phenotype_names)
  if (!is.matrix(cov_matrix) || !is.numeric(cov_matrix)) {
    stop(caller, ": `cov_matrix` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(cov_matrix) != n || ncol(cov_matrix) != n) {
    stop(caller, ": `cov_matrix` must be ", n, " x ", n,
         " to match the phenotype names.", call. = FALSE)
  }
  rn <- rownames(cov_matrix)
  cn <- colnames(cov_matrix)
  if (is.null(rn) || is.null(cn) || !setequal(rn, phenotype_names) ||
      !setequal(cn, phenotype_names) || anyDuplicated(rn) || anyDuplicated(cn)) {
    stop(caller, ": `cov_matrix` dimnames must be exactly ",
         .pvc_set(phenotype_names), ".", call. = FALSE)
  }
  M <- cov_matrix[phenotype_names, phenotype_names, drop = FALSE]
  if (!all(is.finite(M))) {
    stop(caller, ": `cov_matrix` must be finite (no NA, NaN or Inf).",
         call. = FALSE)
  }
  scale <- max(1, max(abs(M)))
  if (n > 1L && max(abs(M - t(M))) > tol * scale) {
    stop(caller, ": `cov_matrix` must be symmetric (max discrepancy: ",
         format(max(abs(M - t(M)))), ").", call. = FALSE)
  }
  if (any(diag(M) < 0)) {
    stop(caller, ": diagonal entries (variances) must be non-negative.",
         call. = FALSE)
  }
  if (n > 1L) {
    ev <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -tol * max(max(ev), .Machine$double.eps)) {
      stop(caller, ": `cov_matrix` is not positive semi-definite ",
           "(smallest eigenvalue ", format(min(ev)), ").", call. = FALSE)
    }
  }
  invisible(M)
}

.pvc_check_complete <- function(effect_name, N, U, rows, caller) {
  omitted <- setdiff(U, N)
  if (length(omitted) == 0L) return(invisible(NULL))

  # For each omitted phenotype, the phenotypes it is already paired with.
  pairs <- vapply(omitted, function(o) {
    with_o <- unique(c(rows$phenotype_name_2[rows$phenotype_name_1 == o],
                       rows$phenotype_name_1[rows$phenotype_name_2 == o]))
    with_o <- sort(setdiff(with_o, o))
    if (length(with_o) == 0L) o
    else paste0(o, " (paired with ", paste(with_o, collapse = ", "), ")")
  }, character(1))

  label <- if (identical(effect_name, "residual")) "residual covariance block"
           else paste0("'", effect_name, "' covariance block")
  stop(
    caller, ": ", .pvc_set(N), " is not a complete ", label, ". ",
    "Stored rows already join ", if (length(N) == 1L) N else "these phenotypes",
    " to ", paste(omitted, collapse = ", "), ", so the block is ", .pvc_set(U),
    " and this call omits: ", paste(pairs, collapse = "; "), ".\n",
    "A block is declared in one call, as a complete matrix. Redeclare it over ",
    "all of ", .pvc_set(U), ":\n\n  ", .pvc_writer_call(effect_name, U), "\n\n",
    "with an explicit 0 for any pair that is uncorrelated.",
    call. = FALSE
  )
}

.pvc_check_strata <- function(N, rows, condition_column, condition_table,
                              condition_level, caller) {
  conditional <- !is.null(condition_column)

  # One (condition_table, condition_column) per block.
  cond_rows <- rows[!is.na(rows$condition_column) & nzchar(rows$condition_column), ,
                    drop = FALSE]
  existing_cols <- unique(paste(cond_rows$condition_table, cond_rows$condition_column,
                                sep = "."))
  new_col <- if (conditional) paste(condition_table, condition_column, sep = ".")
             else character(0)
  all_cols <- unique(c(existing_cols, new_col))
  if (length(all_cols) > 1L) {
    stop(
      caller, ": residual covariance block ", .pvc_set(N), " already has strata ",
      "on ", paste(existing_cols, collapse = ", "), "; a block has at most one ",
      "condition column (this call names ", new_col, "). Model a second ",
      "grouping factor as a named random effect, or clear the block and ",
      "redeclare every stratum on one column.",
      call. = FALSE
    )
  }

  # Every other stratum of the block names exactly N.
  this_key <- if (conditional) {
    .pvc_stratum_key(condition_table, condition_column, condition_level)
  } else ""
  keys <- .pvc_stratum_key(rows$condition_table, rows$condition_column,
                           rows$condition_level)
  for (k in setdiff(unique(keys), this_key)) {
    r   <- rows[keys == k, , drop = FALSE]
    set <- sort(unique(c(r$phenotype_name_1, r$phenotype_name_2)))
    if (!setequal(set, N)) {
      lbl <- if (k == "") "the unconditional stratum" else
        sprintf("stratum (%s.%s = '%s')", r$condition_table[1L],
                r$condition_column[1L], r$condition_level[1L])
      stop(
        caller, ": ", lbl, " of residual covariance block ", .pvc_set(N),
        " is declared over ", .pvc_set(set), " only; every stratum of a block ",
        "names the same phenotypes. Clear the block and redeclare each stratum ",
        "over ", .pvc_set(N), ":\n\n",
        "  pop |> get_table(\"phenotype_var_comp\") |>\n",
        "    filter(effect_name == \"residual\", phenotype_name_1 %in% c(",
        paste0('"', N, '"', collapse = ", "), ")) |>\n",
        "    remove_rows()\n",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

.pvc_check_realization_lock <- function(conn, effect_name, U, caller) {
  in_list <- .pvc_in_list(conn, U)
  names_r <- paste0('c(', paste0('"', U, '"', collapse = ", "), ')')

  if (identical(effect_name, "residual")) {
    n <- DBI::dbGetQuery(conn, sprintf(
      "SELECT COUNT(*) AS n FROM ind_phenotype
       WHERE phenotype_name IN (%s) AND residual_value IS NOT NULL", in_list))$n
    if (n == 0L) return(invisible(NULL))
    stop(
      caller, ": residual covariance block ", .pvc_set(U), " has ",
      format(n, big.mark = ","), " realized draw", if (n != 1L) "s",
      " in ind_phenotype. A block cannot be redefined after realization. ",
      "To redefine it, remove the realizations first:\n\n",
      "  pop |> get_table(\"ind_phenotype\") |>\n",
      "    filter(phenotype_name %in% ", names_r, ", !is.na(residual_value)) |>\n",
      "    remove_rows()\n\n",
      "then call define_residual_cov() again.",
      call. = FALSE
    )
  }

  n <- DBI::dbGetQuery(conn, sprintf(
    "SELECT COUNT(*) AS n FROM phenotype_random_effects
     WHERE effect_name = %s AND phenotype_name IN (%s)",
    DBI::dbQuoteLiteral(conn, effect_name), in_list))$n
  if (n == 0L) return(invisible(NULL))
  stop(
    caller, ": '", effect_name, "' covariance block ", .pvc_set(U), " has ",
    format(n, big.mark = ","), " realized draw", if (n != 1L) "s",
    " in phenotype_random_effects. A block cannot be redefined after ",
    "realization. To redefine it, remove the realizations first:\n\n",
    "  pop |> get_table(\"phenotype_random_effects\") |>\n",
    "    filter(effect_name == \"", effect_name, "\", phenotype_name %in% ",
    names_r, ") |>\n",
    "    remove_rows()\n\n",
    "then call define_effect_cov_matrix() again. Phenotype records already ",
    "computed with those draws are not removed by that call and will no longer ",
    "be explainable from stored state; remove them from ind_phenotype too if ",
    "the population must stay coherent.",
    call. = FALSE
  )
}


#' D6: `condition_change_action` agrees across a residual block
#'
#' @param conn A DBI connection.
#' @param phenotype_names The block members.
#' @param pending Optional `list(phenotype_name =, condition_change_action =)`
#'   for a `phenotype_meta` row about to be written; it replaces any stored
#'   row of the same name in the comparison.
#' @param caller Prefix for the error message.
#' @keywords internal
.check_condition_change_agreement <- function(conn, phenotype_names,
                                              pending = NULL,
                                              caller  = "define_residual_cov()") {
  if (length(phenotype_names) < 2L) return(invisible(NULL))
  meta <- DBI::dbGetQuery(conn, sprintf(
    "SELECT phenotype_name, condition_change_action FROM phenotype_meta
     WHERE phenotype_name IN (%s)", .pvc_in_list(conn, phenotype_names)))
  if (!is.null(pending)) {
    meta <- meta[meta$phenotype_name != pending$phenotype_name, , drop = FALSE]
    meta <- rbind(meta, data.frame(
      phenotype_name          = pending$phenotype_name,
      condition_change_action = pending$condition_change_action,
      stringsAsFactors = FALSE))
  }
  meta$condition_change_action[is.na(meta$condition_change_action)] <- "error"
  if (length(unique(meta$condition_change_action)) <= 1L) return(invisible(NULL))
  meta <- meta[order(meta$phenotype_name), , drop = FALSE]
  stop(
    caller, ": `condition_change_action` must agree across residual covariance ",
    "block ", .pvc_set(sort(phenotype_names)), ": ",
    paste0(meta$phenotype_name, " = '", meta$condition_change_action, "'",
           collapse = ", "),
    ". The value is block-scoped, so set it on the whole block at once with ",
    "define_condition_change_action(pop, '", sort(phenotype_names)[1], "', ",
    "'<action>'); define_phenotype() can only set it while the phenotype is ",
    "still a block of one.",
    call. = FALSE
  )
}


#' §5.6: the `phenotype_effects` rows of a named-effect block are compatible
#'
#' Applies to blocks of two or more phenotypes. Every row for
#' `(effect_name, member)` must be a `random` effect with
#' `distribution = "normal"` and the same `(source_column, source_table)`.
#'
#' @param conn A DBI connection.
#' @param effect_name Character scalar (not `"residual"`).
#' @param phenotype_names The block members.
#' @param pending Optional one-row data frame with columns `phenotype_name`,
#'   `effect_class`, `source_column`, `source_table`, `distribution` for a
#'   `phenotype_effects` row about to be written; it replaces any stored row
#'   of the same phenotype in the comparison.
#' @param caller Prefix for the error message.
#' @keywords internal
validate_named_effect_block <- function(conn, effect_name, phenotype_names,
                                        pending = NULL,
                                        caller  = "define_effect_cov_matrix()") {
  if (length(phenotype_names) < 2L) return(invisible(NULL))
  rows <- DBI::dbGetQuery(conn, sprintf(
    "SELECT phenotype_name, effect_class, source_column, source_table, distribution
     FROM phenotype_effects
     WHERE effect_name = %s AND phenotype_name IN (%s)",
    DBI::dbQuoteLiteral(conn, effect_name), .pvc_in_list(conn, phenotype_names)))
  if (!is.null(pending)) {
    rows <- rows[rows$phenotype_name != pending$phenotype_name, , drop = FALSE]
    rows <- rbind(rows, pending[, names(rows), drop = FALSE])
  }
  if (nrow(rows) == 0L) return(invisible(NULL))
  rows <- rows[order(rows$phenotype_name), , drop = FALSE]
  block_lbl <- paste0("'", effect_name, "' covariance block ",
                      .pvc_set(sort(phenotype_names)))

  not_random <- rows[rows$effect_class != "random", , drop = FALSE]
  if (nrow(not_random) > 0L) {
    stop(
      caller, ": ", block_lbl, " can only join random effects, but '",
      effect_name, "' is a ", not_random$effect_class[1L], " effect on ",
      not_random$phenotype_name[1L], ".",
      call. = FALSE
    )
  }

  dist <- rows$distribution
  dist[is.na(dist)] <- "normal"
  non_normal <- rows[dist != "normal", , drop = FALSE]
  if (nrow(non_normal) > 0L) {
    stop(
      caller, ": ", block_lbl, " requires distribution = \"normal\" on every ",
      "member; ", paste0(non_normal$phenotype_name, " uses \"",
                         dist[dist != "normal"], "\"", collapse = ", "),
      ". A non-normal named effect is supported only as a 1 x 1 block.",
      call. = FALSE
    )
  }

  src <- paste(rows$source_table, rows$source_column, sep = ".")
  if (length(unique(src)) > 1L) {
    stop(
      caller, ": ", block_lbl, " joins effects that read different grouping ",
      "columns: ", paste0(rows$phenotype_name, " reads ", src, collapse = ", "),
      ". Members of one block must share (source_column, source_table).",
      call. = FALSE
    )
  }
  invisible(NULL)
}


# ── The validator and the writer ────────────────────────────────────────────

#' Validate a `phenotype_var_comp` block declaration
#'
#' Runs D1 (completeness, matrix), the residual strata rules, D6, the §5.6
#' named-effect checks and the D3 realization lock against the *current* table
#' state — call it before deleting anything. Returns the block members
#' (`== phenotype_names` when it returns at all).
#'
#' @inheritParams write_phenotype_cov_block
#' @return The sorted block members, invisibly.
#' @keywords internal
validate_phenotype_cov_block <- function(conn, effect_name, phenotype_names,
                                         cov_matrix,
                                         condition_column = NULL,
                                         condition_table  = NULL,
                                         condition_level  = NULL,
                                         caller = "define_residual_cov()",
                                         tol    = 1e-8) {
  N <- sort(unique(as.character(phenotype_names)))
  if (length(N) != length(phenotype_names)) {
    stop(caller, ": `phenotype_names` must not contain duplicates.", call. = FALSE)
  }
  if (is.null(condition_column) != is.null(condition_level)) {
    stop(caller, ": supply both `condition_column` and `condition_level`, ",
         "or neither.", call. = FALSE)
  }
  if (!is.null(condition_column) && !identical(effect_name, "residual")) {
    stop(caller, ": condition strata are supported for the residual effect only.",
         call. = FALSE)
  }

  .pvc_check_matrix(cov_matrix, N, caller, tol = tol)

  U    <- .pvc_block_members(conn, effect_name, N)
  rows <- .pvc_block_rows(conn, effect_name, U)
  .pvc_check_complete(effect_name, N, U, rows, caller)

  if (identical(effect_name, "residual")) {
    .pvc_check_strata(N, rows, condition_column, condition_table,
                      condition_level, caller)
    .check_condition_change_agreement(conn, U, caller = caller)
  } else {
    validate_named_effect_block(conn, effect_name, U, caller = caller)
  }

  .pvc_check_realization_lock(conn, effect_name, U, caller)
  invisible(U)
}


#' Replace one stratum of a block: validate, delete, insert
#'
#' No transaction management — the caller owns the transaction. Writes through
#' `duckdb_register()` + `INSERT`, which does not touch R's RNG.
#'
#' @inheritParams write_phenotype_cov_block
#' @keywords internal
.pvc_write_block <- function(conn, effect_name, phenotype_names, cov_matrix,
                             condition_column = NULL, condition_table = NULL,
                             condition_level = NULL,
                             caller = "define_residual_cov()", tol = 1e-8) {
  N <- validate_phenotype_cov_block(conn, effect_name, phenotype_names,
                                    cov_matrix, condition_column,
                                    condition_table, condition_level,
                                    caller = caller, tol = tol)
  # Symmetrize so the stored (i, j) and (j, i) rows are bit-identical (exact
  # for an already-symmetric input); find_covariance_blocks() relies on it.
  M <- cov_matrix[N, N, drop = FALSE]
  M <- (M + t(M)) / 2
  n <- length(N)

  eff_lit <- DBI::dbQuoteLiteral(conn, effect_name)
  in_list <- .pvc_in_list(conn, N)
  if (is.null(condition_column)) {
    stratum_pred <- "condition_column IS NULL"
  } else {
    stratum_pred <- sprintf(
      "condition_table = %s AND condition_column = %s AND condition_level = %s",
      DBI::dbQuoteLiteral(conn, condition_table),
      DBI::dbQuoteLiteral(conn, condition_column),
      DBI::dbQuoteLiteral(conn, as.character(condition_level)))
  }
  DBI::dbExecute(conn, sprintf(
    "DELETE FROM phenotype_var_comp
     WHERE effect_name = %s AND %s
       AND (phenotype_name_1 IN (%s) OR phenotype_name_2 IN (%s))",
    eff_lit, stratum_pred, in_list, in_list))

  first_id <- next_int_id(conn, "phenotype_var_comp", "id_phenotype_var_comp")
  rows <- data.frame(
    id_phenotype_var_comp = as.integer(first_id + seq_len(n * n) - 1L),
    effect_name           = effect_name,
    phenotype_name_1      = rep(N, each  = n),
    phenotype_name_2      = rep(N, times = n),
    cov_value             = as.numeric(t(M)),
    condition_column      = if (is.null(condition_column)) NA_character_
                            else condition_column,
    condition_table       = if (is.null(condition_column)) NA_character_
                            else condition_table,
    condition_level       = if (is.null(condition_level)) NA_character_
                            else as.character(condition_level),
    stringsAsFactors = FALSE
  )

  tmp <- "__phenotype_cov_block_tmp"
  duckdb::duckdb_register(conn, tmp, rows)
  on.exit(try(duckdb::duckdb_unregister(conn, tmp), silent = TRUE), add = TRUE)
  DBI::dbExecute(conn, sprintf(
    "INSERT INTO phenotype_var_comp
       (id_phenotype_var_comp, effect_name, phenotype_name_1, phenotype_name_2,
        cov_value, condition_column, condition_table, condition_level)
     SELECT id_phenotype_var_comp, effect_name, phenotype_name_1, phenotype_name_2,
            cov_value, condition_column, condition_table, condition_level
     FROM %s", tmp))
  invisible(N)
}


#' Write one stratum of a covariance block in its own transaction
#'
#' The single write path for `phenotype_var_comp`, used by
#' [define_residual_cov()], [define_effect_cov_matrix()] and (inside its own
#' transaction, via `.pvc_write_block()`) [define_effect_random()].
#'
#' @param conn A DBI connection.
#' @param effect_name `"residual"` or a named random effect.
#' @param phenotype_names Character vector; the block being declared.
#' @param cov_matrix Numeric matrix with `dimnames` equal to `phenotype_names`.
#' @param condition_column,condition_table,condition_level The stratum
#'   (`NULL`s for the unconditional stratum). Residual only.
#' @param caller Prefix for error messages.
#' @param tol Relative symmetry / PSD tolerance.
#' @return The sorted block members, invisibly.
#' @keywords internal
write_phenotype_cov_block <- function(conn, effect_name, phenotype_names,
                                      cov_matrix,
                                      condition_column = NULL,
                                      condition_table  = NULL,
                                      condition_level  = NULL,
                                      caller = "define_residual_cov()",
                                      tol    = 1e-8) {
  DBI::dbExecute(conn, "BEGIN TRANSACTION")
  committed <- FALSE
  on.exit(if (!committed) try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE),
          add = TRUE)
  N <- .pvc_write_block(conn, effect_name, phenotype_names, cov_matrix,
                        condition_column, condition_table, condition_level,
                        caller = caller, tol = tol)
  DBI::dbExecute(conn, "COMMIT")
  committed <- TRUE
  invisible(N)
}

# ============================================================================
# formula_helpers.R
#
# Internal helpers for formula-based phenotype specification:
#   - formula_tgv : DSL string for composite genetic-value assembly
#   - formula     : arithmetic derivation over existing ind_phenotype records
# ============================================================================

# ── Constants ────────────────────────────────────────────────────────────────

.FORMULA_MATH_WHITELIST <- c(
  "sqrt", "log", "log2", "log10", "exp", "abs",
  "round", "ceiling", "floor", "sign", "trunc",
  "sin", "cos", "tan", "asin", "acos", "atan"
)

.FORMULA_TGV_DSL_FUNS <- c("self", "dam", "sire", "group_sum", "group_mean")

.FORMULA_ARITH_OPS <- c("+", "-", "*", "/", "^", "(", ")")


# ── Generic AST utilities ────────────────────────────────────────────────────

#' Recursively extract all name nodes from a parsed R expression.
#' Skips the call head (function name) and numeric/logical literals.
#' Used by derived-formula validation and topo sort.
#' @keywords internal
.extract_all_symbols <- function(e) {
  if (is.name(e))   return(as.character(e))
  if (is.call(e)) {
    args <- as.list(e)[-1]  # drop function name slot
    return(unlist(lapply(args, .extract_all_symbols), use.names = FALSE))
  }
  character(0)
}


# ── formula_tgv validation ────────────────────────────────────────────────────

#' Validate a formula_tgv string at define_phenotype() time.
#'
#' Parses the DSL formula and walks it with `.walk_formula_tgv_ast()`, which
#' refuses any call outside the DSL, the arithmetic operators and the math
#' whitelist, and any DSL call with an argument it does not take. Then checks
#' every referenced trait against `trait_meta` (with close-match suggestions
#' via `agrep()`), and every group `table` / `col` against the database: the
#' table must exist and hold `id_ind` and the column.
#'
#' @param conn DBI connection.
#' @param formula_tgv Character. The DSL formula string.
#' @return Invisible NULL on success. Stops on error.
#' @keywords internal
.validate_formula_tgv <- function(conn, formula_tgv) {
  walk_res <- .walk_formula_tgv_ast(.parse_formula_tgv(formula_tgv))

  known_traits <- DBI::dbGetQuery(conn, "SELECT trait_name FROM trait_meta")$trait_name
  ref_traits <- unique(vapply(walk_res$trait_refs, `[[`, character(1), "trait"))
  unknown    <- setdiff(ref_traits, known_traits)

  if (length(unknown) > 0) {
    suggestions <- vapply(unknown, function(u) {
      close <- agrep(u, known_traits, ignore.case = TRUE, value = TRUE,
                     max.distance = list(cost = 2, all = 2))
      if (length(close) > 0)
        paste0("  '", u, "' → did you mean: ", paste(head(close, 3), collapse = ", "), "?")
      else
        paste0("  '", u, "' → not found in trait_meta")
    }, character(1))
    stop(
      "Unknown trait name(s) in `formula_tgv = \"", formula_tgv, "\"`:\n",
      paste(suggestions, collapse = "\n"), "\n",
      "All symbols must be trait names defined with define_trait() in trait_meta.",
      call. = FALSE
    )
  }

  # Group lookups: the table and column must exist now, not at add_phenotype().
  tables <- DBI::dbListTables(conn)
  for (ref in walk_res$trait_refs) {
    if (!ref$type %in% c("group_sum", "group_mean")) next
    where <- paste0("`formula_tgv = \"", formula_tgv, "\"`, ", ref$call)
    if (!ref$table %in% tables) {
      stop(where, ": table '", ref$table, "' does not exist.",
           .case_hint(ref$table, tables), call. = FALSE)
    }
    fields <- DBI::dbListFields(conn, ref$table)
    if (!"id_ind" %in% fields) {
      stop(where, ": table '", ref$table, "' has no 'id_ind' column, so it ",
           "cannot give each individual a group.", call. = FALSE)
    }
    if (!ref$col %in% fields) {
      stop(where, ": column '", ref$col, "' not found in table '", ref$table,
           "'.", .case_hint(ref$col, fields),
           " Add it (e.g. with mutate_table()) before define_phenotype().",
           call. = FALSE)
    }
  }

  invisible(NULL)
}


# " Did you mean 'x'?" when `name` matches one of `choices` only up to case,
# else "". Names are matched exactly everywhere in tidybreed, although DuckDB
# itself ignores case, so a case slip gets a hint rather than a silent match.
.case_hint <- function(name, choices) {
  hit <- choices[tolower(choices) == tolower(name)]
  if (length(hit) == 0L) return("")
  paste0(" Did you mean '", hit[1], "'? Names are case-sensitive.")
}


# Parse a formula_tgv string to one expression.
.parse_formula_tgv <- function(formula_tgv) {
  exprs <- tryCatch(
    parse(text = formula_tgv, keep.source = FALSE),
    error = function(e)
      stop("Could not parse `formula_tgv = \"", formula_tgv, "\"`: ",
           conditionMessage(e), call. = FALSE)
  )
  if (length(exprs) != 1L) {
    stop("`formula_tgv = \"", formula_tgv, "\"` must be a single expression.",
         call. = FALSE)
  }
  exprs[[1]]
}


# ── derived formula validation ────────────────────────────────────────────────

#' Check a derived formula's grammar and return the phenotypes it reads
#'
#' The expression is `eval()`ed at `add_phenotype()` time, so only phenotype
#' names, numbers, the arithmetic operators and the math whitelist are
#' accepted. Any other call or constant is an error naming it.
#'
#' @param expr Parsed R expression.
#' @param formula The formula string, for messages.
#' @return Character vector of the phenotype names referenced (unique).
#' @keywords internal
.check_derived_formula <- function(expr, formula) {
  allowed <- c(.FORMULA_ARITH_OPS, .FORMULA_MATH_WHITELIST)
  where <- paste0("`formula = \"", formula, "\"`: ")
  syms <- character(0)
  walk <- function(e) {
    if (is.numeric(e)) return(invisible())
    if (is.name(e)) {
      nm <- as.character(e)
      if (!nzchar(nm) || nm %in% allowed) {
        stop(where, "`", nm, "` is not a phenotype name.", call. = FALSE)
      }
      syms <<- c(syms, nm)
      return(invisible())
    }
    if (is.call(e)) {
      head <- e[[1]]
      fn <- if (is.name(head)) as.character(head) else ""
      if (!fn %in% allowed) {
        stop(where, "`", paste(deparse(e, width.cutoff = 500L), collapse = " "),
             "` is not allowed. Use phenotype names, numbers, the operators ",
             "+ - * / ^ and the math functions ",
             paste(.FORMULA_MATH_WHITELIST, collapse = ", "), ".",
             call. = FALSE)
      }
      for (child in as.list(e)[-1]) walk(child)
      return(invisible())
    }
    stop(where, "the constant `", paste(deparse(e), collapse = " "),
         "` is not allowed; only numbers and phenotype names may appear.",
         call. = FALSE)
  }
  walk(expr)
  unique(syms)
}


# Parse a derived formula string to one expression.
.parse_derived_formula <- function(formula) {
  exprs <- tryCatch(
    parse(text = formula, keep.source = FALSE),
    error = function(e)
      stop("Could not parse `formula = \"", formula, "\"`: ",
           conditionMessage(e), call. = FALSE)
  )
  if (length(exprs) != 1L) {
    stop("`formula = \"", formula, "\"` must be a single expression.",
         call. = FALSE)
  }
  exprs[[1]]
}


#' Validate a derived formula string at define_phenotype() time.
#'
#' The grammar is checked strictly (`.check_derived_formula()`). Symbols that
#' are not in phenotype_meta generate a warning (not an error), allowing
#' config-first workflows where components are defined before their
#' dependents. A hard error at add_phenotype() time fires if still missing.
#'
#' @param formula Character. The arithmetic formula string.
#' @param known_phenos Character vector of phenotype names from phenotype_meta.
#' @return Invisible NULL on success.
#' @keywords internal
.validate_derived_formula <- function(formula, known_phenos) {
  symbols <- .check_derived_formula(.parse_derived_formula(formula), formula)

  unknown <- setdiff(symbols, known_phenos)
  if (length(unknown) > 0) {
    suggestions <- vapply(unknown, function(u) {
      close <- agrep(u, known_phenos, ignore.case = TRUE, value = TRUE,
                     max.distance = list(cost = 2, all = 2))
      if (length(close) > 0)
        paste0("  '", u, "' → did you mean: ", paste(head(close, 3), collapse = ", "), "?")
      else
        paste0("  '", u, "' → not found in phenotype_meta")
    }, character(1))
    warning(
      "Unknown phenotype name(s) in `formula = \"", formula, "\"`:\n",
      paste(suggestions, collapse = "\n"), "\n",
      "These phenotypes must be defined with define_phenotype() before add_phenotype() is called.",
      call. = FALSE
    )
  }

  invisible(NULL)
}


# ── formula_tgv AST walk ──────────────────────────────────────────────────────

# The arguments each DSL function takes: the positional ones, in order, and
# the optional named ones. Nothing else is accepted.
.FORMULA_TGV_DSL_ARGS <- list(
  self       = list(positional = "trait",          named = "component"),
  dam        = list(positional = "trait",          named = "component"),
  sire       = list(positional = "trait",          named = "component"),
  group_sum  = list(positional = c("trait", "col"), named = c("component", "table")),
  group_mean = list(positional = c("trait", "col"), named = c("component", "table"))
)

#' Walk a formula_tgv expression: validate it, collect its references, and
#' replace each with a placeholder
#'
#' One depth-first pass. Each bare trait symbol or DSL call becomes one
#' reference and is replaced, in the returned expression, by that reference's
#' unique placeholder symbol (`.tgv_1`, `.tgv_2`, ...), so the same trait can
#' appear several times with different contributors, components or tables.
#'
#' The DSL: a bare symbol is `self(trait)`; `self(trait)`, `dam(trait)` and
#' `sire(trait)` take one positional trait; `group_sum(trait, col)` and
#' `group_mean(trait, col)` take a trait and a group column. All five take an
#' optional named `component =` (one of [TGV_COMPONENT_NAMES] or `"total"`,
#' the default); the group functions also take a named `table =` (default
#' `"ind_meta"`). Trait, column and table are symbols or strings; `col` and
#' `table` must be SQL identifiers. Any other argument, any call outside the
#' DSL, the arithmetic operators and the math whitelist, and any constant
#' that is not a number is an error naming the offending call.
#'
#' @param expr Parsed R expression (from `.parse_formula_tgv()`).
#' @return A list:
#'   $trait_refs: list of lists, each with:
#'     - trait:       character trait name
#'     - type:        "self", "dam", "sire", "group_sum", or "group_mean"
#'     - col:         group column name (NA for non-group types)
#'     - table:       group table name  (NA for non-group types)
#'     - component:   "total" or one ind_tgv.component_name
#'     - placeholder: unique R symbol name for the pre-fetched vector
#'     - call:        the reference as written, for messages
#'   $expr: `expr` with every reference replaced by its placeholder
#' @keywords internal
.walk_formula_tgv_ast <- function(expr) {
  trait_refs       <- list()
  allowed_calls    <- c(.FORMULA_ARITH_OPS, .FORMULA_MATH_WHITELIST)

  add_ref <- function(trait, type, col, table, component, call) {
    ph <- paste0(".tgv_", length(trait_refs) + 1L)
    trait_refs[[length(trait_refs) + 1L]] <<- list(
      trait = trait, type = type, col = col, table = table,
      component = component, placeholder = ph, call = call)
    as.name(ph)
  }

  dsl_call <- function(e, fn) {
    shown <- paste(deparse(e, width.cutoff = 500L), collapse = " ")
    bad <- function(...) {
      stop("formula_tgv: ", shown, ": ", ..., call. = FALSE)
    }
    spec <- .FORMULA_TGV_DSL_ARGS[[fn]]
    args <- as.list(e)[-1]
    nms  <- names(args)
    if (is.null(nms)) nms <- rep("", length(args))
    pos  <- args[nms == ""]
    named <- args[nms != ""]
    if (length(pos) != length(spec$positional)) {
      bad(fn, "() takes ", length(spec$positional), " positional argument",
          if (length(spec$positional) > 1L) "s", " (",
          paste(spec$positional, collapse = ", "), "), not ", length(pos),
          if (length(spec$named) > 0L)
            paste0("; ", paste0("`", spec$named, " =`", collapse = " and "),
                   " must be named"), ".")
    }
    unknown <- setdiff(names(named), spec$named)
    if (length(unknown) > 0L) {
      bad(fn, "() has no argument `", unknown[1], "`; it takes ",
          paste0("`", spec$named, " =`", collapse = " and "), ".")
    }
    if (anyDuplicated(names(named))) {
      bad("argument `", names(named)[anyDuplicated(names(named))],
          "` is given twice.")
    }
    value <- function(x, what) {
      if (is.name(x)) return(as.character(x))
      if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) return(x)
      bad("`", what, "` must be a name or a single string.")
    }
    trait <- value(pos[[1]], "trait")
    col   <- if (length(pos) >= 2L) value(pos[[2]], "col") else NA_character_
    component <- if (is.null(named$component)) "total"
                 else value(named$component, "component")
    if (!component %in% c(TGV_COMPONENT_NAMES, "total")) {
      bad("`component = \"", component, "\"` must be one of ",
          paste0("'", c(TGV_COMPONENT_NAMES, "total"), "'", collapse = ", "), ".")
    }
    table <- NA_character_
    if (fn %in% c("group_sum", "group_mean")) {
      table <- if (is.null(named$table)) "ind_meta" else value(named$table, "table")
      tryCatch({
        validate_sql_identifier(col,   what = "group column")
        validate_sql_identifier(table, what = "group table")
      }, error = function(err) bad(conditionMessage(err)))
    }
    add_ref(trait, if (fn == "self") "self" else fn, col, table, component, shown)
  }

  walk <- function(e) {
    # Numbers are ordinary weights and offsets (0.5 * dam(WWM)).
    if (is.numeric(e)) return(e)
    if (is.name(e)) {
      nm <- as.character(e)
      if (nm %in% c(allowed_calls, .FORMULA_TGV_DSL_FUNS)) {
        stop("formula_tgv: `", nm, "` is a function name, not a trait.",
             call. = FALSE)
      }
      return(add_ref(nm, "self", NA_character_, NA_character_, "total", nm))
    }
    if (is.call(e)) {
      head <- e[[1]]
      fn <- if (is.name(head)) as.character(head) else ""
      if (fn %in% .FORMULA_TGV_DSL_FUNS) return(dsl_call(e, fn))
      if (!fn %in% allowed_calls) {
        stop("formula_tgv: `",
             paste(deparse(e, width.cutoff = 500L), collapse = " "),
             "` is not allowed. Use the contributor functions ",
             paste0(.FORMULA_TGV_DSL_FUNS, "()", collapse = ", "),
             ", the operators + - * / ^ and the math functions ",
             paste(.FORMULA_MATH_WHITELIST, collapse = ", "), ".",
             call. = FALSE)
      }
      args <- lapply(as.list(e)[-1], walk)
      return(as.call(c(list(head), args)))
    }
    stop("formula_tgv: the constant `",
         paste(deparse(e), collapse = " "), "` is not allowed; only numbers ",
         "may appear outside a contributor function.", call. = FALSE)
  }

  new_expr <- walk(expr)
  list(trait_refs = trait_refs, expr = new_expr)
}


# ── Genetic-value pre-fetching ────────────────────────────────────────────────

#' Pre-fetch every genetic-value vector a `formula_tgv` expression needs
#'
#' One contributor lookup per reference (see `?contributor_tgv`), reading
#' the reference's `component` (`"total"` by default), returned as a named
#' list ready to be the `eval()` environment.
#'
#' @param trait_refs List from `.walk_formula_tgv_ast()$trait_refs`.
#' @param subset_df The planned `ind_meta` rows.
#' @return Named list: placeholder -> numeric vector (`NA` = missing piece).
#' @keywords internal
.build_tgv_env <- function(pop, trait_refs, subset_df, phenotype_name) {
  conn      <- pop$db_conn
  focal_ids <- as.character(subset_df$id_ind)
  what      <- paste0("formula_tgv for phenotype '", phenotype_name, "'")
  env_list  <- list()
  for (ref in trait_refs) {
    env_list[[ref$placeholder]] <- switch(
      ref$type,
      self       = .tgv_by_id(conn, ref$trait, focal_ids, ref$component),
      dam        = .tgv_by_id(conn, ref$trait, subset_df$id_parent_2,
                              ref$component),
      sire       = .tgv_by_id(conn, ref$trait, subset_df$id_parent_1,
                              ref$component),
      group_sum  = .group_mate_tgv(conn, ref$trait, focal_ids, ref$col,
                                   ref$table, "sum", what, ref$component),
      group_mean = .group_mate_tgv(conn, ref$trait, focal_ids, ref$col,
                                   ref$table, "mean", what, ref$component))
  }
  env_list
}


# ── Top-level formula_tgv evaluator ──────────────────────────────────────────

#' Evaluate a formula_tgv string for a set of individuals.
#'
#' Orchestrates: parse → AST walk (references replaced by placeholders) →
#' genetic-value pre-fetch → eval().
#'
#' @param pop A tidybreed_pop object.
#' @param formula_tgv Character. DSL formula string from phenotype_meta.
#' @param subset_df Data frame: sex-filtered ind_meta rows.
#' @param phenotype_name Character. Used in error messages.
#' @return Named numeric vector (names = id_ind). NA marks excluded individuals
#'         (a missing dam/sire genetic value, or NA group membership).
#' @keywords internal
.eval_formula_tgv <- function(pop, formula_tgv, subset_df, phenotype_name) {
  walk_res <- .walk_formula_tgv_ast(.parse_formula_tgv(formula_tgv))
  tgv_env  <- .build_tgv_env(pop, walk_res$trait_refs, subset_df, phenotype_name)
  result   <- eval(walk_res$expr, envir = list2env(tgv_env, parent = baseenv()))
  stats::setNames(as.numeric(result), as.character(subset_df$id_ind))
}


# ── Derived formula pivot helper ─────────────────────────────────────────────

#' Pivot ind_phenotype rows wide for formula evaluation.
#'
#' Uses pheno_number = 1 by default. If, across the entire `rows` input, not
#' a single row has `pheno_number == 1` (e.g. first records were deleted),
#' falls back to the minimum `pheno_number` per (id_ind, phenotype_name) pair
#' instead. This fallback is all-or-nothing over the whole input, not decided
#' per phenotype or per individual.
#'
#' @param rows Data frame from ind_phenotype query (id_ind, phenotype_name,
#'   pheno_value, pheno_number).
#' @param ids Character vector. Ordered individual IDs.
#' @param pheno_names Character vector. Phenotype column names to create.
#' @return Data frame with one row per id in ids, columns = pheno_names.
#' @keywords internal
.pivot_pheno_wide <- function(rows, ids, pheno_names) {
  if (nrow(rows) == 0) {
    df <- as.data.frame(
      stats::setNames(
        replicate(length(pheno_names), rep(NA_real_, length(ids)), simplify = FALSE),
        pheno_names
      ),
      stringsAsFactors = FALSE
    )
    return(df)
  }

  rows1 <- rows[rows$pheno_number == 1L, , drop = FALSE]
  if (nrow(rows1) == 0) {
    # Fall back to min pheno_number per individual × phenotype
    rows1 <- do.call(rbind, lapply(split(rows, interaction(rows$id_ind, rows$phenotype_name, drop = TRUE)),
      function(sub) sub[which.min(sub$pheno_number), , drop = FALSE]
    ))
  }

  wide <- data.frame(id_ind = ids, stringsAsFactors = FALSE)
  for (pn in pheno_names) {
    sub  <- rows1[rows1$phenotype_name == pn, c("id_ind", "pheno_value"), drop = FALSE]
    names(sub)[2] <- pn
    wide <- merge(wide, sub, by = "id_ind", all.x = TRUE)
  }
  wide[match(ids, wide$id_ind), pheno_names, drop = FALSE]
}


# ── Top-level derived formula evaluator ──────────────────────────────────────

#' Evaluate a derived formula over ind_phenotype records.
#'
#' @param pop A tidybreed_pop object.
#' @param formula Character. Formula string from phenotype_meta.
#' @param ids Character vector. Individual IDs to compute for.
#' @param phenotype_name Character. Name of the derived phenotype (for messages).
#' @param pending Optional data frame of records planned earlier in the same
#'   `add_phenotype()` call (`id_ind`, `phenotype_name`, `pheno_value`,
#'   `pheno_number`) that are not on disk yet. They are treated exactly as if
#'   they had been written, so a derived phenotype can consume a feeder
#'   phenotype from the same call.
#' @return Numeric vector, same length as ids.
#'         NA propagates naturally; Inf/NaN converted to NA with a warning.
#' @keywords internal
.eval_derived_formula <- function(pop, formula, ids, phenotype_name,
                                  pending = NULL) {
  # Re-checked here: the stored string is eval()ed below.
  expr    <- .parse_derived_formula(formula)
  symbols <- .check_derived_formula(expr, formula)

  if (length(ids) == 0) return(numeric(0))

  conn <- pop$db_conn
  tmp  <- "__ap_derived_ids"
  duckdb::duckdb_register(conn, tmp, data.frame(id_ind = unique(ids),
                                                stringsAsFactors = FALSE))
  on.exit(try(duckdb::duckdb_unregister(conn, tmp), silent = TRUE), add = TRUE)
  rows <- DBI::dbGetQuery(conn, paste0(
    "SELECT p.id_ind, p.phenotype_name, p.pheno_value, p.pheno_number ",
    "FROM ind_phenotype AS p JOIN ", tmp, " AS f USING (id_ind) ",
    "WHERE p.phenotype_name IN (", .pvc_in_list(conn, symbols), ")"))

  if (!is.null(pending) && nrow(pending) > 0 &&
      all(c("id_ind", "phenotype_name", "pheno_value", "pheno_number") %in%
          names(pending))) {
    keep <- pending$phenotype_name %in% symbols & pending$id_ind %in% ids
    rows <- rbind(rows, as.data.frame(
      pending[keep, c("id_ind", "phenotype_name", "pheno_value", "pheno_number")]))
  }

  if (nrow(rows) == 0)
    stop(
      "No phenotype records found for component phenotypes (",
      paste(symbols, collapse = ", "), ") needed for derived phenotype '",
      phenotype_name, "'. ",
      "Ensure component phenotypes are simulated with add_phenotype() before '",
      phenotype_name, "'.",
      call. = FALSE
    )

  wide_df <- .pivot_pheno_wide(rows, ids, symbols)

  result <- tryCatch(
    eval(expr, envir = wide_df, enclos = baseenv()),
    error = function(e)
      stop("Error evaluating formula '", formula, "' for phenotype '",
           phenotype_name, "': ", conditionMessage(e), call. = FALSE)
  )

  n_bad <- sum(is.infinite(result) | is.nan(result), na.rm = TRUE)
  if (n_bad > 0) {
    bad_ids <- ids[is.infinite(result) | is.nan(result)]
    warning(
      n_bad, " non-finite value(s) produced by formula '", formula,
      "' for phenotype '", phenotype_name, "'. Converting to NA. ",
      "Example individual IDs: ", paste(head(bad_ids, 5), collapse = ", "),
      if (length(bad_ids) > 5) paste0(" ... +", length(bad_ids) - 5, " more") else "",
      call. = FALSE
    )
    result[is.infinite(result) | is.nan(result)] <- NA_real_
  }

  as.numeric(result)
}


# ── Topological sort ──────────────────────────────────────────────────────────

#' Topologically sort phenotypes for safe evaluation order.
#'
#' Derived formula phenotypes depend on other phenotypes being present in
#' ind_phenotype first. formula_tgv and components phenotypes have no
#' inter-phenotype dependencies (they depend on trait_meta, not phenotype_meta).
#'
#' Uses Kahn's BFS algorithm. Detects cycles and stops with an informative error.
#'
#' @param pheno_meta Data frame from phenotype_meta (all columns).
#' @return Character vector: phenotype_name in safe evaluation order.
#' @keywords internal
.topo_sort_phenotypes <- function(pheno_meta) {
  phenos  <- pheno_meta$phenotype_name
  n       <- length(phenos)
  idx_map <- stats::setNames(seq_len(n), phenos)

  # dep_of[i] = integer indices of phenotypes that i depends on
  dep_of <- vector("list", n)

  for (i in seq_len(n)) {
    formula_col <- pheno_meta$formula[[i]]
    if (is.na(formula_col) || !nzchar(formula_col)) next

    expr <- tryCatch(
      parse(text = formula_col, keep.source = FALSE)[[1]],
      error = function(e) NULL
    )
    if (is.null(expr)) next

    syms <- .extract_all_symbols(expr)
    syms <- setdiff(syms, c(.FORMULA_MATH_WHITELIST, .FORMULA_ARITH_OPS))
    deps <- intersect(syms, phenos)
    dep_of[[i]] <- unname(idx_map[deps])
  }

  # Build successor list (succ[d] = indices that must come after d)
  succ   <- vector("list", n)
  in_deg <- integer(n)
  for (i in seq_len(n)) {
    for (d in dep_of[[i]]) {
      succ[[d]] <- c(succ[[d]], i)
      in_deg[i] <- in_deg[i] + 1L
    }
  }

  # Kahn's BFS
  queue  <- which(in_deg == 0L)
  result <- integer(0)
  while (length(queue) > 0) {
    curr   <- queue[1]
    queue  <- queue[-1]
    result <- c(result, curr)
    for (s in succ[[curr]]) {
      in_deg[s] <- in_deg[s] - 1L
      if (in_deg[s] == 0L) queue <- c(queue, s)
    }
  }

  if (length(result) < n) {
    cycle_phenos <- phenos[setdiff(seq_len(n), result)]
    stop(
      "Circular dependency detected among derived_formula phenotypes: ",
      paste(cycle_phenos, collapse = ", "),
      ". Check your `formula` definitions for cycles.",
      call. = FALSE
    )
  }

  phenos[result]
}

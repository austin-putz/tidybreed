# ============================================================================
# formula_helpers.R
#
# Internal helpers for formula-based phenotype specification:
#   - formula_tbv : DSL string for composite TBV assembly
#   - formula     : arithmetic derivation over existing ind_phenotype records
# ============================================================================

# ── Constants ────────────────────────────────────────────────────────────────

.FORMULA_MATH_WHITELIST <- c(
  "sqrt", "log", "log2", "log10", "exp", "abs",
  "round", "ceiling", "floor", "sign", "trunc",
  "sin", "cos", "tan", "asin", "acos", "atan"
)

.FORMULA_TBV_DSL_FUNS <- c("self", "dam", "sire", "group_sum", "group_mean")

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


# ── formula_tbv validation ────────────────────────────────────────────────────

#' Validate a formula_tbv string at define_phenotype() time.
#'
#' Parses the DSL formula, walks the AST to collect trait references, and
#' checks that every trait name is present in trait_meta. Provides close-match
#' suggestions via agrep() for typos.
#'
#' @param formula_tbv Character. The DSL formula string.
#' @param known_traits Character vector of trait names from trait_meta.
#' @return Invisible NULL on success. Stops on error; warns for scalar constants.
#' @keywords internal
.validate_formula_tbv <- function(formula_tbv, known_traits) {
  expr <- tryCatch(
    parse(text = formula_tbv, keep.source = FALSE)[[1]],
    error = function(e)
      stop("Could not parse `formula_tbv = \"", formula_tbv, "\"`: ",
           conditionMessage(e), call. = FALSE)
  )

  walk_res <- .walk_formula_tbv_ast(expr)

  # Validate each referenced trait
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
      "Unknown trait name(s) in `formula_tbv = \"", formula_tbv, "\"`:\n",
      paste(suggestions, collapse = "\n"), "\n",
      "All symbols must be trait names defined with define_trait() in trait_meta.",
      call. = FALSE
    )
  }

  if (walk_res$has_scalar_constant)
    warning(
      "Scalar arithmetic constant detected in `formula_tbv = \"", formula_tbv, "\"`. ",
      "Adjusting TBVs by a fixed constant is unusual. Proceeding as requested.",
      call. = FALSE
    )

  invisible(NULL)
}


# ── derived formula validation ────────────────────────────────────────────────

#' Validate a derived formula string at define_phenotype() time.
#'
#' Symbols that are not in phenotype_meta generate a warning (not an error),
#' allowing config-first workflows where components are defined before their
#' dependents. A hard error at add_phenotype() time fires if still missing.
#'
#' @param formula Character. The arithmetic formula string.
#' @param known_phenos Character vector of phenotype names from phenotype_meta.
#' @return Invisible NULL on success.
#' @keywords internal
.validate_derived_formula <- function(formula, known_phenos) {
  expr <- tryCatch(
    parse(text = formula, keep.source = FALSE)[[1]],
    error = function(e)
      stop("Could not parse `formula = \"", formula, "\"`: ",
           conditionMessage(e), call. = FALSE)
  )

  symbols <- .extract_all_symbols(expr)
  symbols <- setdiff(symbols,
                     c(.FORMULA_MATH_WHITELIST, .FORMULA_ARITH_OPS))

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


# ── formula_tbv AST walk ──────────────────────────────────────────────────────

#' Walk the AST of a formula_tbv expression, returning structured trait refs.
#'
#' Each occurrence of a bare symbol or DSL function call gets a unique
#' placeholder name (e.g. ".tbv_self_WWD_1", ".tbv_dam_WWM_2") so that the
#' same trait can appear multiple times with distinct pre-fetched vectors.
#'
#' @param expr Parsed R expression (from parse()[[1]]).
#' @return A list:
#'   $trait_refs: list of lists, each with:
#'     - trait:       character trait name
#'     - type:        "self", "dam", "sire", "group_sum", or "group_mean"
#'     - col:         group column name (NA for non-group types)
#'     - table:       group table name  (NA for non-group types)
#'     - placeholder: unique R symbol name for the pre-fetched vector
#'   $has_scalar_constant: logical
#' @keywords internal
.walk_formula_tbv_ast <- function(expr) {
  trait_refs       <- list()
  has_scalar_const <- FALSE
  counter          <- 0L

  recurse <- function(e) {
    if (is.numeric(e) || is.complex(e)) {
      has_scalar_const <<- TRUE
      return()
    }
    if (is.name(e)) {
      nm <- as.character(e)
      if (nm %in% c(.FORMULA_ARITH_OPS, .FORMULA_MATH_WHITELIST,
                    .FORMULA_TBV_DSL_FUNS)) return()
      counter <<- counter + 1L
      ph <- paste0(".tbv_self_", nm, "_", counter)
      trait_refs[[length(trait_refs) + 1L]] <<- list(
        trait = nm, type = "self",
        col = NA_character_, table = NA_character_,
        placeholder = ph
      )
      return()
    }
    if (is.call(e)) {
      fn <- as.character(e[[1]])
      if (fn == "self") {
        tr <- as.character(e[[2]])
        counter <<- counter + 1L
        trait_refs[[length(trait_refs) + 1L]] <<- list(
          trait = tr, type = "self",
          col = NA_character_, table = NA_character_,
          placeholder = paste0(".tbv_self_", tr, "_", counter)
        )
        return()
      }
      if (fn %in% c("dam", "sire")) {
        tr <- as.character(e[[2]])
        counter <<- counter + 1L
        trait_refs[[length(trait_refs) + 1L]] <<- list(
          trait = tr, type = fn,
          col = NA_character_, table = NA_character_,
          placeholder = paste0(".tbv_", fn, "_", tr, "_", counter)
        )
        return()
      }
      if (fn %in% c("group_sum", "group_mean")) {
        tr  <- as.character(e[[2]])
        col <- as.character(e[[3]])
        # Optional table= argument (named or positional 4th)
        args <- as.list(e)[-1]
        tbl  <- "ind_meta"
        if (!is.null(names(args)) && "table" %in% names(args))
          tbl <- as.character(args[["table"]])
        else if (length(args) >= 3)
          tbl <- as.character(args[[3]])
        counter <<- counter + 1L
        trait_refs[[length(trait_refs) + 1L]] <<- list(
          trait = tr, type = fn,
          col = col, table = tbl,
          placeholder = paste0(".tbv_", fn, "_", tr, "_", col, "_", counter)
        )
        return()
      }
      # Arithmetic / other calls: recurse into arguments
      for (child in as.list(e)[-1]) recurse(child)
    }
  }

  recurse(expr)
  list(trait_refs = trait_refs, has_scalar_constant = has_scalar_const)
}


# ── formula_tbv AST substitution ─────────────────────────────────────────────

#' Replace DSL calls and bare symbols in the formula_tbv AST with placeholders.
#'
#' Traverses in the same depth-first order as .walk_formula_tbv_ast(), matching
#' each node to its pre-assigned placeholder from trait_refs (using sequential
#' .used flags so repeated occurrences are handled correctly).
#'
#' @param expr Parsed R expression.
#' @param trait_refs List from .walk_formula_tbv_ast()$trait_refs.
#' @return Modified expression with placeholder symbols.
#' @keywords internal
.substitute_tbv_ast <- function(expr, trait_refs) {
  used <- logical(length(trait_refs))

  transform <- function(e) {
    if (is.name(e)) {
      nm <- as.character(e)
      if (nm %in% c(.FORMULA_ARITH_OPS, .FORMULA_MATH_WHITELIST,
                    .FORMULA_TBV_DSL_FUNS)) return(e)
      for (i in seq_along(trait_refs)) {
        r <- trait_refs[[i]]
        if (!used[i] && r$type == "self" && r$trait == nm) {
          used[i] <<- TRUE
          return(as.name(r$placeholder))
        }
      }
      return(e)
    }
    if (is.call(e)) {
      fn <- as.character(e[[1]])
      if (fn %in% c("self", "dam", "sire")) {
        tr  <- as.character(e[[2]])
        typ <- if (fn == "self") "self" else fn
        for (i in seq_along(trait_refs)) {
          r <- trait_refs[[i]]
          if (!used[i] && r$type == typ && r$trait == tr) {
            used[i] <<- TRUE
            return(as.name(r$placeholder))
          }
        }
      }
      if (fn %in% c("group_sum", "group_mean")) {
        tr  <- as.character(e[[2]])
        col <- as.character(e[[3]])
        for (i in seq_along(trait_refs)) {
          r <- trait_refs[[i]]
          if (!used[i] && r$type == fn && r$trait == tr && r$col == col) {
            used[i] <<- TRUE
            return(as.name(r$placeholder))
          }
        }
      }
      # Recurse: arithmetic and other calls
      new_args <- lapply(as.list(e)[-1], transform)
      as.call(c(list(e[[1]]), new_args))
    } else {
      e  # numeric / logical literals pass through
    }
  }

  transform(expr)
}


# ── TBV vector pre-fetching ───────────────────────────────────────────────────

#' Fetch TBV values for a set of contributor IDs, aligned to focal_ids.
#' Returns NA where contributor is absent or has no TBV record.
#' @keywords internal
.fetch_contributor_tbvs <- function(pop, trait_safe, contr_ids, focal_ids) {
  n       <- length(focal_ids)
  missing <- is.na(contr_ids) | contr_ids == "NA" | !nzchar(as.character(contr_ids))
  notna   <- unique(as.character(contr_ids[!missing]))
  vec     <- rep(NA_real_, n)
  if (length(notna) > 0) {
    ids_sql  <- paste0("'", notna, "'", collapse = ", ")
    tbv_rows <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT id_ind, tbv_value FROM ind_tbv ",
      "WHERE trait_name = '", trait_safe, "' AND id_ind IN (", ids_sql, ")"
    ))
    if (nrow(tbv_rows) > 0) {
      tbv_map <- stats::setNames(tbv_rows$tbv_value, tbv_rows$id_ind)
      for (j in seq_len(n)) {
        if (!missing[j]) {
          val <- tbv_map[as.character(contr_ids[j])]
          vec[j] <- if (length(val) == 1 && !is.null(val)) val else NA_real_
        }
      }
    }
  }
  vec
}


#' Fetch and aggregate group-mate TBVs for a group_sum / group_mean term.
#'
#' Validates that grp_table and grp_col exist (errors with clear message).
#' Singleton pens (no group-mates) receive 0. Individuals with NA group get NA.
#'
#' @keywords internal
.fetch_group_tbvs <- function(pop, trait_safe, focal_ids, grp_col, grp_table,
                               agg_fn, phenotype_name) {
  n <- length(focal_ids)
  vec <- rep(NA_real_, n)

  available_tables <- DBI::dbListTables(pop$db_conn)
  if (!grp_table %in% available_tables)
    stop(
      "Group table '", grp_table, "' not found in database. ",
      "In formula_tbv for phenotype '", phenotype_name, "', ",
      "the `table=` argument of group_sum()/group_mean() must reference ",
      "an existing table. Available tables: ",
      paste(available_tables, collapse = ", "),
      call. = FALSE
    )

  available_cols <- DBI::dbListFields(pop$db_conn, grp_table)
  if (!grp_col %in% available_cols)
    stop(
      "Group column '", grp_col, "' not found in table '", grp_table, "'. ",
      "In formula_tbv for phenotype '", phenotype_name, "', ",
      "check the column name in group_sum()/group_mean(). ",
      "Use `table = \"<table_name>\"` if the column is in a different table. ",
      "Available columns in '", grp_table, "': ",
      paste(available_cols, collapse = ", "),
      call. = FALSE
    )

  focal_sql    <- paste0("'", focal_ids, "'", collapse = ", ")
  focal_grp_df <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, \"", grp_col, "\" AS group_val FROM ", grp_table,
    " WHERE id_ind IN (", focal_sql, ")"
  ))
  focal_grp_map <- stats::setNames(
    as.character(focal_grp_df$group_val), focal_grp_df$id_ind
  )

  non_na_grps <- unique(focal_grp_map[
    !is.na(focal_grp_map) & focal_grp_map != "NA" & nzchar(focal_grp_map)
  ])
  if (length(non_na_grps) == 0) return(vec)

  gv_sql       <- paste0("'", non_na_grps, "'", collapse = ", ")
  all_members  <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, \"", grp_col, "\" AS group_val FROM ", grp_table,
    " WHERE \"", grp_col, "\" IN (", gv_sql, ")"
  ))
  all_members$group_val <- as.character(all_members$group_val)

  all_member_ids <- unique(all_members$id_ind)
  tbv_map        <- stats::setNames(numeric(0), character(0))
  if (length(all_member_ids) > 0) {
    mem_sql  <- paste0("'", all_member_ids, "'", collapse = ", ")
    tbv_rows <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT id_ind, tbv_value FROM ind_tbv WHERE trait_name = '",
      trait_safe, "' AND id_ind IN (", mem_sql, ")"
    ))
    if (nrow(tbv_rows) > 0)
      tbv_map <- stats::setNames(tbv_rows$tbv_value, tbv_rows$id_ind)
  }

  for (j in seq_len(n)) {
    fid     <- focal_ids[j]
    grp_val <- focal_grp_map[fid]
    if (is.na(grp_val) || grp_val == "NA" || !nzchar(grp_val)) {
      vec[j] <- NA_real_
      next
    }
    mate_ids   <- all_members$id_ind[
      all_members$group_val == grp_val & all_members$id_ind != fid
    ]
    valid_tbvs <- tbv_map[mate_ids]
    valid_tbvs <- valid_tbvs[!is.na(valid_tbvs)]
    vec[j]     <- if (length(valid_tbvs) == 0) 0 else agg_fn(valid_tbvs)
  }

  vec
}


#' Pre-fetch all TBV vectors needed by a formula_tbv expression.
#'
#' Calls .fetch_contributor_tbvs() or .fetch_group_tbvs() for each trait_ref,
#' returning a named list ready to be used as an eval() environment.
#'
#' @param pop A tidybreed_pop object.
#' @param trait_refs List from .walk_formula_tbv_ast()$trait_refs.
#' @param subset_df Data frame: sex-filtered ind_meta rows for focal individuals.
#' @param phenotype_name Character. Used in error messages.
#' @return Named list: placeholder → numeric vector (length = nrow(subset_df)).
#' @keywords internal
.build_tbv_env <- function(pop, trait_refs, subset_df, phenotype_name) {
  focal_ids <- as.character(subset_df$id_ind)
  n         <- length(focal_ids)
  env_list  <- list()

  for (ref in trait_refs) {
    ph         <- ref$placeholder
    trait_safe <- gsub("'", "''", ref$trait)

    if (ref$type == "self") {
      ids_sql  <- paste0("'", focal_ids, "'", collapse = ", ")
      tbv_rows <- DBI::dbGetQuery(pop$db_conn, paste0(
        "SELECT id_ind, tbv_value FROM ind_tbv ",
        "WHERE trait_name = '", trait_safe, "' AND id_ind IN (", ids_sql, ")"
      ))
      vec <- stats::setNames(rep(NA_real_, n), focal_ids)
      if (nrow(tbv_rows) > 0) {
        tbv_map <- stats::setNames(tbv_rows$tbv_value, tbv_rows$id_ind)
        for (fid in focal_ids) {
          val <- tbv_map[fid]
          if (length(val) == 1 && !is.null(val) && !is.na(val)) vec[fid] <- val
        }
      }
      env_list[[ph]] <- unname(vec)

    } else if (ref$type == "dam") {
      contr_ids      <- as.character(subset_df$id_parent_2)
      env_list[[ph]] <- .fetch_contributor_tbvs(pop, trait_safe, contr_ids, focal_ids)

    } else if (ref$type == "sire") {
      contr_ids      <- as.character(subset_df$id_parent_1)
      env_list[[ph]] <- .fetch_contributor_tbvs(pop, trait_safe, contr_ids, focal_ids)

    } else if (ref$type %in% c("group_sum", "group_mean")) {
      agg_fn         <- if (ref$type == "group_sum") sum else mean
      env_list[[ph]] <- .fetch_group_tbvs(
        pop            = pop,
        trait_safe     = trait_safe,
        focal_ids      = focal_ids,
        grp_col        = ref$col,
        grp_table      = ref$table,
        agg_fn         = agg_fn,
        phenotype_name = phenotype_name
      )
    }
  }

  env_list
}


# ── Top-level formula_tbv evaluator ──────────────────────────────────────────

#' Evaluate a formula_tbv string for a set of individuals.
#'
#' Orchestrates: parse → AST walk → TBV pre-fetch → AST substitution → eval().
#'
#' @param pop A tidybreed_pop object.
#' @param formula_tbv Character. DSL formula string from phenotype_meta.
#' @param subset_df Data frame: sex-filtered ind_meta rows.
#' @param phenotype_name Character. Used in error messages.
#' @return Named numeric vector (names = id_ind). NA marks excluded individuals
#'         (missing dam/sire TBV, or NA group membership).
#' @keywords internal
.eval_formula_tbv <- function(pop, formula_tbv, subset_df, phenotype_name) {
  expr       <- parse(text = formula_tbv, keep.source = FALSE)[[1]]
  walk_res   <- .walk_formula_tbv_ast(expr)
  trait_refs <- walk_res$trait_refs

  tbv_env       <- .build_tbv_env(pop, trait_refs, subset_df, phenotype_name)
  modified_expr <- .substitute_tbv_ast(expr, trait_refs)

  result <- eval(modified_expr, envir = list2env(tbv_env, parent = baseenv()))
  stats::setNames(as.numeric(result), as.character(subset_df$id_ind))
}


# ── Derived formula pivot helper ─────────────────────────────────────────────

#' Pivot ind_phenotype rows wide for formula evaluation.
#'
#' Uses pheno_number = 1 by default; falls back to min pheno_number if no
#' pheno_number = 1 rows exist.
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

#' Evaluate a derived formula over existing ind_phenotype records.
#'
#' @param pop A tidybreed_pop object.
#' @param formula Character. Formula string from phenotype_meta.
#' @param ids Character vector. Individual IDs to compute for.
#' @param phenotype_name Character. Name of the derived phenotype (for messages).
#' @return Numeric vector, same length as ids.
#'         NA propagates naturally; Inf/NaN converted to NA with a warning.
#' @keywords internal
.eval_derived_formula <- function(pop, formula, ids, phenotype_name) {
  expr    <- parse(text = formula, keep.source = FALSE)[[1]]
  symbols <- .extract_all_symbols(expr)
  symbols <- setdiff(symbols, c(.FORMULA_MATH_WHITELIST, .FORMULA_ARITH_OPS))

  if (length(ids) == 0) return(numeric(0))

  ids_sql    <- paste0("'", ids, "'", collapse = ", ")
  phenos_sql <- paste0("'", gsub("'", "''", symbols), "'", collapse = ", ")
  rows <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT id_ind, phenotype_name, pheno_value, pheno_number FROM ind_phenotype ",
    "WHERE id_ind IN (", ids_sql, ") ",
    "AND phenotype_name IN (", phenos_sql, ")"
  ))

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
    eval(expr, envir = wide_df),
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
#' ind_phenotype first. formula_tbv and components phenotypes have no
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
    row <- pheno_meta[i, ]
    formula_col <- if ("formula" %in% names(row)) row$formula else NA_character_
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

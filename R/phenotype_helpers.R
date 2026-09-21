#' Encode a named numeric vector as a small JSON string
#'
#' @param x Named numeric vector.
#' @return A single character string, or `NA_character_` if `x` is empty.
#' @keywords internal
encode_levels_json <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  stopifnot(!is.null(names(x)))
  keys <- names(x)
  bad <- grepl("[\",:]", keys)
  if (any(bad)) {
    stop("Level names must not contain \" , or : (offending: ",
         paste(keys[bad], collapse = ", "), ")", call. = FALSE)
  }
  vals <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
  paste0("{", paste0("\"", keys, "\":", vals, collapse = ","), "}")
}


#' Decode a level-map JSON string produced by `encode_levels_json()`
#'
#' @param s Character scalar, or `NA`.
#' @return Named numeric vector (possibly length 0), or `NULL` if missing/empty.
#' @keywords internal
decode_levels_json <- function(s) {
  if (is.na(s) || !nzchar(s)) return(NULL)
  inner <- sub("^\\{", "", sub("\\}$", "", s))
  if (!nzchar(inner)) return(stats::setNames(numeric(0), character(0)))
  pairs <- strsplit(inner, ",", fixed = TRUE)[[1]]
  split_on_colon <- strsplit(pairs, ":", fixed = TRUE)
  keys <- vapply(split_on_colon, function(p) gsub("^\"|\"$", "", p[1]), character(1))
  vals <- vapply(split_on_colon, function(p) as.numeric(p[2]), numeric(1))
  stats::setNames(vals, keys)
}


#' Convert a liability vector to ordered integer categories
#'
#' @param liability Numeric vector.
#' @param thresholds Numeric vector of cutpoints (ascending).
#' @return Integer vector of category indices (1-based).
#' @keywords internal
liability_to_categorical <- function(liability, thresholds) {
  thresholds <- sort(thresholds)
  findInterval(liability, thresholds) + 1L
}


#' Clip a count vector to `[min_value, max_value]`
#'
#' @param x Numeric vector.
#' @param min_value,max_value Numeric scalars or `NA`.
#' @return Integer vector (rounded) clipped to bounds.
#' @keywords internal
clip_count <- function(x, min_value = NA_real_, max_value = NA_real_) {
  y <- round(x)
  if (!is.na(min_value)) y <- pmax(y, min_value)
  if (!is.na(max_value)) y <- pmin(y, max_value)
  as.integer(y)
}


#' Stage-1 covariate terms of one phenotype: fixed contributions and the
#' random-effect level every record touches
#'
#' No RNG, no writes. Fixed-class and fixed-covariate effects are evaluated to
#' a per-record contribution, with `null_class_action = "skip"` marking the
#' record `NA`. Random effects are **not** drawn here: the function records
#' which level of the grouping column each record falls in so that Stage 2
#' can resolve the draws once the record list is final (see
#' [add_phenotype_stages]).
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_name Character. Phenotype name to look up in `phenotype_effects`.
#' @param subset_df Data frame: the per-phenotype subset of `ind_meta` (already
#'   sex-filtered). Must contain `id_ind` and any `ind_meta` columns referenced
#'   by effects.
#' @return A list:
#'   - `fixed`: numeric vector of length `nrow(subset_df)`; `NA` for
#'     individuals excluded by `null_class_action = "skip"`.
#'   - `random`: a list, one element per random effect in `effect_name`
#'     order, each `list(effect_name, distribution, level)` where `level` is
#'     a character vector of length `nrow(subset_df)` (`NA` = no level).
#' @keywords internal
.ap_covariate_terms <- function(pop, phenotype_name, subset_df) {
  conn    <- pop$db_conn
  effects <- DBI::dbGetQuery(conn, paste0(
    "SELECT * FROM phenotype_effects WHERE phenotype_name = ",
    DBI::dbQuoteLiteral(conn, phenotype_name), " ORDER BY effect_name"))
  n_ind <- nrow(subset_df)
  if (nrow(effects) == 0) return(list(fixed = rep(0, n_ind), random = list()))

  ids_t     <- subset_df$id_ind
  total     <- rep(0, n_ind)
  skip_mask <- rep(FALSE, n_ind)
  random    <- list()

  for (i in seq_len(nrow(effects))) {
    e <- effects[i, ]

    src_tbl <- if (is.na(e$source_table) || !nzchar(e$source_table)) {
      "ind_meta"
    } else {
      e$source_table
    }

    if (!src_tbl %in% DBI::dbListTables(conn)) {
      stop("Effect '", e$effect_name, "': source table '", src_tbl,
           "' does not exist.", call. = FALSE)
    }
    if (!e$source_column %in% DBI::dbListFields(conn, src_tbl)) {
      stop("Effect '", e$effect_name, "': column '", e$source_column,
           "' not found in table '", src_tbl, "'.", call. = FALSE)
    }
    if (src_tbl == "ind_meta") {
      src_df <- subset_df
    } else {
      src_df <- .ap_read_by_id(conn, src_tbl, ids_t, e$source_column)
      if (any(duplicated(src_df$id_ind))) {
        stop("Effect '", e$effect_name, "': source table '", src_tbl,
             "' has duplicate id_ind rows for the current subset.", call. = FALSE)
      }
      src_df <- src_df[match(ids_t, src_df$id_ind), , drop = FALSE]
    }

    group <- src_df[[e$source_column]]
    ec    <- e$effect_class

    if (ec %in% c("fixed", "fixed_class")) {
      nca <- if (is.na(e$null_class_action)) "skip" else e$null_class_action

      null_mask <- is.na(group) | (as.character(group) == "NA") |
                   (as.character(group) == "")
      if (any(null_mask)) {
        if (nca == "error") {
          stop("Effect '", e$effect_name, "': ", sum(null_mask),
               " individual(s) have NULL in '", e$source_column, "'. ",
               "Set null_class_action = 'skip' or 'zero', or fix the data.",
               call. = FALSE)
        } else if (nca == "skip") {
          skip_mask <- skip_mask | null_mask
          group[null_mask] <- NA_character_
        } else {
          # "zero": treat NULL as reference level (zero contribution)
          group[null_mask] <- NA_character_
        }
      }

      level_map <- decode_levels_json(e$levels_json)
      group_chr <- as.character(group)
      shifts    <- unname(level_map[group_chr])

      # For null_class_action = "zero", NAs get 0
      if (nca == "zero") shifts[is.na(shifts) & null_mask] <- 0

      missing_lvls <- unique(group_chr[!is.na(group_chr) & is.na(shifts)])
      if (length(missing_lvls) > 0) {
        stop("Effect '", e$effect_name, "': levels have no shift defined: ",
             paste(missing_lvls, collapse = ", "),
             ". Update define_effect_fixed_class() to include all levels.",
             call. = FALSE)
      }
      shifts[is.na(shifts)] <- 0
      total <- total + shifts

    } else if (ec == "fixed_cov") {
      vals       <- as.numeric(src_df[[e$source_column]])
      center_val <- if (is.na(e$center)) 0 else e$center
      slope_val  <- if (is.na(e$slope))  0 else e$slope
      p_ord      <- if (is.na(e$poly_order) || e$poly_order < 1L) 1L
                    else as.integer(e$poly_order)
      total <- total + slope_val * (vals - center_val)^p_ord

    } else {
      # Random effect: record the level per individual; Stage 2 draws
      level <- as.character(group)
      level[is.na(group)] <- NA_character_
      random[[length(random) + 1L]] <- list(
        effect_name  = e$effect_name,
        distribution = e$distribution %||% "normal",
        level        = level)
    }
  }

  total[skip_mask] <- NA_real_
  list(fixed = total, random = random)
}


#' Null-coalescing operator for internal use
#' @noRd
`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x


#' Compute the next pheno_number for each individual for a given phenotype
#'
#' RNG-neutral: the ids are registered as a temporary view and joined, never
#' written to a table or pasted into SQL.
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_name Character. Phenotype name.
#' @param ids Character vector of individual IDs.
#' @return Integer vector, same length and order as `ids`.
#' @keywords internal
next_pheno_numbers <- function(pop, phenotype_name, ids) {
  if (length(ids) == 0) return(integer(0))
  conn <- pop$db_conn
  tmp  <- "__ap_pn_ids"
  duckdb::duckdb_register(conn, tmp, data.frame(id_ind = unique(ids),
                                                stringsAsFactors = FALSE))
  on.exit(try(duckdb::duckdb_unregister(conn, tmp), silent = TRUE), add = TRUE)
  res <- DBI::dbGetQuery(conn, paste0(
    "SELECT f.id_ind, COALESCE(MAX(p.pheno_number), 0) + 1 AS next_pn ",
    "FROM ", tmp, " AS f ",
    "LEFT JOIN ind_phenotype AS p ",
    "  ON p.id_ind = f.id_ind AND p.phenotype_name = ",
    DBI::dbQuoteLiteral(conn, phenotype_name), " ",
    "GROUP BY f.id_ind"))
  as.integer(stats::setNames(res$next_pn, res$id_ind)[ids])
}


#' Generate next global integer IDs for ind_phenotype
#'
#' @param pop A `tidybreed_pop` object.
#' @param n Number of IDs to generate.
#' @return Integer vector of length `n`.
#' @keywords internal
next_phenotype_ids <- function(pop, n) {
  if (n == 0L) return(integer(0))
  start <- next_int_id(pop$db_conn, "ind_phenotype", "id_phenotype")
  seq.int(start, start + n - 1L)
}

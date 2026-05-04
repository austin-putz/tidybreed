# Hard-coded descriptions for all core tables
.TABLE_DESCRIPTIONS <- c(
  genome_meta          = "Locus metadata (chr, position, chip flags, QTL flags)",
  genome_haplotype     = "Phased haplotypes — 2 rows per individual",
  genome_genotype      = "Genotypes in 0/1/2 encoding — 1 row per individual",
  founder_haplotypes   = "Founder haplotype pool used by add_founders()",
  ind_meta             = "Individual metadata (IDs, parents, line, sex, user cols)",
  ind_phenotype        = "Phenotype records — long format",
  ind_tbv              = "True breeding values",
  ind_ebv              = "Estimated breeding values (BLUP/GBLUP)",
  trait_meta           = "Trait definitions",
  trait_effects        = "Fixed and random effect configurations per trait",
  trait_effect_cov     = "Variance/covariance matrices (gen_add, residual, random)",
  trait_random_effects = "Sampled random effect levels",
  index_meta           = "Selection index definitions",
  ind_index            = "Computed selection index values"
)


#' Summarize a tidybreed population
#'
#' @description
#' Produces a detailed summary of all tables in a `tidybreed_pop` object.
#' For each table, shows row count, column count, and per-column statistics:
#' frequency tables for low-cardinality categorical columns, 5-number
#' summaries for numeric columns, and date ranges for date/timestamp columns.
#' Wide genome tables (`genome_haplotype`, `genome_genotype`) only summarize
#' their non-locus metadata columns to avoid querying 50,000-column tables.
#'
#' @param object A `tidybreed_pop` object.
#' @param tables Character vector of table names to include. `NULL` (default)
#'   includes all tables registered in `object$tables`.
#' @param max_values Integer. Columns with at most this many distinct values
#'   receive a frequency table display; default `12L`.
#' @param ... Additional arguments (not used; required by the S3 generic).
#'
#' @return A `tidybreed_summary` object (a named list).
#'
#' @seealso [print.tidybreed_pop()]
#'
#' @examples
#' \dontrun{
#' pop <- initialize_genome("MySim", n_loci = 1000, n_chr = 5, chr_len_Mb = 100,
#'                          n_haplotypes = 50)
#' pop <- add_founders(pop, n_males = 50, n_females = 50, line_name = "A", gen = 0L)
#' summary(pop)
#' summary(pop, tables = c("ind_meta", "genome_meta"))
#' summary(pop, max_values = 5L)
#' }
#' @export
summary.tidybreed_pop <- function(object, tables = NULL, max_values = 12L, ...) {
  if (!inherits(object, "tidybreed_pop")) {
    stop("object must be a tidybreed_pop", call. = FALSE)
  }

  all_tables <- object$tables

  if (is.null(tables)) {
    req_tables <- all_tables
  } else {
    unknown <- setdiff(tables, all_tables)
    if (length(unknown) > 0) {
      warning("Tables not found in population: ",
              paste(unknown, collapse = ", "), call. = FALSE)
    }
    req_tables <- intersect(tables, all_tables)
  }

  conn <- object$db_conn

  n_ind <- if ("ind_meta" %in% all_tables) {
    as.integer(DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM ind_meta")$n)
  } else NA_integer_

  n_loci <- if ("genome_meta" %in% all_tables) {
    as.integer(DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM genome_meta")$n)
  } else NA_integer_

  n_traits <- if ("trait_meta" %in% all_tables) {
    as.integer(DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM trait_meta")$n)
  } else NA_integer_

  table_summaries <- lapply(req_tables, function(tbl_name) {
    .tb_summarize_table(conn, tbl_name, as.integer(max_values))
  })
  names(table_summaries) <- req_tables

  structure(
    list(
      pop_name = object$pop_name,
      db_path  = object$db_path,
      n_tables = length(all_tables),
      n_ind    = n_ind,
      n_loci   = n_loci,
      n_traits = n_traits,
      tables   = table_summaries,
      call     = sys.call(0)
    ),
    class = c("tidybreed_summary", "list")
  )
}


#' Print a tidybreed_summary object
#'
#' @param x A `tidybreed_summary` object from [summary.tidybreed_pop()].
#' @param ... Additional arguments (not used).
#' @return `x`, invisibly.
#' @export
print.tidybreed_summary <- function(x, ...) {
  width <- getOption("width", 80L)

  # Header
  .tb_ruler(sprintf("tidybreed Population: %s", x$pop_name), width)
  cat("  Database :", x$db_path, "\n")

  stats <- paste0("  Tables   : ", x$n_tables)
  if (!is.na(x$n_ind))    stats <- paste0(stats, "  │  Individuals: ", format(x$n_ind,  big.mark = ","))
  if (!is.na(x$n_loci))   stats <- paste0(stats, "  │  Loci: ",         format(x$n_loci, big.mark = ","))
  if (!is.na(x$n_traits)) stats <- paste0(stats, "  │  Traits: ",       x$n_traits)
  cat(stats, "\n")

  for (tbl_name in names(x$tables)) {
    cat("\n")
    tbl <- x$tables[[tbl_name]]

    left  <- paste0("── ", tbl_name, " ")
    right <- paste0(" ", format(tbl$n_rows, big.mark = ","),
                    " rows · ", format(tbl$n_cols, big.mark = ","), " cols")
    fill_w <- max(2L, width - nchar(left) - nchar(right))
    cat(left, strrep("─", fill_w), right, "\n", sep = "")

    if (tbl$n_rows == 0L) {
      cat("  (empty)\n")
      next
    }

    cs <- tbl$col_summaries
    if (nrow(cs) > 0) {
      col_w  <- min(30L, max(nchar(cs$col_name)))
      type_w <- max(nchar(cs$col_type))
      uniq_w <- max(nchar(as.character(cs$n_unique)))

      for (i in seq_len(nrow(cs))) {
        cat("  ",
            formatC(cs$col_name[i], width = -col_w),  "  ",
            formatC(cs$col_type[i], width = -type_w), "  ",
            formatC(cs$n_unique[i], width =  uniq_w), "    ",
            cs$display[i], "\n", sep = "")
      }
    }

    if (tbl$is_wide_genome && tbl$n_locus_cols > 0) {
      cat("  [", format(tbl$n_locus_cols, big.mark = ","),
          " locus columns not summarized]\n", sep = "")
    }
  }

  invisible(x)
}


# ── Internal helpers ──────────────────────────────────────────────────────────

.tb_ruler <- function(title, width) {
  left   <- paste0("── ", title, " ")
  fill_w <- max(2L, width - nchar(left))
  cat(left, strrep("─", fill_w), "\n", sep = "")
}


.tb_summarize_table <- function(conn, tbl_name, max_values) {
  desc <- unname(.TABLE_DESCRIPTIONS[tbl_name])
  if (is.na(desc)) desc <- ""

  all_cols   <- DBI::dbListFields(conn, tbl_name)
  locus_cols <- grep("^locus_[0-9]+$", all_cols, value = TRUE)
  is_wide    <- length(locus_cols) > 0
  meta_cols  <- if (is_wide) setdiff(all_cols, locus_cols) else all_cols
  n_cols     <- length(all_cols)

  n_rows <- as.integer(
    DBI::dbGetQuery(conn, paste0("SELECT COUNT(*) AS n FROM ", tbl_name))$n
  )

  empty_cs <- data.frame(
    col_name    = character(),
    col_type    = character(),
    n_unique    = integer(),
    n_nulls_pct = numeric(),
    display     = character(),
    stringsAsFactors = FALSE
  )

  if (n_rows == 0L || length(meta_cols) == 0L) {
    return(list(
      name          = tbl_name,
      description   = desc,
      n_rows        = 0L,
      n_cols        = n_cols,
      is_wide_genome = is_wide,
      n_locus_cols  = length(locus_cols),
      col_summaries = empty_cs
    ))
  }

  col_list <- paste(meta_cols, collapse = ", ")
  summ <- tryCatch(
    DBI::dbGetQuery(conn, paste0("SUMMARIZE SELECT ", col_list, " FROM ", tbl_name)),
    error = function(e) NULL
  )

  if (is.null(summ) || nrow(summ) == 0L) {
    return(list(
      name          = tbl_name,
      description   = desc,
      n_rows        = n_rows,
      n_cols        = n_cols,
      is_wide_genome = is_wide,
      n_locus_cols  = length(locus_cols),
      col_summaries = empty_cs
    ))
  }

  col_summaries <- do.call(rbind, lapply(seq_len(nrow(summ)), function(i) {
    row <- as.list(summ[i, ])
    data.frame(
      col_name    = as.character(row$column_name),
      col_type    = as.character(row$column_type),
      n_unique    = as.integer(row$approx_unique),
      n_nulls_pct = suppressWarnings(as.numeric(row$null_percentage)),
      display     = .tb_format_col(row, conn, tbl_name, max_values, n_rows),
      stringsAsFactors = FALSE
    )
  }))

  list(
    name          = tbl_name,
    description   = desc,
    n_rows        = n_rows,
    n_cols        = n_cols,
    is_wide_genome = is_wide,
    n_locus_cols  = length(locus_cols),
    col_summaries = col_summaries
  )
}


.tb_format_col <- function(row, conn, tbl_name, max_values, n_rows) {
  col_name <- as.character(row$column_name)
  col_type <- as.character(row$column_type)
  n_unique <- as.integer(row$approx_unique)
  null_pct <- suppressWarnings(as.numeric(row$null_percentage))

  n_nas  <- if (!is.na(null_pct) && null_pct > 0) as.integer(round(null_pct / 100 * n_rows)) else 0L
  na_sfx <- if (n_nas > 0L) sprintf("  (%d NAs)", n_nas) else ""

  is_id   <- grepl("^id_|_id$", col_name)
  is_dbl  <- col_type %in% c("DOUBLE", "FLOAT", "REAL", "DECIMAL", "NUMERIC",
                              "FLOAT4", "FLOAT8")
  is_int  <- col_type %in% c("INTEGER", "BIGINT", "SMALLINT", "TINYINT",
                              "UBIGINT", "UINTEGER", "USMALLINT", "UTINYINT",
                              "HUGEINT", "INT4", "INT8", "INT2", "INT1", "INT")
  is_date <- col_type == "DATE"
  is_ts   <- startsWith(col_type, "TIMESTAMP")

  safe_str <- function(x) {
    v <- tryCatch(as.character(x), error = function(e) NA_character_)
    if (length(v) == 0L || is.na(v) || nchar(v) == 0L) NA_character_ else v
  }

  fmt_num <- function(x) {
    v <- suppressWarnings(as.numeric(x))
    if (is.null(v) || length(v) == 0L || is.na(v)) return("NA")
    formatC(v, format = "g", digits = 6)
  }

  if (is_id && !is_dbl) {
    if (n_unique == 0L) {
      return(paste0(sprintf("%d unique", n_unique), na_sfx))
    }
    min_val <- safe_str(row$min)
    max_val <- safe_str(row$max)
    min_str <- if (is.na(min_val)) "?" else min_val
    max_str <- if (is.na(max_val)) "?" else max_val
    sprintf('%d unique  "%s" … "%s"%s', n_unique, min_str, max_str, na_sfx)

  } else if (is_dbl || (is_int && n_unique > max_values)) {
    sprintf("Min: %s  Q1: %s  Med: %s  Q3: %s  Max: %s%s",
            fmt_num(row$min), fmt_num(row$q25), fmt_num(row$q50),
            fmt_num(row$q75), fmt_num(row$max), na_sfx)

  } else if (is_date || is_ts) {
    min_val <- safe_str(row$min)
    max_val <- safe_str(row$max)
    min_str <- if (is.na(min_val)) "?" else min_val
    max_str <- if (is.na(max_val)) "?" else max_val
    sprintf("%s → %s%s", min_str, max_str, na_sfx)

  } else if (n_unique <= max_values) {
    freq_str <- .tb_freq_table(conn, tbl_name, col_name, max_values)
    if (is.null(freq_str)) {
      paste0(sprintf("%d unique", n_unique), na_sfx)
    } else {
      paste0(freq_str, na_sfx)
    }

  } else {
    # High-cardinality VARCHAR or other types
    ex_sql <- paste0(
      "SELECT CAST(", col_name, " AS VARCHAR) AS v FROM ", tbl_name,
      " WHERE ", col_name, " IS NOT NULL LIMIT 2"
    )
    ex <- tryCatch(DBI::dbGetQuery(conn, ex_sql)$v, error = function(e) character(0))
    if (length(ex) >= 2L) {
      sprintf('%d unique  (e.g. "%s", "%s")%s', n_unique, ex[1], ex[2], na_sfx)
    } else if (length(ex) == 1L) {
      sprintf('%d unique  (e.g. "%s")%s', n_unique, ex[1], na_sfx)
    } else {
      sprintf("%d unique%s", n_unique, na_sfx)
    }
  }
}


.tb_freq_table <- function(conn, tbl_name, col_name, max_values) {
  sql <- paste0(
    "SELECT CAST(", col_name, " AS VARCHAR) AS val, COUNT(*) AS n ",
    "FROM ", tbl_name,
    " WHERE ", col_name, " IS NOT NULL ",
    "GROUP BY val ORDER BY n DESC, val ",
    "LIMIT ", max_values + 1L
  )
  freq <- tryCatch(DBI::dbGetQuery(conn, sql), error = function(e) NULL)
  if (is.null(freq) || nrow(freq) == 0L) return(NULL)
  if (nrow(freq) > max_values) return(NULL)
  paste(paste0(freq$val, ": ", freq$n), collapse = ", ")
}

#' Shared SQL-identifier validation utilities
#'
#' Internal helpers used by functions that accept a user-supplied identifier
#' (trait name, effect name, column name, model label) and then interpolate
#' it into SQL. Keeping the rules in one place guarantees every caller
#' enforces the same whitelist and rejects the same reserved keywords.
#'
#' @keywords internal
#' @name sql_utils
NULL


#' DuckDB reserved keywords that must not be used as user identifiers.
#'
#' @keywords internal
SQL_KEYWORDS <- c(
  "SELECT", "FROM", "WHERE", "INSERT", "UPDATE", "DELETE",
  "CREATE", "DROP", "ALTER", "TABLE", "INDEX", "PRIMARY",
  "KEY", "FOREIGN", "REFERENCES", "AND", "OR", "NOT", "NULL",
  "TRUE", "FALSE", "AS", "ON", "JOIN", "LEFT", "RIGHT", "INNER"
)


#' Regex that matches a strictly safe SQL identifier
#' (letter, then letters/digits/underscores).
#'
#' @keywords internal
SQL_IDENTIFIER_RE <- "^[a-zA-Z][a-zA-Z0-9_]*$"


#' Assert that `name` is a valid, non-reserved SQL identifier.
#'
#' Used by every function that interpolates user input into SQL to protect
#' against SQL injection.
#'
#' @param name Character scalar to validate.
#' @param what Human-readable role label used in error messages
#'   (e.g. `"trait name"`, `"effect name"`).
#' @param reserved Optional character vector of additional reserved names
#'   for the calling context (e.g. reserved column names in a table).
#' @return Invisible `NULL` on success; errors otherwise.
#' @keywords internal
validate_sql_identifier <- function(name, what = "identifier",
                                    reserved = character()) {

  if (!is.character(name) || length(name) != 1) {
    stop(what, " must be a single character string", call. = FALSE)
  }
  if (length(reserved) > 0 && name %in% reserved) {
    stop("Cannot modify reserved column '", name, "'. ",
         "Reserved columns: ", paste(reserved, collapse = ", "),
         call. = FALSE)
  }
  if (!grepl(SQL_IDENTIFIER_RE, name)) {
    stop("Invalid ", what, " '", name, "'. ",
         "Must start with a letter and contain only letters, numbers, ",
         "and underscores.",
         call. = FALSE)
  }
  if (toupper(name) %in% SQL_KEYWORDS) {
    stop(what, " '", name, "' is a SQL reserved keyword.", call. = FALSE)
  }
  invisible(NULL)
}


#' Reserved columns per table (cannot be modified by the user)
#'
#' @keywords internal
TABLE_RESERVED_COLS <- list(
  ind_meta         = c("id_ind", "id_parent_1", "id_parent_2", "line_name", "sex"),
  genome_meta      = c("locus_id", "locus_name", "chr", "chr_name", "pos_Mb"),
  ind_phenotype    = c("id_phenotype", "id_ind", "trait_name", "pheno_value", "pheno_number",
                       "liability_value", "cat_name"),
  ind_tbv          = c("id_tbv", "id_ind", "trait_name", "tbv_value"),
  ind_ebv          = c("id_ebv", "id_ind", "trait_name", "model", "ebv_value", "acc", "se", "eval_number"),
  trait_meta       = c("id_trait", "trait_name", "description", "units", "trait_type", "repeatable",
                       "recorded_on", "expressed_sex", "expressed_parent", "target_add_mean",
                       "min_value", "max_value", "prevalence", "thresholds",
                       "index_weight", "economic_value",
                       "cat_values", "cat_names", "store_liability"),
  trait_effects    = c("trait_name", "effect_name", "effect_class", "source_column",
                       "source_table", "distribution", "levels_json", "slope",
                       "center", "value"),
  trait_effect_cov = c("id_trait_effect_cov", "effect_name", "trait_name_1", "trait_name_2", "cov_value"),
  index_meta       = c("id_index_name", "index_name", "trait_name", "index_weight"),
  ind_index        = c("id_index", "id_ind", "index_name", "index_number", "index_value")
)


#' Primary key column per table (used for vector updates and filtered updates)
#'
#' @keywords internal
TABLE_PRIMARY_KEYS <- list(
  ind_meta         = "id_ind",
  genome_meta      = "locus_id",
  ind_phenotype    = "id_phenotype",
  ind_tbv          = "id_tbv",
  ind_ebv          = "id_ebv",
  trait_meta       = "id_trait",
  trait_effect_cov = "id_trait_effect_cov",
  index_meta       = "id_index_name",
  ind_index        = "id_index"
)


#' Composite row keys per table (used for exact single-table DELETE)
#'
#' For tables with a single PK this is length-1. For composite-PK tables all
#' key columns are listed so that `delete_rows()` can target exact rows rather
#' than over-deleting by `id_ind` alone.
#'
#' @keywords internal
TABLE_ROW_KEYS <- list(
  ind_meta         = "id_ind",
  genome_meta      = "locus_id",
  ind_phenotype    = "id_phenotype",
  trait_meta       = "id_trait",
  ind_tbv          = c("id_ind", "trait_name"),
  ind_ebv          = c("id_ind", "trait_name", "model", "eval_number"),
  ind_index        = c("id_ind", "index_name", "index_number"),
  genome_haplotype = c("id_ind", "parent_origin"),
  genome_genotype  = "id_ind",
  trait_effects    = c("trait_name", "effect_name"),
  trait_effect_cov = "id_trait_effect_cov",
  index_meta       = "id_index_name"
)


#' Get the next integer ID for a table's primary-key sequence.
#'
#' Returns MAX(id_col) + 1, or 1 if the table is empty.
#' For multi-row inserts, pre-generate a range:
#'   `seq.int(next_int_id(conn, tbl, col), length.out = n)`
#'
#' @param conn A DBI connection.
#' @param table Character. Table name.
#' @param id_col Character. Integer PK column name.
#' @return Integer scalar.
#' @keywords internal
next_int_id <- function(conn, table, id_col) {
  as.integer(DBI::dbGetQuery(
    conn,
    paste0("SELECT COALESCE(MAX(", id_col, "), 0) + 1 AS nxt FROM ", table)
  )$nxt)
}


#' Tables that contain the `id_ind` column (used for cross-table deletion)
#'
#' @keywords internal
IND_TABLE_ID_IND_COLS <- c(
  "ind_meta", "ind_phenotype", "ind_tbv", "ind_ebv",
  "ind_index", "genome_haplotype", "genome_genotype"
)


#' Infer DuckDB column type from an R value
#'
#' @param value An R value (scalar or vector)
#' @return Character string representing the DuckDB type
#' @keywords internal
infer_duckdb_type <- function(value) {

  if (is.logical(value) && all(is.na(value))) {
    warning(
      "Cannot infer type from NA value. Defaulting to VARCHAR. ",
      "Use typed NA (e.g. NA_real_, NA_integer_, NA_character_) or supply a ",
      "non-NA seed value first to get the correct type.",
      call. = FALSE
    )
    return("VARCHAR")
  }

  if (is.logical(value))                              return("BOOLEAN")
  if (is.integer(value))                              return("INTEGER")
  if (inherits(value, "Date"))                        return("DATE")
  if (inherits(value, "POSIXct") || inherits(value, "POSIXlt")) return("TIMESTAMP")
  if (is.numeric(value))                              return("DOUBLE")
  if (is.character(value))                            return("VARCHAR")

  sample_val <- value[!is.na(value)][1]
  if (is.na(sample_val)) {
    warning(
      "Cannot infer type from NA value. Defaulting to VARCHAR. ",
      "Use typed NA (e.g. NA_real_, NA_integer_, NA_character_) or supply a ",
      "non-NA seed value first to get the correct type.",
      call. = FALSE
    )
    return("VARCHAR")
  }

  stop(
    "Unsupported type: ", class(value)[1], ". ",
    "Supported types: logical, integer, numeric, Date, POSIXct, character",
    call. = FALSE
  )
}


#' Validate, type-infer, and expand extra user columns for pre-insert attachment
#'
#' Called by `add_founders()`, `add_offspring()`, `add_phenotype()`,
#' `add_tbv()`, and `add_ebv()` to process the `...` (extra field) arguments
#' before they are attached to an insertion data frame.
#'
#' Scalars are broadcast to `n_rows`. Vectors must already have length `n_rows`.
#' Any column that does not yet exist in the target table is added via
#' `ALTER TABLE ADD COLUMN` before the caller writes the data frame.
#'
#' @param extra_cols Named list of values (captured from `list(...)`).
#' @param n_rows Integer. Number of rows in the insertion data frame.
#' @param table_name Character. Target table name.
#' @param conn DuckDB connection.
#' @return Named list of expanded value vectors, one per field.
#' @keywords internal
prepare_extra_cols <- function(extra_cols, n_rows, table_name, conn) {
  if (length(extra_cols) == 0) return(list())

  reserved      <- TABLE_RESERVED_COLS[[table_name]]
  if (is.null(reserved)) reserved <- character(0)
  existing_cols <- DBI::dbListFields(conn, table_name)

  result <- list()
  for (field_name in names(extra_cols)) {
    value <- extra_cols[[field_name]]

    validate_sql_identifier(field_name, what = "field name", reserved = reserved)

    db_type <- infer_duckdb_type(value)

    if (length(value) == 1L) {
      value <- rep(value, n_rows)
    } else if (length(value) != n_rows) {
      stop(
        "Vector length (", length(value), ") for field '", field_name,
        "' must equal the number of rows being inserted (", n_rows, ").",
        call. = FALSE
      )
    }

    if (!field_name %in% existing_cols) {
      DBI::dbExecute(conn, paste0(
        "ALTER TABLE ", table_name, " ADD COLUMN ", field_name, " ", db_type
      ))
      existing_cols <- c(existing_cols, field_name)
      message("Added new column '", field_name, "' (", db_type,
              ") to `", table_name, "`")
    }

    result[[field_name]] <- value
  }
  result
}


#' Format a scalar R value as a SQL literal
#'
#' @param value A scalar R value (length-1 vector)
#' @param db_type DuckDB type string from `infer_duckdb_type()`
#' @return A character string safe to embed in a SQL statement
#' @keywords internal
format_sql_value <- function(value, db_type) {
  if (is.na(value))                           return("NULL")
  # Check actual class before db_type — guards against cases where a Date/TIMESTAMP
  # value loses its class (e.g. for-loop iteration stripping class attributes) and
  # infer_duckdb_type() falls through to "DOUBLE", which would embed the date string
  # unquoted in SQL (e.g. 2026-02-26 instead of DATE '2026-02-26').
  if (inherits(value, "POSIXct") || inherits(value, "POSIXlt"))
    return(paste0("TIMESTAMP '", format(value, "%Y-%m-%d %H:%M:%S"), "'"))
  if (inherits(value, "Date"))
    return(paste0("DATE '", as.character(value), "'"))
  if (db_type == "BOOLEAN")                   return(ifelse(value, "TRUE", "FALSE"))
  if (db_type %in% c("INTEGER", "DOUBLE"))    return(as.character(value))
  if (db_type == "DATE")                      return(paste0("DATE '", as.character(value), "'"))
  if (db_type == "TIMESTAMP")                 return(paste0("TIMESTAMP '", format(value, "%Y-%m-%d %H:%M:%S"), "'"))
  if (db_type == "VARCHAR") {
    escaped <- gsub("'", "''", as.character(value))
    return(paste0("'", escaped, "'"))
  }
  stop("Unknown db_type: ", db_type, call. = FALSE)
}

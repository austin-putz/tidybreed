#' Shared SQL-identifier validation utilities
#'
#' Internal helpers used by functions that accept a user-supplied identifier
#' (trait name, effect name, column name, model label) and then interpolate
#' it into SQL. Keeping the rules in one place guarantees every caller
#' enforces the same whitelist and rejects the same reserved keywords.
#'
#' @keywords internal
#' @name sql_utils


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
  ind_meta    = c("id_ind", "id_parent_1", "id_parent_2", "line", "sex"),
  genome_meta = c("locus_id", "locus_name", "chr", "chr_name", "pos_Mb")
)


#' Primary key column per table (used for vector updates and filtered updates)
#'
#' @keywords internal
TABLE_PRIMARY_KEYS <- list(
  ind_meta    = "id_ind",
  genome_meta = "locus_id"
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
  if (is.numeric(value))                              return("DOUBLE")
  if (inherits(value, "Date"))                        return("DATE")
  if (inherits(value, "POSIXct") || inherits(value, "POSIXlt")) return("TIMESTAMP")
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


#' Format a scalar R value as a SQL literal
#'
#' @param value A scalar R value (length-1 vector)
#' @param db_type DuckDB type string from `infer_duckdb_type()`
#' @return A character string safe to embed in a SQL statement
#' @keywords internal
format_sql_value <- function(value, db_type) {
  if (is.na(value))                           return("NULL")
  if (db_type == "BOOLEAN")                   return(ifelse(value, "TRUE", "FALSE"))
  if (db_type %in% c("INTEGER", "DOUBLE"))    return(as.character(value))
  if (db_type == "DATE")                      return(paste0("'", as.character(value), "'"))
  if (db_type == "TIMESTAMP")                 return(paste0("'", format(value, "%Y-%m-%d %H:%M:%S"), "'"))
  if (db_type == "VARCHAR") {
    escaped <- gsub("'", "''", as.character(value))
    return(paste0("'", escaped, "'"))
  }
  stop("Unknown db_type: ", db_type, call. = FALSE)
}

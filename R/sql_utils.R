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

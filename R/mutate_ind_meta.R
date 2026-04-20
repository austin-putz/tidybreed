#' Add or modify columns in individual metadata table
#'
#' @description
#' Adds new columns to the `ind_meta` table or updates existing custom columns.
#' Values can be scalar (applied to all individuals) or vectors (one per individual).
#' The function automatically infers DuckDB column types from R data types.
#'
#' Reserved columns (`id_ind`, `id_parent_1`, `id_parent_2`, `line`, `sex`) cannot be modified.
#'
#' @param pop A `tidybreed_pop` object
#' @param ... Named arguments specifying column names and values.
#'   Each argument should be `field_name = value`, where value is either:
#'   - A scalar (applied to all individuals)
#'   - A vector of length equal to number of individuals
#'
#' @return The `tidybreed_pop` object (invisibly), modified in place
#'
#' @details
#' **Type Inference:**
#' - `logical` → BOOLEAN
#' - `integer` → INTEGER
#' - `numeric` → DOUBLE
#' - `Date` → DATE
#' - `POSIXct`/`POSIXlt` → TIMESTAMP
#' - `character` → VARCHAR
#'
#' **Field Naming Rules:**
#' - Must start with a letter
#' - Can contain letters, numbers, and underscores
#' - Cannot be SQL reserved words
#' - Cannot be reserved columns (id_ind, id_parent_1, id_parent_2, line, sex)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Initialize and add founders
#' pop <- initialize_genome(pop_name = "A", n_loci = 1000, n_chr = 10) %>%
#'   add_founders(n_males = 10, n_females = 100, line_name = "A")
#'
#' # Add scalar fields (same value for all individuals)
#' pop <- pop %>%
#'   mutate_ind_meta(
#'     gen = 0,
#'     farm = "A",
#'     date_birth = Sys.Date()
#'   )
#'
#' # Add vector field (different value per individual)
#' n_ind <- 110
#' pop <- pop %>%
#'   mutate_ind_meta(
#'     weight_kg = rnorm(n_ind, mean = 80, sd = 10),
#'     is_selected = sample(c(TRUE, FALSE), n_ind, replace = TRUE)
#'   )
#'
#' # Update existing field
#' pop <- pop %>%
#'   mutate_ind_meta(gen = 1)
#'
#' # View updated metadata
#' get_table(pop, "ind_meta") %>% collect()
#' }
mutate_ind_meta <- function(pop, ...) {

  # 1. Validate inputs
  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)

  # 2. Check ind_meta table exists
  if (!"ind_meta" %in% pop$tables) {
    stop(
      "ind_meta table does not exist. Call add_founders() first.",
      call. = FALSE
    )
  }

  # 3. Capture named arguments
  field_values <- list(...)

  if (length(field_values) == 0) {
    warning("No fields specified. Returning population unchanged.", call. = FALSE)
    return(invisible(pop))
  }

  # 4. Get current number of individuals
  n_ind <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT COUNT(*) as n FROM ind_meta"
  )$n

  if (n_ind == 0) {
    warning("ind_meta table is empty. No values to update.", call. = FALSE)
    return(invisible(pop))
  }

  # 5. Get existing columns
  existing_cols <- DBI::dbListFields(pop$db_conn, "ind_meta")

  # 6. Process each field
  for (field_name in names(field_values)) {

    value <- field_values[[field_name]]

    # 6a. Validate field name
    validate_field_name(field_name, existing_cols)

    # 6b. Infer DuckDB type (before length check so unsupported types fail early)
    db_type <- infer_duckdb_type(value)

    # 6c. Determine if value is scalar or vector
    if (length(value) == 1) {
      # Scalar - will be replicated for all individuals
      value_to_use <- value
      is_vector <- FALSE
    } else {
      # Vector - must match number of individuals
      if (length(value) != n_ind) {
        stop(
          "Length of '", field_name, "' (", length(value),
          ") does not match number of individuals (", n_ind, ")",
          call. = FALSE
        )
      }
      value_to_use <- value
      is_vector <- TRUE
    }

    # 6d. Add column if it doesn't exist
    if (!field_name %in% existing_cols) {
      add_column_query <- paste0(
        "ALTER TABLE ind_meta ADD COLUMN ",
        field_name, " ", db_type
      )
      DBI::dbExecute(pop$db_conn, add_column_query)
    }

    # 6e. Update values
    if (is_vector) {
      # Vector: Update row by row using id_ind
      update_column_vector(pop$db_conn, field_name, value_to_use, db_type)
    } else {
      # Scalar: Update all rows with same value
      update_column_scalar(pop$db_conn, field_name, value_to_use, db_type)
    }
  }

  # 7. Return modified pop object
  invisible(pop)
}


#' Infer DuckDB column type from R value
#'
#' @param value An R value (scalar or vector)
#' @return Character string representing DuckDB type
#' @keywords internal
infer_duckdb_type <- function(value) {

  # A bare NA (untyped logical NA) carries no useful type information.
  # Warn and default to VARCHAR; callers should use NA_real_, NA_integer_, etc.
  if (is.logical(value) && all(is.na(value))) {
    warning(
      "Cannot infer type from NA value. Defaulting to VARCHAR. ",
      "Use typed NA (e.g. NA_real_, NA_integer_, NA_character_) or supply a ",
      "non-NA seed value first to get the correct type.",
      call. = FALSE
    )
    return("VARCHAR")
  }

  # Check R class first — works correctly even when elements are NA,
  # since class(c(NA, 1.5)) is "numeric", class(c(NA, FALSE)) is "logical", etc.
  if (is.logical(value)) {
    return("BOOLEAN")
  } else if (is.integer(value)) {
    return("INTEGER")
  } else if (is.numeric(value)) {
    return("DOUBLE")
  } else if (inherits(value, "Date")) {
    return("DATE")
  } else if (inherits(value, "POSIXct") || inherits(value, "POSIXlt")) {
    return("TIMESTAMP")
  } else if (is.character(value)) {
    return("VARCHAR")
  }

  # Fallback: try first non-NA value (handles exotic/unclassed vectors)
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


#' Validate field name for SQL safety
#'
#' @param field_name Character string field name
#' @param existing_cols Character vector of existing column names
#' @keywords internal
validate_field_name <- function(field_name, existing_cols) {
  validate_sql_identifier(
    field_name,
    what     = "field name",
    reserved = c("id_ind", "id_parent_1", "id_parent_2", "line", "sex")
  )
}


#' Update column with scalar value
#'
#' @param conn DuckDB connection
#' @param field_name Column name
#' @param value Scalar value
#' @param db_type DuckDB type string
#' @keywords internal
update_column_scalar <- function(conn, field_name, value, db_type) {

  # Format value for SQL based on type
  if (is.na(value)) {
    sql_value <- "NULL"
  } else if (db_type == "BOOLEAN") {
    sql_value <- ifelse(value, "TRUE", "FALSE")
  } else if (db_type %in% c("INTEGER", "DOUBLE")) {
    sql_value <- as.character(value)
  } else if (db_type == "DATE") {
    sql_value <- paste0("'", as.character(value), "'")
  } else if (db_type == "TIMESTAMP") {
    sql_value <- paste0("'", format(value, "%Y-%m-%d %H:%M:%S"), "'")
  } else if (db_type == "VARCHAR") {
    # Escape single quotes
    escaped_value <- gsub("'", "''", as.character(value))
    sql_value <- paste0("'", escaped_value, "'")
  }

  # Execute UPDATE
  update_query <- paste0(
    "UPDATE ind_meta SET ", field_name, " = ", sql_value
  )

  DBI::dbExecute(conn, update_query)

  invisible(NULL)
}


#' Update column with vector of values
#'
#' @param conn DuckDB connection
#' @param field_name Column name
#' @param values Vector of values
#' @param db_type DuckDB type string
#' @keywords internal
update_column_vector <- function(conn, field_name, values, db_type) {

  # Get individual IDs (assumes id_ind exists and is unique)
  ind_ids <- DBI::dbGetQuery(
    conn,
    "SELECT id_ind FROM ind_meta ORDER BY ROWID"
  )$id_ind

  # Create temporary data frame
  update_df <- data.frame(
    id_ind = ind_ids,
    new_value = values,
    stringsAsFactors = FALSE
  )

  # Write to temporary table
  temp_table <- paste0("temp_update_", gsub("[^a-zA-Z0-9]", "", field_name))
  DBI::dbWriteTable(conn, temp_table, update_df, overwrite = TRUE)

  # UPDATE using JOIN
  update_query <- paste0(
    "UPDATE ind_meta ",
    "SET ", field_name, " = temp.new_value ",
    "FROM ", temp_table, " AS temp ",
    "WHERE ind_meta.id_ind = temp.id_ind"
  )

  DBI::dbExecute(conn, update_query)

  # Clean up temporary table
  DBI::dbExecute(conn, paste0("DROP TABLE ", temp_table))

  invisible(NULL)
}

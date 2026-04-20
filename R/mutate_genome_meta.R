#' Add or modify columns in genome metadata table
#'
#' @description
#' Adds new columns to the `genome_meta` table or updates existing custom columns.
#' Values can be scalar (applied to all loci) or vectors (one per locus).
#' The function automatically infers DuckDB column types from R data types.
#'
#' Reserved columns (`locus_id`, `locus_name`, `chr`, `chr_name`, `pos_Mb`) cannot be modified.
#'
#' @param pop A `tidybreed_pop` object
#' @param ... Named arguments specifying column names and values.
#'   Each argument should be `field_name = value`, where value is either:
#'   - A scalar (applied to all loci)
#'   - A vector of length equal to number of loci
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
#' - Cannot be reserved columns (locus_id, locus_name, chr, chr_name, pos_Mb)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Initialize genome
#' pop <- initialize_genome(
#'   pop_name = "test",
#'   n_loci = 1000,
#'   n_chr = 10,
#'   chr_len_Mb = 100
#' )
#'
#' # Add scalar fields (same value for all loci)
#' pop <- pop %>%
#'   mutate_genome_meta(
#'     is_QTL = FALSE,
#'     is_50k = FALSE
#'   )
#'
#' # Add vector field (different value per locus)
#' n_loci <- 1000
#' pop <- pop %>%
#'   mutate_genome_meta(
#'     maf = runif(n_loci, 0.01, 0.5),
#'     is_functional = sample(c(TRUE, FALSE), n_loci, replace = TRUE)
#'   )
#'
#' # Update existing field
#' pop <- pop %>%
#'   mutate_genome_meta(is_QTL = TRUE)
#'
#' # View updated metadata
#' get_table(pop, "genome_meta") %>% collect()
#' }
mutate_genome_meta <- function(pop, ...) {

  # 1. Validate inputs
  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)

  # 2. Check genome_meta table exists
  if (!"genome_meta" %in% pop$tables) {
    stop(
      "genome_meta table does not exist. Call initialize_genome() first.",
      call. = FALSE
    )
  }

  # 3. Capture named arguments
  field_values <- list(...)

  if (length(field_values) == 0) {
    warning("No fields specified. Returning population unchanged.", call. = FALSE)
    return(invisible(pop))
  }

  # 4. Get current number of loci
  n_loci <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT COUNT(*) as n FROM genome_meta"
  )$n

  if (n_loci == 0) {
    warning("genome_meta table is empty. No values to update.", call. = FALSE)
    return(invisible(pop))
  }

  # 5. Get existing columns
  existing_cols <- DBI::dbListFields(pop$db_conn, "genome_meta")

  # 6. Process each field
  for (field_name in names(field_values)) {

    value <- field_values[[field_name]]

    # 6a. Validate field name
    validate_field_name_genome(field_name, existing_cols)

    # 6b. Determine if value is scalar or vector
    if (length(value) == 1) {
      # Scalar - will be replicated for all loci
      value_to_use <- value
      is_vector <- FALSE
    } else {
      # Vector - must match number of loci
      if (length(value) != n_loci) {
        stop(
          "Length of '", field_name, "' (", length(value),
          ") does not match number of loci (", n_loci, ")",
          call. = FALSE
        )
      }
      value_to_use <- value
      is_vector <- TRUE
    }

    # 6c. Infer DuckDB type from R type
    db_type <- infer_duckdb_type(value_to_use)

    # 6d. Add column if it doesn't exist
    if (!field_name %in% existing_cols) {
      add_column_query <- paste0(
        "ALTER TABLE genome_meta ADD COLUMN ",
        field_name, " ", db_type
      )
      DBI::dbExecute(pop$db_conn, add_column_query)
    }

    # 6e. Update values
    if (is_vector) {
      # Vector: Update row by row using locus_id
      update_column_vector_genome(pop$db_conn, field_name, value_to_use, db_type)
    } else {
      # Scalar: Update all rows with same value
      update_column_scalar_genome(pop$db_conn, field_name, value_to_use, db_type)
    }
  }

  # 7. Return modified pop object
  invisible(pop)
}


#' Validate field name for genome_meta table
#'
#' @param field_name Character string field name
#' @param existing_cols Character vector of existing column names
#' @keywords internal
validate_field_name_genome <- function(field_name, existing_cols) {
  validate_sql_identifier(
    field_name,
    what     = "field name",
    reserved = c("locus_id", "locus_name", "chr", "chr_name", "pos_Mb")
  )
}


#' Update column with scalar value in genome_meta
#'
#' @param conn DuckDB connection
#' @param field_name Column name
#' @param value Scalar value
#' @param db_type DuckDB type string
#' @keywords internal
update_column_scalar_genome <- function(conn, field_name, value, db_type) {

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
    "UPDATE genome_meta SET ", field_name, " = ", sql_value
  )

  DBI::dbExecute(conn, update_query)

  invisible(NULL)
}


#' Update column with vector of values in genome_meta
#'
#' @param conn DuckDB connection
#' @param field_name Column name
#' @param values Vector of values
#' @param db_type DuckDB type string
#' @keywords internal
update_column_vector_genome <- function(conn, field_name, values, db_type) {

  # Get locus IDs (assumes locus_id exists and is unique)
  locus_ids <- DBI::dbGetQuery(
    conn,
    "SELECT locus_id FROM genome_meta ORDER BY locus_id"
  )$locus_id

  # Create temporary data frame
  update_df <- data.frame(
    locus_id = locus_ids,
    new_value = values,
    stringsAsFactors = FALSE
  )

  # Write to temporary table
  temp_table <- paste0("temp_update_genome_", gsub("[^a-zA-Z0-9]", "", field_name))
  DBI::dbWriteTable(conn, temp_table, update_df, overwrite = TRUE)

  # UPDATE using JOIN
  update_query <- paste0(
    "UPDATE genome_meta ",
    "SET ", field_name, " = temp.new_value ",
    "FROM ", temp_table, " AS temp ",
    "WHERE genome_meta.locus_id = temp.locus_id"
  )

  DBI::dbExecute(conn, update_query)

  # Clean up temporary table
  DBI::dbExecute(conn, paste0("DROP TABLE ", temp_table))

  invisible(NULL)
}

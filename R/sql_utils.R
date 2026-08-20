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
  "_schema_meta"   = c("id_schema_meta", "object_type", "table_name",
                       "column_name", "description", "notes"),
  ind_meta         = c("id_ind", "id_parent_1", "id_parent_2", "line_name", "sex", "ploidy", "replicate"),
  genome_meta      = c("locus_id", "locus_name", "chr", "chr_name", "pos_bp"),
  genome_map       = c("id_genome_map", "locus_id", "locus_name", "sex",
                       "line_name", "map_name", "pos_cM"),
  genome_effects   = c("id_genome_effect", "locus_name", "line_name", "trait_name",
                       "genome_effect_type", "genome_value", "base_allele_freq"),
  ind_haplotype    = c("id_ind", "parent_origin", "strand", "line_origin",
                       "locus_id", "locus_name", "allele"),
  ind_genotype     = c("id_ind", "locus_id", "locus_name", "dosage_value"),
  ind_crossover    = c("id_crossover", "id_ind", "parent_origin", "chr",
                       "chr_name", "pos_cM"),
  chr_inheritance  = c("chr_name", "offspring_sex", "line_name",
                       "from_parent_1", "from_parent_2"),
  chr_recombination = c("chr_name", "parent_sex", "line_name", "recombines"),
  ind_phenotype    = c("id_phenotype", "id_ind", "phenotype_name", "pheno_value", "pheno_number",
                       "liability_value", "cat_name", "replicate"),
  ind_tbv          = c("id_tbv", "id_ind", "trait_name", "tbv_value", "replicate"),
  ind_ebv          = c("id_ebv", "id_ind", "trait_name", "model", "ebv_value", "acc", "se", "eval_number", "replicate"),
  trait_meta       = c("id_trait", "trait_name", "description", "units",
                       "expressed_parent", "target_add_mean"),
  phenotype_meta   = c("id_phenotype_meta", "phenotype_name", "type", "mean",
                       "expressed_sex", "repeatable", "min_value", "max_value",
                       "prevalence", "thresholds", "cat_values", "cat_names",
                       "store_liability", "missing_component_action",
                       "formula_tbv", "formula"),
  trait_effects    = c("phenotype_name", "effect_name", "effect_class", "source_column",
                       "source_table", "distribution", "levels_json", "slope",
                       "center", "value", "poly_order", "null_class_action"),
  trait_var_comp   = c("id_trait_var_comp", "effect_name", "trait_name_1", "trait_name_2", "cov_value"),
  phenotype_var_comp = c("id_phenotype_var_comp", "effect_name", "phenotype_name_1",
                         "phenotype_name_2", "cov_value", "condition_column",
                         "condition_table", "condition_level", "weight_type", "poly_order"),
  index_meta       = c("id_index_name", "index_name", "trait_name", "index_weight", "economic_weight"),
  ind_index        = c("id_index", "id_ind", "index_name", "index_number", "index_value", "replicate"),
  ind_true_index   = c("id_true_index", "id_ind", "index_name", "weight_type", "true_index_value", "replicate")
)


#' Primary key column per table (used for vector updates and filtered updates)
#'
#' @keywords internal
TABLE_PRIMARY_KEYS <- list(
  ind_meta         = "id_ind",
  genome_meta      = "locus_id",
  genome_map       = "id_genome_map",
  ind_crossover    = "id_crossover",
  genome_effects   = "id_genome_effect",
  ind_phenotype    = "id_phenotype",
  ind_tbv          = "id_tbv",
  ind_ebv          = "id_ebv",
  trait_meta       = "id_trait",
  trait_var_comp     = "id_trait_var_comp",
  phenotype_var_comp = "id_phenotype_var_comp",
  index_meta         = "id_index_name",
  ind_index        = "id_index",
  ind_true_index   = "id_true_index"
)


#' Composite row keys per table (used for exact single-table DELETE)
#'
#' For tables with a single PK this is length-1. For composite-PK tables all
#' key columns are listed so that `remove_rows()` can target exact rows rather
#' than over-deleting by `id_ind` alone.
#'
#' @keywords internal
TABLE_ROW_KEYS <- list(
  ind_meta         = "id_ind",
  genome_meta      = "locus_id",
  genome_map       = "id_genome_map",
  genome_effects   = "id_genome_effect",
  ind_phenotype    = "id_phenotype",
  trait_meta       = "id_trait",
  ind_tbv          = c("id_ind", "trait_name"),
  ind_ebv          = c("id_ind", "trait_name", "model", "eval_number"),
  ind_index        = c("id_ind", "index_name", "index_number"),
  ind_true_index   = c("id_ind", "index_name", "weight_type"),
  ind_haplotype    = c("id_ind", "parent_origin", "strand", "locus_id"),
  ind_genotype     = c("id_ind", "locus_id"),
  ind_crossover    = "id_crossover",
  chr_inheritance   = c("chr_name", "offspring_sex", "line_name"),
  chr_recombination = c("chr_name", "parent_sex", "line_name"),
  trait_effects    = c("phenotype_name", "effect_name"),
  trait_var_comp     = "id_trait_var_comp",
  phenotype_var_comp = "id_phenotype_var_comp",
  index_meta         = "id_index_name"
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


#' Get the next ID for a BIGINT primary-key sequence (overflow-safe).
#'
#' Like [next_int_id()] but returns a **double** rather than an `integer`, so it
#' is safe for `BIGINT` PK columns whose values can exceed `.Machine$integer.max`
#' (`2^31 - 1`). `as.integer()` would return `NA` past that ceiling; a double is
#' exact to `2^53`. Used for `ind_crossover.id_crossover`, which accumulates far
#' more rows than any `INTEGER`-keyed table. DuckDB already returns `BIGINT`
#' `MAX()` as an R double, so no precision is lost on the read.
#'
#' For multi-row inserts, pre-generate a range:
#'   `seq.int(next_row_id(conn, tbl, col), length.out = n)`
#'
#' @param conn A DBI connection.
#' @param table Character. Table name.
#' @param id_col Character. BIGINT PK column name.
#' @return Numeric (double) scalar.
#' @keywords internal
next_row_id <- function(conn, table, id_col) {
  as.numeric(DBI::dbGetQuery(
    conn,
    paste0("SELECT COALESCE(MAX(", id_col, "), 0) + 1 AS nxt FROM ", table)
  )$nxt)
}


#' System-managed tables that users must not overwrite with `define_table()`
#'
#' @keywords internal
SYSTEM_TABLES <- c(
  "_schema_meta",
  "genome_meta", "genome_map", "genome_effects",
  "ind_haplotype", "ind_genotype", "ind_crossover",
  "chr_inheritance", "chr_recombination",
  "founder_haplotypes",
  "ind_meta", "ind_phenotype", "ind_tbv", "ind_ebv",
  "trait_meta", "trait_effects", "trait_var_comp", "trait_random_effects",
  "phenotype_meta", "phenotype_components", "phenotype_var_comp",
  "index_meta", "ind_index", "ind_true_index"
)


#' Tables that contain the `id_ind` column (used for cross-table deletion)
#'
#' @keywords internal
IND_TABLE_ID_IND_COLS <- c(
  "ind_meta", "ind_phenotype", "ind_tbv", "ind_ebv",
  "ind_index", "ind_true_index",
  "ind_haplotype", "ind_genotype", "ind_crossover"
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


#' Quote and escape a character vector for a SQL `IN (...)` list
#'
#' Escapes embedded single quotes (quote-doubling, matching
#' `format_sql_value()`'s VARCHAR branch) and returns the comma-joined,
#' quoted fragment only — callers still write `paste0("... IN (", x, ")")`
#' around the result.
#'
#' @param values Character vector to escape and format for a SQL `IN (...)`
#'   list. Must be non-empty and free of `NA` (a bare `NA` would otherwise
#'   silently become the literal string `'NA'` in SQL).
#' @param what Human-readable role label used in error messages
#'   (e.g. `"locus name"`, `"individual ID"`).
#' @return A single character string: comma-separated, single-quoted,
#'   escaped values.
#' @keywords internal
sql_in_list <- function(values, what = "value") {
  if (!is.character(values) || length(values) == 0L) {
    stop("No ", what, "s provided.", call. = FALSE)
  }
  if (anyNA(values)) {
    stop("NA not allowed in ", what, " list.", call. = FALSE)
  }
  escaped <- gsub("'", "''", values)
  paste0("'", escaped, "'", collapse = ", ")
}

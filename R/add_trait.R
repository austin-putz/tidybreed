#' Add a trait definition to the population
#'
#' @description
#' Creates one row in the `trait_meta` table describing a trait's genetic and
#' phenotypic architecture. The row captures target variance components,
#' trait type (continuous / count / binary / categorical), expression rules
#' (sex-limited, parent-of-origin), and index / economic weights.
#'
#' `add_trait()` does not touch the genome or sample QTL. Those are separate
#' steps handled by [define_qtl()] and [add_additive_effects()]. For the common
#' one-off case, use [add_trait_simple()] which chains all three together.
#'
#' On first call, creates the trait-layer tables (`trait_meta`, `trait_effects`,
#' `trait_random_effects`, `ind_phenotype`, `ind_tbv`, `ind_ebv`).
#' Variance and covariance data is stored in the separate `trait_effect_cov`
#' table (created on first call to [add_effect_cov_matrix()] or when
#' `target_add_var` / `residual_var` are supplied here).
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character. Trait name; used as the primary key and as the
#'   suffix for `is_QTL_{trait_name}` and `add_{trait_name}` columns in `genome_meta`.
#'   Must be a valid SQL identifier.
#' @param description Character. Free-text description of the trait.
#' @param units Character. Measurement units, e.g. `"kg"`, `"count"`.
#' @param trait_type Character. One of `"continuous"`, `"count"`,
#'   `"categorical"`. (`"binary"` has been removed — use `"categorical"` with
#'   one threshold or `prevalence`, and `cat_values = c(0, 1)` for 0/1
#'   encoding.)
#' @param repeatable Logical. Whether an individual can have multiple records
#'   for this trait (e.g. test-day yield).
#' @param recorded_on Character. Which animal records the trait:
#'   `"self"`, `"dam"`, `"sire"`, or `"offspring_mean"`.
#' @param expressed_sex Character. Sex that expresses the trait: `"both"`,
#'   `"M"`, or `"F"`.
#' @param expressed_parent Character. Parent-of-origin expression:
#'   `"both"`, `"parent_1"` (paternal), `"parent_2"` (maternal). Imprinted
#'   traits use only the haplotype from the specified parent.
#' @param target_add_mean Numeric. Target mean additive genetic value for the
#'   base population. Used as the intercept in the phenotype model. Defaults to
#'   0 so that `E[TBV] = 0` when TBV is centered on base allele frequencies.
#' @param target_add_var Numeric. Target additive-genetic variance. When
#'   provided, written to `trait_effect_cov` as a diagonal entry under
#'   `effect_name = "gen_add"`. Used by [add_additive_effects()] to rescale effects.
#'   If already set via [add_effect_cov_matrix()], leave `NULL`.
#' @param residual_var Numeric. Residual variance. When provided, written to
#'   `trait_effect_cov` under `effect_name = "residual"`. Used by
#'   [add_phenotype()]. If already set via [add_effect_cov_matrix()], leave
#'   `NULL`.
#' @param min_value,max_value Numeric. Clipping bounds for count traits. `NA`
#'   means no limit.
#' @param prevalence Numeric between 0 and 1. For categorical traits with a
#'   single threshold (two categories), the fraction expected above the
#'   threshold. The threshold is computed at phenotype time as
#'   `target_add_mean + qnorm(1 - prevalence) * sqrt(VA + VR)`.
#'   Mutually exclusive with `thresholds`.
#' @param thresholds Numeric vector of length K−1 for K ordered categories.
#'   Liability cutpoints in ascending order separating the categories. Stored
#'   as a comma-separated string in `trait_meta`. Mutually exclusive with
#'   `prevalence`.
#' @param cat_values Numeric vector of length K (one per category) giving the
#'   phenotype value stored in `ind_phenotype` for each category. E.g.
#'   `c(0, 1)` for standard binary encoding. Defaults to `1, 2, ..., K` when
#'   `NULL`.
#' @param cat_names Character vector of length K. Human-readable label per
#'   category, e.g. `c("Alive", "Dead")`. Written to the reserved
#'   `cat_name` column in `ind_phenotype`. Must not contain commas.
#' @param store_liability Logical. When `TRUE`, the underlying liability value
#'   is written to the reserved `liability_value` column in `ind_phenotype`
#'   alongside the observed category value.
#' @param index_weight Numeric. Weight in a downstream selection index.
#' @param economic_value Numeric. Economic value per unit of the trait.
#' @param overwrite Logical. If `TRUE` and a trait with the same name already
#'   exists, replace its `trait_meta` row. Associated rows in `trait_effects`
#'   are cleared for that trait.
#'
#' @return The modified `tidybreed_pop` (invisibly). Assign the result back.
#'
#' @seealso [define_qtl()], [add_additive_effects()], [add_effect_fixed_class()],
#'   [add_effect_fixed_cov()], [add_effect_random()], [add_phenotype()],
#'   [add_trait_simple()]
#'
#' @examples
#' \dontrun{
#' pop <- pop |>
#'   add_trait(
#'     trait_name     = "ADG",
#'     trait_type     = "continuous",
#'     units          = "g/day",
#'     target_add_mean = 0,
#'     target_add_var = 0.25,
#'     residual_var   = 0.75
#'   )
#' }
#' @export
add_trait <- function(pop,
                      trait_name,
                      description      = NA_character_,
                      units            = NA_character_,
                      trait_type       = c("continuous", "count", "categorical"),
                      repeatable       = FALSE,
                      recorded_on      = c("self", "dam", "sire",
                                            "offspring_mean"),
                      expressed_sex    = c("both", "M", "F"),
                      expressed_parent = c("both", "parent_1", "parent_2"),
                      target_add_mean  = 0,
                      target_add_var   = NULL,
                      residual_var     = NULL,
                      min_value        = NA_real_,
                      max_value        = NA_real_,
                      prevalence       = NA_real_,
                      thresholds       = NULL,
                      cat_values       = NULL,
                      cat_names        = NULL,
                      store_liability  = FALSE,
                      index_weight     = 0,
                      economic_value   = 0,
                      overwrite        = FALSE) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)

  stopifnot(is.character(trait_name), length(trait_name) == 1, nchar(trait_name) > 0)
  validate_trait_name(trait_name)

  # Catch removed "binary" type before match.arg drops it
  trait_type_raw <- trait_type[1]
  if (identical(trait_type_raw, "binary")) {
    stop(
      "trait_type = 'binary' has been removed. ",
      "Use trait_type = 'categorical' with thresholds = c(<value>) or ",
      "prevalence = <p> (for a 2-category trait). ",
      "To keep 0/1 encoding supply cat_values = c(0, 1).",
      call. = FALSE
    )
  }

  trait_type       <- match.arg(trait_type)
  recorded_on      <- match.arg(recorded_on)
  expressed_sex    <- match.arg(expressed_sex)
  expressed_parent <- match.arg(expressed_parent)

  if (trait_type == "categorical") {
    has_thresholds <- !is.null(thresholds) && length(thresholds) >= 1 &&
                      !all(is.na(thresholds))
    has_prevalence <- !is.na(prevalence) && prevalence > 0 && prevalence < 1

    if (!has_thresholds && !has_prevalence) {
      stop(
        "Categorical traits require `thresholds` (numeric vector of cutpoints) ",
        "OR `prevalence` (for a 2-category trait).",
        call. = FALSE
      )
    }
    if (has_thresholds && has_prevalence) {
      stop("Supply `thresholds` OR `prevalence`, not both.", call. = FALSE)
    }

    n_cats <- if (has_thresholds) length(thresholds) + 1L else 2L

    if (!is.null(cat_values)) {
      if (!is.numeric(cat_values) || length(cat_values) != n_cats) {
        stop(
          "`cat_values` must be a numeric vector of length ", n_cats,
          " (one value per category).",
          call. = FALSE
        )
      }
    }
    if (!is.null(cat_names)) {
      ref_len <- if (!is.null(cat_values)) length(cat_values) else n_cats
      if (!is.character(cat_names) || length(cat_names) != ref_len) {
        stop(
          "`cat_names` must be a character vector of length ", ref_len, ".",
          call. = FALSE
        )
      }
      if (any(grepl(",", cat_names, fixed = TRUE))) {
        stop("`cat_names` entries must not contain commas.", call. = FALSE)
      }
    }
  }

  pop <- ensure_trait_tables(pop)

  # Check for existing trait
  existing <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM trait_meta WHERE trait_name = '",
           trait_name, "'")
  )$n

  if (existing > 0 && !overwrite) {
    stop(
      "Trait '", trait_name, "' already exists. Use `overwrite = TRUE` to replace.",
      call. = FALSE
    )
  }

  if (existing > 0 && overwrite) {
    DBI::dbExecute(pop$db_conn,
      paste0("DELETE FROM trait_meta WHERE trait_name = '", trait_name, "'"))
    DBI::dbExecute(pop$db_conn,
      paste0("DELETE FROM trait_effects WHERE trait_name = '", trait_name, "'"))
  }

  thresholds_str <- if (is.null(thresholds) || all(is.na(thresholds))) {
    NA_character_
  } else {
    paste(thresholds, collapse = ",")
  }

  if (!is.null(target_add_var)) {
    if (!is.numeric(target_add_var) || length(target_add_var) != 1 ||
        is.na(target_add_var) || target_add_var < 0) {
      stop("`target_add_var` must be a non-negative number.", call. = FALSE)
    }
    pop <- write_effect_cov_diagonal(pop, "gen_add", trait_name,
                                     as.numeric(target_add_var))
  }
  if (!is.null(residual_var)) {
    if (!is.numeric(residual_var) || length(residual_var) != 1 ||
        is.na(residual_var) || residual_var < 0) {
      stop("`residual_var` must be a non-negative number.", call. = FALSE)
    }
    pop <- write_effect_cov_diagonal(pop, "residual", trait_name,
                                     as.numeric(residual_var))
  }

  row <- tibble::tibble(
    id_trait         = next_int_id(pop$db_conn, "trait_meta", "id_trait"),
    trait_name       = trait_name,
    description      = as.character(description),
    units            = as.character(units),
    trait_type       = trait_type,
    repeatable       = as.logical(repeatable),
    recorded_on      = recorded_on,
    expressed_sex    = expressed_sex,
    expressed_parent = expressed_parent,
    target_add_mean  = as.numeric(target_add_mean),
    min_value        = as.numeric(min_value),
    max_value        = as.numeric(max_value),
    prevalence       = as.numeric(prevalence),
    thresholds       = thresholds_str,
    index_weight     = as.numeric(index_weight),
    economic_value   = as.numeric(economic_value),
    cat_values       = if (is.null(cat_values)) NA_character_
                       else paste(cat_values, collapse = ","),
    cat_names        = if (is.null(cat_names)) NA_character_
                       else paste(cat_names, collapse = ","),
    store_liability  = as.logical(store_liability)
  )

  DBI::dbWriteTable(pop$db_conn, "trait_meta", row, append = TRUE)

  message("Added trait '", trait_name, "' (type: ", trait_type, ").")

  invisible(pop)
}


#' Ensure the trait-layer tables exist in the database
#'
#' @description
#' Creates `trait_meta`, `trait_effects`, `trait_random_effects`,
#' `ind_phenotype`, `ind_tbv`, and `ind_ebv` if they are not already present.
#' Idempotent: safe to call multiple times. Variance/covariance data lives in
#' `trait_effect_cov`, created lazily by [add_effect_cov_matrix()].
#'
#' @param pop A `tidybreed_pop` object.
#' @return The `tidybreed_pop` object with `$tables` updated.
#' @keywords internal
ensure_trait_tables <- function(pop) {

  existing <- DBI::dbListTables(pop$db_conn)

  ddl <- list(
    trait_meta = "
      CREATE TABLE trait_meta (
        id_trait         INTEGER PRIMARY KEY,
        trait_name       VARCHAR UNIQUE,
        description      VARCHAR,
        units            VARCHAR,
        trait_type       VARCHAR,
        repeatable       BOOLEAN,
        recorded_on      VARCHAR,
        expressed_sex    VARCHAR,
        expressed_parent VARCHAR,
        target_add_mean  DOUBLE,
        min_value        DOUBLE,
        max_value        DOUBLE,
        prevalence       DOUBLE,
        thresholds       VARCHAR,
        index_weight     DOUBLE,
        economic_value   DOUBLE,
        cat_values       VARCHAR,
        cat_names        VARCHAR,
        store_liability  BOOLEAN
      )
    ",
    trait_effects = "
      CREATE TABLE trait_effects (
        trait_name    VARCHAR,
        effect_name   VARCHAR,
        effect_class  VARCHAR,
        source_column VARCHAR,
        source_table  VARCHAR,
        distribution  VARCHAR,
        levels_json   VARCHAR,
        slope         DOUBLE,
        center        DOUBLE,
        value         DOUBLE,
        PRIMARY KEY (trait_name, effect_name)
      )
    ",
    trait_random_effects = "
      CREATE TABLE trait_random_effects (
        trait_name   VARCHAR,
        effect_name  VARCHAR,
        level        VARCHAR,
        draw_value   DOUBLE,
        date_sampled DATE,
        PRIMARY KEY (trait_name, effect_name, level)
      )
    ",
    ind_phenotype = "
      CREATE TABLE ind_phenotype (
        id_phenotype INTEGER PRIMARY KEY,
        id_ind       VARCHAR,
        trait_name   VARCHAR,
        pheno_value  DOUBLE,
        pheno_number INTEGER
      )
    ",
    ind_tbv = "
      CREATE TABLE ind_tbv (
        id_tbv     INTEGER PRIMARY KEY,
        id_ind     VARCHAR,
        trait_name VARCHAR,
        tbv_value  DOUBLE,
        UNIQUE (id_ind, trait_name)
      )
    ",
    ind_ebv = "
      CREATE TABLE ind_ebv (
        id_ebv      INTEGER PRIMARY KEY,
        id_ind      VARCHAR,
        trait_name  VARCHAR,
        model       VARCHAR,
        ebv_value   DOUBLE,
        acc         DOUBLE,
        se          DOUBLE,
        eval_number INTEGER,
        UNIQUE (id_ind, trait_name, model, eval_number)
      )
    ",
    index_meta = "
      CREATE TABLE index_meta (
        id_index_name INTEGER PRIMARY KEY,
        index_name    VARCHAR,
        trait_name    VARCHAR,
        index_weight  DOUBLE,
        UNIQUE (index_name, trait_name)
      )
    ",
    ind_index = "
      CREATE TABLE ind_index (
        id_index     INTEGER PRIMARY KEY,
        id_ind       VARCHAR,
        index_name   VARCHAR,
        index_number INTEGER,
        index_value  DOUBLE
      )
    "
  )

  for (tbl in names(ddl)) {
    if (!tbl %in% existing) {
      DBI::dbExecute(pop$db_conn, ddl[[tbl]])
    }
  }

  # Migrate ind_tbv: drop date_calc if present (no longer used)
  if ("ind_tbv" %in% existing) {
    tbv_cols <- DBI::dbListFields(pop$db_conn, "ind_tbv")
    if ("date_calc" %in% tbv_cols)
      DBI::dbExecute(pop$db_conn, "ALTER TABLE ind_tbv DROP COLUMN date_calc")
  }

  # Migrate ind_ebv: drop date_calc, add eval_number if needed
  if ("ind_ebv" %in% existing) {
    ebv_cols <- DBI::dbListFields(pop$db_conn, "ind_ebv")
    if ("date_calc" %in% ebv_cols)
      DBI::dbExecute(pop$db_conn, "ALTER TABLE ind_ebv DROP COLUMN date_calc")
    if (!"eval_number" %in% ebv_cols)
      DBI::dbExecute(pop$db_conn,
        "ALTER TABLE ind_ebv ADD COLUMN eval_number INTEGER DEFAULT 1")
  }

  # Migrate trait_effects: add new columns if the table already existed
  if ("trait_effects" %in% existing) {
    te_cols <- DBI::dbListFields(pop$db_conn, "trait_effects")
    for (new_col in c("source_table VARCHAR", "slope DOUBLE", "center DOUBLE")) {
      col_name <- strsplit(new_col, " ")[[1]][1]
      if (!col_name %in% te_cols) {
        DBI::dbExecute(pop$db_conn,
          paste0("ALTER TABLE trait_effects ADD COLUMN ", new_col))
      }
    }
  }

  # Migrate trait_meta: add categorical columns if the table already existed
  if ("trait_meta" %in% existing) {
    tm_cols <- DBI::dbListFields(pop$db_conn, "trait_meta")
    for (new_col in c("cat_values VARCHAR", "cat_names VARCHAR",
                      "store_liability BOOLEAN")) {
      col_name <- strsplit(new_col, " ")[[1]][1]
      if (!col_name %in% tm_cols) {
        DBI::dbExecute(pop$db_conn,
          paste0("ALTER TABLE trait_meta ADD COLUMN ", new_col))
      }
    }
  }

  pop$tables <- unique(c(pop$tables, names(ddl)))
  pop
}


#' Validate a trait name for use as SQL identifier and column suffix
#'
#' @param name Character trait name.
#' @keywords internal
validate_trait_name <- function(name) {
  validate_sql_identifier(name, what = "trait name")
}

#' Add a trait definition to the population
#'
#' @description
#' Creates one row in the `trait_meta` table describing a trait's genetic and
#' phenotypic architecture. The row captures target variance components,
#' trait type (continuous / count / binary / categorical), expression rules
#' (sex-limited, parent-of-origin), and index / economic weights.
#'
#' `add_trait()` does not touch the genome or sample QTL. Those are separate
#' steps handled by [define_qtl()] and [set_qtl_effects()]. For the common
#' one-off case, use [add_trait_simple()] which chains all three together.
#'
#' On first call, creates the six trait-layer tables (`trait_meta`,
#' `trait_effects`, `trait_residual_cov`, `ind_phenotype`, `ind_tbv`,
#' `ind_ebv`).
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character. Trait name; used as the primary key and as the
#'   suffix for `is_QTL_{trait_name}` and `add_{trait_name}` columns in `genome_meta`.
#'   Must be a valid SQL identifier.
#' @param description Character. Free-text description of the trait.
#' @param units Character. Measurement units, e.g. `"kg"`, `"count"`.
#' @param trait_type Character. One of `"continuous"`, `"count"`, `"binary"`,
#'   `"categorical"`.
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
#' @param target_add_var Numeric. Target additive-genetic variance used by
#'   [set_qtl_effects()] to rescale sampled effects.
#' @param residual_var Numeric. Residual variance used by [add_phenotype()].
#' @param min_value,max_value Numeric. Clipping bounds for count traits. `NA`
#'   means no limit.
#' @param prevalence Numeric between 0 and 1. For binary traits, the fraction
#'   expected to be affected; determines the liability threshold.
#' @param thresholds Numeric vector. For categorical traits, liability
#'   cutpoints separating ordered levels. Stored as a comma-separated string.
#' @param index_weight Numeric. Weight in a downstream selection index.
#' @param economic_value Numeric. Economic value per unit of the trait.
#' @param overwrite Logical. If `TRUE` and a trait with the same name already
#'   exists, replace its `trait_meta` row. Associated rows in `trait_effects`
#'   and `trait_residual_cov` are cleared for that trait.
#'
#' @return The modified `tidybreed_pop` (invisibly). Assign the result back.
#'
#' @seealso [define_qtl()], [set_qtl_effects()], [add_trait_covariate()],
#'   [add_phenotype()], [add_trait_simple()]
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
                      trait_type       = c("continuous", "count",
                                            "binary", "categorical"),
                      repeatable       = FALSE,
                      recorded_on      = c("self", "dam", "sire",
                                            "offspring_mean"),
                      expressed_sex    = c("both", "M", "F"),
                      expressed_parent = c("both", "parent_1", "parent_2"),
                      target_add_mean  = 0,
                      target_add_var   = 1,
                      residual_var     = 1,
                      min_value        = NA_real_,
                      max_value        = NA_real_,
                      prevalence       = NA_real_,
                      thresholds       = NA_real_,
                      index_weight     = 0,
                      economic_value   = 0,
                      overwrite        = FALSE) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)

  stopifnot(is.character(trait_name), length(trait_name) == 1, nchar(trait_name) > 0)
  validate_trait_name(trait_name)

  trait_type       <- match.arg(trait_type)
  recorded_on      <- match.arg(recorded_on)
  expressed_sex    <- match.arg(expressed_sex)
  expressed_parent <- match.arg(expressed_parent)

  if (trait_type == "binary") {
    if (is.na(prevalence) || prevalence <= 0 || prevalence >= 1) {
      stop(
        "Binary traits require `prevalence` strictly between 0 and 1.",
        call. = FALSE
      )
    }
  }
  if (trait_type == "categorical") {
    if (length(thresholds) < 1 || all(is.na(thresholds))) {
      stop(
        "Categorical traits require at least one value in `thresholds`.",
        call. = FALSE
      )
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
    DBI::dbExecute(pop$db_conn,
      paste0("DELETE FROM trait_residual_cov WHERE trait_1 = '", trait_name,
             "' OR trait_2 = '", trait_name, "'"))
  }

  thresholds_str <- if (all(is.na(thresholds))) {
    NA_character_
  } else {
    paste(thresholds, collapse = ",")
  }

  row <- tibble::tibble(
    trait_name       = trait_name,
    description      = as.character(description),
    units            = as.character(units),
    trait_type       = trait_type,
    repeatable       = as.logical(repeatable),
    recorded_on      = recorded_on,
    expressed_sex    = expressed_sex,
    expressed_parent = expressed_parent,
    target_add_mean  = as.numeric(target_add_mean),
    target_add_var   = as.numeric(target_add_var),
    residual_var     = as.numeric(residual_var),
    min_value        = as.numeric(min_value),
    max_value        = as.numeric(max_value),
    prevalence       = as.numeric(prevalence),
    thresholds       = thresholds_str,
    index_weight     = as.numeric(index_weight),
    economic_value   = as.numeric(economic_value)
  )

  DBI::dbWriteTable(pop$db_conn, "trait_meta", row, append = TRUE)

  message(
    "Added trait '", trait_name, "' (type: ", trait_type,
    ", target_add_var: ", target_add_var,
    ", residual_var: ", residual_var, ")"
  )

  invisible(pop)
}


#' Ensure the trait-layer tables exist in the database
#'
#' @description
#' Creates `trait_meta`, `trait_effects`, `trait_residual_cov`,
#' `ind_phenotype`, `ind_tbv`, and `ind_ebv` if they are not already present.
#' Idempotent: safe to call multiple times.
#'
#' @param pop A `tidybreed_pop` object.
#' @return The `tidybreed_pop` object with `$tables` updated.
#' @keywords internal
ensure_trait_tables <- function(pop) {

  existing <- DBI::dbListTables(pop$db_conn)

  ddl <- list(
    trait_meta = "
      CREATE TABLE trait_meta (
        trait_name       VARCHAR PRIMARY KEY,
        description      VARCHAR,
        units            VARCHAR,
        trait_type       VARCHAR,
        repeatable       BOOLEAN,
        recorded_on      VARCHAR,
        expressed_sex    VARCHAR,
        expressed_parent VARCHAR,
        target_add_mean  DOUBLE,
        target_add_var   DOUBLE,
        residual_var     DOUBLE,
        min_value        DOUBLE,
        max_value        DOUBLE,
        prevalence       DOUBLE,
        thresholds       VARCHAR,
        index_weight     DOUBLE,
        economic_value   DOUBLE
      )
    ",
    trait_effects = "
      CREATE TABLE trait_effects (
        trait_name    VARCHAR,
        effect_name   VARCHAR,
        effect_class  VARCHAR,
        source_column VARCHAR,
        distribution  VARCHAR,
        variance      DOUBLE,
        levels_json   VARCHAR,
        value         DOUBLE,
        PRIMARY KEY (trait_name, effect_name)
      )
    ",
    trait_residual_cov = "
      CREATE TABLE trait_residual_cov (
        trait_1 VARCHAR,
        trait_2 VARCHAR,
        cov     DOUBLE,
        PRIMARY KEY (trait_1, trait_2)
      )
    ",
    ind_phenotype = "
      CREATE TABLE ind_phenotype (
        id_record     VARCHAR PRIMARY KEY,
        id_ind        VARCHAR,
        trait_name    VARCHAR,
        value         DOUBLE,
        env           VARCHAR,
        rep           INTEGER,
        date_measured DATE
      )
    ",
    ind_tbv = "
      CREATE TABLE ind_tbv (
        id_ind     VARCHAR,
        trait_name VARCHAR,
        tbv        DOUBLE,
        date_calc  DATE,
        PRIMARY KEY (id_ind, trait_name)
      )
    ",
    ind_ebv = "
      CREATE TABLE ind_ebv (
        id_ind     VARCHAR,
        trait_name VARCHAR,
        model      VARCHAR,
        ebv        DOUBLE,
        acc        DOUBLE,
        se         DOUBLE,
        date_calc  DATE,
        PRIMARY KEY (id_ind, trait_name, model)
      )
    "
  )

  for (tbl in names(ddl)) {
    if (!tbl %in% existing) {
      DBI::dbExecute(pop$db_conn, ddl[[tbl]])
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

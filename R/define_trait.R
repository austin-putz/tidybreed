#' Define a genetic component trait
#'
#' @description
#' Creates one row in `trait_meta` describing a **genetic component trait**: a
#' quantity with QTL effects in `genome_effects`, TBVs in `ind_tbv`, and
#' additive genetic variance in `trait_effect_cov`. Contains no
#' phenotype-level information.
#'
#' To register the **observed phenotype** that individuals receive records for,
#' call [define_phenotype()] after this function. For the common one-off case
#' use [define_trait_simple()], which chains both steps together.
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character. Unique identifier for this genetic component
#'   trait. Must be a valid SQL identifier.
#' @param target_add_var Numeric. Target additive genetic variance. Written to
#'   `trait_effect_cov` as a diagonal entry under `effect_name = "gen_add"`.
#'   Used by [define_additive_effects()] to rescale effects. If already set via
#'   [define_effect_cov_matrix()], leave `NULL`.
#' @param target_add_mean Numeric. TBV centering mean for the base population.
#'   Default `0`; `E[TBV] = 0` when TBVs are centered on base allele
#'   frequencies. The phenotypic population mean (intercept) is set separately
#'   in [define_phenotype()].
#' @param expressed_parent Character. Parent-of-origin expression: `"both"`
#'   (default), `"parent_1"` (paternal), or `"parent_2"` (maternal). Imprinted
#'   traits use only the haplotype from the specified parent when computing TBVs.
#' @param description Character. Free-text description of the trait.
#' @param units Character. Measurement units, e.g. `"kg"`, `"count"`.
#' @param overwrite Logical. If `TRUE` and a trait with the same name already
#'   exists, replace its `trait_meta` row and clear associated `trait_effects`
#'   rows. Default `FALSE` errors if the trait already exists.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_phenotype()], [define_additive_effects()],
#'   [define_effect_fixed_class()], [define_effect_fixed_cov()],
#'   [define_effect_random()], [add_phenotype()], [define_trait_simple()]
#'
#' @examples
#' \dontrun{
#' # Simple genetic component trait:
#' pop <- pop |>
#'   define_trait("ADG", target_add_var = 100, units = "g/day")
#'
#' # Maternal component traits (no define_phenotype call needed for WWD/WWM):
#' pop <- pop |>
#'   define_trait("WWD", target_add_var = 200) |>
#'   define_trait("WWM", target_add_var = 80)
#'
#' # Then define the observed composite phenotype:
#' pop <- pop |>
#'   define_phenotype("WW", trait_type = "continuous", mean = 230,
#'     residual_var = 180,
#'     components = tibble::tribble(
#'       ~source_trait_name, ~contributor_type,
#'       "WWD", "self",
#'       "WWM", "dam"
#'     ))
#' }
#' @export
define_trait <- function(pop,
                         trait_name,
                         target_add_var   = NULL,
                         target_add_mean  = 0,
                         expressed_parent = c("both", "parent_1", "parent_2"),
                         description      = NULL,
                         units            = NULL,
                         overwrite        = FALSE) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)

  stopifnot(is.character(trait_name), length(trait_name) == 1, nchar(trait_name) > 0)
  validate_trait_name(trait_name)

  expressed_parent <- match.arg(expressed_parent)

  pop <- ensure_trait_tables(pop)

  # Check for existing trait
  existing <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM trait_meta WHERE trait_name = '",
           gsub("'", "''", trait_name), "'")
  )$n

  if (existing > 0 && !overwrite) {
    stop(
      "Trait '", trait_name, "' already exists in trait_meta. ",
      "Use `overwrite = TRUE` to replace.",
      call. = FALSE
    )
  }

  if (existing > 0 && overwrite) {
    DBI::dbExecute(pop$db_conn,
      paste0("DELETE FROM trait_meta WHERE trait_name = '",
             gsub("'", "''", trait_name), "'"))
    DBI::dbExecute(pop$db_conn,
      paste0("DELETE FROM trait_effects WHERE phenotype_name = '",
             gsub("'", "''", trait_name), "'"))
  }

  if (!is.null(target_add_var)) {
    if (!is.numeric(target_add_var) || length(target_add_var) != 1 ||
        is.na(target_add_var) || target_add_var < 0) {
      stop("`target_add_var` must be a non-negative number.", call. = FALSE)
    }
    pop <- write_effect_cov_diagonal(pop, "gen_add", trait_name,
                                     as.numeric(target_add_var))
  }

  row <- tibble::tibble(
    id_trait         = next_int_id(pop$db_conn, "trait_meta", "id_trait"),
    trait_name       = trait_name,
    description      = if (is.null(description)) NA_character_
                       else as.character(description),
    units            = if (is.null(units)) NA_character_
                       else as.character(units),
    expressed_parent = expressed_parent,
    target_add_mean  = as.numeric(target_add_mean)
  )

  DBI::dbWriteTable(pop$db_conn, "trait_meta", row, append = TRUE)

  # Write a global index_meta entry for this trait (index_name = NULL)
  # Allows economic_weight to be assigned later via define_index()
  tn_safe <- gsub("'", "''", trait_name)
  existing_ev <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM index_meta ",
           "WHERE index_name IS NULL AND trait_name = '", tn_safe, "'"))$n

  if (existing_ev == 0) {
    new_id <- next_int_id(pop$db_conn, "index_meta", "id_index_name")
    DBI::dbExecute(pop$db_conn,
      paste0("INSERT INTO index_meta ",
             "(id_index_name, index_name, trait_name, index_weight, economic_weight) ",
             "VALUES (", new_id, ", NULL, '", tn_safe, "', NULL, 0)"))
  } else if (overwrite) {
    DBI::dbExecute(pop$db_conn,
      paste0("UPDATE index_meta SET economic_weight = 0 ",
             "WHERE index_name IS NULL AND trait_name = '", tn_safe, "'"))
  }

  message("Added trait '", trait_name, "'.")

  invisible(pop)
}


#' Ensure the trait-layer tables exist in the database
#'
#' @description
#' Creates all core trait / phenotype / individual tables if they do not already
#' exist, and runs any necessary schema migrations on older databases. Idempotent.
#'
#' @param pop A `tidybreed_pop` object.
#' @return The `tidybreed_pop` object with `$tables` updated.
#' @keywords internal
ensure_trait_tables <- function(pop) {

  existing <- DBI::dbListTables(pop$db_conn)
  con      <- pop$db_conn

  # ── New DDL (only applied if table does not exist) ─────────────────────────

  ddl <- list(

    trait_meta = "
      CREATE TABLE trait_meta (
        id_trait         INTEGER PRIMARY KEY,
        trait_name       VARCHAR UNIQUE NOT NULL,
        description      VARCHAR,
        units            VARCHAR,
        expressed_parent VARCHAR DEFAULT 'both',
        target_add_mean  DOUBLE DEFAULT 0
      )
    ",

    trait_effects = "
      CREATE TABLE trait_effects (
        phenotype_name    VARCHAR NOT NULL,
        effect_name       VARCHAR NOT NULL,
        effect_class      VARCHAR,
        source_column     VARCHAR,
        source_table      VARCHAR,
        distribution      VARCHAR,
        levels_json       VARCHAR,
        slope             DOUBLE,
        center            DOUBLE,
        value             DOUBLE,
        poly_order        INTEGER DEFAULT 1,
        null_class_action VARCHAR DEFAULT 'skip',
        PRIMARY KEY (phenotype_name, effect_name)
      )
    ",

    trait_random_effects = "
      CREATE TABLE trait_random_effects (
        phenotype_name VARCHAR NOT NULL,
        effect_name    VARCHAR NOT NULL,
        level          VARCHAR NOT NULL,
        draw_value     DOUBLE,
        date_sampled   DATE,
        PRIMARY KEY (phenotype_name, effect_name, level)
      )
    ",

    ind_phenotype = "
      CREATE TABLE ind_phenotype (
        id_phenotype   INTEGER PRIMARY KEY,
        id_ind         VARCHAR,
        phenotype_name VARCHAR,
        pheno_value    DOUBLE,
        pheno_number   INTEGER
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
        id_index_name   INTEGER PRIMARY KEY,
        index_name      VARCHAR,
        trait_name      VARCHAR,
        index_weight    DOUBLE,
        economic_weight DOUBLE,
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
    ",

    ind_true_index = "
      CREATE TABLE ind_true_index (
        id_true_index    INTEGER PRIMARY KEY,
        id_ind           VARCHAR NOT NULL,
        index_name       VARCHAR NOT NULL,
        weight_type      VARCHAR NOT NULL,
        true_index_value DOUBLE
      )
    ",

    # ── New tables introduced in v0.31.0 ──────────────────────────────────────

    phenotype_meta = "
      CREATE TABLE phenotype_meta (
        id_phenotype_meta INTEGER PRIMARY KEY,
        phenotype_name    VARCHAR UNIQUE NOT NULL,
        trait_type        VARCHAR,
        mean              DOUBLE DEFAULT 0,
        expressed_sex     VARCHAR DEFAULT 'both',
        repeatable        BOOLEAN DEFAULT FALSE,
        min_value         DOUBLE,
        max_value         DOUBLE,
        prevalence        DOUBLE,
        thresholds        VARCHAR,
        cat_values        VARCHAR,
        cat_names         VARCHAR,
        store_liability   BOOLEAN DEFAULT FALSE
      )
    ",

    phenotype_components = "
      CREATE TABLE phenotype_components (
        id_phenotype_comp   INTEGER PRIMARY KEY,
        phenotype_name      VARCHAR NOT NULL,
        source_trait_name   VARCHAR NOT NULL,
        contributor_type    VARCHAR NOT NULL,
        group_column        VARCHAR,
        group_table         VARCHAR DEFAULT 'ind_meta',
        aggregation         VARCHAR DEFAULT 'sum',
        weight              DOUBLE  DEFAULT 1.0,
        weight_type         VARCHAR DEFAULT 'fixed',
        covariate_name      VARCHAR,
        covariate_table     VARCHAR,
        poly_order          INTEGER,
        poly_scale_min      DOUBLE,
        poly_scale_max      DOUBLE,
        genome_effect_types VARCHAR DEFAULT 'additive',
        missing_action      VARCHAR DEFAULT 'skip',
        contributor_filter  VARCHAR
      )
    ",

    phenotype_residual_cov = "
      CREATE TABLE phenotype_residual_cov (
        id_residual_cov  INTEGER PRIMARY KEY,
        phenotype_name_1 VARCHAR NOT NULL,
        phenotype_name_2 VARCHAR NOT NULL,
        cov_value        DOUBLE  NOT NULL,
        condition_column VARCHAR,
        condition_table  VARCHAR DEFAULT 'ind_meta',
        condition_level  VARCHAR,
        weight_type      VARCHAR DEFAULT 'fixed',
        poly_order       INTEGER
      )
    "
  )

  for (tbl in names(ddl)) {
    if (!tbl %in% existing) {
      DBI::dbExecute(con, ddl[[tbl]])
    }
  }

  # ── Migrations for databases created before v0.31.0 ──────────────────────

  # ind_tbv: drop date_calc if present
  if ("ind_tbv" %in% existing) {
    tbv_cols <- DBI::dbListFields(con, "ind_tbv")
    if ("date_calc" %in% tbv_cols)
      DBI::dbExecute(con, "ALTER TABLE ind_tbv DROP COLUMN date_calc")
  }

  # ind_ebv: drop date_calc, add eval_number if needed
  if ("ind_ebv" %in% existing) {
    ebv_cols <- DBI::dbListFields(con, "ind_ebv")
    if ("date_calc" %in% ebv_cols)
      DBI::dbExecute(con, "ALTER TABLE ind_ebv DROP COLUMN date_calc")
    if (!"eval_number" %in% ebv_cols)
      DBI::dbExecute(con,
        "ALTER TABLE ind_ebv ADD COLUMN eval_number INTEGER DEFAULT 1")
  }

  # index_meta: add economic_weight if missing
  if ("index_meta" %in% existing) {
    idx_cols <- DBI::dbListFields(con, "index_meta")
    if (!"economic_weight" %in% idx_cols)
      DBI::dbExecute(con,
        "ALTER TABLE index_meta ADD COLUMN economic_weight DOUBLE")
  }

  # trait_effects: rename trait_name → phenotype_name; add poly_order + null_class_action
  if ("trait_effects" %in% existing) {
    te_cols <- DBI::dbListFields(con, "trait_effects")
    if ("trait_name" %in% te_cols && !"phenotype_name" %in% te_cols)
      DBI::dbExecute(con,
        "ALTER TABLE trait_effects RENAME COLUMN trait_name TO phenotype_name")
    te_cols <- DBI::dbListFields(con, "trait_effects")   # re-read after rename
    for (new_col in c("source_table VARCHAR", "slope DOUBLE", "center DOUBLE",
                      "poly_order INTEGER", "null_class_action VARCHAR")) {
      col_nm <- strsplit(new_col, " ")[[1]][1]
      if (!col_nm %in% te_cols)
        DBI::dbExecute(con,
          paste0("ALTER TABLE trait_effects ADD COLUMN ", new_col))
    }
  }

  # trait_random_effects: rename trait_name → phenotype_name
  if ("trait_random_effects" %in% existing) {
    tre_cols <- DBI::dbListFields(con, "trait_random_effects")
    if ("trait_name" %in% tre_cols && !"phenotype_name" %in% tre_cols)
      DBI::dbExecute(con,
        "ALTER TABLE trait_random_effects RENAME COLUMN trait_name TO phenotype_name")
  }

  # ind_phenotype: rename trait_name → phenotype_name
  if ("ind_phenotype" %in% existing) {
    ip_cols <- DBI::dbListFields(con, "ind_phenotype")
    if ("trait_name" %in% ip_cols && !"phenotype_name" %in% ip_cols)
      DBI::dbExecute(con,
        "ALTER TABLE ind_phenotype RENAME COLUMN trait_name TO phenotype_name")
  }

  pop$tables <- unique(c(pop$tables, names(ddl)))
  pop
}


#' Validate a trait name for use as SQL identifier
#'
#' @param name Character trait name.
#' @keywords internal
validate_trait_name <- function(name) {
  validate_sql_identifier(name, what = "trait name")
}

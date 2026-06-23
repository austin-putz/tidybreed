#' Define an observed phenotype
#'
#' @description
#' Registers an **observed phenotype**: the quantity written to `ind_phenotype`
#' when [add_phenotype()] is called. For simple (non-composite) traits,
#' `phenotype_name` equals `trait_name` and [define_trait()] must already have
#' been called for that name. For composite traits (e.g. weaning weight
#' assembled from direct and maternal components), no [define_trait()] call is
#' needed for the composite name itself.
#'
#' Writes one row to `phenotype_meta`. If `residual_var` is supplied, also
#' writes an unconditional diagonal entry to `phenotype_residual_cov`. If
#' `components` is supplied, writes one row per component to
#' `phenotype_components`.
#'
#' @param pop A `tidybreed_pop` object.
#' @param phenotype_name Character. Name of the observed phenotype; equals
#'   `trait_name` for simple (non-composite) traits.
#' @param type Character. One of `"continuous"`, `"count"`, `"categorical"`,
#'   or `"derived_formula"`. `"derived_formula"` phenotypes are computed at
#'   [add_phenotype()] time by evaluating the `formula` expression over already-
#'   recorded phenotype values for the same individuals; they have no TBV,
#'   no residual variance, and no QTL of their own.
#' @param mean Numeric. Phenotypic population mean (intercept). Default `0`.
#' @param expressed_sex Character. Who receives a phenotype record: `"both"`
#'   (default), `"M"`, or `"F"`.
#' @param repeatable Logical. Whether an individual can have multiple records
#'   for this phenotype (e.g. repeated test-day or litter-size records).
#'   Default `FALSE`.
#' @param min_value,max_value Numeric. Clipping bounds for count traits.
#'   `NULL` means no limit.
#' @param prevalence Numeric between 0 and 1. For categorical traits with one
#'   threshold (two categories), the fraction expected above the threshold.
#'   Mutually exclusive with `thresholds`.
#' @param thresholds Numeric vector of length K−1 for K ordered categories.
#'   Liability cutpoints in ascending order. Mutually exclusive with
#'   `prevalence`.
#' @param cat_values Numeric vector of length K. Phenotype value stored in
#'   `ind_phenotype` for each category. Defaults to `1, 2, ..., K`.
#' @param cat_names Character vector of length K. Human-readable label per
#'   category, e.g. `c("Alive", "Dead")`. Must not contain commas.
#' @param store_liability Logical. When `TRUE`, the underlying liability value
#'   is written to the reserved `liability_value` column in `ind_phenotype`.
#'   Only meaningful for categorical traits.
#' @param residual_var Numeric. Scalar residual variance. When supplied, writes
#'   one unconditional row to `phenotype_residual_cov`
#'   (`condition_column = NULL`). For heterogeneous residuals or
#'   multi-phenotype correlated residuals, use [define_residual_cov()]
#'   afterwards.
#' @param components A data frame or `tibble` with one row per genetic
#'   component. Columns:
#'   - `source_trait_name` (required): component trait name in `trait_meta`.
#'   - `contributor_type` (required): `"self"`, `"dam"`, `"sire"`, or
#'     `"group"`.
#'   - `weight` (optional, default `1.0`): scalar multiplier.
#'   - `weight_type` (optional, default `"fixed"`): `"fixed"`, `"covariate"`,
#'     `"legendre"`, or `"raw_poly"`.
#'   - `covariate_name` (optional): covariate key.
#'   - `covariate_table` (optional): table containing the covariate column;
#'     `NULL` means value supplied at [add_phenotype()] call time.
#'   - `poly_order` (optional): polynomial basis order.
#'   - `poly_scale_min`, `poly_scale_max` (optional): Legendre scaling bounds.
#'   - `genome_effect_types` (optional, default `"additive"`).
#'   - `group_column` (optional): column defining group membership.
#'   - `group_table` (optional, default `"ind_meta"`): table containing
#'     `group_column`.
#'   - `aggregation` (optional, default `"sum"`): `"sum"` or `"mean"` for
#'     group contributors.
#'
#'   `NULL` (default) → simple single-self trait; `phenotype_components` not
#'   written. Mutually exclusive with `formula_tbv`.
#' @param formula_tbv Character. DSL shorthand for assembling a composite TBV
#'   from component traits already in `trait_meta`. Supports the operators
#'   `+`, `-`, `*` (scalar multiply), `group_sum(trait, col)`, and
#'   `group_mean(trait, col)`. Contributor roles are encoded by suffix:
#'   `_self`, `_dam`, `_sire` (e.g. `"WWD_self + WWM_dam"`). Mutually
#'   exclusive with `components`. Not valid with `type = "derived_formula"`.
#' @param formula Character. Arithmetic expression evaluated over already-
#'   recorded phenotype values to produce a derived phenotype. Phenotype names
#'   reference `ind_phenotype` records produced earlier in the same
#'   [add_phenotype()] call. Operators `+`, `-`, `*`, `/` and parentheses are
#'   supported. Required when `type = "derived_formula"`. Not valid otherwise.
#' @param missing_component_action Character. What to do when an individual is
#'   missing one or more required components (e.g. no group assignment, missing
#'   dam TBV, missing random effect draw). `"skip"` (default) excludes the
#'   individual from `ind_phenotype` and emits a warning with a count.
#'   `"error"` stops with an informative message listing affected individuals.
#'   Stored in `phenotype_meta` so the behaviour is consistent across all
#'   [add_phenotype()] calls for this phenotype.
#' @param overwrite Logical. If `TRUE` and a phenotype with the same name
#'   already exists, replace its rows in `phenotype_meta` and
#'   `phenotype_components`. Default `FALSE` errors on duplicate.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_trait()], [define_residual_cov()], [add_phenotype()],
#'   [define_trait_simple()]
#'
#' @examples
#' \dontrun{
#' # ── Simple continuous trait ──────────────────────────────────────────────
#' pop <- pop |>
#'   define_trait("ADG", target_add_var = 100) |>
#'   get_table("genome_meta") |>
#'   define_additive_effects("ADG") |>
#'   define_phenotype("ADG",
#'     type         = "continuous",
#'     mean         = 850,
#'     residual_var = 120)
#'
#' # ── Count trait with clipping bounds ────────────────────────────────────
#' pop <- pop |>
#'   define_trait("NW", target_add_var = 2) |>
#'   get_table("genome_meta") |>
#'   define_additive_effects("NW") |>
#'   define_phenotype("NW",
#'     type         = "count",
#'     mean         = 10,
#'     min_value    = 1,
#'     max_value    = 30,
#'     residual_var = 8)
#'
#' # ── Categorical trait (binary via prevalence) ────────────────────────────
#' pop <- pop |>
#'   define_trait("mort", target_add_var = 0.05) |>
#'   get_table("genome_meta") |>
#'   define_additive_effects("mort") |>
#'   define_phenotype("mort",
#'     type       = "categorical",
#'     prevalence = 0.05,
#'     cat_names  = c("Alive", "Dead"))
#'
#' # ── Maternal composite via components data frame: WW = WWD (self) + WWM (dam) ──
#' pop <- pop |>
#'   define_phenotype("WW",
#'     type         = "continuous",
#'     mean         = 230,
#'     residual_var = 180,
#'     components   = tibble::tribble(
#'       ~source_trait_name, ~contributor_type,
#'       "WWD",              "self",
#'       "WWM",              "dam"
#'     ))
#'
#' # ── Maternal composite via formula_tbv shorthand (equivalent to above) ──
#' pop <- pop |>
#'   define_phenotype("WW",
#'     type         = "continuous",
#'     mean         = 230,
#'     residual_var = 180,
#'     formula_tbv  = "WWD_self + WWM_dam")
#'
#' # ── SGE (social genetic effects): ADG = direct_self + social group sum ──
#' pop <- pop |>
#'   define_phenotype("ADG_sge",
#'     type         = "continuous",
#'     mean         = 850,
#'     residual_var = 100,
#'     formula_tbv  = "ADG_direct_self + group_sum(ADG_social, pen_id)")
#'
#' # ── Derived formula: FCR computed from already-recorded ADFI and ADG ────
#' # (Define ADFI and ADG first, then derive FCR — no TBV or residual needed)
#' pop <- pop |>
#'   define_phenotype("FCR",
#'     type    = "derived_formula",
#'     formula = "ADFI / ADG")
#' }
#' @export
define_phenotype <- function(pop,
                             phenotype_name,
                             type                     = c("continuous", "count",
                                                          "categorical",
                                                          "derived_formula"),
                             mean                     = 0,
                             expressed_sex            = c("both", "M", "F"),
                             repeatable               = FALSE,
                             min_value                = NULL,
                             max_value                = NULL,
                             prevalence               = NULL,
                             thresholds               = NULL,
                             cat_values               = NULL,
                             cat_names                = NULL,
                             store_liability          = FALSE,
                             residual_var             = NULL,
                             components               = NULL,
                             formula_tbv              = NULL,
                             formula                  = NULL,
                             missing_component_action = c("skip", "error"),
                             overwrite                = FALSE) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)

  stopifnot(is.character(phenotype_name), length(phenotype_name) == 1,
            nchar(phenotype_name) > 0)
  validate_sql_identifier(phenotype_name, what = "phenotype name")

  type                     <- match.arg(type)
  expressed_sex            <- match.arg(expressed_sex)
  missing_component_action <- match.arg(missing_component_action)

  # ── Categorical validation ─────────────────────────────────────────────────

  if (type == "categorical") {
    has_thresholds <- !is.null(thresholds) && length(thresholds) >= 1 &&
                      !all(is.na(thresholds))
    has_prevalence <- !is.null(prevalence) && !is.na(prevalence) &&
                      prevalence > 0 && prevalence < 1

    if (!has_thresholds && !has_prevalence) {
      stop(
        "Categorical phenotypes require `thresholds` (numeric cutpoints vector) ",
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
        stop("`cat_values` must be a numeric vector of length ", n_cats, ".",
             call. = FALSE)
      }
    }
    if (!is.null(cat_names)) {
      ref_len <- if (!is.null(cat_values)) length(cat_values) else n_cats
      if (!is.character(cat_names) || length(cat_names) != ref_len) {
        stop("`cat_names` must be a character vector of length ", ref_len, ".",
             call. = FALSE)
      }
      if (any(grepl(",", cat_names, fixed = TRUE))) {
        stop("`cat_names` entries must not contain commas.", call. = FALSE)
      }
    }
  }

  # ── Formula / derived_formula validation ──────────────────────────────────

  # formula_tbv and components are mutually exclusive
  if (!is.null(formula_tbv) && !is.null(components))
    stop(
      "Supply `formula_tbv` OR `components`, not both. ",
      "`formula_tbv` is the DSL shorthand; `components` is the advanced path ",
      "(covariate weights, Legendre polynomial weights).",
      call. = FALSE
    )

  # formula requires type = "derived_formula"
  if (!is.null(formula) && type != "derived_formula")
    stop(
      "`formula` requires `type = \"derived_formula\"`. ",
      "For composite TBV composition use `formula_tbv` instead.",
      call. = FALSE
    )

  # derived_formula guards
  if (type == "derived_formula") {
    if (!is.null(residual_var))
      stop(
        "`derived_formula` phenotypes have no independent residual variance. ",
        "Do not supply `residual_var` — residual properties are determined by ",
        "the component phenotypes.",
        call. = FALSE
      )
    if (!is.null(components))
      stop(
        "`derived_formula` phenotypes have no TBV components. Do not supply `components`.",
        call. = FALSE
      )
    if (!is.null(formula_tbv))
      stop(
        "`derived_formula` phenotypes use `formula` (not `formula_tbv`). ",
        "`formula_tbv` is for composite TBV assembly, not arithmetic derivation.",
        call. = FALSE
      )
    if (is.null(formula))
      stop(
        "`type = \"derived_formula\"` requires a `formula` string. ",
        "Example: `formula = \"ADFI / ADG\"`.",
        call. = FALSE
      )
  }

  # Validate formula_tbv: parse + trait symbol check (group columns deferred)
  if (!is.null(formula_tbv)) {
    if (!is.character(formula_tbv) || length(formula_tbv) != 1 || !nzchar(formula_tbv))
      stop("`formula_tbv` must be a non-empty character string.", call. = FALSE)
    known_traits <- DBI::dbGetQuery(
      pop$db_conn, "SELECT trait_name FROM trait_meta"
    )$trait_name
    .validate_formula_tbv(formula_tbv, known_traits)
    if (grepl("group_sum|group_mean", formula_tbv))
      message(
        "Note: group column name(s) in formula_tbv are validated at add_phenotype() ",
        "time, not now. A clear error is raised then if a column is missing."
      )
  }

  # Validate derived formula symbols (best-effort — warn if not yet defined)
  if (!is.null(formula)) {
    if (!is.character(formula) || length(formula) != 1 || !nzchar(formula))
      stop("`formula` must be a non-empty character string.", call. = FALSE)
    known_phenos <- DBI::dbGetQuery(
      pop$db_conn, "SELECT phenotype_name FROM phenotype_meta"
    )$phenotype_name
    .validate_derived_formula(formula, known_phenos)
  }

  # ── Ensure tables exist ────────────────────────────────────────────────────

  pop <- ensure_trait_tables(pop)

  # ── Overwrite guard ────────────────────────────────────────────────────────

  pn_safe <- gsub("'", "''", phenotype_name)
  existing_n <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT COUNT(*) AS n FROM phenotype_meta WHERE phenotype_name = '",
           pn_safe, "'"))$n

  if (existing_n > 0 && !overwrite) {
    stop(
      "Phenotype '", phenotype_name, "' already exists in phenotype_meta. ",
      "Use `overwrite = TRUE` to replace.",
      call. = FALSE
    )
  }

  if (existing_n > 0 && overwrite) {
    DBI::dbExecute(pop$db_conn,
      paste0("DELETE FROM phenotype_meta WHERE phenotype_name = '", pn_safe, "'"))
    DBI::dbExecute(pop$db_conn,
      paste0("DELETE FROM phenotype_components WHERE phenotype_name = '", pn_safe, "'"))
    DBI::dbExecute(pop$db_conn,
      paste0("DELETE FROM phenotype_var_comp ",
             "WHERE effect_name = 'residual' AND phenotype_name_1 = '", pn_safe, "'",
             " AND condition_column IS NULL"))
  }

  # ── Serialize categorical fields ──────────────────────────────────────────

  thresholds_str <- if (is.null(thresholds) || all(is.na(thresholds))) {
    NA_character_
  } else {
    paste(thresholds, collapse = ",")
  }

  cat_values_str <- if (is.null(cat_values)) NA_character_
                    else paste(cat_values, collapse = ",")

  cat_names_str  <- if (is.null(cat_names)) NA_character_
                    else paste(cat_names, collapse = ",")

  # ── Insert into phenotype_meta ─────────────────────────────────────────────

  if (is.null(mean) || length(mean) == 0) {
    stop(
      "`mean` for phenotype '", phenotype_name, "' is NULL or empty. ",
      "Check that the config value is defined (e.g. config$general$wean_weight_mean).",
      call. = FALSE
    )
  }

  new_id <- next_int_id(pop$db_conn, "phenotype_meta", "id_phenotype_meta")

  row <- tibble::tibble(
    id_phenotype_meta        = new_id,
    phenotype_name           = phenotype_name,
    type                     = type,
    mean                     = as.numeric(mean),
    expressed_sex            = expressed_sex,
    repeatable               = as.logical(repeatable),
    min_value                = if (is.null(min_value)) NA_real_ else as.numeric(min_value),
    max_value                = if (is.null(max_value)) NA_real_ else as.numeric(max_value),
    prevalence               = if (is.null(prevalence)) NA_real_ else as.numeric(prevalence),
    thresholds               = thresholds_str,
    cat_values               = cat_values_str,
    cat_names                = cat_names_str,
    store_liability          = as.logical(store_liability),
    missing_component_action = missing_component_action,
    formula_tbv              = if (is.null(formula_tbv)) NA_character_ else formula_tbv,
    formula                  = if (is.null(formula))     NA_character_ else formula
  )

  DBI::dbWriteTable(pop$db_conn, "phenotype_meta", row, append = TRUE)

  # ── Residual variance ──────────────────────────────────────────────────────

  if (!is.null(residual_var)) {
    if (!is.numeric(residual_var) || length(residual_var) != 1 ||
        is.na(residual_var) || residual_var < 0) {
      stop("`residual_var` must be a non-negative number.", call. = FALSE)
    }
    pop <- define_residual_cov(
      pop,
      phenotype_names  = phenotype_name,
      cov_matrix       = matrix(as.numeric(residual_var), 1L, 1L,
                                dimnames = list(phenotype_name, phenotype_name)),
      condition_column = NULL
    )
  }

  # ── Components ────────────────────────────────────────────────────────────

  if (!is.null(components)) {
    required_cols <- c("source_trait_name", "contributor_type")
    missing_req <- setdiff(required_cols, names(components))
    if (length(missing_req) > 0) {
      stop("components data frame is missing required column(s): ",
           paste(missing_req, collapse = ", "), call. = FALSE)
    }

    valid_contributor_types <- c("self", "dam", "sire", "group")
    bad_ct <- !components$contributor_type %in% valid_contributor_types
    if (any(bad_ct)) {
      stop("Invalid contributor_type value(s): ",
           paste(unique(components$contributor_type[bad_ct]), collapse = ", "),
           ". Must be one of: ", paste(valid_contributor_types, collapse = ", "),
           call. = FALSE)
    }

    # group contributors require group_column
    grp_rows <- components[components$contributor_type == "group", , drop = FALSE]
    if (nrow(grp_rows) > 0) {
      gcol <- if ("group_column" %in% names(grp_rows)) grp_rows$group_column
              else rep(NA_character_, nrow(grp_rows))
      bad_gc <- is.na(gcol) | !nzchar(as.character(gcol))
      if (any(bad_gc))
        stop("components with contributor_type = 'group' must specify group_column.",
             call. = FALSE)
    }

    # Fill defaults for optional columns
    n_comp <- nrow(components)
    if (!"weight"              %in% names(components)) components$weight              <- 1.0
    if (!"weight_type"         %in% names(components)) components$weight_type         <- "fixed"
    if (!"covariate_name"      %in% names(components)) components$covariate_name      <- NA_character_
    if (!"covariate_table"     %in% names(components)) components$covariate_table     <- NA_character_
    if (!"poly_order"          %in% names(components)) components$poly_order          <- NA_integer_
    if (!"poly_scale_min"      %in% names(components)) components$poly_scale_min      <- NA_real_
    if (!"poly_scale_max"      %in% names(components)) components$poly_scale_max      <- NA_real_
    if (!"genome_effect_types" %in% names(components)) components$genome_effect_types <- "additive"
    if (!"group_column"        %in% names(components)) components$group_column        <- NA_character_
    if (!"group_table"         %in% names(components)) components$group_table         <- "ind_meta"
    if (!"aggregation"         %in% names(components)) components$aggregation         <- "sum"
    if (!"missing_action"      %in% names(components)) components$missing_action      <- "skip"
    if (!"contributor_filter"  %in% names(components)) components$contributor_filter  <- NA_character_

    # Replace NA in group_table with default
    components$group_table[is.na(components$group_table)] <- "ind_meta"
    components$aggregation[is.na(components$aggregation)] <- "sum"

    first_id <- next_int_id(pop$db_conn, "phenotype_components", "id_phenotype_comp")
    comp_rows <- tibble::tibble(
      id_phenotype_comp   = seq(first_id, first_id + n_comp - 1L),
      phenotype_name      = phenotype_name,
      source_trait_name   = as.character(components$source_trait_name),
      contributor_type    = as.character(components$contributor_type),
      group_column        = as.character(components$group_column),
      group_table         = as.character(components$group_table),
      aggregation         = as.character(components$aggregation),
      weight              = as.numeric(components$weight),
      weight_type         = as.character(components$weight_type),
      covariate_name      = as.character(components$covariate_name),
      covariate_table     = as.character(components$covariate_table),
      poly_order          = as.integer(components$poly_order),
      poly_scale_min      = as.numeric(components$poly_scale_min),
      poly_scale_max      = as.numeric(components$poly_scale_max),
      genome_effect_types = as.character(components$genome_effect_types),
      missing_action      = as.character(components$missing_action),
      contributor_filter  = as.character(components$contributor_filter)
    )

    DBI::dbWriteTable(pop$db_conn, "phenotype_components", comp_rows, append = TRUE)
  }

  msg_suffix <- ""
  if (!is.null(formula_tbv)) msg_suffix <- " [formula_tbv DSL]"
  if (!is.null(formula))     msg_suffix <- " [derived_formula]"
  message("Added phenotype '", phenotype_name, "' (type: ", type, ")",
          msg_suffix, ".")

  invisible(pop)
}

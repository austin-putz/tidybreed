#' Schema metadata functions for tidybreed populations
#'
#' @description
#' These functions allow users to view and manage table and column descriptions
#' stored directly in the DuckDB database. Descriptions persist across
#' `restore_pop()` restores and travel with the `.duckdb` file.
#'
#' @keywords internal
#' @name schema
NULL


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Register table/column descriptions in _schema_meta (upsert)
#'
#' @param conn A DBI connection to the DuckDB database.
#' @param entries A data.frame with columns: object_type, table_name,
#'   column_name (NA for table-level), description, notes (NA ok).
#' @keywords internal
register_schema_meta <- function(conn, entries) {
  if (!is.data.frame(entries) || nrow(entries) == 0L) return(invisible(NULL))

  # DELETE existing matching rows (upsert pattern)
  for (i in seq_len(nrow(entries))) {
    row       <- entries[i, ]
    ot_esc    <- gsub("'", "''", as.character(row$object_type))
    tbl_esc   <- gsub("'", "''", as.character(row$table_name))
    col_name  <- row$column_name
    if (is.na(col_name)) {
      del_sql <- paste0(
        "DELETE FROM _schema_meta WHERE object_type = '", ot_esc,
        "' AND table_name = '", tbl_esc, "' AND column_name IS NULL"
      )
    } else {
      col_esc <- gsub("'", "''", as.character(col_name))
      del_sql <- paste0(
        "DELETE FROM _schema_meta WHERE object_type = '", ot_esc,
        "' AND table_name = '", tbl_esc,
        "' AND column_name = '", col_esc, "'"
      )
    }
    DBI::dbExecute(conn, del_sql)
  }

  # INSERT all rows using DBI (handles escaping automatically)
  start_id <- next_int_id(conn, "_schema_meta", "id_schema_meta")
  insert_df <- data.frame(
    id_schema_meta = seq.int(start_id, length.out = nrow(entries)),
    object_type    = as.character(entries$object_type),
    table_name     = as.character(entries$table_name),
    column_name    = entries$column_name,
    description    = as.character(entries$description),
    notes          = entries$notes,
    stringsAsFactors = FALSE
  )
  DBI::dbAppendTable(conn, "_schema_meta", insert_df)
  invisible(NULL)
}


#' Pre-built descriptions for genome table columns (genome_meta, haplotype, genotype)
#' @keywords internal
.genome_table_descriptions <- function() {
  rbind(
    # genome_meta
    .sm_tbl("genome_meta",
            "Locus-level metadata. One row per locus. Stores chromosome assignment, physical position, and chip membership flags added by define_chip()."),
    .sm_col("genome_meta", "locus_id",
            "Sequential integer primary key; 1 to n_loci in creation order"),
    .sm_col("genome_meta", "locus_name",
            "Human-readable locus identifier (e.g. 'Locus_1', 'rs12345')"),
    .sm_col("genome_meta", "chr",
            "Chromosome number (integer index 1 to n_chr)"),
    .sm_col("genome_meta", "chr_name",
            "Chromosome label string; defaults to chr cast to character"),
    .sm_col("genome_meta", "pos_Mb",
            "Physical position in megabases along the chromosome"),
    .sm_col("genome_meta", "founder_allele_freq",
            "Per-locus base allele frequency for founder haplotype sampling; present only when define_founder_haplotypes() has been called"),

    # genome_haplotype
    .sm_tbl("genome_haplotype",
            "Phased haplotypes in wide format. Two rows per individual (parent_origin 1 = paternal, 2 = maternal). Locus alleles stored as UTINYINT columns locus_1, locus_2, etc."),
    .sm_col("genome_haplotype", "id_ind",
            "Individual identifier; FK to ind_meta.id_ind"),
    .sm_col("genome_haplotype", "parent_origin",
            "Haplotype of origin: 1 = paternal (parent_1 side), 2 = maternal (parent_2 side)"),

    # genome_genotype
    .sm_tbl("genome_genotype",
            "Genotypes in 0/1/2 dosage encoding. One row per individual. Computed from genome_haplotype as the sum of both haplotypes at each locus."),
    .sm_col("genome_genotype", "id_ind",
            "Individual identifier; FK to ind_meta.id_ind")
  )
}


#' Pre-built descriptions for core (non-genome) tables created by open_pop()
#' @keywords internal
.core_layer_descriptions <- function() {
  rbind(
    # _schema_meta itself
    .sm_tbl("_schema_meta",
            "System table storing table and column descriptions for all tidybreed database objects."),
    .sm_col("_schema_meta", "id_schema_meta",
            "Auto-incrementing primary key"),
    .sm_col("_schema_meta", "object_type",
            "Entry type: 'table' for table descriptions, 'column' for column descriptions"),
    .sm_col("_schema_meta", "table_name",
            "Name of the table this entry describes"),
    .sm_col("_schema_meta", "column_name",
            "Column name; NULL for table-level entries"),
    .sm_col("_schema_meta", "description",
            "Human-readable description of the table or column"),
    .sm_col("_schema_meta", "notes",
            "Optional supplementary context (NULL if unused)"),

    # ind_meta
    .sm_tbl("ind_meta",
            "Individual-level metadata. One row per individual. Core columns are managed by the system; user-defined columns can be added via mutate_table() or the ... arguments of add_founders() and add_offspring()."),
    .sm_col("ind_meta", "id_ind",
            "Primary key; format '{line_name}_{n}' (e.g. 'A_1', 'Holstein_42')"),
    .sm_col("ind_meta", "id_parent_1",
            "Paternal parent id_ind; NA for founders"),
    .sm_col("ind_meta", "id_parent_2",
            "Maternal parent id_ind; NA for founders"),
    .sm_col("ind_meta", "line_name",
            "Genetic line name (e.g. 'A', 'Holstein')"),
    .sm_col("ind_meta", "sex",
            "Sex of the individual: 'M' for male, 'F' for female"),

    # trait_var_comp
    .sm_tbl("trait_var_comp",
            "Genetic variance component storage. One row per (effect_name, trait_name_1, trait_name_2); both (i,j) and (j,i) stored. Reserved effect_name values: 'gen_add', 'dominance' (future), 'epistasis' (future)."),
    .sm_col("trait_var_comp", "id_trait_var_comp",
            "Auto-incrementing primary key"),
    .sm_col("trait_var_comp", "effect_name",
            "Genetic effect type. Reserved: 'gen_add' (additive genetic G matrix); future: 'dominance', 'epistasis'."),
    .sm_col("trait_var_comp", "trait_name_1",
            "First trait in the (co)variance pair; FK to trait_meta.trait_name"),
    .sm_col("trait_var_comp", "trait_name_2",
            "Second trait; equals trait_name_1 for diagonal (variance) entries"),
    .sm_col("trait_var_comp", "cov_value",
            "Variance (diagonal) or covariance (off-diagonal) value"),

    # genome_effects
    .sm_tbl("genome_effects",
            "QTL effect data in long format. One row per locus x trait x effect type x line. Populated by define_additive_effects(). A locus is a QTL for a trait if it has a matching row here."),
    .sm_col("genome_effects", "id_genome_effect",
            "Auto-incrementing primary key"),
    .sm_col("genome_effects", "locus_name",
            "Locus identifier; FK to genome_meta.locus_name"),
    .sm_col("genome_effects", "line_name",
            "Line-specific effect (NULL = population-wide, shared across all lines)"),
    .sm_col("genome_effects", "trait_name",
            "Genetic component trait; FK to trait_meta.trait_name"),
    .sm_col("genome_effects", "genome_effect_type",
            "Type of genetic effect: 'additive' (current); 'dominance' reserved for future use"),
    .sm_col("genome_effects", "genome_value",
            "Effect size (allele substitution effect) in trait units"),
    .sm_col("genome_effects", "base_allele_freq",
            "Allele frequency used for TBV centering via the Falconer formula"),

    # phenotype_meta
    .sm_tbl("phenotype_meta",
            "Observed phenotype definitions. One row per phenotype. Manages the observation layer: type, mean, sex expression, and distributional parameters. Populated by define_phenotype()."),
    .sm_col("phenotype_meta", "id_phenotype_meta",
            "Auto-incrementing primary key"),
    .sm_col("phenotype_meta", "phenotype_name",
            "Unique phenotype identifier; equals trait_name for simple single-component phenotypes"),
    .sm_col("phenotype_meta", "type",
            "Distribution family: 'continuous', 'count', 'categorical', or 'derived_formula'"),
    .sm_col("phenotype_meta", "mean",
            "Phenotypic population mean (liability intercept for categorical phenotypes); default 0"),
    .sm_col("phenotype_meta", "expressed_sex",
            "Which sex receives records: 'both' (default), 'M', or 'F'"),
    .sm_col("phenotype_meta", "repeatable",
            "TRUE if individuals can receive multiple records for this phenotype; default FALSE"),
    .sm_col("phenotype_meta", "min_value",
            "Lower clipping bound for count-type phenotypes; NULL = no lower bound"),
    .sm_col("phenotype_meta", "max_value",
            "Upper clipping bound for count-type phenotypes; NULL = no upper bound"),
    .sm_col("phenotype_meta", "prevalence",
            "Disease prevalence (proportion) for binary categorical phenotypes"),
    .sm_col("phenotype_meta", "thresholds",
            "Comma-separated liability thresholds for multi-category phenotypes"),
    .sm_col("phenotype_meta", "cat_values",
            "Comma-separated numeric phenotype values assigned to each category level"),
    .sm_col("phenotype_meta", "cat_names",
            "Comma-separated labels for each category level"),
    .sm_col("phenotype_meta", "store_liability",
            "When TRUE, raw liability is written to ind_phenotype.liability_value; default FALSE"),
    .sm_col("phenotype_meta", "missing_component_action",
            "Behavior when a composite component is unresolvable: 'skip' (exclude + warn) or 'error' (stop); default 'skip'"),
    .sm_col("phenotype_meta", "formula_tbv",
            "Optional R formula string for computing TBV from trait columns instead of summing genome effects"),
    .sm_col("phenotype_meta", "formula",
            "Optional R formula string for computing the phenotype directly, bypassing standard component assembly"),

    # phenotype_components
    .sm_tbl("phenotype_components",
            "Component definitions for composite phenotypes. One row per phenotype x component. Enables maternal effects, social genetic effects, and multi-contributor phenotypes. Populated by define_phenotype(..., components = ...)."),
    .sm_col("phenotype_components", "id_phenotype_comp",
            "Auto-incrementing primary key"),
    .sm_col("phenotype_components", "phenotype_name",
            "Composite phenotype this component belongs to; FK to phenotype_meta.phenotype_name"),
    .sm_col("phenotype_components", "source_trait_name",
            "Genetic component trait providing the TBV; FK to trait_meta.trait_name"),
    .sm_col("phenotype_components", "contributor_type",
            "Who contributes this component: 'self', 'dam', 'sire', or 'group'"),
    .sm_col("phenotype_components", "group_column",
            "Column in group_table holding group membership; required for contributor_type = 'group'"),
    .sm_col("phenotype_components", "group_table",
            "Table containing group_column; default 'ind_meta'"),
    .sm_col("phenotype_components", "aggregation",
            "How group members' TBVs are combined: 'sum' (default) or 'mean'"),
    .sm_col("phenotype_components", "weight",
            "Scalar multiplier applied to this component's contribution; default 1.0"),
    .sm_col("phenotype_components", "weight_type",
            "How the weight varies: 'fixed' (default), 'covariate', 'legendre', or 'raw_poly'"),
    .sm_col("phenotype_components", "covariate_name",
            "Column providing the covariate when weight_type = 'covariate'"),
    .sm_col("phenotype_components", "covariate_table",
            "Table containing covariate_name"),
    .sm_col("phenotype_components", "poly_order",
            "Polynomial basis order for legendre or raw_poly weight types"),
    .sm_col("phenotype_components", "poly_scale_min",
            "Lower bound for Legendre polynomial scaling"),
    .sm_col("phenotype_components", "poly_scale_max",
            "Upper bound for Legendre polynomial scaling"),
    .sm_col("phenotype_components", "genome_effect_types",
            "Comma-separated genome_effect_type values to include; default 'additive'"),
    .sm_col("phenotype_components", "missing_action",
            "Per-component fallback; currently unused, governed by phenotype_meta.missing_component_action"),
    .sm_col("phenotype_components", "contributor_filter",
            "Reserved for spatial/neighbourhood lookup; not yet implemented"),

    # phenotype_var_comp
    .sm_tbl("phenotype_var_comp",
            "Phenotype-level variance component storage. One row per (effect_name, phenotype pair, optional condition). Stores residual covariances (effect_name = 'residual') and named random effects (hys, litter, pen, etc.). Populated by define_phenotype(), define_residual_cov(), and define_effect_random()."),
    .sm_col("phenotype_var_comp", "id_phenotype_var_comp",
            "Auto-incrementing primary key"),
    .sm_col("phenotype_var_comp", "effect_name",
            "Effect identifier. Reserved: 'residual' (R matrix). Any other name is a user-defined phenotype-level random effect (e.g. 'hys', 'litter', 'pen')."),
    .sm_col("phenotype_var_comp", "phenotype_name_1",
            "First phenotype in the (co)variance pair; FK to phenotype_meta.phenotype_name"),
    .sm_col("phenotype_var_comp", "phenotype_name_2",
            "Second phenotype; equals phenotype_name_1 for diagonal (variance) entries"),
    .sm_col("phenotype_var_comp", "cov_value",
            "Variance (diagonal) or covariance (off-diagonal) value"),
    .sm_col("phenotype_var_comp", "condition_column",
            "Column in condition_table that partitions residuals; NULL = unconditional. Applies to effect_name = 'residual' only."),
    .sm_col("phenotype_var_comp", "condition_table",
            "Table containing condition_column; default 'ind_meta'"),
    .sm_col("phenotype_var_comp", "condition_level",
            "Value of condition_column this row applies to (e.g. 'M', 'F')"),
    .sm_col("phenotype_var_comp", "weight_type",
            "Weight function type: 'fixed' (default) or 'legendre'"),
    .sm_col("phenotype_var_comp", "poly_order",
            "Polynomial order for legendre weight type")
  )
}


#' Combined genome-layer descriptions (core tables + genome tables)
#'
#' Convenience wrapper used by the deprecated initialize_genome() and
#' migrate_schema_meta(). New code should call .core_layer_descriptions() and
#' .genome_table_descriptions() separately.
#' @keywords internal
.genome_layer_descriptions <- function() {
  rbind(.core_layer_descriptions(), .genome_table_descriptions())
}


#' Pre-built descriptions for trait/phenotype-layer tables
#' @keywords internal
.trait_layer_descriptions <- function() {
  rbind(
    # trait_meta
    .sm_tbl("trait_meta",
            "Genetic component trait definitions. One row per trait. Genetic layer only — no observation-layer metadata. Populated by define_trait()."),
    .sm_col("trait_meta", "id_trait",
            "Auto-incrementing primary key"),
    .sm_col("trait_meta", "trait_name",
            "Unique identifier used across genome_effects, ind_tbv, and trait_var_comp"),
    .sm_col("trait_meta", "description",
            "Free-text description of the biological trait"),
    .sm_col("trait_meta", "units",
            "Measurement units (e.g. 'kg', 'g/day', 'count')"),
    .sm_col("trait_meta", "expressed_parent",
            "Parent-of-origin expression: 'both' (default), 'parent_1' (paternal only), 'parent_2' (maternal only)"),
    .sm_col("trait_meta", "target_add_mean",
            "TBV centering mean for the base population; default 0"),

    # trait_effects
    .sm_tbl("trait_effects",
            "Fixed and random effect configurations for phenotype models. One row per phenotype x effect. Populated by define_effect_fixed_class(), define_effect_fixed_cov(), define_effect_random()."),
    .sm_col("trait_effects", "phenotype_name",
            "Phenotype this effect belongs to; FK to phenotype_meta.phenotype_name"),
    .sm_col("trait_effects", "effect_name",
            "Unique label for this effect within the phenotype (e.g. 'sex', 'gen', 'herd')"),
    .sm_col("trait_effects", "effect_class",
            "Effect type: 'fixed_class', 'fixed_cov', or 'random'"),
    .sm_col("trait_effects", "source_column",
            "Column in source_table used as the grouping or covariate variable"),
    .sm_col("trait_effects", "source_table",
            "Table containing source_column; default 'ind_meta'"),
    .sm_col("trait_effects", "distribution",
            "Sampling distribution for random effects: 'normal', 'gamma', or 'uniform'"),
    .sm_col("trait_effects", "levels_json",
            "JSON map of level to shift value for fixed_class effects (e.g. {\"M\": 30, \"F\": 0})"),
    .sm_col("trait_effects", "slope",
            "Regression coefficient for fixed_cov effects"),
    .sm_col("trait_effects", "center",
            "Centering value subtracted from the covariate before multiplying by slope"),
    .sm_col("trait_effects", "value",
            "Rarely used scalar override"),
    .sm_col("trait_effects", "poly_order",
            "Polynomial order for polynomial covariate effects; default 1 (linear)"),
    .sm_col("trait_effects", "null_class_action",
            "Behavior when grouping column is NULL: 'skip' excludes the individual"),

    # trait_random_effects
    .sm_tbl("trait_random_effects",
            "Sampled random effect levels. One row per phenotype x effect x level. Populated by add_phenotype() on first use; subsequent calls reuse the stored value for consistency."),
    .sm_col("trait_random_effects", "phenotype_name",
            "Phenotype this random effect belongs to; FK to phenotype_meta.phenotype_name"),
    .sm_col("trait_random_effects", "effect_name",
            "Random effect name; FK to trait_effects.effect_name"),
    .sm_col("trait_random_effects", "level",
            "The grouping level (e.g. herd ID, sire ID) as a string"),
    .sm_col("trait_random_effects", "draw_value",
            "The sampled random effect value for this level"),
    .sm_col("trait_random_effects", "date_sampled",
            "Date the value was first sampled"),

    # ind_phenotype
    .sm_tbl("ind_phenotype",
            "Phenotype records in long format. One row per individual x phenotype x record number. Populated by add_phenotype(). User-defined columns can be added via the ... argument of add_phenotype()."),
    .sm_col("ind_phenotype", "id_phenotype",
            "Auto-incrementing primary key"),
    .sm_col("ind_phenotype", "id_ind",
            "Individual identifier; FK to ind_meta.id_ind"),
    .sm_col("ind_phenotype", "phenotype_name",
            "Phenotype name; FK to phenotype_meta.phenotype_name"),
    .sm_col("ind_phenotype", "pheno_value",
            "Observed phenotype value (or category code for categorical phenotypes)"),
    .sm_col("ind_phenotype", "pheno_number",
            "Record number for this individual x phenotype combination; 1 = first record"),

    # ind_tbv
    .sm_tbl("ind_tbv",
            "True breeding values (simulation ground truth). One row per individual x trait. Populated by add_phenotype() and add_tbv(). Values computed from genome effects in genome_effects."),
    .sm_col("ind_tbv", "id_tbv",
            "Auto-incrementing primary key"),
    .sm_col("ind_tbv", "id_ind",
            "Individual identifier; FK to ind_meta.id_ind"),
    .sm_col("ind_tbv", "trait_name",
            "Genetic component trait name; FK to trait_meta.trait_name"),
    .sm_col("ind_tbv", "tbv_value",
            "True breeding value for this individual and trait"),

    # ind_ebv
    .sm_tbl("ind_ebv",
            "Estimated breeding values from external BLUP or GBLUP analyses. One row per individual x trait x model x evaluation number. Populated by add_ebv()."),
    .sm_col("ind_ebv", "id_ebv",
            "Auto-incrementing primary key"),
    .sm_col("ind_ebv", "id_ind",
            "Individual identifier; FK to ind_meta.id_ind"),
    .sm_col("ind_ebv", "trait_name",
            "Trait the EBV was estimated for"),
    .sm_col("ind_ebv", "model",
            "User-supplied label for the estimation model (e.g. 'BLUP_v1', 'ssGBLUP')"),
    .sm_col("ind_ebv", "ebv_value",
            "Estimated breeding value"),
    .sm_col("ind_ebv", "acc",
            "Estimated accuracy of the EBV (optional; NULL if not supplied)"),
    .sm_col("ind_ebv", "se",
            "Standard error of the EBV (optional; NULL if not supplied)"),
    .sm_col("ind_ebv", "eval_number",
            "Auto-incrementing evaluation counter per trait x model; 1 = first evaluation"),

    # index_meta
    .sm_tbl("index_meta",
            "Selection index definitions. One row per index x trait. A special row with index_name = NULL holds the global economic weight per trait written by define_trait(). Named indices hold selection weights from define_index()."),
    .sm_col("index_meta", "id_index_name",
            "Auto-incrementing primary key"),
    .sm_col("index_meta", "index_name",
            "Named selection index (e.g. 'NTI'); NULL for global economic weight rows from define_trait()"),
    .sm_col("index_meta", "trait_name",
            "Trait this index entry applies to; FK to trait_meta.trait_name"),
    .sm_col("index_meta", "index_weight",
            "Selection index weight for this trait; NULL for global economic weight rows"),
    .sm_col("index_meta", "economic_weight",
            "Economic value per unit of this trait"),

    # ind_index
    .sm_tbl("ind_index",
            "Computed selection index values. One row per individual x index x run. Multiple runs are distinguished by index_number. Populated by add_index()."),
    .sm_col("ind_index", "id_index",
            "Auto-incrementing primary key"),
    .sm_col("ind_index", "id_ind",
            "Individual identifier; FK to ind_meta.id_ind"),
    .sm_col("ind_index", "index_name",
            "Named selection index; FK to index_meta.index_name"),
    .sm_col("ind_index", "index_number",
            "Auto-incrementing run counter per individual; 1 = first computation"),
    .sm_col("ind_index", "index_value",
            "Computed selection index value (weighted sum of EBVs or phenotypes)"),

    # ind_true_index
    .sm_tbl("ind_true_index",
            "True selection index values computed from TBVs. One row per individual x index x weight type. Populated by add_tbv() when index_names is supplied."),
    .sm_col("ind_true_index", "id_true_index",
            "Auto-incrementing primary key"),
    .sm_col("ind_true_index", "id_ind",
            "Individual identifier; FK to ind_meta.id_ind"),
    .sm_col("ind_true_index", "index_name",
            "Named selection index; FK to index_meta.index_name"),
    .sm_col("ind_true_index", "weight_type",
            "Which weights were used: 'index' (uses index_weight) or 'economic' (uses economic_weight)"),
    .sm_col("ind_true_index", "true_index_value",
            "True index value: weighted sum of TBVs across all index traits")
  )
}


# ── Tiny row-builder helpers ──────────────────────────────────────────────────

.sm_tbl <- function(table_name, description, notes = NA_character_) {
  data.frame(
    object_type = "table",
    table_name  = table_name,
    column_name = NA_character_,
    description = description,
    notes       = notes,
    stringsAsFactors = FALSE
  )
}

.sm_col <- function(table_name, column_name, description, notes = NA_character_) {
  data.frame(
    object_type = "column",
    table_name  = table_name,
    column_name = column_name,
    description = description,
    notes       = notes,
    stringsAsFactors = FALSE
  )
}


# ── Exported functions ────────────────────────────────────────────────────────

#' View table-level descriptions for a tidybreed population
#'
#' @description
#' Returns a tibble of all user-visible tables with row counts, column counts,
#' and descriptions from `_schema_meta`. Print the result for an aligned
#' overview; use `describe_table(pop, "name")` to drill into a specific table.
#'
#' @param pop A `tidybreed_pop` object.
#'
#' @return A tibble of class `tidybreed_schema` with columns
#'   `table_name`, `n_rows`, `n_cols`, and `description`.
#'   Printed via [print.tidybreed_schema()].
#'
#' @seealso [describe_table()], [define_schema_description()]
#'
#' @examples
#' \dontrun{
#' pop <- initialize_genome("MySim", n_loci = 500, n_chr = 5, chr_len_Mb = 100)
#' schema(pop)
#' }
#' @export
schema <- function(pop) {
  if (inherits(pop, "tidybreed_table")) pop <- pop$pop
  stopifnot(inherits(pop, "tidybreed_pop"))

  user_tables <- setdiff(pop$tables, "_schema_meta")

  has_meta <- "_schema_meta" %in% pop$tables

  # Pull descriptions from _schema_meta when available
  if (has_meta) {
    desc_df <- DBI::dbGetQuery(
      pop$db_conn,
      "SELECT table_name, description FROM _schema_meta WHERE object_type = 'table'"
    )
  } else {
    desc_df <- data.frame(table_name = character(), description = character(),
                          stringsAsFactors = FALSE)
  }

  # Build result row by row
  rows <- lapply(user_tables, function(tbl) {
    n_rows <- tryCatch(
      as.integer(DBI::dbGetQuery(pop$db_conn, paste0("SELECT COUNT(*) AS n FROM ", tbl))$n),
      error = function(e) NA_integer_
    )
    n_cols <- tryCatch(
      length(DBI::dbListFields(pop$db_conn, tbl)),
      error = function(e) NA_integer_
    )
    idx   <- match(tbl, desc_df$table_name)
    desc  <- if (!is.na(idx)) desc_df$description[idx] else ""
    data.frame(table_name  = tbl,
               n_rows      = n_rows,
               n_cols      = n_cols,
               description = desc,
               stringsAsFactors = FALSE)
  })

  result <- do.call(rbind, rows)
  if (is.null(result)) {
    result <- data.frame(table_name  = character(),
                         n_rows      = integer(),
                         n_cols      = integer(),
                         description = character(),
                         stringsAsFactors = FALSE)
  }

  structure(
    tibble::as_tibble(result),
    pop_name = pop$pop_name,
    class    = c("tidybreed_schema", "tbl_df", "tbl", "data.frame")
  )
}


#' Print method for tidybreed_schema
#'
#' @param x A `tidybreed_schema` object from [schema()].
#' @param ... Additional arguments (not used).
#' @return `x`, invisibly.
#' @export
print.tidybreed_schema <- function(x, ...) {
  width    <- getOption("width", 80L)
  pop_name <- attr(x, "pop_name")

  left   <- paste0("── Schema: ", pop_name, " ")
  fill_w <- max(2L, width - nchar(left))
  cat(left, strrep("─", fill_w), "\n", sep = "")
  cat("  Use describe_table(pop, \"name\") for column-level details.\n\n")

  if (nrow(x) == 0L) {
    cat("  (no tables)\n")
    return(invisible(x))
  }

  tbl_w  <- max(nchar("Table"),   max(nchar(x$table_name)))
  rows_w <- max(nchar("Rows"),    nchar(format(max(x$n_rows, na.rm = TRUE), big.mark = ",")))
  cols_w <- max(nchar("Cols"),    nchar(format(max(x$n_cols, na.rm = TRUE), big.mark = ",")))

  hdr <- paste0(
    "  ", formatC("Table",       width = -tbl_w),  "  ",
         formatC("Rows",         width =  rows_w),  "  ",
         formatC("Cols",         width =  cols_w),  "  ",
         "Description"
  )
  cat(hdr, "\n")
  sep_w <- max(2L, width - 2L)
  cat("  ", strrep("─", sep_w), "\n", sep = "")

  for (i in seq_len(nrow(x))) {
    row_rows <- if (is.na(x$n_rows[i])) "?" else format(x$n_rows[i], big.mark = ",")
    row_cols <- if (is.na(x$n_cols[i])) "?" else format(x$n_cols[i], big.mark = ",")
    desc     <- if (is.na(x$description[i]) || x$description[i] == "") "(no description)" else x$description[i]

    avail_w  <- max(10L, width - tbl_w - rows_w - cols_w - 8L)
    desc_trunc <- if (nchar(desc) > avail_w) paste0(substr(desc, 1, avail_w - 3L), "...") else desc

    cat("  ",
        formatC(x$table_name[i], width = -tbl_w),  "  ",
        formatC(row_rows,         width =  rows_w),  "  ",
        formatC(row_cols,         width =  cols_w),  "  ",
        desc_trunc, "\n", sep = "")
  }
  invisible(x)
}


#' View column-level descriptions for a tidybreed table
#'
#' @description
#' Returns a tibble of all columns in `table_name` with their DuckDB types
#' and descriptions from `_schema_meta`. Wide genome tables
#' (`genome_haplotype`, `genome_genotype`, `founder_haplotypes`) only list
#' their non-locus metadata columns to avoid printing thousands of locus
#' columns.
#'
#' @param pop A `tidybreed_pop` object.
#' @param table_name Character. Name of the table to describe.
#'
#' @return A tibble of class `tidybreed_table_desc` (returned invisibly) with
#'   columns `column_name`, `column_type`, `description`, and `notes`.
#'   Printed via [print.tidybreed_table_desc()].
#'
#' @seealso [schema()], [define_schema_description()]
#'
#' @examples
#' \dontrun{
#' pop <- initialize_genome("MySim", n_loci = 500, n_chr = 5, chr_len_Mb = 100)
#' describe_table(pop, "ind_meta")
#' describe_table(pop, "genome_effects")
#' }
#' @export
describe_table <- function(pop, table_name) {
  if (inherits(pop, "tidybreed_table")) pop <- pop$pop
  stopifnot(inherits(pop, "tidybreed_pop"))

  if (!table_name %in% pop$tables) {
    stop(
      "Table '", table_name, "' does not exist in this population. ",
      "Available tables: ", paste(setdiff(pop$tables, "_schema_meta"), collapse = ", "),
      call. = FALSE
    )
  }

  # Get column types from DuckDB DESCRIBE
  desc_rows <- tryCatch(
    DBI::dbGetQuery(pop$db_conn, paste0("DESCRIBE ", table_name)),
    error = function(e) data.frame(column_name = character(), column_type = character(),
                                   stringsAsFactors = FALSE)
  )

  # Detect wide genome tables and restrict to non-locus columns for display
  all_cols   <- desc_rows$column_name
  locus_cols <- grep("^locus_[0-9]+$", all_cols, value = TRUE)
  is_wide    <- length(locus_cols) > 0
  if (is_wide) {
    desc_rows <- desc_rows[!desc_rows$column_name %in% locus_cols, , drop = FALSE]
  }

  # Pull column descriptions from _schema_meta
  has_meta <- "_schema_meta" %in% pop$tables
  if (has_meta) {
    tbl_esc  <- gsub("'", "''", table_name)
    meta_rows <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT column_name, description, notes FROM _schema_meta ",
             "WHERE object_type = 'column' AND table_name = '", tbl_esc, "'")
    )
    tbl_desc_row <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT description FROM _schema_meta ",
             "WHERE object_type = 'table' AND table_name = '", tbl_esc, "' ",
             "AND column_name IS NULL")
    )
    tbl_desc <- if (nrow(tbl_desc_row) > 0) tbl_desc_row$description[1] else ""
  } else {
    meta_rows <- data.frame(column_name = character(), description = character(),
                            notes = character(), stringsAsFactors = FALSE)
    tbl_desc  <- ""
  }

  # Join: all actual columns left-joined with metadata
  result_rows <- lapply(seq_len(nrow(desc_rows)), function(i) {
    col   <- desc_rows$column_name[i]
    ctype <- desc_rows$column_type[i]
    idx   <- match(col, meta_rows$column_name)
    if (!is.na(idx)) {
      cdesc  <- meta_rows$description[idx]
      cnotes <- meta_rows$notes[idx]
    } else {
      cdesc  <- ""
      cnotes <- NA_character_
    }
    data.frame(column_name = col, column_type = ctype,
               description = cdesc, notes = cnotes,
               stringsAsFactors = FALSE)
  })

  result <- if (length(result_rows) > 0) do.call(rbind, result_rows) else {
    data.frame(column_name = character(), column_type = character(),
               description = character(), notes = character(),
               stringsAsFactors = FALSE)
  }

  # Flag stale _schema_meta rows (column described but not in table)
  if (has_meta && nrow(meta_rows) > 0) {
    stale <- setdiff(meta_rows$column_name, all_cols)
    if (length(stale) > 0) {
      stale_rows <- meta_rows[meta_rows$column_name %in% stale, , drop = FALSE]
      stale_rows$column_type  <- NA_character_
      stale_rows$notes        <- paste0("(column not found in table",
                                        ifelse(!is.na(stale_rows$notes) & stale_rows$notes != "",
                                               paste0("; ", stale_rows$notes), ""), ")")
      result <- rbind(result, stale_rows[, c("column_name", "column_type", "description", "notes")])
    }
  }

  n_rows <- tryCatch(
    as.integer(DBI::dbGetQuery(pop$db_conn, paste0("SELECT COUNT(*) AS n FROM ", table_name))$n),
    error = function(e) NA_integer_
  )

  structure(
    tibble::as_tibble(result),
    table_name   = table_name,
    table_desc   = tbl_desc,
    n_rows       = n_rows,
    n_cols       = length(all_cols),
    n_locus_cols = length(locus_cols),
    is_wide      = is_wide,
    class        = c("tidybreed_table_desc", "tbl_df", "tbl", "data.frame")
  )
}


#' Print method for tidybreed_table_desc
#'
#' @param x A `tidybreed_table_desc` object from [describe_table()].
#' @param ... Additional arguments (not used).
#' @return `x`, invisibly.
#' @export
print.tidybreed_table_desc <- function(x, ...) {
  width    <- getOption("width", 80L)
  tbl_name <- attr(x, "table_name")
  tbl_desc <- attr(x, "table_desc")
  n_rows   <- attr(x, "n_rows")
  n_cols   <- attr(x, "n_cols")
  n_locus  <- attr(x, "n_locus_cols")
  is_wide  <- attr(x, "is_wide")

  rows_str <- if (is.na(n_rows)) "? rows" else paste0(format(n_rows, big.mark = ","), " rows")
  cols_str <- paste0(format(n_cols, big.mark = ","), " cols")

  left    <- paste0("── ", tbl_name, " ")
  right   <- paste0(" ", rows_str, " · ", cols_str, " ")
  fill_w  <- max(2L, width - nchar(left) - nchar(right))
  cat(left, strrep("─", fill_w), right, "\n", sep = "")

  if (!is.null(tbl_desc) && nzchar(tbl_desc)) {
    cat("  ", tbl_desc, "\n", sep = "")
  }

  if (nrow(x) == 0L) {
    cat("  (no columns to display)\n")
    return(invisible(x))
  }

  cat("\n")

  col_w  <- max(nchar("Column"),      max(nchar(x$column_name), na.rm = TRUE))
  type_w <- max(nchar("Type"),
                max(nchar(ifelse(is.na(x$column_type), "—", x$column_type)), na.rm = TRUE))

  hdr <- paste0(
    "  ", formatC("Column",      width = -col_w),  "  ",
         formatC("Type",         width = -type_w), "  ",
         "Description"
  )
  cat(hdr, "\n")
  sep_w <- max(2L, width - 2L)
  cat("  ", strrep("─", sep_w), "\n", sep = "")

  for (i in seq_len(nrow(x))) {
    ctype <- if (is.na(x$column_type[i])) "(stale)" else x$column_type[i]
    cdesc <- x$description[i]
    if (is.na(cdesc) || cdesc == "") cdesc <- "(no description)"
    cnote <- x$notes[i]
    if (!is.na(cnote) && nzchar(cnote)) cdesc <- paste0(cdesc, "  [", cnote, "]")

    avail_w   <- max(10L, width - col_w - type_w - 8L)
    cdesc_trunc <- if (nchar(cdesc) > avail_w) paste0(substr(cdesc, 1, avail_w - 3L), "...") else cdesc

    cat("  ",
        formatC(x$column_name[i], width = -col_w),  "  ",
        formatC(ctype,            width = -type_w), "  ",
        cdesc_trunc, "\n", sep = "")
  }

  if (is_wide && n_locus > 0L) {
    cat("\n  [", format(n_locus, big.mark = ","),
        " locus columns not shown]\n", sep = "")
  }

  invisible(x)
}


#' Define or update a description for a table or column
#'
#' @description
#' Upserts a description into `_schema_meta`. Use this to document custom
#' columns added via `mutate_table()` or `...` in `add_founders()` etc.
#' Pipe a `tidybreed_table` from [get_table()] as the first argument.
#'
#' @param tbl A `tidybreed_table` from [get_table()]. The table name is
#'   inferred from `tbl$table_name`.
#' @param column_name Character or `NULL`. When `NULL` (default), the
#'   description applies to the table itself. When a column name is given,
#'   the description applies to that specific column.
#' @param description Character. Human-readable description.
#' @param notes Character or `NULL`. Optional supplementary context.
#'
#' @return The `tidybreed_table`, invisibly, enabling back-to-back chained
#'   calls without repeating `get_table()`. The underlying database is modified
#'   in place via the DBI connection — do not assign the result back to `pop`
#'   (use as a side-effect or extract with `tbl$pop`). `schema()` and
#'   `describe_table()` accept either `tidybreed_pop` or `tidybreed_table`.
#'
#' @seealso [describe_table()], [schema()]
#'
#' @examples
#' \dontrun{
#' # Chain multiple column descriptions — no need to repeat get_table()
#' pop |>
#'   get_table("ind_meta") |>
#'   define_schema_description("id_ind",    "Unique individual identifier") |>
#'   define_schema_description("sex",       "Sex of individual (M or F)")   |>
#'   define_schema_description("line_name", "Genetic line name")
#' describe_table(pop, "ind_meta")  # pop still valid; DBI conn is a reference
#'
#' # Table-level description (column_name = NULL)
#' pop |>
#'   get_table("ind_meta") |>
#'   define_schema_description(description = "Individual metadata table")
#' }
#' @export
define_schema_description <- function(tbl, column_name = NULL, description,
                                      notes = NULL) {
  stopifnot(inherits(tbl, "tidybreed_table"))
  pop        <- tbl$pop
  table_name <- tbl$table_name
  validate_tidybreed_pop(pop)

  stopifnot(is.character(description), length(description) == 1L)

  if (!is.null(column_name)) {
    stopifnot(is.character(column_name), length(column_name) == 1L)
    existing_cols <- DBI::dbListFields(pop$db_conn, table_name)
    if (!column_name %in% existing_cols) {
      stop(
        "Column '", column_name, "' not found in table '", table_name, "'. ",
        "Available columns: ", paste(existing_cols, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (!"_schema_meta" %in% pop$tables) {
    warning(
      "This database does not have a _schema_meta table. ",
      "Run migrate_schema_meta(pop) first.",
      call. = FALSE
    )
    return(invisible(tbl))
  }

  entry <- data.frame(
    object_type = if (is.null(column_name)) "table" else "column",
    table_name  = table_name,
    column_name = if (is.null(column_name)) NA_character_ else column_name,
    description = description,
    notes       = if (is.null(notes)) NA_character_ else as.character(notes),
    stringsAsFactors = FALSE
  )
  register_schema_meta(pop$db_conn, entry)

  invisible(tbl)
}


#' Migrate a legacy database to include _schema_meta
#'
#' @description
#' Creates the `_schema_meta` table and populates it with all standard
#' descriptions for databases created before v0.36.0. Safe to call on a
#' database that already has the table (idempotent).
#'
#' @param pop A `tidybreed_pop` object.
#'
#' @return `pop`, invisibly (with `_schema_meta` added to `pop$tables` if
#'   it was not already present).
#'
#' @seealso [schema()], [describe_table()]
#'
#' @examples
#' \dontrun{
#' pop <- restore_pop("old_sim.duckdb")
#' pop <- migrate_schema_meta(pop)
#' schema(pop)
#' }
#' @export
migrate_schema_meta <- function(pop) {
  stopifnot(inherits(pop, "tidybreed_pop"))

  if ("_schema_meta" %in% pop$tables) {
    message("_schema_meta already present in this database.")
    return(invisible(pop))
  }

  DBI::dbExecute(pop$db_conn, "
    CREATE TABLE _schema_meta (
      id_schema_meta INTEGER PRIMARY KEY,
      object_type    VARCHAR NOT NULL,
      table_name     VARCHAR NOT NULL,
      column_name    VARCHAR,
      description    VARCHAR NOT NULL,
      notes          VARCHAR
    )
  ")

  register_schema_meta(pop$db_conn, .genome_layer_descriptions())

  if ("trait_meta" %in% pop$tables) {
    register_schema_meta(pop$db_conn, .trait_layer_descriptions())
  }

  pop$tables <- unique(c(pop$tables, "_schema_meta"))
  message("Added _schema_meta to database.")

  invisible(pop)
}

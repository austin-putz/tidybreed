# Schema metadata functions for tidybreed populations
#
# These functions allow users to view and manage table and column descriptions
# stored directly in the DuckDB database. Descriptions persist across
# `restore_pop()` restores and travel with the `.duckdb` file.
#
# NOTE: this is a plain comment, not roxygen. An `@name schema` block here would
# collide with the exported schema() topic below and overwrite its documentation.


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


#' Pre-built `_schema_meta` descriptions — Genome group
#'
#' Covers genome_meta, genome_map, genome_effects, and the two per-chromosome rule tables.
#' One function per display group in [schema()]; the group vector in
#' `.schema_table_order()` and these functions must name the same tables.
#' @keywords internal
.genome_descriptions <- function() {
  rbind(
    # genome_meta
    .sm_tbl("genome_meta",
            "Locus-level metadata. One row per locus. Stores chromosome assignment, physical position (pos_bp), and chip membership flags added by define_chip(). Genetic-map positions live in the separate genome_map table."),
    .sm_col("genome_meta", "locus_id",
            "Sequential integer primary key; 1 to n_loci in creation order"),
    .sm_col("genome_meta", "locus_name",
            "Human-readable locus identifier (e.g. 'Locus_1', 'rs12345')"),
    .sm_col("genome_meta", "chr",
            "Chromosome number (integer index 1 to n_chr)"),
    .sm_col("genome_meta", "chr_name",
            "Chromosome label string; defaults to chr cast to character"),
    .sm_col("genome_meta", "pos_bp",
            "Physical position in base pairs along the chromosome (BIGINT, 1-based; VCF/PLINK convention)"),
    .sm_col("genome_meta", "founder_allele_freq",
            "Per-locus base allele frequency for founder haplotype sampling; present only when define_founder_haplotypes() has been called"),
    # genome_map
    .sm_tbl("genome_map",
            "Genetic map in long format. One row per (locus x sex x line x map) with a defined genetic position (pos_cM). sex NULL = both sexes; line_name NULL = all lines. define_genome() writes a single default map (map_name 'default'); sex/line/version-specific maps are added as rows. Logical key (locus_id, sex, line_name, map_name)."),
    .sm_col("genome_map", "id_genome_map",
            "Surrogate integer primary key (assigned via next_int_id)"),
    .sm_col("genome_map", "locus_id",
            "Locus identifier; FK to genome_meta.locus_id; internal join/order key"),
    .sm_col("genome_map", "locus_name",
            "Locus name; FK to genome_meta.locus_name; denormalized"),
    .sm_col("genome_map", "sex",
            "Sex this map applies to: NULL = both sexes, 'M' or 'F' = sex-specific map"),
    .sm_col("genome_map", "line_name",
            "Line/breed this map applies to: NULL = all lines, otherwise line-specific"),
    .sm_col("genome_map", "map_name",
            "Map version/identity; default 'default'"),
    .sm_col("genome_map", "pos_cM",
            "Genetic-map position in centiMorgans along the chromosome"),
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
    # chr_inheritance
    .sm_tbl("chr_inheritance",
            "Per-chromosome copy counts, keyed by offspring sex. One row per (chr_name, offspring_sex, line_name). Seeded default is a diploid autosome (from_parent_1 = from_parent_2 = 1). Non-default rules (sex chromosomes, organelles) are set via define_chromosome()."),
    .sm_col("chr_inheritance", "chr_name",
            "Chromosome label; FK to genome_meta.chr_name"),
    .sm_col("chr_inheritance", "offspring_sex",
            "Offspring sex this rule applies to: NULL (all/default), 'M', or 'F'"),
    .sm_col("chr_inheritance", "line_name",
            "Offspring line this rule applies to: NULL (all lines) or a line name; reserved for crossbreeding"),
    .sm_col("chr_inheritance", "from_parent_1",
            "Absolute number of copies inherited from parent_1 (sire) at ploidy 2"),
    .sm_col("chr_inheritance", "from_parent_2",
            "Absolute number of copies inherited from parent_2 (dam) at ploidy 2"),
    # chr_recombination
    .sm_tbl("chr_recombination",
            "Per-chromosome recombination, keyed by producing-parent sex. One row per (chr_name, parent_sex, line_name). Seeded from define_genome()'s genome-wide recombines_M/recombines_F defaults. Non-default rules (Y, W, achiasmy) are set via define_chromosome()."),
    .sm_col("chr_recombination", "chr_name",
            "Chromosome label; FK to genome_meta.chr_name"),
    .sm_col("chr_recombination", "parent_sex",
            "Producing-parent sex this rule applies to: NULL (both), 'M', or 'F'"),
    .sm_col("chr_recombination", "line_name",
            "Producing-parent line this rule applies to: NULL (all lines) or a line name; reserved for crossbreeding"),
    .sm_col("chr_recombination", "recombines",
            "TRUE if the chromosome recombines in this parent sex's meiosis; FALSE = whole-chromosome achiasmy")
  )
}


#' Pre-built `_schema_meta` descriptions — Founders group
#'
#' Covers the founder haplotype pool.
#' One function per display group in [schema()]; the group vector in
#' `.schema_table_order()` and these functions must name the same tables.
#' @keywords internal
.founder_descriptions <- function() {
  rbind(
    # founder_haplotypes
    .sm_tbl("founder_haplotypes",
            "Pool of founder haplotypes in long format (one row per haplotype x locus) sampled by add_founders() to assign phased alleles."),
    .sm_col("founder_haplotypes", "line_name",
            "Founder line label matching add_founders() line_name. NULL = shared pool for all lines."),
    .sm_col("founder_haplotypes", "haplotype_id",
            "Sequential haplotype identifier within the pool (unique per line_name)"),
    .sm_col("founder_haplotypes", "locus_name",
            "Locus name; FK to genome_meta.locus_name"),
    .sm_col("founder_haplotypes", "allele",
            "Allele on this haplotype at this locus: 0 or 1")
  )
}


#' Pre-built `_schema_meta` descriptions — Individuals group
#'
#' Covers individual metadata and the per-individual genome tables.
#' One function per display group in [schema()]; the group vector in
#' `.schema_table_order()` and these functions must name the same tables.
#' @keywords internal
.individual_descriptions <- function() {
  rbind(
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
    .sm_col("ind_meta", "ploidy",
            "Genome ploidy; declared at add_founders() time (must be 2 in this version), computed at add_offspring() time as the sum of each parent's gamete contribution (own_ploidy / 2 per parent)"),
    # ind_haplotype
    .sm_tbl("ind_haplotype",
            "Phased haplotypes in long format. One row per individual x haplotype x locus. parent_origin (1/2) and strand (1 for diploids) identify the copy; line_origin traces the allele's founding line. No DB primary key (dropped for insert speed); (id_ind, parent_origin, strand, locus_id) is unique by construction, guaranteed R-side."),
    .sm_col("ind_haplotype", "id_ind",
            "Individual identifier; FK to ind_meta.id_ind"),
    .sm_col("ind_haplotype", "parent_origin",
            "Which mating parent contributed this copy: 1 = parent_1 (sire), 2 = parent_2 (dam)"),
    .sm_col("ind_haplotype", "strand",
            "Copy index within a parent's contribution; always 1 for diploids"),
    .sm_col("ind_haplotype", "line_origin",
            "Founding genetic line this allele traces back to; NULL only when genuinely untracked"),
    .sm_col("ind_haplotype", "locus_id",
            "Locus identifier; FK to genome_meta.locus_id; physical sort key"),
    .sm_col("ind_haplotype", "locus_name",
            "Locus name; FK to genome_meta.locus_name; denormalized for direct joins to genome_effects"),
    .sm_col("ind_haplotype", "allele",
            "Phased allele on this strand: 0 (reference) or 1 (alternate)"),
    # ind_genotype
    .sm_tbl("ind_genotype",
            "Genotype dosage cache in long format. One row per individual x locus. NOT auto-populated; filled on demand by add_dosage() from ind_haplotype. May be empty or partial. PRIMARY KEY (id_ind, locus_id)."),
    .sm_col("ind_genotype", "id_ind",
            "Individual identifier; FK to ind_meta.id_ind"),
    .sm_col("ind_genotype", "locus_id",
            "Locus identifier; FK to genome_meta.locus_id"),
    .sm_col("ind_genotype", "locus_name",
            "Locus name; FK to genome_meta.locus_name"),
    .sm_col("ind_genotype", "dosage_value",
            "Sum of alleles across strands (0/1/2 for diploids)"),
    # ind_crossover
    .sm_tbl("ind_crossover",
            "Crossover events in long format, one row per crossover drawn during meiosis. Created empty by define_genome(); populated only when add_offspring(store_crossovers = TRUE) (row writes land with the Stage-2 kernel). Absence of a row for a (id_ind, parent_origin, chr) means that gamete's chromosome did not recombine."),
    .sm_col("ind_crossover", "id_crossover",
            "Primary key (BIGINT) assigned R-side via next_row_id() (the BIGINT-safe allocator); wider than other id_ PKs because crossover events are far more numerous"),
    .sm_col("ind_crossover", "id_ind",
            "The offspring; FK to ind_meta.id_ind"),
    .sm_col("ind_crossover", "parent_origin",
            "1 = event in parent_1's (sire's) meiosis, 2 = parent_2's (dam's); matches ind_haplotype.parent_origin"),
    .sm_col("ind_crossover", "chr",
            "Chromosome number"),
    .sm_col("ind_crossover", "chr_name",
            "Chromosome name; FK to genome_meta.chr_name"),
    .sm_col("ind_crossover", "pos_cM",
            "Crossover genetic location in centiMorgans (the exact drawn position)")
  )
}


#' Pre-built `_schema_meta` descriptions — Genetic model group
#'
#' Covers the trait-keyed configuration tables.
#' One function per display group in [schema()]; the group vector in
#' `.schema_table_order()` and these functions must name the same tables.
#' @keywords internal
.genetic_model_descriptions <- function() {
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
            "Variance (diagonal) or covariance (off-diagonal) value")
  )
}


#' Pre-built `_schema_meta` descriptions — Observation model group
#'
#' Covers the phenotype-keyed configuration tables.
#' One function per display group in [schema()]; the group vector in
#' `.schema_table_order()` and these functions must name the same tables.
#' @keywords internal
.observation_model_descriptions <- function() {
  rbind(
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
            "Polynomial order for legendre weight type"),
    # phenotype_effects
    .sm_tbl("phenotype_effects",
            "Fixed and random effect configurations for phenotype models. One row per phenotype x effect. Populated by define_effect_fixed_class(), define_effect_fixed_cov(), define_effect_random()."),
    .sm_col("phenotype_effects", "phenotype_name",
            "Phenotype this effect belongs to; FK to phenotype_meta.phenotype_name"),
    .sm_col("phenotype_effects", "effect_name",
            "Unique label for this effect within the phenotype (e.g. 'sex', 'gen', 'herd')"),
    .sm_col("phenotype_effects", "effect_class",
            "Effect type: 'fixed_class', 'fixed_cov', or 'random'"),
    .sm_col("phenotype_effects", "source_column",
            "Column in source_table used as the grouping or covariate variable"),
    .sm_col("phenotype_effects", "source_table",
            "Table containing source_column; default 'ind_meta'"),
    .sm_col("phenotype_effects", "distribution",
            "Sampling distribution for random effects: 'normal', 'gamma', or 'uniform'"),
    .sm_col("phenotype_effects", "levels_json",
            "JSON map of level to shift value for fixed_class effects (e.g. {\"M\": 30, \"F\": 0})"),
    .sm_col("phenotype_effects", "slope",
            "Regression coefficient for fixed_cov effects"),
    .sm_col("phenotype_effects", "center",
            "Centering value subtracted from the covariate before multiplying by slope"),
    .sm_col("phenotype_effects", "value",
            "Rarely used scalar override"),
    .sm_col("phenotype_effects", "poly_order",
            "Polynomial order for polynomial covariate effects; default 1 (linear)"),
    .sm_col("phenotype_effects", "null_class_action",
            "Behavior when grouping column is NULL: 'skip' excludes the individual"),
    # phenotype_random_effects
    .sm_tbl("phenotype_random_effects",
            "Sampled random effect levels. One row per phenotype x effect x level. Populated by add_phenotype() on first use; subsequent calls reuse the stored value for consistency."),
    .sm_col("phenotype_random_effects", "phenotype_name",
            "Phenotype this random effect belongs to; FK to phenotype_meta.phenotype_name"),
    .sm_col("phenotype_random_effects", "effect_name",
            "Random effect name; FK to phenotype_effects.effect_name"),
    .sm_col("phenotype_random_effects", "level",
            "The grouping level (e.g. herd ID, sire ID) as a string"),
    .sm_col("phenotype_random_effects", "draw_value",
            "The sampled random effect value for this level"),
    .sm_col("phenotype_random_effects", "date_sampled",
            "Date the value was first sampled")
  )
}


#' Pre-built `_schema_meta` descriptions — Selection group
#'
#' Covers selection index definitions.
#' One function per display group in [schema()]; the group vector in
#' `.schema_table_order()` and these functions must name the same tables.
#' @keywords internal
.selection_descriptions <- function() {
  rbind(
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
            "Economic value per unit of this trait")
  )
}


#' Pre-built `_schema_meta` descriptions — Results group
#'
#' Covers simulation output written by the add_* functions.
#' One function per display group in [schema()]; the group vector in
#' `.schema_table_order()` and these functions must name the same tables.
#' @keywords internal
.results_descriptions <- function() {
  rbind(
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


#' Pre-built `_schema_meta` descriptions — System group
#'
#' Covers the _schema_meta system table itself.
#' One function per display group in [schema()]; the group vector in
#' `.schema_table_order()` and these functions must name the same tables.
#' @keywords internal
.system_descriptions <- function() {
  rbind(
    # _schema_meta
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
            "Optional supplementary context (NULL if unused)")
  )
}



#' All pre-built `_schema_meta` descriptions, in display-group order
#'
#' Registered once by [open_pop()]. Description rows are plain metadata and do
#' not require their table to exist yet, so registering the whole set up front is
#' cheaper than re-registering a subset from every `define_*` entry point — and it
#' removes the window in which a later `define_*` call would clobber a user's
#' [define_schema_description()] override.
#' @keywords internal
.all_schema_descriptions <- function() {
  rbind(
    .genome_descriptions(),
    .founder_descriptions(),
    .individual_descriptions(),
    .genetic_model_descriptions(),
    .observation_model_descriptions(),
    .selection_descriptions(),
    .results_descriptions(),
    .system_descriptions()
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


#' Whole-database size label for the `schema()` header
#'
#' `PRAGMA database_size` returns one row of already-formatted strings, so no
#' byte arithmetic is needed. Two cases:
#'
#' - **File-backed**: report `database_size`, plus `wal_size` when the WAL is
#'   non-empty. The WAL is only folded into the file by a `CHECKPOINT`, so a user
#'   looking at the `.duckdb` on disk sees both; reporting only `database_size`
#'   under-reports. `schema()` must not checkpoint — that is a write.
#' - **In-memory** (`db_name = ":memory:"`): `database_size` is `"0 bytes"` and
#'   `block_size` is `0` because there is no file. Report `memory_usage` and
#'   label it as memory, not disk.
#'
#' @param conn A DBI connection.
#' @param db_path Character. `pop$db_path` — the literal `":memory:"` for an
#'   in-memory population.
#' @return A single string such as `"4.0 MiB on disk"`, or `NA_character_` if
#'   the pragma is unavailable.
#' @keywords internal
.schema_size_label <- function(conn, db_path) {
  ds <- tryCatch(DBI::dbGetQuery(conn, "PRAGMA database_size"),
                 error = function(e) NULL)
  if (is.null(ds) || nrow(ds) < 1L) return(NA_character_)
  ds <- ds[1L, , drop = FALSE]

  in_memory <- identical(db_path, ":memory:") ||
    isTRUE(suppressWarnings(as.numeric(ds$block_size)) == 0)
  if (in_memory) {
    if (is.null(ds$memory_usage) || is.na(ds$memory_usage)) return(NA_character_)
    return(paste0(ds$memory_usage, " in memory"))
  }

  if (is.null(ds$database_size) || is.na(ds$database_size)) return(NA_character_)

  # Both terms are on disk (the WAL is a sibling .wal file), so the suffix goes
  # at the end rather than after the first term. Before the first CHECKPOINT
  # `database_size` can read "0 bytes" with the whole population sitting in the
  # WAL — printing only the first term would be actively wrong there.
  wal <- ds$wal_size
  if (!is.null(wal) && !is.na(wal) && !grepl("^0 +bytes$", wal)) {
    return(paste0(ds$database_size, " + ", wal, " WAL on disk"))
  }
  paste0(ds$database_size, " on disk")
}


# ── Exported functions ────────────────────────────────────────────────────────

#' Display groups for [schema()], in pipeline order
#'
#' Option A from `plans/update_schema_print.md`: a hard-coded vector rather than
#' a name-prefix rule or a `_schema_meta` column. A table added later and not
#' registered here degrades *visibly* — it appears under **User tables** — rather
#' than being silently misfiled by a lexical rule that cannot tell
#' `ind_haplotype` (raw genome data) from `ind_tbv` (simulation output).
#'
#' In-group order is workflow order, not alphabetical. `"User tables"` is an
#' intentionally empty slot: unrecognized tables land there, sorted by name, and
#' `"System"` always prints last.
#'
#' Must name the same tables as the `.\*_descriptions()` helpers above.
#'
#' @return A named list of character vectors, in display order.
#' @keywords internal
.schema_table_order <- function() {
  list(
    "Genome"            = c("genome_meta", "genome_map", "genome_effects",
                            "chr_inheritance", "chr_recombination"),
    "Founders"          = c("founder_haplotypes"),
    "Individuals"       = c("ind_meta", "ind_haplotype", "ind_genotype",
                            "ind_crossover"),
    "Genetic model"     = c("trait_meta", "trait_var_comp"),
    "Observation model" = c("phenotype_meta", "phenotype_components",
                            "phenotype_var_comp", "phenotype_effects",
                            "phenotype_random_effects"),
    "Selection"         = c("index_meta"),
    "Results"           = c("ind_tbv", "ind_phenotype", "ind_ebv", "ind_index",
                            "ind_true_index"),
    "User tables"       = character(0),
    "System"            = c("_schema_meta")
  )
}


#' Assign each table to a display group and a pipeline sort position
#'
#' @param table_names Character vector of table names.
#' @return A data.frame with `table_group` and `sort_key` (integer), aligned to
#'   `table_names`.
#' @keywords internal
.schema_group_of <- function(table_names) {
  groups   <- .schema_table_order()
  g_names  <- names(groups)
  user_pos <- match("User tables", g_names)

  grp  <- character(length(table_names))
  rank <- integer(length(table_names))

  for (i in seq_along(table_names)) {
    hit <- FALSE
    for (gi in seq_along(groups)) {
      wi <- match(table_names[i], groups[[gi]])
      if (!is.na(wi)) {
        grp[i]  <- g_names[gi]
        rank[i] <- wi
        hit     <- TRUE
        break
      }
    }
    if (!hit) {
      grp[i]  <- "User tables"
      rank[i] <- NA_integer_
    }
  }

  # Unrecognized tables sort alphabetically within the User tables group.
  is_user <- grp == "User tables"
  if (any(is_user)) {
    rank[is_user] <- rank(table_names[is_user], ties.method = "first")
  }

  data.frame(
    table_group = factor(grp, levels = g_names),
    sort_key    = match(grp, g_names) * 1000L + rank,
    stringsAsFactors = FALSE
  )
}


#' Abbreviate a row count for display
#'
#' `ind_haplotype` alone otherwise dictates the width of the `Rows` column for
#' every table in the printout.
#'
#' @param n Integer vector of row counts.
#' @return Character vector: `"?"` for `NA`, `"1.10M"` style above a million,
#'   comma-grouped below it.
#' @keywords internal
.schema_format_rows <- function(n) {
  vapply(n, function(v) {
    if (is.na(v)) return("?")
    if (v >= 1e9) return(paste0(formatC(v / 1e9, format = "f", digits = 2), "B"))
    if (v >= 1e6) return(paste0(formatC(v / 1e6, format = "f", digits = 2), "M"))
    format(v, big.mark = ",")
  }, character(1))
}


#' Per-table on-disk size, in bytes
#'
#' DuckDB has no direct per-table byte count. `duckdb_tables().estimated_size`
#' looks right but is estimated *cardinality*, not bytes. The workable route is
#' `PRAGMA storage_info()`, counting the distinct blocks a table's segments
#' occupy and multiplying by the database block size.
#'
#' The caller must `CHECKPOINT` first: before one, recently written data carries
#' `block_id = -1` and every table reports zero. That is a write, which is why
#' `sizes` is opt-in on [schema()] rather than always on.
#'
#' @param conn A DBI connection.
#' @param tables Character vector of table names.
#' @return Numeric vector of bytes, aligned to `tables`; `NA` where the pragma
#'   could not be read.
#' @keywords internal
.schema_table_bytes <- function(conn, tables) {
  block_size <- tryCatch(
    as.numeric(DBI::dbGetQuery(conn, "PRAGMA database_size")$block_size[1]),
    error = function(e) NA_real_
  )
  if (is.na(block_size) || block_size <= 0) {
    return(rep(NA_real_, length(tables)))
  }

  vapply(tables, function(tbl) {
    si <- tryCatch(
      DBI::dbGetQuery(conn, paste0(
        "PRAGMA storage_info(", DBI::dbQuoteString(conn, tbl), ")")),
      error = function(e) NULL
    )
    if (is.null(si) || nrow(si) == 0L) return(0)

    ids <- si$block_id
    if (!is.null(si$additional_block_ids)) {
      ids <- c(ids, unlist(si$additional_block_ids, use.names = FALSE))
    }
    ids <- ids[!is.na(ids) & ids >= 0]
    length(unique(ids)) * block_size
  }, numeric(1), USE.NAMES = FALSE)
}


#' Format a byte count for the `schema()` Size column
#'
#' @param b Numeric vector of bytes.
#' @return Character vector; `"?"` for `NA`.
#' @keywords internal
.schema_format_size <- function(b) {
  vapply(b, function(v) {
    if (is.na(v)) return("?")
    if (v >= 1024^3) return(paste0(formatC(v / 1024^3, format = "f", digits = 2), " GiB"))
    paste0(formatC(v / 1024^2, format = "f", digits = 2), " MiB")
  }, character(1))
}


#' The Size-column footnote, required whenever a Size column prints
#'
#' The column is only defensible because the caveat travels with it, so this has
#' no suppression argument. It must state *both* caveats: the block quantization
#' (why small tables are indistinguishable) and the shortfall against the header
#' total (why the column does not add up). A reader shown only one will read the
#' other as a bug.
#'
#' @param width Integer console width.
#' @param db_size Character. The header's whole-database label, referenced so
#'   the shortfall is self-evidently expected rather than looking like loss.
#' @return `NULL`, invisibly; called for the side effect of printing.
#' @keywords internal
.schema_print_size_footnote <- function(width, db_size) {
  total <- if (is.null(db_size) || is.na(db_size)) "the database" else sub(" on disk$", "", db_size)
  cat("  ", strrep("\u2500", max(2L, width - 2L)), "\n", sep = "")
  txt <- paste0(
    "Size is on-disk storage attributed in whole 256 KiB blocks, so small ",
    "tables are not distinguishable: most read as 0.25 MiB, and a table small ",
    "enough to live in the catalog reads as 0.00 MiB while still holding rows. ",
    "Per-table sizes sum to less than the ", total, " database \u2014 catalog ",
    "and header blocks are not attributed to any table."
  )
  cat(paste0("  ", strwrap(txt, width = max(40L, width - 2L))), sep = "\n")
  cat("\n")
  invisible(NULL)
}


#' View table-level descriptions for a tidybreed population
#'
#' @description
#' Returns a tibble of all user-visible tables with row counts, column counts,
#' and descriptions from `_schema_meta`. Print the result for an aligned
#' overview; use `describe_table(pop, "name")` to drill into a specific table.
#'
#' Tables are ordered by pipeline stage — the order in which a population is
#' built — and printed under group headings. The order does not depend on how
#' the population was opened: [open_pop()] populates `pop$tables` in creation
#' order and [restore_pop()] in DuckDB catalog order, so before this the same
#' `.duckdb` file printed differently depending on which one you used.
#'
#' @param pop A `tidybreed_pop` object, or a `tidybreed_table` from
#'   [get_table()] (its `pop` reference is used).
#' @param order One of `"pipeline"` (default, grouped by build stage with
#'   headings), `"name"` (flat alphabetical), `"rows"` (flat, largest row count
#'   first) or `"size"` (flat, largest on-disk size first). Group headings print
#'   only for `"pipeline"`; `table_group` stays in the returned tibble for all
#'   four. `"size"` requires `sizes = TRUE` and errors otherwise.
#' @param show_empty Logical. `FALSE` (default) collapses each group's zero-row
#'   tables into a single `+ n empty: ...` line — on a freshly built population
#'   most tables are empty and would otherwise bury the ones with data. `TRUE`
#'   prints every table on its own row.
#' @param include_system Logical. `FALSE` (default) hides the `_schema_meta`
#'   system table; `TRUE` lists it under a **System** group.
#' @param sizes Logical. `FALSE` (default). `TRUE` adds a per-table
#'   `size_bytes` column and prints a `Size` column with a mandatory caveat
#'   footnote. Opt-in because collecting it requires a `CHECKPOINT`, which
#'   writes to the database — see **Size reporting**.
#'
#' @return A tibble of class `tidybreed_schema` with columns
#'   `table_name`, `table_group`, `n_rows`, `n_cols`, and `description`.
#'   Printed via [print.tidybreed_schema()]. `table_group` is a factor whose
#'   level order is the display order, so the grouping is available as data for
#'   `filter()` / `split()` and not only as printed text.
#'
#' @section Size reporting:
#' `sizes = TRUE` adds a per-table `Size` column. Three things about it are
#' load-bearing:
#'
#' 1. **It writes.** DuckDB has no direct per-table byte count; the workable
#'    route is `PRAGMA storage_info()`, and before a `CHECKPOINT` recently
#'    written data carries `block_id = -1` so every table reports zero. `schema()`
#'    therefore issues a `CHECKPOINT` when `sizes = TRUE`, and only then. This is
#'    the sole reason the argument is opt-in rather than always on.
#' 2. **It is quantized.** Size is on-disk storage attributed in whole 256 KiB
#'    blocks, so small tables are not distinguishable: most read as 0.25 MiB, and
#'    a table small enough to live in the catalog reads as 0.00 MiB while still
#'    holding rows.
#' 3. **It does not add up.** Per-table sizes sum to less than the database total
#'    in the header, because catalog and header blocks are not attributed to any
#'    table. Blocks are not shared between tables, so the attribution is sound as
#'    far as it goes.
#'
#' Caveats 2 and 3 are printed as a footnote under the table whenever the column
#' appears; the footnote has no suppression argument. For an in-memory population
#' there are no blocks at all, so `sizes = TRUE` warns and omits the column.
#'
#' @section Database size:
#' The printed header reports the size of the whole database, read from
#' `PRAGMA database_size`. For a file-backed population this is the `.duckdb`
#' file size, plus the write-ahead log when it is non-empty (the WAL is only
#' folded into the file by a `CHECKPOINT`, which `schema()` deliberately does
#' not issue — it is a write). For an in-memory population there is no file, so
#' the header reports memory usage instead and says so.
#'
#' @seealso [describe_table()], [define_schema_description()]
#'
#' @examples
#' \dontrun{
#' pop <- open_pop(pop_name = "MySim", db_name = ":memory:") |>
#'   define_genome(n_loci = 500, n_chr = 5, chr_len_Mb = 100)
#' schema(pop)
#'
#' # Expand the collapsed empty tables, and show the system table
#' schema(pop, show_empty = TRUE, include_system = TRUE)
#'
#' # The grouping is data, not just print formatting
#' subset(schema(pop), table_group == "Genome")
#'
#' # Flat orderings answer "what is actually big here?"
#' schema(pop, order = "rows")
#'
#' # On-disk sizes (issues a CHECKPOINT; file-backed populations only)
#' schema(pop, order = "size", sizes = TRUE)
#' }
#' @export
schema <- function(pop,
                   order          = c("pipeline", "name", "rows", "size"),
                   show_empty     = FALSE,
                   include_system = FALSE,
                   sizes          = FALSE) {
  if (inherits(pop, "tidybreed_table")) pop <- pop$pop
  stopifnot(inherits(pop, "tidybreed_pop"))
  stopifnot(is.logical(show_empty),     length(show_empty)     == 1L)
  stopifnot(is.logical(include_system), length(include_system) == 1L)
  stopifnot(is.logical(sizes),          length(sizes)          == 1L)
  order <- match.arg(order)

  # order = "size" needs data that only sizes = TRUE collects. Error rather than
  # silently falling back: the requested ordering is undefined, not merely
  # unavailable, and order = "rows" is the honest proxy the caller can ask for.
  if (identical(order, "size") && !sizes) {
    stop("order = \"size\" requires sizes = TRUE. ",
         "Per-table sizes are opt-in because collecting them issues a ",
         "CHECKPOINT, which writes to the database. ",
         "Use order = \"rows\" for a size proxy that costs nothing.",
         call. = FALSE)
  }

  in_memory <- identical(pop$db_path, ":memory:")
  if (sizes && in_memory) {
    # No file means no blocks, so there is nothing to attribute. Warn and omit
    # the column rather than print zeros under a footnote that does not apply.
    if (identical(order, "size")) {
      stop("order = \"size\" is not available for an in-memory population: ",
           "there are no storage blocks to measure. Use order = \"rows\".",
           call. = FALSE)
    }
    warning("sizes = TRUE ignored: this population is in memory, so there are ",
            "no on-disk blocks to attribute. The header reports memory usage.",
            call. = FALSE)
    sizes <- FALSE
  }

  listed <- pop$tables
  if (!include_system) listed <- setdiff(listed, "_schema_meta")

  # Pull table descriptions from _schema_meta
  desc_df <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT table_name, description FROM _schema_meta WHERE object_type = 'table'"
  )

  # Build result row by row
  rows <- lapply(listed, function(tbl) {
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

  # Pipeline grouping. pop$tables arrives in creation order after open_pop() and
  # in DuckDB catalog order after restore_pop(), so the same .duckdb file used to
  # print differently depending on how it was opened. Sorting here removes that.
  grouped            <- .schema_group_of(result$table_name)
  result$table_group <- grouped$table_group

  keep <- c("table_name", "table_group", "n_rows", "n_cols")
  if (sizes) {
    # storage_info() reports zero for anything still sitting in the WAL, so the
    # checkpoint is mandatory — and is the whole reason `sizes` is opt-in.
    DBI::dbExecute(pop$db_conn, "CHECKPOINT")
    result$size_bytes <- .schema_table_bytes(pop$db_conn, result$table_name)
    keep <- c(keep, "size_bytes")
  }
  keep <- c(keep, "description")

  ord_idx <- switch(
    order,
    pipeline = base::order(grouped$sort_key),
    name     = base::order(result$table_name),
    rows     = base::order(-result$n_rows, result$table_name),
    size     = base::order(-result$size_bytes, result$table_name)
  )
  result           <- result[ord_idx, keep, drop = FALSE]
  rownames(result) <- NULL

  structure(
    tibble::as_tibble(result),
    pop_name   = pop$pop_name,
    db_size    = .schema_size_label(pop$db_conn, pop$db_path),
    show_empty = show_empty,
    order      = order,
    sizes      = sizes,
    class      = c("tidybreed_schema", "tbl_df", "tbl", "data.frame")
  )
}


#' Print method for tidybreed_schema
#'
#' @param x A `tidybreed_schema` object from [schema()].
#' @param ... Additional arguments (not used).
#' @return `x`, invisibly.
#' @export
print.tidybreed_schema <- function(x, ...) {
  width      <- getOption("width", 80L)
  pop_name   <- attr(x, "pop_name")
  show_empty <- isTRUE(attr(x, "show_empty"))
  ord        <- attr(x, "order")
  if (is.null(ord)) ord <- "pipeline"

  left   <- paste0("\u2500\u2500 Schema: ", pop_name, " ")

  # Right-hand header summary: table count, and the whole-database size when
  # PRAGMA database_size was readable. Dropped rather than wrapped if the
  # terminal is too narrow to hold it alongside the population name.
  n_tbl    <- nrow(x)
  summary  <- paste0(n_tbl, " table", if (n_tbl == 1L) "" else "s")
  size_lbl <- attr(x, "db_size")
  if (!is.null(size_lbl) && !is.na(size_lbl)) {
    summary <- paste0(summary, " \u00b7 ", size_lbl)
  }
  right <- paste0(" ", summary, " \u2500\u2500")
  if (nchar(left) + nchar(right) + 2L > width) right <- ""

  fill_w <- max(2L, width - nchar(left) - nchar(right))
  cat(left, strrep("\u2500", fill_w), right, "\n", sep = "")
  cat("  Use describe_table(pop, \"name\") for column-level details.\n")

  if (nrow(x) == 0L) {
    cat("\n  (no tables)\n")
    return(invisible(x))
  }

  rows_txt <- .schema_format_rows(x$n_rows)
  cols_txt <- ifelse(is.na(x$n_cols), "?", format(x$n_cols, trim = TRUE))
  has_size <- "size_bytes" %in% names(x)
  size_txt <- if (has_size) .schema_format_size(x$size_bytes) else NULL

  tbl_w  <- max(nchar("Table"), max(nchar(x$table_name)))
  rows_w <- max(nchar("Rows"),  max(nchar(rows_txt)))
  cols_w <- max(nchar("Cols"),  max(nchar(cols_txt)))
  size_w <- if (has_size) max(nchar("Size"), max(nchar(size_txt))) else 0L
  desc_w <- max(10L, width - (4L + tbl_w + 2L + rows_w + 2L + cols_w + 2L +
                              if (has_size) size_w + 2L else 0L))

  emit_row <- function(i, indent) {
    desc <- x$description[i]
    if (is.na(desc) || desc == "") desc <- "(no description)"
    if (nchar(desc) > desc_w) desc <- paste0(substr(desc, 1L, desc_w - 3L), "...")
    cat(indent,
        formatC(x$table_name[i], width = -tbl_w), "  ",
        formatC(rows_txt[i],     width =  rows_w), "  ",
        formatC(cols_txt[i],     width =  cols_w), "  ",
        if (has_size) paste0(formatC(size_txt[i], width = size_w), "  ") else "",
        desc, "\n", sep = "")
  }

  # Collapsed one-liner for a group's zero-row tables. On a freshly built
  # population most tables are empty and would otherwise bury the signal.
  emit_empty <- function(names_vec, indent) {
    prefix <- paste0(indent, "+ ", length(names_vec), " empty: ")
    line   <- prefix
    cur    <- nchar(prefix)
    for (k in seq_along(names_vec)) {
      piece <- paste0(names_vec[k], if (k < length(names_vec)) ", " else "")
      if (cur + nchar(piece) > width && k > 1L) {
        cat(line, "\n", sep = "")
        line <- paste0(strrep(" ", nchar(prefix)))
        cur  <- nchar(line)
      }
      line <- paste0(line, piece)
      cur  <- cur + nchar(piece)
    }
    cat(line, "\n", sep = "")
  }

  if (identical(ord, "pipeline")) {
    for (g in levels(x$table_group)) {
      idx <- which(x$table_group == g)
      if (!length(idx)) next            # never print a heading for an empty group

      cat("\n  ", g, "\n", sep = "")
      is_empty <- !is.na(x$n_rows[idx]) & x$n_rows[idx] == 0L
      shown    <- if (show_empty) idx else idx[!is_empty]
      for (i in shown) emit_row(i, "    ")
      if (!show_empty && any(is_empty)) {
        emit_empty(x$table_name[idx[is_empty]], "    ")
      }
    }
  } else {
    # Flat orderings print a column header instead of group headings.
    cat("\n  ",
        formatC("Table", width = -tbl_w), "  ",
        formatC("Rows",  width =  rows_w), "  ",
        formatC("Cols",  width =  cols_w), "  ",
        if (has_size) paste0(formatC("Size", width = size_w), "  ") else "",
        "Description\n", sep = "")
    cat("  ", strrep("\u2500", max(2L, width - 2L)), "\n", sep = "")

    is_empty <- !is.na(x$n_rows) & x$n_rows == 0L
    shown    <- if (show_empty) seq_len(nrow(x)) else which(!is_empty)
    for (i in shown) emit_row(i, "  ")
    if (!show_empty && any(is_empty)) emit_empty(x$table_name[is_empty], "  ")
  }

  # Mandatory whenever a Size column prints — see .schema_print_size_footnote().
  if (has_size) {
    cat("\n")
    .schema_print_size_footnote(width, attr(x, "db_size"))
  }

  invisible(x)
}


#' View column-level descriptions for a tidybreed table
#'
#' @description
#' Returns a tibble of all columns in `table_name` with their DuckDB types
#' and descriptions from `_schema_meta`. Any wide table with `locus_<n>`
#' columns only lists its non-locus metadata columns to avoid printing
#' thousands of locus columns.
#'
#' @param pop A `tidybreed_pop` object, or a `tidybreed_table` from
#'   [get_table()] (its `pop` reference is used).
#' @param table_name Character. Name of the table to describe.
#'
#' @return A tibble of class `tidybreed_table_desc` with columns
#'   `column_name`, `column_type`, `description`, and `notes`. Printed via
#'   [print.tidybreed_table_desc()].
#'
#' @seealso [schema()], [define_schema_description()]
#'
#' @examples
#' \dontrun{
#' pop <- open_pop(pop_name = "MySim", db_name = ":memory:") |>
#'   define_genome(n_loci = 500, n_chr = 5, chr_len_Mb = 100)
#' describe_table(pop, "ind_meta")
#' describe_table(pop, "genome_effects")
#'
#' # Also accepts a tidybreed_table, so it chains directly after get_table()
#' pop |> get_table("ind_meta") |> describe_table("ind_meta")
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
  if (nrow(meta_rows) > 0) {
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

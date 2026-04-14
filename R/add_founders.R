#' Add founder individuals to population
#'
#' @description
#' Creates founder individuals by sampling haplotypes from the `founder_haplotypes`
#' table. Each founder receives two randomly sampled haplotypes (with replacement),
#' which are used to populate the `genome_haplotype` and `genome_genotype` tables.
#'
#' The `ind_meta` table is created (if it doesn't exist) or appended to with
#' the new founders. Founder individuals have `NULL` for both parent IDs.
#'
#' @param pop A `tidybreed_pop` object
#' @param n_males Integer. Number of male founders to create
#' @param n_females Integer. Number of female founders to create
#' @param pop_name Character. Population identifier used for individual IDs.
#'   IDs are formatted as `"{pop_name}-{number}"` (e.g., "A-1", "A-2")
#'
#' @return The modified `tidybreed_pop` object (invisibly).
#'   **Important:** Assign the result back to update your object: `pop <- add_founders(pop, ...)`
#'
#' @details
#' **Requirements:**
#' - The `founder_haplotypes` table must exist. Create it by calling
#'   `initialize_genome()` with the `n_haplotypes` parameter.
#'
#' **What it does:**
#' 1. Samples 2 haplotypes per founder from `founder_haplotypes` (with replacement)
#' 2. Creates/updates `ind_meta` table with founder metadata
#' 3. Populates `genome_haplotype` (2 rows per individual)
#' 4. Populates `genome_genotype` (1 row per individual, sum of haplotypes)
#'
#' **ID Format:**
#' - Individual IDs: `"{pop_name}-{number}"` (e.g., "A-1", "A-2", "B-1")
#' - Numbers are sequential within each population
#' - If founders already exist for a population, numbering continues from max ID
#'
#' **Multiple Populations:**
#' - Can be called multiple times to add different populations to same database
#' - Each population has independent ID numbering
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Initialize genome with founder haplotypes
#' pop <- initialize_genome(
#'   pop_name = "test",
#'   n_loci = 1000,
#'   n_chr = 10,
#'   chr_len_Mb = 100,
#'   n_haplotypes = 100
#' )
#'
#' # Add founders
#' pop <- pop %>%
#'   add_founders(n_males = 10, n_females = 100, pop_name = "A")
#'
#' # Add custom metadata
#' pop <- pop %>%
#'   mutate_ind_meta(
#'     gen = 0,
#'     farm = "FarmA",
#'     date_birth = Sys.Date()
#'   )
#'
#' # Add second population to same database
#' pop <- pop %>%
#'   add_founders(n_males = 5, n_females = 50, pop_name = "B")
#'
#' # View founders
#' get_table(pop, "ind_meta") %>% collect()
#' }
add_founders <- function(pop, n_males, n_females, pop_name) {

  # ============================================================================
  # 1. Validate inputs
  # ============================================================================

  # Validate pop object
  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)

  # Validate n_males and n_females
  stopifnot(is.numeric(n_males), length(n_males) == 1, n_males >= 0)
  stopifnot(is.numeric(n_females), length(n_females) == 1, n_females >= 0)

  # Ensure at least one founder
  if ((n_males + n_females) == 0) {
    stop(
      "At least one founder must be specified (n_males + n_females > 0)",
      call. = FALSE
    )
  }

  # Validate pop_name
  stopifnot(is.character(pop_name), length(pop_name) == 1, nchar(pop_name) > 0)

  # Validate field name format for pop_name
  if (!grepl("^[a-zA-Z][a-zA-Z0-9_-]*$", pop_name)) {
    stop(
      "pop_name must start with letter and contain only letters, numbers, underscores, or hyphens",
      call. = FALSE
    )
  }

  # Check founder_haplotypes table exists
  if (!"founder_haplotypes" %in% pop$tables) {
    stop(
      "founder_haplotypes table does not exist. ",
      "Call initialize_genome() with n_haplotypes parameter to create founder haplotypes.",
      call. = FALSE
    )
  }

  # ============================================================================
  # 2. Read founder haplotypes
  # ============================================================================

  # Read all founder haplotypes from database (once, into memory)
  founder_haps_tbl <- get_table(pop, "founder_haplotypes") %>%
    dplyr::collect()

  # Get number of available haplotypes
  n_haplotypes <- nrow(founder_haps_tbl)

  if (n_haplotypes == 0) {
    stop("founder_haplotypes table is empty. Cannot sample haplotypes.", call. = FALSE)
  }

  # Get number of loci from founder_haplotypes columns
  locus_cols <- grep("^locus_", colnames(founder_haps_tbl), value = TRUE)
  n_loci <- length(locus_cols)

  if (n_loci == 0) {
    stop("No locus columns found in founder_haplotypes table.", call. = FALSE)
  }

  # Convert to matrix for efficient access
  hap_data_matrix <- as.matrix(founder_haps_tbl[, locus_cols])

  # ============================================================================
  # 3. Determine ID sequence
  # ============================================================================

  # Check if ind_meta already exists
  ind_meta_exists <- "ind_meta" %in% DBI::dbListTables(pop$db_conn)

  # Determine starting ID number for this population
  if (ind_meta_exists) {
    # Query max ID number for this population
    max_num_query <- paste0(
      "SELECT MAX(CAST(SUBSTRING(ind_id FROM POSITION('-' IN ind_id) + 1) AS INTEGER)) as max_num ",
      "FROM ind_meta ",
      "WHERE ind_id LIKE '", pop_name, "-%'"
    )

    max_num_result <- DBI::dbGetQuery(pop$db_conn, max_num_query)

    if (is.na(max_num_result$max_num)) {
      start_id <- 1
    } else {
      start_id <- max_num_result$max_num + 1
    }
  } else {
    start_id <- 1
  }

  # ============================================================================
  # 4. Sample haplotypes
  # ============================================================================

  # Total number of founders
  n_founders <- n_males + n_females

  # Sample haplotypes: 2 per individual, with replacement
  hap_indices <- matrix(
    sample(1:n_haplotypes, size = n_founders * 2, replace = TRUE),
    nrow = n_founders,
    ncol = 2
  )

  # ============================================================================
  # 5. Create ind_meta data
  # ============================================================================

  # Generate individual IDs
  ind_ids <- paste0(pop_name, "-", seq(start_id, start_id + n_founders - 1))

  # Create sex vector
  sex_vector <- c(rep("M", n_males), rep("F", n_females))

  # Create ind_meta data frame
  ind_meta_df <- tibble::tibble(
    ind_id = ind_ids,
    parent_1 = NA_character_,  # NULL for founders
    parent_2 = NA_character_,  # NULL for founders
    population = pop_name,
    sex = sex_vector
  )

  # ============================================================================
  # 6. Create genome_haplotype data
  # ============================================================================

  # Build haplotype data frame (2 rows per individual)
  hap_list <- vector("list", n_founders)

  for (i in 1:n_founders) {
    # Get the two haplotype indices for this individual
    hap_idx_1 <- hap_indices[i, 1]
    hap_idx_2 <- hap_indices[i, 2]

    # Create 2-row data frame for this individual
    hap_list[[i]] <- tibble::tibble(
      ind_id = rep(ind_ids[i], 2),
      parent_origin = c(1L, 2L)
    )

    # Add locus columns
    for (j in 1:n_loci) {
      locus_name <- paste0("locus_", j)
      hap_list[[i]][[locus_name]] <- as.integer(c(
        hap_data_matrix[hap_idx_1, j],
        hap_data_matrix[hap_idx_2, j]
      ))
    }
  }

  # Combine all individuals into single data frame
  genome_haplotype_df <- dplyr::bind_rows(hap_list)

  # ============================================================================
  # 7. Create genome_genotype data
  # ============================================================================

  # Build genotype data frame (1 row per individual)
  geno_df <- tibble::tibble(ind_id = ind_ids)

  for (j in 1:n_loci) {
    locus_name <- paste0("locus_", j)

    # Extract the two haplotypes for this locus across all individuals
    hap1_vals <- hap_data_matrix[hap_indices[, 1], j]
    hap2_vals <- hap_data_matrix[hap_indices[, 2], j]

    # Sum to get genotype (0, 1, or 2)
    geno_df[[locus_name]] <- as.integer(hap1_vals + hap2_vals)
  }

  genome_genotype_df <- geno_df

  # ============================================================================
  # 8. Write data to database
  # ============================================================================

  # Write ind_meta table
  if (ind_meta_exists) {
    DBI::dbWriteTable(pop$db_conn, "ind_meta", ind_meta_df, append = TRUE)
  } else {
    DBI::dbWriteTable(pop$db_conn, "ind_meta", ind_meta_df, overwrite = FALSE)
    # Refresh tables list from database to keep in sync
    pop$tables <- DBI::dbListTables(pop$db_conn)
  }

  # Append to genome_haplotype
  DBI::dbWriteTable(pop$db_conn, "genome_haplotype", genome_haplotype_df, append = TRUE)

  # Append to genome_genotype
  DBI::dbWriteTable(pop$db_conn, "genome_genotype", genome_genotype_df, append = TRUE)

  # ============================================================================
  # 9. Update and return
  # ============================================================================

  # Update metadata
  if (is.null(pop$metadata$n_individuals)) {
    pop$metadata$n_individuals <- n_founders
  } else {
    pop$metadata$n_individuals <- pop$metadata$n_individuals + n_founders
  }

  message("Added ", n_founders, " founders (", n_males, " males, ", n_females, " females) to population '", pop_name, "'")

  # Return modified pop object
  invisible(pop)
}

#' Define a SNP chip
#'
#' @description
#' Convenience function to mark which loci are on a SNP chip. Creates a logical
#' column in the `genome_meta` table indicating chip membership.
#'
#' Three selection methods are supported:
#' - **By count** (`n_snp`): Select a specified number of SNPs using random, even, or chromosome-proportional spacing
#' - **By ID** (`locus_ids`): Select specific loci by their locus_id values
#' - **By name** (`locus_names`): Select specific loci by their locus_name values
#'
#' @param pop A `tidybreed_pop` object
#' @param name Character. Name of the SNP chip (used in messages and default column name)
#' @param n_snp Integer. Number of SNPs to select (mutually exclusive with locus_ids/locus_names)
#' @param locus_ids Integer vector. Specific locus IDs to select (mutually exclusive with n_snp/locus_names)
#' @param locus_names Character vector. Specific locus names to select (mutually exclusive with n_snp/locus_ids)
#' @param method Character. Selection method when using `n_snp`. One of:
#'   - `"random"`: Randomly sample SNPs across entire genome (default)
#'   - `"even"`: Evenly space SNPs across entire genome
#'   - `"chromosome_even"`: Distribute SNPs proportionally across chromosomes
#' @param col_name Character. Column name to create in genome_meta.
#'   Default: `paste0("is_", name)` (e.g., "50k" → "is_50k")
#'
#' @return The modified `tidybreed_pop` object (invisibly).
#'   **Important:** Assign the result back to update your object: `pop <- define_chip(pop, ...)`
#'
#' @details
#' **Selection Methods:**
#'
#' 1. **Random** (`n_snp` + `method = "random"`):
#'    - Randomly samples n_snp loci without replacement
#'    - Uniform distribution across entire genome
#'
#' 2. **Even spacing** (`n_snp` + `method = "even"`):
#'    - Evenly spaces n_snp loci across genome by position
#'    - Based on sequential locus order in genome_meta
#'
#' 3. **Chromosome-even** (`n_snp` + `method = "chromosome_even"`):
#'    - Distributes n_snp proportionally across chromosomes
#'    - Each chromosome gets n_snp * (chr_loci / total_loci) SNPs
#'    - SNPs are evenly spaced within each chromosome
#'
#' 4. **By locus IDs** (`locus_ids`):
#'    - Select specific loci by their locus_id values (1, 2, 3, ...)
#'    - All specified IDs must exist in genome_meta
#'
#' 5. **By locus names** (`locus_names`):
#'    - Select specific loci by their locus_name values
#'    - All specified names must exist in genome_meta
#'
#' **Column naming:**
#' - Default column name is `is_{name}` (e.g., "50k" → "is_50k")
#' - Can override with `col_name` parameter for custom naming
#' - Column contains logical TRUE/FALSE values
#'
#' **Overwriting:**
#' - If column already exists, it will be overwritten (follows `mutate_*` semantics)
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
#'   chr_len_Mb = 100,
#'   n_haplotypes = 100
#' )
#'
#' # Random selection of 500 SNPs
#' pop <- pop %>%
#'   define_chip(name = "50k", n_snp = 500, method = "random")
#'
#' # Evenly spaced 50,000 SNPs
#' pop <- pop %>%
#'   define_chip(name = "HD", n_snp = 50000, method = "even")
#'
#' # Proportional distribution across chromosomes
#' pop <- pop %>%
#'   define_chip(name = "10k", n_snp = 10000, method = "chromosome_even")
#'
#' # Specific loci by ID
#' pop <- pop %>%
#'   define_chip(name = "custom", locus_ids = c(1, 10, 50, 100))
#'
#' # Specific loci by name (e.g., from external chip manifest)
#' chip_manifest <- c("Locus_1", "Locus_10", "Locus_50")
#' pop <- pop %>%
#'   define_chip(name = "custom", locus_names = chip_manifest)
#'
#' # Custom column name
#' pop <- pop %>%
#'   define_chip(name = "bovine_50k", n_snp = 500, col_name = "SNP_50k")
#'
#' # View chip definition
#' get_table(pop, "genome_meta") %>%
#'   select(locus_id, locus_name, chr, pos_Mb, is_50k) %>%
#'   filter(is_50k == TRUE) %>%
#'   collect()
#' }
define_chip <- function(pop,
                        name,
                        n_snp = NULL,
                        locus_ids = NULL,
                        locus_names = NULL,
                        method = "random",
                        col_name = paste0("is_", name)) {

  # ============================================================================
  # 1. Validate inputs
  # ============================================================================

  # Validate pop object
  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)

  # Validate name
  stopifnot(is.character(name), length(name) == 1, nchar(name) > 0)

  # Validate exactly one selection method provided
  methods_provided <- sum(
    !is.null(n_snp),
    !is.null(locus_ids),
    !is.null(locus_names)
  )

  if (methods_provided == 0) {
    stop(
      "Must specify one selection method: n_snp, locus_ids, or locus_names",
      call. = FALSE
    )
  }

  if (methods_provided > 1) {
    stop(
      "Cannot specify multiple selection methods. ",
      "Choose one of: n_snp, locus_ids, locus_names",
      call. = FALSE
    )
  }

  # Validate method parameter (only relevant for n_snp)
  if (!is.null(n_snp)) {
    valid_methods <- c("random", "even", "chromosome_even")
    if (!method %in% valid_methods) {
      stop(
        "Invalid method: '", method, "'. ",
        "Must be one of: ", paste(valid_methods, collapse = ", "),
        call. = FALSE
      )
    }
  }

  # Validate col_name
  stopifnot(is.character(col_name), length(col_name) == 1, nchar(col_name) > 0)

  # ============================================================================
  # 2. Get genome_meta table
  # ============================================================================

  if (!"genome_meta" %in% pop$tables) {
    stop(
      "genome_meta table does not exist. Call initialize_genome() first.",
      call. = FALSE
    )
  }

  genome <- get_table(pop, "genome_meta") %>%
    dplyr::collect()

  n_loci <- nrow(genome)

  if (n_loci == 0) {
    stop("genome_meta table is empty. Cannot define chip.", call. = FALSE)
  }

  # ============================================================================
  # 3. Generate logical vector based on selection method
  # ============================================================================

  chip_indicator <- rep(FALSE, n_loci)

  if (!is.null(n_snp)) {
    # Selection by n_snp
    chip_indicator <- select_by_n_snp(genome, n_snp, method)
    n_selected <- sum(chip_indicator)
    message(
      "Defined chip '", name, "' with ", n_selected, " SNPs ",
      "(method: ", method, ")"
    )

  } else if (!is.null(locus_ids)) {
    # Selection by locus IDs
    chip_indicator <- select_by_locus_ids(genome, locus_ids)
    n_selected <- sum(chip_indicator)
    message(
      "Defined chip '", name, "' with ", n_selected, " SNPs ",
      "(by locus IDs)"
    )

  } else if (!is.null(locus_names)) {
    # Selection by locus names
    chip_indicator <- select_by_locus_names(genome, locus_names)
    n_selected <- sum(chip_indicator)
    message(
      "Defined chip '", name, "' with ", n_selected, " SNPs ",
      "(by locus names)"
    )
  }

  # ============================================================================
  # 4. Create column using mutate_genome_meta()
  # ============================================================================

  # Create named list for mutate_genome_meta
  args <- setNames(list(chip_indicator), col_name)

  # Call mutate_genome_meta with the chip indicator
  result <- do.call(mutate_genome_meta, c(list(pop = pop), args))

  # Return modified pop object
  invisible(result)
}

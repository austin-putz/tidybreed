#' Define a SNP chip
#'
#' @description
#' Convenience function to mark which loci are on a SNP chip. Creates a logical
#' column in the `genome_meta` table indicating chip membership.
#'
#' Four selection methods are supported:
#' - **By count** (`n`): Select a specified number of SNPs using random, even, or chromosome-proportional spacing
#' - **By logical vector** (`locus_tf`): Pass a TRUE/FALSE vector the same length as the number of loci
#' - **By ID** (`locus_ids`): Select specific loci by their locus_id values
#' - **By name** (`locus_names`): Select specific loci by their locus_name values
#'
#' @param pop A `tidybreed_pop` object
#' @param chip_name Character. Name of the SNP chip (used in messages and default column name)
#' @param n Integer. Number of SNPs to select (mutually exclusive with locus_tf/locus_ids/locus_names)
#' @param locus_tf Logical vector. TRUE/FALSE membership vector with one element per locus, in the
#'   same order as `genome_meta`. Useful when you have already computed membership (e.g., pass
#'   `!pop_50k$is_50k` to select the complement of an existing chip). Mutually exclusive with
#'   n/locus_ids/locus_names.
#' @param locus_ids Integer vector. Specific locus IDs to select (mutually exclusive with n/locus_tf/locus_names)
#' @param locus_names Character vector. Specific locus names to select (mutually exclusive with n/locus_tf/locus_ids)
#' @param method Character. Selection method when using `n`. One of:
#'   - `"random"`: Randomly sample SNPs across entire genome (default)
#'   - `"even"`: Evenly space SNPs across entire genome
#'   - `"chromosome_even"`: Distribute SNPs proportionally across chromosomes
#' @param col_name Character. Column name to create in genome_meta.
#'   Default: `paste0("is_", chip_name)` (e.g., "50k" → "is_50k")
#'
#' @return The modified `tidybreed_pop` object (invisibly).
#'   **Important:** Assign the result back to update your object: `pop <- define_chip(pop, ...)`
#'
#' @details
#' **Selection Methods:**
#'
#' 1. **Random** (`n` + `method = "random"`):
#'    - Randomly samples n loci without replacement
#'    - Uniform distribution across entire genome
#'
#' 2. **Even spacing** (`n` + `method = "even"`):
#'    - Evenly spaces n loci across genome by position
#'    - Based on sequential locus order in genome_meta
#'
#' 3. **Chromosome-even** (`n` + `method = "chromosome_even"`):
#'    - Distributes n proportionally across chromosomes
#'    - Each chromosome gets n * (chr_loci / total_loci) SNPs
#'    - SNPs are evenly spaced within each chromosome
#'
#' 4. **By logical vector** (`locus_tf`):
#'    - Pass a logical vector with one element per locus (same length as genome_meta)
#'    - TRUE marks the locus as on the chip, FALSE marks it as off
#'    - Useful for building chips from existing chip membership, e.g., `!existing_tf` for the complement
#'
#' 5. **By locus IDs** (`locus_ids`):
#'    - Select specific loci by their locus_id values (1, 2, 3, ...)
#'    - All specified IDs must exist in genome_meta
#'
#' 6. **By locus names** (`locus_names`):
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
#'   define_chip(chip_name = "50k", n = 500, method = "random")
#'
#' # Evenly spaced SNPs
#' pop <- pop %>%
#'   define_chip(chip_name = "HD", n = 900, method = "even")
#'
#' # Proportional distribution across chromosomes
#' pop <- pop %>%
#'   define_chip(chip_name = "10k", n = 100, method = "chromosome_even")
#'
#' # Logical vector — complement of an existing chip
#' chip_tf <- get_table(pop, "genome_meta") %>%
#'   dplyr::pull(is_50k)
#' pop <- pop %>%
#'   define_chip(chip_name = "non50k", locus_tf = !chip_tf)
#'
#' # Specific loci by ID
#' pop <- pop %>%
#'   define_chip(chip_name = "custom", locus_ids = c(1, 10, 50, 100))
#'
#' # Specific loci by name (e.g., from external chip manifest)
#' chip_manifest <- c("Locus_1", "Locus_10", "Locus_50")
#' pop <- pop %>%
#'   define_chip(chip_name = "custom", locus_names = chip_manifest)
#'
#' # Custom column name
#' pop <- pop %>%
#'   define_chip(chip_name = "bovine_50k", n = 500, col_name = "SNP_50k")
#'
#' # View chip definition
#' get_table(pop, "genome_meta") %>%
#'   select(locus_id, locus_name, chr, pos_Mb, is_50k) %>%
#'   filter(is_50k == TRUE) %>%
#'   collect()
#' }
define_chip <- function(pop,
                        chip_name,
                        n = NULL,
                        locus_tf = NULL,
                        locus_ids = NULL,
                        locus_names = NULL,
                        method = "random",
                        col_name = paste0("is_", chip_name)) {

  # ============================================================================
  # 1. Validate inputs
  # ============================================================================

  # Validate pop object
  stopifnot(inherits(pop, "tidybreed_pop"))
  validate_tidybreed_pop(pop)

  # Validate chip_name
  stopifnot(is.character(chip_name), length(chip_name) == 1, nchar(chip_name) > 0)

  # Validate exactly one selection method provided
  methods_provided <- sum(
    !is.null(n),
    !is.null(locus_tf),
    !is.null(locus_ids),
    !is.null(locus_names)
  )

  if (methods_provided == 0) {
    stop(
      "Must specify one selection method: n, locus_tf, locus_ids, or locus_names",
      call. = FALSE
    )
  }

  if (methods_provided > 1) {
    stop(
      "Cannot specify multiple selection methods. ",
      "Choose one of: n, locus_tf, locus_ids, locus_names",
      call. = FALSE
    )
  }

  # Validate method parameter (only relevant for n)
  if (!is.null(n)) {
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

  if (!is.null(n)) {
    # Selection by count
    chip_indicator <- select_by_n(genome, n, method)
    n_selected <- sum(chip_indicator)
    message(
      "Defined chip '", chip_name, "' with ", n_selected, " SNPs ",
      "(method: ", method, ")"
    )

  } else if (!is.null(locus_tf)) {
    # Selection by logical vector
    if (!is.logical(locus_tf)) {
      stop("locus_tf must be a logical vector", call. = FALSE)
    }
    if (length(locus_tf) != n_loci) {
      stop(
        "locus_tf length (", length(locus_tf), ") must equal number of loci (",
        n_loci, ")",
        call. = FALSE
      )
    }
    chip_indicator <- locus_tf
    n_selected <- sum(chip_indicator)
    message("Defined chip '", chip_name, "' with ", n_selected, " SNPs (by locus_tf)")

  } else if (!is.null(locus_ids)) {
    # Selection by locus IDs
    chip_indicator <- select_by_locus_ids(genome, locus_ids)
    n_selected <- sum(chip_indicator)
    message(
      "Defined chip '", chip_name, "' with ", n_selected, " SNPs ",
      "(by locus IDs)"
    )

  } else if (!is.null(locus_names)) {
    # Selection by locus names
    chip_indicator <- select_by_locus_names(genome, locus_names)
    n_selected <- sum(chip_indicator)
    message(
      "Defined chip '", chip_name, "' with ", n_selected, " SNPs ",
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

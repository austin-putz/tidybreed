#' Add founder individuals to population
#'
#' @description
#' Creates founder individuals by sampling haplotypes from the `founder_haplotypes`
#' table. Each founder receives two randomly sampled haplotypes (with replacement),
#' which are used to populate the long `ind_haplotype` table (`ind_genotype`
#' dosage is materialized on demand via [add_dosage()]).
#'
#' The `ind_meta` table is created (if it doesn't exist) or appended to with
#' the new founders. Founder individuals have `NULL` for both parent IDs.
#'
#' @param tbl A \code{tidybreed_table} from \code{get_table("founder_haplotypes")}
#'   (optionally piped through \code{dplyr::filter()}). The filtered rows supply
#'   the haplotype pool for sampling; use the filter to select a line-specific
#'   or custom subset.
#' @param n_males Integer. Number of male founders to create
#' @param n_females Integer. Number of female founders to create
#' @param line_name Character. Line identifier used for individual IDs.
#'   IDs are formatted as `"{line_name}_{number}"` (e.g., "A_1", "A_2")
#' @param ... Optional named arguments for custom `ind_meta` columns, e.g.
#'   `gen = 0L`, `farm = "Iowa"`. Scalar values are broadcast to all new
#'   founders; vectors must have length `n_males + n_females`. Column types
#'   are inferred from the R type: use `0L` for INTEGER, `0` for DOUBLE,
#'   `"text"` for VARCHAR, `TRUE`/`FALSE` for BOOLEAN. Reserved column names
#'   (`id_ind`, `sex`, `line`, etc.) are blocked.
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
#' 3. Populates `ind_haplotype` (long: 2 haplotypes x n_loci rows per individual,
#'    with `line_origin` set to the founder's line and `strand = 1`)
#'
#' **ID Format:**
#' - Individual IDs: `"{line_name}_{number}"` (e.g., "A_1", "A_2", "B_1")
#' - Numbers are sequential within each line
#' - If founders already exist for a line, numbering continues from max ID
#'
#' **Multiple Lines:**
#' - Can be called multiple times to add different lines to same database
#' - Each line has independent ID numbering
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pop <- open_pop(pop_name = "test", db_name = ":memory:") |>
#'   define_genome(n_loci = 1000, n_chr = 10, chr_len_Mb = 100) |>
#'   define_founder_haplotypes(n_haplotypes = 100, line_name = "A")
#'
#' # Simple — all haplotypes in the table
#' pop <- pop |>
#'   get_table("founder_haplotypes") |>
#'   add_founders(n_males = 10, n_females = 100, line_name = "A")
#'
#' # Filtered by line
#' pop <- pop |>
#'   get_table("founder_haplotypes") |>
#'   dplyr::filter(line_name == "Yorkshire") |>
#'   add_founders(n_males = 10, n_females = 50, line_name = "Yorkshire")
#'
#' # With custom ind_meta columns via ...
#' pop <- pop |>
#'   get_table("founder_haplotypes") |>
#'   dplyr::filter(line_name == "A") |>
#'   add_founders(n_males = 10, n_females = 100, line_name = "A",
#'                gen = 0L, farm = "FarmA")
#'
#' # Add a second line
#' pop <- pop |>
#'   get_table("founder_haplotypes") |>
#'   dplyr::filter(line_name == "B") |>
#'   add_founders(n_males = 5, n_females = 50, line_name = "B")
#'
#' # View founders
#' get_table(pop, "ind_meta") |> dplyr::collect()
#' }
add_founders <- function(tbl, n_males, n_females, line_name, ...) {

  # ============================================================================
  # 1. Validate inputs
  # ============================================================================

  if (!inherits(tbl, "tidybreed_table")) {
    stop(
      "'tbl' must be a tidybreed_table from get_table('founder_haplotypes') |> filter(...). ",
      "Use: pop |> get_table('founder_haplotypes') |> add_founders()",
      call. = FALSE
    )
  }
  if (tbl$table_name != "founder_haplotypes") {
    stop(
      "'tbl' must be piped from get_table('founder_haplotypes'), not '",
      tbl$table_name, "'.",
      call. = FALSE
    )
  }
  pop <- tbl$pop
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

  # Validate line_name
  stopifnot(is.character(line_name), length(line_name) == 1, nchar(line_name) > 0)

  # Validate field name format for line_name
  if (!grepl("^[a-zA-Z][a-zA-Z0-9_]*$", line_name)) {
    stop(
      "line_name must start with letter and contain only letters, numbers, or underscores",
      call. = FALSE
    )
  }

  # ============================================================================
  # 2. Read founder haplotypes
  # ============================================================================

  founder_haps_long <- dplyr::collect(tbl)
  if (nrow(founder_haps_long) == 0L) {
    stop("No haplotypes selected — your filter returned zero rows.", call. = FALSE)
  }

  # Reconstruct a haplotype x locus allele matrix from the long founder pool.
  # locus_id (from genome_meta) drives column order; haplotypes are keyed by
  # (line_name, haplotype_id) since haplotype_id is only unique within a line.
  gm <- DBI::dbGetQuery(pop$db_conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id")
  n_loci     <- nrow(gm)
  locus_cols <- paste0("locus_", seq_len(n_loci))

  founder_haps_long$locus_id <-
    gm$locus_id[match(founder_haps_long$locus_name, gm$locus_name)]

  key_df <- unique(founder_haps_long[, c("line_name", "haplotype_id")])
  key_df <- key_df[order(key_df$line_name, key_df$haplotype_id), , drop = FALSE]
  uniq_keys <- paste(key_df$line_name, key_df$haplotype_id, sep = "\r")
  n_haplotypes <- length(uniq_keys)

  row_key <- paste(founder_haps_long$line_name, founder_haps_long$haplotype_id,
                   sep = "\r")
  hap_data_matrix <- matrix(0L, nrow = n_haplotypes, ncol = n_loci)
  hap_data_matrix[cbind(match(row_key, uniq_keys), founder_haps_long$locus_id)] <-
    as.integer(founder_haps_long$allele)

  # ============================================================================
  # 3. Determine ID sequence
  # ============================================================================

  # ind_meta always exists (created by initialize_genome()); query current max.
  max_num_result <- DBI::dbGetQuery(pop$db_conn, paste0(
    "SELECT MAX(CAST(SUBSTRING(id_ind FROM POSITION('_' IN id_ind) + 1) AS INTEGER)) as max_num ",
    "FROM ind_meta WHERE LEFT(id_ind, ", nchar(line_name) + 1L, ") = '", line_name, "_'"
  ))
  start_id <- if (is.na(max_num_result$max_num)) 1L else as.integer(max_num_result$max_num) + 1L

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
  ind_ids <- paste0(line_name, "_", seq(start_id, start_id + n_founders - 1))

  # Create sex vector
  sex_vector <- c(rep("M", n_males), rep("F", n_females))

  # Create ind_meta data frame
  ind_meta_df <- tibble::tibble(
    id_ind = ind_ids,
    id_parent_1 = NA_character_,  # NULL for founders
    id_parent_2 = NA_character_,  # NULL for founders
    line_name = line_name,
    sex = sex_vector
  )

  # Attach any user-supplied custom columns from ...
  extra_cols <- list(...)
  if (length(extra_cols) > 0) {
    prepped <- prepare_extra_cols(extra_cols, n_founders, "ind_meta", pop$db_conn)
    for (nm in names(prepped)) ind_meta_df[[nm]] <- prepped[[nm]]
  }

  # ============================================================================
  # 6. Build the wide allele frame (transient; the vehicle for the UNPIVOT into
  #    the long ind_haplotype). It is never written to a wide DB table.
  # ============================================================================

  # Interleave paternal/maternal row indices: [pat1, mat1, pat2, mat2, ...]
  row_idx <- c(rbind(hap_indices[, 1], hap_indices[, 2]))
  hap_matrix <- hap_data_matrix[row_idx, , drop = FALSE]
  storage.mode(hap_matrix) <- "integer"

  hap_wide <- cbind(
    data.frame(
      id_ind        = rep(ind_ids, each = 2),
      parent_origin = rep(c(1L, 2L), times = n_founders),
      stringsAsFactors = FALSE
    ),
    as.data.frame(hap_matrix)
  )
  colnames(hap_wide)[3:(2 + n_loci)] <- locus_cols

  # ============================================================================
  # 7. Write data to database
  # ============================================================================

  # Write ind_meta table (always append; table was created by initialize_genome())
  DBI::dbWriteTable(pop$db_conn, "ind_meta", ind_meta_df, append = TRUE)

  # Write the long ind_haplotype by unpivoting the wide allele frame in DuckDB
  # (so R never materializes the long frame). strand = 1 (diploid); line_origin
  # = the founder's own line for every allele. register + INSERT is RNG-neutral.
  # ind_genotype stays on-demand (add_dosage()); nothing is written to it here.
  duckdb::duckdb_register(pop$db_conn, "__tmp_hap", hap_wide)
  DBI::dbExecute(pop$db_conn, paste0(
    "INSERT INTO ind_haplotype ",
    "(id_ind, parent_origin, strand, line_origin, locus_id, locus_name, allele) ",
    "SELECT u.id_ind, u.parent_origin, 1, '", line_name, "', ",
    "gm.locus_id, gm.locus_name, u.allele ",
    "FROM (UNPIVOT __tmp_hap ON COLUMNS(* EXCLUDE (id_ind, parent_origin)) ",
    "INTO NAME locus_col VALUE allele) u ",
    "JOIN genome_meta gm ON u.locus_col = 'locus_' || gm.locus_id"
  ))
  duckdb::duckdb_unregister(pop$db_conn, "__tmp_hap")

  # ============================================================================
  # 9. Update and return
  # ============================================================================

  message("Added ", n_founders, " founders (", n_males, " males, ", n_females, " females) to line '", line_name, "'")

  # Return modified pop object
  invisible(pop)
}

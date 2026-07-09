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
#' @param ploidy Integer scalar. Genome ploidy for these founders. Must be `2`
#'   in this version of tidybreed (real polyploidy is not yet supported).
#' @param ... Optional named arguments for custom `ind_meta` columns, e.g.
#'   `gen = 0L`, `farm = "Iowa"`. Scalar values are broadcast to all new
#'   founders; vectors must have length `n_males + n_females`. Column types
#'   are inferred from the R type: use `0L` for INTEGER, `0` for DOUBLE,
#'   `"text"` for VARCHAR, `TRUE`/`FALSE` for BOOLEAN. Reserved column names
#'   (`id_ind`, `sex`, `line_name`, `ploidy`, etc.) are blocked.
#' @param batch_size Optional integer. Founders materialized + written per batch.
#'   Bounds peak memory at roughly `batch_size x n_loci` long rows regardless of
#'   the total number of founders. Any value produces byte-identical output for a
#'   fixed seed (only the *write* is batched, never the sampling). Overrides
#'   `max_batch_mem` when both are set.
#' @param max_batch_mem Optional per-batch memory budget as bytes (e.g. `512e6`)
#'   or a string (e.g. `"512MB"`). Used to derive `batch_size` when it is `NULL`.
#'   When both are `NULL` (the default), the batch size is auto-picked from
#'   detected available system memory (conservative fallback if unavailable).
#'
#' @return The modified `tidybreed_pop` object (invisibly).
#'   **Important:** Assign the result back to update your object: `pop <- add_founders(pop, ...)`
#'
#' @details
#' **Requirements:**
#' - The `founder_haplotypes` table must exist. Create it by calling
#'   [define_founder_haplotypes()] (after [open_pop()] and [define_genome()]).
#'
#' **What it does:**
#' 1. Samples 2 haplotypes per founder from `founder_haplotypes` (with replacement)
#' 2. Creates/updates `ind_meta` table with founder metadata
#' 3. Populates `ind_haplotype` (long format; `line_origin` set to the founder's
#'    line and `strand = 1`. Row count per chromosome follows
#'    `chr_meta`'s `copy_mode_M`/`copy_mode_F` for the founder's sex — 2
#'    rows/locus for `"full"` (the default, diploid autosomes), 1 for
#'    `"half"` (e.g. sex chromosomes), 0 for `"none"`; see [define_chr()])
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
#' # Filtered by line (line-specific pools need their own
#' # define_founder_haplotypes() call before they can be filtered to)
#' pop <- pop |>
#'   define_founder_haplotypes(n_haplotypes = 100, line_name = "Yorkshire")
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
#'   define_founder_haplotypes(n_haplotypes = 50, line_name = "B")
#' pop <- pop |>
#'   get_table("founder_haplotypes") |>
#'   dplyr::filter(line_name == "B") |>
#'   add_founders(n_males = 5, n_females = 50, line_name = "B")
#'
#' # View founders
#' get_table(pop, "ind_meta") |> dplyr::collect()
#' }
add_founders <- function(tbl, n_males, n_females, line_name, ploidy = 2L, ...,
                         batch_size = NULL, max_batch_mem = NULL) {

  # ============================================================================
  # 1. Validate inputs
  # ============================================================================

  stopifnot(is.numeric(ploidy), length(ploidy) == 1L)
  if (as.integer(ploidy) != 2L) {
    stop(
      "ploidy must be 2 in this version of tidybreed (sex chromosomes and ",
      "organelles are supported; real polyploidy is not yet supported). ",
      "Got ploidy = ", ploidy, ".",
      call. = FALSE
    )
  }

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
    "SELECT locus_id, locus_name, chr, chr_name FROM genome_meta ORDER BY locus_id")
  n_loci     <- nrow(gm)

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

  # ind_meta always exists (created by open_pop()); query current max.
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
    sex = sex_vector,
    ploidy = as.integer(ploidy)
  )

  # Attach any user-supplied custom columns from ...
  extra_cols <- list(...)
  if (length(extra_cols) > 0) {
    prepped <- prepare_extra_cols(extra_cols, n_founders, "ind_meta", pop$db_conn)
    for (nm in names(prepped)) ind_meta_df[[nm]] <- prepped[[nm]]
  }

  # ============================================================================
  # 6. Batched long write of founder haplotypes (+ ind_meta), one transaction.
  #
  # Founders are whole-haplotype pool draws — no meiosis. The single sample()
  # in step 4 drew every haplotype index up front; batching here only
  # materializes and writes the selected long rows, never the sampling, so
  # seeded pools stay byte-identical for any batch size (the founder RNG trap).
  # Peak memory is bounded by ~B x n_loci long rows, not the full
  # (2 x n_founders) x n_loci selected matrix (which is never built).
  #
  # copy_mode_M/copy_mode_F decides which (sex, parent_origin) slots are written
  # per chromosome — "full" writes both parent_origin slots (diploid-autosome
  # default), "half" writes only the hemi_parent slot, "none" writes nothing —
  # the same rule the old per-chr UNPIVOT applied, now emitted directly in long
  # form. strand = 1 (ploidy 2); line_origin = the founder's own line for every
  # allele. .write_long_haplotypes() is register+INSERT (RNG-neutral); the single
  # dbWriteTable(ind_meta) advances R's RNG exactly once, matching the pre-Stage-1
  # stream. ind_genotype stays on-demand (add_dosage()); nothing written here.
  # ============================================================================

  # hap_data_matrix column c holds locus_id c; founder i's paternal copy is
  # pool row hap_indices[i, 1], maternal is hap_indices[i, 2].
  pat_row <- hap_indices[, 1]
  mat_row <- hap_indices[, 2]

  chr_meta_map  <- get_chr_meta_map(pop$db_conn)
  chr_names_uni <- unique(gm$chr_name)

  B    <- resolve_batch_size(n_founders, n_loci, batch_size, max_batch_mem)
  conn <- pop$db_conn

  DBI::dbExecute(conn, "BEGIN TRANSACTION")
  committed <- FALSE
  on.exit(if (!committed) try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE),
          add = TRUE)

  DBI::dbWriteTable(conn, "ind_meta", ind_meta_df, append = TRUE)

  for (b_start in seq.int(1L, n_founders, by = B)) {
    b_idx <- b_start:min(b_start + B - 1L, n_founders)
    b_ids <- ind_ids[b_idx]
    b_sex <- sex_vector[b_idx]
    b_pat <- pat_row[b_idx]
    b_mat <- mat_row[b_idx]

    # Assemble the batch's long rows by preallocation (no do.call(rbind, parts)):
    # pass 1 enumerates the emitted (chr x parent_origin) blocks and their row
    # counts; pass 2 fills preallocated typed vectors by slice, then builds one
    # data.frame. Byte-identical values/order to the prior parts/rbind build,
    # but without the O(n) rbind copies and per-part data.frame overhead.
    blocks <- list()
    total  <- 0L
    for (cn in chr_names_uni) {
      cols    <- which(gm$chr_name == cn)   # column positions (= locus_id)
      nl      <- length(cols)
      chr_row <- chr_meta_map[[cn]]

      keep_po_for <- function(cm) {
        if (cm == "none") integer(0)
        else if (cm == "full") c(1L, 2L)
        else if (chr_row$hemi_parent == "parent_1") 1L else 2L
      }
      poM <- keep_po_for(chr_row$copy_mode_M)
      poF <- keep_po_for(chr_row$copy_mode_F)

      for (po in c(1L, 2L)) {
        emit <- (b_sex == "M" & po %in% poM) | (b_sex == "F" & po %in% poF)
        sub  <- which(emit)
        if (!length(sub)) next
        n_blk <- length(sub) * nl
        blocks[[length(blocks) + 1L]] <- list(
          po = po, sub = sub, cols = cols, lids = gm$locus_id[cols],
          nl = nl, n = n_blk)
        total <- total + n_blk
      }
    }

    if (total > 0L) {
      id_ind_v        <- character(total)
      parent_origin_v <- integer(total)
      locus_id_v      <- integer(total)
      allele_v        <- integer(total)
      pos <- 0L
      for (blk in blocks) {
        sub <- blk$sub
        ls  <- length(sub)
        idx <- (pos + 1L):(pos + blk$n)
        pool_rows <- if (blk$po == 1L) b_pat[sub] else b_mat[sub]
        block <- hap_data_matrix[pool_rows, blk$cols, drop = FALSE]  # ls x nl
        id_ind_v[idx]        <- rep(b_ids[sub], times = blk$nl)
        parent_origin_v[idx] <- blk$po
        locus_id_v[idx]      <- rep(blk$lids, each = ls)
        allele_v[idx]        <- as.integer(as.vector(block))
        pos <- pos + blk$n
      }
      long_frame <- data.frame(
        id_ind        = id_ind_v,
        parent_origin = parent_origin_v,
        strand        = 1L,
        line_origin   = line_name,
        locus_id      = locus_id_v,
        allele        = allele_v,
        stringsAsFactors = FALSE
      )
      .write_long_haplotypes(conn, long_frame)
      rm(long_frame, id_ind_v, parent_origin_v, locus_id_v, allele_v)
    }
    rm(blocks)
    gc(verbose = FALSE)
  }

  DBI::dbExecute(conn, "COMMIT")
  committed <- TRUE

  # ============================================================================
  # 9. Update and return
  # ============================================================================

  message("Added ", n_founders, " founders (", n_males, " males, ", n_females, " females) to line '", line_name, "'")

  # Return modified pop object
  invisible(pop)
}

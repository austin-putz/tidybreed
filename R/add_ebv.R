#' Add estimated breeding values to a tidybreed population
#'
#' @description
#' Computes EBVs using one of two modes:
#' * `software = "blupf90"` — runs the BLUPF90 suite (`renumf90` + `blupf90+`).
#'   Both binaries must be on PATH. Writes all input/output files to a
#'   timestamped subfolder of `run_dir` for inspection and reproducibility.
#' * `parent_avg = TRUE` — simple average of parent EBVs already stored in
#'   `ind_ebv` for the given `model`. Returns `NA` for any animal whose
#'   parents do not both have EBVs.
#'
#' The first argument is a `tidybreed_table` (from [get_table()] and
#' optionally [filter()]). The filter controls which animals' **phenotypic
#' records** are included in the evaluation; ancestors are traced back
#' automatically for the pedigree via `n_gen_pedigree`.
#'
#' @param tbl A `tidybreed_table` from [get_table()] (optionally filtered).
#'   Must contain an `id_ind` column.
#' @param trait Character vector of trait name(s) to evaluate.
#' @param software Character or `NULL`. Pass `"blupf90"` to run the BLUPF90
#'   suite. Use together with `parent_avg = FALSE`.
#' @param parent_avg Logical. If `TRUE`, compute EBVs as the average of
#'   parent EBVs from `ind_ebv` (for the given `model`). Cannot be combined
#'   with `software`.
#' @param model Character. Label stored in `ind_ebv$model`. Must be a valid
#'   SQL identifier.  Default `"model_1"`.
#' @param overwrite_trait Logical. If `TRUE`, all existing `ind_ebv` rows for
#'   the traits being added are deleted before inserting; new rows receive
#'   `eval_number = 1`. Default `FALSE`.
#' @param delete_all Logical. If `TRUE`, **all** rows in `ind_ebv` are deleted
#'   before inserting; new rows receive `eval_number = 1`. Takes precedence
#'   over `overwrite_trait`. Default `FALSE`.
#' @param n_gen_pedigree Integer. Generations to trace back from the filtered
#'   candidates when building the pedigree file. Default `2`.
#' @param chip_name Character or `NULL`. Name of a SNP chip defined via
#'   [define_chip()] (e.g. `"HD"`). When provided, genotype file is written
#'   and single-step GBLUP (ssGBLUP) is requested via `SNP_FILE` in
#'   `renum.par`. Animals must have `has_{chip_name} = TRUE` in `ind_meta`.
#' @param run_dir Character. Base directory for evaluation sub-folders.
#'   Default `"."` (working directory).
#' @param eval_id Character or `NULL`. Sub-folder name suffix appended to
#'   `"eval_"`. Auto-generated from `model` + timestamp when `NULL`.
#' @param estimate_var Logical. If `TRUE`, request variance component
#'   estimation (`OPTION method VCE`) instead of BLUP-only. Default `FALSE`.
#' @param update_covars Logical. If `estimate_var = TRUE`, write estimated
#'   variance components back to `trait_effect_cov`. Default `FALSE`.
#'
#' @return The modified `tidybreed_pop` (invisibly). EBVs are appended to
#'   `ind_ebv` with an auto-incrementing `eval_number` per trait. Use
#'   `overwrite_trait = TRUE` or `delete_all = TRUE` to clear previous records
#'   instead of accumulating them.
#'
#' @seealso [get_table()], [add_phenotype()], [define_chip()],
#'   [add_effect_cov_matrix()]
#'
#' @examples
#' \dontrun{
#' # Parent average (parents must already have EBVs in ind_ebv for model "pa_v0")
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 2L) |>
#'   add_ebv("ADG", parent_avg = TRUE, model = "pa_v1")
#'
#' # BLUPF90 — renumf90 and blupf90+ must be on PATH
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen >= 3L) |>
#'   add_ebv(c("ADG", "BW"),
#'           software = "blupf90",
#'           model    = "blup_v1",
#'           run_dir  = "./eval_runs",
#'           n_gen_pedigree = 2)
#' }
#' @export
add_ebv <- function(tbl,
                    trait,
                    software       = NULL,
                    parent_avg     = FALSE,
                    model          = "model_1",
                    overwrite_trait = FALSE,
                    delete_all     = FALSE,
                    n_gen_pedigree = 2L,
                    chip_name      = NULL,
                    run_dir        = ".",
                    eval_id        = NULL,
                    estimate_var   = FALSE,
                    update_covars  = FALSE,
                    ...) {

  # ---- Input validation ----
  stopifnot(inherits(tbl, "tidybreed_table"))
  pop <- tbl$pop
  validate_tidybreed_pop(pop)

  stopifnot(is.character(trait), length(trait) >= 1)
  lapply(trait, validate_sql_identifier, what = "trait name")
  validate_sql_identifier(model, what = "model label")

  extra_cols <- list(...)
  if (length(extra_cols) > 0) {
    for (nm in names(extra_cols)) {
      if (length(extra_cols[[nm]]) != 1L) {
        stop("Custom field '", nm, "' in add_ebv() must be a scalar ",
             "(broadcast to all records). Supply per-record vectors with ",
             "mutate_table() after the call.", call. = FALSE)
      }
    }
  }

  use_blupf90    <- identical(software, "blupf90")
  use_parent_avg <- isTRUE(parent_avg)

  if (!is.null(software) && !software %in% "blupf90")
    stop("Unknown software '", software, "'. Currently supported: 'blupf90'.", call. = FALSE)
  if (use_blupf90 && use_parent_avg)
    stop("Specify either software = 'blupf90' OR parent_avg = TRUE, not both.", call. = FALSE)
  if (!use_blupf90 && !use_parent_avg)
    stop("Specify one mode: software = 'blupf90' or parent_avg = TRUE.", call. = FALSE)

  # Validate traits exist
  meta_check <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT trait_name FROM trait_meta WHERE trait_name IN (",
           paste0("'", trait, "'", collapse = ", "), ")")
  )
  missing_t <- setdiff(trait, meta_check$trait_name)
  if (length(missing_t) > 0)
    stop("Traits not found in trait_meta: ", paste(missing_t, collapse = ", "), call. = FALSE)

  # Validate chip_name
  if (!is.null(chip_name)) {
    genome_cols   <- DBI::dbListFields(pop$db_conn, "genome_meta")
    ind_meta_cols <- DBI::dbListFields(pop$db_conn, "ind_meta")
    if (!paste0("is_", chip_name) %in% genome_cols)
      stop("Column 'is_", chip_name, "' not found in genome_meta. ",
           "Call define_chip('", chip_name, "', ...) first.", call. = FALSE)
    if (!paste0("has_", chip_name) %in% ind_meta_cols)
      stop("Column 'has_", chip_name, "' not found in ind_meta. ",
           "Call add_genotypes(..., chip_name = '", chip_name, "') first.", call. = FALSE)
  }

  # ---- Resolve candidate set ----
  if (length(tbl$pending_filter) == 0) {
    subset_ids <- DBI::dbGetQuery(pop$db_conn,
                                  "SELECT DISTINCT id_ind FROM ind_meta")$id_ind
  } else {
    collected <- dplyr::collect(tbl)
    if (!"id_ind" %in% names(collected))
      stop("Filtered table '", tbl$table_name,
           "' must contain 'id_ind' to subset individuals.", call. = FALSE)
    subset_ids <- unique(collected[["id_ind"]])
  }

  if (length(subset_ids) == 0) {
    warning("No individuals matched; no EBVs computed.", call. = FALSE)
    return(invisible(pop))
  }

  # ---- Compute eval_number per trait and optionally delete old rows ----
  if (isTRUE(delete_all)) {
    DBI::dbExecute(pop$db_conn, "DELETE FROM ind_ebv")
    eval_nums <- stats::setNames(rep(1L, length(trait)), trait)
  } else if (isTRUE(overwrite_trait)) {
    trait_sql <- paste0("'", trait, "'", collapse = ", ")
    DBI::dbExecute(pop$db_conn,
      paste0("DELETE FROM ind_ebv WHERE trait_name IN (", trait_sql, ")"))
    eval_nums <- stats::setNames(rep(1L, length(trait)), trait)
  } else {
    eval_nums <- stats::setNames(integer(length(trait)), trait)
    for (t in trait) {
      res <- DBI::dbGetQuery(pop$db_conn,
        paste0("SELECT COALESCE(MAX(eval_number), 0) + 1 AS next_eval ",
               "FROM ind_ebv WHERE trait_name = '", t, "'"))
      eval_nums[[t]] <- as.integer(res$next_eval)
    }
  }

  # ---- Route to mode ----
  if (use_parent_avg) {
    ebv_df <- ebv_parent_avg(pop, subset_ids, trait, model, eval_nums)
  } else {
    ebv_df <- ebv_blupf90(pop, subset_ids, trait, model, eval_nums,
                          n_gen_pedigree = as.integer(n_gen_pedigree),
                          chip_name      = chip_name,
                          run_dir        = run_dir,
                          eval_id        = eval_id,
                          estimate_var   = estimate_var,
                          update_covars  = update_covars)
  }

  # ---- Upsert into ind_ebv ----
  if (!is.null(ebv_df) && nrow(ebv_df) > 0) {
    if (length(extra_cols) > 0) {
      prepped <- prepare_extra_cols(extra_cols, nrow(ebv_df), "ind_ebv",
                                   pop$db_conn)
      for (nm in names(prepped)) ebv_df[[nm]] <- prepped[[nm]]
    }
    upsert_ind_ebv(pop, ebv_df)
    n_written <- sum(!is.na(ebv_df$ebv))
    message("Stored ", nrow(ebv_df), " EBV record(s) in ind_ebv",
            if (n_written < nrow(ebv_df))
              paste0(" (", nrow(ebv_df) - n_written, " NA)"),
            " — model = '", model, "'.")
  }

  invisible(pop)
}


# ---------------------------------------------------------------------------
# Parent average
# ---------------------------------------------------------------------------

ebv_parent_avg <- function(pop, subset_ids, trait, model, eval_nums) {
  id_in    <- paste0("'", subset_ids, "'", collapse = ", ")
  cand_meta <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT id_ind, id_parent_1, id_parent_2 FROM ind_meta ",
           "WHERE id_ind IN (", id_in, ")")
  )

  parent_ids <- unique(c(cand_meta$id_parent_1, cand_meta$id_parent_2))
  parent_ids <- parent_ids[!is.na(parent_ids)]

  if (length(parent_ids) > 0) {
    pid_in   <- paste0("'", parent_ids, "'", collapse = ", ")
    trait_in <- paste0("'", trait, "'", collapse = ", ")
    # Get latest EBV per parent per (trait, model); flag if multiple evals exist
    parent_ebv <- DBI::dbGetQuery(
      pop$db_conn,
      paste0(
        "SELECT e.id_ind, e.trait_name, e.ebv, agg.max_eval, agg.n_evals ",
        "FROM ind_ebv e ",
        "JOIN ( ",
        "  SELECT id_ind, trait_name, ",
        "         MAX(eval_number) AS max_eval, ",
        "         COUNT(DISTINCT eval_number) AS n_evals ",
        "  FROM ind_ebv ",
        "  WHERE model = '", model, "' ",
        "  AND trait_name IN (", trait_in, ") ",
        "  AND id_ind IN (", pid_in, ") ",
        "  GROUP BY id_ind, trait_name ",
        ") agg ON e.id_ind = agg.id_ind ",
        "      AND e.trait_name = agg.trait_name ",
        "      AND e.eval_number = agg.max_eval ",
        "WHERE e.model = '", model, "'"
      )
    )
    # Warn per trait if any parent has multiple evaluations
    if (nrow(parent_ebv) > 0) {
      multi <- unique(parent_ebv$trait_name[parent_ebv$n_evals > 1])
      for (t in multi) {
        max_ev <- max(parent_ebv$max_eval[parent_ebv$trait_name == t], na.rm = TRUE)
        warning(
          "Parent EBVs for trait '", t, "': multiple evaluations found. ",
          "Using the latest (eval_number = ", max_ev, ").",
          call. = FALSE
        )
      }
    }
  } else {
    warning("No parents found for the candidate animals. All EBVs will be NA.",
            call. = FALSE)
    parent_ebv <- data.frame(id_ind = character(), trait_name = character(),
                             ebv = numeric(), max_eval = integer(),
                             n_evals = integer(), stringsAsFactors = FALSE)
  }

  n  <- length(subset_ids) * length(trait)
  id_out      <- character(n)
  trait_out   <- character(n)
  ebv_out     <- numeric(n)
  eval_out    <- integer(n)
  k <- 1L

  for (id in subset_ids) {
    row <- cand_meta[cand_meta$id_ind == id, , drop = FALSE]
    p1  <- if (nrow(row) > 0 && !is.na(row$id_parent_1)) row$id_parent_1 else NA_character_
    p2  <- if (nrow(row) > 0 && !is.na(row$id_parent_2)) row$id_parent_2 else NA_character_

    for (t in trait) {
      p1_ebv <- if (!is.na(p1)) {
        v <- parent_ebv$ebv[parent_ebv$id_ind == p1 & parent_ebv$trait_name == t]
        if (length(v) == 1L) v else NA_real_
      } else NA_real_

      p2_ebv <- if (!is.na(p2)) {
        v <- parent_ebv$ebv[parent_ebv$id_ind == p2 & parent_ebv$trait_name == t]
        if (length(v) == 1L) v else NA_real_
      } else NA_real_

      if (is.na(p1_ebv) || is.na(p2_ebv)) {
        warning("Individual '", id, "' trait '", t,
                "': one or both parent EBVs missing — EBV set to NA.",
                call. = FALSE)
        pa_val <- NA_real_
      } else {
        pa_val <- mean(c(p1_ebv, p2_ebv))
      }

      id_out[k]    <- id
      trait_out[k] <- t
      ebv_out[k]   <- pa_val
      eval_out[k]  <- eval_nums[[t]]
      k <- k + 1L
    }
  }

  tibble::tibble(
    id_ind      = id_out,
    trait_name  = trait_out,
    model       = model,
    ebv         = ebv_out,
    acc         = NA_real_,
    se          = NA_real_,
    eval_number = eval_out
  )
}


# ---------------------------------------------------------------------------
# BLUPF90 path
# ---------------------------------------------------------------------------

ebv_blupf90 <- function(pop, subset_ids, trait, model, eval_nums,
                         n_gen_pedigree, chip_name, run_dir, eval_id,
                         estimate_var, update_covars) {

  renumf90_path <- find_blupf90_binary("renumf90")
  blupf90_path  <- find_blupf90_binary("blupf90+")

  if (is.null(eval_id))
    eval_id <- paste0(model, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  eval_dir <- file.path(normalizePath(run_dir, mustWork = FALSE),
                        paste0("eval_", eval_id))
  dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)
  message("Evaluation folder: ", eval_dir)

  # Trace pedigree
  message("Tracing pedigree (", n_gen_pedigree, " generation(s) back)...")
  ped_df      <- trace_pedigree(pop, subset_ids, n_gen_pedigree)
  all_ped_ids <- unique(ped_df$id_ind)
  message("Pedigree: ", nrow(ped_df), " animals (",
          length(subset_ids), " candidates + ",
          length(all_ped_ids) - length(subset_ids), " ancestors).")

  # Write pedigree file
  write.table(ped_df,
              file      = file.path(eval_dir, "pedigree.txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = " ")

  # Build and write data file
  message("Building data file...")
  data_res <- build_data_file(pop, subset_ids, trait, eval_dir)

  # Write genotype file
  geno_res <- NULL
  if (!is.null(chip_name)) {
    message("Writing genotype file (chip: ", chip_name, ")...")
    geno_res <- write_geno_file(pop, all_ped_ids, chip_name, eval_dir)
  }

  # Determine alpha_size (max ID length, rounded up to multiple of 5)
  alpha_size <- ceiling((max(nchar(all_ped_ids)) + 5) / 5) * 5
  if (!is.null(geno_res) && geno_res$id_width > alpha_size)
    alpha_size <- geno_res$id_width

  # Write renum.par
  message("Writing renum.par...")
  animal_effect_num <- write_renum_par(
    eval_dir         = eval_dir,
    col_map          = data_res$col_map,
    distinct_effects = data_res$distinct_effects,
    effects_df       = data_res$effects_df,
    trait            = trait,
    pop              = pop,
    chip_name        = chip_name,
    estimate_var     = estimate_var,
    alpha_size       = alpha_size
  )

  # Write meta.txt
  write_meta_file(
    eval_dir          = eval_dir,
    eval_id           = eval_id,
    col_map           = data_res$col_map,
    distinct_effects  = data_res$distinct_effects,
    trait             = trait,
    effects_df        = data_res$effects_df,
    chip_name         = chip_name,
    n_loci            = if (!is.null(geno_res)) geno_res$n_loci  else 0L,
    id_width          = if (!is.null(geno_res)) geno_res$id_width else 0L,
    animal_effect_num = animal_effect_num
  )

  # Run renumf90
  message("Running renumf90...")
  run_renumf90(eval_dir, renumf90_path)
  message("renumf90 completed.")

  # Run blupf90+
  message("Running blupf90+...")
  run_blupf90_plus(eval_dir, blupf90_path)
  message("blupf90+ completed.")

  # Parse solutions
  message("Parsing solutions (animal effect = ", animal_effect_num, ")...")
  ebv_df <- parse_blupf90_solutions(
    eval_dir          = eval_dir,
    trait             = trait,
    animal_effect_num = animal_effect_num,
    all_ped_ids       = all_ped_ids,
    model             = model,
    eval_nums         = eval_nums
  )

  # Optional VCE writeback
  if (isTRUE(estimate_var) && isTRUE(update_covars))
    update_covars_from_blupf90(pop, eval_dir, trait)

  ebv_df
}


# ---------------------------------------------------------------------------
# Upsert ind_ebv
# ---------------------------------------------------------------------------

upsert_ind_ebv <- function(pop, ebv_df) {
  if (nrow(ebv_df) == 0) return(invisible(NULL))

  tmp <- paste0("_ebv_tmp_", as.character(round(as.numeric(Sys.time()) * 1000)))
  duckdb::duckdb_register(pop$db_conn, tmp, as.data.frame(ebv_df))
  on.exit(duckdb::duckdb_unregister(pop$db_conn, tmp), add = TRUE)

  cols        <- names(ebv_df)
  key_cols    <- c("id_ind", "trait_name", "model", "eval_number")
  update_cols <- setdiff(cols, key_cols)

  col_list   <- paste(cols, collapse = ", ")
  update_set <- paste(
    paste0(update_cols, " = EXCLUDED.", update_cols),
    collapse = ", "
  )
  DBI::dbExecute(pop$db_conn, paste0(
    "INSERT INTO ind_ebv (", col_list, ") ",
    "SELECT ", col_list, " FROM ", tmp, " ",
    "ON CONFLICT (id_ind, trait_name, model, eval_number) DO UPDATE SET ", update_set
  ))
  invisible(NULL)
}

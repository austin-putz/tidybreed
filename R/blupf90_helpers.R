#' Internal helpers for BLUPF90 evaluation
#'
#' @keywords internal
#' @name blupf90_helpers
NULL


#' Find a BLUPF90 binary on PATH
#' @keywords internal
find_blupf90_binary <- function(name) {
  path <- Sys.which(name)
  if (!nzchar(path))
    stop("'", name, "' not found in PATH. ",
         "Install the BLUPF90 suite and ensure the binary is on PATH.",
         call. = FALSE)
  path
}


#' Trace pedigree n generations back from a set of candidate IDs
#'
#' @param pop tidybreed_pop
#' @param subset_ids character vector of starting animal IDs
#' @param n_gen integer number of generations to trace back
#' @return data.frame with columns id_ind, id_parent_1, id_parent_2
#'   (unknown parents replaced with "0")
#' @keywords internal
trace_pedigree <- function(pop, subset_ids, n_gen) {
  all_ids <- subset_ids
  for (gen in seq_len(n_gen)) {
    id_in <- paste0("'", all_ids, "'", collapse = ", ")
    new_parents <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT DISTINCT id_parent_1, id_parent_2 FROM ind_meta ",
             "WHERE id_ind IN (", id_in, ")")
    )
    new_parent_ids <- unique(c(new_parents$id_parent_1, new_parents$id_parent_2))
    new_parent_ids <- new_parent_ids[!is.na(new_parent_ids) & !new_parent_ids %in% all_ids]
    if (length(new_parent_ids) == 0) break
    all_ids <- c(all_ids, new_parent_ids)
  }
  id_in <- paste0("'", all_ids, "'", collapse = ", ")
  ped <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT id_ind, id_parent_1, id_parent_2 FROM ind_meta ",
           "WHERE id_ind IN (", id_in, ")")
  )
  ped$id_parent_1[is.na(ped$id_parent_1)] <- "0"
  ped$id_parent_2[is.na(ped$id_parent_2)] <- "0"
  ped
}


#' Build the BLUPF90 data file
#'
#' Column layout: mu(1), id_ind(2), fixed_class effects, fixed_cov effects,
#' trait observations (missing → 0).
#'
#' @param pop tidybreed_pop
#' @param subset_ids character vector of animal IDs whose phenotypes to include
#' @param trait character vector of trait names
#' @param eval_dir path to evaluation folder; writes data.txt there
#' @return list with data, col_map, distinct_effects, effects_df,
#'   n_fixed_effects, trait_cols
#' @keywords internal
build_data_file <- function(pop, subset_ids, trait, eval_dir) {
  n_traits <- length(trait)

  # Pull fixed effects for these traits (fixed_class and fixed_cov only)
  effects_df <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT trait_name, effect_name, effect_class, source_column, source_table ",
           "FROM trait_effects ",
           "WHERE trait_name IN (", paste0("'", trait, "'", collapse = ", "), ") ",
           "AND effect_class IN ('fixed_class', 'fixed_cov') ",
           "ORDER BY effect_class, effect_name")
  )

  # Distinct effects (by effect_name; warn if source_column differs across traits)
  if (nrow(effects_df) > 0) {
    de <- effects_df[!duplicated(effects_df$effect_name),
                     c("effect_name", "effect_class", "source_column", "source_table"),
                     drop = FALSE]
    de <- de[order(factor(de$effect_class, levels = c("fixed_class", "fixed_cov"))), ]
    rownames(de) <- NULL
  } else {
    de <- data.frame(effect_name = character(), effect_class = character(),
                     source_column = character(), source_table = character(),
                     stringsAsFactors = FALSE)
  }
  n_eff <- nrow(de)

  # Column numbers
  mu_col    <- 1L
  id_col    <- 2L
  eff_cols  <- if (n_eff > 0) seq(3L, 2L + n_eff) else integer(0)
  trt_start <- 3L + n_eff
  trt_cols  <- seq(trt_start, trt_start + n_traits - 1L)

  col_map <- c(mu = mu_col, id_ind = id_col)
  if (n_eff > 0) col_map <- c(col_map, setNames(eff_cols, de$effect_name))
  col_map <- c(col_map, setNames(trt_cols, trait))

  # Pull phenotypic data: average per id_ind per trait
  id_in <- paste0("'", subset_ids, "'", collapse = ", ")
  pheno_long <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT id_ind, trait_name, AVG(value) AS value FROM ind_phenotype ",
           "WHERE trait_name IN (", paste0("'", trait, "'", collapse = ", "), ") ",
           "AND id_ind IN (", id_in, ") ",
           "GROUP BY id_ind, trait_name")
  )
  if (nrow(pheno_long) == 0)
    stop("No phenotypic records found for the requested individuals and traits.",
         call. = FALSE)

  # Pivot wider manually
  all_ids_with_pheno <- unique(pheno_long$id_ind)
  pheno_wide <- data.frame(id_ind = all_ids_with_pheno, stringsAsFactors = FALSE)
  for (t in trait) {
    t_rows <- pheno_long[pheno_long$trait_name == t, c("id_ind", "value"), drop = FALSE]
    pheno_wide <- merge(pheno_wide, t_rows, by = "id_ind", all.x = TRUE)
    names(pheno_wide)[names(pheno_wide) == "value"] <- t
  }
  # Replace NA with 0 (BLUPF90 missing indicator)
  for (t in trait) pheno_wide[[t]][is.na(pheno_wide[[t]])] <- 0

  # Build data frame
  data_df <- data.frame(mu = 1L, id_ind = pheno_wide$id_ind,
                        stringsAsFactors = FALSE)

  # Add fixed-effect columns from source tables
  if (n_eff > 0) {
    unique_src_tables <- unique(de$source_table)
    for (src_tbl in unique_src_tables) {
      eff_sub <- de[de$source_table == src_tbl, , drop = FALSE]
      src_cols <- unique(eff_sub$source_column)
      id_in_df <- paste0("'", pheno_wide$id_ind, "'", collapse = ", ")
      src_data <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT id_ind, ", paste(src_cols, collapse = ", "),
               " FROM ", src_tbl, " WHERE id_ind IN (", id_in_df, ")")
      )
      for (eff_row in seq_len(nrow(eff_sub))) {
        eff  <- eff_sub[eff_row, ]
        col_vals <- src_data[[eff$source_column]][
          match(pheno_wide$id_ind, src_data$id_ind)]
        if (is.character(col_vals)) {
          col_vals <- gsub(" ", "_", col_vals)
          col_vals[is.na(col_vals)] <- "0"
        }
        data_df[[eff$effect_name]] <- col_vals
      }
    }
  }

  # Add trait columns
  for (t in trait) data_df[[t]] <- pheno_wide[[t]]

  # Final column order
  col_order <- c("mu", "id_ind",
                 if (n_eff > 0) de$effect_name else character(0),
                 trait)
  data_df <- data_df[, col_order, drop = FALSE]

  # Write data.txt (header row + SKIP_HEADER 1 in renum.par)
  write.table(data_df, file = file.path(eval_dir, "data.txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = " ")

  list(
    data             = data_df,
    col_map          = col_map,
    distinct_effects = de,
    effects_df       = effects_df,
    n_fixed_effects  = n_eff,
    trait_cols       = trt_cols
  )
}


#' Write the BLUPF90 genotype file (fixed-width ID + concatenated 0/1/2 string)
#'
#' @param pop tidybreed_pop
#' @param all_ped_ids character vector of all pedigree animals
#' @param chip_name character chip name (e.g. "HD")
#' @param eval_dir path to evaluation folder
#' @return list with geno_file path, n_loci, id_width, geno_ids
#' @keywords internal
write_geno_file <- function(pop, all_ped_ids, chip_name, eval_dir) {
  has_col <- paste0("has_", chip_name)
  is_col  <- paste0("is_", chip_name)
  id_in   <- paste0("'", all_ped_ids, "'", collapse = ", ")

  geno_ids <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT id_ind FROM ind_meta WHERE id_ind IN (", id_in, ") ",
           "AND ", has_col, " = TRUE")
  )$id_ind

  if (length(geno_ids) == 0) {
    message("No genotyped animals found in the evaluation set for chip '", chip_name, "'.")
    return(list(geno_file = NULL, n_loci = 0L, id_width = 0L, geno_ids = character(0)))
  }

  chip_loci <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT locus_id FROM genome_meta WHERE ", is_col, " = TRUE ORDER BY locus_id")
  )$locus_id
  n_loci    <- length(chip_loci)
  locus_cols <- paste0("locus_", chip_loci)

  id_in_geno <- paste0("'", geno_ids, "'", collapse = ", ")
  geno_data  <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT id_ind, ", paste(locus_cols, collapse = ", "),
           " FROM genome_genotype WHERE id_ind IN (", id_in_geno, ")")
  )

  id_width <- ceiling((max(nchar(geno_data$id_ind)) + 5) / 5) * 5

  geno_path <- file.path(eval_dir, "genotype.txt")
  lines <- character(nrow(geno_data))
  for (i in seq_len(nrow(geno_data))) {
    geno_str <- paste0(unlist(geno_data[i, locus_cols, drop = FALSE]), collapse = "")
    lines[i] <- sprintf("%-*s%s", id_width, geno_data$id_ind[i], geno_str)
  }
  writeLines(lines, con = geno_path)

  message("Wrote genotype file: ", n_loci, " loci, ", nrow(geno_data), " animals.")
  list(geno_file = geno_path, n_loci = n_loci, id_width = id_width, geno_ids = geno_ids)
}


#' Write the renumf90 parameter file (renum.par)
#'
#' @return integer animal_effect_num (effect number to extract EBVs from solutions)
#' @keywords internal
write_renum_par <- function(eval_dir, col_map, distinct_effects, effects_df,
                             trait, pop, chip_name, estimate_var, alpha_size) {
  n_traits     <- length(trait)
  n_fixed_effs <- nrow(distinct_effects)

  # Load variance components
  R_mat <- load_effect_cov(pop, "residual", trait)
  G_mat <- load_effect_cov(pop, "gen_add",  trait)

  if (is.null(R_mat))
    stop("Residual covariance matrix not found for traits: ",
         paste(trait, collapse = ", "),
         ". Call add_effect_cov_matrix(pop, 'residual', ...) first.", call. = FALSE)
  if (is.null(G_mat))
    stop("Additive genetic covariance matrix not found for traits: ",
         paste(trait, collapse = ", "),
         ". Call add_effect_cov_matrix(pop, 'gen_add', ...) first.", call. = FALSE)

  # Format matrix rows: one line per row
  fmt_mat <- function(mat) {
    apply(mat, 1, function(row) {
      paste(format(as.numeric(row), scientific = FALSE, digits = 9, trim = TRUE),
            collapse = " ")
    })
  }

  # animal_effect_num: 1 (mu) + n_fixed_effs + 1 (RANDOM animal)
  animal_effect_num <- 1L + n_fixed_effs + 1L

  # TRAITS line: column numbers of observations
  trait_cols_str <- paste(col_map[trait], collapse = " ")

  # Build EFFECT blocks
  effect_blocks <- character(0)

  # mu: always effect 1 (column 1, all traits, covariate — always 1)
  mu_cols_str  <- paste(rep(1L, n_traits), collapse = " ")
  effect_blocks <- c(effect_blocks, "EFFECT", paste0(mu_cols_str, " cov"))

  # Fixed effects
  if (n_fixed_effs > 0) {
    for (i in seq_len(n_fixed_effs)) {
      eff <- distinct_effects[i, ]
      # Per-trait column: use the effect's column if trait has it, else 0
      cols_per_trait <- sapply(trait, function(t) {
        has_it <- effects_df$trait_name == t & effects_df$effect_name == eff$effect_name
        if (any(has_it)) col_map[[eff$effect_name]] else 0L
      })
      cols_str <- paste(cols_per_trait, collapse = " ")
      type_str <- if (eff$effect_class == "fixed_cov") "cov" else "cross alpha"
      effect_blocks <- c(effect_blocks, "EFFECT", paste0(cols_str, " ", type_str))
    }
  }

  # RANDOM animal block — must be preceded by its own EFFECT block (column 2 = id_ind)
  id_cols_str  <- paste(rep(col_map[["id_ind"]], n_traits), collapse = " ")
  random_block <- c(
    "EFFECT", paste0(id_cols_str, " cross alpha"),
    "RANDOM", "animal",
    "FILE",   "pedigree.txt",
    "SKIP_HEADER", "1",
    "FILE_POS", "1 2 3 0 0"
  )
  if (!is.null(chip_name)) {
    random_block <- c(random_block, "SNP_FILE", "genotype.txt")
  }
  random_block <- c(random_block,
    "PED_DEPTH", "0",
    "INBREEDING", "pedigree",
    "(CO)VARIANCES", fmt_mat(G_mat)
  )

  alpha_sz  <- max(20L, as.integer(ceiling(alpha_size / 5) * 5))
  method_str <- if (isTRUE(estimate_var)) "VCE" else "BLUP"

  par_lines <- c(
    "DATAFILE",         "data.txt",
    "SKIP_HEADER",      "1",
    "TRAITS",           trait_cols_str,
    "FIELDS_PASSED TO OUTPUT", "2",
    "WEIGHT(S)",        "",
    "RESIDUAL_VARIANCE", fmt_mat(R_mat),
    effect_blocks,
    random_block,
    paste0("OPTION alpha_size ", alpha_sz),
    "OPTION missing 0",
    paste0("OPTION method ", method_str),
    "OPTION origID",
    "OPTION remove_all_missing",
    ""
  )

  writeLines(par_lines, file.path(eval_dir, "renum.par"))
  animal_effect_num
}


#' Write the human-readable meta.txt for the evaluation folder
#' @keywords internal
write_meta_file <- function(eval_dir, eval_id, col_map, distinct_effects,
                             trait, effects_df, chip_name, n_loci, id_width,
                             animal_effect_num) {
  lines <- c(
    paste0("=== tidybreed add_ebv() Evaluation: ", eval_id, " ==="),
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "Data file: data.txt",
    "  Column 1: mu           [intercept, value = 1 for all rows]",
    "  Column 2: id_ind       [animal ID; FIELDS_PASSED TO OUTPUT = 2]"
  )
  if (nrow(distinct_effects) > 0) {
    for (i in seq_len(nrow(distinct_effects))) {
      eff     <- distinct_effects[i, ]
      col_num <- col_map[[eff$effect_name]]
      lines   <- c(lines,
        sprintf("  Column %d: %-15s [%s, source: %s.%s]",
                col_num, eff$effect_name, eff$effect_class,
                eff$source_table, eff$source_column))
    }
  }
  for (t in trait) {
    col_num   <- col_map[[t]]
    trait_idx <- which(trait == t)
    lines <- c(lines,
      sprintf("  Column %d: %-15s [observation, trait %d]", col_num, t, trait_idx))
  }

  lines <- c(lines, "",
    "Pedigree file: pedigree.txt",
    "  Columns: id_ind id_parent_1 id_parent_2 (unknown parent = 0)",
    "  FILE_POS: 1 2 3 0 0")

  if (!is.null(chip_name)) {
    lines <- c(lines, "",
      "Genotype file: genotype.txt",
      sprintf("  Chip: %s | Loci: %d | ID width: %d chars, then genotype string",
              chip_name, n_loci, id_width))
  }

  eff_labels <- paste0(seq_len(1L + nrow(distinct_effects)),
                       c("=mu", paste0("=", distinct_effects$effect_name)))
  eff_labels <- c(eff_labels, paste0(animal_effect_num, "=animal (EBVs)"))

  lines <- c(lines, "",
    paste0("Animal effect number in solutions file: ", animal_effect_num),
    paste0("  Effect sequence: ", paste(eff_labels, collapse = ", ")),
    "",
    "Solutions file format (OPTION origID, 5 columns):",
    "  trait_num  effect_num  level  original_id  solution",
    paste0("  Filter to effect_num == ", animal_effect_num, " for EBVs"))

  writeLines(lines, file.path(eval_dir, "meta.txt"))
  invisible(NULL)
}


#' Run renumf90 in the eval directory
#' @keywords internal
run_renumf90 <- function(eval_dir, renumf90_path) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(eval_dir)
  result <- system2(renumf90_path,
                    args   = "renum.par",
                    stdout = "renumf90.out",
                    stderr = "renumf90.err")
  if (result != 0)
    stop("renumf90 failed (exit code ", result, "). See: ",
         file.path(eval_dir, "renumf90.out"), call. = FALSE)
  invisible(NULL)
}


#' Run blupf90+ in the eval directory (reads renf90.par written by renumf90)
#' @keywords internal
run_blupf90_plus <- function(eval_dir, blupf90_path) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(eval_dir)
  result <- system2(blupf90_path,
                    args   = "renf90.par",
                    stdout = "blupf90.out",
                    stderr = "blupf90.err")
  if (result != 0)
    stop("blupf90+ failed (exit code ", result, "). See: ",
         file.path(eval_dir, "blupf90.out"), call. = FALSE)
  invisible(NULL)
}


#' Parse the blupf90+ solutions file and return EBV tibble
#'
#' With OPTION origID the solutions file has 5 columns:
#'   trait_num  effect_num  level  original_id  solution
#'
#' @param eval_dir path to evaluation folder
#' @param trait character vector of trait names (in model order)
#' @param animal_effect_num integer effect number for the animal random effect
#' @param all_ped_ids character vector of all pedigree animal IDs (used to
#'   distinguish animal solutions from other random effect solutions)
#' @param model character model label
#' @param date_calc Date
#' @return tibble: id_ind, trait_name, model, ebv, acc, se, date_calc
#' @keywords internal
parse_blupf90_solutions <- function(eval_dir, trait, animal_effect_num,
                                     all_ped_ids, model, date_calc) {
  sols_path <- file.path(eval_dir, "solutions.orig")
  if (!file.exists(sols_path))
    stop("Solutions file not found: ", sols_path, call. = FALSE)

  raw_lines  <- readLines(sols_path)
  data_lines <- raw_lines[grepl("^\\s*[0-9]", raw_lines)]
  if (length(data_lines) == 0)
    stop("Solutions file appears empty or has unexpected format: ", sols_path,
         call. = FALSE)

  sols <- read.table(text   = paste(data_lines, collapse = "\n"),
                     header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  if (ncol(sols) < 5)
    stop("Solutions file has ", ncol(sols), " columns; expected 5 (with OPTION origID).",
         call. = FALSE)
  colnames(sols)[1:5] <- c("trait_num", "effect_num", "level", "orig_id", "solution")
  sols$orig_id <- as.character(sols$orig_id)

  # Filter to animal effect rows that match pedigree IDs
  anim_rows <- sols[sols$effect_num == animal_effect_num &
                    sols$orig_id %in% all_ped_ids, , drop = FALSE]

  if (nrow(anim_rows) == 0) {
    warning("No animal EBVs matched pedigree IDs. ",
            "Check that OPTION origID is active and the solutions file is correct.",
            call. = FALSE)
    return(tibble::tibble(
      id_ind     = character(),
      trait_name = character(),
      model      = character(),
      ebv        = numeric(),
      acc        = numeric(),
      se         = numeric(),
      date_calc  = as.Date(character())
    ))
  }

  anim_rows$trait_name <- trait[as.integer(anim_rows$trait_num)]

  result <- tibble::tibble(
    id_ind     = anim_rows$orig_id,
    trait_name = anim_rows$trait_name,
    model      = model,
    ebv        = as.numeric(anim_rows$solution),
    acc        = NA_real_,
    se         = NA_real_,
    date_calc  = as.Date(date_calc)
  )

  # Summary message
  for (t in unique(result$trait_name)) {
    vals <- result$ebv[result$trait_name == t]
    message(sprintf("  EBV [%s]: n=%d, mean=%.4f, sd=%.4f",
                    t, length(vals), mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE)))
  }
  result
}


#' Stub for VCE writeback — parse blupf90.out and update trait_effect_cov
#' @keywords internal
update_covars_from_blupf90 <- function(pop, eval_dir, trait) {
  message("VCE writeback: automated parsing not yet implemented. ",
          "Inspect blupf90.out and update trait_effect_cov manually via ",
          "add_effect_cov_matrix().")
  invisible(NULL)
}

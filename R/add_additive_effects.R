#' Sample or assign additive QTL effects for a trait
#'
#' @description
#' Selects QTL from a filtered `genome_meta` table and writes additive effects
#' to the `genome_effects` table. This function absorbs the role previously
#' played by `define_qtl()`: the filtered `tidybreed_table` determines which
#' loci become QTL, and the effects are written in a single step.
#'
#' Two modes:
#'
#' * **Manual**: pass `effects`, a numeric vector of length `n_qtl` (number of
#'   filtered loci) in ascending `locus_id` order.
#' * **Sampled**: draw effects from `distribution` (`"normal"` or `"gamma"`).
#'   If `scale_to_target = TRUE`, effects are rescaled using the Falconer
#'   formula so the expected additive variance in the base population equals
#'   the `target_add_var` stored for this trait.
#'
#' The `base` argument controls which allele frequencies are used:
#'
#' * `"founder_haplotypes"` (default) — uses `founder_allele_freq` from
#'   `genome_meta` (requires `initialize_genome()` was called with
#'   `n_haplotypes`).
#' * `"current_pop"` — computes allele frequencies from the current
#'   `genome_haplotype` table. Pass a filtered `tidybreed_table` via
#'   `base_tbl` to restrict which individuals define the base population.
#'
#' Calling this function again for the same `(trait_name, genome_effect_type,
#' line_name)` replaces the existing rows in `genome_effects`.
#'
#' @param tbl A `tidybreed_table` from [get_table()]`("genome_meta")` (with an
#'   optional [filter()]). The filtered rows determine which loci are QTL.
#' @param trait_name Character. Name of an existing trait in `trait_meta`.
#' @param effects Optional numeric vector of length `n_qtl` (manual mode), in
#'   ascending `locus_id` order.
#' @param distribution Character. `"normal"` (default) or `"gamma"`, used when
#'   `effects` is `NULL`.
#' @param base Character. `"founder_haplotypes"` (default) or `"current_pop"`.
#' @param base_tbl Optional `tidybreed_table` (from [get_table()] on any table
#'   with an `id_ind` column) used when `base = "current_pop"` to restrict
#'   which individuals define the base allele frequencies. When `NULL`, all
#'   individuals in `genome_haplotype` are used.
#' @param line_name Optional character. When set, effects are tagged to this
#'   genetic line (for future line-specific TBV). `NULL` (default) means
#'   population-wide effects.
#' @param scale_to_target Logical. If `TRUE`, rescale sampled effects using the
#'   Falconer formula so the expected additive variance equals the stored
#'   `target_add_var`.
#' @param seed Optional integer for reproducibility.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [set_qtl_effects_multi()] for correlated multi-trait effects.
#'
#' @examples
#' \dontrun{
#' # Default: all loci on chr 1-5 become QTL; scale to target variance
#' pop <- pop |>
#'   add_trait("ADG", target_add_var = 0.25, residual_var = 0.75) |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(chr %in% 1:5) |>
#'   add_additive_effects("ADG", distribution = "normal")
#'
#' # current_pop: use generation-0 individuals to define base allele frequencies
#' gen0_tbl <- get_table(pop, "ind_meta") |> dplyr::filter(gen == 0L)
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(chr %in% 1:5) |>
#'   add_additive_effects("ADG", base = "current_pop", base_tbl = gen0_tbl)
#' }
#' @export
add_additive_effects <- function(tbl,
                                  trait_name,
                                  effects         = NULL,
                                  distribution    = c("normal", "gamma"),
                                  base            = c("founder_haplotypes", "current_pop"),
                                  base_tbl        = NULL,
                                  line_name       = NULL,
                                  scale_to_target = TRUE,
                                  seed            = NULL) {

  if (!inherits(tbl, "tidybreed_table")) {
    stop(
      "'tbl' must be a tidybreed_table from get_table('genome_meta') |> filter(...). ",
      "Use: pop |> get_table('genome_meta') |> add_additive_effects()",
      call. = FALSE
    )
  }
  if (tbl$table_name != "genome_meta") {
    stop(
      "'tbl' must be piped from get_table('genome_meta'), not '", tbl$table_name, "'.",
      call. = FALSE
    )
  }

  pop <- tbl$pop
  validate_tidybreed_pop(pop)
  validate_sql_identifier(trait_name, what = "trait name")
  distribution <- match.arg(distribution)
  base         <- match.arg(base)

  if (!is.null(seed)) set.seed(seed)

  if (!is.null(line_name) &&
      (!is.character(line_name) || length(line_name) != 1L)) {
    stop("'line_name' must be a single character string or NULL.", call. = FALSE)
  }

  if (!DBI::dbExistsTable(pop$db_conn, "trait_meta") ||
      nrow(DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT 1 FROM trait_meta WHERE trait_name = '", trait_name, "'")
      )) == 0L) {
    stop("Trait '", trait_name, "' not found in trait_meta.", call. = FALSE)
  }
  target_add_var <- get_effect_var(pop, "gen_add", trait_name)

  # Collect filtered loci; sort by locus_id for consistent effect ordering
  loci_df <- dplyr::collect(tbl)
  if (!"locus_name" %in% names(loci_df)) {
    stop(
      "The filtered table must contain 'locus_name'. ",
      "Pipe get_table('genome_meta') into add_additive_effects().",
      call. = FALSE
    )
  }
  if (nrow(loci_df) == 0L) {
    stop("No loci selected — your filter returned zero rows.", call. = FALSE)
  }
  if ("locus_id" %in% names(loci_df)) {
    loci_df <- loci_df[order(loci_df$locus_id), ]
  }
  selected_locus_names <- loci_df$locus_name
  n_qtl <- length(selected_locus_names)

  # Full-genome locus order for Falconer formula (p_base is in locus_id order)
  genome_order <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id"
  )
  qtl_tf <- genome_order$locus_name %in% selected_locus_names

  # Resolve base_ids from base_tbl when base = "current_pop"
  base_ids <- if (!is.null(base_tbl)) {
    if (!inherits(base_tbl, "tidybreed_table")) {
      stop("'base_tbl' must be a tidybreed_table.", call. = FALSE)
    }
    b <- dplyr::collect(base_tbl$tbl)
    if (!"id_ind" %in% names(b)) {
      stop("'base_tbl' must contain an 'id_ind' column.", call. = FALSE)
    }
    unique(b[["id_ind"]])
  } else {
    NULL
  }

  p_base <- compute_base_allele_freq(pop, base, base_ids)

  if (!is.null(effects)) {
    if (!is.numeric(effects)) stop("`effects` must be numeric.", call. = FALSE)
    if (length(effects) != n_qtl) {
      stop("`effects` length (", length(effects), ") must equal number of selected loci (",
           n_qtl, ").", call. = FALSE)
    }
    qtl_effects <- as.numeric(effects)
  } else {
    qtl_effects <- switch(
      distribution,
      normal = stats::rnorm(n_qtl),
      gamma  = stats::rgamma(n_qtl, shape = 0.4, rate = 1.66) *
                 sample(c(-1, 1), n_qtl, replace = TRUE)
    )
    if (scale_to_target) {
      if (is.na(target_add_var)) {
        stop(
          "No additive genetic variance stored for trait '", trait_name, "'. ",
          "Call add_effect_cov_matrix(pop, 'gen_add', ...) or ",
          "add_trait(pop, '", trait_name, "', target_add_var = ...) first.",
          call. = FALSE
        )
      }
      qtl_effects <- rescale_effects_to_target(qtl_tf, qtl_effects, target_add_var, p_base)
    }
  }

  # Base allele frequencies at the selected loci (in locus_id order)
  p_base_qtl <- as.numeric(p_base[qtl_tf])

  # Replace any existing rows for this (trait, effect type, line) combination
  line_sql <- if (is.null(line_name)) {
    "line_name IS NULL"
  } else {
    paste0("line_name = '", line_name, "'")
  }
  DBI::dbExecute(
    pop$db_conn,
    paste0(
      "DELETE FROM genome_effects ",
      "WHERE trait_name = '", trait_name, "' ",
      "AND genome_effect_type = 'additive' ",
      "AND ", line_sql
    )
  )

  start_id <- next_int_id(pop$db_conn, "genome_effects", "id_genome_effect")
  effects_df <- tibble::tibble(
    id_genome_effect   = seq.int(start_id, length.out = n_qtl),
    locus_name         = selected_locus_names,
    line_name          = line_name,
    trait_name         = trait_name,
    genome_effect_type = "additive",
    genome_value       = as.numeric(qtl_effects),
    base_allele_freq   = p_base_qtl
  )
  DBI::dbWriteTable(pop$db_conn, "genome_effects", effects_df, append = TRUE)

  message("Set additive effects for ", n_qtl, " QTL on trait '", trait_name,
          "' (base: ", base, ").")
  invisible(pop)
}


#' Sample correlated additive QTL effects across multiple traits
#'
#' @description
#' Draws additive effects for multiple traits from a multivariate normal
#' distribution keyed by the additive-genetic covariance matrix `G`, and writes
#' them to `genome_effects`. Two selection methods:
#'
#' * `method = "shared"` — effects are drawn at the loci supplied via `tbl`
#'   (used as the shared QTL set for all traits). Loci that are QTL for only a
#'   subset of traits in `genome_effects` also receive independent draws (with
#'   the diagonal variance of `G` for that trait).
#' * `method = "union"` — the loci in `tbl` form the union of QTL across all
#'   traits. Existing per-trait QTL sets (from prior `add_additive_effects()`
#'   calls) in `genome_effects` are used to determine which loci are active for
#'   each trait within that union; loci not yet assigned to a trait receive `NA`
#'   and are excluded from the TBV for that trait.
#'
#' If `scale_to_target = TRUE`, each trait's effects are rescaled independently
#' using the Falconer formula.
#'
#' @param tbl A `tidybreed_table` from [get_table()]`("genome_meta")` (with an
#'   optional [filter()]) that defines the loci to use. For `method = "shared"`
#'   these loci become the QTL for **all** `trait_names`. For `method = "union"`
#'   these loci form the candidate pool; per-trait membership is determined from
#'   existing rows in `genome_effects`.
#' @param trait_names Character vector of trait names (length >= 2). All must
#'   exist in `trait_meta`.
#' @param G Optional numeric matrix of additive-genetic (co)variances. Must be
#'   square and symmetric with side length `length(trait_names)`. When supplied,
#'   stored to `trait_effect_cov` under `"gen_add"`. When `NULL`, read from
#'   `trait_effect_cov`.
#' @param method Character. `"shared"` (default) or `"union"`.
#' @param base Character. `"founder_haplotypes"` (default) or `"current_pop"`.
#' @param base_tbl Optional `tidybreed_table` used when `base = "current_pop"`
#'   to restrict which individuals define base allele frequencies.
#' @param line_name Optional character. Tag effects to a specific line (NULL =
#'   population-wide).
#' @param scale_to_target Logical. Rescale each trait's effects to its
#'   `target_add_var` using the Falconer formula.
#' @param seed Optional integer for reproducibility.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @examples
#' \dontrun{
#' G <- matrix(c(0.25, 0.10, 0.10, 0.30), 2, 2,
#'             dimnames = list(c("ADG", "BW"), c("ADG", "BW")))
#' pop <- pop |>
#'   add_effect_cov_matrix("gen_add", G) |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(chr %in% 1:5) |>
#'   set_qtl_effects_multi(trait_names = c("ADG", "BW"), G = G)
#' }
#' @export
set_qtl_effects_multi <- function(tbl,
                                  trait_names,
                                  G               = NULL,
                                  method          = c("shared", "union"),
                                  base            = c("founder_haplotypes", "current_pop"),
                                  base_tbl        = NULL,
                                  line_name       = NULL,
                                  scale_to_target = TRUE,
                                  seed            = NULL) {

  if (!inherits(tbl, "tidybreed_table")) {
    stop(
      "'tbl' must be a tidybreed_table from get_table('genome_meta') |> filter(...). ",
      "Use: pop |> get_table('genome_meta') |> set_qtl_effects_multi()",
      call. = FALSE
    )
  }
  if (tbl$table_name != "genome_meta") {
    stop(
      "'tbl' must be piped from get_table('genome_meta'), not '", tbl$table_name, "'.",
      call. = FALSE
    )
  }

  pop <- tbl$pop
  validate_tidybreed_pop(pop)
  base   <- match.arg(base)
  method <- match.arg(method)
  stopifnot(is.character(trait_names), length(trait_names) >= 2)
  lapply(trait_names, validate_sql_identifier, what = "trait name")

  if (!is.null(line_name) &&
      (!is.character(line_name) || length(line_name) != 1L)) {
    stop("'line_name' must be a single character string or NULL.", call. = FALSE)
  }

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for set_qtl_effects_multi(). ",
         "Install with install.packages('MASS').", call. = FALSE)
  }

  # Resolve G: if not supplied, read from trait_effect_cov
  if (!is.null(G)) {
    if (!is.matrix(G) || nrow(G) != length(trait_names) || ncol(G) != length(trait_names)) {
      stop("`G` must be a square matrix with side = length(trait_names).", call. = FALSE)
    }
    if (!isSymmetric(unname(G))) stop("`G` must be symmetric.", call. = FALSE)
    g_named <- G
    dimnames(g_named) <- list(trait_names, trait_names)
    pop <- add_effect_cov_matrix(pop, "gen_add", g_named)
  } else {
    G_stored <- load_effect_cov(pop, "gen_add", trait_names)
    if (is.null(G_stored)) {
      stop("No 'gen_add' covariance matrix found for traits: ",
           paste(trait_names, collapse = ", "),
           ". Call add_effect_cov_matrix(pop, 'gen_add', G) first or pass G directly.",
           call. = FALSE)
    }
    G <- G_stored
  }

  if (!is.null(seed)) set.seed(seed)

  # Validate traits exist
  trait_meta_rows <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT trait_name FROM trait_meta WHERE trait_name IN (",
           paste0("'", trait_names, "'", collapse = ", "), ")")
  )
  missing_traits <- setdiff(trait_names, trait_meta_rows$trait_name)
  if (length(missing_traits) > 0) {
    stop("Traits not found in trait_meta: ", paste(missing_traits, collapse = ", "),
         call. = FALSE)
  }

  target_var <- stats::setNames(
    vapply(trait_names, function(t) get_effect_var(pop, "gen_add", t), numeric(1)),
    trait_names
  )

  # Collect candidate loci from tbl (sorted by locus_id)
  loci_df <- dplyr::collect(tbl)
  if (!"locus_name" %in% names(loci_df)) {
    stop("The filtered table must contain 'locus_name'.", call. = FALSE)
  }
  if (nrow(loci_df) == 0L) stop("No loci selected — filter returned zero rows.", call. = FALSE)
  if ("locus_id" %in% names(loci_df)) loci_df <- loci_df[order(loci_df$locus_id), ]
  candidate_locus_names <- loci_df$locus_name

  # Full-genome order (for Falconer formula)
  genome_order <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id"
  )
  n_loci <- nrow(genome_order)

  # Build per-trait QTL membership mask (length n_loci, in locus_id order)
  qtl_tf_mat <- matrix(FALSE, nrow = n_loci, ncol = length(trait_names),
                       dimnames = list(NULL, trait_names))

  line_sql <- if (is.null(line_name)) "line_name IS NULL" else
    paste0("line_name = '", line_name, "'")

  for (t in trait_names) {
    if (method == "shared") {
      qtl_tf_mat[, t] <- genome_order$locus_name %in% candidate_locus_names
    } else {
      # union: per-trait QTL from existing genome_effects rows, restricted to candidate set
      existing <- DBI::dbGetQuery(
        pop$db_conn,
        paste0("SELECT locus_name FROM genome_effects ",
               "WHERE trait_name = '", t, "' AND genome_effect_type = 'additive' ",
               "AND ", line_sql)
      )$locus_name
      active <- intersect(existing, candidate_locus_names)
      if (length(active) == 0) {
        warning("Trait '", t, "' has no existing additive effects in genome_effects ",
                "within the candidate loci; it will receive no effects from this call.",
                call. = FALSE)
      }
      qtl_tf_mat[, t] <- genome_order$locus_name %in% active
    }
  }

  effects_mat <- matrix(NA_real_, nrow = n_loci, ncol = length(trait_names),
                        dimnames = list(NULL, trait_names))

  if (method == "shared") {
    shared <- apply(qtl_tf_mat, 1, all)
    n_shared <- sum(shared)
    if (n_shared == 0) {
      warning("No loci are QTL for all traits; using union fallback.", call. = FALSE)
    } else {
      draws <- MASS::mvrnorm(n = n_shared, mu = rep(0, length(trait_names)), Sigma = G)
      if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1)
      effects_mat[shared, ] <- draws
    }
    # Independent draws for trait-specific (non-shared) QTL
    for (k in seq_along(trait_names)) {
      t    <- trait_names[k]
      solo <- qtl_tf_mat[, t] & !shared
      if (sum(solo) > 0) {
        effects_mat[solo, t] <- stats::rnorm(sum(solo), sd = sqrt(G[k, k]))
      }
    }
  } else {  # union
    any_qtl <- apply(qtl_tf_mat, 1, any)
    n_any   <- sum(any_qtl)
    if (n_any == 0) stop("No QTL found across any trait in the candidate loci.", call. = FALSE)
    draws <- MASS::mvrnorm(n = n_any, mu = rep(0, length(trait_names)), Sigma = G)
    if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1)
    block_mask <- qtl_tf_mat[any_qtl, , drop = FALSE]
    draws[!block_mask] <- NA_real_
    effects_mat[any_qtl, ] <- draws
  }

  # Resolve base_ids from base_tbl
  base_ids <- if (!is.null(base_tbl)) {
    if (!inherits(base_tbl, "tidybreed_table")) stop("'base_tbl' must be a tidybreed_table.", call. = FALSE)
    b <- dplyr::collect(base_tbl$tbl)
    if (!"id_ind" %in% names(b)) stop("'base_tbl' must contain 'id_ind'.", call. = FALSE)
    unique(b[["id_ind"]])
  } else {
    NULL
  }

  p_base <- compute_base_allele_freq(pop, base, base_ids)

  if (scale_to_target) {
    na_targets <- trait_names[is.na(target_var)]
    if (length(na_targets) > 0) {
      stop("No additive genetic variance stored for trait(s): ",
           paste(na_targets, collapse = ", "),
           ". Supply `G` or call add_effect_cov_matrix() first.", call. = FALSE)
    }
    for (k in seq_along(trait_names)) {
      t        <- trait_names[k]
      qtl_tf_k <- qtl_tf_mat[, t]
      a_qtl    <- effects_mat[qtl_tf_k, t]
      if (any(!is.na(a_qtl))) {
        effects_mat[qtl_tf_k, t] <- rescale_effects_to_target(
          qtl_tf_k, a_qtl, target_var[[t]], p_base
        )
      }
    }
  }

  # Write to genome_effects: one row per (locus, trait) that has an effect
  for (t in trait_names) {
    qtl_tf_t        <- qtl_tf_mat[, t]
    locus_names_t   <- genome_order$locus_name[qtl_tf_t]
    effects_t       <- effects_mat[qtl_tf_t, t]
    non_na          <- !is.na(effects_t)
    locus_names_t   <- locus_names_t[non_na]
    effects_t       <- effects_t[non_na]
    p_base_qtl_t    <- as.numeric(p_base[qtl_tf_t][non_na])
    n_t             <- length(locus_names_t)
    if (n_t == 0L) next

    # Delete existing rows for this trait+type+line
    DBI::dbExecute(
      pop$db_conn,
      paste0(
        "DELETE FROM genome_effects ",
        "WHERE trait_name = '", t, "' ",
        "AND genome_effect_type = 'additive' ",
        "AND ", line_sql
      )
    )

    start_id <- next_int_id(pop$db_conn, "genome_effects", "id_genome_effect")
    df_t <- tibble::tibble(
      id_genome_effect   = seq.int(start_id, length.out = n_t),
      locus_name         = locus_names_t,
      line_name          = line_name,
      trait_name         = t,
      genome_effect_type = "additive",
      genome_value       = as.numeric(effects_t),
      base_allele_freq   = p_base_qtl_t
    )
    DBI::dbWriteTable(pop$db_conn, "genome_effects", df_t, append = TRUE)
  }

  message("Set correlated additive effects for traits: ",
          paste(trait_names, collapse = ", "), " (method: ", method, ")")
  invisible(pop)
}


#' Rescale QTL effects to hit a target additive variance using Falconer formula
#'
#' @param qtl_tf Logical mask of QTL loci, length `n_loci`.
#' @param qtl_effects Numeric effects at QTL loci, length `sum(qtl_tf)`.
#' @param target_add_var Target additive variance.
#' @param p_base Numeric vector of base allele frequencies, length `n_loci`.
#' @return Rescaled `qtl_effects` vector.
#' @keywords internal
rescale_effects_to_target <- function(qtl_tf, qtl_effects, target_add_var, p_base) {
  p_qtl <- p_base[qtl_tf]
  V_A   <- sum(2 * p_qtl * (1 - p_qtl) * qtl_effects^2)
  if (V_A <= 0 || !is.finite(V_A)) {
    warning("Falconer V_A is zero or infinite; effects returned unchanged.", call. = FALSE)
    return(qtl_effects)
  }
  qtl_effects * sqrt(target_add_var / V_A)
}


#' Compute per-locus allele frequencies from the base population
#'
#' @param pop A `tidybreed_pop` object.
#' @param base Character. `"founder_haplotypes"` or `"current_pop"`.
#' @param base_ids Optional character vector of `id_ind` for `"current_pop"`.
#' @return Numeric vector of allele frequencies, length `n_loci`, in
#'   `locus_id` order.
#' @keywords internal
compute_base_allele_freq <- function(pop, base, base_ids = NULL) {
  if (base == "founder_haplotypes") {
    if (!"founder_haplotypes" %in% pop$tables) {
      stop(
        "founder_haplotypes table not found. ",
        "Did you call initialize_genome() with n_haplotypes? ",
        "Use base = 'current_pop' instead.",
        call. = FALSE
      )
    }
    df <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM founder_haplotypes")
    if (nrow(df) == 0) {
      stop("founder_haplotypes table is empty.", call. = FALSE)
    }
    locus_cols  <- setdiff(names(df), "hap_id")
    locus_order <- order(as.integer(sub("^locus_", "", locus_cols)))
    hap_mat     <- as.matrix(df[, locus_cols[locus_order], drop = FALSE])
    return(colMeans(hap_mat))
  }

  # current_pop: compute from genome_haplotype rows
  if (!"genome_haplotype" %in% pop$tables) {
    stop("genome_haplotype table does not exist.", call. = FALSE)
  }
  if (is.null(base_ids)) {
    df <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_haplotype")
  } else {
    if (length(base_ids) == 0) stop("base_ids is empty.", call. = FALSE)
    ids_sql <- paste0("'", base_ids, "'", collapse = ", ")
    df <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT * FROM genome_haplotype WHERE id_ind IN (", ids_sql, ")")
    )
  }
  if (nrow(df) == 0) {
    stop("No haplotype rows found for the base population.", call. = FALSE)
  }
  locus_cols  <- setdiff(names(df), c("id_ind", "parent_origin"))
  locus_order <- order(as.integer(sub("^locus_", "", locus_cols)))
  hap_mat     <- as.matrix(df[, locus_cols[locus_order], drop = FALSE])
  colMeans(hap_mat)
}


#' Pull the genotype matrix (individuals x loci) into memory
#'
#' @param pop A `tidybreed_pop` object.
#' @param subset_ids Optional character vector of `id_ind`. When supplied,
#'   only those rows are pulled from DuckDB (push-down filter) — important
#'   for large genomes where the full genotype matrix wouldn't fit in R.
#' @return Numeric matrix with rows in `id_ind` order and columns in
#'   `locus_id` order.
#' @keywords internal
get_genotype_matrix <- function(pop, subset_ids = NULL) {
  if (!"genome_genotype" %in% pop$tables) {
    stop("genome_genotype table does not exist. Cannot rescale effects.",
         call. = FALSE)
  }
  if (is.null(subset_ids)) {
    df <- DBI::dbGetQuery(pop$db_conn, "SELECT * FROM genome_genotype")
  } else {
    if (length(subset_ids) == 0) {
      stop("subset_ids is empty.", call. = FALSE)
    }
    ids_sql <- paste0("'", subset_ids, "'", collapse = ", ")
    df <- DBI::dbGetQuery(
      pop$db_conn,
      paste0("SELECT * FROM genome_genotype WHERE id_ind IN (", ids_sql, ")")
    )
  }
  if (nrow(df) == 0) {
    stop("No genotype rows available for the requested subset.",
         call. = FALSE)
  }
  locus_cols <- setdiff(names(df), "id_ind")
  locus_order <- order(as.integer(sub("^locus_", "", locus_cols)))
  mat <- as.matrix(df[, locus_cols[locus_order], drop = FALSE])
  rownames(mat) <- df$id_ind
  mat
}

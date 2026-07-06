#' Define additive QTL effects for one or more traits
#'
#' @description
#' Selects QTL from a filtered `genome_meta` table and writes additive effects
#' to the `genome_effects` table.
#'
#' **Single trait** (`trait_name` length 1) — two modes:
#'
#' * **Manual**: pass `effects`, a numeric vector of length `n_qtl` (number of
#'   filtered loci) in ascending `locus_id` order.
#' * **Sampled**: draw effects from `distribution` (`"normal"` or `"gamma"`).
#'   If `scale_to_target = TRUE`, effects are rescaled using the Falconer
#'   formula so the expected additive variance in the base population equals
#'   the `target_add_var` stored for this trait.
#'
#' **Multiple traits** (`trait_name` length >= 2) — effects are drawn jointly
#' from a multivariate normal distribution keyed by the additive-genetic
#' covariance matrix `G`. Two locus-selection methods:
#'
#' * `method = "shared"` — the loci in `tbl` become the shared QTL set for
#'   all traits. Loci that are QTL for only a subset of traits in
#'   `genome_effects` also receive independent draws (with the diagonal
#'   variance of `G` for that trait).
#' * `method = "union"` — the loci in `tbl` form the candidate pool; per-trait
#'   membership is determined from existing rows in `genome_effects`.
#'
#' The `base` argument controls which allele frequencies are used:
#'
#' * `"founder_haplotypes"` (default) — uses `founder_allele_freq` from
#'   `genome_meta` (requires `initialize_genome()` was called with
#'   `n_haplotypes`).
#' * `"current_pop"` — computes allele frequencies from the current
#'   `ind_haplotype` table. Pass a filtered `tidybreed_table` via
#'   `base_tbl` to restrict which individuals define the base population.
#'
#' Calling this function again for the same `(trait_name, genome_effect_type,
#' line_name)` replaces the existing rows in `genome_effects`.
#'
#' @param tbl A `tidybreed_table` from [get_table()]`("genome_meta")` (with an
#'   optional [dplyr::filter()]). The filtered rows determine which loci are QTL.
#' @param trait_name Character scalar **or** vector. Name(s) of existing traits
#'   in `trait_meta`. When length >= 2, effects are drawn jointly from
#'   `MVN(0, G)` and `G` / `method` become active.
#' @param effects Optional numeric vector of length `n_qtl` (manual mode, single
#'   trait only), in ascending `locus_id` order. Error if `length(trait_name) > 1`.
#' @param distribution Character. `"normal"` (default) or `"gamma"`, used when
#'   `effects` is `NULL` and `length(trait_name) == 1`. Ignored for multi-trait.
#' @param G Optional numeric matrix of additive-genetic (co)variances (multi-trait
#'   only). Must be square and symmetric with side length `length(trait_name)`.
#'   When supplied, stored to `trait_var_comp` under `"gen_add"`. When `NULL`,
#'   read from `trait_var_comp`.
#' @param method Character. `"shared"` (default) or `"union"`. Multi-trait only.
#'   `"shared"` — all listed traits use the filtered loci as their shared QTL
#'   set. `"union"` — per-trait QTL sets are read from existing `genome_effects`
#'   rows, restricted to the filtered loci.
#' @param base Character. `"founder_haplotypes"` (default) or `"current_pop"`.
#' @param base_tbl Optional `tidybreed_table` (from [get_table()] on any table
#'   with an `id_ind` column) used when `base = "current_pop"` to restrict
#'   which individuals define the base allele frequencies. When `NULL`, all
#'   individuals in `ind_haplotype` are used.
#' @param line_name Optional character. When set, effects are tagged to this
#'   genetic line (for future line-specific TBV). `NULL` (default) means
#'   population-wide effects.
#' @param scale_to_target Logical. If `TRUE`, rescale effects using the Falconer
#'   formula so the expected additive variance equals the stored `target_add_var`.
#' @param seed Optional integer for reproducibility.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_trait()], [define_effect_cov_matrix()].
#'
#' @examples
#' \dontrun{
#' # Single trait — all loci on chr 1-5 become QTL; scale to target variance
#' pop <- pop |>
#'   define_trait("ADG", target_add_var = 0.25) |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(chr %in% 1:5) |>
#'   define_additive_effects("ADG", distribution = "normal")
#'
#' # Multiple correlated traits — shared QTL set, joint MVN draw
#' G <- matrix(c(0.25, 0.10, 0.10, 0.30), 2, 2,
#'             dimnames = list(c("ADG", "BW"), c("ADG", "BW")))
#' pop <- pop |>
#'   define_effect_cov_matrix("gen_add", G) |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(chr %in% 1:5) |>
#'   define_additive_effects(c("ADG", "BW"), G = G)
#'
#' # current_pop: use generation-0 individuals to define base allele frequencies
#' gen0_tbl <- get_table(pop, "ind_meta") |> dplyr::filter(gen == 0L)
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(chr %in% 1:5) |>
#'   define_additive_effects("ADG", base = "current_pop", base_tbl = gen0_tbl)
#' }
#' @export
define_additive_effects <- function(tbl,
                                    trait_name,
                                    effects         = NULL,
                                    distribution    = c("normal", "gamma"),
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
      "Use: pop |> get_table('genome_meta') |> define_additive_effects()",
      call. = FALSE
    )
  }
  if (tbl$table_name != "genome_meta") {
    stop(
      "'tbl' must be piped from get_table('genome_meta'), not '", tbl$table_name, "'.",
      call. = FALSE
    )
  }

  stopifnot(is.character(trait_name), length(trait_name) >= 1L)
  pop          <- tbl$pop
  validate_tidybreed_pop(pop)
  distribution <- match.arg(distribution)
  base         <- match.arg(base)
  method       <- match.arg(method)

  if (!is.null(line_name) &&
      (!is.character(line_name) || length(line_name) != 1L)) {
    stop("'line_name' must be a single character string or NULL.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  # ------------------------------------------------------------------ #
  #  Single-trait path                                                   #
  # ------------------------------------------------------------------ #
  if (length(trait_name) == 1L) {

    validate_sql_identifier(trait_name, what = "trait name")

    if (!DBI::dbExistsTable(pop$db_conn, "trait_meta") ||
        nrow(DBI::dbGetQuery(
          pop$db_conn,
          paste0("SELECT 1 FROM trait_meta WHERE trait_name = '", trait_name, "'")
        )) == 0L) {
      stop("Trait '", trait_name, "' not found in trait_meta.", call. = FALSE)
    }
    target_add_var <- get_trait_var(pop, "gen_add", trait_name)

    loci_df <- dplyr::collect(tbl)
    if (!"locus_name" %in% names(loci_df)) {
      stop(
        "The filtered table must contain 'locus_name'. ",
        "Pipe get_table('genome_meta') into define_additive_effects().",
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

    genome_order <- DBI::dbGetQuery(
      pop$db_conn,
      "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id"
    )
    qtl_tf <- genome_order$locus_name %in% selected_locus_names

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
            "Call define_effect_cov_matrix(pop, 'gen_add', ...) or ",
            "define_trait(pop, '", trait_name, "', target_add_var = ...) first.",
            call. = FALSE
          )
        }
        assert_qtl_autosomal(pop$db_conn, selected_locus_names)
        qtl_effects <- rescale_effects_to_target(qtl_tf, qtl_effects, target_add_var, p_base)
      }
    }

    p_base_qtl <- as.numeric(p_base[qtl_tf])

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
    return(invisible(pop))
  }

  # ------------------------------------------------------------------ #
  #  Multi-trait path (length(trait_name) >= 2)                         #
  # ------------------------------------------------------------------ #
  if (!is.null(effects)) {
    stop("'effects' cannot be used when 'trait_name' has length > 1.", call. = FALSE)
  }
  if (distribution != "normal") {
    warning("'distribution' is ignored for multi-trait; effects are always drawn from MVN.",
            call. = FALSE)
  }

  lapply(trait_name, validate_sql_identifier, what = "trait name")

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for multi-trait effect sampling. ",
         "Install with install.packages('MASS').", call. = FALSE)
  }

  # Resolve G: if not supplied, read from trait_var_comp
  if (!is.null(G)) {
    if (!is.matrix(G) || nrow(G) != length(trait_name) || ncol(G) != length(trait_name)) {
      stop("`G` must be a square matrix with side = length(trait_name).", call. = FALSE)
    }
    if (!isSymmetric(unname(G))) stop("`G` must be symmetric.", call. = FALSE)
    g_named <- G
    dimnames(g_named) <- list(trait_name, trait_name)
    pop <- define_effect_cov_matrix(pop, "gen_add", g_named)
  } else {
    G_stored <- load_trait_cov(pop, "gen_add", trait_name)
    if (is.null(G_stored)) {
      stop("No 'gen_add' covariance matrix found for traits: ",
           paste(trait_name, collapse = ", "),
           ". Call define_effect_cov_matrix(pop, 'gen_add', G) first or pass G directly.",
           call. = FALSE)
    }
    G <- G_stored
  }

  # Validate traits exist
  trait_meta_rows <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT trait_name FROM trait_meta WHERE trait_name IN (",
           paste0("'", trait_name, "'", collapse = ", "), ")")
  )
  missing_traits <- setdiff(trait_name, trait_meta_rows$trait_name)
  if (length(missing_traits) > 0) {
    stop("Traits not found in trait_meta: ", paste(missing_traits, collapse = ", "),
         call. = FALSE)
  }

  target_var <- stats::setNames(
    vapply(trait_name, function(t) get_trait_var(pop, "gen_add", t), numeric(1)),
    trait_name
  )

  loci_df <- dplyr::collect(tbl)
  if (!"locus_name" %in% names(loci_df)) {
    stop("The filtered table must contain 'locus_name'.", call. = FALSE)
  }
  if (nrow(loci_df) == 0L) stop("No loci selected — filter returned zero rows.", call. = FALSE)
  if ("locus_id" %in% names(loci_df)) loci_df <- loci_df[order(loci_df$locus_id), ]
  candidate_locus_names <- loci_df$locus_name

  genome_order <- DBI::dbGetQuery(
    pop$db_conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id"
  )
  n_loci <- nrow(genome_order)

  qtl_tf_mat <- matrix(FALSE, nrow = n_loci, ncol = length(trait_name),
                       dimnames = list(NULL, trait_name))

  line_sql <- if (is.null(line_name)) "line_name IS NULL" else
    paste0("line_name = '", line_name, "'")

  for (t in trait_name) {
    if (method == "shared") {
      qtl_tf_mat[, t] <- genome_order$locus_name %in% candidate_locus_names
    } else {
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

  effects_mat <- matrix(NA_real_, nrow = n_loci, ncol = length(trait_name),
                        dimnames = list(NULL, trait_name))

  if (method == "shared") {
    shared   <- apply(qtl_tf_mat, 1, all)
    n_shared <- sum(shared)
    if (n_shared == 0) {
      warning("No loci are QTL for all traits; using union fallback.", call. = FALSE)
    } else {
      draws <- MASS::mvrnorm(n = n_shared, mu = rep(0, length(trait_name)), Sigma = G)
      if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1)
      effects_mat[shared, ] <- draws
    }
    for (k in seq_along(trait_name)) {
      t    <- trait_name[k]
      solo <- qtl_tf_mat[, t] & !shared
      if (sum(solo) > 0) {
        effects_mat[solo, t] <- stats::rnorm(sum(solo), sd = sqrt(G[k, k]))
      }
    }
  } else {
    any_qtl <- apply(qtl_tf_mat, 1, any)
    n_any   <- sum(any_qtl)
    if (n_any == 0) stop("No QTL found across any trait in the candidate loci.", call. = FALSE)
    draws <- MASS::mvrnorm(n = n_any, mu = rep(0, length(trait_name)), Sigma = G)
    if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1)
    block_mask <- qtl_tf_mat[any_qtl, , drop = FALSE]
    draws[!block_mask] <- NA_real_
    effects_mat[any_qtl, ] <- draws
  }

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
    na_targets <- trait_name[is.na(target_var)]
    if (length(na_targets) > 0) {
      stop("No additive genetic variance stored for trait(s): ",
           paste(na_targets, collapse = ", "),
           ". Supply `G` or call define_effect_cov_matrix() first.", call. = FALSE)
    }
    for (k in seq_along(trait_name)) {
      t        <- trait_name[k]
      qtl_tf_k <- qtl_tf_mat[, t]
      a_qtl    <- effects_mat[qtl_tf_k, t]
      if (any(!is.na(a_qtl))) {
        assert_qtl_autosomal(pop$db_conn, genome_order$locus_name[qtl_tf_k])
        effects_mat[qtl_tf_k, t] <- rescale_effects_to_target(
          qtl_tf_k, a_qtl, target_var[[t]], p_base
        )
      }
    }
  }

  for (t in trait_name) {
    qtl_tf_t      <- qtl_tf_mat[, t]
    locus_names_t <- genome_order$locus_name[qtl_tf_t]
    effects_t     <- effects_mat[qtl_tf_t, t]
    non_na        <- !is.na(effects_t)
    locus_names_t <- locus_names_t[non_na]
    effects_t     <- effects_t[non_na]
    p_base_qtl_t  <- as.numeric(p_base[qtl_tf_t][non_na])
    n_t           <- length(locus_names_t)
    if (n_t == 0L) next

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
          paste(trait_name, collapse = ", "), " (method: ", method, ")")
  invisible(pop)
}


#' Sample correlated additive QTL effects across multiple traits
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated. Use [define_additive_effects()] with a character
#' vector for `trait_name` instead:
#'
#' ```r
#' # old
#' tbl |> set_qtl_effects_multi(trait_names = c("ADG", "BW"), G = G)
#'
#' # new
#' tbl |> define_additive_effects(trait_name = c("ADG", "BW"), G = G)
#' ```
#'
#' @param tbl Passed to [define_additive_effects()].
#' @param trait_names Passed to [define_additive_effects()] as `trait_name`.
#' @param G Passed to [define_additive_effects()].
#' @param method Passed to [define_additive_effects()].
#' @param base Passed to [define_additive_effects()].
#' @param base_tbl Passed to [define_additive_effects()].
#' @param line_name Passed to [define_additive_effects()].
#' @param scale_to_target Passed to [define_additive_effects()].
#' @param seed Passed to [define_additive_effects()].
#'
#' @return The modified `tidybreed_pop` (invisibly).
#' @keywords internal
set_qtl_effects_multi <- function(tbl,
                                  trait_names,
                                  G               = NULL,
                                  method          = c("shared", "union"),
                                  base            = c("founder_haplotypes", "current_pop"),
                                  base_tbl        = NULL,
                                  line_name       = NULL,
                                  scale_to_target = TRUE,
                                  seed            = NULL) {
  .Deprecated(
    "define_additive_effects",
    msg = paste0(
      "'set_qtl_effects_multi()' is deprecated. ",
      "Use: tbl |> define_additive_effects(trait_name = c(...), G = G)"
    )
  )
  define_additive_effects(
    tbl             = tbl,
    trait_name      = trait_names,
    G               = G,
    method          = method,
    base            = base,
    base_tbl        = base_tbl,
    line_name       = line_name,
    scale_to_target = scale_to_target,
    seed            = seed
  )
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
        "Did you call define_founder_haplotypes()? ",
        "Use base = 'current_pop' instead.",
        call. = FALSE
      )
    }
    freq <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT gm.locus_id, AVG(CAST(fh.allele AS DOUBLE)) AS f ",
      "FROM founder_haplotypes fh ",
      "JOIN genome_meta gm ON fh.locus_name = gm.locus_name ",
      "GROUP BY gm.locus_id ORDER BY gm.locus_id"))
    if (nrow(freq) == 0) {
      stop("founder_haplotypes table is empty.", call. = FALSE)
    }
    n_loci <- DBI::dbGetQuery(pop$db_conn,
      "SELECT COUNT(*) AS n FROM genome_meta")$n
    out <- numeric(n_loci)
    out[freq$locus_id] <- freq$f
    return(out)
  }

  if (!"ind_haplotype" %in% pop$tables) {
    stop("ind_haplotype table does not exist.", call. = FALSE)
  }
  where <- ""
  if (!is.null(base_ids)) {
    if (length(base_ids) == 0) stop("base_ids is empty.", call. = FALSE)
    ids_sql <- paste0("'", base_ids, "'", collapse = ", ")
    where   <- paste0("WHERE id_ind IN (", ids_sql, ") ")
  }
  # Per-locus base allele frequency = mean allele over all haplotype rows.
  freq <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT locus_id, AVG(CAST(allele AS DOUBLE)) AS f ",
           "FROM ind_haplotype ", where, "GROUP BY locus_id ORDER BY locus_id")
  )
  if (nrow(freq) == 0) {
    stop("No haplotype rows found for the base population.", call. = FALSE)
  }
  n_loci <- DBI::dbGetQuery(pop$db_conn,
    "SELECT COUNT(*) AS n FROM genome_meta")$n
  out <- numeric(n_loci)
  out[freq$locus_id] <- freq$f
  out
}


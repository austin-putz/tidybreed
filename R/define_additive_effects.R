#' Define additive QTL effects for one or more traits
#'
#' @description
#' Selects QTL from a filtered `genome_meta` table and writes one order-one
#' `additive` term per locus through [define_genome_effects()], under the
#' reserved effect owner `"generated_additive_tbv"`. [add_tbv()] reads order-one
#' `additive` variants from that owner and nothing else, so effects written here
#' and effects a user writes with [define_genome_effects()] can never be
#' confused for one another.
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
#'   membership is read from the terms already stored at this scope.
#'
#' The `base` argument controls which allele frequencies are used:
#'
#' * `"founder_haplotypes"` (default) — computes allele frequencies directly
#'   from the `founder_haplotypes` table (requires [define_founder_haplotypes()]
#'   was called). Restrict to one founder pool with `base_line_name`. (This does
#'   **not** read `genome_meta.founder_allele_freq`, which is informational only.)
#' * `"current_pop"` — computes allele frequencies from the current
#'   `ind_haplotype` table. Pass a filtered `tidybreed_table` via
#'   `base_tbl` to restrict which individuals define the base population.
#'
#' @section Which population centers the effects:
#' Base allele frequencies center the true breeding value (the Falconer
#' `allele - p` term) and set the `2pq` denominator used by `scale_to_target`.
#' By default they come from **the population the effect applies to**:
#' `base_line_name` inherits `line_name`, so a line-specific effect is centered
#' on that line's own founder pool and a population-wide effect
#' (`line_name = NULL`) on the whole founder base.
#'
#' This matters because pooling divergent lines overstates within-line
#' heterozygosity — the Wahlund effect. Two lines fixed for opposite alleles
#' each have zero within-line variance, but pool to `p = 0.5` and an apparent
#' `2pq = 0.5`; the inflated denominator then makes `scale_to_target` **under**-scale
#' the effects, and realized within-line additive variance falls short of
#' `target_add_var`. Pass `base_line_name = NULL` explicitly to force pooling
#' anyway.
#'
#' The centering constant is stored per member as
#' `genome_effect_members.center_value` and travels with its `genome_value`, so
#' evaluation applies each allele copy's own line's centering — a crossbred
#' animal's line-A alleles are centered on line A and its line-B alleles on
#' line B.
#'
#' @section Scope, and what a re-run replaces:
#' `line_name` and `parent_origin` compose into the single origin row an
#' `additive` member is allowed:
#'
#' | `line_name` | `parent_origin` | Stored scope |
#' |---|---|---|
#' | `NULL` | `NULL` | no origin rows (the common scope) |
#' | `"A"` | `NULL` | `('exact', 'A', parent NULL, copy_count = 1)` |
#' | `NULL` | `1` / `2` | `('any', NULL, parent, copy_count = 1)` |
#' | `"A"` | `1` / `2` | `('exact', 'A', parent, copy_count = 1)` |
#'
#' Re-running replaces **only the variant at the same scope**
#' (`mode = "replace_scope"`), so successive common / line-A / line-B calls each
#' keep the others: the per-copy fallback that makes crossbred breeding values
#' correct depends on all of them standing. It also means changing
#' `parent_origin` on a re-run **adds** a variant rather than replacing one —
#' the two are in a containment relation and both apply, to different copies.
#' That is legal and rarely intended, so the function warns on exactly that
#' case.
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
#'   individuals in `ind_haplotype` are used. Ignored (with a warning) when
#'   `base = "founder_haplotypes"` — use `base_line_name` there.
#' @param base_line_name Optional character, `base = "founder_haplotypes"` only.
#'   Which founder pool defines the base allele frequencies. **Defaults to
#'   `line_name`**, so line-specific effects are centered on their own line;
#'   pass `NULL` explicitly to pool every line instead. Errors if no
#'   `founder_haplotypes` rows carry that line. See *Which population centers
#'   the effects* above.
#' @param line_name Optional character. When set, effects are scoped to allele
#'   copies of this genetic line: a copy whose `line_origin` matches takes these
#'   values, and falls back per copy to the common variant where no
#'   line-specific one exists. Also becomes the default for `base_line_name`.
#'   `NULL` (default) means the common scope, matching every copy.
#' @param parent_origin Optional `1` (sire / parent_1) or `2` (dam / parent_2) —
#'   imprinting, restricting the term to copies inherited from that parent.
#'   `NULL` (default) means both parents' copies. **Per trait**: a scalar is
#'   recycled, a vector must match `trait_name` positionally, or name its
#'   entries by trait. A call mixing origins across traits while supplying `G`
#'   is rejected — under random mating the paternal and maternal copies at a
#'   locus are independent, so the requested genetic covariance between a
#'   paternal-only and a maternal-only trait is zero and cannot be realized.
#'   This replaces the removed `trait_meta.expressed_parent` flag, which could
#'   only be set trait-wide; [define_genome_effects()] expresses the per-locus
#'   case the flag never could.
#' @param scale_to_target Logical. If `TRUE`, rescale effects so the expected
#'   additive variance equals the stored `target_add_var`:
#'   `V_A = sum_j n_eligible,j * p_j q_j a_j^2`, where `n_eligible` is 2 for an
#'   unparented term and 1 for a parent-qualified one.
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
#'
#' # Crossbreeding: each line's effects centered on its own founder pool.
#' # base_line_name inherits line_name, so nothing extra is needed.
#' pop <- pop |>
#'   get_table("genome_meta") |> dplyr::filter(chr %in% 1:5) |>
#'   define_additive_effects("ADG", line_name = "Duroc")
#' pop <- pop |>
#'   get_table("genome_meta") |> dplyr::filter(chr %in% 1:5) |>
#'   define_additive_effects("ADG", line_name = "Landrace")
#'
#' # Line-specific effects, but deliberately centered on the pooled base
#' pop <- pop |>
#'   get_table("genome_meta") |> dplyr::filter(chr %in% 1:5) |>
#'   define_additive_effects("ADG", line_name = "Duroc", base_line_name = NULL)
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
                                    base_line_name  = NULL,
                                    line_name       = NULL,
                                    parent_origin   = NULL,
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
  .ge_require_effect_tables(pop)

  if (!is.null(line_name) &&
      (!is.character(line_name) || length(line_name) != 1L)) {
    stop("'line_name' must be a single character string or NULL.", call. = FALSE)
  }
  if (!is.null(line_name)) validate_sql_identifier(line_name, what = "line name")
  po <- .dae_parent_origin(parent_origin, trait_name)

  # Centre on whatever population the effect applies to: a line-specific effect
  # defaults to its own line's founder pool, a population-wide one to the pooled
  # base. missing() distinguishes "not supplied" (inherit line_name) from an
  # explicit base_line_name = NULL (force pooling even for a line-specific effect).
  base_line_name_supplied <- !missing(base_line_name)
  if (!base_line_name_supplied) {
    base_line_name <- line_name
  } else if (!is.null(base_line_name) &&
             (!is.character(base_line_name) || length(base_line_name) != 1L)) {
    stop("'base_line_name' must be a single character string or NULL.",
         call. = FALSE)
  }
  if (base != "founder_haplotypes") {
    # Only an *explicit* base_line_name is an error here. A value merely
    # inherited from line_name must not fire: `base = "current_pop"` with both
    # line_name and a line-filtered base_tbl is a legitimate, documented call.
    if (base_line_name_supplied && !is.null(base_line_name)) {
      stop("'base_line_name' applies only to base = \"founder_haplotypes\". ",
           "For base = \"current_pop\", restrict the base population with a ",
           "line-filtered 'base_tbl' instead.", call. = FALSE)
    }
    base_line_name <- NULL
  }
  if (!is.null(base_tbl) && base == "founder_haplotypes") {
    warning("'base_tbl' is ignored when base = \"founder_haplotypes\". ",
            "Use 'base_line_name' to restrict the founder pool by line.",
            call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  # ------------------------------------------------------------------ #
  #  Single-trait path                                                   #
  # ------------------------------------------------------------------ #
  if (length(trait_name) == 1L) {

    validate_sql_identifier(trait_name, what = "trait name")
    .ge_require_trait(pop$db_conn, trait_name)
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

    p_base <- compute_base_allele_freq(pop, base, base_ids, base_line_name)

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
        qtl_effects <- rescale_effects_to_target(
          qtl_tf, qtl_effects, target_add_var, p_base,
          n_eligible = .dae_n_eligible(po[[trait_name]])
        )
      }
    }

    scope <- .dae_scope(line_name, po[[trait_name]])
    built <- .dae_build(pop$db_conn, trait_name, selected_locus_names,
                        qtl_effects, as.numeric(p_base[qtl_tf]), scope)
    model <- .ge_read_model(pop$db_conn)
    drop  <- .ge_resolve_deletes(model, trait_name, GE_ADDITIVE_OWNER,
                                 "replace_scope", .ge_scope_from_origin(
                                   scope, "replace_scope"), TRUE)
    .ge_commit(pop$db_conn, drop, built)
    .dae_warn_parent_only(pop$db_conn, trait_name)

    message("Set additive effects for ", n_qtl, " QTL on trait '", trait_name,
            "' (base: ", base, "; scope: ", .dae_scope_label(line_name,
            po[[trait_name]]), ").")
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

  # A correlated draw across traits that express from different parents cannot
  # deliver the requested off-diagonal: under random mating the paternal and
  # maternal copies at a locus are independent, so the genetic covariance
  # between a paternal-only and a maternal-only trait is exactly zero however
  # strongly the sampled coefficients correlate. Unobtainable, not approximate.
  if (length(unique(vapply(trait_name, function(t) .dae_po_key(po[[t]]),
                           character(1)))) > 1L) {
    stop("'parent_origin' differs across traits in one correlated call: ",
         paste0(trait_name, " = ",
                vapply(trait_name, function(t) .dae_po_key(po[[t]]),
                       character(1)), collapse = ", "),
         ". Under random mating a locus's paternal and maternal copies are ",
         "independent, so the genetic covariance between a paternal-only and ",
         "a maternal-only trait is zero — the off-diagonal of G cannot be ",
         "realized, whatever the sampled coefficients. Use one ",
         "parent_origin for the whole call, or define the traits separately.",
         call. = FALSE)
  }

  # Validate traits exist
  trait_meta_rows <- DBI::dbGetQuery(
    pop$db_conn,
    paste0("SELECT trait_name FROM trait_meta WHERE trait_name IN (",
           sql_in_list(trait_name, what = "trait name"), ")")
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

  model <- .ge_read_model(pop$db_conn)

  for (t in trait_name) {
    if (method == "shared") {
      qtl_tf_mat[, t] <- genome_order$locus_name %in% candidate_locus_names
    } else {
      existing <- .dae_existing_loci(model, t, .dae_scope(line_name, po[[t]]))
      active   <- intersect(genome_order$locus_name[
        genome_order$locus_id %in% existing], candidate_locus_names)
      if (length(active) == 0) {
        warning("Trait '", t, "' has no existing generated additive effects at ",
                "this scope within the candidate loci; it will receive no ",
                "effects from this call.", call. = FALSE)
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

  p_base <- compute_base_allele_freq(pop, base, base_ids, base_line_name)

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
          qtl_tf_k, a_qtl, target_var[[t]], p_base,
          n_eligible = .dae_n_eligible(po[[t]])
        )
      }
    }
  }

  # One transaction for the whole call: a partially-written correlated set is
  # not a weaker version of the requested G, it is a different model.
  built <- NULL
  drop  <- integer(0)
  for (t in trait_name) {
    qtl_tf_t      <- qtl_tf_mat[, t]
    locus_names_t <- genome_order$locus_name[qtl_tf_t]
    effects_t     <- effects_mat[qtl_tf_t, t]
    non_na        <- !is.na(effects_t)
    locus_names_t <- locus_names_t[non_na]
    effects_t     <- effects_t[non_na]
    p_base_qtl_t  <- as.numeric(p_base[qtl_tf_t][non_na])
    if (length(locus_names_t) == 0L) next

    scope_t <- .dae_scope(line_name, po[[t]])
    drop <- c(drop, .ge_resolve_deletes(
      model, t, GE_ADDITIVE_OWNER, "replace_scope",
      .ge_scope_from_origin(scope_t, "replace_scope"), TRUE))
    built <- .dae_stack(built, .dae_build(pop$db_conn, t, locus_names_t,
                                          effects_t, p_base_qtl_t, scope_t))
  }
  if (is.null(built)) {
    stop("No trait received any effect from this call.", call. = FALSE)
  }
  .ge_commit(pop$db_conn, unique(drop), built)
  .dae_warn_parent_only(pop$db_conn, trait_name)

  message("Set correlated additive effects for traits: ",
          paste(trait_name, collapse = ", "), " (method: ", method,
          "; scope: ", .dae_scope_label(line_name, po[[trait_name[1]]]), ")")
  invisible(pop)
}


# ── define_additive_effects() internals ─────────────────────────────────────

#' The effect owner `define_additive_effects()` writes under
#'
#' Reserved: `add_tbv()` reads order-one `additive` variants from this owner and
#' nothing else, and the general writer refuses to touch it. Two distinct
#' defaults are the point — v4.1 shared one, so rerunning the generator in
#' replace mode would have deleted a user's own terms.
#'
#' @keywords internal
#' @noRd
GE_ADDITIVE_OWNER <- "generated_additive_tbv"

#' Resolve `parent_origin` to one value per trait
#'
#' Per trait because the deleted `trait_meta.expressed_parent` was per trait; a
#' scalar-only argument would lose exactly the expressiveness its removal was
#' meant to preserve. Accepts a scalar (recycled), a vector matching
#' `trait_name` positionally, or a vector named by trait.
#'
#' @keywords internal
#' @noRd
.dae_parent_origin <- function(parent_origin, trait_name) {
  if (is.null(parent_origin)) {
    return(stats::setNames(rep(list(NULL), length(trait_name)), trait_name))
  }
  if (!is.numeric(parent_origin) && !is.integer(parent_origin)) {
    stop("'parent_origin' must be 1 (sire / parent_1), 2 (dam / parent_2), or ",
         "NULL (both parents' copies).", call. = FALSE)
  }
  v <- stats::setNames(as.integer(parent_origin), names(parent_origin))
  if (!is.null(names(parent_origin))) {
    unknown <- setdiff(names(parent_origin), trait_name)
    if (length(unknown) > 0L) {
      stop("'parent_origin' names trait(s) not in this call: ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    out <- stats::setNames(rep(list(NULL), length(trait_name)), trait_name)
    for (nm in names(parent_origin)) out[[nm]] <- v[[nm]]
  } else if (length(v) == 1L) {
    out <- stats::setNames(rep(list(v), length(trait_name)), trait_name)
  } else if (length(v) == length(trait_name)) {
    out <- stats::setNames(as.list(v), trait_name)
  } else {
    stop("'parent_origin' must have length 1, length(trait_name) (",
         length(trait_name), "), or be named by trait.", call. = FALSE)
  }
  bad <- vapply(out, function(x) !is.null(x) && (is.na(x) || !x %in% c(1L, 2L)),
                logical(1))
  if (any(bad)) {
    stop("'parent_origin' must be 1 (sire / parent_1) or 2 (dam / parent_2); ",
         "use NULL for both.", call. = FALSE)
  }
  out
}

#' @keywords internal
#' @noRd
.dae_po_key <- function(po) if (is.null(po)) "NULL" else as.character(po)

#' Copies contributing to an additive value at one locus
#'
#' `V_A = sum_j n_eligible,j * p_j q_j a_j^2`. Two for an unparented term, one
#' for a parent-qualified one — a parent-qualified term reads a single copy, so
#' an imprinted model asked for `target_add_var = V` would otherwise land at
#' `V/2`. The enumeration is complete only because `assert_qtl_autosomal()`
#' refuses `scale_to_target = TRUE` at any locus that is not `(1,1)` for both
#' offspring sexes; without that guard a hemizygous locus would need a third
#' value and a sex ratio.
#'
#' @keywords internal
#' @noRd
.dae_n_eligible <- function(po) if (is.null(po)) 2 else 1

#' The `origin` argument for one (line_name, parent_origin) combination
#'
#' One row in every scoped case, which is all an additive member may carry.
#' `'any'` is used only for (no line, one parent) — stamping it unconditionally
#' would erase the line dimension, so two line-specific imprinted calls would
#' land on one identical scope and the second would replace the first.
#'
#' @keywords internal
#' @noRd
.dae_scope <- function(line_name, po) {
  if (is.null(line_name) && is.null(po)) return(NULL)
  if (is.null(line_name)) {
    return(list(line_match_type = "any", parent_origin = po))
  }
  c(list(line_match_type = "exact", line_name = line_name),
    if (is.null(po)) NULL else list(parent_origin = po))
}

#' @keywords internal
#' @noRd
.dae_scope_label <- function(line_name, po) {
  paste0(if (is.null(line_name)) "all lines" else paste0("line ", line_name),
         ", ",
         if (is.null(po)) "both parents' copies"
         else paste0("parent_origin ", po, " only"))
}

#' Build one trait's additive terms in the locally-indexed candidate form
#'
#' @keywords internal
#' @noRd
.dae_build <- function(conn, trait_name, locus_names, effects, centers, scope) {
  .ge_build(conn, trait_name,
            terms = data.frame(term_id       = locus_names,
                               locus_name    = locus_names,
                               contrast_name = "additive",
                               center_value  = as.numeric(centers),
                               genome_value  = as.numeric(effects),
                               stringsAsFactors = FALSE),
            origin = scope, effect_owner = GE_ADDITIVE_OWNER)
}

#' Concatenate two candidate builds, renumbering the second's local ids
#'
#' @keywords internal
#' @noRd
.dae_stack <- function(a, b) {
  if (is.null(a)) return(b)
  off <- max(a$terms$id_genome_effect)
  b$terms$id_genome_effect   <- b$terms$id_genome_effect + off
  b$members$id_genome_effect <- b$members$id_genome_effect + off
  if (nrow(b$origins) > 0L) {
    b$origins$id_genome_effect <- b$origins$id_genome_effect + off
  }
  names(b$labels) <- as.character(as.integer(names(b$labels)) + off)
  list(terms   = rbind(a$terms, b$terms),
       members = rbind(a$members, b$members),
       origins = rbind(a$origins, b$origins),
       labels  = c(a$labels, b$labels))
}

#' Loci already carrying a generated additive term for this trait at this scope
#'
#' `method = "union"` reads per-trait QTL membership from what is already
#' stored, which since v0.66.0 lives in the term/member/origin tables.
#'
#' @keywords internal
#' @noRd
.dae_existing_loci <- function(model, trait_name, scope) {
  ids <- model$terms$id_genome_effect[
    model$terms$trait_name == trait_name &
      model$terms$effect_owner == GE_ADDITIVE_OWNER]
  if (length(ids) == 0L) return(integer(0))
  sc  <- .ge_scope_from_origin(scope, "replace_scope")
  hit <- vapply(ids, function(id) {
    .ge_term_at_scope(model$members[model$members$id_genome_effect == id, ,
                                    drop = FALSE],
                      model$origins[model$origins$id_genome_effect == id, ,
                                    drop = FALSE],
                      sc)
  }, logical(1))
  unique(model$members$locus_id[model$members$id_genome_effect %in% ids[hit]])
}

#' Warn when a write leaves two variants differing only in the parent dimension
#'
#' `replace_scope` keys on origin-predicate equality and `parent_origin` is part
#' of the predicate, so `define_additive_effects(line_name = "A")` followed by
#' the same call with `parent_origin = 1` leaves **both** variants standing.
#' They are in a containment relation, so the result is a legal fallback pair —
#' paternal copies take the imprinted value, maternal copies fall back to the
#' generic one — which is correct by the rules and almost certainly not what a
#' user re-running the call intended. It is also a behaviour change from the
#' deleted `expressed_parent`, which was one trait-wide flag that could only be
#' overwritten.
#'
#' Fires only on that case: the members must agree on **every** line predicate
#' and differ only by one carrying a parent the other leaves open. The
#' common/line-A/line-B fallback differs in *line*, and a reciprocal
#' `(exact A, 1)` / `(exact A, 2)` pair has disjoint parents rather than nested
#' ones, so neither is ever reported.
#'
#' @keywords internal
#' @noRd
.dae_warn_parent_only <- function(conn, trait_names) {
  model <- .ge_read_model(conn)
  t <- model$terms[model$terms$trait_name %in% trait_names &
                     model$terms$effect_owner == GE_ADDITIVE_OWNER, ,
                   drop = FALSE]
  if (nrow(t) < 2L) return(invisible(NULL))
  keys <- vapply(t$id_genome_effect, function(id) {
    .ge_family_key(t$trait_name[t$id_genome_effect == id], GE_ADDITIVE_OWNER,
                   model$members[model$members$id_genome_effect == id, ,
                                 drop = FALSE])
  }, character(1))
  hits <- character(0)
  for (fam in split(t$id_genome_effect, keys)) {
    if (length(fam) < 2L) next
    preds <- lapply(fam, function(id) {
      .ge_predicate(model$members[model$members$id_genome_effect == id, ,
                                  drop = FALSE],
                    model$origins[model$origins$id_genome_effect == id, ,
                                  drop = FALSE])
    })
    for (i in seq_along(fam)) for (j in seq_along(fam)) {
      if (j == i) next
      if (!.ge_pred_leq(preds[[i]], preds[[j]]) ||
          .ge_pred_leq(preds[[j]], preds[[i]])) next
      if (!.dae_parent_only_pair(preds[[i]], preds[[j]])) next
      hits <- c(hits, paste0(
        .ge_scope_label(preds[[i]]), " and ", .ge_scope_label(preds[[j]]),
        " on trait '", t$trait_name[t$id_genome_effect == fam[i]], "'"))
    }
  }
  if (length(hits) == 0L) return(invisible(NULL))
  warning("This call left two additive variants that differ only in the ",
          "parent dimension: ", paste(unique(hits), collapse = "; "),
          ". Both stand — the qualified one applies to its parent's copies and ",
          "the other falls back for the rest — which is a legal fallback pair ",
          "but is rarely what re-running the same call with a new ",
          "parent_origin was meant to do. Use mode replace_owner via ",
          "define_genome_effects(), or remove the unwanted variant, if you ",
          "meant to replace it.", call. = FALSE)
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.dae_parent_only_pair <- function(P, Q) {
  if (length(P) != length(Q)) return(FALSE)
  all(vapply(seq_along(P), function(i) {
    p <- P[[i]]; q <- Q[[i]]
    if (p$kind != "additive" || q$kind != "additive") return(FALSE)
    pl <- if (identical(p$scope, "any")) "any" else p$line
    ql <- if (identical(q$scope, "any")) "any" else q$line
    identical(pl, ql)
  }, logical(1)))
}



#' Rescale QTL effects to hit a target additive variance
#'
#' `V_A = sum_j n_eligible,j * p_j q_j a_j^2`. The familiar Falconer `2pq a^2`
#' is the `n_eligible = 2` case; a parent-qualified term reads one copy, so
#' fixing the 2 would land an imprinted model at `V/2`.
#'
#' @param qtl_tf Logical mask of QTL loci, length `n_loci`.
#' @param qtl_effects Numeric effects at QTL loci, length `sum(qtl_tf)`.
#' @param target_add_var Target additive variance.
#' @param p_base Numeric vector of base allele frequencies, length `n_loci`.
#' @param n_eligible Copies contributing to the additive value at a locus: `2`
#'   for an unparented term, `1` for a parent-qualified one. The enumeration is
#'   complete only because `assert_qtl_autosomal()` refuses `scale_to_target`
#'   at any locus that is not `(1,1)` for both offspring sexes.
#' @return Rescaled `qtl_effects` vector.
#' @keywords internal
rescale_effects_to_target <- function(qtl_tf, qtl_effects, target_add_var,
                                      p_base, n_eligible = 2) {
  p_qtl <- p_base[qtl_tf]
  V_A   <- sum(n_eligible * p_qtl * (1 - p_qtl) * qtl_effects^2)
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
compute_base_allele_freq <- function(pop, base, base_ids = NULL,
                                     line_name = NULL) {
  if (base == "founder_haplotypes") {
    if (!"founder_haplotypes" %in% pop$tables) {
      stop(
        "founder_haplotypes table not found. ",
        "Did you call define_founder_haplotypes()? ",
        "Use base = 'current_pop' instead.",
        call. = FALSE
      )
    }

    if (is.null(line_name)) {
      # Pooling every line's haplotypes overstates within-line heterozygosity
      # (Wahlund), so 2pq is too large and scale_to_target under-scales the
      # effects. That is correct for a population-wide effect, which is centred
      # on the whole founder base; it is wrong for a line-specific one, so warn
      # only when no line is in play.
      n_lines <- DBI::dbGetQuery(pop$db_conn,
        "SELECT COUNT(DISTINCT line_name) AS n FROM founder_haplotypes")$n
      if (isTRUE(n_lines > 1L)) {
        warning(
          "founder_haplotypes contains ", n_lines, " lines; base allele ",
          "frequencies for this population-wide effect are pooled across all ",
          "of them, which overstates within-line heterozygosity. Set ",
          "line_name (or base_line_name) for per-line centering.",
          call. = FALSE
        )
      }
      params <- NULL
      line_filter <- ""
    } else {
      # Parameterized rather than interpolated. line_name reaches SQL in several
      # places in this file via raw paste0(); this one is bound.
      validate_sql_identifier(line_name, what = "line name")
      params <- list(line_name)
      line_filter <- "WHERE fh.line_name = ? "
    }

    freq <- DBI::dbGetQuery(pop$db_conn, paste0(
      "SELECT gm.locus_id, AVG(CAST(fh.allele AS DOUBLE)) AS f ",
      "FROM founder_haplotypes fh ",
      "JOIN genome_meta gm ON fh.locus_name = gm.locus_name ",
      line_filter,
      "GROUP BY gm.locus_id ORDER BY gm.locus_id"), params = params)
    if (nrow(freq) == 0) {
      if (is.null(line_name)) {
        stop("founder_haplotypes table is empty.", call. = FALSE)
      }
      # Must be loud: `out` below is zero-initialised, so silently returning it
      # would centre every allele at 0 and contribute nothing to the Falconer
      # V_A -- a typo'd line name would produce plausible-looking output.
      avail <- DBI::dbGetQuery(pop$db_conn,
        "SELECT DISTINCT line_name FROM founder_haplotypes ORDER BY line_name")$line_name
      stop(
        "No founder_haplotypes rows for line '", line_name, "'. Available: ",
        if (all(is.na(avail))) "only an unnamed (line_name = NULL) pool"
          else paste0("'", stats::na.omit(avail), "'", collapse = ", "),
        ".", call. = FALSE
      )
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
    where <- paste0("WHERE id_ind IN (",
                    sql_in_list(base_ids, what = "individual id"), ") ")
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


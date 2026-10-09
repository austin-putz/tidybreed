#' Define additive, dominance and epistatic QTL effects calibrated to targets
#'
#' @description
#' Selects QTL from a filtered `genome_meta` table, samples additive,
#' dominance and additive-by-additive (A x A) effects, and **calibrates** them
#' so that the stored targets `G_A`, `G_D` and `G_AA` are delivered exactly
#' under a named anchor. Writes them under the reserved effect owner
#' `"generated"`, replacing the traits' whole generated model. The
#' non-additive sibling of [define_additive_effects()]: the same pipe subject,
#' target rules, base resolution and owner.
#'
#' Use [define_additive_effects()] for line-scoped or parent-of-origin
#' additive models (crossbreeding variants); use this function when dominance
#' or epistasis may be part of the model. With an additive-only target the two
#' write identical rows for the same seed.
#'
#' @details
#' **The model.** In functional coding the genotypic value is
#' `g = sum_j a_j (x_j - 1) + sum_j d_j 1[x_j = 1] +
#' sum_(k,l) e_kl (x_k - 1)(x_l - 1)`, `x` the dosage of allele 1. Its exact
#' statistical (NOIA) re-expression has additive effects
#' `alpha_j = a_j + b_j d_j + sum_l e_jl c_l` (`c_l = 2 p_l - 1`), dominance
#' effects `d` and A x A effects `e`. The three targets are covariances of
#' those three components.
#'
#' **When the result is exact.** Every present block (`G_A`, and `G_D` /
#' `G_AA` when present) is delivered under the named anchor to a relative
#' tolerance of `1e-8` on the correlation scale. That is checked, never
#' assumed: a calibration that misses it is an error before anything is
#' written.
#'
#' **How.** Highest order first: the A x A architecture is calibrated to
#' `G_AA`; the dominance architecture, `d = h |a|` with dominance degrees
#' `h ~ N(dominance_degree_mean, dominance_degree_sd)`, is calibrated to
#' `G_D`; then, with the coupling `b d + sum e c` the first two induce now
#' fixed, the additive architecture is solved so the statistical additive
#' effects deliver `G_A`. Each stage is a `k x k` right factor of its drawn
#' architecture, so the same seed gives the same architecture.
#'
#' **The additive floor.** Dominance and epistasis already induce additive
#' variance (`b d` and `e c` are part of `alpha`). Part of it the additive
#' architecture can cancel, part it cannot; `G_A` below what it cannot cancel
#' is refused, and the error names that minimum per trait. The floor belongs
#' to the **sampled architecture**, not to `G_D` / `G_AA` alone: another seed
#' gives another floor, and with every base frequency at `0.5` under
#' `"genic"` it is zero however large `G_D` is.
#'
#' **The anchor** is the reference covariance the calibration is exact for:
#'
#' * `"genic"` (default): the random-mating (HWE + linkage-equilibrium) limit
#'   at the base allele frequencies. Use it for multi-generation studies.
#' * `"realised"`: the covariances of the individuals `base_tbl` selects,
#'   linkage disequilibrium included. The effects are still **stored** with
#'   HWE-referenced contrasts centred at the base frequencies, so the stored
#'   `additive` / `dominance` split is not the cohort's; the total is exact,
#'   and [extract_genetic_variance()]`(anchor = "realised")` on the same
#'   individuals gives the targets back. `base_tbl` must select individuals
#'   with complete genotypes, and the in-memory design arrays have a size
#'   limit.
#'
#' Under linkage disequilibrium the `additive` component of a realised cohort
#' is a contrast component, not a joint regression: calibration is exact for
#' the anchor and nothing more.
#'
#' After calibration the delivered blocks are compared with what another
#' population sees: the genic limit under `"realised"`; the base individuals'
#' own covariances under `"genic"` when `base_tbl` selects individuals; the
#' founder pool's expectation of the additive block only (dominance and A x A
#' would need the pool's multi-locus LD). A departure outside `warn_bounds`
#' warns (the pool comparison is a message). Nothing is stored.
#'
#' @section Targets:
#' Each block (`additive`, `dominance`, `additive_by_additive`) comes from
#' exactly one source:
#'
#' * a passed matrix (`G_A`, `G_D`, `G_AA`; a number for one trait), written
#'   to `trait_var_comp` in the same transaction as the effects. It is refused
#'   when that block is already stored for any of the traits, even as an
#'   identical matrix, and whether or not `trait_var_comp_tbl` shows it;
#' * otherwise the rows of `trait_var_comp_tbl` (a filtered
#'   `get_table(pop, "trait_var_comp")`), or, when it is `NULL`, the stored
#'   population-wide rows.
#'
#' A block that comes from neither is **absent** from the model; a zero matrix
#' is an explicit exact zero, and writes no terms. An `additive` block is
#' required. A model with no non-zero dominance or A x A block takes the
#' additive-only route: the same draw, calibration and rows as
#' [define_additive_effects()]. To add a block later, store or pass it and
#' re-run: the additive architecture draw is unchanged by the new block (same
#' seed).
#'
#' Non-zero dominance or A x A targets need a positive-definite `G_A`: this
#' release solves the additive stage for a full-rank factor. A correlation of
#' +/-1 or a zero additive variance is accepted for additive-only models.
#' Singular targets are reported with a message, since a typed `1` meant as
#' `0.99` looks the same.
#'
#' @section Pairs:
#' A x A acts on pairs of QTL, used only with an `additive_by_additive` block:
#'
#' * `pairs`, a data frame `locus_name_1`, `locus_name_2`: your design. A
#'   locus may appear in several pairs (a hub gene); unknown loci, loci
#'   outside the filter, self-pairs and repeated pairs (in either order) are
#'   errors.
#' * `pairs = NULL`: a random matching, each QTL in at most one pair, as
#'   AlphaSimR does. `n_pairs = NULL` pairs every QTL once
#'   (`floor(m / 2)` pairs), and a message says so.
#'
#' Within a pair the loci are ordered by name (C locale), and pairs are
#' written sorted by locus id, so the stored order does not depend on the
#' draw.
#'
#' @section Crossbreeding and scope:
#' Effects are common to every line (no `line_name` or `parent_origin`; line
#' or parent-specific non-additive effects are a later release). Common
#' effects still give heterosis, through allele-frequency differences between
#' lines. The default base pools the founder table, so the targets then hold
#' at the **pooled** frequencies (with a Wahlund warning); pass `base_tbl` to
#' calibrate on one reference line.
#'
#' @section What a re-run replaces:
#' Every `"generated"` term of the traits, at every scope: line-scoped
#' additive variants from [define_additive_effects()] included (the message
#' counts them). Terms of other owners are untouched. Afterwards
#' [define_additive_effects()] refuses these traits while their generated
#' model has dominance or A x A terms; re-run this function with an
#' additive-only target, or use [remove_generated_effects()].
#'
#' @section Generated means calibrated:
#' There is no option to write fixed or unscaled effects. Exact coefficients
#' go through [define_genome_effect_terms()] with [ad_terms()] /
#' [aa_terms()] under a user owner, for example:
#'
#' ```
#' pop |>
#'   define_genome_effect_terms(
#'     trait_name = "ADG",
#'     terms = rbind(
#'       ad_terms(locus_name = c("Locus_10", "Locus_44"),
#'                a = c(0.30, -0.12), d = c(0.10, 0.05), p = c(0.35, 0.60),
#'                coding = "functional"),
#'       aa_terms(locus_name_1 = "Locus_10", locus_name_2 = "Locus_44",
#'                e = 0.08, p_1 = 0.35, p_2 = 0.60, coding = "functional")))
#' ```
#'
#' @section Inbreeding depression:
#' `inbreeding_depression` (the drop in the mean per unit of inbreeding `F`,
#' positive = depression, `sum 2pq d`) sets each named trait's mean
#' dominance degree. It is exact for one trait. With several traits each
#' requested trait's mean is solved as if it stood alone, and the joint
#' calibration then mixes the columns, so the delivered depression is
#' approximate, with no closeness guarantee; a message gives requested and
#' delivered values. It is not stored.
#'
#' @param tbl A `tidybreed_table` from [get_table()]`("genome_meta")`,
#'   optionally filtered. Its loci are the QTL of every trait.
#' @param trait_name Character vector of existing traits.
#' @param G_A,G_D,G_AA Optional `k x k` targets (a number for one trait),
#'   named in `trait_name` order or unnamed. Written with the effects; never
#'   over a stored block. `NULL` reads the block from `trait_var_comp_tbl`.
#' @param trait_var_comp_tbl Optional filtered
#'   `get_table(pop, "trait_var_comp")`: the stored rows to use for blocks not
#'   passed. `NULL` reads the stored population-wide rows.
#' @param pairs Optional data frame `locus_name_1`, `locus_name_2`.
#' @param n_pairs Optional number of random pairs, `1..floor(m / 2)`.
#' @param anchor `"genic"` (default) or `"realised"`.
#' @param dominance_degree_mean,dominance_degree_sd Mean and standard
#'   deviation of the dominance degrees that shape the dominance
#'   architecture. Defaults `0.19` and `0.097`.
#' @param inbreeding_depression Optional numeric vector named by trait.
#' @param base_tbl Optional `tidybreed_table` selecting the base population,
#'   as in [define_additive_effects()]. Required under `"realised"`.
#' @param warn_bounds `c(lower, upper)` for the comparison with another
#'   population, or `NULL` to turn it off.
#'
#' @return The `tidybreed_pop`, invisibly.
#'
#' @seealso [define_additive_effects()], [extract_genetic_variance()],
#'   [define_effect_cov_matrix()], [remove_generated_effects()],
#'   [aa_terms()].
#'
#' @examples
#' \dontrun{
#' pop <- pop |> define_trait("ADG")
#' set.seed(1)
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(chr_name %in% c("1", "2")) |>
#'   define_genome_effects("ADG", G_A = 0.3, G_D = 0.1, G_AA = 0.05)
#'
#' # Additive-only first, dominance added later: the additive draw is the same
#' set.seed(2)
#' pop <- pop |> get_table("genome_meta") |> define_genome_effects("BF", G_A = 1)
#' pop <- pop |> define_effect_cov_matrix("dominance", 0.2, trait_name = "BF")
#' set.seed(2)
#' pop <- pop |> get_table("genome_meta") |> define_genome_effects("BF")
#' }
#' @export
define_genome_effects <- function(tbl,
                                  trait_name,
                                  G_A                   = NULL,
                                  G_D                   = NULL,
                                  G_AA                  = NULL,
                                  trait_var_comp_tbl    = NULL,
                                  pairs                 = NULL,
                                  n_pairs               = NULL,
                                  anchor                = c("genic", "realised"),
                                  dominance_degree_mean = 0.19,
                                  dominance_degree_sd   = 0.097,
                                  inbreeding_depression = NULL,
                                  base_tbl              = NULL,
                                  warn_bounds           = c(0.8, 1.25)) {

  # ── 1. Arguments. No RNG and no write until the draw (step 8). ─────────────
  if (!inherits(tbl, "tidybreed_table")) {
    stop("'tbl' must be a tidybreed_table from get_table('genome_meta') |> ",
         "filter(...). Use: pop |> get_table('genome_meta') |> ",
         "define_genome_effects()", call. = FALSE)
  }
  if (tbl$table_name != "genome_meta") {
    stop("'tbl' must be piped from get_table('genome_meta'), not '",
         tbl$table_name, "'.", call. = FALSE)
  }
  if (!is.character(trait_name) || length(trait_name) == 0L ||
      anyNA(trait_name)) {
    stop("`trait_name` must be a character vector of trait names.",
         call. = FALSE)
  }
  pop <- tbl$pop
  validate_tidybreed_pop(pop)
  anchor <- match.arg(anchor)
  conn   <- pop$db_conn
  .ge_require_effect_tables(pop)
  lapply(trait_name, validate_sql_identifier, what = "trait name")
  if (anyDuplicated(trait_name)) {
    stop("`trait_name` must not contain duplicates.", call. = FALSE)
  }
  k <- length(trait_name)
  .qtl_validate_scalar(dominance_degree_mean, "dominance_degree_mean")
  .qtl_validate_scalar(dominance_degree_sd, "dominance_degree_sd", min = 0)
  .qtl_validate_warn_bounds(warn_bounds)
  inbreeding_depression <- .dge_check_inbreeding(inbreeding_depression,
                                                 trait_name, dominance_degree_sd)
  pairs <- .dge_check_pairs_arg(pairs)
  if (!is.null(n_pairs)) {
    .qtl_validate_scalar(n_pairs, "n_pairs", integral = TRUE, min = 1)
    n_pairs <- as.integer(round(n_pairs))
  }
  if (!is.null(pairs) && !is.null(n_pairs)) {
    stop("Pass `pairs` (a chosen design) or `n_pairs` (random pairs), not ",
         "both.", call. = FALSE)
  }
  if (k > 1L && !requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for multi-trait effect sampling. ",
         "Install with install.packages('MASS').", call. = FALSE)
  }
  for (t in trait_name) .ge_require_trait(conn, t)
  if (anchor == "realised") .dae_check_realised(base_tbl, NULL, list())

  # ── 2. Targets, once per block (§6C, D8). ─────────────────────────────────
  tg <- .dge_resolve_targets(conn, trait_name, list(
    additive = G_A, dominance = G_D, additive_by_additive = G_AA),
    trait_var_comp_tbl)
  G   <- tg$G
  std <- tg$std
  has_d  <- !is.null(G$dominance)
  has_aa <- !is.null(G$additive_by_additive)
  live_d  <- has_d  && std$dominance$rank > 0L
  live_aa <- has_aa && std$additive_by_additive$rank > 0L
  if (!has_aa && (!is.null(pairs) || !is.null(n_pairs))) {
    stop("`", if (!is.null(pairs)) "pairs" else "n_pairs", "` is used only ",
         "with an 'additive_by_additive' target (`G_AA`, or a stored block), ",
         "and this call has none.", call. = FALSE)
  }
  if (!is.null(inbreeding_depression)) {
    if (!has_d) {
      stop("`inbreeding_depression` needs a 'dominance' target (`G_D`, or a ",
           "stored block): inbreeding depression is a property of the ",
           "dominance effects.", call. = FALSE)
    }
    zero_d <- names(inbreeding_depression)[
      diag(G$dominance)[names(inbreeding_depression)] == 0]
    if (length(zero_d)) {
      stop("`inbreeding_depression` names ", paste(zero_d, collapse = ", "),
           ", whose dominance variance is 0: with no dominance effects there ",
           "is no inbreeding depression to set.", call. = FALSE)
    }
  }
  route <- if (live_d || live_aa) "nonadditive" else "additive"
  if (route == "nonadditive" &&
      (!all(std$additive$pos) || std$additive$rank < k)) {
    stop("This release calibrates dominance or epistasis only with a ",
         "positive-definite `G_A` (the sampled additive architecture must ",
         "have full rank). The 'additive' target for ",
         paste(trait_name, collapse = ", "), " is singular (rank ",
         std$additive$rank, " of ", k, "). A correlation of +/-1 or a zero ",
         "additive variance is accepted for additive-only models (no ",
         "non-zero `G_D` / `G_AA`).", call. = FALSE)
  }
  if (live_d && dominance_degree_mean == 0 && dominance_degree_sd == 0) {
    stop("`dominance_degree_mean` and `dominance_degree_sd` are both 0, so ",
         "every dominance effect drawn is 0 and no calibration can scale it ",
         "to the non-zero 'dominance' target. This restricts the ",
         "architecture parameters, not the target: set a non-zero mean or ",
         "sd.", call. = FALSE)
  }

  # ── 3. What the call replaces (§5). ───────────────────────────────────────
  replaced <- .dge_generated_counts(conn, trait_name)

  # ── 4. Loci and base, lightweight. ────────────────────────────────────────
  loci_df <- dplyr::collect(tbl)
  if (!"locus_name" %in% names(loci_df)) {
    stop("The filtered table must contain 'locus_name'. Pipe ",
         "get_table('genome_meta') into define_genome_effects().",
         call. = FALSE)
  }
  if (nrow(loci_df) == 0L) {
    stop("No loci selected — your filter returned zero rows.", call. = FALSE)
  }
  genome_order <- DBI::dbGetQuery(conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id")
  candidate <- genome_order$locus_name %in% loci_df$locus_name
  mask <- matrix(candidate, nrow = nrow(genome_order), ncol = k,
                 dimnames = list(NULL, trait_name))
  qtl_name <- genome_order$locus_name[candidate]
  qtl_id   <- genome_order$locus_id[candidate]
  m <- length(qtl_name)
  base   <- .dae_resolve_base(pop, base_tbl, NULL)
  p_base <- base$p_base
  .dae_require_base_at(p_base, candidate, genome_order$locus_name)
  assert_qtl_autosomal(conn, qtl_name)
  n_ind <- if (anchor == "realised") {
    DBI::dbGetQuery(conn, paste0(
      "SELECT COUNT(*) AS n FROM (SELECT DISTINCT id_ind FROM (",
      as.character(dbplyr::sql_render(base$tbl$tbl)), ") b)"))$n
  }

  # ── 5. Pairs, and every shortage knowable before the draw. ────────────────
  pair_idx <- NULL
  r <- 0L
  if (has_aa) {
    if (!is.null(pairs)) {
      pair_idx <- .dge_supplied_pairs(pairs, genome_order, candidate)
      pair_idx <- .dge_canonical_pairs(pair_idx, qtl_name, qtl_id)
      r <- nrow(pair_idx)
    } else {
      if (m < 2L) {
        stop("Random A x A pairs need at least 2 QTL; the filter selects ",
             m, ". Select more loci, or pass `pairs`.", call. = FALSE)
      }
      max_pairs <- m %/% 2L
      if (is.null(n_pairs)) n_pairs <- max_pairs
      if (n_pairs > max_pairs) {
        stop("`n_pairs` = ", n_pairs, " exceeds the ", max_pairs, " pairs a ",
             "random matching of ", m, " QTL can hold (each QTL in at most ",
             "one pair, floor(m / 2)). For designs where a locus is in ",
             "several pairs, pass `pairs`.", call. = FALSE)
      }
      r <- n_pairs
    }
    if (live_aa && std$additive_by_additive$rank > r) {
      stop("The 'additive_by_additive' target has rank ",
           std$additive_by_additive$rank, " but there ",
           if (r == 1L) "is 1 pair" else paste0("are ", r, " pairs"),
           ": ", r, " pair column", if (r != 1L) "s", " cannot carry a rank ",
           "above ", r, " under any anchor. Add pairs.", call. = FALSE)
    }
  }
  if (anchor == "realised") {
    for (e in names(G)) {
      if (!is.null(G[[e]]) && std[[e]]$rank > n_ind - 1L) {
        stop("Under anchor = \"realised\" a centred cohort of ", n_ind,
             " individuals has at most ", n_ind - 1L, " independent ",
             "directions, but the '", e, "' target has rank ", std[[e]]$rank,
             ". Select more base individuals, or use anchor = \"genic\".",
             call. = FALSE)
      }
    }
  }

  # ── 6. Dosage guard, then the designs (realised only). ────────────────────
  design <- NULL
  if (anchor == "realised") {
    if (route == "nonadditive") {
      .dge_dosage_guard(n_ind, m, if (live_d) m else 0L, if (live_aa) r else 0L)
    }
    design <- .dae_collect_dosages(pop, base$tbl, qtl_id)
  }

  # ── 7. Anchor feasibility per block. ──────────────────────────────────────
  w_all <- 2 * p_base * (1 - p_base)
  anc_A <- .dae_anchor_at(anchor, candidate, w_all, design,
                          genome_order$locus_id)
  .qtl_anchor_rank_check(anc_A$rank(1e-10), std$additive$rank)
  na <- NULL
  if (route == "nonadditive") {
    na <- .na_anchors(anchor, p = p_base[candidate], X = design$X)
    if (live_d) {
      .qtl_anchor_rank_check(na$D$rank(1e-10), std$dominance$rank, "dominance")
    }
    if (live_aa && !is.null(pair_idx)) {
      na <- .na_aa_anchor(na, pair_idx)
      .qtl_anchor_rank_check(na$AA$rank(1e-10), std$additive_by_additive$rank,
                             "additive_by_additive")
    }
  }

  # ── 8. Draw (D2 order). The only RNG use. ─────────────────────────────────
  dr <- .dge_draw(mask, G$additive, has_d, has_aa, pair_idx, n_pairs,
                  qtl_name, qtl_id)

  # ── 9. Calibrate. ─────────────────────────────────────────────────────────
  if (route == "additive") {
    B <- .dae_calibrate_shared(dr$B, candidate, std$additive, anc_A)
    cal <- list(route = "additive", B_alpha = B[candidate, , drop = FALSE],
                delivered = list(A = .dge_named(anc_A$cov(B[candidate, ,
                                                            drop = FALSE]),
                                                trait_name)))
  } else {
    if (live_aa && is.null(pair_idx)) {
      na <- .na_aa_anchor(na, dr$pairs)
      rk <- na$AA$rank(1e-10)
      if (rk < std$additive_by_additive$rank) {
        stop("The drawn random pair design is rank-deficient under the \"",
             anchor, "\" anchor: its ", nrow(dr$pairs), " pairs have rank ",
             rk, " but the 'additive_by_additive' target has rank ",
             std$additive_by_additive$rank, ". Select more segregating loci, ",
             "pass more `n_pairs`, or pass `pairs`. Nothing was written.",
             call. = FALSE)
      }
    } else if (has_aa && is.null(na$AA)) {
      na <- .na_aa_anchor(na, dr$pairs)
    }
    cal <- .na_calibrate(na, G$additive, G$dominance, G$additive_by_additive,
                         B_a = dr$B_a, z = dr$z, B_aa = dr$B_aa,
                         pairs = dr$pairs,
                         dominance_degree_mean = dominance_degree_mean,
                         dominance_degree_sd = dominance_degree_sd,
                         inbreeding_depression = inbreeding_depression)
  }

  # ── 10. Diagnostics: computed before the commit, reported after. ──────────
  diag_res <- if (!is.null(warn_bounds)) {
    if (route == "additive") {
      .dae_diagnostics(pop, base$tbl, anchor, genome_order, candidate,
                       cal$B_alpha, std$additive, p_base)
    } else {
      .dge_diagnostics(pop, base$tbl, anchor, cal, na, std, qtl_id,
                       genome_order, candidate, p_base, dr$pairs)
    }
  }

  # ── 11. Build, and commit terms and passed targets in one transaction. ────
  model <- .ge_read_model(conn)
  drop  <- integer(0)
  for (t in trait_name) {
    drop <- c(drop, .ge_resolve_deletes(model, t, GE_GENERATED_OWNER,
                                        "replace_owner", NULL, TRUE))
  }
  built <- if (route == "additive") {
    .dae_build_traits(conn, trait_name, B, mask, genome_order$locus_name,
                      p_base, function(t) NULL)
  } else {
    .dge_build_nonadditive(conn, trait_name, cal, qtl_name, p_base[candidate],
                           dr$pairs)
  }
  write <- names(G)[tg$write]
  write_targets <- if (length(write)) {
    G_write <- G[write]
    function(conn) {
      for (e in names(G_write)) .tvc_write_block(conn, e, G_write[[e]], NULL)
    }
  }
  .ge_commit(conn, unique(drop), built, before_commit = write_targets)

  # ── 12. Messages. Reported, never stored. ─────────────────────────────────
  .dge_messages(trait_name, replaced, has_aa, is.null(pairs), dr$pairs, m,
                anchor, base$label, cal, G, std, dominance_degree_mean,
                dominance_degree_sd, inbreeding_depression, na)
  if (!is.null(diag_res)) {
    if (route == "additive") {
      .dae_report_diagnostics(diag_res, anchor, warn_bounds)
    } else {
      .dge_report_diagnostics(diag_res, anchor, warn_bounds)
    }
  }
  for (e in names(G)) if (!is.null(G[[e]])) .qtl_rank_note(G[[e]], e)
  invisible(pop)
}


# ── define_genome_effects() internals ───────────────────────────────────────

#' @keywords internal
#' @noRd
.dge_named <- function(x, traits) {
  x <- (x + t(x)) / 2
  dimnames(x) <- list(traits, traits)
  x
}

#' Validate `inbreeding_depression`: named by call traits, aligned by name
#' @keywords internal
#' @noRd
.dge_check_inbreeding <- function(x, trait_name, sd) {
  if (is.null(x)) return(NULL)
  nm <- names(x)
  if (!is.numeric(x) || length(x) == 0L || anyNA(x) || any(!is.finite(x))) {
    stop("`inbreeding_depression` must be a finite numeric vector named by ",
         "trait.", call. = FALSE)
  }
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
    stop("`inbreeding_depression` must name each of its traits once, e.g. ",
         "c(ADG = 1.2); it is aligned by name, never by position.",
         call. = FALSE)
  }
  unknown <- setdiff(nm, trait_name)
  if (length(unknown)) {
    stop("`inbreeding_depression` names trait(s) not in this call: ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  if (sd <= 0) {
    stop("`inbreeding_depression` needs `dominance_degree_sd` > 0: with a ",
         "constant dominance degree the ratio of inbreeding depression to ",
         "dominance standard deviation is fixed by the additive architecture ",
         "and cannot be targeted.", call. = FALSE)
  }
  # Ordered as the call's traits, so the result does not depend on the order
  # the names were typed in.
  x[intersect(trait_name, nm)]
}

#' Validate the `pairs` argument's shape (its loci are checked in step 5)
#' @keywords internal
#' @noRd
.dge_check_pairs_arg <- function(pairs) {
  if (is.null(pairs)) return(NULL)
  if (!is.data.frame(pairs) ||
      !setequal(names(pairs), c("locus_name_1", "locus_name_2")) ||
      ncol(pairs) != 2L) {
    stop("`pairs` must be a data frame with exactly the columns ",
         "locus_name_1 and locus_name_2.", call. = FALSE)
  }
  if (nrow(pairs) == 0L) {
    stop("`pairs` has no rows. Leave it NULL for random pairs.", call. = FALSE)
  }
  out <- data.frame(locus_name_1 = as.character(pairs$locus_name_1),
                    locus_name_2 = as.character(pairs$locus_name_2),
                    stringsAsFactors = FALSE)
  if (anyNA(out) || any(!nzchar(as.matrix(out)))) {
    stop("`pairs` must name a locus in every cell.", call. = FALSE)
  }
  out
}

#' Resolve the three target blocks, each from exactly one source (D8)
#'
#' @param passed Named list of passed matrices (`NULL` = not passed).
#' @return list(G = named list of named k x k blocks or `NULL` (absent),
#'   std = their `.qtl_target_std()`, write = logical, passed per block).
#' @keywords internal
#' @noRd
.dge_resolve_targets <- function(conn, trait_name, passed, trait_var_comp_tbl) {
  explicit <- !is.null(trait_var_comp_tbl)
  rows_all <- if (explicit) {
    .tvc_collect_explicit(trait_var_comp_tbl, conn)
  } else {
    .dae_default_target_rows(conn, trait_name, NULL)
  }
  rows_all <- rows_all[rows_all$trait_name_1 %in% trait_name |
                         rows_all$trait_name_2 %in% trait_name, , drop = FALSE]
  arg <- c(additive = "G_A", dominance = "G_D", additive_by_additive = "G_AA")
  G <- stats::setNames(vector("list", 3L), names(arg))
  std <- G
  write <- stats::setNames(logical(3L), names(arg))
  for (e in names(arg)) {
    rows_e <- rows_all[rows_all$effect_name == e, , drop = FALSE]
    if (!is.null(passed[[e]])) {
      M <- .check_cov_dimnames(passed[[e]], trait_name, arg[[e]])
      if (anyNA(M) || any(!is.finite(M))) {
        stop("`", arg[[e]], "` must contain only finite values.", call. = FALSE)
      }
      .qtl_target_std(unname(M), name = arg[[e]])
      if (explicit && nrow(rows_e) > 0L) {
        stop("The '", e, "' block comes from two sources: `", arg[[e]],
             "` and rows of `trait_var_comp_tbl`. Each block comes from one; ",
             "drop `", arg[[e]], "`, or filter the block out of ",
             "`trait_var_comp_tbl`.", call. = FALSE)
      }
      block <- .tvc_block_traits(conn, e, NULL, trait_name)
      if (length(block)) {
        stop("A '", e, "' block is already stored for ",
             paste(block, collapse = ", "), " (population-wide). Stored ",
             "targets are never overwritten, even by an identical matrix. ",
             "Drop `", arg[[e]], "` to calibrate to the stored block, or ",
             "remove it and re-run this call, which writes the new target and ",
             "re-draws the terms in one transaction. To remove it first:\n",
             .tvc_removal_call(e, NULL, block), call. = FALSE)
      }
      G[[e]] <- (M + t(M)) / 2
      write[[e]] <- TRUE
    } else {
      G[e] <- list(.tvc_block_from_rows(conn, e, trait_name, rows_e, explicit,
                                        NULL, "all lines, both parents' copies"))
    }
    if (!is.null(G[[e]])) std[[e]] <- .qtl_target_std(unname(G[[e]]), arg[[e]])
  }
  if (is.null(G$additive)) {
    stop("No 'additive' target for ", paste(trait_name, collapse = ", "),
         if (explicit) " in `trait_var_comp_tbl`" else " is stored",
         ". Every model has an additive block: pass `G_A =`, or store one ",
         "with define_effect_cov_matrix(pop, \"additive\", ...).",
         call. = FALSE)
  }
  list(G = G, std = std, write = write)
}

#' Count the generated terms a call replaces, and the line-scoped ones
#' @return data frame trait_name, n, n_line.
#' @keywords internal
#' @noRd
.dge_generated_counts <- function(conn, trait_name) {
  gm <- .gev_read_model(conn, trait_name, effect_owner = GE_GENERATED_OWNER)
  tl <- if (nrow(gm$terms)) .gev_term_line(gm) else character(0)
  data.frame(
    trait_name = trait_name,
    n      = vapply(trait_name, function(t) sum(gm$terms$trait_name == t),
                    integer(1)),
    n_line = vapply(trait_name, function(t)
      sum(gm$terms$trait_name == t & !is.na(tl)), integer(1)),
    stringsAsFactors = FALSE, row.names = NULL)
}

#' Map supplied pairs to QTL positions, refusing bad keys
#' @return r x 2 integer matrix of positions among the QTL (`locus_id` order).
#' @keywords internal
#' @noRd
.dge_supplied_pairs <- function(pairs, genome_order, candidate) {
  l1 <- pairs$locus_name_1; l2 <- pairs$locus_name_2
  show <- function(x) paste0("'", utils::head(unique(x), 10L), "'",
                             collapse = ", ")
  unknown <- setdiff(c(l1, l2), genome_order$locus_name)
  if (length(unknown)) {
    stop("`pairs` names loci not in genome_meta: ", show(unknown), ".",
         call. = FALSE)
  }
  qtl <- genome_order$locus_name[candidate]
  outside <- setdiff(c(l1, l2), qtl)
  if (length(outside)) {
    stop("`pairs` names loci outside the filtered QTL set: ", show(outside),
         ". Every pair locus must be a QTL of this call.", call. = FALSE)
  }
  self <- l1 == l2
  if (any(self)) {
    stop("`pairs` pairs a locus with itself: ", show(l1[self]), ".",
         call. = FALSE)
  }
  key <- paste(pmin(l1, l2), pmax(l1, l2), sep = "\r")
  if (anyDuplicated(key)) {
    d <- duplicated(key)
    stop("`pairs` repeats a pair (in either order): ",
         paste0("(", l1[d], ", ", l2[d], ")", collapse = ", "),
         ". Give each pair once.", call. = FALSE)
  }
  cbind(match(l1, qtl), match(l2, qtl))
}

#' Canonical pair order: by name within a pair, by locus id across pairs
#'
#' Within a pair the loci take [aa_terms()]'s C-locale order; pairs are then
#' sorted by `(locus_id_1, locus_id_2)`, so the written order depends on
#' neither the draw nor the platform.
#' @keywords internal
#' @noRd
.dge_canonical_pairs <- function(P, qtl_name, qtl_id) {
  rk <- integer(length(qtl_name))
  rk[order(qtl_name, method = "radix")] <- seq_along(qtl_name)
  swap <- rk[P[, 1]] > rk[P[, 2]]
  P[swap, ] <- P[swap, 2:1, drop = FALSE]
  P <- P[order(qtl_id[P[, 1]], qtl_id[P[, 2]]), , drop = FALSE]
  storage.mode(P) <- "integer"
  P
}

#' Refuse realised design arrays above the cell limit
#'
#' Counts the columns the non-additive route keeps: `m` additive, `m`
#' dominance with a non-zero dominance target, `r` pair columns with a
#' non-zero A x A target. The additive-only route keeps Part A's `n x m`
#' guard (in `.dae_collect_dosages()`).
#' @keywords internal
#' @noRd
.dge_dosage_guard <- function(n, m, m_D, r) {
  cells <- as.numeric(n) * (m + m_D + r)
  if (cells > QTL_REALISED_MAX_CELLS) {
    stop("define_genome_effects(anchor = \"realised\") would keep design ",
         "arrays of ", n, " individuals x (", m, " additive + ", m_D,
         " dominance + ", r, " pair columns) = ",
         format(cells, big.mark = ",", scientific = FALSE), " cells, above ",
         "the limit of ", format(QTL_REALISED_MAX_CELLS, big.mark = ",",
                                 scientific = FALSE), ". The limit bounds ",
         "the retained design arrays, not peak memory (sweeps and pair ",
         "products allocate temporaries). Use anchor = \"genic\", or a ",
         "smaller base_tbl, QTL set or pair set.", call. = FALSE)
  }
  invisible(NULL)
}

#' Draw every architecture, in the D2 order
#'
#' 1. The additive architecture, `.draw_additive_architecture(mask, "normal",
#'    G_A)` (shared with [define_additive_effects()], gate C4).
#' 2. With a dominance block (zero or not): standard-normal degree deviations
#'    `z` (m x k); the degrees are `mean + sd z`.
#' 3. With an A x A block (zero or not) and no supplied pairs: one `sample()`
#'    permutation of the QTL in `locus_id` order, paired off consecutively,
#'    the first `n_pairs` kept, then put in canonical order.
#' 4. With an A x A block: the architecture `B_aa` (r x k).
#'
#' Each block's draws follow those of every block that is always present, so
#' adding A x A changes neither the additive nor the dominance draw (gate
#' C17).
#'
#' @return list(B = the full n_loci x k draw, B_a = its QTL rows, z, pairs,
#'   B_aa).
#' @keywords internal
#' @noRd
.dge_draw <- function(mask, G_A, has_d, has_aa, pairs, n_pairs, qtl_name,
                      qtl_id) {
  B   <- .draw_additive_architecture(mask, "normal", G_A)
  B_a <- B[mask[, 1L], , drop = FALSE]
  m <- nrow(B_a); k <- ncol(B_a)
  z <- if (has_d) matrix(stats::rnorm(m * k), m, k)
  if (has_aa && is.null(pairs)) {
    perm  <- sample(m)
    pairs <- matrix(perm[seq_len(2L * n_pairs)], ncol = 2L, byrow = TRUE)
    pairs <- .dge_canonical_pairs(pairs, qtl_name, qtl_id)
  }
  B_aa <- if (has_aa) matrix(stats::rnorm(nrow(pairs) * k), nrow(pairs), k)
  list(B = B, B_a = B_a, z = z, pairs = if (has_aa) pairs, B_aa = B_aa)
}

#' Store a non-additive model: Cockerham terms at the base frequencies
#'
#' Per trait, the functional `(a, d, e)` become statistical (NOIA) terms at
#' the base `p` (`.noia_to_stored()`, `.noia_terms()`): `ad_terms()` +
#' `aa_terms()` rows. Exact zeros are dropped by the builders; a zero
#' dominance or A x A column writes no terms of that kind.
#' @keywords internal
#' @noRd
.dge_build_nonadditive <- function(conn, trait_name, cal, qtl_name, p, pairs) {
  built <- NULL
  pk <- stats::setNames(p, qtl_name)
  for (j in seq_along(trait_name)) {
    a <- stats::setNames(cal$B_a[, j], qtl_name)
    d <- if (!is.null(cal$B_d)) stats::setNames(cal$B_d[, j], qtl_name) else
      numeric(0)
    pr <- if (!is.null(cal$B_aa)) {
      data.frame(locus_1 = qtl_name[pairs[, 1]], locus_2 = qtl_name[pairs[, 2]],
                 e = cal$B_aa[, j], stringsAsFactors = FALSE)
    } else {
      data.frame(locus_1 = character(0), locus_2 = character(0),
                 e = numeric(0), stringsAsFactors = FALSE)
    }
    stat  <- .noia_to_stored(a, d, pr, pk)
    built <- .dae_stack(built, .ge_build(conn, trait_name[j], .noia_terms(stat),
                                         NULL, GE_GENERATED_OWNER))
  }
  built
}

#' Compare each delivered block with what another population sees (D6)
#'
#' * `"realised"`: each block's genic limit at the base individuals'
#'   frequencies (closed forms).
#' * `"genic"`, base selects individuals: each block realised on them
#'   (A x A accumulated in pair chunks).
#' * `"genic"`, founder pool: the additive block's pool expectation only.
#' * `"genic"`, `ind_haplotype` copies: none.
#'
#' Zero blocks are skipped. Nothing is stored.
#' @return NULL, or list(spec = named list of spectra, label, pool, n_h,
#'   note).
#' @keywords internal
#' @noRd
.dge_diagnostics <- function(pop, base_tbl, anchor, cal, na, std, qtl_id,
                             genome_order, candidate, p_base, pairs) {
  live <- function(e) !is.null(std[[e]]) && std[[e]]$rank > 0L
  blocks <- function(A, D, AA) {
    out <- list()
    if (live("additive")) out$additive <- .qtl_relative_spectrum_std(std$additive, A)
    if (live("dominance")) out$dominance <- .qtl_relative_spectrum_std(std$dominance, D)
    if (live("additive_by_additive")) {
      out$additive_by_additive <- .qtl_relative_spectrum_std(
        std$additive_by_additive, AA)
    }
    out
  }
  k <- ncol(cal$B_a)
  genic_at <- function(p) {
    w <- 2 * p * (1 - p)
    C <- .na_coupling(nrow(cal$B_a), k, (1 - p) - p, 2 * p - 1, cal$B_d,
                      cal$B_aa, pairs)
    alpha <- cal$B_a + C
    list(A = crossprod(alpha, w * alpha),
         D = if (!is.null(cal$B_d)) crossprod(cal$B_d, w^2 * cal$B_d),
         AA = if (!is.null(cal$B_aa))
           crossprod(cal$B_aa, (w[pairs[, 1]] * w[pairs[, 2]]) * cal$B_aa))
  }
  if (anchor == "realised") {
    g <- genic_at(na$p)
    return(list(spec = blocks(g$A, g$D, g$AA),
                label = "genic limit (expectation)", pool = FALSE))
  }
  if (base_tbl$table_name == "founder_haplotypes") {
    res <- .dae_diagnostics(pop, base_tbl, "genic", genome_order, candidate,
                            cal$B_alpha, std$additive, p_base)
    if (!is.null(res) && is.null(res$note)) {
      res$spec <- list(additive = res$spec)
      res$label <- paste0("founder pool (pool expectation under random ",
                          "pairing; additive block only, dominance and ",
                          "additive-by-additive not compared)")
    }
    return(res)
  }
  if (base_tbl$table_name == "ind_haplotype") return(NULL)
  X <- tryCatch(.dae_collect_dosages(pop, base_tbl, qtl_id)$X,
                error = function(e) NULL)
  if (is.null(X)) {
    return(list(note = paste0(
      "Skipped the observed comparison: the base individuals' genotypes ",
      "could not be collected (partial or above the in-memory limit).")))
  }
  ob <- .na_anchors("realised", X = X)
  n  <- nrow(X)
  C  <- .na_coupling(nrow(cal$B_a), k, ob$b, ob$cc, cal$B_d, cal$B_aa, pairs)
  A  <- crossprod(ob$Z_A %*% (cal$B_a + C)) / (n - 1)
  D  <- if (!is.null(cal$B_d)) ob$D$cov(cal$B_d)
  AA <- if (!is.null(cal$B_aa)) {
    v <- .egv_aa_values(ob$Z_A, pairs, cal$B_aa, .egv_pair_chunk(n))
    crossprod(v) / (n - 1)
  }
  list(spec = blocks(A, D, AA), label = "base individuals (observed)",
       pool = FALSE)
}

#' Report a `.dge_diagnostics()` result: one warning naming each block outside
#' `warn_bounds`; the founder pool's comparison is a message
#' @keywords internal
#' @noRd
.dge_report_diagnostics <- function(res, anchor, warn_bounds) {
  if (!is.null(res$note)) {
    message(res$note)
    return(invisible(NULL))
  }
  if (length(res$spec) == 0L) return(invisible(NULL))
  if (isTRUE(res$pool)) {
    one <- res
    one$spec <- res$spec$additive
    return(.dae_report_diagnostics(one, anchor, warn_bounds))
  }
  out <- vapply(res$spec, function(sp) {
    length(sp) > 0L && (min(sp) < warn_bounds[1] || max(sp) > warn_bounds[2])
  }, logical(1))
  if (any(out)) {
    warning("The calibration is exact under the ", anchor, " anchor, but the ",
            res$label, " sees covariances departing from the targets: ",
            paste0(names(res$spec)[out], " relative spectrum ",
                   vapply(res$spec[out], .qtl_spectrum_text, character(1)),
                   collapse = "; "),
            ", outside warn_bounds [", warn_bounds[1], ", ", warn_bounds[2],
            "]. warn_bounds = NULL turns this off.", call. = FALSE)
  }
  invisible(res$spec)
}

#' The closing messages (phase-5 plan 5b.7)
#' @keywords internal
#' @noRd
.dge_messages <- function(trait_name, replaced, has_aa, random_pairs, pairs, m,
                          anchor, base_label, cal, G, std, dd_mean, dd_sd,
                          inbreeding_depression, na) {
  for (i in which(replaced$n > 0L)) {
    message("Replaced ", replaced$n[i], " generated term",
            if (replaced$n[i] != 1L) "s", " of trait '",
            replaced$trait_name[i], "'",
            if (replaced$n_line[i] > 0L) paste0(
              " (", replaced$n_line[i], " of them line-scoped: crossbred ",
              "variants from define_additive_effects())"), ".")
  }
  if (has_aa && random_pairs) {
    r <- nrow(pairs)
    message("Drew ", r, " random A x A pair", if (r != 1L) "s",
            if (r == m %/% 2L) paste0(": every QTL paired once (floor(", m,
                                      " / 2)), as AlphaSimR does") else
              paste0(" (at most floor(", m, " / 2) = ", m %/% 2L, ")"),
            ". Pass n_pairs for fewer pairs, or pairs for a chosen design (a ",
            "locus may then appear in several pairs).")
  }
  lbl <- c(additive = "additive", dominance = "dominance",
           additive_by_additive = "additive-by-additive")
  key <- c(additive = "A", dominance = "D", additive_by_additive = "AA")
  parts <- character(0)
  for (e in names(lbl)) {
    if (is.null(G[[e]])) next
    parts <- c(parts, if (std[[e]]$rank == 0L) {
      paste0(lbl[[e]], ": zero target, no terms")
    } else {
      paste0(lbl[[e]], " ", .dae_format_cov(cal$delivered[[key[[e]]]]))
    })
  }
  message("Set genome effects for ", m, " QTL",
          if (has_aa) paste0(" and ", nrow(pairs), " pair",
                             if (nrow(pairs) != 1L) "s"),
          " on trait", if (length(trait_name) > 1L) "s", " ",
          paste(trait_name, collapse = ", "), " (base: ", base_label,
          "): exact under the ", anchor, " anchor; delivered ",
          paste(parts, collapse = "; "), ".")
  if (cal$route == "additive") return(invisible(NULL))

  fl <- diag(cal$floor)
  message("Additive floor for this sampled architecture (the least additive ",
          "variance its dominance and A x A effects allow; another draw gives ",
          "another floor): ", paste0(trait_name, " = ", signif(fl, 4),
                                     collapse = ", "),
          ". G_A - floor has smallest eigenvalue ",
          signif(cal$G_tilde_min, 4), " on the correlation scale.")

  if (!is.null(cal$B_d)) {
    ib <- cal$inbreeding
    head <- paste0("Inbreeding depression (drop in the mean per unit of F, ",
                   "positive = depression; sum 2pq d under the ", anchor,
                   " anchor's frequencies): ")
    if (is.null(inbreeding_depression)) {
      message(head, "implied ", paste0(trait_name, " = ",
                                       signif(ib$delivered, 4),
                                       collapse = ", "), ".")
    } else {
      nm <- names(inbreeding_depression)
      message(head, paste0(nm, " requested ", signif(inbreeding_depression, 6),
                           ", delivered ", signif(ib$delivered[nm], 6),
                           collapse = "; "),
              if (length(trait_name) > 1L) paste0(
                ". Approximate for several traits: each requested trait's ",
                "degree mean is solved as if it stood alone and the joint ",
                "calibration mixes the columns, so there is no closeness ",
                "guarantee"), ".")
    }
  }

  solved <- cal$inbreeding$mean[!is.na(cal$inbreeding$mean)]
  per <- vapply(seq_along(trait_name), function(j) {
    paste0(trait_name[j], ": ", sum(cal$B_a[, j] != 0), " QTL with a != 0",
           if (!is.null(cal$B_d)) paste0(", ", sum(cal$B_d[, j] != 0),
                                         " with d != 0"),
           if (!is.null(cal$B_aa)) paste0(", ", sum(cal$B_aa[, j] != 0),
                                          " pairs with e != 0"))
  }, character(1))
  message("Functional effects (sampling: dominance_degree_mean = ", dd_mean,
          ", dominance_degree_sd = ", dd_sd,
          if (length(solved)) paste0("; solved degree mean ",
                                     paste0(names(solved), " = ",
                                            signif(solved, 6),
                                            collapse = ", ")),
          "): ", paste(per, collapse = "; "), ".")
  invisible(NULL)
}

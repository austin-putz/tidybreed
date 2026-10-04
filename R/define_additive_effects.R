#' Define additive QTL effects for one or more traits
#'
#' @description
#' Selects QTL from a filtered `genome_meta` table, samples an effect
#' architecture, **calibrates** it to the stored additive target, and writes
#' one order-one `additive` term per locus and trait through the same engine as
#' [define_genome_effect_terms()], under the reserved effect owner
#' `"generated"`. [add_tbv()] reads order-one `additive` variants from that
#' owner and nothing else, so effects written here and effects a user writes
#' with [define_genome_effect_terms()] can never be confused for one another.
#'
#' @details
#' **When the result is exact.** The requested covariance `G` is delivered
#' exactly, `B' M B = G` to machine precision, **only** for sampled effects
#' with `method = "shared"` and `scale_to_target = TRUE`, when the rank is
#' feasible: `rank(G) <= rank(M)` (the anchor has enough independent
#' segregating directions at the selected loci) and `rank(G) <=
#' rank(B0' M B0)` (the drawn architecture does too). Each infeasibility is its
#' own error. Rank and positive semidefiniteness are judged on the target's
#' **correlation** scale, so a trait recorded in small units is never
#' truncated away. "Exact" is checked, not assumed: the delivered `B' M B` is
#' compared with `G` as stored, entry by entry, to a relative tolerance of
#' `1e-8` on the correlation scale; a calibration that misses it (a
#' numerically ill-conditioned architecture) is an error before anything is
#' written. The closing message says "exact" or "approximate" and gives the
#' delivered covariance under the anchor.
#'
#' **How.** Effects are drawn as today (one draw per QTL for one trait; joint
#' `MVN(0, G)` rows for several), and that draw is only the *architecture*
#' `B0`. It is then right-multiplied by a `k x k` matrix `A` so that
#' `B = B0 A` satisfies `B' M B = G` exactly (the congruence of
#' Proposition 2 in the source method). For one trait this is exactly the
#' scalar rescale `b0 * sqrt(G / sum(w b0^2))`. For two or more it also fixes
#' the genetic **correlations**, which a per-trait rescale cannot: at 200 QTL
#' and a target correlation of 0.4 a scalar rescale delivers anything from
#' about 0.18 to 0.60. The same seed gives the same `B0`.
#'
#' **The anchor `M`** is the reference-population genotype covariance the
#' calibration is exact for:
#'
#' * `anchor = "genic"` (default): `M = diag(n_eligible * p * q)` at the base
#'   allele frequencies, the random-mating (HWE + linkage-equilibrium) limit.
#'   Use it for multi-generation studies: under random mating the realised
#'   covariance converges to it. The manuscript's "reference" anchor is this
#'   plus a non-default `base_tbl`, whose frequencies the weights use.
#' * `anchor = "realised"`: `M = Cov(X)` of the individuals `base_tbl` selects,
#'   linkage disequilibrium included. The single-generation / clonal option.
#'   `base_tbl` must select **individuals** (a table with `id_ind`, not
#'   `founder_haplotypes` or `ind_haplotype`) with complete genotypes at the
#'   QTL, and the call must use the common scope (no `line_name`, no
#'   `parent_origin`). The genotype matrix is collected into memory, so there
#'   is a size limit; above it the call errors and suggests `"genic"`.
#'
#' Under `method = "union"` each trait keeps its own QTL set and is scaled by
#' its own variance only, so the variances are exact and the covariances are
#' **approximate** -- including a zero target covariance, which overlapping QTL
#' sets do not deliver. A warning gives the delivered covariance and
#' correlation and names `method = "shared"` as the exact option. A trait with
#' a positive target variance and no QTL in the call is an error.
#'
#' After calibration the delivered covariance is compared with what another
#' population sees (§7.4 of the plan): the **pool expectation** `2 Cov(H)`
#' when the base is the founder pool, the **observed** `Cov(X)` when it selects
#' individuals, and the genic limit under `anchor = "realised"`. The founder
#' pool's comparison is always a message: a small pool's departure is its
#' sampling LD, not a mistake in the call. To get the target exactly in the
#' founders, add them first and calibrate with `anchor = "realised"` on them.
#' The observed and genic-limit comparisons warn when the relative spectrum
#' leaves `warn_bounds`. Nothing is stored.
#'
#' @section Targets:
#' The target is the population-wide (or line) `additive` block of
#' `trait_var_comp`, the single source of generation targets:
#'
#' * pass `G` (a `k x k` matrix, or a number for one trait) to write it **and**
#'   calibrate to it, in the same transaction as the effects. If a block is
#'   already stored for any of the traits at that `line_name` the call errors,
#'   even for an identical matrix, and gives the [remove_rows()] call;
#' * or leave `G = NULL` to use the stored rows, optionally chosen with
#'   `trait_var_comp_tbl = get_table(pop, "trait_var_comp") |> filter(...)`.
#'   The rows for the call's traits must form one complete symmetric block. A
#'   stored block that pairs one of the traits with a trait outside the call
#'   is an error (calibrating one trait alone would break the stored
#'   covariance), as is a stored `dominance` or `additive_by_additive` block
#'   for the traits: this generator calibrates the additive block only and
#'   never silently ignores a stored target. Filter them away with
#'   `trait_var_comp_tbl` to say so explicitly.
#'
#' With `line_name = "C"` the default reads line C's block when one exists and
#' otherwise the population-wide one.
#'
#' Manual `effects` are written unchanged and take no target; so do sampled
#' effects with `scale_to_target = FALSE` (which for several traits still use
#' the stored `G` as the draw's covariance). `G` with either is an error:
#' nothing would be calibrated to it.
#'
#' @section Which population centers the effects:
#' Base allele frequencies center the true breeding value (the Falconer
#' `allele - p` term) and set the genic weights `n_eligible * p q`. They come
#' from `base_tbl`, a filtered `tidybreed_table` whose identity says *what kind
#' of thing* is selected and whose [dplyr::filter()] says *which* (see
#' [extract_allele_freq()] for the three accepted shapes: the founder pool,
#' allele copies in `ind_haplotype`, or individuals from any table with
#' `id_ind`).
#'
#' When `base_tbl = NULL` the base is **the population the effect applies to**,
#' resolved with the same `line -> NULL` precedence as `resolve_genome_map()`:
#' a line-specific effect (`line_name = "A"`) centers on line A's own founder
#' pool, or on the shared (`line_name = NULL`) pool when no named pool exists;
#' a population-wide effect (`line_name = NULL`) centers on the whole founder
#' table. Only that last case warns when the founder table holds more than one
#' pool: pooling divergent lines overstates within-line heterozygosity — the
#' Wahlund effect — so the calibration **under**-scales the effects and the
#' realised within-line additive variance falls short of the target. An
#' explicit `base_tbl` is an intentional selection and never warns — pass
#' `base_tbl = get_table(pop, "founder_haplotypes")` to pool on purpose, which
#' is how the common fallback variant of a crossbreeding model is defined.
#'
#' A selected QTL locus with no allele copies in the base is an error, never
#' silently centered at `p = 0`.
#'
#' The centering constant is stored per member as
#' `genome_effect_members.center_value` and travels with its `genome_value`, so
#' evaluation applies each allele copy's own line's centering — a crossbred
#' animal's line-A alleles are centered on line A and its line-B alleles on
#' line B.
#'
#' @section Lines and line means:
#' Line-specific effects (`line_name = "A"`) are centred on line A's own base,
#' so they add **no** difference between line means. Differences between lines
#' come from allele-frequency differences at QTL whose effects are shared:
#' use common effects centred on one reference line, e.g.
#' `base_tbl = get_table(pop, "founder_haplotypes") |>
#' filter(line_name == "Terminal")`. The calibration then hits the target
#' within that line, and the other lines get whatever variance their
#' frequencies give.
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
#'   in `trait_meta`. When length >= 2, the architecture is drawn jointly from
#'   `MVN(0, G)` and `method` becomes active.
#' @param effects Optional numeric vector of length `n_qtl` (manual mode, single
#'   trait only), in ascending `locus_id` order, written unchanged. Error if
#'   `length(trait_name) > 1` or with `G`.
#' @param distribution Character. `"normal"` (default) or `"gamma"`, the
#'   single-trait architecture. Ignored for multi-trait (always MVN).
#' @param G Optional additive-genetic (co)variance target: a `k x k` matrix
#'   (named in `trait_name` order, or unnamed), or a single number for one
#'   trait. Written to `trait_var_comp` (with the call's `line_name`) in the
#'   same transaction as the effects, and never over a stored block. `NULL`
#'   reads the stored target (see *Targets*).
#' @param trait_var_comp_tbl Optional filtered
#'   `get_table(pop, "trait_var_comp")`: the stored rows to calibrate to. Use it
#'   to pick a block explicitly, e.g. to leave a stored non-additive block out
#'   or to calibrate one trait of a stored block alone. Not with `G`.
#' @param anchor Character. `"genic"` (default) or `"realised"`: the
#'   reference covariance the calibration is exact for. See *Details*.
#' @param method Character. `"shared"` (default) or `"union"`. Multi-trait only.
#'   `"shared"` — all listed traits use the filtered loci as their shared QTL
#'   set; exact. `"union"` — per-trait QTL sets are read from existing generated
#'   terms at this scope, restricted to the filtered loci; per-trait scaling,
#'   approximate for a non-zero covariance.
#' @param base_tbl Optional `tidybreed_table` from [get_table()] (optionally
#'   filtered) selecting the allele copies that define base allele
#'   frequencies: `founder_haplotypes`, `ind_haplotype`, or any table with an
#'   `id_ind` column. Must come from the same `pop` as `tbl`. `NULL` (default)
#'   resolves to the founder pool of the line the effect applies to — see
#'   *Which population centers the effects*. Under `anchor = "realised"` it is
#'   required and names the individuals whose `Cov(X)` is the anchor.
#' @param line_name Optional character. When set, effects are scoped to allele
#'   copies of this genetic line: a copy whose `line_origin` matches takes these
#'   values, and falls back per copy to the common variant where no
#'   line-specific one exists. Also selects the default `base_tbl` and the
#'   line's own target block. `NULL` (default) means the common scope.
#' @param parent_origin Optional `1` (sire / parent_1) or `2` (dam / parent_2) —
#'   imprinting, restricting the term to copies inherited from that parent.
#'   `NULL` (default) means both parents' copies. **Per trait**: a scalar is
#'   recycled, a vector must match `trait_name` positionally, or name its
#'   entries by trait; each value must be exactly `1` or `2`. One call must
#'   use one origin for every trait: it calibrates against one reference
#'   covariance, defined for one set of inherited copies (a mixed-scope anchor
#'   is not supported yet). Define differently-scoped traits in separate
#'   calls. For imprinting that varies locus by locus, write the terms with
#'   [define_genome_effect_terms()].
#' @param scale_to_target Logical. `TRUE` (default) calibrates sampled effects
#'   to the target. `FALSE` writes the draw unscaled and takes no target.
#' @param warn_bounds Numeric length 2, `c(lower, upper)` with
#'   `0 < lower <= upper`, or `NULL` to turn the comparison off. Outside
#'   the bounds, an observed or genic-limit comparison warns; a founder-pool
#'   comparison adds the realised-anchor hint to its message.
#'   Default `c(0.8, 1.25)` (±25%, multiplicative).
#' @param seed Optional integer, applied with [set.seed()] immediately before
#'   the draw, after every input check and after the anchor's feasibility
#'   check (`rank(G) <= rank(M)`), so those refusals never touch the RNG. A
#'   failure that depends on the draw itself -- the drawn architecture's rank,
#'   or a calibration that fails verification -- comes after the seed and has
#'   consumed RNG draws.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_trait()], [define_effect_cov_matrix()],
#'   [define_genome_effect_terms()], [extract_allele_freq()].
#'
#' @examples
#' \dontrun{
#' # Single trait — all loci on chr 1-5 become QTL; write the target and
#' # calibrate to it in one call
#' pop <- pop |> define_trait("ADG")
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   dplyr::filter(chr_name %in% as.character(1:5)) |>
#'   define_additive_effects("ADG", G = 0.25)
#'
#' # Correlated traits — shared QTL set, exact G (variances and correlation)
#' pop <- pop |> define_trait("BW")
#' G <- matrix(c(0.25, 0.10, 0.10, 0.30), 2, 2,
#'             dimnames = list(c("ADG", "BW"), c("ADG", "BW")))
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   define_additive_effects(c("ADG", "BW"), G = G)
#'
#' # A target stored beforehand is read back: no G
#' pop <- pop |> define_trait("FCR") |>
#'   define_effect_cov_matrix("additive", 0.1, trait_name = "FCR")
#' pop <- pop |> get_table("genome_meta") |> define_additive_effects("FCR")
#'
#' # Exact in the generation-0 individuals themselves, LD included
#' pop <- pop |>
#'   get_table("genome_meta") |>
#'   define_additive_effects("FCR", anchor = "realised",
#'     base_tbl = get_table(pop, "ind_meta") |> dplyr::filter(gen == 0L))
#'
#' # Crossbreeding: three variants. The common fallback names its base to say
#' # "yes, pool"; each line's variant centers on its own founder pool by default.
#' gm <- pop |> get_table("genome_meta") |> dplyr::filter(chr_name == "1")
#' pop <- gm |> define_additive_effects("ADG",
#'                base_tbl = get_table(pop, "founder_haplotypes"))
#' pop <- gm |> define_additive_effects("ADG", line_name = "Duroc")
#' pop <- gm |> define_additive_effects("ADG", line_name = "Landrace")
#' }
#' @export
define_additive_effects <- function(tbl,
                                    trait_name,
                                    effects            = NULL,
                                    distribution       = c("normal", "gamma"),
                                    G                  = NULL,
                                    trait_var_comp_tbl = NULL,
                                    anchor             = c("genic", "realised"),
                                    method             = c("shared", "union"),
                                    base_tbl           = NULL,
                                    line_name          = NULL,
                                    parent_origin      = NULL,
                                    scale_to_target    = TRUE,
                                    warn_bounds        = c(0.8, 1.25),
                                    seed               = NULL) {

  # ── 1. Arguments. Nothing below draws from the RNG until step 4. ──────────
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
  anchor       <- match.arg(anchor)
  method       <- match.arg(method)
  conn         <- pop$db_conn
  .ge_require_effect_tables(pop)

  lapply(trait_name, validate_sql_identifier, what = "trait name")
  if (anyDuplicated(trait_name)) {
    stop("`trait_name` must not contain duplicates.", call. = FALSE)
  }
  k <- length(trait_name)
  if (!is.null(line_name) &&
      (!is.character(line_name) || length(line_name) != 1L || is.na(line_name))) {
    stop("'line_name' must be a single character string or NULL.", call. = FALSE)
  }
  if (!is.null(line_name)) validate_sql_identifier(line_name, what = "line name")
  po <- .dae_parent_origin(parent_origin, trait_name)
  .qtl_validate_warn_bounds(warn_bounds)
  if (!is.logical(scale_to_target) || length(scale_to_target) != 1L ||
      is.na(scale_to_target)) {
    stop("`scale_to_target` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(seed)) .qtl_validate_scalar(seed, "seed", integral = TRUE)

  if (!is.null(G) && (!is.null(effects) || !scale_to_target)) {
    stop("`G` is a calibration target for sampled effects only: ",
         if (!is.null(effects)) "manual `effects` are written unchanged"
         else "`scale_to_target = FALSE` writes the draw unscaled",
         ", so nothing would be calibrated to it and trait_var_comp would ",
         "record a target the model does not deliver. Such effects take no ",
         "target: drop `G`.", call. = FALSE)
  }
  if (!is.null(G) && !is.null(trait_var_comp_tbl)) {
    stop("Pass `G` (a new target) or `trait_var_comp_tbl` (stored rows), ",
         "not both.", call. = FALSE)
  }
  if (k > 1L && !is.null(effects)) {
    stop("'effects' cannot be used when 'trait_name' has length > 1.", call. = FALSE)
  }
  if (k > 1L && distribution != "normal") {
    warning("'distribution' is ignored for multi-trait; effects are always drawn from MVN.",
            call. = FALSE)
  }
  if (k > 1L && !requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for multi-trait effect sampling. ",
         "Install with install.packages('MASS').", call. = FALSE)
  }
  for (t in trait_name) .ge_require_trait(conn, t)

  # One call calibrates against one anchor, and the anchor's weights
  # n_eligible * p q depend on which copies a term reads (2 for both parents,
  # 1 for one). Traits reading different copies need a block anchor over the
  # paternal and maternal copies separately, which does not exist yet. The
  # covariance is not always zero: paternal-only vs maternal-only is zero under
  # random mating, but both-parents vs paternal-only is pq per locus.
  po_keys <- vapply(trait_name, function(t) .dae_po_key(po[[t]]), character(1))
  if (length(unique(po_keys)) > 1L) {
    stop("'parent_origin' differs across traits in one call: ",
         paste0(trait_name, " = ", po_keys, collapse = ", "),
         ". One call calibrates all its traits against one reference ",
         "covariance, which is defined for one set of inherited copies; a ",
         "mixed-scope anchor is not supported yet. Use one parent_origin for ",
         "the whole call, or define the traits in separate calls (their ",
         "genetic covariance is then whatever the effects give: zero in ",
         "expectation for a paternal-only vs a maternal-only trait under ",
         "random mating, but not in general when one trait reads both copies).",
         call. = FALSE)
  }
  n_elig <- .dae_n_eligible(po[[trait_name[1L]]])

  if (anchor == "realised") .dae_check_realised(base_tbl, line_name, po,
                                                effects, scale_to_target)

  # ── 2. The target (§6C). ──────────────────────────────────────────────────
  need <- if (!is.null(effects)) "none"
          else if (!scale_to_target) (if (k > 1L) "sigma" else "none")
          else "target"
  tgt <- .dae_resolve_target(pop, trait_name, line_name, G,
                             trait_var_comp_tbl, need)
  G <- tgt$G

  # ── 3. Loci, base, and the anchor's data. ─────────────────────────────────
  loci_df <- dplyr::collect(tbl)
  if (!"locus_name" %in% names(loci_df)) {
    stop("The filtered table must contain 'locus_name'. ",
         "Pipe get_table('genome_meta') into define_additive_effects().",
         call. = FALSE)
  }
  if (nrow(loci_df) == 0L) {
    stop("No loci selected — your filter returned zero rows.", call. = FALSE)
  }
  # Everything downstream (effects, p_base, the written members) is in
  # locus_id order, so take the names from genome_meta rather than from the
  # collected filter, whose order is whatever the projection left.
  genome_order <- DBI::dbGetQuery(conn,
    "SELECT locus_id, locus_name FROM genome_meta ORDER BY locus_id")
  n_loci    <- nrow(genome_order)
  candidate <- genome_order$locus_name %in% loci_df$locus_name

  model <- .ge_read_model(conn)
  mask  <- matrix(candidate, nrow = n_loci, ncol = k,
                  dimnames = list(NULL, trait_name))
  if (k > 1L && method == "union") {
    for (t in trait_name) {
      existing <- .dae_existing_loci(model, t, .dae_scope(line_name, po[[t]]))
      active   <- candidate & genome_order$locus_id %in% existing
      if (!any(active)) {
        # A positive target with no QTL to carry it would be stored (or kept)
        # as a target the model does not deliver.
        if (need == "target" && G[t, t] > 0) {
          stop("Trait '", t, "' has no existing generated additive effects at ",
               "this scope within the candidate loci, so method = \"union\" ",
               "cannot deliver its target variance ", format(G[t, t]), ". ",
               "Define its QTL first, or use method = \"shared\".",
               call. = FALSE)
        }
        if (need == "target") {
          message("Trait '", t, "' has no QTL in this call; its target ",
                  "variance is 0, so it receives no effects.")
        } else {
          warning("Trait '", t, "' has no existing generated additive effects ",
                  "at this scope within the candidate loci; it will receive no ",
                  "effects from this call.", call. = FALSE)
        }
      }
      mask[, t] <- active
    }
    if (!any(mask)) {
      stop("No QTL found across any trait in the candidate loci.", call. = FALSE)
    }
  }
  any_qtl <- apply(mask, 1L, any)

  base       <- .dae_resolve_base(pop, base_tbl, line_name)
  p_base     <- base$p_base
  for (t in trait_name) {
    .dae_require_base_at(p_base, mask[, t], genome_order$locus_name,
                         trait = if (k > 1L) t)
  }
  if (!is.null(effects)) {
    if (!is.numeric(effects)) stop("`effects` must be numeric.", call. = FALSE)
    if (length(effects) != sum(candidate)) {
      stop("`effects` length (", length(effects), ") must equal number of selected loci (",
           sum(candidate), ").", call. = FALSE)
    }
  }
  if (need == "target") {
    assert_qtl_autosomal(conn, genome_order$locus_name[any_qtl])
  }
  design <- if (anchor == "realised")
    .dae_collect_dosages(pop, base$tbl, genome_order$locus_id[any_qtl])

  w_all <- n_elig * p_base * (1 - p_base)
  anchor_at <- function(rows) {
    if (anchor == "genic") return(.qtl_anchor_diag(w_all[rows]))
    cols <- match(genome_order$locus_id[rows], design$locus_id)
    Xc   <- sweep(design$X[, cols, drop = FALSE], 2L,
                  colMeans(design$X[, cols, drop = FALSE]), "-")
    .qtl_anchor_design(Xc, nrow(Xc) - 1L)
  }

  # Anchor feasibility depends on the base alone, so it is refused before the
  # seed: only an architecture-rank failure (step 5) depends on the draw.
  std <- if (need == "target") .qtl_target_std(G)
  if (need == "target") {
    if (k == 1L || method == "shared") {
      .qtl_anchor_rank_check(anchor_at(mask[, 1L])$rank(1e-10), std$rank)
    } else {
      for (t in trait_name) {
        if (G[t, t] > 0) .qtl_anchor_rank_check(anchor_at(mask[, t])$rank(1e-10), 1L)
      }
    }
  }

  # ── 4. Seed, then draw. No RNG use before this line. ──────────────────────
  if (!is.null(seed)) set.seed(seed)
  B <- if (!is.null(effects)) {
    out <- matrix(NA_real_, n_loci, 1L, dimnames = list(NULL, trait_name))
    out[candidate, 1L] <- as.numeric(effects)
    out
  } else {
    .draw_additive_architecture(mask, distribution, G)
  }

  # ── 5. Calibrate, verified against the target as stored. ─────────────────
  exact <- NA
  if (need == "target") {
    if (k == 1L || method == "shared") {
      rows <- mask[, 1L]
      B[rows, ] <- .qtl_calibrate(B[rows, , drop = FALSE], std,
                                  anchor_at(rows))$B
    } else {
      for (t in trait_name) {
        rows <- mask[, t]
        if (!any(rows)) next
        B[rows, t] <- .qtl_calibrate(B[rows, t, drop = FALSE],
                                     .qtl_target_std(G[t, t, drop = FALSE]),
                                     anchor_at(rows))$B
      }
    }
  }

  # Delivered covariance under the anchor, on every locus any trait uses.
  # "exact" is a check of this against the target, never assumed: under
  # "union" even a zero target covariance is missed where the QTL sets overlap.
  B_any     <- B[any_qtl, , drop = FALSE]
  B_any[is.na(B_any)] <- 0
  delivered <- anchor_at(any_qtl)$cov(B_any)
  dimnames(delivered) <- list(trait_name, trait_name)
  if (need == "target") {
    exact <- .qtl_target_error(delivered, std) <= QTL_CALIBRATION_TOL
  }
  if (need == "target" && isFALSE(exact)) {
    warning("method = \"union\" is approximate: each trait keeps its own QTL ",
            "set and is scaled by its own variance only, so the covariances ",
            "are not calibrated. Delivered covariance ", .dae_format_cov(delivered),
            " (correlation ", .dae_format_cor(delivered), ") against target ",
            .dae_format_cov(G), ". method = \"shared\" is exact.", call. = FALSE)
  }

  # ── 6. Diagnostics: what another population sees (§7.4). Computed before
  # the commit, so a failure here leaves the database untouched; reported
  # after it. Nothing stored.
  diag_res <- if (need == "target" && !is.null(warn_bounds) &&
                  is.null(line_name) && is.null(po[[trait_name[1L]]])) {
    .dae_diagnostics(pop, base$tbl, anchor, genome_order, any_qtl, B_any, std,
                     p_base)
  }

  # ── 7. Build, and commit terms with any new target in one transaction. ────
  built <- NULL
  drop  <- integer(0)
  for (t in trait_name) {
    rows <- mask[, t] & !is.na(B[, t])
    if (!any(rows)) next
    scope_t <- .dae_scope(line_name, po[[t]])
    drop <- c(drop, .ge_resolve_deletes(
      model, t, GE_GENERATED_OWNER, "replace_scope",
      .ge_scope_from_origin(scope_t, "replace_scope"), TRUE))
    built <- .dae_stack(built, .dae_build(conn, t, genome_order$locus_name[rows],
                                          B[rows, t], p_base[rows], scope_t))
  }
  write_target <- if (tgt$write) {
    G_write <- G
    function(conn) .tvc_write_block(conn, "additive", G_write, line_name)
  }
  .ge_commit(conn, unique(drop), built, before_commit = write_target)
  .dae_warn_parent_only(conn, trait_name)
  if (!is.null(diag_res)) .dae_report_diagnostics(diag_res, anchor, warn_bounds)

  # ── 8. Messages. ──────────────────────────────────────────────────────────
  calib <- if (need != "target") {
    "effects written unchanged (not calibrated)"
  } else {
    paste0(if (isTRUE(exact)) "exact" else "approximate", " under the ",
           anchor, " anchor; delivered ", .dae_format_cov(delivered))
  }
  scope_lbl <- .dae_scope_label(line_name, po[[trait_name[1L]]])
  if (k == 1L) {
    message("Set additive effects for ", sum(mask[, 1L]), " QTL on trait '",
            trait_name, "' (base: ", base$label, "; scope: ", scope_lbl,
            "): ", calib, ".")
  } else {
    message("Set correlated additive effects for traits: ",
            paste(trait_name, collapse = ", "), " (method: ", method,
            "; base: ", base$label, "; scope: ", scope_lbl, "): ", calib, ".")
  }
  if (!is.null(line_name)) {
    message("Line-specific effects are centred on line ", line_name,
            "'s own base, so they add no difference between line means. For ",
            "differences between lines, use common effects centred on one ",
            "reference line (see ?define_additive_effects).")
  }
  invisible(pop)
}


# ── define_additive_effects() internals ─────────────────────────────────────

#' The effect owner `define_additive_effects()` writes under
#'
#' Reserved: `add_tbv()` reads order-one `additive` variants from this owner and
#' nothing else, and the general writer refuses to touch it. Keeping it distinct
#' from the writer's `"custom"` default is what stops a rerun of the generator in
#' replace mode from deleting a user's own terms.
#'
#' @keywords internal
#' @noRd
GE_GENERATED_OWNER <- "generated"

#' Resolve `parent_origin` to one value per trait
#'
#' Per trait, because imprinting is a property of a trait and one correlated
#' call may define several. Accepts a scalar (recycled), a vector matching
#' `trait_name` positionally, or a vector named by trait.
#'
#' @keywords internal
#' @noRd
.dae_parent_origin <- function(parent_origin, trait_name) {
  if (is.null(parent_origin)) {
    return(stats::setNames(rep(list(NULL), length(trait_name)), trait_name))
  }
  # Checked on the values as given, before any coercion: as.integer(1.9) is 1,
  # which would silently turn a typo into a paternal-only effect.
  if (!is.numeric(parent_origin) || length(parent_origin) == 0L ||
      anyNA(parent_origin) || !all(parent_origin %in% c(1, 2))) {
    stop("'parent_origin' must be 1 (sire / parent_1) or 2 (dam / parent_2) ",
         "for every trait, or NULL (both parents' copies); got ",
         paste(format(parent_origin), collapse = ", "), ".", call. = FALSE)
  }
  nms <- names(parent_origin)
  if (!is.null(nms) && (anyNA(nms) || any(nms == "") || anyDuplicated(nms))) {
    stop("A named 'parent_origin' must name each trait once, with no empty ",
         "or duplicated names.", call. = FALSE)
  }
  v <- stats::setNames(as.integer(parent_origin), nms)
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
  out
}

#' @keywords internal
#' @noRd
.dae_po_key <- function(po) if (is.null(po)) "NULL" else as.character(po)

#' Copies contributing to an additive value at one locus
#'
#' `V_A = sum_j n_eligible,j * p_j q_j a_j^2`. Two for an unparented term, one
#' for a parent-qualified one — a parent-qualified term reads a single copy, so
#' an imprinted model asked for an additive target `V` would otherwise land at
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
            origin = scope, effect_owner = GE_GENERATED_OWNER)
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
#' stored in the term/member/origin tables.
#'
#' @keywords internal
#' @noRd
.dae_existing_loci <- function(model, trait_name, scope) {
  ids <- model$terms$id_genome_effect[
    model$terms$trait_name == trait_name &
      model$terms$effect_owner == GE_GENERATED_OWNER]
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
#' user re-running the call intended.
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
                     model$terms$effect_owner == GE_GENERATED_OWNER, ,
                   drop = FALSE]
  if (nrow(t) < 2L) return(invisible(NULL))
  keys <- .ge_family_keys(t, model$members[
    model$members$id_genome_effect %in% t$id_genome_effect, , drop = FALSE])
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
          "define_genome_effect_terms(), or remove the unwanted variant, if you ",
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


#' Largest genotype matrix (individuals x QTL cells) collected into memory
#'
#' `anchor = "realised"` forms `Cov(X)` in R from the collected matrix, which
#' defeats the larger-than-RAM design; the first release refuses above this
#' limit (gate A21). The same limit bounds the pool and observed comparison
#' matrices of the diagnostics, which are skipped above it.
#' @keywords internal
#' @noRd
QTL_REALISED_MAX_CELLS <- 2e7

#' Refuse what `anchor = "realised"` does not support (gate A6)
#' @keywords internal
#' @noRd
.dae_check_realised <- function(base_tbl, line_name, po, effects,
                                scale_to_target) {
  if (!is.null(effects) || !scale_to_target) {
    stop("`anchor` chooses what the calibration is exact for; manual ",
         "`effects` and `scale_to_target = FALSE` are not calibrated. Use the ",
         "default anchor.", call. = FALSE)
  }
  if (is.null(base_tbl)) {
    stop("anchor = \"realised\" needs `base_tbl` selecting the individuals ",
         "whose genotype covariance is the anchor, e.g. ",
         "base_tbl = get_table(pop, \"ind_meta\") |> filter(gen == 0L).",
         call. = FALSE)
  }
  if (base_tbl$table_name %in% c("founder_haplotypes", "ind_haplotype")) {
    stop("anchor = \"realised\" needs `base_tbl` to select individuals (a ",
         "table with id_ind such as ind_meta), not '", base_tbl$table_name,
         "'. ", if (base_tbl$table_name == "founder_haplotypes")
           "The founder pool has no individuals yet. "
         else "ind_haplotype selects allele copies, and a copy filter such as ",
         if (base_tbl$table_name == "ind_haplotype")
           "line_origin leaves crossbreds with partial genotypes, whose Cov(X) is not a population covariance. ",
         "Use anchor = \"genic\" for a pool or copy selection.", call. = FALSE)
  }
  if (!is.null(line_name) || any(!vapply(po, is.null, logical(1)))) {
    stop("anchor = \"realised\" supports the common scope only (no ",
         "`line_name`, no `parent_origin`): a scoped anchor needs the ",
         "covariance of the eligible copies, not of dosages. Use ",
         "anchor = \"genic\".", call. = FALSE)
  }
  invisible(TRUE)
}

#' Resolve the call's additive target (§6C)
#'
#' @param need `"target"` (calibrate), `"sigma"` (k >= 2 unscaled draw: the
#'   stored block is the draw covariance), or `"none"` (manual or unscaled
#'   single-trait effects: no target is read or required, and the block
#'   checks are skipped).
#' @return list(G = named k x k matrix or NULL, write = TRUE when `G` was
#'   passed and must be written with the terms).
#' @keywords internal
#' @noRd
.dae_resolve_target <- function(pop, trait_name, line_name, G,
                                trait_var_comp_tbl, need) {
  conn <- pop$db_conn
  if (need == "none") {
    if (!is.null(trait_var_comp_tbl)) {
      stop("`trait_var_comp_tbl` chooses a target to calibrate to, but these ",
           "effects are not calibrated.", call. = FALSE)
    }
    return(list(G = NULL, write = FALSE))
  }

  if (!is.null(G)) {
    G <- .check_cov_dimnames(G, trait_name, "G")
    if (anyNA(G) || any(!is.finite(G))) {
      stop("`G` must contain only finite values.", call. = FALSE)
    }
    .qtl_target_std(unname(G), name = "G")
    # A stored non-additive target is never silently ignored, whichever way
    # the additive target arrives. With `G` there is no trait_var_comp_tbl to
    # leave it out explicitly, so the explicit route is: store `G` first, then
    # select the additive rows.
    stored <- .dae_default_target_rows(conn, trait_name, line_name)
    other  <- setdiff(unique(stored$effect_name[
      stored$trait_name_1 %in% trait_name | stored$trait_name_2 %in% trait_name]),
      "additive")
    if (length(other)) {
      stop("A stored '", paste(other, collapse = "', '"), "' target exists for ",
           paste(trait_name, collapse = ", "), ". define_additive_effects() ",
           "calibrates the additive block only and never silently ignores a ",
           "stored target. To generate additive effects only, store `G` first ",
           "and select it explicitly:\n",
           "  define_effect_cov_matrix(pop, \"additive\", G, trait_name = c(",
           paste0('"', trait_name, '"', collapse = ", "), ")",
           if (!is.null(line_name)) paste0(', line_name = "', line_name, '"'),
           ")\n  ... |> define_additive_effects(..., trait_var_comp_tbl = ",
           'get_table(pop, "trait_var_comp") |>\n',
           '      dplyr::filter(effect_name == "additive", ',
           if (is.null(line_name)) "is.na(line_name)"
           else paste0('line_name == "', line_name, '"'), "))\n",
           "or remove the stored block with remove_rows().", call. = FALSE)
    }
    block <- .tvc_block_traits(conn, "additive", line_name, trait_name)
    if (length(block)) {
      stop("An 'additive' block is already stored for ",
           paste(block, collapse = ", "),
           if (is.null(line_name)) " (population-wide)"
           else paste0(" (line '", line_name, "')"),
           ". Stored targets are never overwritten, even by an identical ",
           "matrix. Drop `G` to calibrate to the stored block, or remove it ",
           "first:\n", .tvc_removal_call("additive", line_name, block),
           call. = FALSE)
    }
    return(list(G = (G + t(G)) / 2, write = TRUE))
  }

  explicit <- !is.null(trait_var_comp_tbl)
  rows <- if (explicit) {
    if (!inherits(trait_var_comp_tbl, "tidybreed_table") ||
        trait_var_comp_tbl$table_name != "trait_var_comp") {
      stop("`trait_var_comp_tbl` must be get_table(pop, \"trait_var_comp\") ",
           "|> filter(...).", call. = FALSE)
    }
    if (!identical(trait_var_comp_tbl$pop$db_conn, conn)) {
      stop("`trait_var_comp_tbl` must come from the same pop as `tbl`.",
           call. = FALSE)
    }
    dplyr::collect(trait_var_comp_tbl)
  } else {
    .dae_default_target_rows(conn, trait_name, line_name)
  }
  rows <- rows[rows$trait_name_1 %in% trait_name |
               rows$trait_name_2 %in% trait_name, , drop = FALSE]
  hint_filter <- paste0(
    '  trait_var_comp_tbl = get_table(pop, "trait_var_comp") |>\n',
    '    dplyr::filter(effect_name == "additive"',
    if (is.null(line_name)) ", is.na(line_name)"
    else paste0(', line_name == "', line_name, '"'),
    ', trait_name_1 %in% c(', paste0('"', trait_name, '"', collapse = ", "),
    '),\n                  trait_name_2 %in% c(',
    paste0('"', trait_name, '"', collapse = ", "), '))')

  # A stored target is never silently ignored.
  other <- setdiff(unique(rows$effect_name), "additive")
  if (length(other)) {
    stop("A stored '", paste(other, collapse = "', '"), "' target exists for ",
         paste(trait_name, collapse = ", "), ". define_additive_effects() ",
         "calibrates the additive block only and never silently ignores a ",
         "stored target. To generate additive effects only, say so with\n",
         hint_filter, "\nor remove the stored block with remove_rows().",
         call. = FALSE)
  }
  outside <- setdiff(unique(c(rows$trait_name_1, rows$trait_name_2)), trait_name)
  if (length(outside)) {
    block <- sort(unique(c(trait_name, outside)))
    stop("The stored 'additive' block links ",
         paste(trait_name, collapse = ", "), " with ",
         paste(outside, collapse = ", "), ". Calibrating ",
         paste(trait_name, collapse = ", "), " alone would break the stored ",
         "covariance. Either pass all of them, trait_name = c(",
         paste0('"', block, '"', collapse = ", "), "), or calibrate ",
         paste(trait_name, collapse = ", "), " alone on purpose with\n",
         hint_filter, call. = FALSE)
  }
  if (explicit && length(unique(rows$line_name)) > 1L) {
    stop("`trait_var_comp_tbl` holds two candidate 'additive' blocks for ",
         paste(trait_name, collapse = ", "), " (line_name ",
         paste(ifelse(is.na(unique(rows$line_name)), "NULL",
                      unique(rows$line_name)), collapse = " and "),
         "). Filter it to one.", call. = FALSE)
  }
  if (nrow(rows) == 0L) {
    stop("No 'additive' target is stored for ", paste(trait_name, collapse = ", "),
         if (!is.null(line_name)) paste0(" (line '", line_name, "' or population-wide)"),
         if (explicit) " in `trait_var_comp_tbl`", ". Pass `G =`, or store one with ",
         "define_effect_cov_matrix(pop, \"additive\", ...).", call. = FALSE)
  }
  k <- length(trait_name)
  M <- matrix(NA_real_, k, k, dimnames = list(trait_name, trait_name))
  key <- paste(rows$trait_name_1, rows$trait_name_2)
  if (anyDuplicated(key)) {
    stop("The 'additive' rows for ", paste(trait_name, collapse = ", "),
         " hold a (trait_name_1, trait_name_2) pair twice. Filter to one block.",
         call. = FALSE)
  }
  M[cbind(rows$trait_name_1, rows$trait_name_2)] <- rows$cov_value
  if (anyNA(M)) {
    miss <- which(is.na(M), arr.ind = TRUE)
    stop("The 'additive' rows for ", paste(trait_name, collapse = ", "),
         " do not form one complete k x k block: missing (",
         paste0(trait_name[miss[, 1]], ", ", trait_name[miss[, 2]],
                collapse = "), ("), "). Both (i, j) and (j, i) are needed.",
         call. = FALSE)
  }
  if (!isSymmetric(unname(M), tol = 1e-12)) {
    stop("The stored 'additive' block for ", paste(trait_name, collapse = ", "),
         " is not symmetric.", call. = FALSE)
  }
  .qtl_target_std(unname(M), name = "the stored 'additive' block")
  list(G = M, write = FALSE)
}

#' Default target rows: the call's line block, else population-wide
#'
#' Every genetic `effect_name` is resolved on its own: line C can have its own
#' `additive` block and share the population-wide `dominance` one.
#' @keywords internal
#' @noRd
.dae_default_target_rows <- function(conn, trait_name, line_name) {
  out <- NULL
  for (e in GENETIC_EFFECT_NAMES) {
    ln <- .tvc_resolve_line(conn, e, trait_name, line_name)
    r <- DBI::dbGetQuery(conn, paste0(
      "SELECT effect_name, line_name, trait_name_1, trait_name_2, cov_value ",
      "FROM trait_var_comp WHERE effect_name = ", DBI::dbQuoteLiteral(conn, e),
      " AND ", .tvc_line_sql(conn, ln)))
    out <- rbind(out, r)
  }
  out
}

#' Draw the effect architecture B0 (n_loci x k, NA off each trait's mask)
#'
#' The one sampler both generators call first (§7.2, gate C4): today's draws,
#' unchanged. One trait: `rnorm()` or the signed gamma at the masked loci in
#' `locus_id` order. Several: one `MVN(0, G)` row per locus any trait uses,
#' then masked per trait (`method = "union"`).
#' @keywords internal
#' @noRd
.draw_additive_architecture <- function(mask, distribution, G) {
  n_loci <- nrow(mask); k <- ncol(mask)
  B <- matrix(NA_real_, n_loci, k, dimnames = dimnames(mask))
  if (k == 1L) {
    n <- sum(mask[, 1L])
    B[mask[, 1L], 1L] <- switch(distribution,
      normal = stats::rnorm(n),
      gamma  = stats::rgamma(n, shape = 0.4, rate = 1.66) *
                 sample(c(-1, 1), n, replace = TRUE))
    return(B)
  }
  any_qtl <- apply(mask, 1L, any)
  draws <- MASS::mvrnorm(n = sum(any_qtl), mu = rep(0, k), Sigma = G)
  if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1L)
  draws[!mask[any_qtl, , drop = FALSE]] <- NA_real_
  B[any_qtl, ] <- draws
  B
}

#' Collect the base individuals' dosages at the QTL, `id_ind` then `locus_id`
#'
#' Integer sums only (exact). Refuses an empty or single-individual base, a
#' base above `QTL_REALISED_MAX_CELLS`, and any individual without both
#' copies at every QTL: `extract_genotypes()` would read a missing copy as 0,
#' and `Cov(X)` of partial genotypes is not a population covariance.
#' @return list(X = n x m integer-valued matrix, id_ind, locus_id).
#' @keywords internal
#' @noRd
.dae_collect_dosages <- function(pop, base_tbl, locus_ids) {
  conn <- pop$db_conn
  sub  <- as.character(dbplyr::sql_render(base_tbl$tbl))
  ids_sql <- paste0("SELECT DISTINCT id_ind FROM (", sub, ") b")
  n <- DBI::dbGetQuery(conn, paste0("SELECT COUNT(*) AS n FROM (", ids_sql, ")"))$n
  m <- length(locus_ids)
  if (n < 2L) {
    stop("anchor = \"realised\" needs at least 2 individuals in `base_tbl`; ",
         "it selects ", n, ".", call. = FALSE)
  }
  if (as.numeric(n) * m > QTL_REALISED_MAX_CELLS) {
    stop("anchor = \"realised\" would collect a ", n, " x ", m,
         " genotype matrix (", format(as.numeric(n) * m, big.mark = ","),
         " cells), above the limit of ",
         format(QTL_REALISED_MAX_CELLS, big.mark = ",", scientific = FALSE),
         ". Use anchor = \"genic\", or a smaller `base_tbl` or QTL set.",
         call. = FALSE)
  }
  lst <- paste(as.integer(locus_ids), collapse = ", ")
  d <- DBI::dbGetQuery(conn, paste0(
    "SELECT h.id_ind, h.locus_id, CAST(SUM(h.allele) AS INTEGER) AS dosage, ",
    "COUNT(*) AS n_copies FROM ind_haplotype h ",
    "JOIN (", ids_sql, ") ids USING (id_ind) ",
    "WHERE h.locus_id IN (", lst, ") ",
    "GROUP BY h.id_ind, h.locus_id ORDER BY h.id_ind, h.locus_id"))
  per_ind <- table(factor(d$id_ind))
  ids <- sort(unique(d$id_ind))
  bad <- length(ids) < n || any(per_ind != m) || any(d$n_copies != 2L)
  if (bad) {
    stop("anchor = \"realised\": some individuals in `base_tbl` lack both ",
         "allele copies at every QTL (e.g. crossbreds selected by line, or ",
         "individuals without haplotypes). Cov(X) of partial genotypes is not ",
         "a population covariance. Use anchor = \"genic\".", call. = FALSE)
  }
  X <- matrix(as.numeric(d$dosage), nrow = length(ids), ncol = m, byrow = TRUE,
              dimnames = list(ids, NULL))
  list(X = X, id_ind = ids, locus_id = sort(as.integer(locus_ids)))
}

#' Compare the delivered covariance with what another population sees (§7.4)
#'
#' Exact under the chosen anchor does not mean exact everywhere. Computes the
#' relative spectrum of what the **other** population sees against the target:
#' * genic anchor, founder-pool base: the pool expectation of a founder's
#'   dosage covariance, `2 Cov(H)` with the **population** divisor `n_h`.
#'   `add_founders()` draws each haplotype independently and with replacement
#'   from the pool, so a founder's two copies are iid draws from the empirical
#'   haplotype distribution, whose covariance divides by `n_h`, not `n_h - 1`
#'   (the pool keeps its LD; the founders actually drawn scatter around this);
#' * genic anchor, individuals base: their observed `Cov(X)` (`n - 1`);
#' * realised anchor: the genic (random-mating) limit at their frequencies.
#' An `ind_haplotype` copy selection has no pairing to compare with and is
#' skipped, as is anything above `QTL_REALISED_MAX_CELLS` (with a message).
#'
#' Runs **before** the commit, so its queries cannot fail after the target
#' and terms are written; `.dae_report_diagnostics()` emits the result after
#' the commit. The pool's identity columns (`line_name`, `haplotype_id`) are
#' read from the physical table with `base_tbl`'s filters re-applied, so a
#' `select()`ed selection works. Nothing is stored.
#' @return NULL (nothing to report) or list(spec, label, pool, n_h, note).
#' @keywords internal
#' @noRd
.dae_diagnostics <- function(pop, base_tbl, anchor, genome_order, rows,
                             B, std, p_base) {
  if (std$rank == 0L) return(NULL)
  conn <- pop$db_conn
  locus_ids <- genome_order$locus_id[rows]
  pool <- FALSE
  n_h  <- NA_integer_
  if (anchor == "realised") {
    p <- p_base[rows]
    cand  <- crossprod(B * sqrt(2 * p * (1 - p)))
    label <- "genic limit (expectation)"
  } else if (base_tbl$table_name == "founder_haplotypes") {
    lazy <- dplyr::tbl(conn, "founder_haplotypes")
    if (length(base_tbl$pending_filter) > 0L) {
      lazy <- dplyr::filter(lazy, !!!base_tbl$pending_filter)
    }
    sub <- as.character(dbplyr::sql_render(lazy))
    n_h <- DBI::dbGetQuery(conn, paste0(
      "SELECT COUNT(*) AS n FROM (SELECT DISTINCT line_name, haplotype_id ",
      "FROM (", sub, ") b)"))$n
    if (n_h < 2L) return(NULL)
    if (as.numeric(n_h) * length(locus_ids) > QTL_REALISED_MAX_CELLS) {
      return(list(note = paste0(
        "Skipped the pool-expectation comparison: ", n_h, " haplotypes x ",
        length(locus_ids), " QTL is above the in-memory limit.")))
    }
    h <- DBI::dbGetQuery(conn, paste0(
      "SELECT COALESCE(b.line_name, '') || '|' || CAST(b.haplotype_id AS VARCHAR) ",
      "AS hap, gm.locus_id, CAST(b.allele AS DOUBLE) AS allele FROM (", sub,
      ") b JOIN genome_meta gm ON b.locus_name = gm.locus_name WHERE gm.locus_id IN (",
      paste(as.integer(locus_ids), collapse = ", "), ") ORDER BY hap, gm.locus_id"))
    haps <- sort(unique(h$hap))
    H <- matrix(NA_real_, length(haps), length(locus_ids))
    H[cbind(match(h$hap, haps), match(h$locus_id, locus_ids))] <- h$allele
    if (anyNA(H)) return(NULL)
    Hc <- sweep(H, 2L, colMeans(H), "-")
    cand  <- 2 * crossprod(Hc %*% B) / nrow(H)
    label <- "founder pool (pool expectation under random pairing)"
    pool  <- TRUE
  } else if (base_tbl$table_name == "ind_haplotype") {
    return(NULL)
  } else {
    X <- tryCatch(.dae_collect_dosages(pop, base_tbl, locus_ids)$X,
                  error = function(e) NULL)
    if (is.null(X)) {
      return(list(note = paste0(
        "Skipped the observed comparison: the base individuals' genotypes ",
        "could not be collected (partial or above the in-memory limit).")))
    }
    Xc <- sweep(X, 2L, colMeans(X), "-")
    cand  <- crossprod(Xc %*% B) / (nrow(X) - 1)
    label <- "base individuals (observed)"
  }
  list(spec = .qtl_relative_spectrum_std(std, cand), label = label,
       pool = pool, n_h = n_h)
}

#' Report a `.dae_diagnostics()` result
#'
#' A pool's departure is its sampling LD, not a mistake in the call (Q22):
#' always a message, never a warning. The observed and genic-limit
#' comparisons warn outside `warn_bounds`.
#' @keywords internal
#' @noRd
.dae_report_diagnostics <- function(res, anchor, warn_bounds) {
  if (!is.null(res$note)) {
    message(res$note)
    return(invisible(NULL))
  }
  spec <- res$spec
  outside <- min(spec) < warn_bounds[1] || max(spec) > warn_bounds[2]
  if (res$pool) {
    message("The ", res$label, " sees relative spectrum ",
            .qtl_spectrum_text(spec), " of the target (sampling LD of ",
            res$n_h, " haplotypes).",
            if (outside) paste0(
              " For the target exactly in the founders, add them first and ",
              "call again with anchor = \"realised\" and base_tbl = those ",
              "individuals."))
  } else if (outside) {
    warning("The calibration is exact under the ", anchor, " anchor, but the ",
            res$label, " sees a covariance departing from the target: relative ",
            "spectrum ", .qtl_spectrum_text(spec), " outside warn_bounds [",
            warn_bounds[1], ", ", warn_bounds[2], "]. warn_bounds = NULL ",
            "turns this off.", call. = FALSE)
  }
  invisible(spec)
}

#' Compact text for a covariance matrix in a message
#' @keywords internal
#' @noRd
.dae_format_cov <- function(M) {
  nm <- rownames(M)
  if (length(nm) == 1L) return(sprintf("variance %.6g", M[1, 1]))
  ut <- which(upper.tri(M, diag = TRUE), arr.ind = TRUE)
  paste0("[", paste(sprintf("%s,%s = %.6g", nm[ut[, 1]], nm[ut[, 2]], M[ut]),
                    collapse = "; "), "]")
}

#' Compact text for the correlations of a covariance matrix
#' @keywords internal
#' @noRd
.dae_format_cor <- function(M) {
  nm <- rownames(M)
  d  <- sqrt(pmax(diag(M), 0))
  ut <- which(upper.tri(M), arr.ind = TRUE)
  r  <- M[ut] / (d[ut[, 1]] * d[ut[, 2]])
  paste(sprintf("%s,%s = %.3f", nm[ut[, 1]], nm[ut[, 2]], r), collapse = "; ")
}


#' Resolve `base_tbl` (default or explicit) into `p_base`
#'
#' Called after argument validation in each path, so an argument error is
#' reported before the base is touched or warned about. The default is the
#' population the effect applies to (`.dae_default_base()`) and may warn; an
#' explicit `base_tbl` is an intentional selection and never warns. `p_base`
#' may hold `NA` at loci the base has no copies for -- `.dae_require_base_at()`
#' refuses those at the QTL before anything folds `p` into `V_A`.
#'
#' @keywords internal
#' @noRd
.dae_resolve_base <- function(pop, base_tbl, line_name) {
  if (is.null(base_tbl)) {
    base_tbl <- .dae_default_base(pop, line_name)
  } else {
    .validate_base_tbl(base_tbl, pop)
  }
  list(p_base = extract_allele_freq(base_tbl)$allele_freq,
       label  = .dae_base_label(base_tbl),
       tbl    = base_tbl)
}

#' Resolve the default base population for `define_additive_effects()`
#'
#' The founder pool of the line the effect applies to, with the same
#' `line -> NULL` precedence as `resolve_genome_map()` and
#' `resolve_chr_inheritance()`: a named pool wins, the shared (`NULL`) pool is
#' the fallback, and only when neither exists is it an error. This is
#' pool-level fallback, never per-locus stitching -- a named pool that is
#' missing loci is the base as it stands, and `.dae_require_base_at()` refuses
#' the gap at the QTL.
#'
#' Warns only for a population-wide effect on a founder table holding more
#' than one pool (Wahlund; see the roxygen). A line-scoped effect that falls
#' back to the shared pool does not warn: a shared pool is one population by
#' construction, so nothing is pooled.
#'
#' @param pop A `tidybreed_pop`.
#' @param line_name The effect's `line_name`, or `NULL`.
#' @return A `tidybreed_table` on `founder_haplotypes`, filtered as resolved.
#' @keywords internal
#' @noRd
.dae_default_base <- function(pop, line_name) {
  if (!"founder_haplotypes" %in% pop$tables) {
    stop("founder_haplotypes table not found. Did you call ",
         "define_founder_haplotypes()? Otherwise pass base_tbl explicitly, ",
         "e.g. base_tbl = get_table(pop, \"ind_meta\").", call. = FALSE)
  }
  fh   <- get_table(pop, "founder_haplotypes")
  conn <- pop$db_conn
  if (is.null(line_name)) {
    # NULL counts as its own pool; COUNT(DISTINCT line_name) would ignore it.
    n_pools <- DBI::dbGetQuery(conn,
      "SELECT COUNT(DISTINCT COALESCE(line_name, '')) AS n FROM founder_haplotypes")$n
    if (isTRUE(n_pools > 1L)) {
      warning(
        "founder_haplotypes holds ", n_pools, " pools; base allele frequencies ",
        "for this population-wide effect are pooled across all of them, which ",
        "overstates within-line heterozygosity. Set line_name for per-line ",
        "centering, or pass base_tbl explicitly to pool on purpose.",
        call. = FALSE)
    }
    return(fh)
  }
  has_line <- DBI::dbGetQuery(conn,
    "SELECT EXISTS(SELECT 1 FROM founder_haplotypes WHERE line_name = ?) AS ok",
    params = list(line_name))$ok
  if (isTRUE(has_line)) {
    return(dplyr::filter(fh, .data$line_name == .env$line_name))
  }
  has_shared <- DBI::dbGetQuery(conn,
    "SELECT EXISTS(SELECT 1 FROM founder_haplotypes WHERE line_name IS NULL) AS ok")$ok
  if (isTRUE(has_shared)) {
    return(dplyr::filter(fh, is.na(.data$line_name)))
  }
  .founder_base_empty_error(conn, paste0("line '", line_name, "'"))
}

#' Refuse a QTL locus the base has no allele copies for
#'
#' `p_base` may hold `NA` where `extract_allele_freq()` found no copies. The
#' genic anchor folds `p` into its weights `n_eligible * p q`, so an `NA`
#' reaching it would make every effect `NA`, and treating a missing `p` as `0`
#' would centre the locus silently. Either is worse than stopping here and naming the loci.
#'
#' @param p_base Numeric, length `n_loci`, `locus_id` order; may hold `NA`.
#' @param written_tf Logical mask of the loci this call will write.
#' @param locus_names Character, length `n_loci`, `locus_id` order.
#' @param trait Optional trait name for the message (multi-trait path).
#' @keywords internal
#' @noRd
.dae_require_base_at <- function(p_base, written_tf, locus_names, trait = NULL) {
  bad <- written_tf & is.na(p_base)
  if (!any(bad)) return(invisible(TRUE))
  shown <- locus_names[bad]
  if (length(shown) > 5L) {
    shown <- c(shown[1:5], paste0("... (", length(locus_names[bad]), " total)"))
  }
  stop("base_tbl has no allele copies at ", sum(bad), " selected QTL ",
       if (!is.null(trait)) paste0("for trait '", trait, "' ") else "",
       "(", paste(shown, collapse = ", "), "). ",
       "Widen the base selection or drop those loci from the QTL filter.",
       call. = FALSE)
}

#' Human-readable base description for the completion message
#'
#' @keywords internal
#' @noRd
.dae_base_label <- function(base_tbl) {
  n_filters <- length(base_tbl$pending_filter)
  paste0(base_tbl$table_name,
         if (n_filters > 0L) paste0(" [", n_filters, " filter",
                                    if (n_filters > 1L) "s", "]") else "")
}

#' Generate phenotype records for a subset of individuals
#'
#' @description
#' Simulates phenotype values for one or more phenotypes and writes them to
#' `ind_phenotype`. Also computes and stores the underlying true breeding
#' value (TBV) per individual per trait in `ind_tbv`.
#'
#' **Model** (per phenotype, on the liability / continuous scale):
#' \preformatted{
#'   y_i = mean + sum(fixed_shifts) + sum(random_shifts) + TBV_i + e_i
#' }
#'
#' * `mean` comes from `phenotype_meta.mean`.
#' * Fixed shifts come from the `fixed_class` / `fixed_cov` rows of
#'   `phenotype_effects`.
#' * Random shifts come from its `random` rows (see [define_effect_random()]):
#'   one draw per distinct level of the effect's `source_column`, realized
#'   the first time any record touches the level and stored in
#'   `phenotype_random_effects`, then reused by every later record with that
#'   level — in this call or any later one. A level is persistent: a pen
#'   that is re-realized per batch is a different level (`pen_batch`), not a
#'   different feature. When the effect is correlated across phenotypes
#'   ([define_effect_cov_matrix()] with the effect's name), a level's draw
#'   for one phenotype is conditional on the draws it already has stored for
#'   the block's other phenotypes — `add_phenotype("ADG")` today and
#'   `add_phenotype("BF")` next season gives pen `P1` a `(ADG, BF)` pair with
#'   the declared covariance, whichever came first. A record whose level is
#'   `NULL` gets no draw and a shift of `0`.
#' * For **simple** phenotypes (`phenotype_name == trait_name`), `TBV_i` is the
#'   standard additive TBV from `genome_effects` (computed via [add_tbv()],
#'   which this function calls internally for every source trait it needs).
#' * For **composite** phenotypes (rows in `phenotype_components`, written by
#'   `define_phenotype(..., components = ...)`), `TBV_i` is the weighted sum
#'   of contributor TBVs (self, dam, sire, or group) — see
#'   `.assemble_composite_tbv()`.
#' * For **`formula_tbv`** composite phenotypes (`phenotype_meta.formula_tbv`
#'   set, written by `define_phenotype(..., formula_tbv = ...)`), `TBV_i` is
#'   evaluated from a small DSL expression referencing self/dam/sire/group
#'   TBVs instead of a `phenotype_components` data frame.
#' * `e_i` is the residual, drawn from the phenotype's residual covariance
#'   block in `phenotype_var_comp` (see [define_residual_cov()]). Within a
#'   block, residuals are **correlated across phenotypes and sequential in
#'   time**: each record is drawn from the block's multivariate normal
#'   conditional on every residual the same individual has already realized
#'   for the other phenotypes of the block at the same `pheno_number` —
#'   whether those were drawn in this call, in an earlier call, or supplied
#'   through `user_residual`. So `add_phenotype("A")` today and
#'   `add_phenotype("B")` after culling gives the survivors' `B` residual
#'   the stored correlation with their `A` residual, with no `B` record for
#'   the culled. Heterogeneous residual variance (strata by
#'   `condition_column`) is applied per record: each record draws from the
#'   stratum its condition value selects, falling back to the unconditional
#'   `R` (with `residual_condition_level` stored as `NULL`) when the value
#'   is `NULL` or matches no stratum, and erroring if there is no
#'   unconditional stratum to fall back on. A stored residual drawn under a
#'   different stratum than the current record resolves to is an error, or
#'   is dropped from the conditioning set with a warning when the block's
#'   phenotypes have `condition_change_action = "independent"`. The realized
#'   residual and its stratum are written to `ind_phenotype.residual_value`
#'   / `residual_condition_level` (liability scale; `NULL` for `user_values`
#'   and `derived_formula` records). `pheno_number` pairs records
#'   ordinally, not by simulated time; repeated records of the *same*
#'   phenotype have independent residuals (use a permanent-environment
#'   random effect for within-animal covariance).
#'
#' **`derived_formula` phenotypes are the one exception to the model above.**
#' When `phenotype_meta.type == "derived_formula"` (`phenotype_meta.formula`
#' set), the phenotype value is computed directly as an arithmetic expression
#' over other individuals' already-written `ind_phenotype` records — there is
#' no TBV, no mean/fixed/random contribution, and no residual draw for that
#' phenotype. When a call mixes `derived_formula` phenotypes with others that
#' feed them, the phenotypes are topologically sorted first so dependencies
#' are written before the formulas that consume them.
#'
#' **Subset selection**: pipe a `tidybreed_table` (from [get_table()] and
#' optionally [dplyr::filter()]) as the first argument.
#'
#' **How a call runs.** Every record is planned first — sex expression, the
#' repeatable guard, fixed-effect skips, missing-component exclusions and
#' `pheno_number` are all decided before a single random number is drawn —
#' then every draw is made in memory, then everything is written in one
#' transaction. An individual that ends up without a record therefore never
#' consumes RNG, and seeded output does not depend on physical row order:
#' records are planned in `id_ind` order within each phenotype. See
#' `?add_phenotype_stages`.
#'
#' **If a call fails**, nothing is written: `ind_phenotype` and
#' `phenotype_random_effects` are exactly as they were, whether the error
#' came from validation, from sampling (a missing variance, a residual
#' stratum change under `condition_change_action = "error"`) or from the
#' write itself — down to the schema, so a new column named in `...` is
#' rolled back with the rows it was added for. The random-number stream is **not** rewound: `.Random.seed`
#' stays advanced by the draws made before the error, as after any other
#' failed R call, so re-running the call draws different values. Pass `seed`
#' (or call `set.seed()`) again if the retry must reproduce the failed call.
#' The TBVs the call materialized through [add_tbv()] remain; they do not
#' depend on the RNG, and the retry rewrites them.
#'
#' **Escape hatches**:
#' * `user_values`: skip model computation and write these values as phenotype
#'   records for the subset.
#' * `user_residual`: supply the residuals of some or all phenotypes instead
#'   of drawing them; the rest are drawn conditional on the supplied values.
#'
#' @param tbl A `tidybreed_table` from [get_table()], optionally piped through
#'   [dplyr::filter()]. Any table with an `id_ind` column is accepted; the
#'   individuals acted on are the distinct `id_ind` values present in the
#'   (filtered) table. An unfiltered `ind_meta` selects every individual; an
#'   unfiltered `ind_ebv`, `ind_index`, `ind_genotype`, ... selects only the
#'   individuals that have rows there. A table without `id_ind` is an error.
#' @param phenotype_name Character vector of phenotype name(s). When `NULL`
#'   (default), all phenotypes in `phenotype_meta` are used in
#'   `id_phenotype_meta` order.
#' @param user_residual Optional residuals to use instead of drawing them
#'   (the mean, covariate and TBV contributions are still computed and
#'   added). When exactly one phenotype in the call is generated from the
#'   model, a plain numeric vector matched **by position** to that
#'   phenotype's planned records — sorted `id_ind` order after sex
#'   expression, the repeatable guard and any exclusion, so its length must
#'   equal the planned record count, which may be smaller than the filtered
#'   `tbl`. Otherwise a named list keyed by `phenotype_name` that may name
#'   **any subset** of the model-generated phenotypes, each element following
#'   the same positional rule; the phenotypes not named are drawn conditional
#'   on the supplied values. Supplied residuals are stored in
#'   `residual_value` like drawn ones and condition later calls; a phenotype
#'   whose residuals are all supplied needs no residual variance declared.
#'   A value outside the support of a singular covariance (e.g. non-zero
#'   for a zero-variance phenotype) is an error. Named (per-`id_ind`) vectors are
#'   **not** supported here; cannot be combined with `user_values`.
#' @param user_values Optional override for the full phenotype value —
#'   skips the model entirely (mean, covariates, and residual are not
#'   evaluated), though TBVs are still computed and stored in `ind_tbv`. For a
#'   single `phenotype_name`: a plain numeric vector matching the planned
#'   records by position (sorted `id_ind` order, after sex expression and the
#'   repeatable guard), or a named numeric vector (e.g.
#'   `c(id_1 = 555, id_2 = 560)`) to match by `id_ind` regardless of order —
#'   every name must be an individual in that planned set, each once, and
#'   only the named individuals receive a record. For multiple phenotypes: a
#'   named list keyed by `phenotype_name`, each element following the same
#'   (positional-or-named) rule.
#' @param seed Optional integer for reproducibility.
#' @param ... Optional scalar extra columns written to `ind_phenotype`
#'   (broadcast to all records). Supply per-record vectors with
#'   [mutate_table()] after the call.
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_phenotype()], [define_trait()], [define_additive_effects()],
#'   [define_residual_cov()], [define_effect_cov_matrix()],
#'   [define_effect_fixed_class()], [define_effect_fixed_cov()],
#'   [define_effect_random()], [add_tbv()]
#'
#' @examples
#' \dontrun{
#' # All individuals — all phenotypes
#' pop <- pop |> get_table("ind_meta") |> add_phenotype()
#'
#' # Named phenotype, filtered subset
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(sex == "F", gen == 1L) |>
#'   add_phenotype("ADG")
#'
#' # Any table with id_ind chooses the individuals. Marker-assisted
#' # pre-selection: only carriers of two copies at Locus_10 (run add_dosage()
#' # first, since ind_genotype is an on-demand cache)
#' pop <- pop |>
#'   get_table("ind_genotype") |>
#'   dplyr::filter(locus_name == "Locus_10", dosage_value == 2L) |>
#'   add_phenotype("ADG")
#'
#' # EBV-based: only animals above an EBV threshold in the latest evaluation
#' pop <- pop |>
#'   get_table("ind_ebv") |>
#'   dplyr::filter(trait_name == "ADG", eval_number == 3L, ebv_value > 0.5) |>
#'   add_phenotype("ADG")
#'
#' # Unfiltered ind_ebv means "every animal that has an EBV", not everyone
#' pop <- pop |> get_table("ind_ebv") |> add_phenotype("ADG")
#'
#' # Composite (maternal) phenotype: WW = direct (self) + maternal (dam) TBV,
#' # registered once via define_phenotype(components = ...)
#' pop <- pop |>
#'   define_phenotype("WW", type = "continuous", mean = 230, residual_var = 180,
#'     components = tibble::tribble(
#'       ~source_trait_name, ~contributor_type,
#'       "WWD",              "self",
#'       "WWM",              "dam"
#'     ))
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 1L) |>
#'   add_phenotype("WW")
#'
#' # Culling between records. A and B share a residual covariance block, so
#' # B's residual is drawn conditional on each survivor's stored A residual.
#' # The culled animals simply get no B record -- and because a residual is
#' # realized only when a record is planned, they consume no draw either.
#' pop <- pop |>
#'   define_phenotype("A", type = "continuous", mean = 100) |>
#'   define_phenotype("B", type = "continuous", mean = 250) |>
#'   define_residual_cov(c("A", "B"),
#'     matrix(c(40, 18, 18, 30), 2, 2,
#'            dimnames = list(c("A", "B"), c("A", "B"))))
#'
#' # 1. Record A on everyone
#' pop <- pop |> get_table("ind_meta") |> add_phenotype("A", seed = 1)
#'
#' # 2. Cull on the realized A: keep the top half
#' cut <- pop |> get_table("ind_phenotype") |>
#'   dplyr::filter(phenotype_name == "A") |> dplyr::pull(pheno_value) |>
#'   stats::median()
#'
#' # 3. Record B on the survivors only. Selecting from ind_phenotype means
#' #    "the animals with this A record", not everyone.
#' pop <- pop |>
#'   get_table("ind_phenotype") |>
#'   dplyr::filter(phenotype_name == "A", pheno_value >= cut) |>
#'   add_phenotype("B", seed = 2)
#'
#' # Checking the result: do NOT expect the correlation between the two
#' # stored residual_value columns to equal the declared 18/sqrt(40*30) =
#' # 0.52. The survivors were selected on A, so their A residuals are
#' # range-restricted and the observed correlation is attenuated (~0.35 for
#' # a top-half cull). What selection does *not* change is the conditional
#' # slope, so that is the quantity to check:
#' #
#' #   coef(lm(residual_B ~ residual_A))[2]  ==  18 / 40  ==  0.45
#'
#' # Escape hatch: supply phenotype values directly (skips the model, but
#' # still computes and stores TBVs); named vector matches by id_ind
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 0L, sex == "M") |>
#'   add_phenotype("ADG", user_values = c(A_1 = 555, A_2 = 560))
#' }
#' @export
add_phenotype <- function(tbl,
                          phenotype_name = NULL,
                          user_residual  = NULL,
                          user_values   = NULL,
                          seed          = NULL,
                          ...) {

  stopifnot(inherits(tbl, "tidybreed_table"))
  pop <- tbl$pop
  validate_tidybreed_pop(pop)

  # Resolve phenotype names from phenotype_meta
  if (is.null(phenotype_name)) {
    phenotype_name <- DBI::dbGetQuery(
      pop$db_conn,
      "SELECT phenotype_name FROM phenotype_meta ORDER BY id_phenotype_meta"
    )$phenotype_name
    if (length(phenotype_name) == 0L)
      stop("No phenotypes found in phenotype_meta. ",
           "Call define_phenotype() first.", call. = FALSE)
  }

  stopifnot(is.character(phenotype_name), length(phenotype_name) >= 1)
  lapply(phenotype_name, validate_sql_identifier, what = "phenotype name")
  phenotype_name <- unique(phenotype_name)

  extra_cols <- list(...)
  if (length(extra_cols) > 0) {
    # Check for common misspellings of phenotype_name
    misspelled <- c("trait", "phenotype", "trait_name", "pheno_name")
    for (mis in misspelled) {
      if (mis %in% names(extra_cols)) {
        stop("Did you mean 'phenotype_name'? Found argument '",
             mis, "' in ...", call. = FALSE)
      }
    }

    for (nm in names(extra_cols)) {
      if (length(extra_cols[[nm]]) != 1L) {
        stop("Custom field '", nm,
             "' in add_phenotype() must be a scalar ",
             "(broadcast to all records). Supply per-record vectors ",
             "with mutate_table() after the call.", call. = FALSE)
      }
    }
  }

  if (!is.null(user_values) && !is.null(user_residual)) {
    stop("user_values and user_residual cannot be combined: user_values ",
         "bypasses the model, so there is no residual to fix.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  # Three stages (see ?add_phenotype_stages): PLAN decides every record with
  # no RNG and no writes; RESOLVE draws everything in memory; COMMIT writes
  # once, in one transaction, without touching the RNG.
  plan <- .ap_plan(tbl, phenotype_name, user_values = user_values)
  if (is.null(plan)) return(invisible(pop))
  pop <- plan$pop

  resolved <- .ap_resolve(plan, user_residual = user_residual)
  .ap_commit(pop, plan, resolved, extra_cols)

  invisible(pop)
}


# ── Internal helpers ───────────────────────────────────────────────────────────

#' Assemble the composite TBV of one phenotype from `phenotype_components`
#'
#' Sums `weight * contributor TBV` over the phenotype's component rows, one
#' contributor lookup per row (see `?contributor_tbv`). `ind_tbv` must
#' already hold the source traits for every contributor
#' (`.ap_materialize_tbvs()`). A missing piece — a `NULL` dam or sire, a
#' contributor with no TBV, a `NULL` group value, a `NULL` covariate —
#' makes the individual's composite `NA`.
#'
#' @param comp_rows The phenotype's `phenotype_components` rows.
#' @param subset_df The planned `ind_meta` rows (needs `id_ind`,
#'   `id_parent_1`, `id_parent_2`).
#' @param missing_component_action `"skip"` warns and returns `NA` for the
#'   excluded individuals; `"error"` stops. Both name the count and up to
#'   five ids.
#' @return Numeric vector named by `id_ind`; `NA` marks an excluded
#'   individual.
#' @keywords internal
.assemble_composite_tbv <- function(pop, phenotype_name, comp_rows, subset_df,
                                    missing_component_action) {
  conn      <- pop$db_conn
  focal_ids <- as.character(subset_df$id_ind)
  n         <- length(focal_ids)
  composite <- rep(0, n)

  for (i in seq_len(nrow(comp_rows))) {
    comp  <- comp_rows[i, , drop = FALSE]
    trait <- comp$source_trait_name
    what  <- paste0("Phenotype '", phenotype_name, "', component '", trait,
                    "' (", comp$contributor_type, ")")

    weight <- if (is.na(comp$weight)) 1 else comp$weight
    wt <- switch(
      comp$weight_type,
      fixed = rep(weight, n),
      covariate = {
        if (is.na(comp$covariate_name)) {
          stop(what, ": weight_type = 'covariate' needs covariate_name.",
               call. = FALSE)
        }
        cov_tbl <- if (is.na(comp$covariate_table)) "ind_meta" else comp$covariate_table
        weight * as.numeric(.read_one_per_id(conn, cov_tbl, comp$covariate_name,
                                             focal_ids, what))
      },
      stop(what, ": weight_type '", comp$weight_type, "' is not implemented; ",
           "use 'fixed' or 'covariate'.", call. = FALSE))

    vec <- switch(
      comp$contributor_type,
      self  = .tbv_by_id(conn, trait, focal_ids),
      dam   = .tbv_by_id(conn, trait, subset_df$id_parent_2),
      sire  = .tbv_by_id(conn, trait, subset_df$id_parent_1),
      group = .group_mate_tbv(conn, trait, focal_ids, comp$group_column,
                              comp$group_table, comp$aggregation, what))
    composite <- composite + wt * vec
  }

  n_missing <- sum(is.na(composite))
  if (n_missing > 0L) {
    missing_ids <- focal_ids[is.na(composite)]
    msg <- paste0(
      n_missing, " individual(s) had one or more missing components for ",
      "phenotype '", phenotype_name, "' and were excluded. (IDs: ",
      paste(utils::head(missing_ids, 5L), collapse = ", "),
      if (n_missing > 5L) paste0(" ... +", n_missing - 5L, " more") else "", ")")
    if (missing_component_action == "error") stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
  }
  stats::setNames(composite, focal_ids)
}

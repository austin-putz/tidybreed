#' Define a trait with QTL and sampled effects in one call
#'
#' @description
#' Convenience wrapper that chains [define_trait()], [define_additive_effects()],
#' and [define_phenotype()] for a single uncorrelated trait using random QTL
#' placement. For correlated multi-trait simulations or custom QTL placement,
#' use the functions individually with
#' `get_table("genome_meta") |> filter(...) |> define_additive_effects()`.
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character. Trait name; must be a valid SQL identifier.
#' @param n_qtl Integer. Number of QTL to randomly select from `genome_meta`.
#' @param target_add_var Numeric. Additive genetic variance target. Passed to
#'   [define_trait()].
#' @param mean Numeric. Phenotypic population mean. Passed to
#'   [define_phenotype()] as the `mean` argument. Default `0`.
#' @param residual_var Numeric. Residual variance. Passed to [define_phenotype()].
#' @param type Character. One of `"continuous"` (default), `"count"`,
#'   `"categorical"`. Passed to [define_phenotype()].
#' @param expressed_sex Character. `"both"` (default), `"M"`, or `"F"`.
#'   Passed to [define_phenotype()].
#' @param repeatable Logical. Whether the trait allows repeated records.
#'   Passed to [define_phenotype()]. Default `FALSE`.
#' @param effect_distribution Character. `"normal"` (default) or `"gamma"`.
#'   Passed to [define_additive_effects()].
#' @param scale_to_target Logical. Passed to [define_additive_effects()].
#'   Default `TRUE`.
#' @param seed Optional integer for reproducibility (applied before QTL draw
#'   and effect sampling).
#' @param ... Additional arguments forwarded to [define_phenotype()] (e.g.
#'   `prevalence`, `thresholds`, `cat_values`, `cat_names`, `min_value`,
#'   `max_value`, `store_liability`).
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @seealso [define_trait()], [define_additive_effects()], [define_phenotype()]
#'
#' @examples
#' \dontrun{
#' pop <- pop |>
#'   define_trait_simple(
#'     trait_name     = "ADG",
#'     n_qtl          = 100,
#'     target_add_var = 100,
#'     mean           = 850,
#'     residual_var   = 120
#'   )
#' }
#' @export
define_trait_simple <- function(pop,
                                trait_name,
                                n_qtl,
                                target_add_var,
                                mean                = 0,
                                residual_var,
                                type                = "continuous",
                                expressed_sex       = "both",
                                repeatable          = FALSE,
                                effect_distribution = "normal",
                                scale_to_target     = TRUE,
                                seed                = NULL,
                                ...) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  stopifnot(is.numeric(n_qtl), length(n_qtl) == 1, n_qtl > 0)

  if (!is.null(seed)) set.seed(seed)

  pop <- define_trait(pop, trait_name = trait_name,
                      target_add_var = target_add_var)

  sel <- pop |>
    get_table("genome_meta") |>
    dplyr::collect() |>
    dplyr::slice_sample(n = as.integer(n_qtl)) |>
    dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_additive_effects(trait_name      = trait_name,
                            distribution    = effect_distribution,
                            scale_to_target = scale_to_target,
                            seed            = NULL)

  pop <- define_phenotype(pop,
                          phenotype_name = trait_name,
                          type           = type,
                          mean           = mean,
                          expressed_sex  = expressed_sex,
                          repeatable     = repeatable,
                          residual_var   = residual_var,
                          ...)
  invisible(pop)
}

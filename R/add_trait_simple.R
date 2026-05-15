#' Add a trait with QTL and sampled effects in one call
#'
#' @description
#' Convenience wrapper that chains [add_trait()], [define_qtl()], and
#' [set_qtl_effects()] for a single uncorrelated trait using random QTL
#' placement. For correlated multi-trait simulations, or to control QTL
#' placement (e.g. by chromosome or position), use the three functions
#' individually with `get_table("genome_meta") |> filter(...) |> define_qtl()`.
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Character. Trait name; must be a valid SQL identifier.
#' @param n_qtl Integer. Number of QTL to randomly select from `genome_meta`.
#' @param effect_distribution Character. `"normal"` or `"gamma"`.
#' @param scale_to_target Logical. Passed to [set_qtl_effects()].
#' @param seed Optional integer for reproducibility (applied before both the
#'   QTL random draw and the effect sampling).
#' @param ... Additional arguments forwarded to [add_trait()] (e.g.
#'   `trait_type`, `target_add_var`, `residual_var`, `target_add_mean`,
#'   `prevalence`, `expressed_sex`).
#'
#' @return The modified `tidybreed_pop` (invisibly).
#'
#' @examples
#' \dontrun{
#' pop <- pop |>
#'   add_trait_simple(
#'     trait_name          = "ADG",
#'     trait_type          = "continuous",
#'     n_qtl               = 100,
#'     target_add_var      = 0.25,
#'     residual_var        = 0.75,
#'     target_add_mean     = 850,
#'     effect_distribution = "normal"
#'   )
#' }
#' @export
add_trait_simple <- function(pop,
                             trait_name,
                             n_qtl,
                             effect_distribution = "normal",
                             scale_to_target     = TRUE,
                             seed                = NULL,
                             ...) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  stopifnot(is.numeric(n_qtl), length(n_qtl) == 1, n_qtl > 0)

  if (!is.null(seed)) set.seed(seed)

  pop <- add_trait(pop, trait_name = trait_name, ...)

  sel <- pop |>
    get_table("genome_meta") |>
    dplyr::collect() |>
    dplyr::slice_sample(n = as.integer(n_qtl)) |>
    dplyr::pull(locus_name)

  pop <- pop |>
    get_table("genome_meta") |>
    dplyr::filter(locus_name %in% sel) |>
    define_qtl(trait_name)

  pop <- set_qtl_effects(pop, trait_name = trait_name,
                         distribution    = effect_distribution,
                         scale_to_target = scale_to_target,
                         seed            = NULL)
  invisible(pop)
}

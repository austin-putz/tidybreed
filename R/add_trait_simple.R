#' Add a trait with QTL and sampled effects in one call
#'
#' @description
#' Convenience wrapper that chains [add_trait()], [define_qtl()], and
#' [set_qtl_effects()] for a single uncorrelated trait. For correlated
#' multi-trait simulations use the three low-level functions directly plus
#' [set_qtl_effects_multi()].
#'
#' @param pop A `tidybreed_pop` object.
#' @param trait_name Trait name.
#' @param n_qtl Integer. Number of QTL to select.
#' @param qtl_method Character. Selection method: `"random"`, `"even"`, or
#'   `"chromosome_even"`.
#' @param effect_distribution Character. `"normal"` or `"gamma"`.
#' @param scale_to_target Logical. Passed to [set_qtl_effects()].
#' @param seed Optional integer for reproducibility.
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
#'     qtl_method          = "chromosome_even",
#'     effect_distribution = "normal"
#'   )
#' }
#' @export
add_trait_simple <- function(pop,
                             trait_name,
                             n_qtl,
                             qtl_method          = "random",
                             effect_distribution = "normal",
                             scale_to_target     = TRUE,
                             seed                = NULL,
                             ...) {

  stopifnot(inherits(pop, "tidybreed_pop"))
  stopifnot(is.numeric(n_qtl), length(n_qtl) == 1, n_qtl > 0)

  pop <- add_trait(pop, trait_name = trait_name, ...)
  pop <- define_qtl(pop, trait_name = trait_name, n = as.integer(n_qtl),
                    method = qtl_method)
  pop <- set_qtl_effects(pop, trait_name = trait_name,
                         distribution    = effect_distribution,
                         scale_to_target = scale_to_target,
                         seed            = seed)
  invisible(pop)
}

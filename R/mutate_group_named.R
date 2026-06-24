#' Assign named group labels to filtered rows by count or proportion
#'
#' @description
#' Randomly assigns each filtered row to one of the supplied `group_names`.
#' Group membership is determined either by explicit `counts` (how many animals
#' per group) or by `proportions` (relative fractions). The column is always
#' `VARCHAR`. Useful for farm assignment, treatment groups, line splitting, and
#' any situation where groups have user-defined string labels.
#'
#' Pipe `get_table("ind_meta")` — optionally through `dplyr::filter()` — into
#' this function. Filtered primary keys are sorted ascending before shuffling,
#' giving a stable, reproducible starting point.
#'
#' @section Counts vs proportions:
#' Exactly one of `counts` or `proportions` must be supplied.
#'
#' - **`counts`**: integer vector, same length as `group_names`. When
#'   `sum(counts) > n_filtered`, `underfill_action` decides what happens.
#'   When `sum(counts) < n_filtered`, remaining rows are left `NULL` unless
#'   `include_leftover = TRUE` (which adds them to the last group).
#' - **`proportions`**: numeric vector, same length as `group_names`. Values
#'   must be positive and are normalized to sum to 1 within
#'   `sqrt(.Machine$double.eps)` tolerance. With `prop_method = "exact"` (default),
#'   the largest-remainder method converts proportions to integer counts,
#'   tie-breaking left-to-right. With `prop_method = "random"`, each animal is
#'   independently drawn from the distribution.
#'
#' @section Underfill behavior (counts mode only):
#' When `sum(counts) > n_filtered` and `underfill_action`:
#' - `"error"` (default): stop before sampling.
#' - `"proportional"`: rescale counts to fit `n_filtered` using the
#'   largest-remainder method.
#' - `"sequential"`: fill groups in order until animals run out; later groups
#'   may be partially or fully empty.
#'
#' @param tbl A `tidybreed_table` from `get_table()` (optionally filtered).
#' @param col_name Character scalar. Name of the group column to create or
#'   update. Must be a valid SQL identifier and not a reserved column.
#' @param group_names Character vector. Group labels. Must be the same length
#'   as `counts` or `proportions`.
#' @param counts Integer vector of animals per group. Mutually exclusive with
#'   `proportions`.
#' @param proportions Numeric vector of relative group fractions. Mutually
#'   exclusive with `counts`.
#' @param prop_method `"exact"` (default) uses the largest-remainder method for
#'   deterministic integer rounding; `"random"` draws each label independently.
#' @param include_leftover Logical. When `sum(counts) < n_filtered` (or when
#'   proportions sum to < n due to rounding), add the unassigned rows to the
#'   last group. Default `FALSE` leaves them `NULL`.
#' @param underfill_action `"error"` (default), `"proportional"`, or
#'   `"sequential"`. Only applies when `counts` is supplied and
#'   `sum(counts) > n_filtered`.
#' @param seed Integer or `NULL`. If supplied, sets the RNG seed before
#'   shuffling. Validation happens before the seed is set.
#' @param overwrite Logical. If `FALSE` (default), error when any filtered row
#'   already has a non-`NULL` value in `col_name`.
#' @param quiet Logical. If `TRUE`, suppress informational messages.
#'
#' @return The modified `tidybreed_pop` object (invisibly).
#'
#' @seealso [mutate_group_seq()], [mutate_group_concatenate()], [mutate_table()]
#'
#' @examples
#' \dontrun{
#' # Assign females to two farms by count
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(sex == "F") |>
#'   mutate_group_named(
#'     col_name    = "farm",
#'     group_names = c("farm_A", "farm_B"),
#'     counts      = c(60L, 40L)
#'   )
#'
#' # Assign all individuals by proportion (60 / 40 split)
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   mutate_group_named(
#'     col_name     = "farm",
#'     group_names  = c("farm_A", "farm_B"),
#'     proportions  = c(0.6, 0.4)
#'   )
#'
#' # If available animals are fewer than requested counts, rescale proportionally
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 1L) |>
#'   mutate_group_named(
#'     col_name         = "trt",
#'     group_names      = c("control", "high", "low"),
#'     counts           = c(50L, 50L, 50L),
#'     underfill_action = "proportional"
#'   )
#' }
#' @export
mutate_group_named <- function(tbl,
                                col_name,
                                group_names,
                                counts           = NULL,
                                proportions      = NULL,
                                prop_method      = "exact",
                                include_leftover = FALSE,
                                underfill_action = "error",
                                seed             = NULL,
                                overwrite        = FALSE,
                                quiet            = FALSE) {

  # --- Validate structure (before any RNG) ---
  info <- .resolve_group_target(tbl, col_name, "VARCHAR", overwrite)
  if (is.null(info)) return(invisible(tbl$pop))

  pop        <- info$pop
  conn       <- info$conn
  table_name <- info$table_name
  pk_col     <- info$pk_col
  pks        <- info$pks
  n_total    <- info$n_total

  if (!is.character(group_names) || length(group_names) < 1L || any(is.na(group_names))) {
    stop("'group_names' must be a non-empty character vector with no NA values.", call. = FALSE)
  }
  n_groups <- length(group_names)

  if (is.null(counts) == is.null(proportions)) {
    stop("Exactly one of 'counts' or 'proportions' must be specified.", call. = FALSE)
  }

  prop_method      <- match.arg(prop_method,      c("exact", "random"))
  underfill_action <- match.arg(underfill_action,  c("error", "proportional", "sequential"))

  # --- Validate counts ---
  if (!is.null(counts)) {
    if (!is.numeric(counts) || length(counts) != n_groups || any(is.na(counts)) || any(counts < 0L)) {
      stop(
        "'counts' must be a non-negative integer vector of length ", n_groups,
        " (same length as 'group_names').",
        call. = FALSE
      )
    }
    counts      <- as.integer(counts)
    n_required  <- sum(counts)

    if (n_total < n_required) {
      if (underfill_action == "error") {
        stop(
          "Available filtered rows (", n_total, ") are fewer than the total ",
          "requested (", n_required, "). ",
          "Set underfill_action = 'proportional' or 'sequential' to handle this.",
          call. = FALSE
        )
      }
      if (underfill_action == "proportional") {
        counts <- .largest_remainder(counts / n_required, n_total)
      }
      if (underfill_action == "sequential") {
        caps   <- pmax(0L, n_total - c(0L, cumsum(as.integer(counts[-n_groups]))))
        counts <- pmin(counts, caps)
      }
      n_required <- sum(counts)
    }

    if (n_total > n_required) {
      n_extra <- n_total - n_required
      if (include_leftover) {
        counts[n_groups] <- counts[n_groups] + n_extra
        n_required       <- n_total
      } else {
        warning(
          n_extra, " row", if (n_extra != 1L) "s" else "",
          " will be left unassigned (NULL). ",
          "Set include_leftover = TRUE to add them to the last group.",
          call. = FALSE
        )
      }
    }

    n_assigned <- sum(counts)
    label_vec  <- rep(group_names, times = counts)
  }

  # --- Validate proportions ---
  if (!is.null(proportions)) {
    if (!is.numeric(proportions) || length(proportions) != n_groups ||
        any(is.na(proportions)) || any(proportions <= 0)) {
      stop(
        "'proportions' must be a positive numeric vector of length ", n_groups,
        " (same length as 'group_names').",
        call. = FALSE
      )
    }
    prop_sum <- sum(proportions)
    if (abs(prop_sum - 1) > sqrt(.Machine$double.eps) * 100) {
      stop(
        "'proportions' must sum to approximately 1 (sum = ", round(prop_sum, 6), ").",
        call. = FALSE
      )
    }
    proportions <- proportions / prop_sum  # normalize
  }

  # --- Generate assignments ---
  if (!is.null(seed)) set.seed(seed)

  if (!is.null(counts)) {
    shuffled_pks <- sample(pks)
    values       <- rep(NA_character_, n_total)
    if (n_assigned > 0L) {
      idx_in_sorted          <- match(shuffled_pks[seq_len(n_assigned)], pks)
      values[idx_in_sorted]  <- label_vec
    }

  } else {
    if (prop_method == "exact") {
      counts_exact <- .largest_remainder(proportions, n_total)
      n_assigned   <- sum(counts_exact)
      label_vec    <- rep(group_names, times = counts_exact)
      shuffled_pks <- sample(pks)
      values       <- rep(NA_character_, n_total)
      if (n_assigned > 0L) {
        idx_in_sorted         <- match(shuffled_pks[seq_len(n_assigned)], pks)
        values[idx_in_sorted] <- label_vec
      }

    } else {
      # "random": each animal drawn independently
      values <- sample(group_names, n_total, replace = TRUE, prob = proportions)
    }
  }

  # --- Write and return ---
  .ensure_col_exists(conn, table_name, col_name, "VARCHAR")
  .write_group_values(conn, table_name, col_name, pk_col, pks, values)
  .summarize_group_assignment(values, n_total, table_name, col_name, quiet)

  invisible(pop)
}

#' Assign sequential integer group IDs to filtered rows
#'
#' @description
#' Randomly shuffles the filtered rows and assigns them to consecutive integer
#' group labels (pen 1, 2, 3 …). The column is always `INTEGER`. Useful for
#' pen assignment, contemporary group numbering, and any situation where each
#' group is identified by a unique integer.
#'
#' Pipe `get_table("ind_meta")` — optionally through `dplyr::filter()` — into
#' this function. Filtered primary keys are sorted ascending before shuffling,
#' giving a stable, reproducible starting point.
#'
#' @section Group-size modes:
#' Exactly one of `n_per_group` or `n_groups` must be supplied.
#'
#' - **`n_per_group` scalar**: every group has exactly `n_per_group` animals.
#'   Remainder animals are left `NULL` unless `include_leftover = TRUE`.
#' - **`n_per_group` length-2 vector `c(min, max)`**: each group is sized by
#'   drawing uniformly from `min:max`. When remaining animals fall below `min`,
#'   they are left `NULL` (or included if `include_leftover = TRUE`).
#' - **`n_groups`**: splits the filtered set into `n_groups` groups of size
#'   `floor(n / n_groups)`. When `include_leftover = TRUE`, the first
#'   `n %% n_groups` groups each receive one extra animal (balanced split).
#'   When `include_leftover = FALSE`, the remainder animals are left `NULL`.
#'
#' @section Start numbering:
#' When `start = NULL` (default), the next label is `MAX(col_name) + 1` if the
#' column already exists, or `1` for a new column. Supplying `start = 1L` resets
#' numbering — useful when assigning pens independently per generation.
#' Setting `start` does **not** require `overwrite = TRUE`; overwrite only
#' applies when **filtered rows** already hold non-`NULL` values.
#'
#' @param tbl A `tidybreed_table` from `get_table()` (optionally filtered).
#' @param col_name Character scalar. Name of the group column to create or
#'   update. Must be a valid SQL identifier and not a reserved column.
#' @param n_per_group Integer scalar or length-2 integer vector `c(min, max)`.
#'   Mutually exclusive with `n_groups`.
#' @param n_groups Positive integer scalar. Split filtered rows into this many
#'   groups. Mutually exclusive with `n_per_group`.
#' @param start Integer scalar. First group label. Default: `MAX(col_name) + 1`
#'   if the column exists, else `1`.
#' @param include_leftover Logical. If `TRUE`, remainder animals are assigned to
#'   a final (smaller) group. If `FALSE` (default), they are left `NULL`.
#' @param seed Integer or `NULL`. If supplied, sets the RNG seed before
#'   shuffling so results are reproducible. Validation happens before the seed
#'   is set, so a failed call does not consume RNG state.
#' @param overwrite Logical. If `FALSE` (default), error when any filtered row
#'   already has a non-`NULL` value in `col_name`.
#' @param quiet Logical. If `TRUE`, suppress informational messages.
#'
#' @return The modified `tidybreed_pop` object (invisibly).
#'
#' @seealso [mutate_group_named()], [mutate_group_concatenate()], [mutate_table()]
#'
#' @examples
#' \dontrun{
#' # Assign males to pens of 10 (remainder left NULL)
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(sex == "M") |>
#'   mutate_group_seq(col_name = "pen", n_per_group = 10L)
#'
#' # Assign all individuals to 5 balanced groups
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   mutate_group_seq(col_name = "pen", n_groups = 5L, include_leftover = TRUE)
#'
#' # Restart pen numbering from 1 for generation 2 (no overwrite = TRUE needed)
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   dplyr::filter(gen == 2L) |>
#'   mutate_group_seq(col_name = "pen", n_per_group = 10L, start = 1L)
#'
#' # Reproducible assignment
#' pop <- pop |>
#'   get_table("ind_meta") |>
#'   mutate_group_seq(col_name = "pen", n_per_group = 10L, seed = 42L)
#' }
#' @export
mutate_group_seq <- function(tbl,
                              col_name,
                              n_per_group      = NULL,
                              n_groups         = NULL,
                              start            = NULL,
                              include_leftover = FALSE,
                              seed             = NULL,
                              overwrite        = FALSE,
                              quiet            = FALSE) {

  # --- Validate structure (before any RNG) ---
  info <- .resolve_group_target(tbl, col_name, "INTEGER", overwrite)
  if (is.null(info)) return(invisible(tbl$pop))

  pop        <- info$pop
  conn       <- info$conn
  table_name <- info$table_name
  pk_col     <- info$pk_col
  pks        <- info$pks
  n_total    <- info$n_total

  if (is.null(n_per_group) == is.null(n_groups)) {
    stop("Exactly one of 'n_per_group' or 'n_groups' must be specified.", call. = FALSE)
  }

  if (!is.null(n_groups)) {
    if (!is.numeric(n_groups) || length(n_groups) != 1L ||
        is.na(n_groups) || n_groups < 1L) {
      stop("'n_groups' must be a positive integer scalar.", call. = FALSE)
    }
    n_groups <- as.integer(n_groups)
    if (n_groups > n_total) {
      stop(
        "'n_groups' (", n_groups, ") exceeds the number of filtered rows (",
        n_total, ").",
        call. = FALSE
      )
    }
  }

  if (!is.null(n_per_group)) {
    if (!is.numeric(n_per_group) || any(is.na(n_per_group)) || any(n_per_group < 1L)) {
      stop(
        "'n_per_group' must be a positive integer scalar or a length-2 vector [min, max].",
        call. = FALSE
      )
    }
    if (length(n_per_group) > 2L) {
      stop("'n_per_group' must be a scalar or a length-2 vector [min, max].", call. = FALSE)
    }
    if (length(n_per_group) == 2L && n_per_group[1L] > n_per_group[2L]) {
      stop("'n_per_group[1]' must be <= 'n_per_group[2]'.", call. = FALSE)
    }
    n_per_group <- as.integer(n_per_group)
  }

  # Determine start_num before sampling
  if (!is.null(start)) {
    start_num <- as.integer(start)
    if (is.na(start_num)) stop("'start' must be a finite integer.", call. = FALSE)
  } else {
    if (info$col_exists) {
      max_val   <- DBI::dbGetQuery(conn, paste0("SELECT MAX(", col_name, ") FROM ", table_name))[[1L]]
      start_num <- if (is.na(max_val)) 1L else as.integer(max_val) + 1L
    } else {
      start_num <- 1L
    }
  }

  # --- Generate assignments ---
  if (!is.null(seed)) set.seed(seed)

  if (!is.null(n_groups)) {
    group_size      <- as.integer(floor(n_total / n_groups))
    n_leftover_rows <- n_total %% n_groups

    if (include_leftover) {
      sizes <- rep(group_size, n_groups)
      if (n_leftover_rows > 0L) {
        sizes[seq_len(n_leftover_rows)] <- sizes[seq_len(n_leftover_rows)] + 1L
      }
    } else {
      sizes <- rep(group_size, n_groups)
    }

    n_assigned   <- sum(sizes)
    group_labels <- c(
      unlist(lapply(seq_along(sizes), function(g) rep(start_num + g - 1L, sizes[g]))),
      rep(NA_integer_, n_total - n_assigned)
    )

    shuffled_pks <- sample(pks)
    values       <- group_labels[match(pks, shuffled_pks)]

  } else {
    shuffled_pks <- sample(pks)

    if (length(n_per_group) == 1L) {
      group_size      <- n_per_group[1L]
      n_full          <- n_total %/% group_size
      n_leftover_rows <- n_total %% group_size

      labels <- rep(seq_len(n_full) + start_num - 1L, each = group_size)
      if (n_leftover_rows > 0L && include_leftover) {
        labels <- c(labels, rep(start_num + n_full, n_leftover_rows))
      }
      group_labels <- c(labels, rep(NA_integer_, n_total - length(labels)))

    } else {
      # Variable size: c(min, max)
      min_sz      <- n_per_group[1L]
      max_sz      <- n_per_group[2L]
      remaining   <- n_total
      g           <- 0L
      label_parts <- list()

      while (remaining > 0L) {
        size_sampled <- if (min_sz == max_sz) min_sz else sample(min_sz:max_sz, 1L)
        size         <- min(size_sampled, remaining)

        if (size < min_sz) {
          # Remaining animals can't fill a minimum-size group
          if (include_leftover) {
            label_parts <- c(label_parts, list(rep(start_num + g, size)))
          } else {
            label_parts <- c(label_parts, list(rep(NA_integer_, size)))
          }
          break
        }

        label_parts <- c(label_parts, list(rep(start_num + g, size)))
        g         <- g + 1L
        remaining <- remaining - size
      }

      group_labels <- unlist(label_parts)
      if (length(group_labels) < n_total) {
        group_labels <- c(group_labels, rep(NA_integer_, n_total - length(group_labels)))
      }
    }

    values <- group_labels[match(pks, shuffled_pks)]
  }

  # Warn about unassigned rows
  n_unassigned <- sum(is.na(values))
  if (!include_leftover && n_unassigned > 0L) {
    warning(
      n_unassigned, " row", if (n_unassigned != 1L) "s" else "",
      " left unassigned (NULL) in '", col_name, "'. ",
      "Set include_leftover = TRUE to include them.",
      call. = FALSE
    )
  }

  # --- Write and return ---
  .ensure_col_exists(conn, table_name, col_name, "INTEGER")
  .write_group_values(conn, table_name, col_name, pk_col, pks, values)
  .summarize_group_assignment(values, n_total, table_name, col_name, quiet)

  invisible(pop)
}

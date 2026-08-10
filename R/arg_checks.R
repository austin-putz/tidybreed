# Shared argument-validation helpers. Internal; none are exported.
#
# These are deliberately generic (no simulation semantics) so any exported
# function can validate its scalar arguments the same way and produce the same
# error wording.

# ---------------------------------------------------------------------------
# Scalar argument validation
# ---------------------------------------------------------------------------

# Validate one numeric scalar argument and return it. Every check is opt-in so
# each call site reads as a single line, and every failure names the argument
# and the value the user actually passed -- the previous `stopifnot()` calls
# reported expressions like "n_templates >= 2L is not TRUE" after coercion had
# already changed the value.
#
# `min`/`max` are inclusive unless `exclusive = TRUE`. `whole = TRUE` requires
# an exact whole number (no silent truncation) and returns an integer.
.check_scalar <- function(x, name, min = NULL, max = NULL, exclusive = FALSE,
                          whole = FALSE, finite = TRUE) {
  if (!is.numeric(x)) {
    stop("`", name, "` must be numeric, not ", class(x)[1], ".", call. = FALSE)
  }
  if (length(x) != 1L) {
    stop("`", name, "` must be a single value, got length ", length(x), ".",
         call. = FALSE)
  }
  if (is.na(x)) {
    stop("`", name, "` must not be NA.", call. = FALSE)
  }
  if (finite && !is.finite(x)) {
    stop("`", name, "` must be finite, got ", x, ".", call. = FALSE)
  }
  if (!is.null(min)) {
    ok <- if (exclusive) x > min else x >= min
    if (!ok) {
      stop("`", name, "` must be ", if (exclusive) "greater than " else "at least ",
           min, ", got ", x, ".", call. = FALSE)
    }
  }
  if (!is.null(max)) {
    ok <- if (exclusive) x < max else x <= max
    if (!ok) {
      stop("`", name, "` must be ", if (exclusive) "less than " else "at most ",
           max, ", got ", x, ".", call. = FALSE)
    }
  }
  if (whole) {
    if (x != trunc(x)) {
      stop("`", name, "` must be a whole number, got ", x, ".", call. = FALSE)
    }
    if (x > .Machine$integer.max) {
      stop("`", name, "` is too large to represent as an integer, got ", x, ".",
           call. = FALSE)
    }
    return(as.integer(x))
  }
  as.numeric(x)
}

# Validate a logical scalar (e.g. `exact_freq`, `recombines_M`).
.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

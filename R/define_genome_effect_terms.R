#' Reserved effect owners
#'
#' Owner names the package writes itself. The exported writer refuses to
#' create or delete rows under these, with no override, so that rerunning
#' [define_additive_effects()] can never remove a user's own terms, a user's
#' call can never remove the generated ones, and "every active term is
#' `generated`" proves the terms were calibrated to the stored target (Q21,
#' Q23). Only the package's internal `.ge_write_terms()` may write them. Derived from
#' `GE_GENERATED_OWNER` rather than repeating the literal, so the reserved list
#' and the writer that owns it cannot drift apart.
#'
#' @keywords internal
#' @noRd
GE_RESERVED_OWNERS <- GE_GENERATED_OWNER

#' Define genome effects as terms over one or more loci
#'
#' @description
#' Writes genome effects in the term / member / origin form: one coefficient
#' (`genome_effects`) over one or more loci (`genome_effect_members`), each
#' locus optionally scoped to allele copies of a given line and/or parent of
#' origin (`genome_effect_member_origins`).
#'
#' `terms` is a **long data frame, one row per (term x locus)** — the same shape
#' as every table [get_table()] returns. Scope is supplied separately in
#' `origin` so the common case stays flat.
#'
#' `define_genome_effect_terms()` writes any effect you supply; `define_*_effects()`
#' functions such as [define_additive_effects()] sample effects of one shape
#' and write them through the same path.
#'
#' @section The `terms` data frame:
#'
#' | Column | Required | Meaning |
#' |---|---|---|
#' | `term_id` | yes, unless one term | Groups rows into one term. **User-facing only** — never stored; the writer replaces it with an `id_genome_effect`. Any atomic type |
#' | `genome_value` | yes | The term's coefficient. Constant within a `term_id` |
#' | `effect_name` | no | Per-term label. Constant within a `term_id`; no mathematical meaning |
#' | `locus_name` | yes | Resolved to `locus_id`; each locus at most once per term |
#' | `contrast_name` | yes | `"additive"`, `"dominance"` or `"indicator"` |
#' | `center_value` | non-indicator | `p` for Cockerham coding, `0.5` for functional |
#' | `copy_count_value`, `dosage_value` | indicator | The local genotype state. `copy_count_value` is inferred at ordinary diploid-autosomal loci |
#'
#' A single-term call may omit `term_id` entirely.
#'
#' @section Scope (`origin`):
#'
#' * `NULL` (default) — the common scope: no origin rows, matches every allele
#'   copy.
#' * A **named scalar list**, e.g. `list(line_name = "Duroc", parent_origin = 1)`
#'   — one scope applied to every member of every term. Accepted names are
#'   `line_match_type`, `line_name`, `parent_origin` and `copy_count`.
#' * A **data frame** with columns `term_id`, `locus_name`, `line_match_type`,
#'   `line_name`, `parent_origin`, `copy_count` — per-member scopes, needed for
#'   the exact multisets a `dominance` or `indicator` member takes. Keyed by
#'   `locus_name`, so you never touch canonical slot order.
#'
#' An `additive` member takes at most one origin row; a genotype member takes an
#' exact multiset whose `copy_count`s sum to the state's copy count (2 for
#' `dominance`, `copy_count_value` for `indicator`). A scalar list therefore
#' scopes `additive` members but is usually not enough for a genotype member —
#' the validator says so, naming the locus.
#'
#' @section Replacement modes:
#'
#' * `"append"` — insert; a term duplicating an existing family + scope identity
#'   is rejected.
#' * `"replace_scope"` — delete only the variants in `(trait_name, effect_owner)`
#'   whose origin predicate **equals** the supplied scope, then insert. This is
#'   what [define_additive_effects()] uses, so successive common / line-A /
#'   line-B calls each replace only their own variant. Requires `origin` to be
#'   `NULL` or a scalar list (a per-member `origin` data frame has no single
#'   scope to key on).
#' * `"replace_owner"` — replace everything under `(trait_name, effect_owner)`.
#' * `"replace_trait"` — clear every owner's terms for the trait. Never implied.
#'
#' @param pop A `tidybreed_pop`.
#' @param trait_name Character scalar. Must exist in `trait_meta`.
#' @param terms Long data frame; see **The `terms` data frame**.
#' @param effect_owner Character scalar naming the writer that owns these rows,
#'   **for replacement only**. Owners always sum and are never selected between,
#'   because `effect_owner` is part of the fallback-family signature. Default
#'   `"custom"`.
#' @param mode One of `"append"`, `"replace_scope"`, `"replace_owner"`,
#'   `"replace_trait"`.
#' @param origin Scope; see **Scope (`origin`)**.
#' @param base_tbl Optional `tidybreed_table` from [get_table()] (optionally
#'   filtered) selecting the allele copies whose frequencies fill
#'   `center_value` on any `additive` or `dominance` member that has none:
#'   `founder_haplotypes`, `ind_haplotype`, or any table with an `id_ind`
#'   column, from the same `pop` — see [extract_allele_freq()]. An explicit
#'   `center_value` is never overwritten; `indicator` members are never
#'   touched; the base is queried only if some centre is actually missing.
#'   `NULL` (default) fills nothing, and a missing centre is an error. The
#'   fill is Cockerham `p` only — functional coding writes `0.5` explicitly
#'   ([ad_terms()] does). One `base_tbl` gives one `p` per locus per call, so
#'   line-scoped surfaces are written one line at a time with
#'   `mode = "append"`.
#' @param require_complete Logical. When `TRUE`, an indicator surface must name
#'   every reachable `(copy_count, dosage)` state on every member — including
#'   `copy_count_value = 0` where a chromosome can be absent. Default `FALSE`
#'   (sparse: a cell you do not write contributes zero).
#'
#' @return The `tidybreed_pop`, invisibly.
#'
#' @seealso [ad_terms()] and [genotype_terms()] build `terms`;
#'   [define_additive_effects()] for generated additive QTL effects.
#'
#' @examples
#' \dontrun{
#' # One dominance term, Cockerham coding at p = 0.3
#' pop <- pop |> define_genome_effect_terms(
#'   trait_name = "ADG",
#'   terms = data.frame(locus_name    = "Locus_10",
#'                      contrast_name = "dominance",
#'                      center_value  = 0.3,
#'                      genome_value  = 0.8)
#' )
#'
#' # A hand-entered 3x3 A x A surface: nine cells, nine terms, two members each
#' cells <- expand.grid(g1 = 0:2, g2 = 0:2)
#' cells$value <- c(0, 0, 0,  0, 1.4, 2.1,  0, 2.1, 3.6)
#' surface <- rbind(
#'   data.frame(term_id = seq_len(nrow(cells)), locus_name = "Locus_10",
#'              contrast_name = "indicator", dosage_value = cells$g1,
#'              genome_value  = cells$value),
#'   data.frame(term_id = seq_len(nrow(cells)), locus_name = "Locus_44",
#'              contrast_name = "indicator", dosage_value = cells$g2,
#'              genome_value  = cells$value)
#' )
#' pop <- pop |> define_genome_effect_terms("ADG", surface[surface$genome_value != 0, ],
#'                                          effect_owner = "epistasis_AxA")
#'
#' # Reciprocal dominance: the F1 value depends on which parent gave which line
#' pop <- pop |> define_genome_effect_terms(
#'   "ADG",
#'   terms = data.frame(term_id = 1L, locus_name = "Locus_10",
#'                      contrast_name = "dominance",
#'                      center_value = 0.3, genome_value = 1.2),
#'   origin = data.frame(term_id = 1L, locus_name = "Locus_10",
#'                       line_match_type = "exact",
#'                       line_name     = c("Duroc", "Landrace"),
#'                       parent_origin = c(1L, 2L),
#'                       copy_count    = c(1L, 1L)),
#'   effect_owner = "reciprocal"
#' )
#'
#' # Let the writer fill Cockerham p from a base population: leave
#' # center_value out and pass base_tbl (see extract_allele_freq()).
#' pop <- pop |> define_genome_effect_terms(
#'   "ADG",
#'   data.frame(locus_name = "Locus_10", contrast_name = "dominance",
#'              genome_value = 0.8),
#'   base_tbl = get_table(pop, "founder_haplotypes") |>
#'     dplyr::filter(line_name == "A")
#' )
#' }
#' @export
define_genome_effect_terms <- function(pop,
                                       trait_name,
                                       terms,
                                       effect_owner        = "custom",
                                       mode                = c("append",
                                                               "replace_scope",
                                                               "replace_owner",
                                                               "replace_trait"),
                                       origin              = NULL,
                                       base_tbl            = NULL,
                                       require_complete    = FALSE) {
  mode <- match.arg(mode)
  .ge_write_terms(pop, trait_name, terms, effect_owner, mode, origin,
                  base_tbl, require_complete, allow_reserved_owner = FALSE)
}

#' The writer behind [define_genome_effect_terms()]
#'
#' `allow_reserved_owner = TRUE` is internal only: the exported function always
#' passes `FALSE`, so no exported path can write coefficients under the
#' reserved owner (Q23). Tests use it to show the generator is sugar over this
#' engine.
#' @keywords internal
#' @noRd
.ge_write_terms <- function(pop, trait_name, terms, effect_owner = "custom",
                            mode = "append", origin = NULL, base_tbl = NULL,
                            require_complete = FALSE,
                            allow_reserved_owner = FALSE) {
  validate_tidybreed_pop(pop)
  conn <- pop$db_conn

  # base_tbl is validated whenever supplied (a wrong object should error
  # whether or not anything needs filling) but queried only when an additive
  # or dominance row is missing its centre. An omitted center_value column is
  # the documented way to ask for a fill, so normalise it the way .ge_build()
  # will before deciding -- terms$center_value is NULL then, and
  # any(NULL & ...) is FALSE.
  base_freq <- NULL
  if (!is.null(base_tbl)) {
    .validate_base_tbl(base_tbl, pop)
    if (is.data.frame(terms) && nrow(terms) > 0L &&
        "contrast_name" %in% names(terms)) {
      centres <- if ("center_value" %in% names(terms)) terms$center_value
                 else rep(NA_real_, nrow(terms))
      needs_fill <- any(is.na(centres) &
                        terms$contrast_name %in% c("additive", "dominance"))
      if (isTRUE(needs_fill)) base_freq <- extract_allele_freq(base_tbl)
    }
  }

  validate_sql_identifier(trait_name, what = "trait name")
  validate_sql_identifier(effect_owner, what = "effect owner")
  if (!isTRUE(allow_reserved_owner) && effect_owner %in% GE_RESERVED_OWNERS) {
    stop("'", effect_owner, "' is a reserved effect owner written by the ",
         "package itself (", paste(GE_RESERVED_OWNERS, collapse = ", "), "). ",
         "Choose your own effect_owner — owners always sum, so a separate one ",
         "costs nothing and keeps define_additive_effects() from deleting your ",
         "rows.", call. = FALSE)
  }
  .ge_require_effect_tables(pop)
  .ge_require_trait(conn, trait_name)

  built  <- .ge_build(conn, trait_name, terms, origin, effect_owner,
                      base_freq = base_freq)
  scope  <- .ge_scope_from_origin(origin, mode)
  model  <- .ge_read_model(conn)
  drop   <- .ge_resolve_deletes(model, trait_name, effect_owner, mode, scope,
                                allow_reserved_owner)

  if (isTRUE(require_complete)) .ge_assert_complete(conn, built)

  .ge_commit(conn, drop, built)
  message("Wrote ", nrow(built$terms), " genome-effect term",
          if (nrow(built$terms) == 1L) "" else "s", " (",
          nrow(built$members), " member rows) for trait '", trait_name,
          "' under owner '", effect_owner, "' [", mode, "].")
  invisible(pop)
}


# ── Input parsing ───────────────────────────────────────────────────────────

GE_TERMS_COLS  <- c("term_id", "genome_value", "effect_name", "locus_name",
                    "contrast_name", "center_value", "copy_count_value",
                    "dosage_value")
GE_ORIGIN_COLS <- c("term_id", "locus_name", "line_match_type", "line_name",
                    "parent_origin", "copy_count")

#' Turn the user's `terms` (+ `origin`) into canonical candidate rows
#'
#' Returns `list(terms, members, origins, labels)` where the three frames match
#' the three tables except that `id_genome_effect` is a local 1..n index, and
#' `labels` maps that index to the `term_id` the user typed. Every error raised
#' here names the user's `term_id`.
#'
#' @keywords internal
#' @noRd
.ge_build <- function(conn, trait_name, terms, origin, effect_owner,
                      base_freq = NULL) {
  if (!is.data.frame(terms) || nrow(terms) == 0L) {
    stop("'terms' must be a data frame with at least one row.", call. = FALSE)
  }
  unknown <- setdiff(names(terms), GE_TERMS_COLS)
  if (length(unknown) > 0L) {
    stop("Unknown column", if (length(unknown) == 1L) "" else "s",
         " in 'terms': ", paste0("'", unknown, "'", collapse = ", "),
         ". Accepted: ", paste(GE_TERMS_COLS, collapse = ", "), ".",
         call. = FALSE)
  }
  for (need in c("locus_name", "contrast_name", "genome_value")) {
    if (!need %in% names(terms)) {
      stop("'terms' must have a '", need, "' column.", call. = FALSE)
    }
  }
  if (!"term_id" %in% names(terms)) {
    if (length(unique(terms$locus_name)) != nrow(terms)) {
      stop("'terms' has no 'term_id' column, so it is read as a single term, ",
           "but it names the same locus more than once. Add 'term_id'.",
           call. = FALSE)
    }
    terms$term_id <- 1L
  }
  if (anyNA(terms$term_id)) stop("'term_id' must not be NA.", call. = FALSE)
  if (!is.atomic(terms$term_id)) {
    stop("'term_id' must be an atomic vector.", call. = FALSE)
  }
  keys   <- as.character(terms$term_id)
  labels <- unique(keys)                      # user order, not sorted
  idx    <- match(keys, labels)

  loci <- DBI::dbGetQuery(conn, "SELECT locus_id, locus_name FROM genome_meta")
  miss <- setdiff(unique(terms$locus_name), loci$locus_name)
  if (length(miss) > 0L) {
    stop("Locus name", if (length(miss) == 1L) "" else "s",
         " not in genome_meta: ", paste0("'", head(miss, 5), "'", collapse = ", "),
         if (length(miss) > 5L) paste0(" (", length(miss), " total)") else "",
         ".", call. = FALSE)
  }
  terms$locus_id <- loci$locus_id[match(terms$locus_name, loci$locus_name)]

  bad_contrast <- setdiff(terms$contrast_name,
                          c("additive", "dominance", "indicator"))
  if (length(bad_contrast) > 0L) {
    stop("Unknown contrast_name: ",
         paste0("'", unique(bad_contrast), "'", collapse = ", "),
         ". Use 'additive', 'dominance' or 'indicator'.", call. = FALSE)
  }

  for (col in c("center_value", "copy_count_value", "dosage_value",
                "effect_name")) {
    if (!col %in% names(terms)) {
      terms[[col]] <- if (col == "effect_name") NA_character_ else NA_real_
    }
  }
  # Coercion happens below via as.numeric()/as.integer(), which turn a
  # character column into silent NAs and a fractional state into a truncated
  # one. Neither should reach storage as a plausible-looking value.
  for (col in c("genome_value", "center_value")) {
    if (!is.numeric(terms[[col]]) && !all(is.na(terms[[col]]))) {
      stop("'", col, "' must be numeric (got ", class(terms[[col]])[1], ").",
           call. = FALSE)
    }
  }
  for (col in c("copy_count_value", "dosage_value")) {
    x <- terms[[col]]
    if (!is.numeric(x) && !all(is.na(x))) {
      stop("'", col, "' must be a whole number (got ", class(x)[1], ").",
           call. = FALSE)
    }
    if (is.numeric(x) && any(!is.na(x) & (x != trunc(x) | x < 0))) {
      stop("'", col, "' must be a non-negative whole number.", call. = FALSE)
    }
  }

  # Per-term scalars must actually be scalar within the term. Each rule is
  # decided for every term at once (a per-term subset is quadratic), and the
  # error is the first violation by (term in user order, rule): the term
  # comes first, so a later rule failing on an earlier term wins.
  n_terms <- length(labels)
  first   <- match(seq_len(n_terms), idx)              # first row of each term
  val1    <- terms$genome_value[first]
  nm1     <- terms$effect_name[first]
  gv      <- terms$genome_value
  en      <- terms$effect_name
  bad_value <- is.na(val1) |
    .ge_any_by(idx, is.na(gv) | gv != val1[idx], n_terms)
  bad_name  <- .ge_any_by(idx, ifelse(is.na(en) | is.na(nm1[idx]),
                                      is.na(en) != is.na(nm1[idx]),
                                      en != nm1[idx]), n_terms)
  dup_row   <- duplicated(.ge_pair_key(idx, terms$locus_id))
  bad_locus <- .ge_any_by(idx, dup_row, n_terms)
  bad <- which(bad_value | bad_name | bad_locus)
  if (length(bad) > 0L) {
    k   <- bad[1]
    tag <- paste0("term_id '", labels[k], "'")
    sub <- terms[idx == k, , drop = FALSE]
    if (bad_value[k]) {
      stop(tag, ": 'genome_value' must be one non-missing value for the whole ",
           "term (got ", paste(unique(sub$genome_value), collapse = ", "),
           "). A term is one coefficient over its members.", call. = FALSE)
    }
    if (bad_name[k]) {
      stop(tag, ": 'effect_name' must be constant within a term (got ",
           paste0("'", unique(sub$effect_name), "'", collapse = ", "), ").",
           call. = FALSE)
    }
    dup <- unique(sub$locus_name[duplicated(sub$locus_id)])
    stop(tag, ": locus ", paste0("'", dup, "'", collapse = ", "),
         " appears more than once in one term. Each locus contributes at ",
         "most one member.", call. = FALSE)
  }
  term_frame <- data.frame(
    id_genome_effect = seq_len(n_terms),
    trait_name       = rep(trait_name, n_terms),
    effect_owner     = rep(effect_owner, n_terms),
    effect_name      = as.character(nm1),
    genome_value     = as.numeric(val1),
    stringsAsFactors = FALSE
  )

  # Members, canonicalized by ascending locus_id within each term.
  o   <- order(idx, terms$locus_id)
  sub <- terms[o, , drop = FALSE]
  members <- data.frame(
    id_genome_effect = idx[o],
    member_slot      = sequence(tabulate(idx, n_terms)),
    locus_id         = as.integer(sub$locus_id),
    locus_name       = as.character(sub$locus_name),
    contrast_name    = as.character(sub$contrast_name),
    copy_count_value = as.integer(sub$copy_count_value),
    dosage_value     = as.integer(sub$dosage_value),
    center_value     = as.numeric(sub$center_value),
    stringsAsFactors = FALSE
  )

  # Fill Cockerham centres from the base where the user left them out, on
  # additive/dominance members only. An explicit centre always wins; an
  # indicator has no centre and is never touched. fill_failed exists only so
  # the validator can say why a centre is still missing.
  members$fill_failed <- FALSE
  if (!is.null(base_freq)) {
    fill <- is.na(members$center_value) &
            members$contrast_name %in% c("additive", "dominance")
    members$center_value[fill] <-
      base_freq$allele_freq[match(members$locus_id[fill], base_freq$locus_id)]
    members$fill_failed <- fill & is.na(members$center_value)
  }

  members <- .ge_infer_copy_counts(conn, members, labels)
  .ge_check_member_fields(members, labels)
  members$fill_failed <- NULL

  origins <- .ge_build_origins(origin, members, labels)
  .ge_check_origin_fields(origins, members, labels)
  list(terms   = term_frame,
       members = members,
       origins = origins,
       labels  = stats::setNames(labels, seq_along(labels)))
}

#' Per group: does any row of the group satisfy `flag`?
#'
#' `idx` holds group numbers `1..n`. An `NA` flag counts as `TRUE`, so a rule
#' that cannot be decided is reported rather than passed.
#'
#' @keywords internal
#' @noRd
.ge_any_by <- function(idx, flag, n) {
  flag[is.na(flag)] <- TRUE
  tabulate(idx[flag], n) > 0L
}

#' One exact numeric key per (group, integer value) pair
#'
#' For `duplicated()` / `match()` over two integer columns without pasting
#' strings. The key `x * radix + y`, with `radix = max(y) + 1`, is below
#' `(max(x) + 1) * radix`; it is used only while that bound is at most 2^53, where
#' every key is an exact double and distinct pairs get distinct keys. Rounding
#' is monotonic, so the bound computed in doubles never passes a key that is
#' really larger. Otherwise -- a bound too large, or a missing, negative or
#' non-whole value -- the key is pasted strings, which compare `NA` as a value
#' the way `duplicated()` does on one column.
#'
#' @keywords internal
#' @noRd
.ge_pair_key <- function(x, y) {
  whole <- function(v) !anyNA(v) && all(v >= 0) && all(v == trunc(v))
  if (!whole(x) || !whole(y)) return(paste(x, y))
  radix <- max(c(0, y)) + 1
  if ((max(c(0, x)) + 1) * radix > 2^53) return(paste(x, y))
  as.numeric(x) * radix + as.numeric(y)
}

#' Infer `copy_count_value` for indicator members at diploid-autosomal loci
#'
#' Dosage alone is not a genotype state, so the column is required — but at an
#' ordinary autosome there is exactly one copy count it can be, and making every
#' user type `copy_count_value = 2` for a 3x3 surface buys nothing.
#'
#' @keywords internal
#' @noRd
.ge_infer_copy_counts <- function(conn, members, labels) {
  need <- members$contrast_name == "indicator" & is.na(members$copy_count_value)
  if (!any(need)) return(members)

  counts <- .ge_reachable_copy_counts(conn)
  chr_of <- DBI::dbGetQuery(conn, "SELECT locus_id, chr_name FROM genome_meta")
  rows <- which(need)
  chr  <- chr_of$chr_name[match(members$locus_id[rows], chr_of$locus_id)]
  # One copy count per chromosome, or the first member (in member order) on a
  # chromosome where dosage alone is not a state.
  single <- vapply(counts, length, integer(1)) == 1L
  bad <- which(!single[chr] %in% TRUE)
  if (length(bad) > 0L) {
    i  <- rows[bad[1]]
    cc <- counts[[chr[bad[1]]]]
    stop("term_id '", labels[members$id_genome_effect[i]], "' locus '",
         members$locus_name[i], "': copy_count_value cannot be inferred — ",
         "chromosome '", chr[bad[1]], "' can carry ",
         paste(cc, collapse = " or "),
         " copies, so dosage alone is not a genotype state here. Supply ",
         "copy_count_value explicitly.", call. = FALSE)
  }
  members$copy_count_value[rows] <- unlist(counts[chr], use.names = FALSE)
  members
}

#' Reachable total copy counts per chromosome, over both offspring sexes
#'
#' Line-agnostic: an effect definition is population-level.
#'
#' @keywords internal
#' @noRd
.ge_reachable_copy_counts <- function(conn) {
  inh <- lapply(c("M", "F"), function(sx) resolve_chr_inheritance(conn, sx))
  chrs <- inh[[1]]$chr_name
  stats::setNames(lapply(seq_along(chrs), function(i) {
    sort(unique(vapply(inh, function(d) {
      as.integer(d$from_parent_1[i] + d$from_parent_2[i])
    }, integer(1))))
  }), chrs)
}

#' Row-local `contrast_name` <-> state/center agreement, before any SQL
#'
#' Duplicates the SQL CHECK deliberately: the constraint fires with a DuckDB
#' message naming a surrogate id the user has never seen, this one names the
#' `term_id` they typed.
#'
#' @keywords internal
#' @noRd
.ge_check_member_fields <- function(members, labels) {
  # Each member yields at most two messages, in this order: the state/centre
  # rule, then the "must not carry" rule. Built for every row at once and
  # interleaved row by row, so the list reads exactly as a per-row loop would.
  tag <- paste0("term_id '", labels[members$id_genome_effect], "' locus '",
                members$locus_name, "'")
  kind <- members$contrast_name
  cc   <- members$copy_count_value
  dv   <- members$dosage_value
  ctr  <- members$center_value
  ind  <- kind == "indicator"
  first <- second <- rep(NA_character_, nrow(members))

  w <- ind & is.na(dv)
  first[w] <- paste(tag[w], "is an indicator and needs 'dosage_value'")
  w <- ind & !is.na(dv) & !is.na(cc) & dv > cc
  first[w] <- paste0(tag[w], ": dosage_value (", dv[w],
                     ") exceeds copy_count_value (", cc[w], ")")
  w <- ind & !is.na(ctr)
  second[w] <- paste(tag[w], "is an indicator and must not carry",
                     "'center_value' (an indicator is a state, not a",
                     "centred contrast)")

  fail <- if (is.null(members$fill_failed)) logical(nrow(members))
          else members$fill_failed %in% TRUE
  w <- !ind & is.na(ctr)
  first[w] <- paste0(tag[w], " is '", kind[w],
                     "' and needs 'center_value' (p under Cockerham",
                     " coding, 0.5 under functional)",
                     ifelse(fail[w],
                            " -- and base_tbl has no allele copies at this locus",
                            ""))
  w <- !ind & !is.na(ctr) & (ctr < 0 | ctr > 1)
  first[w] <- paste0(tag[w], ": center_value must be between 0 and 1 (got ",
                     ctr[w], ")")
  w <- !ind & (!is.na(cc) | !is.na(dv))
  second[w] <- paste0(tag[w], " is '", kind[w], "' and must ",
                      "not carry an indicator state")

  v <- c(rbind(first, second))
  v <- v[!is.na(v)]
  if (length(v) > 0L) {
    stop("Invalid 'terms':\n  - ", paste(unique(v), collapse = "\n  - "),
         call. = FALSE)
  }
  invisible(NULL)
}

#' Row-local origin rules, before any SQL
#'
#' These are declared `CHECK` constraints on the table, but a constraint that
#' fires mid-transaction reports a DuckDB expression and aborts the transaction
#' in a state the caller then has to unwind. Checking here means the user gets a
#' sentence about their own call and the write never starts.
#'
#' @keywords internal
#' @noRd
.ge_check_origin_fields <- function(origins, members, labels) {
  if (nrow(origins) == 0L) return(invisible(NULL))
  locus <- members$locus_name[match(paste(origins$id_genome_effect,
                                          origins$member_slot),
                                    paste(members$id_genome_effect,
                                          members$member_slot))]
  tag <- paste0("term_id '", labels[origins$id_genome_effect], "' locus '",
                locus, "'")
  lmt <- origins$line_match_type
  ln  <- origins$line_name
  po  <- origins$parent_origin
  # One column per rule, in rule order; read row by row. An unknown match type
  # reports only itself, as nothing else about the row is meaningful.
  n <- nrow(origins)
  msg <- matrix(NA_character_, nrow = 5L, ncol = n)
  known <- lmt %in% c("exact", "unknown", "any")
  w <- !known
  msg[1, w] <- paste0(tag[w], ": line_match_type must be 'exact', 'unknown' ",
                      "or 'any' (got '", lmt[w], "')")
  w <- known & lmt == "exact" & is.na(ln)
  msg[1, w] <- paste(tag[w], "is line_match_type 'exact' and needs a line_name")
  w <- known & lmt != "exact" & !is.na(ln)
  msg[2, w] <- paste0(tag[w], ": line_match_type '", lmt[w], "' takes no ",
                      "line_name (got '", ln[w], "')")
  w <- known & lmt == "any" & is.na(po)
  msg[3, w] <- paste(tag[w], "is line_match_type 'any' with no parent_origin,",
                     "which constrains nothing and is not a second spelling",
                     "of the common scope. Use origin = NULL for that, or",
                     "give a parent_origin")
  w <- known & !is.na(po) & !po %in% c(1L, 2L)
  msg[4, w] <- paste0(tag[w], ": parent_origin must be 1 (sire) or 2 (dam), ",
                      "or NA for either (got ", po[w], ")")
  w <- known & (is.na(origins$copy_count) | origins$copy_count < 1L)
  msg[5, w] <- paste(tag[w], "needs copy_count >= 1")
  v <- c(msg)
  v <- v[!is.na(v)]
  if (length(v) > 0L) {
    stop("Invalid 'origin':\n  - ", paste(unique(v), collapse = "\n  - "),
         call. = FALSE)
  }
  invisible(NULL)
}

#' Expand `origin` into canonical origin rows
#'
#' @keywords internal
#' @noRd
.ge_build_origins <- function(origin, members, labels) {
  empty <- data.frame(id_genome_effect = integer(0), member_slot = integer(0),
                      origin_slot = integer(0), line_match_type = character(0),
                      line_name = character(0), parent_origin = integer(0),
                      copy_count = integer(0), stringsAsFactors = FALSE)
  if (is.null(origin)) return(empty)

  if (is.list(origin) && !is.data.frame(origin)) {
    one <- .ge_origin_from_list(origin)
    raw <- data.frame(
      id_genome_effect = members$id_genome_effect,
      member_slot      = members$member_slot,
      line_match_type  = one$line_match_type,
      line_name        = one$line_name,
      parent_origin    = one$parent_origin,
      copy_count       = one$copy_count,
      stringsAsFactors = FALSE
    )
    return(.ge_canonical_origins(raw))
  }

  if (!is.data.frame(origin)) {
    stop("'origin' must be NULL, a named scalar list, or a data frame.",
         call. = FALSE)
  }
  unknown <- setdiff(names(origin), GE_ORIGIN_COLS)
  if (length(unknown) > 0L) {
    stop("Unknown column", if (length(unknown) == 1L) "" else "s",
         " in 'origin': ", paste0("'", unknown, "'", collapse = ", "),
         ". Accepted: ", paste(GE_ORIGIN_COLS, collapse = ", "), ".",
         call. = FALSE)
  }
  if (!"locus_name" %in% names(origin)) {
    stop("'origin' data frames are keyed by 'locus_name' (never member_slot), ",
         "so a 'locus_name' column is required.", call. = FALSE)
  }
  if (!"term_id" %in% names(origin)) {
    if (length(labels) != 1L) {
      stop("'origin' needs a 'term_id' column when 'terms' defines more than ",
           "one term.", call. = FALSE)
    }
    origin$term_id <- labels
  }
  for (col in c("line_match_type", "line_name", "parent_origin", "copy_count")) {
    if (!col %in% names(origin)) {
      origin[[col]] <- switch(col,
                              line_match_type = NA_character_,
                              line_name       = NA_character_,
                              parent_origin   = NA_integer_,
                              copy_count      = 1L)
    }
  }

  k <- match(as.character(origin$term_id), labels)
  if (anyNA(k)) {
    stop("'origin' names term_id ",
         paste0("'", unique(as.character(origin$term_id)[is.na(k)]), "'",
                collapse = ", "),
         ", which does not appear in 'terms'.", call. = FALSE)
  }
  slot <- members$member_slot[match(paste(k, origin$locus_name),
                                    paste(members$id_genome_effect,
                                          members$locus_name))]
  if (anyNA(slot)) {
    i <- which(is.na(slot))[1]
    stop("'origin' scopes locus '", origin$locus_name[i], "' on term_id '",
         as.character(origin$term_id)[i],
         "', which has no member at that locus.", call. = FALSE)
  }

  lmt <- as.character(origin$line_match_type)
  ln  <- as.character(origin$line_name)
  po  <- suppressWarnings(as.integer(origin$parent_origin))
  lmt[is.na(lmt)] <- ifelse(!is.na(ln[is.na(lmt)]), "exact",
                            ifelse(!is.na(po[is.na(lmt)]), "any", NA_character_))
  if (anyNA(lmt)) {
    stop("'origin' row ", which(is.na(lmt))[1], " constrains nothing: give a ",
         "'line_name', a 'parent_origin', or an explicit 'line_match_type'.",
         call. = FALSE)
  }
  for (nm in stats::na.omit(unique(ln))) {
    validate_sql_identifier(nm, what = "line name")
  }
  .ge_canonical_origins(data.frame(
    id_genome_effect = k,
    member_slot      = as.integer(slot),
    line_match_type  = lmt,
    line_name        = ln,
    parent_origin    = po,
    copy_count       = as.integer(origin$copy_count),
    stringsAsFactors = FALSE
  ))
}

#' Normalize a scalar-list scope into one origin row's fields
#'
#' @keywords internal
#' @noRd
.ge_origin_from_list <- function(origin) {
  unknown <- setdiff(names(origin),
                     c("line_match_type", "line_name", "parent_origin",
                       "copy_count"))
  if (length(unknown) > 0L || is.null(names(origin)) || any(names(origin) == "")) {
    stop("A scalar-list 'origin' takes only named entries line_match_type, ",
         "line_name, parent_origin, copy_count",
         if (length(unknown) > 0L)
           paste0(" (got ", paste0("'", unknown, "'", collapse = ", "), ")")
         else "", ".", call. = FALSE)
  }
  if (any(lengths(origin) != 1L)) {
    stop("Every entry of a scalar-list 'origin' must have length 1; ",
         paste0("'", names(origin)[lengths(origin) != 1L], "'",
                collapse = ", "), " is longer. A scalar list is one scope ",
         "applied to every member. For per-member scopes use the data-frame ",
         "form; for \"line A or line B\" write one call per line, which is ",
         "what makes them separate variants that can be replaced ",
         "independently.", call. = FALSE)
  }
  ln  <- if (is.null(origin$line_name)) NA_character_ else
    as.character(origin$line_name)
  if (!is.na(ln)) validate_sql_identifier(ln, what = "line name")
  po  <- if (is.null(origin$parent_origin)) NA_integer_ else
    as.integer(origin$parent_origin)
  lmt <- if (!is.null(origin$line_match_type)) {
    as.character(origin$line_match_type)
  } else if (!is.na(ln)) "exact" else if (!is.na(po)) "any" else NA_character_
  if (is.na(lmt)) {
    stop("'origin' constrains nothing: give a 'line_name', a 'parent_origin', ",
         "or an explicit 'line_match_type'. Use origin = NULL for the common ",
         "scope.", call. = FALSE)
  }
  list(line_match_type = lmt, line_name = ln, parent_origin = po,
       copy_count = if (is.null(origin$copy_count)) 1L
                    else as.integer(origin$copy_count))
}

#' Canonical origin ordering: `origin_slot` runs 1..n over the sorted tuple
#'
#' Sorted by `line_match_type`, then `line_name`, then `parent_origin`, then
#' `copy_count`, **NULLs last** on each key, so a scope written in either order
#' stores identically and `replace_scope` can compare stored rows to fresh
#' input. A NULL `line_name` is in fact unreachable as a tie-breaker (it occurs
#' only under `unknown`/`any`, which differ already on the primary key), but the
#' sentinel is written anyway so the code says what the rule is rather than what
#' the current lattice happens to allow.
#'
#' @keywords internal
#' @noRd
.ge_canonical_origins <- function(raw) {
  if (nrow(raw) == 0L) return(raw)
  ord <- order(raw$id_genome_effect, raw$member_slot, raw$line_match_type,
               ifelse(is.na(raw$line_name), "\uFFFF", raw$line_name),
               ifelse(is.na(raw$parent_origin), .Machine$integer.max,
                      raw$parent_origin),
               raw$copy_count)
  raw <- raw[ord, , drop = FALSE]
  key <- paste(raw$id_genome_effect, raw$member_slot)
  raw$origin_slot <- as.integer(stats::ave(seq_len(nrow(raw)), key,
                                           FUN = seq_along))
  rownames(raw) <- NULL
  raw[, c("id_genome_effect", "member_slot", "origin_slot", "line_match_type",
          "line_name", "parent_origin", "copy_count"), drop = FALSE]
}


# ── Completeness ────────────────────────────────────────────────────────────

#' `require_complete`: every reachable `(copy_count, dosage)` state must appear
#'
#' Coverage is over the **state grid**, not over dosage values: at a locus whose
#' chromosome can be absent, `copy_count_value = 0` is a reachable state and a
#' surface that omits it is not complete.
#'
#' @keywords internal
#' @noRd
.ge_assert_complete <- function(conn, built) {
  members <- built$members
  if (any(members$contrast_name != "indicator")) {
    stop("require_complete = TRUE applies to indicator surfaces; this call ",
         "also contains ",
         paste0("'", unique(members$contrast_name[
           members$contrast_name != "indicator"]), "'", collapse = " and "),
         " members, whose contrasts are continuous functions of the genotype ",
         "and have no state grid to cover.", call. = FALSE)
  }
  counts <- .ge_reachable_copy_counts(conn)
  chr_of <- DBI::dbGetQuery(conn, "SELECT locus_id, chr_name FROM genome_meta")

  grp <- vapply(split(members, members$id_genome_effect), function(m) {
    paste(m$locus_id[order(m$locus_id)], collapse = ",")
  }, character(1))
  for (g in unique(grp)) {
    ids <- as.integer(names(grp)[grp == g])
    mm  <- members[members$id_genome_effect %in% ids, , drop = FALSE]
    loci <- sort(unique(mm$locus_id))
    want <- lapply(loci, function(lid) {
      chr <- chr_of$chr_name[match(lid, chr_of$locus_id)]
      do.call(rbind, lapply(counts[[chr]], function(cc) {
        data.frame(cc = cc, dv = seq.int(0L, cc))
      }))
    })
    grid <- expand.grid(lapply(want, function(w) paste(w$cc, w$dv, sep = "/")),
                        stringsAsFactors = FALSE)
    need <- apply(grid, 1, paste, collapse = " & ")
    have <- vapply(ids, function(id) {
      m <- mm[mm$id_genome_effect == id, , drop = FALSE]
      m <- m[order(m$locus_id), , drop = FALSE]
      paste(paste(m$copy_count_value, m$dosage_value, sep = "/"),
            collapse = " & ")
    }, character(1))
    missing <- setdiff(need, have)
    if (length(missing) > 0L) {
      nms <- mm$locus_name[match(loci, mm$locus_id)]
      stop("require_complete = TRUE: the surface over ",
           paste0("'", nms, "'", collapse = " x "), " is missing ",
           length(missing), " of ", length(need),
           " reachable (copy_count/dosage) state",
           if (length(need) == 1L) "" else "s", ": ",
           paste(head(missing, 5), collapse = "; "),
           if (length(missing) > 5L) ", ..." else "", ".", call. = FALSE)
    }
  }
  invisible(NULL)
}


# ── Replacement ─────────────────────────────────────────────────────────────

#' The single scope a `replace_scope` call keys on
#'
#' @keywords internal
#' @noRd
.ge_scope_from_origin <- function(origin, mode) {
  if (mode != "replace_scope") return(NULL)
  if (is.null(origin)) return(list(scope = "any"))
  if (is.data.frame(origin)) {
    stop("mode = \"replace_scope\" keys on one origin predicate, but 'origin' ",
         "is a per-member data frame with no single scope. Use ",
         "mode = \"replace_owner\", or pass the scope as a named scalar list.",
         call. = FALSE)
  }
  one <- .ge_origin_from_list(origin)
  list(scope = "scoped", one = one)
}

#' Does a stored term sit at exactly the supplied scope?
#'
#' Every member must carry that scope: `replace_scope` replaces a *uniform*
#' variant, which is what the scalar-list `origin` produces.
#'
#' @keywords internal
#' @noRd
.ge_term_at_scope <- function(mm, oo, scope) {
  for (i in seq_len(nrow(mm))) {
    so <- oo[oo$member_slot == mm$member_slot[i], , drop = FALSE]
    if (identical(scope$scope, "any")) {
      if (nrow(so) > 0L) return(FALSE)
      next
    }
    if (nrow(so) != 1L) return(FALSE)
    one <- scope$one
    if (!identical(so$line_match_type[1], one$line_match_type)) return(FALSE)
    if (!identical(.ge_line_token(so$line_match_type[1], so$line_name[1]),
                   .ge_line_token(one$line_match_type, one$line_name))) {
      return(FALSE)
    }
    if (!identical(is.na(so$parent_origin[1]), is.na(one$parent_origin))) {
      return(FALSE)
    }
    if (!is.na(one$parent_origin) &&
        so$parent_origin[1] != one$parent_origin) return(FALSE)
    if (so$copy_count[1] != one$copy_count) return(FALSE)
  }
  TRUE
}

#' Which stored term ids this call deletes
#'
#' @keywords internal
#' @noRd
.ge_resolve_deletes <- function(model, trait_name, effect_owner, mode, scope,
                                allow_reserved_owner) {
  if (mode == "append" || nrow(model$terms) == 0L) return(integer(0))
  t <- model$terms
  hit <- switch(
    mode,
    replace_trait = t$trait_name == trait_name,
    t$trait_name == trait_name & t$effect_owner == effect_owner
  )
  if (mode == "replace_trait" && !isTRUE(allow_reserved_owner)) {
    reserved <- hit & t$effect_owner %in% GE_RESERVED_OWNERS
    if (any(reserved)) {
      stop("mode = \"replace_trait\" would delete ", sum(reserved),
           " term(s) under the reserved owner '",
           paste(unique(t$effect_owner[reserved]), collapse = "', '"),
           "', which the generators define_additive_effects() and ",
           "define_genome_effects() own. Generated terms are replaced only by ",
           "re-running a generator (or removed with remove_generated_effects()); ",
           "use replace_owner or replace_scope to replace your own terms.",
           call. = FALSE)
    }
  }
  ids <- t$id_genome_effect[hit]
  if (mode != "replace_scope" || length(ids) == 0L) return(ids)

  ids[.ge_terms_at_scope(model, ids, scope)]
}

#' `.ge_term_at_scope()` for many stored terms, splitting the model once
#'
#' A per-id subset of the member and origin frames is quadratic in the size of
#' the stored model.
#'
#' @return Logical, one per `ids`.
#' @keywords internal
#' @noRd
.ge_terms_at_scope <- function(model, ids, scope) {
  if (length(ids) == 0L) return(logical(0))
  rows <- .ge_rows_by_id(model, ids)
  vapply(seq_along(ids), function(i) {
    .ge_term_at_scope(model$members[rows$members[[i]], , drop = FALSE],
                      model$origins[rows$origins[[i]], , drop = FALSE], scope)
  }, logical(1))
}

#' Member and origin row numbers of each of `ids`, from one pass over each
#'
#' @return `list(members, origins)`, each a list parallel to `ids`.
#' @keywords internal
#' @noRd
.ge_rows_by_id <- function(model, ids) {
  u   <- unique(ids)
  lev <- seq_along(u)
  by  <- function(x) {
    split(seq_along(x), factor(match(x, u), levels = lev))[match(ids, u)]
  }
  list(members = by(model$members$id_genome_effect),
       origins = by(model$origins$id_genome_effect))
}


# ── Transaction ─────────────────────────────────────────────────────────────

#' @keywords internal
#' @noRd
.ge_read_model <- function(conn) {
  list(terms   = DBI::dbGetQuery(conn, "SELECT * FROM genome_effects"),
       members = DBI::dbGetQuery(conn, "SELECT * FROM genome_effect_members"),
       origins = DBI::dbGetQuery(conn,
                                 "SELECT * FROM genome_effect_member_origins"))
}

#' @keywords internal
#' @noRd
.ge_require_effect_tables <- function(pop) {
  if (!all(c("genome_effects", "genome_effect_members",
             "genome_effect_member_origins") %in% pop$tables)) {
    stop("The genome-effect tables do not exist yet. Call define_genome() ",
         "before defining effects — genome_effect_members declares a foreign ",
         "key to genome_meta.locus_id.", call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.ge_require_trait <- function(conn, trait_name) {
  ok <- DBI::dbExistsTable(conn, "trait_meta") &&
    nrow(DBI::dbGetQuery(conn, paste0(
      "SELECT 1 FROM trait_meta WHERE trait_name IN (",
      sql_in_list(trait_name, what = "trait name"), ") LIMIT 1"))) > 0L
  if (!ok) {
    stop("Trait '", trait_name, "' not found in trait_meta. Call ",
         "define_trait() first.", call. = FALSE)
  }
  invisible(NULL)
}

#' Delete terms, children before parents
#'
#' DuckDB does not cascade, so the origin rows go first, then the members, then
#' the term — the reverse of the foreign-key direction.
#'
#' @keywords internal
#' @noRd
.ge_delete_terms <- function(conn, ids) {
  if (length(ids) == 0L) return(invisible(NULL))
  lst <- paste(as.integer(ids), collapse = ", ")
  for (tbl in c("genome_effect_member_origins", "genome_effect_members",
                "genome_effects")) {
    DBI::dbExecute(conn, paste0("DELETE FROM ", tbl,
                                " WHERE id_genome_effect IN (", lst, ")"))
  }
  invisible(NULL)
}

#' Apply one genome-effect write in a single transaction
#'
#' Deletes, then inserts, then validates the **whole** table set before COMMIT,
#' so a family conflict with rows written by an earlier call is caught here and
#' rolls the write back rather than being discovered at evaluation time.
#'
#' `built` may carry several traits' rows at once (the multi-trait
#' [define_additive_effects()] path), which is why ids are assigned here.
#'
#' `before_commit`, when given, is called with `conn` inside the transaction
#' after the whole-table validation and before `COMMIT`: the hook through which
#' a generator writes its target (`.tvc_write_block()`) atomically with the
#' terms. DuckDB has no nested `BEGIN`, so the target insert cannot open its
#' own.
#'
#' @keywords internal
#' @noRd
.ge_commit <- function(conn, delete_ids, built, before_commit = NULL) {
  n  <- nrow(built$terms)
  # Locally-indexed ids become real ones only now, so a build is reusable and
  # nothing depends on a MAX() read taken before the deletes.
  DBI::dbExecute(conn, "BEGIN TRANSACTION")
  ok <- FALSE
  # Every row-local rule is checked in R before this point, so a constraint
  # error here would be a package bug rather than bad input -- but roll back
  # regardless, or the connection is left aborted and every later query fails.
  on.exit({
    if (!ok) try(DBI::dbExecute(conn, "ROLLBACK"), silent = TRUE)
  }, add = TRUE)

  .ge_delete_terms(conn, delete_ids)
  start <- next_int_id(conn, "genome_effects", "id_genome_effect")
  remap <- stats::setNames(seq.int(start, length.out = n),
                           built$terms$id_genome_effect)
  labels <- if (is.null(built$labels)) NULL else
    stats::setNames(unname(built$labels), as.character(remap[names(built$labels)]))

  terms   <- built$terms
  members <- built$members
  origins <- built$origins
  terms$id_genome_effect   <- unname(remap[as.character(terms$id_genome_effect)])
  members$id_genome_effect <- unname(remap[as.character(members$id_genome_effect)])
  if (nrow(origins) > 0L) {
    origins$id_genome_effect <-
      unname(remap[as.character(origins$id_genome_effect)])
  }

  # Pre-SQL structural check, so the user's own term_id survives into the
  # message for everything the frames can decide on their own.
  v <- .ge_validate_frames(terms, members[, setdiff(names(members),
                                                    "locus_name")], origins,
                           labels = labels)
  if (length(v) > 0L) {
    stop("Invalid genome effects:\n  - ", paste(v, collapse = "\n  - "),
         call. = FALSE)
  }

  DBI::dbWriteTable(conn, "genome_effects", terms, append = TRUE)
  DBI::dbWriteTable(conn, "genome_effect_members",
                    members[, setdiff(names(members), "locus_name")],
                    append = TRUE)
  if (nrow(origins) > 0L) {
    DBI::dbWriteTable(conn, "genome_effect_member_origins", origins,
                      append = TRUE)
  }

  validate_genome_effects(conn, labels = labels)
  # A generator's new target (define_additive_effects(G = )) commits with its
  # terms or not at all: a failed call must not leave a target the model does
  # not deliver, nor terms calibrated to a target that was never stored.
  if (!is.null(before_commit)) before_commit(conn)
  DBI::dbExecute(conn, "COMMIT")
  ok <- TRUE
  invisible(unname(remap))
}

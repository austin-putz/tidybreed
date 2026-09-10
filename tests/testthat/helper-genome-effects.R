# ---------------------------------------------------------------------------
# Phase A fixtures for the term/member genome-effects model.
# Spec: plans/update_genome_effects_v4.md (v4.9). Progress notes and any
# design changes this phase forced: plans/update_genome_effects_phase_A.md.
#
# Nothing here touches the database. Phase A exists to answer two questions
# before a single table is created:
#
#   1. Is every case in the origin truth table *representable* in the three
#      proposed tables, and does it survive the row-local and cross-row rules?
#   2. Is its value *derivable from the stored rows alone*, and does that
#      value match a number computed by hand?
#
# Two independent evaluators are provided. `gefx_eval_naive()` implements the
# semantic definition literally -- Cartesian product over evaluation units,
# one variant selected per tuple (see "Contribution and resource guard").
# `gefx_eval_grouped()` implements the label-vector factorization that
# "Evaluation strategy" claims is the same function and is what Phase D will
# translate into SQL. Every fixture asserts naive == hand-computed == grouped,
# so the factorization is proven before any SQL depends on it.
#
# Loci are named, not id'd: locus_name -> locus_id resolution is the writer's
# job in Phase C and has no bearing on the semantics under test here.
# ---------------------------------------------------------------------------

# --- Genotype context ------------------------------------------------------

# One row per realized allele copy. An absent copy is an absent row, which is
# how the zero-copy state arises (fem_noY has no LY rows at all).
gefx_copies <- function() {
  spec <- c(
    # id_ind,  locus, parent_origin, line_origin, allele
    "pure_A,  L1, 1, A,  1",
    "pure_A,  L1, 2, A,  0",
    "pure_A,  L2, 1, A,  1",
    "pure_A,  L2, 2, A,  1",
    "pure_A,  L3, 1, A,  0",
    "pure_A,  L3, 2, A,  0",
    "pure_A,  L4, 1, A,  1",
    "pure_A,  L4, 2, A,  0",
    "pure_A,  LX, 1, A,  1",
    "pure_A,  LX, 2, A,  1",

    "f1_AB,   L1, 1, A,  1",
    "f1_AB,   L1, 2, B,  1",
    "f1_AB,   L2, 1, A,  0",
    "f1_AB,   L2, 2, B,  1",
    "f1_AB,   L3, 1, A,  1",
    "f1_AB,   L3, 2, B,  0",
    "f1_AB,   L4, 1, A,  1",
    "f1_AB,   L4, 2, B,  1",

    "f1_BA,   L1, 1, B,  1",
    "f1_BA,   L1, 2, A,  1",
    "f1_BA,   L2, 1, B,  1",
    "f1_BA,   L2, 2, A,  0",
    "f1_BA,   L3, 1, B,  0",
    "f1_BA,   L3, 2, A,  1",
    "f1_BA,   L4, 1, B,  1",
    "f1_BA,   L4, 2, A,  1",

    # One copy of known line, one copy whose founding line is unknown (NULL).
    "unk_A,   L1, 1, A,  1",
    "unk_A,   L1, 2, NA, 0",
    "unk_A,   L2, 1, A,  1",
    "unk_A,   L2, 2, NA, 1",
    "unk_A,   L3, 1, A,  1",
    "unk_A,   L3, 2, NA, 1",

    # Hemizygous carriers: one LX copy, from the dam (parent_origin 2).
    "male_X,  L1, 1, A,  1",
    "male_X,  L1, 2, A,  1",
    "male_X,  LX, 2, A,  1",
    "male_X0, L1, 1, A,  0",
    "male_X0, L1, 2, A,  0",
    "male_X0, LX, 2, A,  0",

    # Two LX copies, both allele 0 -> dosage 0 at copy count 2.
    "fem_XX0, L1, 1, A,  0",
    "fem_XX0, L1, 2, A,  0",
    "fem_XX0, LX, 1, A,  0",
    "fem_XX0, LX, 2, A,  0",

    # Carries no LY copy at all -> the synthesized zero-copy state.
    "fem_noY, L1, 1, A,  0",
    "fem_noY, L1, 2, A,  0"
  )
  parts <- do.call(rbind, lapply(strsplit(spec, ",", fixed = TRUE), trimws))
  data.frame(
    id_ind        = parts[, 1],
    locus_name    = parts[, 2],
    parent_origin = as.integer(parts[, 3]),
    line_origin   = ifelse(parts[, 4] == "NA", NA_character_, parts[, 4]),
    allele        = as.integer(parts[, 5]),
    stringsAsFactors = FALSE
  )
}

gefx_loci <- function() c("L1", "L2", "L3", "L4", "LX", "LY")

gefx_individuals <- function() unique(gefx_copies()$id_ind)

# --- Row constructors: the three proposed tables ---------------------------

gefx_term <- function(id, value, trait = "ADG", owner = "custom",
                      effect_name = NA_character_) {
  data.frame(id_genome_effect = as.integer(id), trait_name = trait,
             effect_owner = owner, effect_name = effect_name,
             genome_value = as.numeric(value), stringsAsFactors = FALSE)
}

gefx_member <- function(id, slot, locus, contrast, center = NA_real_,
                        copy_count_value = NA_integer_, dosage_value = NA_integer_) {
  data.frame(id_genome_effect = as.integer(id), member_slot = as.integer(slot),
             locus_name = locus, contrast_name = contrast,
             copy_count_value = as.integer(copy_count_value),
             dosage_value = as.integer(dosage_value),
             center_value = as.numeric(center), stringsAsFactors = FALSE)
}

gefx_origin <- function(id, slot, origin_slot, line_match_type,
                        line_name = NA_character_, parent_origin = NA_integer_,
                        copy_count = 1L) {
  data.frame(id_genome_effect = as.integer(id), member_slot = as.integer(slot),
             origin_slot = as.integer(origin_slot),
             line_match_type = line_match_type, line_name = line_name,
             parent_origin = as.integer(parent_origin),
             copy_count = as.integer(copy_count), stringsAsFactors = FALSE)
}

gefx_model <- function(terms, members, origins = NULL) {
  if (is.null(origins)) {
    origins <- gefx_origin(0L, 0L, 0L, "exact", "x", 1L)[0, , drop = FALSE]
  }
  list(terms = terms, members = members, origins = origins)
}

# --- Validation: the row-local CHECKs plus the R-enforced cross-row rules ---

# Returns a character vector of violations; empty means the model is writable.
gefx_validate <- function(model) {
  v <- character(0)
  mem <- model$members
  org <- model$origins

  for (i in seq_len(nrow(mem))) {
    m <- mem[i, ]
    tag <- paste0("member(", m$id_genome_effect, ",", m$member_slot, ")")
    if (!m$contrast_name %in% c("additive", "dominance", "indicator")) {
      v <- c(v, paste(tag, "contrast_name not in closed set"))
      next
    }
    if (m$contrast_name == "indicator") {
      if (is.na(m$copy_count_value) || is.na(m$dosage_value)) {
        v <- c(v, paste(tag, "indicator requires copy_count_value and dosage_value"))
      } else if (m$dosage_value > m$copy_count_value) {
        v <- c(v, paste(tag, "dosage_value > copy_count_value"))
      }
      if (!is.na(m$center_value)) v <- c(v, paste(tag, "indicator must not carry center_value"))
    } else {
      if (!is.na(m$copy_count_value) || !is.na(m$dosage_value)) {
        v <- c(v, paste(tag, "non-indicator must not carry a state"))
      }
      if (is.na(m$center_value)) {
        v <- c(v, paste(tag, "center_value IS NOT NULL required"))
      } else if (m$center_value < 0 || m$center_value > 1) {
        v <- c(v, paste(tag, "center_value outside [0, 1]"))
      }
    }
  }

  for (i in seq_len(nrow(org))) {
    o <- org[i, ]
    tag <- paste0("origin(", o$id_genome_effect, ",", o$member_slot, ",", o$origin_slot, ")")
    if (!o$line_match_type %in% c("exact", "unknown", "any")) {
      v <- c(v, paste(tag, "line_match_type not in closed set"))
      next
    }
    if (is.na(o$copy_count) || o$copy_count <= 0) v <- c(v, paste(tag, "copy_count must be > 0"))
    if (!is.na(o$parent_origin) && !o$parent_origin %in% c(1L, 2L)) {
      v <- c(v, paste(tag, "parent_origin not in (1, 2)"))
    }
    if (o$line_match_type == "exact" && is.na(o$line_name)) {
      v <- c(v, paste(tag, "'exact' requires line_name"))
    }
    if (o$line_match_type != "exact" && !is.na(o$line_name)) {
      v <- c(v, paste(tag, "only 'exact' may carry line_name"))
    }
    if (o$line_match_type == "any" && is.na(o$parent_origin)) {
      v <- c(v, paste(tag, "'any' requires a non-NULL parent_origin"))
    }
    owner_member <- mem[mem$id_genome_effect == o$id_genome_effect &
                          mem$member_slot == o$member_slot, , drop = FALSE]
    if (nrow(owner_member) != 1L) {
      v <- c(v, paste(tag, "orphan origin row"))
      next
    }
    if (o$line_match_type == "any" && owner_member$contrast_name != "additive") {
      v <- c(v, paste(tag, "'any' permitted only on additive members"))
    }
  }

  for (id in unique(mem$id_genome_effect)) {
    mm <- mem[mem$id_genome_effect == id, , drop = FALSE]
    if (anyDuplicated(mm$locus_name)) {
      v <- c(v, paste0("term(", id, ") repeats a locus"))
    }
    if (!nrow(mm)) v <- c(v, paste0("term(", id, ") has no members"))
    for (slot in mm$member_slot) {
      oo <- org[org$id_genome_effect == id & org$member_slot == slot, , drop = FALSE]
      kind <- mm$contrast_name[mm$member_slot == slot]
      if (kind == "additive" && nrow(oo) > 1L) {
        v <- c(v, paste0("member(", id, ",", slot, ") additive members take at most one origin row"))
      }
      if (kind == "additive" && nrow(oo) == 1L && oo$copy_count != 1L) {
        v <- c(v, paste0("member(", id, ",", slot, ") additive origin requires copy_count = 1"))
      }
    }
  }

  fam <- gefx_families(model)
  for (f in fam) {
    preds <- lapply(f$ids, function(id) gefx_variant_pred(model, id))
    if (length(f$ids) > 1L) {
      for (i in seq_along(f$ids)) for (j in seq_along(f$ids)) {
        if (j <= i) next
        le <- gefx_pred_leq(preds[[i]], preds[[j]])
        ge <- gefx_pred_leq(preds[[j]], preds[[i]])
        if (le && ge) {
          v <- c(v, paste0("duplicate family + scope identity: terms ",
                           f$ids[i], " and ", f$ids[j]))
        } else if (!le && !ge && gefx_pred_overlap(preds[[i]], preds[[j]])) {
          v <- c(v, paste0("overlapping but incomparable scopes in one family: terms ",
                           f$ids[i], " and ", f$ids[j]))
        }
      }
    }
  }
  unique(v)
}

# --- Families --------------------------------------------------------------

# Signature: (trait_name, effect_owner, ordered [(locus, contrast, cc, dosage)]).
# center_value, origin rows and genome_value are deliberately excluded.
gefx_family_key <- function(model, id) {
  tt <- model$terms[model$terms$id_genome_effect == id, ]
  mm <- model$members[model$members$id_genome_effect == id, , drop = FALSE]
  mm <- mm[order(mm$locus_name), , drop = FALSE]
  paste0(tt$trait_name, "|", tt$effect_owner, "|",
         paste(mm$locus_name, mm$contrast_name,
               ifelse(is.na(mm$copy_count_value), "-", mm$copy_count_value),
               ifelse(is.na(mm$dosage_value), "-", mm$dosage_value),
               sep = ":", collapse = "&"))
}

gefx_families <- function(model, trait = NULL) {
  ids <- model$terms$id_genome_effect
  if (!is.null(trait)) ids <- ids[model$terms$trait_name == trait]
  keys <- vapply(ids, function(id) gefx_family_key(model, id), character(1))
  lapply(split(ids, keys), function(g) list(key = keys[match(g[1], ids)], ids = g))
}

# --- Predicates ------------------------------------------------------------

# A variant's predicate is one per-member predicate, in member_slot order.
gefx_variant_pred <- function(model, id) {
  mm <- model$members[model$members$id_genome_effect == id, , drop = FALSE]
  mm <- mm[order(mm$member_slot), , drop = FALSE]
  lapply(seq_len(nrow(mm)), function(i) {
    slot <- mm$member_slot[i]
    kind <- if (mm$contrast_name[i] == "additive") "additive" else "genotype"
    oo <- model$origins[model$origins$id_genome_effect == id &
                          model$origins$member_slot == slot, , drop = FALSE]
    if (!nrow(oo)) return(list(kind = kind, scope = "any"))
    if (kind == "additive") {
      list(kind = kind, scope = "scoped",
           line = if (oo$line_match_type[1] == "exact") oo$line_name[1] else oo$line_match_type[1],
           parent = oo$parent_origin[1])
    } else {
      items <- do.call(rbind, lapply(seq_len(nrow(oo)), function(k) {
        data.frame(line = rep(oo$line_name[k], oo$copy_count[k]),
                   parent = rep(oo$parent_origin[k], oo$copy_count[k]),
                   stringsAsFactors = FALSE)
      }))
      list(kind = kind, scope = "scoped", items = items)
    }
  })
}

# Does a per-member predicate match a unit's label?
gefx_pred_match_unit <- function(p, unit) {
  if (identical(p$scope, "any")) return(TRUE)
  if (p$kind == "additive") {
    line_ok <- if (identical(p$line, "any")) TRUE
               else if (identical(p$line, "unknown")) is.na(unit$line)
               else !is.na(unit$line) && unit$line == p$line
    parent_ok <- is.na(p$parent) || unit$parent == p$parent
    return(line_ok && parent_ok)
  }
  gefx_bijection(p$items, unit$items)
}

# Exact-multiset satisfaction: every demand consumed by a distinct copy.
# Full search rather than greedy -- with a demand carrying an ANY parent
# alongside a parent-qualified one, greedy can fail a satisfiable multiset.
gefx_bijection <- function(demands, copies) {
  if (nrow(demands) != nrow(copies)) return(FALSE)
  if (!nrow(demands)) return(TRUE)
  d <- demands[1, ]
  for (i in seq_len(nrow(copies))) {
    cc <- copies[i, ]
    line_ok <- !is.na(cc$line) && !is.na(d$line) && cc$line == d$line
    parent_ok <- is.na(d$parent) || cc$parent == d$parent
    if (line_ok && parent_ok &&
        gefx_bijection(demands[-1, , drop = FALSE], copies[-i, , drop = FALSE])) {
      return(TRUE)
    }
  }
  FALSE
}

# Containment, componentwise across members. P <= Q means "P is at least as
# specific as Q", so the winner among matching variants is the minimum.
gefx_pred_leq <- function(P, Q) {
  if (length(P) != length(Q)) return(FALSE)
  all(vapply(seq_along(P), function(i) gefx_member_leq(P[[i]], Q[[i]]), logical(1)))
}

gefx_member_leq <- function(p, q) {
  if (identical(q$scope, "any")) return(TRUE)
  if (identical(p$scope, "any")) return(FALSE)
  if (p$kind == "additive") {
    line_ok <- identical(q$line, "any") || identical(p$line, q$line)
    parent_ok <- is.na(q$parent) || (!is.na(p$parent) && p$parent == q$parent)
    return(line_ok && parent_ok)
  }
  # Same line multiset, and P's parent assignments refine Q's.
  if (nrow(p$items) != nrow(q$items)) return(FALSE)
  if (!identical(sort(p$items$line), sort(q$items$line))) return(FALSE)
  gefx_refines(p$items, q$items)
}

gefx_refines <- function(pi, qi) {
  if (!nrow(pi)) return(TRUE)
  d <- pi[1, ]
  for (i in seq_len(nrow(qi))) {
    q <- qi[i, ]
    if (identical(q$line, d$line) &&
        (is.na(q$parent) || (!is.na(d$parent) && d$parent == q$parent)) &&
        gefx_refines(pi[-1, , drop = FALSE], qi[-i, , drop = FALSE])) {
      return(TRUE)
    }
  }
  FALSE
}

# Two predicates overlap if some reachable unit label satisfies both. Decided
# by enumeration over the fixture universe, which is what makes it decidable
# here without a general satisfiability argument.
gefx_pred_overlap <- function(P, Q) {
  if (length(P) != length(Q)) return(FALSE)
  copies <- gefx_copies()
  all(vapply(seq_along(P), function(i) {
    p <- P[[i]]; q <- Q[[i]]
    units <- gefx_universe_units(copies, p$kind)
    any(vapply(units, function(u) gefx_pred_match_unit(p, u) && gefx_pred_match_unit(q, u),
               logical(1)))
  }, logical(1)))
}

gefx_universe_units <- function(copies, kind) {
  if (kind == "additive") {
    lab <- unique(copies[, c("line_origin", "parent_origin")])
    return(lapply(seq_len(nrow(lab)), function(i) {
      list(line = lab$line_origin[i], parent = lab$parent_origin[i])
    }))
  }
  keys <- unique(paste(copies$id_ind, copies$locus_name))
  out <- lapply(keys, function(k) {
    parts <- strsplit(k, " ", fixed = TRUE)[[1]]
    rows <- copies[copies$id_ind == parts[1] & copies$locus_name == parts[2], , drop = FALSE]
    list(items = data.frame(line = rows$line_origin, parent = rows$parent_origin,
                            stringsAsFactors = FALSE))
  })
  c(out, list(list(items = data.frame(line = character(0), parent = integer(0),
                                      stringsAsFactors = FALSE))))
}

# --- Evaluation units ------------------------------------------------------

gefx_units <- function(copies, id_ind, locus, kind) {
  rows <- copies[copies$id_ind == id_ind & copies$locus_name == locus, , drop = FALSE]
  rows <- rows[order(rows$parent_origin), , drop = FALSE]
  if (kind == "additive") {
    return(lapply(seq_len(nrow(rows)), function(i) {
      list(label = paste0(ifelse(is.na(rows$line_origin[i]), "NA", rows$line_origin[i]),
                          "@", rows$parent_origin[i]),
           line = rows$line_origin[i], parent = rows$parent_origin[i],
           allele = rows$allele[i])
    }))
  }
  # Genotype grain: exactly one unit per locus, including the zero-copy state,
  # which has no row to join and must be materialized.
  items <- data.frame(line = rows$line_origin, parent = rows$parent_origin,
                      stringsAsFactors = FALSE)
  list(list(
    label = paste(sort(paste0(ifelse(is.na(items$line), "NA", items$line),
                              "@", items$parent)), collapse = ","),
    items = items,
    copy_count = nrow(rows),
    dosage = sum(rows$allele)
  ))
}

gefx_contrast_value <- function(unit, member) {
  switch(member$contrast_name,
    additive  = unit$allele - member$center_value,
    dominance = {
      p <- member$center_value; q <- 1 - p
      if (unit$copy_count != 2L) {
        stop("dominance evaluated at a non-diploid state -- rejected at write time")
      }
      c(-2 * p^2, 2 * p * q, -2 * q^2)[unit$dosage + 1L]
    },
    indicator = as.numeric(unit$copy_count == member$copy_count_value &&
                             unit$dosage == member$dosage_value),
    stop("unknown contrast")
  )
}

gefx_member_kind <- function(contrast) if (contrast == "additive") "additive" else "genotype"

# --- Evaluator 1: the semantic definition, per tuple -----------------------

gefx_eval_naive <- function(model, id_ind, trait = "ADG", copies = gefx_copies()) {
  total <- 0
  for (fam in gefx_families(model, trait)) {
    ref <- model$members[model$members$id_genome_effect == fam$ids[1], , drop = FALSE]
    ref <- ref[order(ref$member_slot), , drop = FALSE]
    units <- lapply(seq_len(nrow(ref)), function(i) {
      gefx_units(copies, id_ind, ref$locus_name[i], gefx_member_kind(ref$contrast_name[i]))
    })
    if (any(lengths(units) == 0L)) next          # no eligible copy anywhere
    preds <- lapply(fam$ids, function(id) gefx_variant_pred(model, id))
    grid <- expand.grid(lapply(units, seq_along), KEEP.OUT.ATTRS = FALSE)
    for (r in seq_len(nrow(grid))) {
      tuple <- lapply(seq_along(units), function(m) units[[m]][[grid[r, m]]])
      pick <- gefx_select_variant(preds, tuple, fam$ids)
      if (is.na(pick)) next                      # unmatched tuple contributes 0
      mm <- model$members[model$members$id_genome_effect == pick, , drop = FALSE]
      mm <- mm[order(mm$member_slot), , drop = FALSE]
      prod <- 1
      for (m in seq_along(tuple)) prod <- prod * gefx_contrast_value(tuple[[m]], mm[m, ])
      total <- total + model$terms$genome_value[model$terms$id_genome_effect == pick] * prod
    }
  }
  total
}

gefx_select_variant <- function(preds, tuple, ids) {
  hit <- which(vapply(preds, function(P) {
    all(vapply(seq_along(tuple), function(m) gefx_pred_match_unit(P[[m]], tuple[[m]]),
               logical(1)))
  }, logical(1)))
  if (!length(hit)) return(NA_integer_)
  minimal <- hit[vapply(hit, function(i) {
    !any(vapply(setdiff(hit, i), function(j) gefx_pred_leq(preds[[j]], preds[[i]]), logical(1)))
  }, logical(1))]
  if (length(minimal) != 1L) {
    stop("ambiguous specificity among terms ", paste(ids[minimal], collapse = ", "),
         " -- must have been rejected at write time")
  }
  ids[minimal]
}

# --- Evaluator 2: the label-vector factorization ---------------------------

gefx_eval_grouped <- function(model, id_ind, trait = "ADG", copies = gefx_copies()) {
  total <- 0
  for (fam in gefx_families(model, trait)) {
    ref <- model$members[model$members$id_genome_effect == fam$ids[1], , drop = FALSE]
    ref <- ref[order(ref$member_slot), , drop = FALSE]
    units <- lapply(seq_len(nrow(ref)), function(i) {
      gefx_units(copies, id_ind, ref$locus_name[i], gefx_member_kind(ref$contrast_name[i]))
    })
    if (any(lengths(units) == 0L)) next
    labels <- lapply(units, function(u) unique(vapply(u, function(x) x$label, character(1))))
    preds <- lapply(fam$ids, function(id) gefx_variant_pred(model, id))
    grid <- expand.grid(lapply(labels, seq_along), KEEP.OUT.ATTRS = FALSE)
    for (r in seq_len(nrow(grid))) {
      rep_unit <- lapply(seq_along(units), function(m) {
        lab <- labels[[m]][grid[r, m]]
        Filter(function(x) x$label == lab, units[[m]])[[1]]
      })
      pick <- gefx_select_variant(preds, rep_unit, fam$ids)
      if (is.na(pick)) next
      mm <- model$members[model$members$id_genome_effect == pick, , drop = FALSE]
      mm <- mm[order(mm$member_slot), , drop = FALSE]
      prod <- 1
      for (m in seq_along(units)) {
        lab <- labels[[m]][grid[r, m]]
        same <- Filter(function(x) x$label == lab, units[[m]])
        prod <- prod * sum(vapply(same, function(u) gefx_contrast_value(u, mm[m, ]), numeric(1)))
      }
      total <- total + model$terms$genome_value[model$terms$id_genome_effect == pick] * prod
    }
  }
  total
}

# --- Fixture registry ------------------------------------------------------
#
# Each fixture carries:
#   description        what the case proves
#   model              the stored rows (the three proposed tables)
#   expected           hand-computed value per individual, NULL when the
#                      fixture is a rejection case
#   expect_invalid     regex the validator's complaint must match
#   expect_eval_error  regex an evaluation-time stop must match
#   note               anything the case forced or exposed
#
# Hand computations live in plans/update_genome_effects_phase_A.md, one
# worked line per number below. Nothing here is generated from an
# implementation; that is the point of doing it before the DDL exists.

gefx_fixtures <- function() list(

  F01 = list(
    description = "common vs A-specific additive: per-copy line fallback (today's engine)",
    model = gefx_model(
      rbind(gefx_term(1, 2.0), gefx_term(2, 3.0)),
      rbind(gefx_member(1, 1, "L1", "additive", center = 0.5),
            gefx_member(2, 1, "L1", "additive", center = 0.4)),
      gefx_origin(2, 1, 1, "exact", "A")
    ),
    expected = c(pure_A = 0.6, f1_AB = 2.8, f1_BA = 2.8, unk_A = 0.8)
  ),

  F02 = list(
    description = "generic A vs paternal A: strict subset wins, no equal-row-count tie",
    model = gefx_model(
      rbind(gefx_term(1, 1.0), gefx_term(2, 5.0)),
      rbind(gefx_member(1, 1, "L2", "additive", center = 0.4),
            gefx_member(2, 1, "L2", "additive", center = 0.4)),
      rbind(gefx_origin(1, 1, 1, "exact", "A"),
            gefx_origin(2, 1, 1, "exact", "A", parent_origin = 1L))
    ),
    expected = c(pure_A = 3.6, f1_AB = -2.0, f1_BA = -0.4)
  ),

  F03 = list(
    description = "common vs A/B dominance: exact multiset beats the common scope",
    model = gefx_model(
      rbind(gefx_term(1, 4.0), gefx_term(2, 10.0)),
      rbind(gefx_member(1, 1, "L1", "dominance", center = 0.5),
            gefx_member(2, 1, "L1", "dominance", center = 0.5)),
      rbind(gefx_origin(2, 1, 1, "exact", "A"),
            gefx_origin(2, 1, 2, "exact", "B"))
    ),
    expected = c(pure_A = 2.0, f1_AB = -5.0, f1_BA = -5.0, unk_A = 2.0)
  ),

  F04 = list(
    description = "generic A/B vs one reciprocal: the other direction falls back",
    model = gefx_model(
      rbind(gefx_term(1, 2.0), gefx_term(2, 7.0)),
      rbind(gefx_member(1, 1, "L3", "dominance", center = 0.5),
            gefx_member(2, 1, "L3", "dominance", center = 0.5)),
      rbind(gefx_origin(1, 1, 1, "exact", "A"),
            gefx_origin(1, 1, 2, "exact", "B"),
            gefx_origin(2, 1, 1, "exact", "A", parent_origin = 1L),
            gefx_origin(2, 1, 2, "exact", "B", parent_origin = 2L))
    ),
    expected = c(f1_AB = 3.5, f1_BA = 1.0, pure_A = 0.0, unk_A = 0.0)
  ),

  F05 = list(
    description = "overlapping but incomparable scopes in one family are rejected",
    model = gefx_model(
      rbind(gefx_term(1, 1.0), gefx_term(2, 2.0)),
      rbind(gefx_member(1, 1, "L1", "additive", center = 0.5),
            gefx_member(2, 1, "L1", "additive", center = 0.5)),
      rbind(gefx_origin(1, 1, 1, "exact", "A"),
            gefx_origin(2, 1, 1, "any", parent_origin = 1L))
    ),
    expect_invalid = "overlapping but incomparable"
  ),

  F06 = list(
    description = "common vs origin-specific A x A: the selected variant varies by tuple",
    model = gefx_model(
      rbind(gefx_term(1, 1.0), gefx_term(2, 4.0)),
      rbind(gefx_member(1, 1, "L1", "additive", center = 0.5),
            gefx_member(1, 2, "L2", "additive", center = 0.5),
            gefx_member(2, 1, "L1", "additive", center = 0.4),
            gefx_member(2, 2, "L2", "additive", center = 0.4)),
      rbind(gefx_origin(2, 1, 1, "exact", "A"),
            gefx_origin(2, 2, 1, "exact", "A"))
    ),
    expected = c(f1_AB = -0.71, pure_A = 0.96),
    note = paste("The reason the sum does not factor into a plain product of",
                 "per-member sums: different tuples select different variants.")
  ),

  F07 = list(
    description = "partial specificity at one member of two",
    model = gefx_model(
      rbind(gefx_term(1, 1.0), gefx_term(2, 3.0)),
      rbind(gefx_member(1, 1, "L1", "additive", center = 0.5),
            gefx_member(1, 2, "L4", "additive", center = 0.5),
            gefx_member(2, 1, "L1", "additive", center = 0.5),
            gefx_member(2, 2, "L4", "additive", center = 0.5)),
      gefx_origin(2, 1, 1, "exact", "A")
    ),
    expected = c(f1_AB = 2.0),
    note = "All-common would be 1.0 and all-scoped 3.0, so 2.0 separates all three."
  ),

  F08 = list(
    description = "two disjoint specific scopes both contribute, on different copies",
    model = gefx_model(
      rbind(gefx_term(1, 2.0), gefx_term(2, 5.0)),
      rbind(gefx_member(1, 1, "L1", "additive", center = 0.5),
            gefx_member(2, 1, "L1", "additive", center = 0.5)),
      rbind(gefx_origin(1, 1, 1, "exact", "A", parent_origin = 1L),
            gefx_origin(2, 1, 1, "exact", "B", parent_origin = 2L))
    ),
    expected = c(f1_AB = 3.5, f1_BA = 0.0, pure_A = 1.0),
    note = paste("f1_BA scores 0 with every copy unmatched -- the case the",
                 "trait-level no-match diagnostic must still catch.")
  ),

  F09 = list(
    description = "indicator states separate copy counts 0, 1 and 2 at dosage zero",
    model = gefx_model(
      rbind(gefx_term(1, 10.0), gefx_term(2, 20.0), gefx_term(3, 30.0)),
      rbind(gefx_member(1, 1, "LX", "indicator", copy_count_value = 0L, dosage_value = 0L),
            gefx_member(2, 1, "LX", "indicator", copy_count_value = 1L, dosage_value = 0L),
            gefx_member(3, 1, "LX", "indicator", copy_count_value = 2L, dosage_value = 0L))
    ),
    expected = c(fem_noY = 10.0, male_X0 = 20.0, fem_XX0 = 30.0,
                 male_X = 0.0, pure_A = 0.0),
    note = paste("Dosage alone conflates all three. fem_noY has no LX row at",
                 "all, so the zero-copy state has to be materialized.")
  ),

  F10 = list(
    description = "imprinting as an origin predicate: ('any', parent 1)",
    model = gefx_model(
      gefx_term(1, 2.0),
      gefx_member(1, 1, "L2", "additive", center = 0.5),
      gefx_origin(1, 1, 1, "any", parent_origin = 1L)
    ),
    expected = c(pure_A = 1.0, f1_AB = -1.0, f1_BA = 1.0),
    note = "Numerically the removed expressed_parent = 'parent_1' filter."
  ),

  F11 = list(
    description = "line-specific imprinting: the wrapper must compose line with parent",
    model = gefx_model(
      rbind(gefx_term(1, 1.0), gefx_term(2, 4.0)),
      rbind(gefx_member(1, 1, "L1", "additive", center = 0.4),
            gefx_member(2, 1, "L1", "additive", center = 0.6)),
      rbind(gefx_origin(1, 1, 1, "exact", "A", parent_origin = 1L),
            gefx_origin(2, 1, 1, "exact", "B", parent_origin = 1L))
    ),
    expected = c(f1_AB = 0.6, f1_BA = 1.6),
    note = paste("Stamping 'any' instead of composing would collapse both",
                 "variants onto one scope and lose a per-line centre.")
  ),

  F12 = list(
    description = "{A, NULL}: exact, unknown and common are three disjoint-or-nested scopes",
    model = gefx_model(
      rbind(gefx_term(1, 2.0), gefx_term(2, 3.0), gefx_term(3, 7.0)),
      rbind(gefx_member(1, 1, "L1", "additive", center = 0.5),
            gefx_member(2, 1, "L1", "additive", center = 0.4),
            gefx_member(3, 1, "L1", "additive", center = 0.5)),
      rbind(gefx_origin(2, 1, 1, "exact", "A"),
            gefx_origin(3, 1, 1, "unknown"))
    ),
    expected = c(unk_A = -1.7, f1_AB = 2.8, pure_A = 0.6)
  ),

  F13 = list(
    description = "hemizygous additive centring uses the realized copy count",
    model = gefx_model(
      gefx_term(1, 2.0),
      gefx_member(1, 1, "LX", "additive", center = 0.5)
    ),
    expected = c(male_X = 1.0, male_X0 = -1.0, pure_A = 2.0, fem_noY = 0.0),
    note = "male_X is dosage - c, not dosage - 2c. fem_noY has no eligible copy."
  ),

  F14 = list(
    description = "dominance at a non-diploid state stops rather than scoring",
    model = gefx_model(
      gefx_term(1, 3.0),
      gefx_member(1, 1, "LX", "dominance", center = 0.5)
    ),
    eval_ind = "male_X",
    expect_eval_error = "non-diploid",
    note = paste("In the real system this never reaches evaluation: the write",
                 "is refused against chr_inheritance. Phase A has no such table,",
                 "so the case is pinned at the evaluator instead.")
  ),

  F15 = list(
    description = "duplicate family + scope identity is rejected",
    model = gefx_model(
      rbind(gefx_term(1, 2.0), gefx_term(2, 1.0)),
      rbind(gefx_member(1, 1, "L1", "additive", center = 0.5),
            gefx_member(2, 1, "L1", "additive", center = 0.5))
    ),
    expect_invalid = "duplicate family \\+ scope identity"
  ),

  F16 = list(
    description = "mixed centring at one locus collides as a duplicate (Q2)",
    model = gefx_model(
      rbind(gefx_term(1, 2.0), gefx_term(2, 1.0)),
      rbind(gefx_member(1, 1, "L1", "additive", center = 0.5),
            gefx_member(2, 1, "L1", "additive", center = 0.2))
    ),
    expect_invalid = "duplicate family \\+ scope identity",
    combined = c(genome_value = 3.0, center_value = 0.4),
    note = "center_value is outside the family signature, so these collide."
  ),

  F17 = list(
    description = "'any' without a parent is the common scope written the long way",
    model = gefx_model(
      gefx_term(1, 2.0),
      gefx_member(1, 1, "L1", "additive", center = 0.5),
      gefx_origin(1, 1, 1, "any")
    ),
    expect_invalid = "'any' requires a non-NULL parent_origin"
  ),

  F18 = list(
    description = "'any' on a genotype member constrains nothing and is rejected",
    model = gefx_model(
      gefx_term(1, 2.0),
      gefx_member(1, 1, "L1", "dominance", center = 0.5),
      gefx_origin(1, 1, 1, "any", parent_origin = 1L)
    ),
    expect_invalid = "permitted only on additive members"
  ),

  F19 = list(
    description = "additive main, dominance main and an interaction at one locus all fire",
    model = gefx_model(
      rbind(gefx_term(1, 2.0), gefx_term(2, 4.0), gefx_term(3, 1.0)),
      rbind(gefx_member(1, 1, "L1", "additive", center = 0.5),
            gefx_member(2, 1, "L1", "dominance", center = 0.5),
            gefx_member(3, 1, "L1", "additive", center = 0.5),
            gefx_member(3, 2, "L4", "additive", center = 0.5))
    ),
    expected = c(f1_AB = 1.0),
    note = "2.0 additive + (-2.0) dominance + 1.0 interaction; three families, never competing."
  )
)

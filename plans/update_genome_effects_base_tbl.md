# Base population as a filtered table — `base_tbl` for the genome-effect writers

**Status:** **IMPLEMENTED** — v0.69.0, 2026-09-20, commits `48ebce2` (step 1),
`6a7d3a5` (step 2), `845d1cd` (step 3) + a docs commit. Deviations and the
as-executed call-site inventory are in
`plans/update_genome_effects_base_tbl_implementation.md`; `> Implemented:`
notes below mark where the code differs from the sketch. **Created:** 2026-09-20.
**Revision:** v5 (2026-09-20) — final pre-implementation pass. **Approved for
implementation** by the Codex v4 review
(`plans/update_genome_effects_base_tbl_codex_v1.md`, "Update after v4"); Q10 and
Q11 accepted; two-function split endorsed.
v3 — multi-line / crossbreeding walkthrough (§4). v2 — after the first Codex review.
**Baseline:** v0.68.4 on `feat/genome-effects-v49`.
**Relates to:** `plans/update_genome_effects.md` (the v2 storage architecture) and
the v4 / v4.9 series. This plan touches only how the writers *choose the population
that defines allele frequencies*; storage and evaluation are unchanged.

## Revision history

**v4 → v5 (final pass).** Codex's v4 review approved implementation and the
two-function split with one test clarification, accepted: test 22 now **omits
`center_value`** so the writer half must fill from `base_tbl` — copying the
stored centres would leave the base validated but unqueried and prove only
storage-path parity. A full read-through of the plan then fixed six things an
implementer would have tripped on: the §4.5 dominance example lacked
`copy_count = 2L` (a scalar `origin` list scopes additive members but is refused
for a genotype member, `R/define_genome_effects.R:52-56`);
`.founder_base_empty_error()` was sketched with two different signatures;
the `NA`-in-`p_base` check was not pinned to run *before* the two
`rescale_effects_to_target()` call sites that would fold `NA` into `V_A`;
§2.5 did not say that `base_tbl = NULL` deliberately differs between the verbs;
the §4.1 comment and test 16 still said "no base arguments" for a script with
one explicit base; and test 22 said "fresh population" where "same population,
`replace_scope`" is both simpler and guarantees an identical base.

**v4 addendum (same day).** Merging the two writers was considered and rejected;
§2.5 records why, states that `base_tbl` means the same on both, adds a
one-sentence family statement to their roxygen, and pins "generator ≡ writer"
as acceptance test 22. The scope/`mode` vocabulary mismatch between the verbs is
logged in §3.5 as a separate change, not folded in here.

**v3 → v4.** The Codex re-review approved the design and Q10/Q11, with three
amendments — all verified against the tree and accepted:

1. **`line_name` dropped from the `founder_haplotypes` column requirement** (§2.3,
   §3.1). Once the Wahlund check moved onto the default path (v3), nothing in
   `extract_allele_freq()` reads `line_name` from the rendered subquery; v3 kept the
   requirement by inertia and would have rejected a valid
   `select(locus_name, allele)` base.
2. **The `needs_fill` test in §3.4 was wrong for the plan's own showcase example.**
   An omitted `center_value` column makes `terms$center_value` `NULL`, and
   `any(NULL & …)` is `FALSE` — no query, then `.ge_build()` (which creates the
   column as `NA_real_`, `R/define_genome_effects.R:245-249`) raises the old
   error. The decision now normalises the column exactly as `.ge_build()` does,
   and test 18 covers both an explicit `NA` and an omitted column.
3. **Acceptance test 17 was non-deterministic.** Realized `line_origin == "Duroc"`
   copies are a finite transmitted sample of the Duroc pool and equal its
   frequency only in expectation. The test now compares against a hand-computed
   `AVG(allele)` over exactly the selected realized copies, on a fixed fixture.

Three wording fixes from the same review: the §2.2 contract states both halves
(partial coverage → per-locus `NA`; no coverage → error); §Q10 says "at most two
`EXISTS` queries", not one; test 16 no longer calls the explicit pooled common
write a "default-base call". One Codex note recorded for the implementer under
Q10: the fallback is **pool-level, never per-locus** — an incomplete named pool is
not stitched from the shared pool; the QTL-level missing-copy check (§3.3) is
what catches that.

**v2 → v3.** A walkthrough of every multi-line / crossbreeding configuration the
package supports (§4) found the design sound but three details wrong or unsaid:

1. **The default base now falls back to the shared pool.** A shared
   `founder_haplotypes` pool (`line_name = NULL`) feeding several named founder
   lines is a legitimate design — one base, lines diverge by selection — and v2's
   default `filter(line_name == "A")` would error on it (today's code does too).
   The package's own precedence rule (`resolve_genome_map()`,
   `resolve_chr_inheritance()`: `line = L → line = NULL`) is applied to the
   default instead (§3.3, Q10).
2. **The Wahlund warning moves out of `extract_allele_freq()` onto the default
   path**, and counts the `NULL` pool as a group. An explicit
   `base_tbl = get_table(pop, "founder_haplotypes")` is the user *asking* to pool —
   it is how the common fallback variant is defined in every crossbreeding
   program — and must not warn on every call (§3.2, §3.3, Q11).
3. **`ind_meta.line_name` vs `ind_haplotype.line_origin` are different things**
   and the plan now says which base to reach for when (§2.3, §4.2).

**v1 → v2.** The Codex review approved the direction and every open-question
recommendation, and was accepted in full on its five required revisions. Each was
checked against the tree before acceptance:

1. **One shared `base_tbl` validator** (§3.1). `select.tidybreed_table()` swaps the
   lazy projection and keeps the wrapper, so dispatching on `table_name` alone can
   hand the user a raw DuckDB error. Required columns are now part of the contract,
   and the same-connection check applies to both writers.
2. **The `define_genome_effects()` fill lives inside `.ge_build()`** (§3.4).
   `.ge_build()` builds `members` and calls `.ge_check_member_fields()` with no seam
   between (`R/define_genome_effects.R:304-321`); v1's "fill before validation"
   could not be bolted on from outside.
3. **Two overclaims corrected.** v1 said the generic empty check preserved the
   "no rows for line X, available: …" diagnostic — it cannot; the founder-specific
   diagnostic is now retained explicitly (§3.2). v1 counted 26 call sites; a wider
   search finds 73, including **six executable** `base = "current_pop"` calls in the
   swine vignette that would silently flip to the founder default if the argument
   were merely deleted (§3.6).
4. **SQL precision** (§3.2): `sql_render()`, a semi-join on distinct ids, the
   all-`NA` check named for what it proves, and honesty that the Wahlund diagnostic
   is a second statement.
5. **Acceptance tests expanded** (§3.8) for repeated individual rows, columns
   removed by `select()`, cross-population rejection in both writers, no query
   when every centre is explicit, and `0`/`1` frequencies distinct from `NA`.

Also from the review: the Q8 rationale for keeping `genome_meta.founder_allele_freq`
was wrong and is replaced (it is last-call-wins in a multi-line population, not a
snapshot); Q6 no longer implies `ad_terms(p = NA)`, which `.ad_recycle()` rejects;
Q9's claim that rendered SQL is "simply wrong" on a second connection is softened.
One refinement added beyond the review: `base_tbl` is **validated always, queried
lazily** — a wrong object errors whether or not any centre happens to need filling.

---

## Summary of decisions

1. **`define_additive_effects()` loses `base`, `base_tbl`-as-current-pop-only, and
   `base_line_name`, and gains one argument: `base_tbl`, a filtered
   `tidybreed_table`.** The table's identity says *what kind of thing* is being
   selected (a founder pool, allele copies, or individuals); the user's `filter()`
   says *which ones*. There is no enum.
2. **The pipe subject stays `genome_meta`.** The first argument already answers
   "which loci become QTL"; the population that centres them is an orthogonal
   selection and comes in by name. This is the same two-table shape `add_ebv()`
   already uses (`tbl` = which animals, `phenotype =` = which records).
3. **Three table kinds are accepted**, dispatched on `base_tbl$table_name` **and**
   checked for the columns the generated SQL touches: `founder_haplotypes` (the
   pool), `ind_haplotype` (these allele copies), and any other table carrying
   `id_ind` (these individuals, resolved through `ind_haplotype`).
4. **The default is the founder pool of the line the effect applies to, with the
   package's standard `line → NULL` fallback.** `base_tbl = NULL` means
   `founder_haplotypes |> filter(line_name == <line_name>)` when that pool exists;
   the shared pool (`line_name IS NULL`) when the effect is line-scoped but only a
   shared pool exists; and the whole founder table (with the Wahlund warning) when
   the effect is population-wide. What changes from today is how you *override*
   it: "force pooling for a line-A effect" becomes the explicit
   `base_tbl = get_table(pop, "founder_haplotypes")` — which does **not** warn,
   because pooling was asked for — instead of the sentinel `base_line_name = NULL`
   and a `missing()` check.
5. **One shared, exported helper: `extract_allele_freq(tbl)`.** It is the single
   place per-locus allele frequency is computed from a `tidybreed_table`, used by
   both writers internally and by users to obtain `p` for `ad_terms()`. It follows
   the `extract_` contract — returns analysis data, changes no state.
6. **`define_genome_effects()` gains the same `base_tbl`,** used only to fill
   `center_value` where the user left it `NA` on an `additive` or `dominance`
   member. An explicit `center_value` always wins; `indicator` members are never
   touched. The fill happens inside `.ge_build()`, before member validation.
7. **The frequency is computed in one SQL statement with the user's filter as a
   subquery.** Nothing is collected into R except the per-locus frequencies. This
   deliberately improves on the current `current_pop` path (collect ids → paste an
   `IN (...)` list) and on `add_ebv()`'s `phenotype` handling (collect the whole
   filtered table to pull out ids). The default-path Wahlund check and the
   founder empty-selection diagnostic are separate, cheap statements.

---

## 1. Why

### 1.1 What exists today

`define_additive_effects()` answers one question — *which allele copies define
`p`?* — with three arguments that constrain each other:

| Argument | Honoured when | Ignored / errors when |
|---|---|---|
| `base = c("founder_haplotypes", "current_pop")` | always | — |
| `base_tbl` (a `tidybreed_table` with `id_ind`) | `base = "current_pop"` | `founder_haplotypes` → warning, ignored |
| `base_line_name` | `base = "founder_haplotypes"` | `current_pop` → error if explicitly supplied |

`R/define_additive_effects.R:223-248` is ~25 lines whose only job is to keep the
three from contradicting each other, including a `missing()` check to distinguish
"`base_line_name` not supplied → inherit `line_name`" from "`base_line_name =
NULL` supplied → force pooling". `compute_base_allele_freq()`
(`R/define_additive_effects.R:820-915`) is two hand-rolled SQL paths; the
`current_pop` path collects `id_ind` into R and pastes it back as an `IN (...)`
list via `sql_in_list()`. Its zero-initialised output vector silently centres any
locus absent from the base at `p = 0`.

`define_genome_effects()` has no base argument at all: `center_value` is supplied
in `terms`, and `ad_terms()` refuses to run without an explicit `p`. There is no
sanctioned way to *get* `p` from the database — users must write their own
`AVG(allele)` query.

### 1.2 The package already has the pattern

`add_founders()` takes `get_table("founder_haplotypes") |> filter(line_name ==
"A")` as its pipe subject. `add_ebv()` takes a second `tidybreed_table` by name:

```r
add_ebv(tbl,                 # which animals enter the evaluation
        trait_name,
        ...,
        phenotype = get_table(pop, "ind_phenotype") |> filter(...))   # which records
```

validated at `R/add_ebv.R:172-192` with `inherits(x, "tidybreed_table")` plus a
`table_name` check. `define_additive_effects()` is the one writer still expressing a
population selection as an enum plus side-arguments.

### 1.3 What the table-typed argument buys beyond parity

Piping `ind_haplotype` directly gives a base that has **no expression today**:

```r
base_tbl = get_table(pop, "ind_haplotype") |> filter(line_origin == "Duroc")
```

centres on Duroc allele copies *wherever they sit*, including inside crossbreds.
For a line-scoped effect that is the population it actually applies to, and it
falls out of the design for free rather than needing a fourth argument.

---

## 2. Surface

### 2.1 `define_additive_effects()`

```r
define_additive_effects(tbl,                       # genome_meta, filtered: which loci
                        trait_name,
                        effects         = NULL,
                        distribution    = c("normal", "gamma"),
                        G               = NULL,
                        method          = c("shared", "union"),
                        base_tbl        = NULL,    # <- the only base argument
                        line_name       = NULL,
                        parent_origin   = NULL,
                        scale_to_target = TRUE,
                        seed            = NULL)
```

Removed: `base`, `base_line_name`. `base_tbl` keeps its name but changes meaning —
it is now honoured always and accepts three table kinds (§2.3).

```r
# Default: founder pool of the effect's own line
pop |> get_table("genome_meta") |> filter(chr %in% 1:5) |>
  define_additive_effects("ADG", line_name = "Duroc")
#  -> base_tbl = founder_haplotypes |> filter(line_name == "Duroc")

# Generation-0 individuals define p
pop |> get_table("genome_meta") |> filter(chr %in% 1:5) |>
  define_additive_effects("ADG",
    base_tbl = get_table(pop, "ind_meta") |> filter(gen == 0L))

# Force pooling for a line-specific effect (was: base_line_name = NULL)
pop |> get_table("genome_meta") |>
  define_additive_effects("ADG", line_name = "Duroc",
    base_tbl = get_table(pop, "founder_haplotypes"))

# Duroc allele copies wherever they sit, including inside crossbreds
pop |> get_table("genome_meta") |>
  define_additive_effects("ADG", line_name = "Duroc",
    base_tbl = get_table(pop, "ind_haplotype") |> filter(line_origin == "Duroc"))
```

### 2.2 `extract_allele_freq()` — new, exported

```r
extract_allele_freq(tbl)
#> # A tibble: n_loci × 3
#>   locus_id locus_name allele_freq
#>      <int> <chr>            <dbl>
```

**Contract.** Exactly one row per `genome_meta` locus, ascending `locus_id`;
`locus_id` integer, `locus_name` character, `allele_freq` double; **`NA_real_` at
each locus with no selected copies when the base is partly covered, an error when
no locus is covered at all**; every non-`NA` value finite and in `[0, 1]`
(checked — the haplotype tables should guarantee it, and the public helper fails
clearly if corruption breaks it); no database writes, no population mutation.

`allele_freq` is the frequency of allele 1 (the quantity used as Cockerham
`center_value`). `NA`, never `0`, for an absent locus: silently centring at 0 is
the failure mode the current zero-initialised vector invites. An observed
frequency of exactly `0` or `1` is a real value and is retained.

Usage with `ad_terms()`:

```r
p <- pop |> get_table("founder_haplotypes") |> filter(line_name == "A") |>
  extract_allele_freq()
loci <- c("Locus_10", "Locus_44")
tt <- ad_terms(loci, a = c(0.4, 0.1), d = c(0.2, 0.05),
               p = p$allele_freq[match(loci, p$locus_name)], coding = "cockerham")
pop |> define_genome_effects("ADG", tt, effect_owner = "custom")
```

### 2.3 Accepted `base_tbl` kinds

| `base_tbl$table_name` | Meaning | Required columns after any `select()` | `p` per locus |
|---|---|---|---|
| `founder_haplotypes` | the founder pool (filtered or not) | `locus_name`, `allele` | `AVG(allele)` over the filtered rows, joined to `genome_meta` on `locus_name` for `locus_id` |
| `ind_haplotype` | these allele copies (filtered or not) | `locus_id`, `allele` | `AVG(allele)` over the filtered rows directly |
| any other table with `id_ind` (`ind_meta`, `ind_phenotype`, `ind_tbv`, …) | these individuals | `id_ind` | `AVG(allele)` over `ind_haplotype` semi-joined to `SELECT DISTINCT id_ind FROM (<filtered>)` |
| anything else | error naming the table and listing the three accepted shapes | — | — |

The third row is the same "any table that has an `id_ind` column" convention the
action functions already document for candidate sets. **Selection semantics, not
row multiplicity:** for an `id_ind` table the frequency depends only on the
distinct selected individuals — `ind_phenotype` with five records per animal gives
the same `p` as `ind_meta` for the same ids. The `DISTINCT` is what guarantees
that and is tested (§3.8).

`line_name` is **not** required for `founder_haplotypes`: nothing in the helper
reads it from the rendered subquery. The default-path Wahlund check
(`.dae_default_base()`) and the empty-selection diagnostic
(`.founder_base_empty_error()`) both query the physical table, and the implicit
default always receives the unprojected table. So
`get_table(pop, "founder_haplotypes") |> select(locus_name, allele) |>
extract_allele_freq()` is a valid call. (v3 required it; v4 corrects that.)

**Two notions of "line".** `ind_meta.line_name` is a pedigree label — an F1 is
whatever the user called it (`"F1"`, `"DL"`); `ind_haplotype.line_origin` is the
founding line each allele copy traces to (`"Duroc"`, `"Landrace"`), propagated
through every `add_offspring()` call. They coincide for purebreds and diverge for
crosses. So:

- `ind_meta |> filter(line_name == "Duroc")` — *these animals*; for purebreds the
  same copies as the `line_origin` filter, for a crossbred label it pools every
  founding line present in those animals.
- `ind_haplotype |> filter(line_origin == "Duroc")` — *these copies*, wherever
  they sit; the correct base for a line-scoped effect in a crossbreeding program
  at any cross depth (F1, backcross, three-way, rotational).

**Warning boundary.** The Wahlund warning is raised only on the **default** path
in `define_additive_effects()` (a population-wide effect with `base_tbl = NULL`
on a multi-pool founder table). Any explicit `base_tbl` — including an explicit
whole `founder_haplotypes` — is an intentional selection and is not diagnosed for
pooling. `extract_allele_freq()` itself never warns: it is a computation, not a
modelling check.

### 2.4 `define_genome_effects()`

```r
define_genome_effects(pop, trait_name, terms,
                      effect_owner = "custom",
                      mode = c("append", "replace_scope", "replace_owner", "replace_trait"),
                      origin = NULL,
                      base_tbl = NULL,          # <- new
                      require_complete = FALSE,
                      allow_reserved_owner = FALSE)
```

When `base_tbl` is supplied, every `additive`/`dominance` member whose
`center_value` is `NA` receives `allele_freq` for its locus from
`extract_allele_freq(base_tbl)`. Rules:

- An explicit non-`NA` `center_value` is never overwritten.
- `indicator` members are never touched (a state has no centre); supplying
  `center_value` on one stays the error it is today.
- `base_tbl` is validated on every call it is supplied to (§3.1), but the
  frequency query runs **only if at least one additive/dominance input row has a
  missing `center_value`**. Supplying a base when every centre is explicit costs
  nothing.
- If a member needs a fill and its locus has `NA` frequency in the base, the
  existing row-specific "needs `center_value`" validation error fires — still
  naming the `term_id` and `locus_name` — extended with *"and `base_tbl` has no
  allele copies at this locus"*. It is not replaced by a generic list of loci.
- Without `base_tbl`, a missing centre is the error it is today. `NA` means
  "fill me" only in the presence of an explicit opt-in.
- Functional coding (`center_value = 0.5`) is not a fill target: `ad_terms()`
  writes `0.5` explicitly, and a user who wants functional coding writes it
  explicitly. The fill is Cockerham `p` only.

Two workflows, stated plainly so neither implies the other:

- Building terms through `ad_terms()`: get `p` from `extract_allele_freq()` and
  pass it. `ad_terms()` stays pure and rejects `NA` in `p` (`.ad_recycle()`).
- Hand-writing a terms data frame: leave `center_value` out (or `NA`) and supply
  `base_tbl` to the writer.

```r
# A hand-written dominance term with no p lookup step
pop |> define_genome_effects("ADG",
  data.frame(locus_name = "Locus_10", contrast_name = "dominance", genome_value = 0.8),
  base_tbl = get_table(pop, "founder_haplotypes") |> filter(line_name == "A"))
```

### 2.5 How the two writers relate

Both functions stay; they are **one engine, two verbs**, and `base_tbl` means the
same thing on each. This was considered and settled on 2026-09-20 (merge rejected:
mutually exclusive argument sets, a pipe subject whose type would change by mode,
RNG use that would depend on arguments, and a reserved-owner permission that
would become conditional — the same reasoning as the `add_offspring()` /
`add_doubled_haploids()` split in `plans/doubled_haploids.md`).

| | `define_genome_effects()` | `define_additive_effects()` |
|---|---|---|
| Altitude | **writer** — you supply coefficients | **generator** — samples coefficients to a target variance, then writes |
| Input | `terms` data frame | filtered `genome_meta` + trait spec |
| RNG | none | yes (`seed`) |
| `effect_owner` | any non-reserved | always `generated_additive_tbv` (what `add_tbv()` reads) |
| Replace | explicit `mode` | `replace_scope`, implicit |
| `base_tbl` | fills `NA` `center_value` | supplies `p` for centring **and** `scale_to_target` |
| Storage path | `.ge_build → .ge_read_model → .ge_resolve_deletes → .ge_commit` | **the same** (`R/define_additive_effects.R:333-337, 657-667`) |

`base_tbl = NULL` deliberately means different things: the generator has a
domain default (the effect's own founder pool, §3.3) because it knows what it is
defining; the writer has no basis for inventing a population, so `NULL` there
means "no fill" and a missing centre stays an error. Same argument, same
population semantics when supplied, different defaults — and the roxygen for
each says which.

The generator is sugar over the writer, and that must be **provable**, not just
asserted — see acceptance test 22. Both share `extract_allele_freq()` for `p`,
so a base selection means the same population in either call.

The roxygen for each function carries one sentence naming the family:
*"`define_genome_effects()` writes any effect you supply; `define_*_effects()`
functions sample effects of one shape and write them through the same path."*
This is what a future `define_dominance_effects()` inherits.

Two vocabulary mismatches remain between the verbs and are deliberately **not**
fixed here (§3.5): scope is `line_name` + `parent_origin` on the generator but
`origin = list(...)` on the writer; replacement is implicit on the generator but
`mode =` on the writer.

---

## 3. Implementation sketch

### 3.1 Shared validator — `.validate_base_tbl(base_tbl, pop)`

One internal helper, called by both writers before `extract_allele_freq()` (and by
`extract_allele_freq()` itself with `pop = base_tbl$pop`, so the public helper
gets the same column checks):

```r
.validate_base_tbl <- function(base_tbl, pop, arg = "base_tbl") {
  if (!inherits(base_tbl, "tidybreed_table"))
    stop("'", arg, "' must be a tidybreed_table from get_table() |> filter(...).")
  validate_tidybreed_pop(base_tbl$pop)
  if (!identical(base_tbl$pop$db_conn, pop$db_conn))
    stop("'", arg, "' must be piped from get_table() on the same pop as 'tbl'.")

  cols <- colnames(base_tbl$tbl)            # the *projected* columns, after select()
  need <- switch(base_tbl$table_name,
    founder_haplotypes = c("locus_name", "allele"),
    ind_haplotype      = c("locus_id", "allele"),
    "id_ind")
  miss <- setdiff(need, cols)
  if (length(miss))
    stop("'", arg, "' (", base_tbl$table_name, ") is missing column(s) ",
         toString(miss), " needed to compute allele frequencies. ",
         "Accepted shapes: founder_haplotypes, ind_haplotype, or any table with id_ind.")
  invisible(base_tbl)
}
```

The reason the column check exists: `select.tidybreed_table()` replaces
`.data$tbl` with the projected lazy query and returns the same wrapper, so a
`tidybreed_table` whose `table_name` is `"ind_haplotype"` may no longer expose
`allele`. Without this check the failure surfaces as a DuckDB binder error inside
generated SQL.

> **Implemented** as sketched, in `R/extract_allele_freq.R` (the lower layer,
> so the public helper can call it with `pop = NULL`). `.founder_base_empty_error()`
> lives in the same file for the same reason.

### 3.2 `extract_allele_freq()` (`R/extract_allele_freq.R`)

```r
extract_allele_freq <- function(tbl) {
  .validate_base_tbl(tbl, tbl$pop, arg = "tbl")
  conn <- tbl$pop$db_conn
  sub  <- dbplyr::sql_render(tbl$tbl)       # the user's filter/select, rendered to SQL

  freq_sql <- switch(tbl$table_name,
    founder_haplotypes = paste0(
      "SELECT gm.locus_id, AVG(CAST(b.allele AS DOUBLE)) AS allele_freq ",
      "FROM (", sub, ") b JOIN genome_meta gm ON b.locus_name = gm.locus_name ",
      "GROUP BY gm.locus_id"),
    ind_haplotype = paste0(
      "SELECT b.locus_id, AVG(CAST(b.allele AS DOUBLE)) AS allele_freq ",
      "FROM (", sub, ") b GROUP BY b.locus_id"),
    # individuals: semi-join on the distinct selected ids
    paste0(
      "SELECT h.locus_id, AVG(CAST(h.allele AS DOUBLE)) AS allele_freq ",
      "FROM ind_haplotype h ",
      "JOIN (SELECT DISTINCT id_ind FROM (", sub, ") b) ids USING (id_ind) ",
      "GROUP BY h.locus_id"))

  # LEFT JOIN to genome_meta is what guarantees one row per locus and lets
  # "no copies" (NA) be told apart from an observed frequency of 0.
  out <- DBI::dbGetQuery(conn, paste0(
    "SELECT gm.locus_id, gm.locus_name, f.allele_freq ",
    "FROM genome_meta gm LEFT JOIN (", freq_sql, ") f USING (locus_id) ",
    "ORDER BY gm.locus_id"))

  if (all(is.na(out$allele_freq))) {
    # This proves no usable locus matched -- not that the source had zero rows
    # (a founder pool whose locus_name keys do not match genome_meta looks the same).
    if (tbl$table_name == "founder_haplotypes") .founder_base_empty_error(conn, "this selection")
    stop("The filtered base contains no allele copies matching genome_meta loci.")
  }
  f <- out$allele_freq[!is.na(out$allele_freq)]
  if (any(!is.finite(f) | f < 0 | f > 1)) stop("allele frequencies outside [0, 1] -- corrupt haplotype rows?")

  tibble::as_tibble(out)
}
```

One founder-only diagnostic is carried over from today's
`compute_base_allele_freq()`: `.founder_base_empty_error(conn, what)` runs
`SELECT DISTINCT line_name FROM founder_haplotypes` (the physical table) and
raises the old *"No founder_haplotypes rows for <what>. Available: 'A', 'B', and
an unnamed (NULL) pool"* message. `what` is `"this selection"` from the helper
and `"line 'A'"` from the default resolver (§3.3); same function, one signature. v1 claimed a generic empty check preserved this; it does
not, so it is kept as an explicit branch.

The Wahlund warning is **not** here (v2 had it here). `extract_allele_freq()` is a
pure computation; whether pooling is a modelling mistake depends on what the
caller is defining, which only `define_additive_effects()` knows (§3.3).

So: one statement for the frequencies, at most one more on the empty-selection
error path.

Should `pop$tables` not contain `founder_haplotypes` at all, `get_table()` already
errors before any of this runs; the old "Did you call define_founder_haplotypes()?"
hint moves onto that path (or is dropped if `get_table()`'s message is judged
sufficient).

> **Implemented** as sketched. The "Did you call `define_founder_haplotypes()`?"
> hint was **kept**, on `.dae_default_base()`, and now suggests
> `base_tbl = get_table(pop, "ind_meta")` as the alternative — more useful than
> `get_table()`'s generic table-missing message.

### 3.3 `define_additive_effects()` rewiring

```r
  if (is.null(base_tbl)) {
    base_tbl <- .dae_default_base(pop, line_name)          # may warn (Wahlund)
  } else {
    .validate_base_tbl(base_tbl, pop)                       # explicit: never warns
  }
  p_base <- extract_allele_freq(base_tbl)$allele_freq     # locus_id order; may hold NA
```

The default resolver mirrors `resolve_genome_map()` / `resolve_chr_inheritance()`:
a line-specific pool wins, the shared (`NULL`) pool is the fallback, and only when
neither exists is it an error.

```r
.dae_default_base <- function(pop, line_name) {
  fh <- get_table(pop, "founder_haplotypes")   # errors if the table does not exist
  conn <- pop$db_conn
  if (is.null(line_name)) {
    # Population-wide effect: the whole founder table. Pooling several pools
    # overstates within-line heterozygosity (Wahlund) -- warn. NULL counts as
    # its own pool; COUNT(DISTINCT line_name) would ignore it.
    n_pools <- DBI::dbGetQuery(conn,
      "SELECT COUNT(DISTINCT COALESCE(line_name, '')) AS n FROM founder_haplotypes")$n
    if (n_pools > 1L) warning("founder_haplotypes holds ", n_pools, " pools; base allele ",
      "frequencies for this population-wide effect pool all of them ... ",
      "Pass base_tbl explicitly to silence this.")
    return(fh)
  }
  has_line <- DBI::dbGetQuery(conn,
    "SELECT EXISTS(SELECT 1 FROM founder_haplotypes WHERE line_name = ?) AS ok",
    params = list(line_name))$ok
  if (has_line)
    return(dplyr::filter(fh, .data$line_name == .env$line_name))   # .data/.env: no ambiguity
  has_shared <- DBI::dbGetQuery(conn,
    "SELECT EXISTS(SELECT 1 FROM founder_haplotypes WHERE line_name IS NULL) AS ok")$ok
  if (has_shared)
    return(dplyr::filter(fh, is.na(.data$line_name)))               # shared pool: line -> NULL
  .founder_base_empty_error(conn, paste0("line '", line_name, "'"))   # "... Available: ..."
}
```

A line-scoped effect that falls back to the shared pool does **not** warn: a
shared pool is one population by construction, so there is nothing to pool.

> **Implemented** with one placement change. The sketch above resolves the base
> once at the top of the function; in practice that made the Wahlund warning
> fire *before* the multi-trait argument checks (the mixed-`parent_origin`
> error in `gate 44` acquired a warning it never had). The resolution now runs
> **after argument validation, in each path**, at the point the old
> `compute_base_allele_freq()` call sat, through one helper
> `.dae_resolve_base(pop, base_tbl, line_name)` → `list(p_base, label)`. The
> completion message reports the base as `founder_haplotypes [1 filter]` /
> `ind_meta` rather than an enum value.

Back in the writer, after `p_base` is in hand and **before anything consumes
it**:

```r
  # Loud at the QTL that will actually be written -- per trait under method = "union"
  .dae_require_base_at(p_base, qtl_tf, tbl_loci)    # stop("base_tbl has no allele copies at QTL loci: Locus_12, Locus_40 ...")
```

Placement matters: `p_base` can now hold `NA`, and both
`rescale_effects_to_target()` call sites (`R/define_additive_effects.R:324`,
`:511`) fold it into `V_A = Σ n_eligible · p q a²`. An `NA` reaching either one
would make the scale factor `NA` and every effect `NA` — or, worse, be silently
dropped by an `na.rm` — so the check runs first, on the same per-trait
`qtl_tf` (or `qtl_tf_t[non_na]` in the union path, line 529) that the write will
use.

Everything downstream (`.dae_sample_effects()`, `scale_to_target`, the origin
composition, the writer call) is untouched — `p_base` is the same `n_loci`-length
numeric it is today. `compute_base_allele_freq()` is deleted.

**Parity claim, stated precisely.** For valid, complete founder pools the
line-aware default is numerically unchanged. Diagnostics for empty or malformed
selections change wording (and the "Did you call `define_founder_haplotypes()`?"
hint may move); this is behavioural parity on valid inputs, not identical failure
modes.

### 3.4 `define_genome_effects()` fill — inside `.ge_build()`

`.ge_build()` (`R/define_genome_effects.R:195`) resolves `locus_name → locus_id`,
assembles `members` (line 304), infers copy counts (320) and validates member
fields (321) in one pass with no seam. The fill goes **into** `.ge_build()` via a
new argument, the smaller of the two shapes the review offered:

```r
# define_genome_effects()
  if (!is.null(base_tbl)) .validate_base_tbl(base_tbl, pop)      # always
  # An omitted center_value column is the documented way to ask for a fill
  # (§2.4), and terms$center_value is then NULL -- any(NULL & ...) is FALSE, so the
  # naive test would skip the query and .ge_build() would raise the old error.
  # Normalise exactly as .ge_build() does (l.245-249) before deciding.
  centres <- if ("center_value" %in% names(terms)) terms$center_value
             else rep(NA_real_, nrow(terms))
  needs_fill <- !is.null(base_tbl) && "contrast_name" %in% names(terms) &&
    any(is.na(centres) & terms$contrast_name %in% c("additive", "dominance"))
  base_freq <- if (needs_fill) extract_allele_freq(base_tbl) else NULL   # lazily
  built <- .ge_build(conn, trait_name, terms, origin, effect_owner, base_freq = base_freq)

# .ge_build(), between members assembly (l.304) and .ge_check_member_fields() (l.321)
  if (!is.null(base_freq)) {
    fill <- is.na(members$center_value) &
            members$contrast_name %in% c("additive", "dominance")
    members$center_value[fill] <-
      base_freq$allele_freq[match(members$locus_id[fill], base_freq$locus_id)]
    members$fill_failed <- fill & is.na(members$center_value)   # for the message only
  }
  members <- .ge_infer_copy_counts(conn, members, labels)
  .ge_check_member_fields(members, labels)   # existing per-row error, now with the
                                             # "base_tbl has no copies at this locus" suffix
                                             # when fill_failed is TRUE
```

`.ge_check_member_fields()` keeps its existing per-row message (term_id from
`labels`, `locus_name`, contrast); only the suffix is new, and `fill_failed` is
dropped before anything is written.

> **Implemented** exactly as sketched (including the v4 `needs_fill`
> normalisation). The roxygen `@param base_tbl` also states the one-`p`-per-call
> rule from §4.5.

### 3.5 `plans/TODO.md` entries (not implemented here)

- `add_ebv()`: rename `phenotype` → `phenotype_tbl` for `*_tbl` consistency, and
  replace `collect()`-the-table-to-get-ids with the same rendered-subquery
  approach.
- `genome_meta.founder_allele_freq`: decide between dropping it, making it
  line-keyed, or documenting its last-call-wins meaning (see Q8).
- **Unify the two writers' vocabulary** (§2.5): one scope spelling — either the
  generator accepts `origin =` with `line_name`/`parent_origin` kept as sugar
  that composes into one origin row, or the writer's `origin` accepts the short
  form `list(line_name = , parent_origin = )` — and `mode =` exposed on the
  generator instead of a hard-wired `replace_scope`. Separate change; it touches
  argument surfaces this plan does not.

### 3.6 Call-site inventory

A repository-wide search — `grep -rn 'base = "\|base_tbl\|base_line_name\|base = c('
R/ tests/testthat/ vignettes/` — finds **73** lines (v1 said 26; that count used a
narrower pattern). Executable `base = "current_pop"` calls, which are the ones
that would **silently change meaning** if the argument were merely removed:

| File | `current_pop` call sites | Migration |
|---|---|---|
| `vignettes/swine/swine-time-based-age-at-puberty-sex-semen.R` | 6 executable (lines 917, 1157, 1239, 1312, 1440, 1509) + 1 comment (313) | each becomes an explicit `base_tbl = get_table(pop, "ind_meta")` (whole current population), or the line-filtered `ind_meta` where the comment says so |
| `tests/testthat/test-add_tbv.R` | 3 | explicit `base_tbl = get_table(pop, "ind_meta") \|> filter(...)` |
| `tests/testthat/test-define_additive_effects.R` | 2 | same |
| `R/define_additive_effects.R` | 3 (roxygen) | rewritten examples |
| `R/define_founder_haplotypes.R` | 1 (roxygen) | rewritten cross-reference |

The remaining references are `base = "founder_haplotypes"` (drop — it is the
default), `base_line_name = "A"` (drop when `line_name = "A"` is also set, else
`base_tbl = founder_haplotypes |> filter(line_name == "A")`), explicit
`base_line_name = NULL` (→ `base_tbl = get_table(pop, "founder_haplotypes")`),
and `base = "current_pop", base_tbl = x` (→ `base_tbl = x`). Roxygen in
`R/open_pop.R` and `R/restore_pop.R` mentions the old arguments in prose.

**Rule:** the implementer regenerates this inventory with the grep above at the
time of the change and migrates every executable `current_pop` call to an
explicit whole-population `base_tbl`. Never just delete the argument.

> **Implemented.** The as-executed inventory is in the implementation summary.
> Beyond the table above, twelve `suppressWarnings()` wrappers around
> pooled-default calls in `test-add_tbv.R` and `test-genome-effects-eval.R`
> became an explicit `base_tbl = get_table(pop, "founder_haplotypes")` — the
> Q11 "tell" is gone — and the ten `base_line_name = "A"` in
> `test-genome-effects-writer.R` became `base_tbl = gew_base_A(pop)`. The grep
> returns nothing after commit 2.

Per the pre-1.0 policy the old names are removed outright — no aliases, no
deprecation messages, no tests that exercise them.

### 3.7 Order of work

1. `.validate_base_tbl()` + `extract_allele_freq()` with focused SQL/contract
   tests.
2. Rewire `define_additive_effects()`; delete `base`, `base_line_name`,
   `compute_base_allele_freq()`; rewrite its tests; migrate every
   `current_pop` call site from §3.6.
3. Add `base_freq` to `.ge_build()`, expose `base_tbl` on
   `define_genome_effects()`; tests.
4. Only after both writer paths pass the full suite: roxygen (including the
   one-sentence family statement from §2.5 on both functions), `NAMESPACE`,
   `man/`, `CLAUDE.md` (the `define_additive_effects()` section and its
   `base = "current_pop"` example), `NEWS.md`, `DESCRIPTION` version.
5. Log the §3.5 items in `plans/TODO.md`.

### 3.8 Acceptance tests

`test-extract_allele_freq.R` (new):

1. Each of the three kinds against a hand-computed frequency.
2. `ind_meta`, repeated-record `ind_phenotype`, and direct `ind_haplotype`
   selections **agree** when they select the same allele copies — individuals are
   not weighted by record count.
3. `founder_haplotypes |> filter(line_name == "A")` agrees with
   `ind_haplotype |> filter(line_origin == "A")` on founders-only data.
4. A base with some but not all loci returns `NA` only at absent loci; an
   observed `0` or `1` is retained and is not `NA`.
5. `select()` removing a required column → package error naming the table and
   the column, not a DuckDB binder error.
6. A table without `id_ind` → error listing the three accepted shapes.
7. Empty filtered founder base → the "Available: …" diagnostic; empty
   non-founder base → the generic message.
8. `extract_allele_freq()` never warns — not on a pooled `founder_haplotypes`,
   not on any other shape.
9. Output is `locus_id`-ordered, one row per locus, typed as in §2.2; the
   `_schema_meta` row counts and every table are unchanged afterwards
   (read-only).

`test-define_additive_effects.R`:

10. Existing base-population tests rewritten to the new surface.
11. `base_tbl` from another population/connection errors.
12. A selected QTL locus with no base copies errors naming the loci; under
    `method = "union"` the check is per trait, on the loci actually written.
13. Same `seed` reproduces itself (within-code reproducibility; **not** a golden
    comparison against pre-change output).
14. **Default resolution** (§3.3): (a) line pool exists → used; (b) only a
    shared pool → used, silently; (c) both exist → line pool wins; (d) neither →
    error listing available pools, naming the `NULL` pool as such.
15. **Wahlund warning boundary**: fires for a population-wide effect with
    `base_tbl = NULL` on a two-pool table; counts a `NULL` pool + one named pool
    as two; does **not** fire for the same effect with an explicit
    `base_tbl = get_table(pop, "founder_haplotypes")`; does not fire for a
    line-scoped effect on any path.
16. **Crossbreeding end to end** (the §4.1 script): two default-base line calls
    plus one explicit pooled common call produce three variants with distinct `center_value`s, and `add_tbv()` on
    Duroc, Landrace, F1 and backcross animals matches a hand-computed per-copy
    TBV. This is the existing test at `test-add_tbv.R:220` rewritten without
    `base =` / `base_line_name`.
17. On a fixed fixture containing F1s: `ind_haplotype |> filter(line_origin ==
    "Duroc")` equals a hand-computed `AVG(allele)` over exactly those realized
    Duroc-origin copies (**not** the Duroc founder-pool frequency — realized
    copies are a finite transmitted sample and equal the pool only in
    expectation); `ind_meta |> filter(line_name == "F1")` equals the
    hand-computed frequency over all copies of those F1 animals and includes
    both origins; and, with the fixture chosen so they differ, the two are not
    equal.

`test-genome-effects-writer.R`:

18. `base_tbl` fills `center_value` on additive/dominance members both when the
    column is present with `NA` and when the column is **omitted entirely**;
    leaves explicit values alone; never touches indicators.
19. No frequency query when every centre is explicit (assert via a mocked
    `extract_allele_freq()` or a `DBI` call counter).
    > **Implemented** without a mock or counter: the test passes a base whose
    > *query would fail* (no copies at any locus) with every centre explicit —
    > if it were queried the call would error, and it does not; the same base
    > with a missing centre does error. Sharper, and no tracing of `DBI`.
20. A failed fill keeps the `term_id` + `locus_name` in the error and appends the
    "no copies at this locus" suffix.
21. `base_tbl` from another population/connection errors here too.
22. **Generator ≡ writer (§2.5).** With a fixed seed, capture the
    `genome_value`s `define_additive_effects()` sampled (from
    `genome_effect_loci`), snapshot the three effect tables and `ind_tbv`, then
    **on the same population** (so the base is identical by construction) write
    those values with
    `define_genome_effects(terms = <locus_name, contrast_name = "additive",
    genome_value>, origin = <the composed scope>,
    effect_owner = "generated_additive_tbv", mode = "replace_scope",
    allow_reserved_owner = TRUE, base_tbl = <the same base>)` — with
    **`center_value` omitted**, so the writer must fill it from `base_tbl`.
    Copying the stored centres would leave `base_tbl` validated but unqueried
    and prove only storage-path parity; omitting them proves the shared base
    semantics too. `mode = "replace_scope"` replaces the generator's variant
    with the writer's; the three effect tables must then equal the snapshot
    (surrogate ids aside) and a re-run `add_tbv()` must equal the snapshot
    exactly. Run once for the common
    scope and once for a `line_name` + `parent_origin` scope. If this cannot be
    made to pass, the two verbs have drifted and that is a bug in this change,
    not in the test.

`vignettes/swine/…`: runs end to end after migration.
> **Implemented:** run through the last `define_additive_effects()` call (line
> 1520) on the full 10k-locus / 2000-founder config with `load_all()`; all six
> migrated calls and the WWD/WWM correlated call ran with `base: ind_meta`,
> exit 0. The breeding loop beyond that point needs the installed package and
> BLUPF90 and does not touch the changed code.

---

## 4. Multi-line and crossbreeding walkthrough

The design's job in a multi-line program is to give each **variant** of an
additive term (common / line A / line B …) the centring that matches the copies
it will be applied to. `add_tbv()` picks, per allele copy, the most specific
variant whose origin matches the copy's `(line_origin, parent_origin)` label and
falls back to the common variant; each variant carries its own `center_value`.
So the question for every configuration is: *does the default, or an easy
explicit `base_tbl`, hand each variant the right `p`?*

### 4.1 Two purebred lines and their F1 — the canonical case

```r
pop <- open_pop(...) |> define_genome(...) |>
  define_founder_haplotypes(n_haplotypes = 200, line_name = "Duroc",    ...) |>
  define_founder_haplotypes(n_haplotypes = 200, line_name = "Landrace", ...)

pop <- pop |> get_table("founder_haplotypes") |> filter(line_name == "Duroc") |>
  add_founders(n_males = 20, n_females = 200, line_name = "Duroc")
pop <- pop |> get_table("founder_haplotypes") |> filter(line_name == "Landrace") |>
  add_founders(n_males = 20, n_females = 200, line_name = "Landrace")

pop <- define_trait(pop, "ADG", target_add_var = 0.25)
gm  <- pop |> get_table("genome_meta") |> filter(chr %in% 1:5)

# Three variants, three calls. The common one names its base to say "yes, pool";
# the line-scoped ones take the default:
pop <- gm |> define_additive_effects("ADG",
         base_tbl = get_table(pop, "founder_haplotypes"))   # common fallback: pooled *on purpose*, no warning
pop <- gm |> define_additive_effects("ADG", line_name = "Duroc")      # default -> Duroc pool
pop <- gm |> define_additive_effects("ADG", line_name = "Landrace")   # default -> Landrace pool

# F1s and backcrosses: every copy is centred on its own founding line's p
pop <- add_offspring(pop, crosses)          # line_origin travels with each copy
pop <- pop |> get_table("ind_meta") |> add_tbv("ADG")
```

Today the same three calls need `base = "current_pop", base_tbl = duroc_tbl`
style arguments (`tests/testthat/test-add_tbv.R:236-241`); after the change the
line-scoped ones need nothing. The common call is written with an explicit whole
`founder_haplotypes` to say "yes, pool" — with `base_tbl = NULL` it would produce
the same numbers plus a Wahlund warning, which is correct: an *unintended*
population-wide effect on a two-line base is exactly what the warning exists for.

### 4.2 A shared founder pool feeding several lines

```r
pop <- pop |> define_founder_haplotypes(n_haplotypes = 400)          # line_name = NULL
pop <- pop |> get_table("founder_haplotypes") |> add_founders(..., line_name = "A")
pop <- pop |> get_table("founder_haplotypes") |> add_founders(..., line_name = "B")
```

One base population, lines diverge later by selection and drift. A line-scoped
effect here should centre on the shared pool — there is no other founder `p` —
and the v3 default does so via the `line → NULL` fallback (§3.3). v2's default
(and today's code) errors with "no rows for line 'A'". If, generations later, the
user wants line-specific centring reflecting drift, the explicit forms are one
line each:

```r
base_tbl = get_table(pop, "ind_haplotype") |> filter(line_origin == "A")   # A's copies, any depth
base_tbl = get_table(pop, "ind_meta") |> filter(line_name == "A", gen == 0L) # A's founders
```

### 4.3 Crossbred base populations

`ind_meta |> filter(line_name == "F1")` is a valid base — it is "these animals" —
but its `p` is a 50:50 mixture of Duroc and Landrace copies. That is the right
base for something like a terminal-cross commercial population evaluated as its
own entity, and the wrong one for a line-scoped variant. The plan does not try
to detect this (§2.3 warning boundary): `line_origin` is the tool for line-scoped
centring, and `ind_meta.line_name` for animal-set centring, and the roxygen
states the distinction with this example.

### 4.4 Three-way, rotational and backcross programs

Nothing new: `line_origin` on every copy is the founding line however many
generations of crossing sit between. `filter(line_origin == "Duroc")` is exact at
any depth, and the evaluator's per-copy fallback already handles copies whose
founding line has no variant of its own. Line-specific pools are optional —
§4.2's shared-pool default covers programs that start from one base.

### 4.5 Line-scoped terms in `define_genome_effects()`

One `base_tbl` per call gives one `p` per locus per call. That is the same
one-scope-per-call shape `define_additive_effects()` has, and the writer already
supports `mode = "append"`, so line-scoped surfaces are written one line at a
time:

```r
# copy_count = 2: a dominance member takes an exact multiset summing to the
# state's copy count; a scalar list without it is refused for genotype members
pop |> define_genome_effects("ADG", dom_terms,
  origin   = list(line_match_type = "exact", line_name = "Duroc", copy_count = 2L),
  base_tbl = get_table(pop, "founder_haplotypes") |> filter(line_name == "Duroc"))
pop |> define_genome_effects("ADG", dom_terms,
  origin   = list(line_match_type = "exact", line_name = "Landrace", copy_count = 2L),
  base_tbl = get_table(pop, "founder_haplotypes") |> filter(line_name == "Landrace"))
```

A reciprocal member (one term whose genotype member names two lines, as in the
`CLAUDE.md` example) has exactly one `center_value` by storage design; which `p`
belongs there is a modelling choice the user makes, and the fill supplies
whichever base was passed. This is a pre-existing property of the storage model,
not something the plan changes, and the roxygen for `base_tbl` says so.

### 4.6 Things checked and found unaffected

- **Sex chromosomes / hemizygous copies.** `AVG(allele)` over however many rows a
  locus has (2, 1 or 0 per individual) is the copy frequency in every case; the
  same arithmetic as today.
- **`method = "union"` per-trait QTL sets and the multi-trait `G` path.** They
  consume `p_base` exactly as before; the per-trait "no base copies at a written
  QTL" check (§3.3) is the only touch point.
- **`parent_origin` (imprinting) variants.** Orthogonal to the base; the line
  default composes with `parent_origin` exactly as `base_line_name` did.
- **Founders not yet added.** The founder-pool default works before any
  `add_founders()` call, as today; an `ind_meta`/`ind_haplotype` base naturally
  requires individuals and errors with the empty-selection message otherwise.
- **`define_founder_haplotypes()` roxygen** (lines 19–23) currently tells
  multi-line users to reach for `base = "current_pop"`; it is rewritten to point
  at the default's per-line centring and the `line_origin` form.

---

## 5. Open questions — for review

Each with a recommendation. Codex agreed with Q1–Q9 in its first review and with
Q10–Q11 in its re-review of v3; the verdicts below record that and the
qualifications attached. **All eleven are settled** unless a reviewer reopens one.

### Q1. Should `base_tbl` be required rather than default to the effect's own line?

**Options.** (a) `base_tbl = NULL` → founder pool of `line_name`, pooled when
`line_name = NULL`. (b) Always required.

**Recommendation: (a). Codex: agree.** The default is not a convenience shortcut;
it encodes the Wahlund argument already documented in the roxygen ("centre on the
population the effect applies to"). Making it required would force every
single-line simulation to write `base_tbl = get_table(pop, "founder_haplotypes")`
for no information gain, and would make the wrong thing (pooling across lines for
a line-scoped effect) exactly as easy to type as the right thing. Keep the
default; keep the Wahlund warning for the pooled case.

### Q2. Is the name `base_tbl` right?

**Options.** `base_tbl` (keeps the existing name, matches the `tbl` first-argument
convention). `base` alone (shorter, but reads like an enum — which is what we are
removing). `base_pop_tbl` (more descriptive).

**Recommendation: `base_tbl`. Codex: agree.** The `_tbl` suffix is the package's
signal that a `tidybreed_table` is expected; `add_ebv()` broke that with
`phenotype =` and this plan should not add a second exception (§3.5 logs the
`add_ebv()` fix).

### Q3. Accept three table kinds, or pin to one like `add_ebv()` pins `phenotype` to `ind_phenotype`?

**Options.** (a) Three kinds dispatched on `table_name` + required columns (§2.3).
(b) Only tables with `id_ind` — founders would be selected as `ind_meta |>
filter(gen == 0L)`. (c) Only `founder_haplotypes` and `ind_haplotype`, the two
tables that hold alleles.

**Recommendation: (a). Codex: agree, provided the required-column and
same-connection checks are part of the contract** — they now are (§3.1). (b) loses
the founder *pool* as a base — the pool is larger than any sampled founder set and
is the theoretically correct base for `scale_to_target`; it also cannot express
"Duroc copies inside crossbreds". (c) loses the "these individuals" selection,
which is how every other action function takes a candidate set.

### Q4. Should a QTL locus with no copies in the base be an error or a warning + `NA` skip?

**Options.** (a) Error naming the loci. (b) Warn and drop those loci from the
QTL set. (c) Warn and centre at `p = 0.5`.

**Recommendation: (a). Codex: agree, with the check limited to loci actually
written per trait** (matters under `method = "union"`; §3.3). The current code's
zero-initialised vector silently centres at `p = 0` — exactly the bug this plan
removes. (b) changes the QTL set behind the user's back, which violates "no
implicit selection". (c) invents a frequency.

### Q5. Should `define_genome_effects()` get `base_tbl` at all, or is `extract_allele_freq()` enough?

**Options.** (a) Both: the helper for `ad_terms()`, plus the `NA`-fill in the
writer. (b) Helper only; `center_value` stays fully user-supplied.

**Recommendation: (a). Codex: agree, with the `.ge_build()` integration (§3.4),
lazy querying, explicit-wins, and never-fill-indicators.** The concern with (a)
is that "`NA` means fill me" is a second meaning for `NA` in a column that today
means "missing, error". Mitigation, now explicit in §2.4: the fill runs only when
`base_tbl` is supplied; without it, `NA` is the error it is today. (b) remains a
legitimate stopping point if reviewers still find the dual meaning uncomfortable.

### Q6. Should `ad_terms()` accept a `base_tbl` too, so `p` can be omitted?

**Options.** (a) No — `ad_terms()` is a pure builder with no database access;
users call `extract_allele_freq()` and pass `p`. (b) Yes — `ad_terms(..., p =
NULL, base_tbl = )`.

**Recommendation: (a). Codex: agree.** `ad_terms()` reports the implied mean
`μ = a(p − q) + 2pq·d`, is testable without a population, and rejects `NA` in `p`
via `.ad_recycle()` — none of that should change. The two workflows are stated in
§2.4: `ad_terms()` users fetch `p` with the helper; users who want the writer to
fill centres hand-write the terms data frame instead.

### Q7. Should `add_ebv()`'s `phenotype` argument be renamed / re-implemented for consistency?

**Options.** (a) Out of scope; leave it. (b) Rename to `phenotype_tbl` and switch
its `collect()`-the-whole-table id extraction to the same subquery approach.

**Recommendation: (a) for this plan, (b) as a follow-up. Codex: agree.** It is a
real inconsistency and a real inefficiency, but it has nothing to do with genome
effects and the BLUPF90 paths have no CI coverage. Logged in §3.5.

### Q8. Should `genome_meta.founder_allele_freq` be dropped now that `extract_allele_freq()` exists?

**Options.** (a) Keep for now, decide separately. (b) Drop in this change.

**Recommendation: (a), for scope reasons only. Codex: agree on the outcome,
disagree with v1's rationale — accepted.** v1 called the column "a snapshot of
the pool as sampled". It is not: `define_founder_haplotypes()` documents that in a
multi-line setup it "describes only the pool written last"
(`R/define_founder_haplotypes.R:19`). It is last-call-wins, cannot say which line
it describes, and is read by no writer. The honest reasons to leave it alone
*here* are that nothing depends on it and that its fate (drop / make line-keyed /
document the limitation) deserves its own small decision rather than riding on
this plan. Logged in §3.5.

### Q9. What if `base_tbl` and `tbl` come from different `tidybreed_pop` objects that point at the same file?

**Options.** (a) Compare `db_conn` identity — strict, errors on a `restore_pop()`
copy. (b) Compare the database path — permissive.

**Recommendation: (a), applied and tested in both writers. Codex: agree.** Two
connections to the same DuckDB file are not a supported composition boundary
anywhere else in the package, and accepting them complicates connection-lifetime
reasoning. (v1 said rendered SQL would be "simply wrong" on a second connection;
that overstates it — the SQL text itself is connection-independent — but the
strict provenance rule is still the safer public contract.) The error message
says "pipe `get_table()` from the same `pop`".

### Q10. Should the default base for a line-scoped effect fall back to the shared (`NULL`) pool?

**Options.** (a) Yes — `line = L → line = NULL`, mirroring `resolve_genome_map()`
and `resolve_chr_inheritance()`; error only when neither pool exists. (b) No —
a line-scoped effect requires a same-named pool; a shared-pool program passes
`base_tbl = get_table(pop, "founder_haplotypes")` explicitly on every line call.

**Recommendation: (a). Codex: agree.** `founder_haplotypes.line_name = NULL` already *means*
"the pool every line samples from" (`add_founders()` reads it that way), and the
package has one precedence rule for nullable line dimensions everywhere else.
Making the base the single place that rule does not apply would be a surprise,
and (b) errors on a correct, common design (§4.2) with a message whose only
possible fix is boilerplate. The fallback adds at most two `EXISTS` queries
(`has_line`, then `has_shared`) and cannot
mis-centre: if a same-named pool exists it always wins, and the shared pool is one
population by construction so nothing is pooled by falling back to it.

Codex's qualification, adopted: this is **pool-level fallback, not per-locus
stitching**. Once a named pool exists it is the base, complete or not; a named pool
missing some loci is not silently completed from the shared pool, and the
QTL-level "no base copies at a written locus" error (§3.3) is the correct failure
for that malformed case.

### Q11. Where does the Wahlund warning live, and when does it fire?

**Options.** (a) Inside `extract_allele_freq()` for any multi-pool
`founder_haplotypes` selection (v2). (b) Only on the `base_tbl = NULL` path of
`define_additive_effects()` for a population-wide effect; never for an explicit
`base_tbl`; never in the helper.

**Recommendation: (b). Codex: agree.** The warning exists to catch a user who *did not notice*
they were centring a population-wide effect on a pooled base. In every
crossbreeding program the common fallback variant is defined by deliberately
pooling (§4.1), and under (a) that call warns forever — the test suite already
wraps it in `suppressWarnings()` three times (`test-add_tbv.R:130,161,203`),
which is what a mis-placed warning looks like. An explicit `base_tbl` is the
user's answer to the warning's question. (b) also keeps `extract_allele_freq()`
a pure computation, which is what an `extract_` function should be. Count pools
as `COUNT(DISTINCT COALESCE(line_name, ''))` so the shared pool counts as one.

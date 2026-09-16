#------------------------------------------------------------------------------#
# tidybreed — package summary generator
#------------------------------------------------------------------------------#
#
# Recomputes every package metric from the working tree and writes two files:
#
#   package_summary.md                       (repo root; the Markdown version)
#   dev/package_summary/package_summary.html (styled page, from template.html)
#
# Run from the package root:
#
#   Rscript dev/package_summary/render_package_summary.R
#
# The HTML page is meant to be published as a Claude artifact (or opened in a
# browser / printed to PDF). Nothing here is loaded by the package; `dev/` is
# in .Rbuildignore.
#
# The "Database tables" section is read from a live in-memory population
# (open_pop() |> define_genome() |> ... |> define_trait() materializes every
# table and view), so column counts, types and descriptions come from the DDL
# the package actually runs. That needs pkgload::load_all() on the working
# tree; everything else is base R.
#
# Counting method (kept stable so snapshots are comparable across versions):
#   - "function definitions" = lines in R/*.R matching
#       ^name <- function  or  ^name = function   (top level, dotted names count)
#   - "lines" = raw line counts (wc -l), comments and blanks included
#   - test_that() / expect_*() are counted as occurrences in tests/testthat/test-*.R
#------------------------------------------------------------------------------#

root <- normalizePath(".")
if (!file.exists(file.path(root, "DESCRIPTION")))
  stop("Run from the package root (DESCRIPTION not found in ", root, ")")

out_dir  <- file.path(root, "dev", "package_summary")
template <- file.path(out_dir, "template.html")
out_html <- file.path(out_dir, "package_summary.html")
out_md   <- file.path(root, "package_summary.md")

#------------------------------------------------------------------------------#
# helpers
#------------------------------------------------------------------------------#

n_lines <- function(files) sum(vapply(files, function(f) length(readLines(f, warn = FALSE)), 1L))
count_pattern <- function(files, pattern) {
  sum(vapply(files, function(f) sum(grepl(pattern, readLines(f, warn = FALSE), perl = TRUE)), 1L))
}
count_matches <- function(files, pattern) {
  sum(vapply(files, function(f) {
    m <- gregexpr(pattern, readLines(f, warn = FALSE), perl = TRUE)
    sum(vapply(m, function(x) sum(x > 0), 1L))
  }, 1L))
}
fmt  <- function(x) format(x, big.mark = ",", trim = TRUE)
k    <- function(x) paste0(format(round(x / 1000, 1), nsmall = 1), "k")
git  <- function(...) trimws(system2("git", c(...), stdout = TRUE))
desc_field <- function(field) {
  d <- read.dcf(file.path(root, "DESCRIPTION"))
  if (!field %in% colnames(d)) return(character(0))
  v <- trimws(strsplit(d[1, field], ",")[[1]])
  v[nzchar(v)]
}
html_esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

#------------------------------------------------------------------------------#
# metrics
#------------------------------------------------------------------------------#

r_files    <- list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)
cpp_files  <- list.files(file.path(root, "src"), "\\.cpp$", full.names = TRUE)
man_files  <- list.files(file.path(root, "man"), "\\.Rd$", full.names = TRUE)
test_files <- list.files(file.path(root, "tests", "testthat"), "^test-.*\\.R$", full.names = TRUE)
help_files <- list.files(file.path(root, "tests", "testthat"), "^helper-.*\\.R$", full.names = TRUE)
vign_files <- list.files(file.path(root, "vignettes"), "\\.(Rmd|qmd)$", full.names = TRUE)
plan_files <- list.files(file.path(root, "plans"), "\\.md$", full.names = TRUE)

ns <- readLines(file.path(root, "NAMESPACE"))
exports    <- sort(sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE)))
s3_methods <- sum(grepl("^S3method\\(", ns))

m <- list(
  version   = unname(read.dcf(file.path(root, "DESCRIPTION"))[1, "Version"]),
  date      = format(Sys.Date()),
  branch    = git("branch", "--show-current"),
  sha       = git("rev-parse", "--short", "HEAD"),
  license   = unname(read.dcf(file.path(root, "DESCRIPTION"))[1, "License"]),
  n_exports = length(exports),
  n_s3      = s3_methods,
  n_fns     = count_pattern(r_files, "^[.A-Za-z_][A-Za-z0-9._]* *(<-|=) *function"),
  n_r_files = length(r_files),
  r_loc     = n_lines(r_files),
  roxy_loc  = count_pattern(r_files, "^#'"),
  n_cpp     = length(cpp_files),
  cpp_loc   = n_lines(cpp_files),
  n_man     = length(man_files),
  n_tests   = length(test_files),
  test_loc  = n_lines(test_files),
  n_help    = length(help_files),
  help_loc  = n_lines(help_files),
  n_test_that = count_matches(test_files, "test_that\\("),
  n_expect  = count_matches(test_files, "expect_[a-z_]*\\("),
  n_vign    = length(vign_files),
  vign_desc = paste(sprintf("`%s`, %s lines", basename(vign_files),
                            fmt(vapply(vign_files, n_lines, 1L))), collapse = "; "),
  news_loc  = n_lines(file.path(root, "NEWS.md")),
  news_ver  = count_pattern(file.path(root, "NEWS.md"), "^# tidybreed"),
  readme_loc = n_lines(file.path(root, "README.md")),
  n_plans   = length(plan_files),
  plan_loc  = n_lines(plan_files),
  commits   = as.integer(git("rev-list", "--count", "HEAD")),
  imports   = desc_field("Imports"),
  linkingto = desc_field("LinkingTo"),
  suggests  = desc_field("Suggests"),
  sysreq    = desc_field("SystemRequirements")
)
m$ratio <- round(m$test_loc / m$r_loc, 2)
# DuckDB stack first, then the rest in DESCRIPTION order
core_imports    <- c("duckdb", "DBI", "dbplyr")
imports_ordered <- c(core_imports[core_imports %in% m$imports], setdiff(m$imports, core_imports))
m$license <- sub(" \\+ file LICENSE$", "", m$license)

# ---- database tables --------------------------------------------------------
# Materialize every table/view in memory and read the registries the package
# itself uses (schema(), describe_table(), TABLE_ROW_KEYS, archive_replicate()).
old_opt <- options(tidybreed.quiet = TRUE)     # silences the load banner too
suppressMessages(pkgload::load_all(root, quiet = TRUE))
pop <- suppressMessages(
  open_pop(pop_name = "summary", db_name = ":memory:") |>
    define_genome(n_loci = 20, n_chr = 1, chr_len_Mb = 10) |>
    define_founder_haplotypes(n_haplotypes = 4) |>
    get_table("founder_haplotypes") |>
    add_founders(n_males = 1, n_females = 1, line_name = "A") |>
    define_trait("T")
)
sch <- as.data.frame(schema(pop, show_empty = TRUE, include_system = TRUE))
kinds <- DBI::dbGetQuery(pop$db_conn,
  "SELECT table_name, table_type FROM information_schema.tables")
cols <- DBI::dbGetQuery(pop$db_conn,
  "SELECT table_name, column_name, data_type FROM information_schema.columns
   ORDER BY table_name, ordinal_position")
cons <- DBI::dbGetQuery(pop$db_conn,
  "SELECT table_name, constraint_type, constraint_column_names AS cols
     FROM duckdb_constraints()
    WHERE constraint_type IN ('PRIMARY KEY', 'UNIQUE', 'FOREIGN KEY')")
n_col_described <- DBI::dbGetQuery(pop$db_conn,
  "SELECT count(*) AS n FROM _schema_meta WHERE object_type = 'column'")$n
close_pop(pop)

# Key shown per table: the SQL PRIMARY KEY when one is declared, otherwise the
# logical row key the package enforces in R (TABLE_ROW_KEYS) -- marked as such.
key_of <- function(tbl) {
  pk <- cons$cols[cons$table_name == tbl & cons$constraint_type == "PRIMARY KEY"]
  if (length(pk)) return(list(cols = pk[[1]], logical = FALSE))
  k <- TABLE_ROW_KEYS[[tbl]]
  if (is.null(k)) list(cols = character(0), logical = FALSE) else list(cols = k, logical = TRUE)
}
options(old_opt)

# which function creates each table: first CREATE TABLE/VIEW site in R/, in
# pipeline order of the files that own DDL
ddl_files <- c("open_pop.R", "define_genome.R", "define_founder_haplotypes.R",
               "define_trait.R", "define_effect_cov_matrix.R", "genome_effects_helpers.R")
# tables written with dbWriteTable() rather than CREATE TABLE
creator_override <- c(founder_haplotypes = "define_founder_haplotypes()")
creator_of <- function(tbl) {
  if (tbl %in% names(creator_override)) return(creator_override[[tbl]])
  pat <- sprintf('CREATE (TABLE|VIEW)( IF NOT EXISTS)? "?%s"?\\b', tbl)
  for (f in ddl_files) {
    if (any(grepl(pat, readLines(file.path(root, "R", f), warn = FALSE), perl = TRUE))) {
      if (f == "genome_effects_helpers.R")           # views built by helpers
        return(if (grepl("^genome_effect_", tbl)) "define_genome()" else "define_trait()")
      return(paste0(sub("\\.R$", "", f), "()"))
    }
  }
  ""
}

# archive_replicate() treatment, read from its formals
arch <- formals(archive_replicate)
archive_of <- function(tbl) {
  if (tbl %in% eval(arch$store_and_reset)) "per replicate"
  else if (tbl %in% eval(arch$store_once)) "once"
  else if (tbl %in% eval(arch$reset_only)) "reset only"
  else "kept"
}

tables <- data.frame(
  table_name  = sch$table_name,
  table_group = as.character(sch$table_group),
  kind        = ifelse(kinds$table_type[match(sch$table_name, kinds$table_name)] == "VIEW", "view", "table"),
  n_cols      = sch$n_cols,
  row_key     = vapply(sch$table_name, function(t) paste(key_of(t)$cols, collapse = ", "), ""),
  key_logical = vapply(sch$table_name, function(t) key_of(t)$logical, TRUE),
  created_by  = vapply(sch$table_name, creator_of, ""),
  archive     = vapply(sch$table_name, archive_of, ""),
  description = sch$description,
  stringsAsFactors = FALSE
)
tables$table_group[tables$table_name == "_schema_meta"] <- "System"
type_tally <- sort(table(cols$data_type), decreasing = TRUE)
m$n_tables  <- sum(tables$kind == "table")
m$n_views   <- sum(tables$kind == "view")
m$n_columns <- sum(tables$n_cols)
m$n_groups  <- length(unique(tables$table_group))
m$n_described <- n_col_described
m$n_pk      <- sum(!tables$key_logical & nzchar(tables$row_key))
m$n_logical <- sum(tables$key_logical)

# ---- largest files -----------------------------------------------------------
# Known purposes; anything else falls back to the file's first roxygen title.
purposes <- c(
  "schema.R"                  = "Table registry, descriptions, `schema()` / `describe_table()`",
  "add_phenotype.R"           = "Phenotype simulation (composite, SGE, fixed/random effects, residuals)",
  "define_genome_effects.R"   = "General genome-effect writer: terms, members, origins, replace modes",
  "define_additive_effects.R" = "QTL effect sampling, Falconer rescale, multi-trait MVN, line/parent-origin scope",
  "genome_effects_eval.R"     = "The one evaluator behind `add_tbv()` / `add_tgv()`",
  "add_offspring.R"           = "Mating, gamete formation, offspring haplotype writes",
  "mutate_table.R"            = "Generic typed column add/update on any table",
  "formula_helpers.R"         = "Formula-based derived phenotypes",
  "add_ebv.R"                 = "EBV import and BLUPF90 wrapper",
  "add_tbv.R"                 = "True breeding values from the reserved additive terms",
  "define_founder_haplotypes.R" = "Founder haplotype pools (no-LD and LD methods)",
  "genome_effects_helpers.R"  = "Validation and canonical forms for genome-effect rows"
)
first_roxygen_title <- function(f) {
  l <- readLines(f, warn = FALSE)
  t <- grep("^#' \\S", l, value = TRUE)
  if (length(t)) sub("^#' ", "", t[[1]]) else ""
}
top_files <- function(files, n = 6) {
  loc <- vapply(files, n_lines, 1L)
  o <- order(-loc)[seq_len(min(n, length(files)))]
  data.frame(file = basename(files[o]), loc = unname(loc[o]), stringsAsFactors = FALSE)
}
src_top  <- top_files(r_files)
src_top$purpose <- ifelse(src_top$file %in% names(purposes),
                          purposes[src_top$file],
                          vapply(file.path(root, "R", src_top$file), first_roxygen_title, ""))
test_top <- top_files(test_files)

# ---- API by prefix ---------------------------------------------------------
api_group <- function(fn) {
  if (grepl("^(open|restore|close)_", fn)) "session"
  else if (grepl("^define_", fn))          "define_"
  else if (grepl("^add_", fn))             "add_"
  else if (grepl("^mutate_", fn))          "mutate_"
  else if (grepl("^(extract|remove|archive)_", fn)) "extract_"
  else if (grepl("_terms$", fn))           "terms"
  else                                     "inspection"
}
api_meta <- list(
  session    = c(label = "`open_` / `restore_` / `close_`",       html = "open_ / restore_ / close_",     verb = "session"),
  define_    = c(label = "`define_`",                             html = "define_",                       verb = "write model configuration"),
  add_       = c(label = "`add_`",                                html = "add_",                          verb = "insert simulation output"),
  mutate_    = c(label = "`mutate_`",                             html = "mutate_",                       verb = "add or update columns"),
  extract_   = c(label = "`extract_` / `remove_` / `archive_`",  html = "extract_ / remove_ / archive_", verb = "read out, delete, stamp replicates"),
  terms      = c(label = "Term builders",                         html = "term builders",                 verb = "rows for define_genome_effects()"),
  inspection = c(label = "Inspection",                            html = "inspection",                    verb = "lazy dplyr access")
)
api <- split(exports, vapply(exports, api_group, ""))
api <- api[names(api_meta)[names(api_meta) %in% names(api)]]
# Within a group, list functions in workflow order (the order a simulation
# calls them), not alphabetically; anything new lands at the end, sorted.
pipeline_order <- c(
  "open_pop", "restore_pop", "close_pop",
  "define_genome", "define_chromosome", "define_founder_haplotypes", "define_chip",
  "define_trait", "define_trait_simple", "define_additive_effects", "define_genome_effects",
  "define_phenotype", "define_residual_cov", "define_effect_cov_matrix", "define_effect_random",
  "define_effect_fixed_class", "define_effect_fixed_cov", "define_effect_intercept",
  "define_index", "define_table", "define_schema_description",
  "add_founders", "add_offspring", "add_phenotype", "add_tbv", "add_tgv",
  "add_ebv", "add_index", "add_dosage", "add_genotypes",
  "mutate_table", "mutate_derived", "mutate_group_seq", "mutate_group_named", "mutate_group_concatenate",
  "extract_genotypes", "remove_rows", "archive_replicate",
  "ad_terms", "genotype_terms",
  "get_table", "schema", "describe_table"
)
api <- lapply(api, function(v) c(pipeline_order[pipeline_order %in% v], sort(setdiff(v, pipeline_order))))

#------------------------------------------------------------------------------#
# Markdown
#------------------------------------------------------------------------------#

md <- c(
  "# tidybreed — Package Summary", "",
  sprintf("**Version:** %s · **Snapshot date:** %s · **Branch:** `%s` (`%s`)",
          m$version, m$date, m$branch, m$sha), "",
  "A database-first (DuckDB) breeding-program simulator in R. All genomic and",
  "individual data lives in a file-based DuckDB database rather than R memory,",
  "enabling simulations larger than RAM, resumable runs, and replicate archiving.",
  "A small Rcpp kernel handles meiosis/recombination.", "",
  "## Highlights", "",
  "| | |", "|---|---|",
  sprintf("| Exported functions | %d (+%d S3 methods; %d functions total incl. internals) |", m$n_exports, m$n_s3, m$n_fns),
  sprintf("| R source | %d files, ~%s lines (~%s of which are roxygen docs) |", m$n_r_files, fmt(round(m$r_loc, -2)), fmt(round(m$roxy_loc, -2))),
  sprintf("| C++ (Rcpp) | %d files, ~%s lines (gamete/recombination kernel) |", m$n_cpp, fmt(round(m$cpp_loc, -1))),
  sprintf("| Documentation | %d man pages, %d vignette%s, %s-line README |", m$n_man, m$n_vign, if (m$n_vign == 1) "" else "s", fmt(round(m$readme_loc, -1))),
  sprintf("| Tests | %d testthat files, ~%s lines, %d tests, ~%s assertions |", m$n_tests, fmt(round(m$test_loc, -2)), m$n_test_that, fmt(round(m$n_expect, -1))),
  sprintf("| Test : source ratio | %.2f : 1 |", m$ratio),
  sprintf("| History | %d commits, %d released versions in NEWS.md |", m$commits, m$news_ver),
  sprintf("| Dependencies | %s |", paste(imports_ordered[imports_ordered != "stats"], collapse = ", ")),
  "",
  "## Detailed Counts", "",
  "| Metric | Count |", "|---|---:|",
  sprintf("| Exported functions (`NAMESPACE`) | %d |", m$n_exports),
  sprintf("| S3 methods registered | %d |", m$n_s3),
  sprintf("| Total R function definitions | %d |", m$n_fns),
  sprintf("| R source files (`R/`) | %d |", m$n_r_files),
  sprintf("| R lines of code | %s |", fmt(m$r_loc)),
  sprintf("| Roxygen doc lines (`#'`) in `R/` | %s |", fmt(m$roxy_loc)),
  sprintf("| C++ source files (`src/`) | %d |", m$n_cpp),
  sprintf("| C++ lines of code | %s |", fmt(m$cpp_loc)),
  sprintf("| Man pages (`man/*.Rd`) | %d |", m$n_man),
  sprintf("| testthat test files | %d |", m$n_tests),
  sprintf("| testthat lines of code | %s (+ %d helper files, %s lines) |", fmt(m$test_loc), m$n_help, fmt(m$help_loc)),
  sprintf("| `test_that()` blocks | %d |", m$n_test_that),
  sprintf("| `expect_*()` assertions | %s |", fmt(m$n_expect)),
  sprintf("| Vignettes | %d (%s) |", m$n_vign, m$vign_desc),
  sprintf("| `NEWS.md` | %s lines, %d version headings |", fmt(m$news_loc), m$news_ver),
  sprintf("| `README.md` | %s lines |", fmt(m$readme_loc)),
  sprintf("| Design docs (`plans/`) | %d files, %s lines |", m$n_plans, fmt(m$plan_loc)),
  sprintf("| Git commits | %d |", m$commits),
  "",
  sprintf("## Exported API (%d functions)", m$n_exports), "",
  "| Prefix | Functions |", "|---|---|",
  vapply(names(api), function(g) sprintf("| %s | %s |", api_meta[[g]][["label"]],
                                         paste0("`", api[[g]], "`", collapse = ", ")), ""),
  "",
  "## Largest Source Files", "",
  "| File | LOC | Purpose |", "|---|---:|---|",
  sprintf("| `%s` | %s | %s |", src_top$file, fmt(src_top$loc), src_top$purpose),
  "",
  "## Largest Test Files", "",
  "| Test file | LOC |", "|---|---:|",
  sprintf("| `%s` | %s |", test_top$file, fmt(test_top$loc)),
  "",
  "## Database Tables", "",
  sprintf("%d tables and %d views in %d groups, %d columns in total, %d of them described in `_schema_meta` (`schema()`, `describe_table()`).",
          m$n_tables, m$n_views, m$n_groups, m$n_columns, m$n_described), "",
  sprintf("Column types: %s.", paste(sprintf("%s ×%d", names(type_tally), type_tally), collapse = ", ")), "",
  sprintf("Keys: %d tables declare a SQL `PRIMARY KEY`; %d use a logical key enforced in R (`TABLE_ROW_KEYS`), by design where a DuckDB constraint would block bulk inserts or transactional replacement. `Archive` is how `archive_replicate()` treats the table: copied and stamped *per replicate*, copied *once*, *reset only*, or *kept* in the working database.",
          m$n_pk, m$n_logical), "",
  unlist(lapply(unique(tables$table_group), function(g) {
    tg <- tables[tables$table_group == g, ]
    c(sprintf("### %s (%d)", g, nrow(tg)), "",
      "| Table | Kind | Cols | Key | Created by | Archive | Description |",
      "|---|---|---:|---|---|---|---|",
      sprintf("| `%s` | %s | %d | %s | %s | %s | %s |",
              tg$table_name, tg$kind, tg$n_cols,
              ifelse(nzchar(tg$row_key),
                     paste0("`", gsub(", ", "`, `", tg$row_key), "`", ifelse(tg$key_logical, " (logical)", "")),
                     "—"),
              ifelse(nzchar(tg$created_by), paste0("`", tg$created_by, "`"), "—"),
              tg$archive, gsub("\\|", "\\\\|", tg$description)),
      "")
  })),
  "## Dependencies", "",
  sprintf("- **Imports:** %s", paste(m$imports, collapse = ", ")),
  sprintf("- **LinkingTo:** %s", paste(m$linkingto, collapse = ", ")),
  sprintf("- **Suggests:** %s", paste(gsub(">=", "≥", m$suggests), collapse = ", ")),
  if (length(m$sysreq)) sprintf("- **SystemRequirements:** %s", paste(m$sysreq, collapse = ", ")),
  sprintf("- **License:** %s", m$license)
)
writeLines(md, out_md)
message("Wrote ", out_md)

#------------------------------------------------------------------------------#
# HTML (fill template.html)
#------------------------------------------------------------------------------#

md_code_to_html <- function(x) html_esc(gsub("`([^`]+)`", "\\1", x)) |>
  (\(s) gsub("([^]+)", "<code>\\1</code>", s))()

stat <- function(n, label, sub, cls = "") sprintf(
  '      <div class="stat%s">\n        <div class="n">%s</div>\n        <div class="l">%s</div>\n        <div class="sub">%s</div>\n      </div>',
  if (nzchar(cls)) paste0(" ", cls) else "", n, label, sub)

ledger <- paste(c(
  stat(m$n_exports, "Exported functions", sprintf("+%d S3 methods · %d total", m$n_s3, m$n_fns), "accent"),
  stat(sub("k$", "<small>k</small>", k(m$r_loc)), "Lines of R", sprintf("%d files · %s roxygen", m$n_r_files, k(m$roxy_loc))),
  stat(fmt(m$cpp_loc), "Lines of C++", sprintf("%d files · gamete kernel", m$n_cpp)),
  stat(m$n_man, "Man pages", sprintf("%d vignette%s · %s-line README", m$n_vign, if (m$n_vign == 1) "" else "s", fmt(m$readme_loc))),
  stat(m$n_test_that, "Tests", sprintf("%d files · %s assertions", m$n_tests, fmt(m$n_expect)), "gold"),
  stat(sprintf("%.2f<small> : 1</small>", m$ratio), "Test : source ratio", sprintf("%s test lines / %s source", k(m$test_loc), k(m$r_loc)), "gold"),
  stat(m$commits, "Commits", sprintf("%d released versions", m$news_ver)),
  stat(length(m$imports), "Imports", paste(c(intersect(c("duckdb", "dplyr", "Rcpp", "dqrng"), m$imports), "…"), collapse = " · "))
), collapse = "\n")

row <- function(label, n, group = FALSE, note = NULL) sprintf(
  '          <tr%s><td>%s%s</td><td class="num">%s</td></tr>',
  if (group) ' class="group"' else "", md_code_to_html(label),
  if (is.null(note)) "" else sprintf(' <span class="note">(%s)</span>', md_code_to_html(note)), n)

counts_rows <- paste(c(
  row("Exported functions (`NAMESPACE`)", m$n_exports, TRUE),
  row("S3 methods registered", m$n_s3),
  row("Total R function definitions", m$n_fns),
  row("R source files (`R/`)", m$n_r_files, TRUE),
  row("R lines of code", fmt(m$r_loc)),
  row("Roxygen doc lines (`#'`) in `R/`", fmt(m$roxy_loc)),
  row("C++ source files (`src/`)", m$n_cpp, TRUE),
  row("C++ lines of code", fmt(m$cpp_loc)),
  row("Man pages (`man/*.Rd`)", m$n_man, TRUE),
  row("Vignettes", m$n_vign, note = m$vign_desc),
  row("`README.md` lines", fmt(m$readme_loc)),
  row("`NEWS.md` lines", fmt(m$news_loc), note = sprintf("%d version headings", m$news_ver)),
  row("Design docs (`plans/`)", fmt(m$plan_loc), note = sprintf("%d files", m$n_plans)),
  row("testthat test files", m$n_tests, TRUE),
  row("testthat lines of code", fmt(m$test_loc), note = sprintf("+ %d helper files, %s lines", m$n_help, fmt(m$help_loc))),
  row("`test_that()` blocks", m$n_test_that),
  row("`expect_*()` assertions", fmt(m$n_expect)),
  row("Git commits", m$commits, TRUE)
), collapse = "\n")

api_rows <- paste(vapply(names(api), function(g) sprintf(
  '      <div class="api-row">\n        <div class="prefix">%s<span class="verb">%s</span></div>\n        <div class="fns">%s</div>\n      </div>',
  api_meta[[g]][["html"]], html_esc(api_meta[[g]][["verb"]]),
  paste0("<code>", api[[g]], "</code>", collapse = "")), ""), collapse = "\n")

bar_rows <- function(df, cls = "", purpose = NULL) {
  w <- 100 * df$loc / max(df$loc)
  paste(sprintf(
    '      <div class="bar-row%s">\n        <div class="name">%s%s</div>\n        <div class="track"><div class="fill" style="width:%.1f%%"></div></div>\n        <div class="loc">%s</div>\n      </div>',
    if (nzchar(cls)) paste0(" ", cls) else "", df$file,
    if (is.null(purpose)) "" else sprintf('<span class="purpose">%s</span>', md_code_to_html(purpose)),
    w, fmt(df$loc)), collapse = "\n")
}

dep <- function(key, codes, core = character(0)) sprintf(
  '      <div class="dep">\n        <div class="k">%s</div>\n        <div class="v">%s</div>\n      </div>',
  key, paste0(sprintf('<code%s>%s</code>', ifelse(codes %in% core, ' class="core"', ""), html_esc(codes)), collapse = ""))
dep_rows <- paste(c(
  dep("Imports", imports_ordered, core = core_imports),
  dep("LinkingTo", m$linkingto),
  dep("Suggests", gsub(">=", "≥", m$suggests)),
  if (length(m$sysreq)) dep("System", m$sysreq),
  dep("License", m$license)
), collapse = "\n")

tables_ledger <- paste(c(
  stat(m$n_tables, "Tables", sprintf("in %d groups", m$n_groups), "accent"),
  stat(m$n_views, "Views", "derived, never stored"),
  stat(m$n_columns, "Columns", sprintf("%d described in _schema_meta", m$n_described)),
  stat(m$n_pk, "SQL primary keys", sprintf("%d logical keys enforced in R", m$n_logical))
), collapse = "\n")

tables_note <- paste(
  "Key is the declared SQL <code>PRIMARY KEY</code>, or the logical row key enforced in R",
  "where a DuckDB constraint would block bulk inserts or transactional replacement.",
  "Archive is how <code>archive_replicate()</code> treats the table: copied and stamped",
  "<b>per replicate</b>, copied <b>once</b>, <b>reset only</b>, or <b>kept</b> in the working database.")

tables_types <- paste(sprintf('<span class="type"><code>%s</code><b>%d</b></span>',
                              names(type_tally), type_tally), collapse = "\n")

tables_groups <- paste(unlist(lapply(unique(tables$table_group), function(g) {
  tg <- tables[tables$table_group == g, ]
  rows <- sprintf(
    '          <tr>\n            <td><span class="tname">%s</span><span class="tdesc">%s</span></td>\n            <td class="kind %s">%s</td>\n            <td class="num">%d</td>\n            <td class="key">%s</td>\n            <td class="key">%s</td>\n            <td class="arch">%s</td>\n          </tr>',
    tg$table_name, html_esc(tg$description), tg$kind, tg$kind, tg$n_cols,
    ifelse(nzchar(tg$row_key),
           paste0("<code>", gsub(", ", "</code> <code>", tg$row_key), "</code>",
                  ifelse(tg$key_logical, ' <span class="lk">logical</span>', "")),
           "—"),
    ifelse(nzchar(tg$created_by), paste0("<code>", tg$created_by, "</code>"), "—"),
    html_esc(tg$archive))
  c(sprintf('    <div class="tgroup">\n      <h3>%s <span class="count">%d %s</span></h3>\n      <div class="tablewrap">\n        <table class="tables">\n          <thead><tr><th>Table</th><th>Kind</th><th class="num">Cols</th><th>Key</th><th>Created by</th><th>Archive</th></tr></thead>\n          <tbody>',
            html_esc(g), nrow(tg), if (nrow(tg) == 1) "table" else "tables"),
    rows,
    '          </tbody>\n        </table>\n      </div>\n    </div>')
})), collapse = "\n")

fill <- list(
  version = m$version, date = m$date, branch = m$branch, sha = m$sha, license = m$license,
  n_exports = m$n_exports, ledger = ledger, counts_rows = counts_rows, api_rows = api_rows,
  src_bars = bar_rows(src_top, purpose = src_top$purpose),
  src_max_name = src_top$file[1], src_max_loc = fmt(src_top$loc[1]),
  test_bars = bar_rows(test_top, cls = "test"),
  test_max_name = test_top$file[1], test_max_loc = fmt(test_top$loc[1]),
  dep_rows = dep_rows,
  tables_ledger = tables_ledger, tables_types = tables_types, tables_groups = tables_groups,
  tables_note = tables_note
)
html <- paste(readLines(template, warn = FALSE), collapse = "\n")
for (key in names(fill)) html <- gsub(paste0("{{", key, "}}"), fill[[key]], html, fixed = TRUE)
left <- regmatches(html, gregexpr("\\{\\{[a-z_]+\\}\\}", html))[[1]]
if (length(left)) stop("Unfilled placeholders: ", paste(left, collapse = ", "))
writeLines(html, out_html)
message("Wrote ", out_html)

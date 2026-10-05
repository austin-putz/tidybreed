# SQL text recording, for the hard rule that individual ids never appear in
# SQL text (CLAUDE.md). DuckDB's dbSendQuery method carries every statement,
# dbplyr's included.

# Every SQL statement DuckDB receives while `code` runs.
record_sql <- function(code) {
  rec <- new.env()
  rec$sql <- character()
  # do.call() passes the built expression: trace() quotes its `tracer`.
  suppressMessages(do.call(trace, list("dbSendQuery",
    signature = c("duckdb_connection", "character"),
    tracer = bquote(assign("sql", c(get("sql", envir = .(rec)), statement),
                           envir = .(rec))),
    where = asNamespace("DBI"), print = FALSE), quote = TRUE))
  on.exit(suppressMessages(untrace("dbSendQuery",
    signature = c("duckdb_connection", "character"),
    where = asNamespace("DBI"))), add = TRUE)
  force(code)
  rec$sql
}

# The ids that appear quoted ('id') in any of the statements.
leaked_ids <- function(sql, ids) {
  quoted <- paste0("'", ids, "'")
  ids[vapply(quoted, function(q) any(grepl(q, sql, fixed = TRUE)), logical(1))]
}

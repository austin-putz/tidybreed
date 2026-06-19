library(devtools)
results <- test(pkg = ".", reporter = "summary")
if (!all(results$failed == 0)) quit(status = 1)

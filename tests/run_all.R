# Run from any directory: Rscript /path/to/project/tests/run_all.R
argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(argument) != 1L) stop("Run this file with Rscript.")
test_dir <- dirname(normalizePath(sub("^--file=", "", argument)))
for (script in c("validate.R", "validate_preparations.R", "validate_project.R")) {
  source(file.path(test_dir, script), local = new.env(parent = globalenv()))
}
cat("\nAll validation suites passed.\n")

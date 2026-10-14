# Integration checks for project paths, evaluation, summaries, and cached data.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
source(file.path(project_dir, "R", "load_all.R"))

passed <- 0L
check <- function(label, expr) {
  force(expr)
  passed <<- passed + 1L
  cat("PASS:", label, "\n")
}
expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}

check("all active R files parse", {
  directories <- file.path(project_dir, c("R", "scripts", "simulations", "applications", "tests"))
  files <- unlist(lapply(directories, list.files, pattern = "\\.R$", full.names = TRUE))
  for (path in files) parse(path)
})

check("entry scripts only define functions when sourced from another directory", {
  previous <- getwd()
  setwd(tempdir())
  env <- new.env(parent = globalenv())
  entries <- c("simulations/01_grid_x_pattern.R", "simulations/02_knn_graph.R",
                 "simulations/03_sbm_graph.R", "simulations/05_oracle_vs_cv.R",
                 "scripts/grid_signal.R", "scripts/plot_summary_curves.R",
                 "applications/chicago_comparison.R", "applications/india_preparation.R")
  for (entry in entries) source(file.path(project_dir, entry), local = env)
  stopifnot(identical(env$gfm_project_root(), project_dir),
            exists("run_chicago_comparison", envir = env, inherits = FALSE))
  setwd(previous)
})

check("raw proportions preserve zero-total rows rather than inventing values", {
  Y <- rbind(c(2, 1), c(0, 0))
  P <- multinomial_mle(Y)
  stopifnot(all(P[1, ] == c(2 / 3, 1 / 3)), all(is.na(P[2, ])))
})

check("scores agree with individual outcomes and correctly exclude empty nodes", {
  Y <- matrix(c(3, 1, 2, 4, 0, 0), 3, 2, byrow = TRUE,
                dimnames = list(c("a", "b", "c"), c("one", "two")))
  P <- matrix(c(0.8, 0.2, 0.1, 0.9, 0.5, 0.5), 3, 2, byrow = TRUE, dimnames = dimnames(Y))
  scored <- score_multinomial(Y, P)
  expected_mse <- mean((matrix(c(0.75, 0.25, 1 / 3, 2 / 3), 2, 2, byrow = TRUE) - P[1:2, ])^2)
  expected_loss <- -(3 * log(0.8) + log(0.2) + 2 * log(0.1) + 4 * log(0.9)) / 10
  stopifnot(scored$tested_nodes == 2, scored$test_population == 10,
            abs(scored$mse - expected_mse) < 1e-12,
            abs(scored$log_loss - expected_loss) < 1e-12,
            abs(scored$brier - 0.484) < 1e-12,
            identical(scored, score_multinomial(Y, P[3:1, 2:1])))
  P[1, ] <- c(0, 1)
  stopifnot(is.infinite(score_multinomial(Y, P)$log_loss))
})

check("MCAR fitted-count conversion and isolated-node rejection", {
  Y <- matrix(c(5, 5, 3, 7), 2, 2, byrow = TRUE,
                dimnames = list(c("a", "b"), c("one", "two")))
  fake <- list(fitted.values = matrix(c(8, 2, 1, 9), 2, 2, byrow = TRUE))
  P <- gfm_mcar_probabilities(fake, Y)
  stopifnot(all(P == matrix(c(0.8, 0.2, 0.1, 0.9), 2, 2, byrow = TRUE)),
            identical(dimnames(P), dimnames(Y)))
  expect_error(fit_mcar_multinomial(Y, matrix(0, 2, 2)))
})

check("summary means, axes, and stored uncertainty are unchanged", {
  path <- file.path(project_dir, "results", "mse_grid_n.csv")
  table <- utils::read.csv(path)
  original <- table
  long <- gfm_summary_long(table, "n")
  stopifnot(identical(table, original), identical(long$x[1:5], table$n),
            identical(long$value[1:5], table$gfm), identical(long$sd[1:5], table$gfm_sig))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    plot <- plot_summary_curves(setNames(rep(list(table), 3), c("Grid", "kNN", "SBM")))
    stopifnot(inherits(plot, "patchwork"))
  }
})

check("supplied Chicago cache, CSVs, graph, and holdout are aligned", {
  cache <- file.path(project_dir, "data", "chicago_crime_data")
  chicago <- readRDS(file.path(cache, "chicago_2024_grid1000m_min1_all.rds"))
  counts <- utils::read.csv(file.path(cache, "chicago_2024_grid1000m_min1_all_counts.csv"),
                            check.names = FALSE)
  stopifnot(identical(dim(chicago$Y), c(652L, 6L)), sum(chicago$Y) == 256625,
            nrow(chicago$D_G) == 1176, all(Matrix::rowSums(chicago$W) > 0),
            length(unique(chicago$component)) == 2,
            all(chicago$Y_train + chicago$Y_test == chicago$Y),
            sum(chicago$Y_train) == 205051, sum(chicago$Y_test) == 51574,
            all(as.matrix(counts[, colnames(chicago$Y)]) == chicago$Y))
  connected <- chicago_largest_component(chicago)
  repeated <- chicago_largest_component(connected)
  stopifnot(nrow(connected$Y) == 650, connected$analysis_selection$excluded_records == 5,
            identical(connected$Y_train,
                      chicago$Y_train[rownames(connected$Y), , drop = FALSE]),
            identical(connected$analysis_selection, repeated$analysis_selection),
            is.null(connected$output_file))
})

cat("\nAll", passed, "project validation groups passed.\n")

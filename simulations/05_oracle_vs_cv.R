# Load the project without changing the working directory.
local({
  target_env <- parent.env(environment())
  source_files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1L))
  source_files <- source_files[nzchar(source_files)]
  if (length(source_files) > 0L) {
    entry <- tail(source_files, 1L)
  } else {
    cli_file <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(cli_file) != 1L) stop("Source this script or run it with Rscript.")
    entry <- sub("^--file=", "", cli_file)
  }
  root <- dirname(dirname(normalizePath(entry, mustWork = TRUE)))
  source(file.path(root, "scripts", "bootstrap.R"), local = target_env)
})

#' One dataset per target size: truth-selected versus CV-selected lambda.
#' This illustration is not a replicated Monte Carlo validation of CV.
run_oracle_vs_cv <- function(n_targets = c(100, 500, 1000, 2000), D = 4L, seed = 1L,
                              lambda_grid = c(0.5, 1, 2.5, 4, 6), K = 5L,
                              T_max = 4000L, tol = 1e-6,
                              project_root = gfm_project_root()) {
  if (!is.numeric(n_targets) || !length(n_targets) || any(!is.finite(n_targets)) ||
      any(n_targets < 4 | n_targets != floor(n_targets))) {
    stop("n_targets must contain integer sizes of at least four.", call. = FALSE)
  }
  rows <- .gfm_with_seed(seed, lapply(n_targets, function(target) {
    side <- round(sqrt(target))
    g <- build_grid_graph(side, side)
    D_G <- build_incidence_matrix(g)
    truth <- generate_true_beta_X(g, D, side, side)
    Y <- simulate_data(truth$P_true, sample(150:300, igraph::vcount(g), replace = TRUE))
    mse <- numeric(length(lambda_grid))
    warm <- NULL
    for (index in order(lambda_grid)) {
      fit <- fit_graph_fused_multinomial(Y, D_G, lambda_grid[index], T_max = T_max,
                                         tol = tol, beta_init = warm, verbose = FALSE)
      warm <- fit$beta_hat
      mse[index] <- mean((fit$P_hat - truth$P_true)^2)
    }
    best <- which.min(mse)
    cv <- cv_graph_fused_multinomial(Y, D_G, lambda_grid, K = K, seed = 42L,
                                     T_max = T_max, tol = tol, verbose = FALSE)
    fitted <- fit_graph_fused_multinomial(Y, D_G, cv$lambda_min, T_max = T_max,
                                          tol = tol, verbose = FALSE)
    mse_cv <- mean((fitted$P_hat - truth$P_true)^2)
    data.frame(n_target = target, grid_dim = paste0(side, "x", side),
                 n_actual = igraph::vcount(g), n_edges = nrow(D_G),
                 lambda_oracle = lambda_grid[best], mse_oracle = mse[best],
                 lambda_cv = cv$lambda_min, mse_cv = mse_cv,
                 mse_diff = mse_cv - mse[best],
                 pct_increase = if (mse[best] > 0) 100 * (mse_cv - mse[best]) / mse[best] else NA_real_,
                 cv_fit_converged = fitted$converged)
  }))
  table <- do.call(rbind, rows)
  result <- list(metrics = table, settings = list(n_targets = n_targets, seed = seed, D = D, K = K,
                                                  lambda_grid = lambda_grid, cv_seed = 42L,
                                                  T_max = T_max, tol = tol))
  paths <- gfm_save_workflow(result, "oracle_vs_cv", project_root)
  utils::write.csv(table, file.path(paths$results, "oracle_vs_cv_table.csv"), row.names = FALSE)
  if (requireNamespace("xtable", quietly = TRUE)) {
    latex <- xtable::xtable(table[, c("n_actual", "n_edges", "lambda_oracle", "mse_oracle",
                                      "lambda_cv", "mse_cv", "mse_diff", "pct_increase")],
                            caption = "Truth-selected and cross-validated lambda on the grid.",
                            label = "tab:oracle_vs_cv", digits = c(0, 0, 0, 2, 6, 2, 6, 6, 2))
    text <- capture.output(print(latex, include.rownames = FALSE))
    writeLines(text, file.path(paths$results, "oracle_vs_cv_table.tex"))
  }
  print(table, row.names = FALSE)
  invisible(result)
}

if (sys.nframe() == 0L) {
  flags <- gfm_cli_flags(allowed = "--quick")
  if (flags$quick) run_oracle_vs_cv(c(25, 49), K = 2, lambda_grid = c(0.5, 2.5),
                                   T_max = 4000, tol = 1e-4)
  else run_oracle_vs_cv()
}


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

#' Generate the representative category-by-method grid signal figure.
run_grid_signal <- function(side = 20L, D = 4L, seed = 1L, lambda = 2.5,
                             include_mcar = TRUE, T_max = 5000L, tol = 1e-6,
                             burnin = 5000L, n_sample = 20000L, thin = 10L,
                             project_root = gfm_project_root()) {
  gfm_require_packages("ggplot2")
  .gfm_check_flag(include_mcar, "include_mcar")
  g <- build_grid_graph(side, side)
  truth <- generate_true_beta_X(g, D, side, side)
  Y <- .gfm_with_seed(seed, simulate_data(truth$P_true,
                                         sample(150:300, igraph::vcount(g), replace = TRUE)))
  category_names <- paste("Category", seq_len(D))
  colnames(Y) <- colnames(truth$P_true) <- category_names
  D_G <- build_incidence_matrix(g)
  mle_time <- system.time(P_mle <- multinomial_mle(Y))
  gfm_time <- system.time(fit <- fit_graph_fused_multinomial(
    Y, D_G, lambda, T_max = T_max, tol = tol, verbose = FALSE
  ))
  probabilities <- list(True = truth$P_true, MLE = P_mle, GFM = fit$P_hat)
  timings <- c(MLE = unname(mle_time["elapsed"]), GFM = unname(gfm_time["elapsed"]))
  mcar <- NULL
  if (include_mcar) {
    mcar <- fit_mcar_multinomial(Y, igraph::as_adjacency_matrix(g, sparse = TRUE),
                                burnin, n_sample, thin, seed = seed + 1L)
    probabilities$MCAR <- mcar$P_hat
    timings["MCAR"] <- mcar$elapsed_seconds
  }
  estimates <- probabilities[-1L]
  metrics <- data.frame(method = names(estimates),
                         mse = vapply(estimates, function(P) mean((P - truth$P_true)^2), numeric(1L)),
                         relative_error = vapply(estimates, relative_probability_error,
                                                   P_true = truth$P_true, numeric(1L)),
                         elapsed_seconds = unname(timings[names(estimates)]), row.names = NULL)
  labels <- c("True", sprintf("%s\n(RE = %.2f%%)", metrics$method, 100 * metrics$relative_error))
  plot <- plot_probability_surfaces(cbind(igraph::V(g)$col, igraph::V(g)$row),
                                    probabilities, labels = labels)
  result <- list(Y = Y, truth = truth, g = g, fit = fit, mcar = mcar,
                   probabilities = probabilities, metrics = metrics,
                   settings = list(side = side, D = D, seed = seed, lambda = lambda,
                                     include_mcar = include_mcar, T_max = T_max, tol = tol,
                                     burnin = burnin, n_sample = n_sample, thin = thin,
                                     mcar_seed = seed + 1L))
  paths <- gfm_save_workflow(result, "grid_signal", project_root)
  gfm_save_plot(plot, file.path(paths$figures, "signal.pdf"), 7.2, 6.4)
  print(metrics)
  invisible(result)
}

if (sys.nframe() == 0L) {
  flags <- gfm_cli_flags()
  if (flags$quick) run_grid_signal(side = 5, include_mcar = flags$include_mcar,
                                   T_max = 2000, tol = 1e-5)
  else run_grid_signal(include_mcar = flags$include_mcar)
}


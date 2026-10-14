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

#' Illustrate X-pattern recovery on rook and king grids.
#' CV selects lambda on the rook graph; the same value is used on the king graph.
run_grid_demo <- function(nrow_grid = 50L, ncol_grid = 50L, D = 4L, seed = 1L,
                           lambda_grid = c(0, 0.05, 0.1, 0.25, 0.5, 1, 2, 4),
                           K = 5L, T_max = 4000L, tol = 1e-6,
                           project_root = gfm_project_root()) {
  gfm_require_packages("ggplot2")
  g_rook <- build_grid_graph(nrow_grid, ncol_grid)
  g_king <- build_king_grid_graph(nrow_grid, ncol_grid)
  D_rook <- build_incidence_matrix(g_rook)
  D_king <- build_incidence_matrix(g_king)
  truth <- generate_true_beta_X(g_rook, D, nrow_grid, ncol_grid)
  Y <- .gfm_with_seed(seed, simulate_data(truth$P_true,
                                         sample(150:300, igraph::vcount(g_rook), replace = TRUE)))
  oracle <- numeric(length(lambda_grid))
  warm <- NULL
  for (index in order(lambda_grid)) {
    fit <- fit_graph_fused_multinomial(Y, D_rook, lambda_grid[index],
                                       T_max = T_max, tol = tol, beta_init = warm, verbose = FALSE)
    warm <- fit$beta_hat
    oracle[index] <- mean((fit$P_hat - truth$P_true)^2)
  }
  cv <- cv_graph_fused_multinomial(Y, D_rook, lambda_grid, K = K, seed = 42L,
                                   T_max = T_max, tol = tol, verbose = FALSE)
  fit_rook <- fit_graph_fused_multinomial(Y, D_rook, cv$lambda_min,
                                         T_max = T_max + 1000L, tol = tol, verbose = FALSE)
  fit_king <- fit_graph_fused_multinomial(Y, D_king, cv$lambda_min,
                                         T_max = T_max + 1000L, tol = tol, verbose = FALSE)
  metrics <- data.frame(method = c("GFM rook", "GFM king"),
                         mse = c(mean((fit_rook$P_hat - truth$P_true)^2),
                                  mean((fit_king$P_hat - truth$P_true)^2)),
                         lambda = cv$lambda_min,
                         converged = c(fit_rook$converged, fit_king$converged))
  result <- list(Y = Y, truth = truth, g_rook = g_rook, g_king = g_king,
                   fit_rook = fit_rook, fit_king = fit_king, cv = cv, metrics = metrics,
                   oracle = data.frame(lambda = lambda_grid, mse = oracle),
                   settings = list(nrow_grid = nrow_grid, ncol_grid = ncol_grid,
                                     seed = seed, D = D, lambda_grid = lambda_grid,
                                     K = K, cv_seed = 42L, T_max = T_max, tol = tol))
  paths <- gfm_save_workflow(result, "grid_demo", project_root)
  coords <- cbind(igraph::V(g_rook)$col, igraph::V(g_rook)$row)
  surface <- plot_probability_surfaces(coords, list(True = truth$P_true,
                                                    Rook = fit_rook$P_hat, King = fit_king$P_hat),
                                        all_categories = FALSE)
  gfm_save_plot(surface, file.path(paths$figures, "x_recovery.png"), 10, 4)
  convergence <- data.frame(iteration = seq_along(fit_rook$resid_hist),
                             residual = pmax(fit_rook$resid_hist, .Machine$double.xmin))
  curve <- ggplot2::ggplot(convergence, ggplot2::aes(iteration, residual)) +
    ggplot2::geom_line(color = "steelblue") + ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Iteration", y = "Relative iterate change") + ggplot2::theme_minimal()
  gfm_save_plot(curve, file.path(paths$figures, "convergence.png"), 7, 5)
  oracle_plot <- ggplot2::ggplot(result$oracle, ggplot2::aes(lambda, mse)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::theme_minimal()
  gfm_save_plot(oracle_plot, file.path(paths$figures, "oracle_curve.png"), 7, 5)
  cv_data <- data.frame(lambda = lambda_grid, value = cv$cv_mean, se = cv$cv_se)
  cv_plot <- ggplot2::ggplot(cv_data, ggplot2::aes(lambda, value)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = value - se, ymax = value + se), width = 0.05) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::labs(y = "Mean held-out negative log-likelihood per trial") + ggplot2::theme_minimal()
  gfm_save_plot(cv_plot, file.path(paths$figures, "cv_curve.png"), 7, 5)
  print(metrics)
  invisible(result)
}

if (sys.nframe() == 0L) {
  flags <- gfm_cli_flags(allowed = "--quick")
  if (flags$quick) run_grid_demo(5, 5, K = 2, lambda_grid = c(0.5, 2.5), T_max = 4000, tol = 1e-4)
  else run_grid_demo()
}


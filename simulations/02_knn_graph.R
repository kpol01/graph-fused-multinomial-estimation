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

#' Compare two union kNN graphs on the same data and a small-k selected lambda.
#' This is the original fixed-k comparison, not joint selection of (lambda, k).
run_knn_demo <- function(n_pts = 2000L, D = 4L, k_small = 5L, k_large = 15L,
                          seed = 7L, lambda_grid = c(0, 0.05, 0.1, 0.25, 0.5, 1, 2, 4),
                          K = 5L, T_max = 4000L, tol = 1e-6,
                          project_root = gfm_project_root()) {
  gfm_require_packages("ggplot2")
  .gfm_check_scalar(n_pts, "n_pts", lower = max(k_small, k_large) + 1L, integer = TRUE)
  data <- .gfm_with_seed(seed, {
    coords <- cbind(x = stats::runif(n_pts), y = stats::runif(n_pts))
    truth <- generate_true_beta_X_coords(coords, D)
    Y <- simulate_data(truth$P_true, sample(150:300, n_pts, replace = TRUE))
    list(coords = coords, truth = truth, Y = Y)
  })
  g_small <- build_knn_graph(data$coords, k_small)
  g_large <- build_knn_graph(data$coords, k_large)
  D_small <- build_incidence_matrix(g_small)
  D_large <- build_incidence_matrix(g_large)
  cv <- cv_graph_fused_multinomial(data$Y, D_small, lambda_grid, K = K, seed = 11L,
                                   T_max = T_max, tol = tol, verbose = FALSE)
  small <- fit_graph_fused_multinomial(data$Y, D_small, cv$lambda_min,
                                      T_max = T_max + 1000L, tol = tol, verbose = FALSE)
  large <- fit_graph_fused_multinomial(data$Y, D_large, cv$lambda_min,
                                      T_max = T_max + 1000L, tol = tol, verbose = FALSE)
  metrics <- data.frame(k = c(k_small, k_large), lambda = cv$lambda_min,
                         mse = c(mean((small$P_hat - data$truth$P_true)^2),
                                  mean((large$P_hat - data$truth$P_true)^2)),
                         converged = c(small$converged, large$converged))
  result <- c(data, list(g_small = g_small, g_large = g_large, cv = cv,
                         fit_small = small, fit_large = large, metrics = metrics,
                         settings = list(n_pts = n_pts, k_small = k_small, k_large = k_large,
                                           seed = seed, D = D, lambda_grid = lambda_grid,
                                           K = K, cv_seed = 11L, T_max = T_max, tol = tol)))
  paths <- gfm_save_workflow(result, "knn_demo", project_root)
  recovery <- plot_probability_surfaces(data$coords,
    list(True = data$truth$P_true, Small = small$P_hat, Large = large$P_hat),
    tile = FALSE, all_categories = FALSE,
    labels = c("True", paste("k =", k_small), paste("k =", k_large)))
  gfm_save_plot(recovery, file.path(paths$figures, "x_recovery.png"), 10, 4)
  edges <- igraph::as_edgelist(g_small, names = FALSE)
  graph_data <- data.frame(x = data$coords[edges[, 1], 1], y = data$coords[edges[, 1], 2],
                            xend = data$coords[edges[, 2], 1], yend = data$coords[edges[, 2], 2])
  graph_plot <- ggplot2::ggplot(graph_data, ggplot2::aes(x, y, xend = xend, yend = yend)) +
    ggplot2::geom_segment(color = "grey70", linewidth = 0.15) +
    ggplot2::geom_point(data = as.data.frame(data$coords), ggplot2::aes(x, y),
                        inherit.aes = FALSE, size = 0.5) +
    ggplot2::coord_equal() + ggplot2::theme_minimal()
  gfm_save_plot(graph_plot, file.path(paths$figures, "graph.png"), 7, 6)
  print(metrics)
  invisible(result)
}

if (sys.nframe() == 0L) {
  flags <- gfm_cli_flags(allowed = "--quick")
  if (flags$quick) run_knn_demo(40, k_small = 3, k_large = 6, K = 2,
                               lambda_grid = c(0.5, 2.5), T_max = 4000, tol = 1e-4)
  else run_knn_demo()
}


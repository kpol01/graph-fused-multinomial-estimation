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

#' Illustrate the repaired five-block SBM with the original noisy block signal.
run_sbm_demo <- function(block_sizes = c(300, 250, 400, 350, 300), D = 4L, seed = 21L,
                          p1_by_block = c(0.75, 0.10, 0.55, 0.20, 0.40),
                          lambda_grid = c(0, 0.05, 0.1, 0.25, 0.5, 1, 2, 4),
                          K = 5L, T_max = 4000L, tol = 1e-6,
                          project_root = gfm_project_root()) {
  gfm_require_packages("ggplot2")
  data <- .gfm_with_seed(seed, {
    g <- build_sbm_graph(block_sizes, p_in = 0.03, p_out = 0.001)
    truth <- generate_true_beta_blocks(g, D, p1_by_block)
    Y <- simulate_data(truth$P_true, sample(150:300, igraph::vcount(g), replace = TRUE))
    list(g = g, truth = truth, Y = Y)
  })
  D_G <- build_incidence_matrix(data$g)
  cv <- cv_graph_fused_multinomial(data$Y, D_G, lambda_grid, K = K, seed = 5L,
                                   T_max = T_max, tol = tol, verbose = FALSE)
  fused <- fit_graph_fused_multinomial(data$Y, D_G, cv$lambda_min,
                                      T_max = T_max + 1000L, tol = tol, verbose = FALSE)
  unfused <- fit_graph_fused_multinomial(data$Y, D_G, 0,
                                        T_max = T_max + 1000L, tol = tol, verbose = FALSE)
  metrics <- data.frame(method = c("Unregularized", "GFM"),
                         lambda = c(0, cv$lambda_min),
                         mse = c(mean((unfused$P_hat - data$truth$P_true)^2),
                                  mean((fused$P_hat - data$truth$P_true)^2)),
                         converged = c(unfused$converged, fused$converged))
  layout <- .gfm_with_seed(99L, igraph::layout_with_fr(data$g))
  result <- c(data, list(cv = cv, fit = fused, fit_unregularized = unfused,
                         layout = layout, metrics = metrics,
                         settings = list(block_sizes = block_sizes, p1_by_block = p1_by_block,
                                           p_in = 0.03, p_out = 0.001, within_block_sd = 0.35,
                                           seed = seed, D = D, lambda_grid = lambda_grid,
                                           K = K, cv_seed = 5L, layout_seed = 99L,
                                           T_max = T_max, tol = tol)))
  paths <- gfm_save_workflow(result, "sbm_demo", project_root)
  recovery <- plot_probability_surfaces(layout,
    list(True = data$truth$P_true, Unregularized = unfused$P_hat, GFM = fused$P_hat),
    tile = FALSE, all_categories = FALSE)
  gfm_save_plot(recovery, file.path(paths$figures, "recovery.png"), 10, 4)
  community <- data.frame(x = layout[, 1], y = layout[, 2],
                           block = factor(igraph::V(data$g)$block))
  community_plot <- ggplot2::ggplot(community, ggplot2::aes(x, y, color = block)) +
    ggplot2::geom_point(size = 0.8) + ggplot2::coord_equal() + ggplot2::theme_minimal()
  gfm_save_plot(community_plot, file.path(paths$figures, "communities.png"), 7, 6)
  print(metrics)
  invisible(result)
}

if (sys.nframe() == 0L) {
  flags <- gfm_cli_flags(allowed = "--quick")
  if (flags$quick) run_sbm_demo(rep(8L, 5), K = 2, lambda_grid = c(0.5, 2.5),
                               T_max = 4000, tol = 1e-4)
  else run_sbm_demo()
}


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

#' Fit training counts and score against raw held-out proportions.
#' Default data and MCAR settings match the supplied comparison script.
run_chicago_comparison <- function(lambda = 1, train_fraction = 0.8, split_seed = 123L,
                                    include_mcar = TRUE, T_max = 5000L, tol = 1e-6,
                                    burnin = 5000L, n_sample = 20000L, thin = 10L,
                                    largest_component = FALSE,
                                    project_root = gfm_project_root()) {
  gfm_require_packages(c("sf", "ggplot2"))
  .gfm_check_flag(include_mcar, "include_mcar")
  .gfm_check_flag(largest_component, "largest_component")
  cache <- file.path(project_root, "data", "chicago_crime_data")
  prepared <- file.path(cache, "chicago_2024_grid1000m_min1_all.rds")
  if (file.exists(prepared)) {
    chicago <- split_chicago_crimes(readRDS(prepared), train_fraction, split_seed)
  } else {
    chicago <- prepare_chicago_crimes(output_dir = cache,
                                      train_fraction = train_fraction, split_seed = split_seed)
  }
  if (largest_component) chicago <- chicago_largest_component(chicago)
  training_totals <- rowsum(chicago$Y_train, chicago$component)
  if (any(training_totals == 0)) {
    warning("A Chicago component has zero training counts in a category; ",
              "a finite unconstrained logit minimizer need not exist. ",
              "Inspect component totals or explicitly choose largest_component = TRUE.",
              call. = FALSE)
  }
  gfm_time <- system.time(fit <- fit_graph_fused_multinomial(
    chicago$Y_train, chicago$D_G, lambda, T_max = T_max, tol = tol, verbose = FALSE
  ))
  if (!fit$converged) {
    warning("GFM reached T_max. Saved metrics and maps describe the final iterate; ",
              "inspect residual and convergence before interpreting the fit.", call. = FALSE)
  }
  probabilities <- list(GFM = fit$P_hat)
  elapsed <- c(GFM = unname(gfm_time["elapsed"]))
  mcar <- NULL
  if (include_mcar) {
    mcar <- fit_mcar_multinomial(chicago$Y_train, chicago$W, burnin, n_sample, thin,
                                seed = split_seed + 1L)
    probabilities$MCAR <- mcar$P_hat
    elapsed["MCAR"] <- mcar$elapsed_seconds
  }
  metrics <- do.call(rbind, lapply(names(probabilities), function(method) {
    data.frame(method = method, score_multinomial(chicago$Y_test, probabilities[[method]]),
                 elapsed_seconds = unname(elapsed[method]),
                 converged = if (method == "GFM") fit$converged else NA)
  }))
  result <- list(chicago = chicago, fit = fit, mcar = mcar,
                   probabilities = probabilities, metrics = metrics,
                   settings = list(lambda = lambda, train_fraction = train_fraction,
                                     split_seed = split_seed, include_mcar = include_mcar,
                                     largest_component = largest_component,
                                     T_max = T_max, tol = tol,
                                     burnin = burnin, n_sample = n_sample, thin = thin,
                                     mcar_seed = split_seed + 1L))
  name <- if (largest_component) "chicago_largest_component" else "chicago_comparison"
  paths <- gfm_save_workflow(result, name, project_root)
  plot <- plot_chicago_comparison(chicago, probabilities)
  gfm_save_plot(plot, file.path(paths$figures, "chicago_probabilities.pdf"), 6.69, 9)
  print(metrics, row.names = FALSE)
  invisible(result)
}

if (sys.nframe() == 0L) {
  flags <- gfm_cli_flags(allowed = c("--quick", "--no-mcar", "--largest-component"))
  largest <- "--largest-component" %in% commandArgs(trailingOnly = TRUE)
  if (flags$quick) run_chicago_comparison(include_mcar = flags$include_mcar,
                                          T_max = 2000, tol = 1e-4, largest_component = largest)
  else run_chicago_comparison(include_mcar = flags$include_mcar, largest_component = largest)
}


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

#' Plot the supplied summary tables without rescaling means, axes, or SDs.
run_summary_plots <- function(project_root = gfm_project_root(), uncertainty_multiplier = 1) {
  gfm_require_packages(c("ggplot2", "patchwork"))
  paths <- gfm_output_dirs(project_root, "summary_curves")
  regimes <- list(n = list(suffix = "n", x = "n"),
                   categories = list(suffix = "D", x = "D"),
                   joint = list(suffix = "n_D", x = "n"))
  output <- list()
  for (name in names(regimes)) {
    regime <- regimes[[name]]
    tables <- setNames(lapply(c("grid", "knn", "sbm"), function(graph) {
      utils::read.csv(file.path(project_root, "results",
                                paste0("mse_", graph, "_", regime$suffix, ".csv")),
                      check.names = FALSE)
    }), c("Grid graph", "kNN graph", "SBM graph"))
    plot <- plot_summary_curves(tables, x_col = regime$x, y_label = "Reported error",
                                uncertainty_multiplier = uncertainty_multiplier)
    path <- file.path(paths$figures, paste0("error_", name, ".pdf"))
    gfm_save_plot(plot, path, 17, 14, units = "cm")
    output[[name]] <- path
  }
  timing <- setNames(lapply(c("grid", "knn", "sbm"), function(graph) {
    utils::read.csv(file.path(project_root, "results", paste0("t_", graph, "_n.csv")),
                    check.names = FALSE)
  }), c("Grid graph", "kNN graph", "SBM graph"))
  plot <- plot_summary_curves(timing, x_col = "n", y_label = "Time (seconds)", log_y = TRUE,
                              uncertainty_multiplier = uncertainty_multiplier)
  output$timing <- file.path(paths$figures, "time_n.pdf")
  gfm_save_plot(plot, output$timing, 17, 14, units = "cm")
  writeLines(capture.output(utils::sessionInfo()), file.path(paths$results, "session_info.txt"))
  message("Created ", length(output), " summary figures in ", paths$figures)
  invisible(output)
}

if (sys.nframe() == 0L) {
  gfm_cli_flags(allowed = character())
  run_summary_plots()
}


# Source this file once to load all reusable GFM functions.
# Paths are resolved relative to this file, independently of getwd().
# Packages are checked but not attached, and temporary loader variables
# stay out of the calling environment. No packages are installed here.

local({
  target_env <- parent.env(environment())
  source_paths <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1L))
  source_paths <- source_paths[nzchar(source_paths)]
  if (length(source_paths) == 0L) {
    stop("Load the project with source('/path/to/R/load_all.R').", call. = FALSE)
  }
  loader_path <- tail(source_paths, 1L)
  if (!file.exists(loader_path) && file.exists(basename(loader_path))) {
    loader_path <- basename(loader_path)  # source(..., chdir = TRUE)
  }
  script_dir <- dirname(normalizePath(loader_path, mustWork = TRUE))
  assign(".gfm_project_path", dirname(script_dir), envir = target_env)
  required_packages <- c("igraph", "Matrix", "FNN")
  installed <- vapply(required_packages, requireNamespace, logical(1L),
                      quietly = TRUE)
  if (any(!installed)) {
    missing <- required_packages[!installed]
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         ".\nInstall them with install.packages(c(",
         paste(sprintf('"%s"', missing), collapse = ", "), ")).", call. = FALSE)
  }
  scripts <- c("model_utils.R", "project_utils.R", "graph_utils.R",
               "simulation_utils.R", "fit_graph_fused_multinomial.R",
               "cross_validation.R", "evaluation_utils.R", "plotting_utils.R",
               "prepare_chicago_crimes.R", "prepare_india_mother_tongue.R")
  for (script in scripts) {
    source(file.path(script_dir, script), local = target_env)
  }
})

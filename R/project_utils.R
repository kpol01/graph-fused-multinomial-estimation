# Project paths, optional dependencies, and workflow output helpers.

gfm_require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         ". Install the corresponding dependency group; see README.md.", call. = FALSE)
  }
  invisible(packages)
}

#' Locate the project from a directory or a file inside its tree.
#' @param path Optional starting path; defaults to the path recorded by load_all.R.
#' @return The absolute project root, recognized by R/load_all.R and simulations/.
gfm_project_root <- function(path = NULL) {
  if (is.null(path)) {
    path <- get0(".gfm_project_path", envir = environment(gfm_project_root),
                 inherits = FALSE, ifnotfound = getwd())
  }
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!dir.exists(path)) path <- dirname(path)
  repeat {
    if (file.exists(file.path(path, "R", "load_all.R")) &&
        dir.exists(file.path(path, "simulations"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  stop("Cannot locate the GFM project root from the supplied path.", call. = FALSE)
}

gfm_output_dirs <- function(project_root = gfm_project_root(), name) {
  .gfm_check_scalar(nchar(name), "name length", lower = 1)
  if (length(name) != 1L || !grepl("^[A-Za-z0-9_-]+$", name)) {
    stop("Use a simple workflow name containing letters, numbers, '_' or '-'.",
         call. = FALSE)
  }
  paths <- list(results = file.path(project_root, "results", "generated", name),
                figures = file.path(project_root, "figures", "generated", name))
  for (path in paths) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(path)) stop("Cannot create output directory: ", path, call. = FALSE)
  }
  paths
}

gfm_save_workflow <- function(result, name, project_root = gfm_project_root()) {
  paths <- gfm_output_dirs(project_root, name)
  saveRDS(result, file.path(paths$results, "analysis.rds"))
  if (!is.null(result$metrics)) {
    utils::write.csv(result$metrics, file.path(paths$results, "metrics.csv"), row.names = FALSE)
  }
  writeLines(capture.output(utils::sessionInfo()), file.path(paths$results, "session_info.txt"))
  invisible(paths)
}

gfm_cli_flags <- function(allowed = c("--quick", "--no-mcar")) {
  args <- commandArgs(trailingOnly = TRUE)
  unknown <- setdiff(args, allowed)
  if (length(unknown) > 0L) {
    stop("Unknown arguments: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  list(quick = "--quick" %in% args, include_mcar = !"--no-mcar" %in% args)
}

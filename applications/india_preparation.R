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

#' Prepare the cached full 2011 district dataset, with optional explicit selection.
#' Preparation and geometry checks do not fit either statistical model.
run_india_preparation <- function(largest_component = FALSE, train_fraction = 0.8,
                                   split_seed = 123L, project_root = gfm_project_root()) {
  .gfm_check_flag(largest_component, "largest_component")
  india <- prepare_india_mother_tongue(
    cache_dir = file.path(project_root, "data", "india_mother_tongue_cache"),
    train_fraction = train_fraction, split_seed = split_seed
  )
  if (largest_component) india <- india_largest_component(india)
  name <- if (largest_component) "india_largest_component" else "india_preparation"
  paths <- gfm_output_dirs(project_root, name)
  saveRDS(india, file.path(paths$results, "prepared_india.rds"))
  utils::write.csv(india$component_summary, file.path(paths$results, "component_summary.csv"),
                   row.names = FALSE)
  counts <- data.frame(india$districts, india$Y, check.names = FALSE)
  utils::write.csv(counts, file.path(paths$results, "district_counts.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(paths$results, "session_info.txt"))
  message("Saved prepared India data to ", paths$results)
  invisible(india)
}

if (sys.nframe() == 0L) {
  gfm_cli_flags(allowed = "--largest-component")
  run_india_preparation(largest_component = "--largest-component" %in% commandArgs(trailingOnly = TRUE))
}


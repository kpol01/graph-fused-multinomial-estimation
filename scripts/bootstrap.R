# Internal bootstrap for entry scripts, independent of the working directory.
# This file is sourced by an entry script; it does not fit models.
local({
  target_env <- parent.env(environment())
  paths <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1L))
  paths <- paths[nzchar(paths)]
  current <- tail(paths, 1L)
  root <- dirname(dirname(normalizePath(current, mustWork = TRUE)))
  source(file.path(root, "R", "load_all.R"), local = target_env)
})

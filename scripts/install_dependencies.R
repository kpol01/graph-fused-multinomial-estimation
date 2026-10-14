# Explicit, optional installation step. Analysis scripts never install packages.
# Usage: Rscript scripts/install_dependencies.R core plots data
# Add "mcar" for CARBayes or "latex" for optional LaTeX table export.
groups <- list(core = c("Matrix", "igraph", "FNN"),
                 plots = c("ggplot2", "patchwork"), data = c("sf", "readxl"),
                 mcar = "CARBayes", latex = "xtable")
chosen <- commandArgs(trailingOnly = TRUE)
if (length(chosen) == 0L) chosen <- c("core", "plots", "data")
if (any(!chosen %in% names(groups))) {
  stop("Choose dependency groups: ", paste(names(groups), collapse = ", "), call. = FALSE)
}
packages <- unique(unlist(groups[chosen], use.names = FALSE))
missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing) > 0L) {
  utils::install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All requested packages are already installed.")
}

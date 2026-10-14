# Shared publication plots. Optional plotting packages load only when called.

gfm_save_plot <- function(plot, path, width, height, units = "in") {
  gfm_require_packages("ggplot2")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height,
                  units = units, dpi = 300, bg = "white")
  invisible(path)
}

plot_probability_surfaces <- function(coords, probabilities, tile = TRUE,
                                       all_categories = TRUE, labels = names(probabilities)) {
  gfm_require_packages("ggplot2")
  .gfm_check_matrix(coords, "coords", min_cols = 2L)
  if (!is.list(probabilities) || length(probabilities) < 1L ||
      length(labels) != length(probabilities) || anyNA(labels) || anyDuplicated(labels)) {
    stop("Supply probability matrices and distinct panel labels.", call. = FALSE)
  }
  template <- probabilities[[1L]]
  category_names <- colnames(template)
  if (is.null(category_names) || anyNA(category_names) || any(!nzchar(category_names)) ||
      anyDuplicated(category_names)) category_names <- paste("Category", seq_len(ncol(template)))
  categories <- if (all_categories) seq_len(ncol(template)) else 1L
  data <- do.call(rbind, lapply(seq_along(probabilities), function(index) {
    P <- probabilities[[index]]
    .gfm_check_matrix(P, "probability surface")
    if (!identical(dim(P), dim(template)) || nrow(P) != nrow(coords) ||
        any(P < 0 | P > 1) || any(abs(rowSums(P) - 1) > 1e-6)) {
      stop("All probability surfaces must match the coordinates and categories.", call. = FALSE)
    }
    data.frame(x = rep(coords[, 1L], length(categories)),
                y = rep(coords[, 2L], length(categories)),
                value = as.vector(P[, categories, drop = FALSE]),
                category = rep(category_names[categories], each = nrow(P)),
                method = labels[index])
  }))
  data$category <- factor(data$category, levels = category_names[categories])
  data$method <- factor(data$method, levels = labels)
  if (tile) {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x, y, fill = value)) +
      ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(option = "C", limits = c(0, 1),
                                                         name = "Probability") +
      ggplot2::scale_y_reverse()
  } else {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x, y, color = value)) +
      ggplot2::geom_point(size = 0.8) +
      ggplot2::scale_color_viridis_c(option = "C", limits = c(0, 1), name = "Probability")
  }
  plot + ggplot2::facet_grid(category ~ method, switch = "y") +
    ggplot2::coord_equal() + ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                    strip.text.x = ggplot2::element_text(face = "bold", size = 9),
                    strip.text.y.left = ggplot2::element_text(angle = 0),
                    axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                    legend.position = "bottom", strip.placement = "outside")
}

#' Plot held-out observed Chicago proportions and the supplied training fits.
#' Labels use freshly calculated MSE; clipping affects display geometry only.
plot_chicago_comparison <- function(chicago, probabilities) {
  gfm_require_packages(c("sf", "ggplot2"))
  if (!is.list(probabilities) || is.null(names(probabilities)) ||
      any(!nzchar(names(probabilities))) || anyDuplicated(names(probabilities))) {
    stop("Supply a named list of fitted probability matrices.", call. = FALSE)
  }
  Y <- chicago$Y_test
  .gfm_check_counts(Y)
  if (!identical(chicago$grid_sf$cell_key, rownames(Y))) {
    stop("Chicago geometry and count rows must follow the same cell order.", call. = FALSE)
  }
  active <- rowSums(Y) > 0
  if (!any(active)) stop("There are no test crimes to plot.", call. = FALSE)
  fitted <- lapply(probabilities, gfm_align_probabilities, Y = Y)
  mse <- vapply(fitted, function(P) score_multinomial(Y, P)$mse, numeric(1L))
  matrices <- c(list(Observed = multinomial_mle(Y)), fitted)
  panel_labels <- c("Observed (test)", sprintf("%s\n(MSE = %.5f)", names(fitted), mse))
  city <- sf::st_make_valid(sf::st_transform(chicago$city_boundary, 26916))
  grid <- sf::st_transform(chicago$grid_sf, sf::st_crs(city))
  grid$node_index <- seq_len(nrow(grid))
  sf::st_agr(grid) <- "constant"
  clipped <- suppressWarnings(sf::st_intersection(sf::st_make_valid(grid),
                                                   sf::st_geometry(city)))
  if (any(sf::st_geometry_type(clipped) == "GEOMETRYCOLLECTION")) {
    clipped <- sf::st_collection_extract(clipped, "POLYGON", warn = FALSE)
  }
  clipped <- clipped[!sf::st_is_empty(clipped) & sf::st_dimension(clipped) == 2L, ]
  if (nrow(clipped) == 0L) stop("No grid polygons overlap the city.", call. = FALSE)
  categories <- colnames(Y)
  layers <- list()
  for (j in seq_along(categories)) for (k in seq_along(matrices)) {
    values <- matrices[[k]][, j]
    values[!active] <- NA_real_
    layers[[length(layers) + 1L]] <- sf::st_sf(
      category = categories[j], method = panel_labels[k],
      probability = values[clipped$node_index], geometry = sf::st_geometry(clipped)
    )
  }
  data <- do.call(rbind, layers)
  data$category <- factor(data$category, levels = categories)
  data$method <- factor(data$method, levels = panel_labels)
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = city, fill = "grey90", color = NA) +
    ggplot2::geom_sf(data = data, ggplot2::aes(fill = probability), color = NA) +
    ggplot2::geom_sf(data = city, fill = NA, color = "grey30", linewidth = 0.2) +
    ggplot2::facet_grid(category ~ method) +
    ggplot2::scale_fill_viridis_c(option = "C", limits = c(0, 1), na.value = "grey90",
                                  name = "Probability") +
    ggplot2::coord_sf(datum = NA) + ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", strip.text = ggplot2::element_text(size = 9),
                    panel.spacing = grid::unit(0.15, "cm"),
                    plot.margin = ggplot2::margin(5, 5, 5, 5))
}

gfm_summary_long <- function(table, x_col) {
  methods <- c(gfm = "GFM", mcar = "MCAR", mmle = "MLE")
  required <- c(x_col, names(methods), paste0(names(methods), "_sig"))
  if (!is.data.frame(table) || !all(required %in% names(table)) || nrow(table) < 1L ||
      any(!vapply(table[required], is.numeric, logical(1L))) ||
      any(!is.finite(as.matrix(table[required])))) {
    stop("Summary table schema must include ", paste(required, collapse = ", "), ".",
         call. = FALSE)
  }
  if (any(as.matrix(table[paste0(names(methods), "_sig")]) < 0)) {
    stop("Stored uncertainty values must be nonnegative.", call. = FALSE)
  }
  do.call(rbind, lapply(names(methods), function(method) {
    data.frame(x = table[[x_col]], value = table[[method]],
                sd = table[[paste0(method, "_sig")]], method = methods[[method]])
  }))
}

#' Plot summary means and their stored uncertainty without changing CSV values.
#' @param tables Named list of grid, kNN, and SBM summary tables.
#' @param x_col Literal CSV axis name, either n or D.
#' @param uncertainty_multiplier Explicit multiplier; default is one.
plot_summary_curves <- function(tables, x_col = "n", y_label = "Reported error",
                                 log_y = FALSE, uncertainty_multiplier = 1) {
  gfm_require_packages(c("ggplot2", "patchwork"))
  .gfm_check_scalar(uncertainty_multiplier, "uncertainty_multiplier", lower = 0)
  long <- lapply(tables, gfm_summary_long, x_col = x_col)
  range_all <- unlist(lapply(long, function(data) {
    c(data$value - uncertainty_multiplier * data$sd,
      data$value + uncertainty_multiplier * data$sd)
  }))
  if (log_y && any(range_all <= 0)) {
    stop("Log-scale mean +/- uncertainty intervals must be positive.", call. = FALSE)
  }
  limits <- if (log_y) range(range_all) * c(0.8, 1.2) else c(0, max(range_all) * 1.05)
  plots <- lapply(seq_along(long), function(index) {
    data <- long[[index]]
    data$method <- factor(data$method, levels = c("GFM", "MCAR", "MLE"))
    plot <- ggplot2::ggplot(data, ggplot2::aes(x, value, color = method,
                                              shape = method, group = method)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = value - uncertainty_multiplier * sd,
                                         ymax = value + uncertainty_multiplier * sd),
                              width = max(diff(range(data$x)), 1) * 0.015, linewidth = 0.35) +
      ggplot2::geom_line(linewidth = 0.45) + ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = c(GFM = "#0072B2", MCAR = "#D55E00", MLE = "#009E73")) +
      ggplot2::scale_shape_manual(values = c(GFM = 16, MCAR = 17, MLE = 15)) +
      ggplot2::scale_x_continuous(breaks = sort(unique(data$x))) +
      ggplot2::labs(title = names(tables)[index], x = x_col, y = y_label,
                      color = NULL, shape = NULL) + ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
    if (log_y) plot <- plot + ggplot2::scale_y_log10()
    plot + ggplot2::coord_cartesian(ylim = limits)
  })
  patchwork::wrap_plots(plots, design = "AABB\n#CC#", guides = "collect") +
    patchwork::plot_annotation(theme = ggplot2::theme(legend.position = "bottom"))
}

# Chicago crime counts and a common graph for GFM/MCAR analyses.
# Source this file to define the functions; downloads occur only on preparation.
# Raw zeros and the OTHER reference category are preserved. Distances use
# NAD83 / UTM zone 16N (EPSG:26916); adjacency is rook contiguity.
# Source records: https://data.cityofchicago.org/d/ijzp-q8t2
# City boundary: https://data.cityofchicago.org/d/qqq8-j68g
# See PREPARATION_GUIDE.md for setup, examples, and interpretation.

.chicago_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Install these R packages first: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
}

.chicago_check_scalar <- function(x, name, lower, upper = Inf, integer = FALSE) {
  valid <- is.numeric(x) && !is.complex(x) && length(x) == 1L && is.finite(x) &&
    x >= lower && x <= upper && (!integer || x == floor(x))
  if (!valid) {
    stop(name, " must be a finite ", if (integer) "integer" else "number",
         " in [", lower, ", ", upper, "].", call. = FALSE)
  }
  invisible(x)
}

.chicago_check_string <- function(x, name, allow_empty = FALSE) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      (!allow_empty && !nzchar(trimws(x)))) {
    stop(name, " must be a single ", if (!allow_empty) "nonempty ",
         "character string.", call. = FALSE)
  }
  invisible(x)
}

.chicago_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(x)
}

.chicago_check_split <- function(train_fraction, split_seed) {
  .chicago_check_scalar(train_fraction, "train_fraction", lower = 0, upper = 1)
  if (train_fraction == 0 || train_fraction == 1) {
    stop("Use 0 < train_fraction < 1.", call. = FALSE)
  }
  .chicago_check_scalar(split_seed, "split_seed", lower = 0,
                        upper = .Machine$integer.max, integer = TRUE)
}

.chicago_check_counts <- function(Y) {
  if (!is.matrix(Y) || !is.numeric(Y) || is.complex(Y) || nrow(Y) < 1L ||
      ncol(Y) < 2L || any(!is.finite(Y)) || any(Y < 0 | Y != floor(Y))) {
    stop("Y must be a nonempty numeric matrix of nonnegative integer counts, ",
         "with at least two categories.", call. = FALSE)
  }
  invisible(Y)
}

.chicago_require_columns <- function(data, columns, label) {
  if (!is.data.frame(data) || !all(columns %in% names(data))) {
    stop(label, " must contain: ", paste(columns, collapse = ", "), ".",
         call. = FALSE)
  }
}

.chicago_check_snapshot <- function(snapshot, period, api) {
  required <- c("data", "period", "api", "expected_records", "downloaded_at")
  if (!is.list(snapshot) || !all(required %in% names(snapshot)) ||
      !identical(snapshot$period, period) || !identical(snapshot$api, api)) {
    stop("The cached snapshot does not match the requested year/source. ",
         "Use refresh = TRUE to replace it.", call. = FALSE)
  }
  .chicago_require_columns(snapshot$data,
                           c("id", "date", "primary_type", "latitude", "longitude"),
                           "The record snapshot")
  .chicago_check_scalar(snapshot$expected_records, "snapshot$expected_records",
                        lower = 1, integer = TRUE)
  if (nrow(snapshot$data) != snapshot$expected_records ||
      anyNA(snapshot$data$id) || anyDuplicated(snapshot$data$id)) {
    stop("The record snapshot has a count mismatch or missing/duplicate IDs. ",
         "Use refresh = TRUE to replace it.", call. = FALSE)
  }
}

# HTTP and CSV helpers use built-in libcurl, without curl/httr2 packages.
.chicago_request <- function(url, query = list(), app_token = "") {
  .chicago_check_string(url, "url")
  .chicago_check_string(app_token, "app_token", allow_empty = TRUE)
  if (length(query)) {
    encode <- function(x) utils::URLencode(as.character(x), reserved = TRUE)
    keys <- vapply(names(query), encode, character(1), USE.NAMES = FALSE)
    values <- vapply(query, encode, character(1), USE.NAMES = FALSE)
    url <- paste0(url, "?", paste0(keys, "=", values, collapse = "&"))
  }
  old_timeout <- getOption("timeout", 60)
  options(timeout = max(120, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  path <- tempfile(fileext = ".download")
  on.exit(unlink(path), add = TRUE)
  headers <- if (nzchar(app_token)) c(`X-App-Token` = app_token) else NULL
  last_error <- "Empty download."
  for (attempt in seq_len(4L)) {
    unlink(path)
    status <- tryCatch(
      suppressWarnings(utils::download.file(url, path, method = "libcurl",
                                            quiet = TRUE, mode = "wb",
                                            headers = headers)),
      error = function(e) e
    )
    if (!inherits(status, "error") && isTRUE(status == 0L) &&
        file.exists(path) && file.info(path)$size > 0) {
      return(readBin(path, what = "raw", n = file.info(path)$size))
    }
    last_error <- if (inherits(status, "error")) conditionMessage(status) else
      paste("Download returned status", status, "or an empty response.")
    if (attempt < 4L) {
      message("Download failed; retrying (attempt ", attempt + 1L, " of 4).")
      Sys.sleep(2 * attempt)
    }
  }
  stop("Download failed after 4 attempts: ", last_error, call. = FALSE)
}
.chicago_fetch_table <- function(api, query, app_token = "") {
  txt <- rawToChar(.chicago_request(api, query, app_token))
  utils::read.csv(text = txt, colClasses = "character", check.names = FALSE,
                  na.strings = "", stringsAsFactors = FALSE)
}

.chicago_cell_key <- function(x, y) paste(x, y, sep = "_")

.chicago_rook_graph <- function(grid) {
  required <- c("cell_key", "grid_x", "grid_y")
  .chicago_require_columns(grid, required, "grid")
  keys <- grid$cell_key
  if (nrow(grid) < 1L || anyNA(keys) || anyDuplicated(keys) ||
      any(!nzchar(keys)) || any(!is.finite(grid$grid_x)) ||
      any(!is.finite(grid$grid_y))) {
    stop("grid must contain distinct cell keys and finite grid coordinates.",
         call. = FALSE)
  }
  east <- match(.chicago_cell_key(grid$grid_x + 1L, grid$grid_y), keys)
  north <- match(.chicago_cell_key(grid$grid_x, grid$grid_y + 1L), keys)
  from <- rep(seq_len(nrow(grid)), 2L)
  to <- c(east, north)
  ok <- !is.na(to)
  edge_data <- data.frame(from = keys[from[ok]], to = keys[to[ok]])
  igraph::graph_from_data_frame(edge_data, directed = FALSE,
                                vertices = data.frame(name = keys))
}

#' Split observed crimes without replacement within each retained cell.
#' @param chicago Prepared list containing the raw count matrix Y.
#' @param train_fraction Requested training fraction, strictly between 0 and 1.
#' @param split_seed Nonnegative integer seed, local to this split.
#' @return The input list with Y_train, Y_test, their row totals, and split_info.
#'   Singleton cells stay in training. Larger cells have at least one crime in
#'   each set. Graphs, geometries, and full counts keep the same node order.
split_chicago_crimes <- function(chicago, train_fraction = 0.8,
                                 split_seed = 123L) {
  if (!is.list(chicago) || is.null(chicago$Y)) {
    stop("chicago must be a prepared object containing the count matrix $Y.",
         call. = FALSE)
  }
  Y <- as.matrix(chicago$Y)
  .chicago_check_counts(Y)
  .chicago_check_split(train_fraction, split_seed)
  n_i <- rowSums(Y)
  if (any(n_i == 0)) {
    stop("Every retained cell must contain at least one recorded crime.",
         call. = FALSE)
  }

  # Restore the caller's RNG state after generating the reproducible split.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(split_seed)

  # Select the training total from the cell total alone, not its category mix.
  # Rounding and small cells can change the achieved overall fraction.
  train_total <- ifelse(n_i == 1, 1,
                        pmax(1, pmin(n_i - 1, floor(train_fraction * n_i))))
  Y_train <- matrix(0, nrow(Y), ncol(Y), dimnames = dimnames(Y))
  remaining_total <- n_i
  remaining_draws <- train_total
  for (j in seq_len(ncol(Y) - 1L)) {
    active <- remaining_draws > 0
    Y_train[active, j] <- stats::rhyper(
      nn = sum(active), m = Y[active, j],
      n = remaining_total[active] - Y[active, j],
      k = remaining_draws[active]
    )
    remaining_draws <- remaining_draws - Y_train[, j]
    remaining_total <- remaining_total - Y[, j]
  }
  Y_train[, ncol(Y)] <- remaining_draws
  if (is.integer(Y)) storage.mode(Y_train) <- "integer"
  Y_test <- Y - Y_train
  stopifnot(all(Y_train >= 0), all(Y_test >= 0),
            all(Y_train + Y_test == Y), all(rowSums(Y_train) == train_total),
            identical(dimnames(Y_train), dimnames(Y)),
            identical(dimnames(Y_test), dimnames(Y)))

  chicago$Y_train <- Y_train
  chicago$Y_test <- Y_test
  chicago$n_i_train <- rowSums(Y_train)
  chicago$n_i_test <- rowSums(Y_test)
  chicago$split_info <- list(
    method = "Within-cell sampling without replacement",
    train_fraction_requested = train_fraction, split_seed = as.integer(split_seed),
    train_fraction_actual = sum(Y_train) / sum(Y),
    total_crimes = sum(Y), train_crimes = sum(Y_train), test_crimes = sum(Y_test),
    total_cells = nrow(Y), tested_cells = sum(chicago$n_i_test > 0),
    singleton_cells_training_only = sum(n_i == 1)
  )
  chicago
}

#' Download/cache Chicago crimes and construct grid counts and a rook graph.
#' @param year Calendar year, at least 2001.
#' @param cellsize_m Positive square-cell side length in metres.
#' @param min_count Minimum full-data crime count for a retained cell.
#' @param largest_component Whether to retain only the largest graph component.
#' @param output_dir Directory for snapshots, prepared RDS, and CSV outputs.
#' @param refresh Whether to replace cached records and the city boundary.
#' @param app_token Optional Socrata application token; never stored in results.
#' @param train_fraction,split_seed Arguments passed to split_chicago_crimes().
#' @return A list with full/training/test counts, aligned graph matrices and
#'   geometries, exclusion and download audits, and the saved output path.
#'   No models are fitted. The boundary is the public version at download time.
prepare_chicago_crimes <- function(year = 2024L, cellsize_m = 1000,
                                   min_count = 1L, largest_component = FALSE,
                                   output_dir = "chicago_crime_data",
                                   refresh = FALSE,
                                   app_token = Sys.getenv("SOCRATA_APP_TOKEN"),
                                   train_fraction = 0.8, split_seed = 123L) {
  .chicago_require(c("sf", "Matrix", "igraph"))
  .chicago_check_scalar(year, "year", lower = 2001,
                        upper = .Machine$integer.max - 1, integer = TRUE)
  .chicago_check_scalar(cellsize_m, "cellsize_m", lower = .Machine$double.xmin)
  .chicago_check_scalar(min_count, "min_count", lower = 1,
                        upper = .Machine$integer.max, integer = TRUE)
  .chicago_check_flag(largest_component, "largest_component")
  .chicago_check_flag(refresh, "refresh")
  .chicago_check_string(output_dir, "output_dir")
  .chicago_check_string(app_token, "app_token", allow_empty = TRUE)
  .chicago_check_split(train_fraction, split_seed)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) {
    stop("Could not create output_dir: ", output_dir, call. = FALSE)
  }
  api <- "https://data.cityofchicago.org/resource/ijzp-q8t2.csv"
  boundary_url <- "https://data.cityofchicago.org/resource/qqq8-j68g.geojson"
  categories <- c("THEFT", "BATTERY", "CRIMINAL DAMAGE", "ASSAULT", "BURGLARY", "OTHER")
  period <- sprintf("date >= '%d-01-01T00:00:00' AND date < '%d-01-01T00:00:00'",
                    as.integer(year), as.integer(year + 1L))

  fetch_table <- function(query) .chicago_fetch_table(api, query, app_token)

  # Page by unique ID, retaining the original fixed upper-ID/count audit.
  raw_path <- file.path(output_dir, sprintf("raw_chicago_%d.rds", as.integer(year)))
  if (file.exists(raw_path) && !refresh) {
    snapshot <- readRDS(raw_path)
    .chicago_check_snapshot(snapshot, period, api)
    message("Using cached records downloaded on ", snapshot$downloaded_at)
  } else {
    size <- fetch_table(list(`$select` = "count(*) AS n, max(id) AS max_id",
                             `$where` = period))
    .chicago_require_columns(size, c("n", "max_id"), "The API summary")
    expected <- as.numeric(size$n[1L])
    max_id <- as.numeric(size$max_id[1L])
    if (!is.finite(expected) || expected < 1 || expected != floor(expected) ||
        !is.finite(max_id) || max_id != floor(max_id)) {
      stop("The API returned no records for this year.", call. = FALSE)
    }
    upper <- format(max_id, scientific = FALSE, trim = TRUE)
    pages <- list()
    last_id <- -1
    downloaded <- 0L
    repeat {
      where <- paste0(period, " AND id <= ", upper, " AND id > ",
                      format(last_id, scientific = FALSE, trim = TRUE))
      page <- fetch_table(list(`$select` = "id,date,primary_type,latitude,longitude",
                               `$where` = where, `$order` = "id ASC", `$limit` = 50000))
      .chicago_require_columns(page,
                               c("id", "date", "primary_type", "latitude", "longitude"),
                               "The API record page")
      if (!nrow(page)) break
      ids <- suppressWarnings(as.numeric(page$id))
      if (any(!is.finite(ids)) || any(ids != floor(ids)) || any(ids > max_id) ||
          any(ids <= last_id) || anyDuplicated(ids) || is.unsorted(ids)) {
        stop("Unexpected record IDs in the API response.", call. = FALSE)
      }
      pages[[length(pages) + 1L]] <- page
      last_id <- max(ids)
      downloaded <- downloaded + nrow(page)
      message("Downloaded ", format(downloaded, big.mark = ","), " / ",
              format(expected, big.mark = ","), " records")
      if (downloaded >= expected) break
    }
    if (downloaded != expected) {
      stop("Download count differs from the API count. The live data may have ",
           "changed during retrieval; rerun the function.", call. = FALSE)
    }
    raw <- do.call(rbind, pages)
    if (anyDuplicated(raw$id)) stop("Duplicate record IDs in download.", call. = FALSE)
    snapshot <- list(data = raw, api = api, period = period,
                     expected_records = expected, max_id = max_id,
                     downloaded_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
    saveRDS(snapshot, raw_path)
  }
  .chicago_check_snapshot(snapshot, period, api)
  raw <- snapshot$data
  n_downloaded <- nrow(raw)

  # Project into NAD83 / UTM zone 16N: all grid distances below are in METRES.
  # A one-degree grid in longitude/latitude would not be a one-kilometre grid.
  boundary_path <- file.path(output_dir, "chicago_city_boundary.geojson")
  if (!file.exists(boundary_path) || file.info(boundary_path)$size == 0 || refresh) {
    response <- .chicago_request(boundary_url, list(`$limit` = 50000), app_token)
    writeBin(response, boundary_path)
  }
  boundary <- sf::st_read(boundary_path, quiet = TRUE)
  if (!nrow(boundary) || all(sf::st_is_empty(boundary))) {
    stop("The downloaded city boundary has no geometry.", call. = FALSE)
  }
  boundary <- sf::st_make_valid(sf::st_transform(boundary, 26916))
  city <- sf::st_sf(geometry = sf::st_union(sf::st_geometry(boundary)))

  raw$longitude <- suppressWarnings(as.numeric(raw$longitude))
  raw$latitude <- suppressWarnings(as.numeric(raw$latitude))
  valid_xy <- is.finite(raw$longitude) & is.finite(raw$latitude) &
    raw$longitude >= -180 & raw$longitude <= 180 &
    raw$latitude >= -90 & raw$latitude <= 90
  n_invalid_xy <- sum(!valid_xy)
  raw <- raw[valid_xy, , drop = FALSE]
  raw$primary_type <- toupper(trimws(raw$primary_type))
  valid_type <- !is.na(raw$primary_type) & nzchar(raw$primary_type)
  n_missing_type <- sum(!valid_type)
  raw <- raw[valid_type, , drop = FALSE]
  if (!nrow(raw)) stop("No records have valid coordinates and categories.", call. = FALSE)

  # A coarse city bounding-box filter also rejects impossible far-away geocodes
  # before projection. The exact polygon check follows in the projected CRS.
  city_ll_bbox <- sf::st_bbox(sf::st_transform(city, 4326))
  in_bbox <- raw$longitude >= city_ll_bbox["xmin"] &
    raw$longitude <= city_ll_bbox["xmax"] &
    raw$latitude >= city_ll_bbox["ymin"] &
    raw$latitude <= city_ll_bbox["ymax"]
  n_outside <- sum(!in_bbox)
  raw <- raw[in_bbox, , drop = FALSE]
  if (!nrow(raw)) stop("No records fall within the city bounding box.", call. = FALSE)
  points <- sf::st_as_sf(raw, coords = c("longitude", "latitude"), crs = 4326,
                         remove = FALSE)
  points <- sf::st_transform(points, 26916)
  inside <- lengths(sf::st_intersects(points, city)) > 0L
  n_outside <- n_outside + sum(!inside)
  points <- points[inside, , drop = FALSE]
  if (!nrow(points)) stop("No records fall within the city boundary.", call. = FALSE)
  points$category <- ifelse(points$primary_type %in% categories[1:5],
                            points$primary_type, "OTHER")

  # Anchor the lattice at multiples of cellsize_m in the projected CRS.
  # Integer grid coordinates assign boundary points to exactly one cell.
  bb <- sf::st_bbox(city)
  origin <- floor(unname(bb[c("xmin", "ymin")]) / cellsize_m) * cellsize_m
  upper <- (floor(unname(bb[c("xmax", "ymax")]) / cellsize_m) + 1) * cellsize_m
  shape <- as.integer(round((upper - origin) / cellsize_m))
  cells <- sf::st_make_grid(city, cellsize = cellsize_m, offset = origin,
                            n = shape, square = TRUE)
  cells <- cells[lengths(sf::st_intersects(cells, city)) > 0L]
  centers <- sf::st_coordinates(sf::st_centroid(cells))
  grid_x <- as.integer(floor(centers[, 1L] / cellsize_m))
  grid_y <- as.integer(floor(centers[, 2L] / cellsize_m))
  grid <- sf::st_sf(cell_key = .chicago_cell_key(grid_x, grid_y), grid_x = grid_x,
                    grid_y = grid_y, geometry = cells)
  grid <- grid[order(grid$grid_y, grid$grid_x), ]
  xy <- sf::st_coordinates(points)
  point_key <- .chicago_cell_key(floor(xy[, 1L] / cellsize_m), floor(xy[, 2L] / cellsize_m))
  if (anyNA(match(point_key, grid$cell_key))) {
    stop("Some records could not be assigned to a grid cell.", call. = FALSE)
  }
  counts <- table(factor(point_key, levels = grid$cell_key),
                  factor(points$category, levels = categories))
  Y <- matrix(as.integer(counts), nrow = nrow(grid), ncol = length(categories),
              dimnames = list(grid$cell_key, categories))
  stopifnot(sum(Y) == nrow(points))
  audit <- sf::st_drop_geometry(grid)
  audit$n_i <- as.integer(rowSums(Y))
  eligible <- audit$n_i >= min_count
  n_below_min <- sum(audit$n_i[!eligible])
  if (!any(eligible)) stop("No cells meet min_count; lower it or enlarge the cells.",
                           call. = FALSE)
  grid <- grid[eligible, ]
  Y <- Y[eligible, , drop = FALSE]

  # Rook edges only: horizontal/vertical neighbors, no diagonal edges.
  # Explicit vertex names retain isolated nodes and enforce matrix row order.
  g <- .chicago_rook_graph(grid)
  components_before <- igraph::components(g)
  n_excluded_component <- 0L
  if (largest_component && components_before$no > 1L) {
    keep <- components_before$membership == which.max(components_before$csize)
    n_excluded_component <- sum(Y[!keep, , drop = FALSE])
    grid <- grid[keep, ]
    Y <- Y[keep, , drop = FALSE]
    g <- .chicago_rook_graph(grid)
  }
  n <- nrow(Y)
  d <- ncol(Y)
  n_i <- as.integer(rowSums(Y))
  grid$node_id <- seq_len(n)
  grid$n_i <- n_i
  audit$retained <- audit$cell_key %in% grid$cell_key
  coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(grid)))
  rownames(coords) <- grid$cell_key
  colnames(coords) <- c("x_m", "y_m")
  edge_index <- igraph::as_edgelist(g, names = FALSE)
  n_edges <- nrow(edge_index)
  W <- igraph::as_adjacency_matrix(g, sparse = TRUE)
  # D_G has dimensions |E| x n. Each row is +1 at one endpoint, -1 at the other.
  D_G <- Matrix::sparseMatrix(i = rep(seq_len(n_edges), 2L),
                              j = c(edge_index[, 1L], edge_index[, 2L]),
                              x = c(rep(1, n_edges), rep(-1, n_edges)),
                              dims = c(n_edges, n),
                              dimnames = list(NULL, grid$cell_key))
  component <- igraph::components(g)$membership
  component_totals <- rowsum(Y, group = component, reorder = FALSE)
  if (any(component_totals == 0L)) {
    warning("At least one graph component has no records in a category. ",
            "An unbounded logit fit can diverge there; inspect component_totals. ",
            "The raw counts have been preserved.", call. = FALSE)
  }
  # Matrix is loaded with requireNamespace(), not attached with library().
  # Explicit Matrix generics dispatch correctly on sparse S4 matrices; base
  # diag(W) can otherwise try to construct an enormous diagonal matrix.
  stopifnot(identical(rownames(Y), grid$cell_key),
            identical(igraph::V(g)$name, grid$cell_key),
            identical(colnames(D_G), grid$cell_key),
            identical(rownames(W), grid$cell_key),
            all(n_i >= min_count), all(Matrix::rowSums(W) <= 4),
            all(Matrix::diag(W) == 0), isTRUE(Matrix::isSymmetric(W)),
            all(as.numeric(D_G %*% rep(1, n)) == 0),
            sum(n_i) + n_below_min + n_excluded_component == nrow(points))

  prep_info <- list(year = year, cellsize_m = cellsize_m, crs_epsg = 26916L,
                    min_count = min_count, largest_component = largest_component,
                    downloaded_records = n_downloaded,
                    invalid_coordinate_records = n_invalid_xy,
                    missing_category_records = n_missing_type,
                    outside_city_records = n_outside,
                    zero_count_city_cells = sum(audit$n_i == 0),
                    records_below_min_count = n_below_min,
                    records_excluded_by_component = n_excluded_component,
                    retained_records = sum(n_i), nodes = n, edges = n_edges,
                    components = igraph::components(g)$no,
                    isolated_nodes = sum(igraph::degree(g) == 0L))
  result <- list(Y = Y, n_i = n_i, n = n, D = d, categories = categories,
                 reference_category = tail(categories, 1L), g = g, D_G = D_G, W = W,
                 coords = coords, grid_sf = grid, city_boundary = city,
                 edges = data.frame(from = edge_index[, 1L], to = edge_index[, 2L]),
                 component = component, component_totals = component_totals,
                 cell_audit = audit, prep_info = prep_info,
                 download_info = snapshot[setdiff(names(snapshot), "data")],
                 session_info = utils::sessionInfo())
  result <- split_chicago_crimes(result, train_fraction = train_fraction,
                                 split_seed = split_seed)
  stem <- sprintf("chicago_%d_grid%sm_min%d_%s", as.integer(year),
                  format(cellsize_m, scientific = FALSE, trim = TRUE),
                  as.integer(min_count), if (largest_component) "largest" else "all")
  result$output_file <- file.path(output_dir, paste0(stem, ".rds"))
  saveRDS(result, result$output_file)
  node_table <- cbind(sf::st_drop_geometry(grid), as.data.frame(coords),
                      as.data.frame(Y, check.names = FALSE))
  utils::write.csv(node_table, file.path(output_dir, paste0(stem, "_counts.csv")),
                   row.names = FALSE)
  utils::write.csv(result$edges, file.path(output_dir, paste0(stem, "_edges.csv")),
                   row.names = FALSE)
  message("Prepared ", n, " nodes, ", n_edges, " edges, ", d, " categories and ",
          format(sum(n_i), big.mark = ","), " retained records.")
  message("Graph components: ", prep_info$components, "; isolated nodes: ",
          prep_info$isolated_nodes, ". See $prep_info for all exclusions.")
  message("Train/test split: ", format(result$split_info$train_crimes, big.mark = ","),
          " training crimes; ", format(result$split_info$test_crimes, big.mark = ","),
          " test crimes. See $Y_train, $Y_test and $split_info.")
  message("Saved: ", result$output_file)
  result
}

#' Restrict an existing Chicago object without changing its retained holdout.
#' @param chicago Prepared object containing counts, splits, and grid geometry.
#' @return A largest-component object with rebuilt graph matrices and an
#'   exclusion audit. The source cache and city boundary are unchanged.
chicago_largest_component <- function(chicago) {
  .chicago_require(c("sf", "Matrix", "igraph"))
  required <- c("Y", "Y_train", "Y_test", "grid_sf", "coords", "prep_info", "split_info")
  if (!is.list(chicago) || !all(required %in% names(chicago))) {
    stop("Supply a complete object from prepare_chicago_crimes().", call. = FALSE)
  }
  for (field in c("Y", "Y_train", "Y_test")) .chicago_check_counts(chicago[[field]])
  if (!identical(dimnames(chicago$Y), dimnames(chicago$Y_train)) ||
      !identical(dimnames(chicago$Y), dimnames(chicago$Y_test)) ||
      any(chicago$Y_train + chicago$Y_test != chicago$Y) ||
      !identical(rownames(chicago$Y), chicago$grid_sf$cell_key)) {
    stop("Chicago counts, splits, and grid rows must be aligned.", call. = FALSE)
  }
  graph <- .chicago_rook_graph(chicago$grid_sf)
  components <- igraph::components(graph)
  keep <- which(components$membership == which.max(components$csize))
  if (length(keep) == nrow(chicago$Y)) return(chicago)
  result <- chicago
  excluded <- setdiff(seq_len(nrow(chicago$Y)), keep)
  result$analysis_selection <- list(
    rule = "Largest connected component of the occupied rook grid",
    original_n_cells = nrow(chicago$Y), retained_n_cells = length(keep),
    excluded_cell_keys = rownames(chicago$Y)[excluded],
    excluded_records = sum(chicago$Y[excluded, , drop = FALSE])
  )
  result$excluded_cells <- sf::st_drop_geometry(chicago$grid_sf[excluded, ])
  for (field in c("Y", "Y_train", "Y_test")) {
    result[[field]] <- chicago[[field]][keep, , drop = FALSE]
  }
  result$n_i <- as.integer(rowSums(result$Y))
  result$n_i_train <- rowSums(result$Y_train)
  result$n_i_test <- rowSums(result$Y_test)
  result$n <- length(keep)
  result$grid_sf <- chicago$grid_sf[keep, ]
  result$grid_sf$node_id <- seq_len(result$n)
  result$coords <- chicago$coords[keep, , drop = FALSE]
  result$g <- .chicago_rook_graph(result$grid_sf)
  result$W <- igraph::as_adjacency_matrix(result$g, sparse = TRUE)
  edges <- igraph::as_edgelist(result$g, names = FALSE)
  m <- nrow(edges)
  result$D_G <- Matrix::sparseMatrix(
    i = rep(seq_len(m), 2L), j = c(edges[, 1L], edges[, 2L]),
    x = c(rep(1, m), rep(-1, m)), dims = c(m, result$n),
    dimnames = list(NULL, rownames(result$Y))
  )
  result$edges <- data.frame(from = edges[, 1L], to = edges[, 2L])
  result$component <- igraph::components(result$g)$membership
  result$component_totals <- rowsum(result$Y, result$component, reorder = FALSE)
  result$cell_audit$retained <- result$cell_audit$cell_key %in% rownames(result$Y)
  result$prep_info$largest_component <- TRUE
  result$prep_info$records_excluded_by_component <-
    result$prep_info$records_excluded_by_component + result$analysis_selection$excluded_records
  result$prep_info$retained_records <- sum(result$Y)
  result$prep_info$nodes <- result$n
  result$prep_info$edges <- m
  result$prep_info$components <- 1L
  result$prep_info$isolated_nodes <- sum(igraph::degree(result$g) == 0)
  result$split_info$train_fraction_actual <- sum(result$Y_train) / sum(result$Y)
  result$split_info$total_crimes <- sum(result$Y)
  result$split_info$train_crimes <- sum(result$Y_train)
  result$split_info$test_crimes <- sum(result$Y_test)
  result$split_info$total_cells <- result$n
  result$split_info$tested_cells <- sum(result$n_i_test > 0)
  result$split_info$singleton_cells_training_only <- sum(result$n_i == 1)
  result$source_output_file <- chicago$output_file
  result$output_file <- NULL  # Save this derived object to a new, explicit path.
  result
}


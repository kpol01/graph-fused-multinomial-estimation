# Run with: Rscript tests/validate_preparations.R
# All fixtures are synthetic. No live datasets are downloaded.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
search_before <- search()
for (script in c("prepare_chicago_crimes.R", "prepare_india_mother_tongue.R")) {
  source(file.path(project_dir, "R", script))
}
stopifnot(identical(search(), search_before))
for (package in c("sf", "Matrix", "igraph", "readxl", "ggplot2")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Install ", package)
}

passed <- 0L
check <- function(label, expr) {
  force(expr)
  passed <<- passed + 1L
  cat("PASS:", label, "\n")
}
expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
quiet <- function(expr) {
  suppressMessages(withCallingHandlers(force(expr), warning = function(w) {
    invokeRestart("muffleWarning")
  }))
}
rectangle <- function(x, y, width = 1000, height = 1000) {
  sf::st_polygon(list(rbind(c(x, y), c(x + width, y),
                            c(x + width, y + height), c(x, y + height), c(x, y))))
}
file_url <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  paste0(if (startsWith(path, "/")) "file://" else "file:///", path)
}

# When original uploads are available alongside this project, also compare
# their signatures and outputs. The downloadable tests do not need originals.
original_dir <- file.path(dirname(project_dir), "originals")
legacy <- NULL
scripts <- c("prepare_chicago_crimes.R", "prepare_india_mother_tongue.R")
if (all(file.exists(file.path(original_dir, scripts)))) {
  legacy <- new.env(parent = globalenv())
  for (script in scripts) source(file.path(original_dir, script), local = legacy)
}

check("standalone sourcing and public interface compatibility", {
  public <- c("split_chicago_crimes", "prepare_chicago_crimes",
              "split_mother_tongue_counts", "prepare_india_mother_tongue",
              "india_largest_component", "plot_india_mother_tongue",
              "score_india_mother_tongue")
  if (!is.null(legacy)) for (name in public) {
    stopifnot(identical(formals(get(name)), formals(get(name, envir = legacy))))
  }
  isolated <- new.env(parent = globalenv())
  source(file.path(project_dir, "R", scripts[1]), local = isolated)
  source(file.path(project_dir, "R", scripts[2]), local = isolated)
  stopifnot(all(vapply(public, exists, logical(1), envir = isolated, inherits = FALSE)),
            identical(search(), search_before))
})

Y <- rbind(c(1L, 0L, 0L), c(0L, 2L, 0L), c(0L, 0L, 3L), c(5L, 10L, 4L),
           c(500000L, 700000L, 200000L))
dimnames(Y) <- list(sprintf("%03d", 1:5), c("A", "B", "OTHER"))
check("Chicago splits retain exact assignments, integer counts, and cell rules", {
  chicago <- list(Y = Y, marker = "graph metadata stays fixed")
  for (fraction in c(0.2, 0.8)) {
    split <- split_chicago_crimes(chicago, fraction, 25L)
    stopifnot(all(split$Y_train + split$Y_test == Y),
              all(split$n_i_train >= 1), all(split$n_i_test[-1] >= 1),
              split$n_i_test[1] == 0, is.integer(split$Y_train),
              identical(split$marker, chicago$marker),
              identical(dimnames(split$Y_test), dimnames(Y)))
    if (!is.null(legacy)) {
      stopifnot(identical(split, legacy$split_chicago_crimes(chicago, fraction, 25L)))
    }
  }
})

check("India splits retain exact assignments and preserve zero-total rows", {
  counts <- rbind(Y, c(0L, 0L, 0L))
  rownames(counts)[6] <- "006"
  for (fraction in c(0.2, 0.8)) {
    split <- split_mother_tongue_counts(counts, fraction, 25L)
    stopifnot(all(split$Y_train + split$Y_test == counts),
              all(split$Y_train[6, ] == 0), all(split$Y_test[6, ] == 0),
              split$n_i_test[1] == 0, all(split$n_i_test[2:5] >= 1))
    if (!is.null(legacy)) {
      stopifnot(identical(split, legacy$split_mother_tongue_counts(counts, fraction, 25L)))
    }
  }
})

check("both splitters preserve the caller's RNG state", {
  set.seed(40)
  previous <- .Random.seed
  invisible(split_chicago_crimes(list(Y = Y)))
  invisible(split_mother_tongue_counts(Y))
  stopifnot(identical(.Random.seed, previous))
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(split_chicago_crimes(list(Y = Y)))
  invisible(split_mother_tongue_counts(Y))
  stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

check("rook adjacency preserves isolated cells and excludes diagonals", {
  grid <- data.frame(cell_key = c("0_0", "1_0", "2_1", "-3_-3"),
                      grid_x = c(0, 1, 2, -3), grid_y = c(0, 0, 1, -3))
  graph <- .chicago_rook_graph(grid)
  stopifnot(igraph::vcount(graph) == 4, igraph::ecount(graph) == 1,
            identical(igraph::V(graph)$name, grid$cell_key),
            all(igraph::degree(graph) == c(1, 1, 0, 0)))
})

india_crs <- "+proj=laea +lat_0=22 +lon_0=80 +datum=WGS84 +units=m +no_defs"
geometries <- sf::st_sfc(rectangle(0, 0), rectangle(1000, 1000),
                          rectangle(4000, 0), crs = india_crs)
district_sf <- sf::st_sf(district_code = c("001", "002", "003"), geometry = geometries)
check("queen graph, incidence, components, and Laplacian", {
  graph <- .india_graph(district_sf)
  stopifnot(nrow(graph$D_G) == 1, all(graph$degree == c(1, 1, 0)),
            all(graph$component == c(1, 1, 2)),
            all(rowSums(as.matrix(graph$D_G)) == 0),
            max(abs(as.matrix(graph$L - Matrix::crossprod(graph$D_G)))) < 1e-12)
  isolated <- .india_graph_from_neighbors(list(integer(), integer()), c("001", "002"))
  stopifnot(identical(dim(isolated$D_G), c(0L, 2L)), all(isolated$degree == 0))
  expect_error(.india_graph_from_neighbors(list(2L, integer()), c("001", "002")))
  expect_error(.india_graph_from_neighbors(list(1L), "001"))
})

fixture_dir <- file.path(project_dir, "tests", "fixtures")
cache <- tempfile("india-fixture-")
dir.create(cache)
for (state in c("00", "01")) {
  stopifnot(file.copy(file.path(fixture_dir, paste0("synthetic_C16_", state, ".xlsx")),
                     file.path(cache, paste0("C16_", state, ".xlsx"))))
}
shape <- sf::st_sf(ST_CEN_CD = rep(1L, 4), censuscode = c(1L, 2L, 3L, 0L),
                   geometry = c(geometries, sf::st_sfc(rectangle(6000, 0), crs = india_crs)))
sf::st_write(shape, file.path(cache, "2011_Dist.shp"), quiet = TRUE)

check("C-16 code padding, hierarchy, and population reconciliation", {
  table <- .india_read_c16(file.path(cache, "C16_01.xlsx"))
  stopifnot(identical(.india_code(c(1, 35), 2), c("01", "35")),
            sum(table$persons[table$district_code == "000" & table$is_language_total]) == 95,
            sum(table$persons[table$district_code != "000" & table$is_language_total]) == 95)
  expect_error(.india_code("1000", 3))
  expect_error(.india_code("-1", 3))
  expect_error(.india_read_c16(readxl::readxl_example("datasets.xlsx")))
})

india <- NULL
check("India preparation from real workbook-format and shapefile fixtures", {
  india <- quiet(prepare_india_mother_tongue(1, c(Hindi = "006240"), cache_dir = cache))
  stopifnot(identical(rownames(india$Y), c("001", "002", "003")),
            identical(colnames(india$Y), c("Hindi", "OTHER")),
            all(india$Y[, 1] == c(30, 10, 4)), all(india$n_i == c(50, 30, 15)),
            all(india$Y_train + india$Y_test == india$Y),
            identical(india$isolated_districts, "003"))
  if (!is.null(legacy)) {
    old <- quiet(legacy$prepare_india_mother_tongue(1, c(Hindi = "006240"), cache_dir = cache))
    for (field in c("Y", "Y_train", "Y_test", "W", "D_G", "component", "degree")) {
      stopifnot(identical(india[[field]], old[[field]]))
    }
  }
})

connected <- NULL
check("largest-component selection keeps splits and a stable exclusion audit", {
  stale <- india
  stale$component <- rep(99L, 3)
  connected <- quiet(india_largest_component(stale))
  repeated <- quiet(india_largest_component(connected))
  stopifnot(identical(rownames(connected$Y), c("001", "002")),
            identical(connected$Y_train, india$Y_train[1:2, , drop = FALSE]),
            identical(connected$boundary_sf, india$boundary_sf),
            identical(connected$analysis_selection, repeated$analysis_selection),
            identical(connected$excluded_districts, repeated$excluded_districts),
            all(connected$degree > 0))
  if (!is.null(legacy)) {
    old <- quiet(legacy$india_largest_component(india))
    stopifnot(identical(connected$Y, old$Y), identical(connected$D_G, old$D_G))
  }
})

check("probability alignment and independently calculated test scores", {
  object <- list(Y = matrix(c(5, 3, 3, 6), 2, 2, byrow = TRUE,
                             dimnames = list(c("001", "002"), c("Hindi", "OTHER"))))
  object$Y_test <- matrix(c(3, 1, 2, 4), 2, 2, byrow = TRUE, dimnames = dimnames(object$Y))
  P <- matrix(c(0.8, 0.2, 0.1, 0.9), 2, 2, byrow = TRUE, dimnames = dimnames(object$Y))
  scored <- score_india_mother_tongue(object, P)
  expected <- -(3 * log(0.8) + log(0.2) + 2 * log(0.1) + 4 * log(0.9)) / 10
  stopifnot(abs(scored$log_loss - expected) < 1e-12, abs(scored$brier - 0.484) < 1e-12,
            identical(scored, score_india_mother_tongue(object, P[2:1, 2:1])))
  unlabeled <- unname(P[, 2:1])
  stopifnot(identical(scored, score_india_mother_tongue(object, unlabeled,
                                                     category_order = c("OTHER", "Hindi"))))
  if (!is.null(legacy)) stopifnot(identical(scored, legacy$score_india_mother_tongue(object, P)))
  P[1, ] <- c(0, 1)
  stopifnot(is.infinite(score_india_mother_tongue(object, P)$log_loss))
  object$Y_test <- object$Y_test[, 2:1]
  expect_error(score_india_mother_tongue(object, P))
})

# Cache a small city with two rook-adjacent cells and a separate occupied cell.
chicago_cache <- tempfile("chicago-fixture-")
dir.create(chicago_cache)
city <- sf::st_sf(geometry = sf::st_sfc(rectangle(500000, 4600000, 2000, 1000),
                                       rectangle(506000, 4606000), crs = 26916))
sf::st_write(city, file.path(chicago_cache, "chicago_city_boundary.geojson"), quiet = TRUE)
types <- c("THEFT", "BATTERY", "CRIMINAL DAMAGE", "ASSAULT", "BURGLARY", "OTHER")
xy <- rbind(matrix(c(500500, 4600500), 7, 2, byrow = TRUE),
             matrix(c(501500, 4600500), 6, 2, byrow = TRUE), c(506500, 4606500),
             c(500500, 4600500), c(500500, 4600500),
             c(510500, 4610500), c(503500, 4603500))
points <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"), crs = 26916)
ll <- sf::st_coordinates(sf::st_transform(points, 4326))
raw <- data.frame(id = as.character(seq_len(nrow(xy))), date = "2024-06-01T12:00:00",
                   primary_type = c(types, "THEFT", types, "THEFT", "THEFT", "", "THEFT", "THEFT"),
                   latitude = as.character(ll[, 2]), longitude = as.character(ll[, 1]))
raw$longitude[15] <- "181"
api <- "https://data.cityofchicago.org/resource/ijzp-q8t2.csv"
period <- "date >= '2024-01-01T00:00:00' AND date < '2025-01-01T00:00:00'"
snapshot <- list(data = raw, api = api, period = period, expected_records = nrow(raw),
                  max_id = nrow(raw), downloaded_at = "Synthetic fixture")
saveRDS(snapshot, file.path(chicago_cache, "raw_chicago_2024.rds"))
chicago <- NULL
check("Chicago preparation and exclusions from a cached spatial fixture", {
  chicago <- quiet(prepare_chicago_crimes(output_dir = chicago_cache))
  stopifnot(chicago$n == 3, nrow(chicago$D_G) == 1, sum(chicago$Y) == 14,
            chicago$prep_info$invalid_coordinate_records == 1,
            chicago$prep_info$missing_category_records == 1,
            chicago$prep_info$outside_city_records == 2,
            all(chicago$Y_train + chicago$Y_test == chicago$Y),
            identical(rownames(chicago$Y), chicago$grid_sf$cell_key),
            identical(colnames(chicago$D_G), rownames(chicago$Y)),
            file.exists(chicago$output_file))
  if (!is.null(legacy)) {
    old <- quiet(legacy$prepare_chicago_crimes(output_dir = chicago_cache))
    for (field in c("Y", "Y_train", "Y_test", "D_G", "W", "prep_info", "split_info")) {
      stopifnot(identical(chicago[[field]], old[[field]]))
    }
  }
  subset <- quiet(prepare_chicago_crimes(output_dir = chicago_cache, largest_component = TRUE))
  filtered <- quiet(prepare_chicago_crimes(output_dir = chicago_cache, min_count = 2))
  stopifnot(subset$n == 2, subset$prep_info$records_excluded_by_component == 1,
            filtered$n == 2, filtered$prep_info$records_below_min_count == 1)
})

check("ID pagination and API schema checking with a mocked fetcher", {
  mock_env <- new.env(parent = globalenv())
  source(file.path(project_dir, "R", scripts[1]), local = mock_env)
  calls <- 0L
  mock_env$.chicago_fetch_table <- function(api, query, app_token = "") {
    calls <<- calls + 1L
    if (calls == 1L) return(data.frame(n = nrow(raw), max_id = nrow(raw)))
    if (calls == 2L) return(raw[1:9, ])
    stopifnot(grepl("id > 9", query$`$where`, fixed = TRUE))
    raw[10:18, ]
  }
  fresh <- tempfile("chicago-paging-")
  dir.create(fresh)
  file.copy(file.path(chicago_cache, "chicago_city_boundary.geojson"), fresh)
  result <- quiet(mock_env$prepare_chicago_crimes(output_dir = fresh))
  stopifnot(calls == 3L, identical(result$Y, chicago$Y))
  mock_env$.chicago_fetch_table <- function(...) data.frame(unexpected = "schema")
  expect_error(mock_env$prepare_chicago_crimes(output_dir = tempfile("invalid-paging-")))
})

check("local transfer helpers, empty-file recovery, and timeout restoration", {
  path <- tempfile(fileext = ".txt")
  writeLines("Synthetic test download", path)
  previous <- getOption("timeout")
  bytes <- .chicago_request(file_url(path))
  stopifnot(identical(bytes, readBin(path, "raw", n = file.info(path)$size)),
            identical(getOption("timeout"), previous))
  target <- tempfile()
  file.create(target)
  .india_download(file_url(path), target)
  stopifnot(file.info(target)$size > 0, identical(getOption("timeout"), previous),
            identical(.india_download("invalid://not-used", target), target))
  expect_error(.india_download(file_url(path), tempdir()))
})

check("map panels, probability scale, and row-order checks", {
  P <- connected$Y_train / connected$n_i_train
  plot <- plot_india_mother_tongue(connected, P, P, j = "Hindi")
  built <- ggplot2::ggplot_build(plot)
  stopifnot(nrow(built$layout$layout) == 3,
            identical(plot$scales$get_scales("fill")$limits, c(0, 1)))
  ggplot2::ggsave(file.path(tempdir(), "preparation_map_check.png"), plot,
                  width = 10, height = 4, dpi = 100)
  bad <- connected
  bad$district_sf <- bad$district_sf[2:1, ]
  expect_error(plot_india_mother_tongue(bad, P, P))
})

check("invalid arguments and malformed snapshots fail before preparation", {
  expect_error(split_chicago_crimes(list(Y = matrix(0, 1, 2))))
  expect_error(split_mother_tongue_counts(matrix(c(1.5, 2), 1)))
  expect_error(prepare_india_mother_tongue(split_seed = "bad", cache_dir = tempfile()))
  expect_error(prepare_chicago_crimes(train_fraction = 1, output_dir = tempfile()))
  expect_error(.chicago_check_snapshot(list(data = raw), period, api))
  broken <- snapshot
  broken$data$id[2] <- broken$data$id[1]
  expect_error(.chicago_check_snapshot(broken, period, api))
})

cat("\nAll", passed, "validation groups passed.\n")
cat(R.version.string, "\n")
for (package in c("sf", "Matrix", "igraph", "readxl", "ggplot2")) {
  cat(package, as.character(utils::packageVersion(package)), "\n")
}

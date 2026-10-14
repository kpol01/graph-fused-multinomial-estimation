# Census of India 2011 mother-tongue counts, district graph, maps, and scores.
# Source this file to define functions; preparation downloads/cache reads are
# explicit. Categories are published mother tongues, with OTHER as reference.
# Queen contiguity uses polygon contacts; isolated districts are preserved.
# See PREPARATION_GUIDE.md for sources, setup, and examples.

.india_check_scalar <- function(x, name, lower, upper = Inf, integer = FALSE) {
  valid <- is.numeric(x) && !is.complex(x) && length(x) == 1L && is.finite(x) &&
    x >= lower && x <= upper && (!integer || x == floor(x))
  if (!valid) {
    stop(name, " must be a finite ", if (integer) "integer" else "number",
         " in [", lower, ", ", upper, "].", call. = FALSE)
  }
  invisible(x)
}

.india_check_split <- function(train_fraction, split_seed) {
  .india_check_scalar(train_fraction, "train_fraction", lower = 0, upper = 1)
  if (train_fraction == 0 || train_fraction == 1) {
    stop("Use 0 < train_fraction < 1.", call. = FALSE)
  }
  .india_check_scalar(split_seed, "split_seed", lower = 0,
                      upper = .Machine$integer.max, integer = TRUE)
}

.india_check_counts <- function(Y, label = "Y") {
  if (!is.matrix(Y) || !is.numeric(Y) || is.complex(Y) || nrow(Y) < 1L ||
      ncol(Y) < 2L || any(!is.finite(Y)) || any(Y < 0 | Y != floor(Y))) {
    stop(label, " must be a nonempty numeric matrix of nonnegative integer ",
         "counts, with at least two categories.", call. = FALSE)
  }
  invisible(Y)
}

.india_check_named_counts <- function(india) {
  if (!is.list(india) || is.null(india$Y)) {
    stop("Supply an object containing the raw count matrix india$Y.", call. = FALSE)
  }
  .india_check_counts(india$Y, "india$Y")
  ids <- rownames(india$Y)
  categories <- colnames(india$Y)
  if (is.null(ids) || anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids)) ||
      is.null(categories) || anyNA(categories) || anyDuplicated(categories) ||
      any(!nzchar(categories))) {
    stop("india$Y must have unique, nonempty district and category names.",
         call. = FALSE)
  }
  invisible(india)
}

.india_check_count_order <- function(Y, reference, label) {
  .india_check_counts(Y, label)
  if (!identical(dim(Y), dim(reference)) ||
      (!is.null(rownames(Y)) && !identical(rownames(Y), rownames(reference))) ||
      (!is.null(colnames(Y)) && !identical(colnames(Y), colnames(reference)))) {
    stop(label, " must have the same dimensions and named row/category order ",
         "as india$Y.", call. = FALSE)
  }
}

.india_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install these R packages first: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
}

.india_download <- function(url, path) {
  if (!is.character(url) || length(url) != 1L || is.na(url) || !nzchar(url) ||
      !is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("url and path must be single nonempty character strings.", call. = FALSE)
  }
  if (dir.exists(path)) stop("Download path is a directory: ", path, call. = FALSE)
  if (file.exists(path) && file.info(path)$size > 0) return(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(dirname(path))) {
    stop("Could not create the cache directory: ", dirname(path), call. = FALSE)
  }
  temporary <- tempfile(pattern = paste0(basename(path), "."),
                        tmpdir = dirname(path), fileext = ".part")
  old_options <- options(timeout = max(300, getOption("timeout", 60)))
  on.exit({ options(old_options); unlink(temporary) }, add = TRUE)
  status <- tryCatch(
    utils::download.file(url, temporary, mode = "wb", method = "libcurl", quiet = TRUE),
    error = function(e) stop("Could not download ", url, "\n", conditionMessage(e),
                             call. = FALSE)
  )
  if (status != 0 || !file.exists(temporary) || file.info(temporary)$size == 0) {
    stop("Download failed: ", url, call. = FALSE)
  }
  # A zero-byte cached file is not reusable; remove it before the final rename.
  if (file.exists(path) && file.info(path)$size == 0) unlink(path)
  if (!file.rename(temporary, path)) {
    stop("Could not save download: ", path, call. = FALSE)
  }
  path
}

.india_code <- function(x, width) {
  .india_check_scalar(width, "width", lower = 1, upper = 9, integer = TRUE)
  x <- trimws(as.character(x))
  if (anyNA(x) || any(!grepl("^[0-9]+$", x))) {
    stop("Census codes must contain digits only.", call. = FALSE)
  }
  value <- suppressWarnings(as.numeric(x))
  if (any(!is.finite(value)) || any(value > 10^width - 1)) {
    stop("A Census code exceeds its ", width, "-digit field.", call. = FALSE)
  }
  sprintf(paste0("%0", width, "d"), as.integer(value))
}

.india_read_c16 <- function(path) {
  raw <- readxl::read_excel(path, sheet = 1, col_names = FALSE,
                            col_types = "text", .name_repair = "minimal")
  if (nrow(raw) < 1L || ncol(raw) < 16L ||
      !isTRUE(grepl("C-16", raw[[1L]][1L], fixed = TRUE))) {
    stop("Unexpected C-16 workbook structure: ", path, call. = FALSE)
  }
  # Keep state/district totals. Never sum subdistrict rows into their parents.
  keep <- which(trimws(raw[[1]]) == "C0116" &
                  trimws(raw[[4]]) %in% c("0", "00000"))
  if (!length(keep)) stop("No state/district C-16 rows in ", path)
  raw <- raw[keep, , drop = FALSE]
  numeric_count <- function(j) {
    x <- suppressWarnings(as.numeric(gsub(",", "", raw[[j]], fixed = TRUE)))
    if (anyNA(x) || any(!is.finite(x) | x < 0 | x != floor(x))) {
      stop("Missing or invalid population count in ", path, ", column ", j)
    }
    x
  }
  persons <- numeric_count(8)
  if (any(persons != numeric_count(9) + numeric_count(10)) ||
      any(persons != numeric_count(11) + numeric_count(14))) {
    stop("Persons/sex or total/rural/urban counts do not reconcile: ", path)
  }
  d <- data.frame(
    state_code = .india_code(raw[[2]], 2),
    district_code = .india_code(raw[[3]], 3),
    area_name = trimws(raw[[5]]),
    mother_tongue_code = .india_code(raw[[6]], 6),
    mother_tongue = trimws(raw[[7]]), persons = persons,
    stringsAsFactors = FALSE
  )
  d$is_language_total <- grepl("000$", d$mother_tongue_code)
  d$language_code <- paste0(substr(d$mother_tongue_code, 1, 3), "000")
  if (anyDuplicated(d[c("state_code", "district_code", "mother_tongue_code")])) {
    stop("Duplicate area/category rows in ", path)
  }
  # Omitted mother-tongue rows can be zero only if published child rows fully
  # reconcile with each parent total. The standalone 124000 OTHERS has no children.
  children <- d[!d$is_language_total, ]
  child_sums <- if (nrow(children) > 0L) {
    stats::aggregate(persons ~ state_code + district_code + language_code,
                     children, sum)
  } else {
    children[, c("state_code", "district_code", "language_code", "persons")]
  }
  parents <- d[d$is_language_total,
               c("state_code", "district_code", "language_code", "persons")]
  checked <- merge(child_sums, parents,
                   by = c("state_code", "district_code", "language_code"),
                   all = TRUE, suffixes = c("_children", "_parent"))
  exempt <- checked$language_code == "124000" & is.na(checked$persons_children)
  bad <- !exempt & (is.na(checked$persons_children) |
                      is.na(checked$persons_parent) |
                      checked$persons_children != checked$persons_parent)
  if (any(bad)) stop("Mother-tongue hierarchy does not reconcile: ", path)
  d
}

#' Split aggregate counts without replacement within each district.
#' @param Y Nonnegative integer district-by-category matrix, including raw zeros.
#' @param train_fraction Requested training fraction, strictly between 0 and 1.
#' @param split_seed Nonnegative integer seed, local to this split.
#' @return Y_train, Y_test, their row totals, and split_info. Larger districts
#'   have at least one observation in each set; singletons stay in training.
#'   Zero-total input rows are preserved. The caller's RNG state is restored.
split_mother_tongue_counts <- function(Y, train_fraction = 0.8, split_seed = 123L) {
  Y <- as.matrix(Y)
  .india_check_counts(Y)
  .india_check_split(train_fraction, split_seed)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(split_seed)
  n <- rowSums(Y)
  n_train <- ifelse(n <= 1, n, pmax(1, pmin(n - 1, floor(train_fraction * n))))
  Y_train <- matrix(0, nrow(Y), ncol(Y), dimnames = dimnames(Y))
  # Sequential multivariate hypergeometric draw; no expansion into individuals.
  for (i in seq_len(nrow(Y))) {
    remaining_total <- n[i]
    remaining_train <- n_train[i]
    for (j in seq_len(ncol(Y) - 1L)) {
      draw <- if (remaining_train == 0 || Y[i, j] == 0) 0 else
        stats::rhyper(1, m = Y[i, j], n = remaining_total - Y[i, j],
                      k = remaining_train)
      Y_train[i, j] <- draw
      remaining_total <- remaining_total - Y[i, j]
      remaining_train <- remaining_train - draw
    }
    Y_train[i, ncol(Y)] <- remaining_train
  }
  Y_test <- Y - Y_train
  stopifnot(all(Y_train + Y_test == Y), all(Y_train >= 0), all(Y_test >= 0),
            all(rowSums(Y_train) == n_train),
            identical(dimnames(Y_train), dimnames(Y)),
            identical(dimnames(Y_test), dimnames(Y)))
  list(Y_train = Y_train, Y_test = Y_test,
       n_i_train = rowSums(Y_train), n_i_test = rowSums(Y_test),
       split_info = list(train_fraction = train_fraction, split_seed = split_seed,
                         scheme = "Without replacement within each district",
                         singleton_rule = "A singleton goes to training"))
}

.india_graph <- function(district_sf) {
  # Queen contiguity: touching at an edge or a point. No distance-based bridges.
  nb <- sf::st_touches(district_sf)
  .india_graph_from_neighbors(nb, district_sf$district_code)
}

.india_graph_from_neighbors <- function(nb, ids) {
  n <- length(nb)
  if (!is.list(nb) || n < 1L || length(ids) != n || anyNA(ids) ||
      anyDuplicated(ids) || any(!nzchar(as.character(ids)))) {
    stop("Supply a nonempty neighbor list and one unique district ID per node.",
         call. = FALSE)
  }
  for (i in seq_len(n)) {
    neighbors <- nb[[i]]
    if (!is.numeric(neighbors) || is.complex(neighbors) ||
        any(!is.finite(neighbors)) ||
        any(neighbors < 1 | neighbors > n | neighbors != floor(neighbors)) ||
        anyDuplicated(neighbors) || i %in% neighbors) {
      stop("Neighbor lists must contain distinct, valid indices without self links.",
           call. = FALSE)
    }
  }
  reciprocal <- vapply(seq_len(n), function(i) {
    all(vapply(nb[[i]], function(j) i %in% nb[[j]], logical(1L)))
  }, logical(1L))
  if (!all(reciprocal)) {
    stop("Neighbor lists must describe a symmetric undirected graph.", call. = FALSE)
  }
  ids <- as.character(ids)
  if (sum(lengths(nb)) == 0L) {
    edges <- matrix(integer(), nrow = 0L, ncol = 2L)
  } else {
    ij <- cbind(rep(seq_len(n), lengths(nb)), unlist(nb, use.names = FALSE))
    edges <- ij[ij[, 1] < ij[, 2], , drop = FALSE]
  }
  colnames(edges) <- c("from", "to")
  m <- nrow(edges)
  W <- Matrix::sparseMatrix(
    i = c(edges[, 1], edges[, 2]), j = c(edges[, 2], edges[, 1]),
    x = rep(1, 2 * m), dims = c(n, n), dimnames = list(ids, ids))
  D_G <- Matrix::sparseMatrix(
    i = rep(seq_len(m), 2), j = c(edges[, 1], edges[, 2]),
    x = c(rep(1, m), rep(-1, m)), dims = c(m, n),
    dimnames = list(paste(ids[edges[, 1]], ids[edges[, 2]], sep = "--"), ids))
  degree <- as.numeric(Matrix::rowSums(W))
  component <- integer(n)
  nc <- 0L
  for (root in seq_len(n)) {
    if (component[root] != 0L) next
    nc <- nc + 1L
    component[root] <- nc
    queue <- integer(n)
    queue[1L] <- root
    head <- 1L
    tail <- 1L
    while (head <= tail) {
      i <- queue[head]
      head <- head + 1L
      fresh <- nb[[i]][component[nb[[i]]] == 0L]
      component[fresh] <- nc
      if (length(fresh) > 0L) {
        queue[tail + seq_along(fresh)] <- fresh
        tail <- tail + length(fresh)
      }
    }
  }
  names(degree) <- names(component) <- ids
  stopifnot(all(Matrix::diag(W) == 0), Matrix::isSymmetric(W))
  list(W = W, D_G = D_G, degree = degree, component = component,
       edges = edges, neighbors = nb,
       L = Matrix::Diagonal(x = degree) - W)
}

#' Prepare 2011 district counts and a queen-contiguity graph.
#' @param state_codes 2011 Census state/union-territory codes, 01 through 35.
#' @param categories Named vector of distinct leaf mother-tongue codes.
#'   OTHER is added last; language totals ending in 000 are not categories.
#' @param train_fraction,split_seed Arguments for the within-district holdout.
#' @param cache_dir Directory for pinned C-16 workbooks and district boundaries.
#' @return Full/training/test counts, district geometries, graph matrices,
#'   a category catalog, source metadata, and component audits. The full
#'   requested district set is retained; no bridges or pseudocounts are added.
prepare_india_mother_tongue <- function(
    state_codes = 1:35,
    categories = c(Hindi = "006240", Bengali = "002007", Marathi = "013071",
                   Telugu = "021046", Tamil = "020027"),
    train_fraction = 0.8, split_seed = 123L,
    cache_dir = "india_mother_tongue_cache") {
  .india_require(c("readxl", "sf", "Matrix"))
  .india_check_split(train_fraction, split_seed)
  if (!is.character(cache_dir) || length(cache_dir) != 1L || is.na(cache_dir) ||
      !nzchar(trimws(cache_dir))) {
    stop("cache_dir must be a single nonempty character string.", call. = FALSE)
  }
  states <- sort(unique(.india_code(state_codes, 2)))
  if (!length(states) || any(!states %in% sprintf("%02d", 1:35))) {
    stop("state_codes must be 2011 Census codes from 01 through 35.")
  }
  if (is.null(names(categories)) || !length(categories) ||
      anyNA(names(categories)) || any(!nzchar(trimws(names(categories)))) ||
      anyDuplicated(names(categories)) || any(names(categories) == "OTHER")) {
    stop("categories must be uniquely named codes; OTHER is added automatically.")
  }
  labels <- names(categories)
  codes <- .india_code(unname(categories), 6)
  if (anyDuplicated(codes) || any(grepl("000$", codes))) {
    stop("Choose distinct leaf mother-tongue codes, not language totals ending 000.")
  }
  census_revision <- "d70c2ae4c18f603e165bf3adb4b1c06566848e4a"
  boundary_revision <- "b3fbbde595310b397a55d718e0958ce249a4fa1f"
  census_base <- paste0("https://raw.githubusercontent.com/",
                        "sunil-dhaka/census-language-analysis/", census_revision, "/datasets/C-16/")
  boundary_base <- paste0("https://raw.githubusercontent.com/datameet/maps/",
                          boundary_revision, "/Districts/Census_2011/")
  workbook <- function(st) .india_download(
    paste0(census_base, st, ".xlsx"), file.path(cache_dir, paste0("C16_", st, ".xlsx")))
  national <- .india_read_c16(workbook("00"))
  country <- national[national$state_code == "00", ]
  catalog <- country[!country$is_language_total,
                     c("mother_tongue_code", "mother_tongue", "language_code")]
  parent_names <- country[country$is_language_total,
                          c("mother_tongue_code", "mother_tongue")]
  parent_index <- match(catalog$language_code, parent_names$mother_tongue_code)
  catalog$language_group <- sub("^[0-9]+[[:space:]]+", "",
                                parent_names$mother_tongue[parent_index])
  if (any(!codes %in% catalog$mother_tongue_code)) {
    stop("Unknown mother-tongue code(s): ", paste(setdiff(codes, catalog$mother_tongue_code),
                                                  collapse = ", "))
  }
  district_rows <- lapply(states, function(st) {
    message("Reading Census C-16: state ", st)
    d <- .india_read_c16(workbook(st))
    if (any(d$state_code != st)) stop("Wrong state in workbook ", st)
    d <- d[d$district_code != "000", ]
    expected <- sum(national$persons[national$state_code == st &
                                       national$is_language_total])
    if (!nrow(d) || sum(d$persons[d$is_language_total]) != expected) {
      stop("District counts do not reconcile with the national workbook for state ", st)
    }
    d
  })
  d <- do.call(rbind, district_rows)
  districts <- unique(d[c("state_code", "district_code", "area_name")])
  districts <- districts[order(districts$district_code), ]
  rownames(districts) <- NULL
  names(districts)[names(districts) == "area_name"] <- "district_name"
  ids <- districts$district_code
  if (anyDuplicated(ids)) stop("District codes are not unique.")
  totals <- tapply(d$persons[d$is_language_total],
                   d$district_code[d$is_language_total], sum)
  Y <- matrix(0, length(ids), length(codes) + 1L,
              dimnames = list(ids, c(labels, "OTHER")))
  for (j in seq_along(codes)) {
    z <- d[d$mother_tongue_code == codes[j], ]
    Y[match(z$district_code, ids), j] <- z$persons
  }
  Y[, ncol(Y)] <- as.numeric(totals[ids]) - rowSums(Y)
  stopifnot(all(Y >= 0), all(rowSums(Y) == totals[ids]))
  districts$n_i <- rowSums(Y)
  for (ext in c("shp", "shx", "dbf", "prj")) {
    .india_download(paste0(boundary_base, "2011_Dist.", ext),
                    file.path(cache_dir, paste0("2011_Dist.", ext)))
  }
  raw_shapes <- sf::st_read(file.path(cache_dir, "2011_Dist.shp"), quiet = TRUE)
  if (!all(c("ST_CEN_CD", "censuscode") %in% names(raw_shapes))) {
    stop("Unexpected boundary fields; require ST_CEN_CD and censuscode.")
  }
  # censuscode is the all-India district code. DT_CEN_CD is NOT that code.
  shape_ids <- .india_code(raw_shapes$censuscode, 3)
  shape_states <- .india_code(raw_shapes$ST_CEN_CD, 2)
  available <- shape_ids[shape_states %in% states & shape_ids != "000"]
  if (anyDuplicated(available) || !setequal(ids, available)) {
    stop("Census and boundary district codes do not match.\nMissing boundaries: ",
         paste(setdiff(ids, available), collapse = ", "),
         "\nMissing counts: ", paste(setdiff(available, ids), collapse = ", "))
  }
  index <- match(ids, shape_ids)
  if (any(shape_states[index] != districts$state_code)) stop("State-code mismatch.")
  shapes <- raw_shapes[index, ]
  # A metre-based India-centred equal-area projection; no change to counts.
  map_crs <- "+proj=laea +lat_0=22 +lon_0=80 +datum=WGS84 +units=m +no_defs"
  district_sf <- sf::st_sf(districts, geometry = sf::st_geometry(shapes))
  district_sf <- sf::st_make_valid(sf::st_transform(district_sf, map_crs))
  if (any(sf::st_is_empty(district_sf)) || any(sf::st_dimension(district_sf) != 2)) {
    stop("An empty or non-polygon district geometry was encountered.")
  }
  graph <- .india_graph(district_sf)
  split <- split_mother_tongue_counts(Y, train_fraction, split_seed)
  component_totals <- rowsum(split$Y_train, graph$component, reorder = TRUE)
  audit <- data.frame(
    component = seq_len(max(graph$component)),
    n_districts = as.integer(table(graph$component)),
    training_population = rowSums(component_totals),
    missing_categories = apply(component_totals, 1, function(x)
      paste(colnames(Y)[x == 0], collapse = ", ")),
    row.names = NULL
  )
  if (any(component_totals == 0)) warning(
    "A graph component has zero training counts in a category. Inspect ",
    "component_summary/component_totals; an unconstrained logit fit may diverge. ",
    "Raw zeros have been preserved.", call. = FALSE)
  result <- c(list(
    Y = Y, n_i = rowSums(Y), categories = colnames(Y),
    category_codes = setNames(c(codes, NA_character_), colnames(Y)),
    reference_category = "OTHER", mother_tongue_catalog = catalog,
    districts = districts, district_sf = district_sf,
    boundary_sf = sf::st_sf(geometry = sf::st_union(sf::st_geometry(district_sf))),
    component_totals = component_totals, component_summary = audit,
    isolated_districts = ids[graph$degree == 0],
    source_info = list(census_year = 2011L, state_codes = states,
                       census_base = census_base, boundary_base = boundary_base,
                       boundary_license = "DataMeet, CC BY 2.5 India",
                       national_c16_population = sum(country$persons[country$is_language_total]),
                       unmatched_boundary_features = sf::st_drop_geometry(
                         raw_shapes[shape_ids == "000", ]))
  ), graph, split)
  result$train <- list(Y = result$Y_train, n_i = result$n_i_train)
  result$test <- list(Y = result$Y_test, n_i = result$n_i_test)
  message("Prepared ", nrow(Y), " districts, ", ncol(Y), " categories, ",
          nrow(graph$D_G), " edges, ", nrow(audit), " graph components.")
  result
}

#' Restrict a prepared object to its largest connected component.
#' @param india Object returned by prepare_india_mother_tongue().
#' @return A new object with aligned subset counts, geometry, and rebuilt
#'   graph matrices. Splits stay fixed. The full boundary and exclusion audit
#'   are retained, including on repeated calls. Ties use the first component.
india_largest_component <- function(india) {
  # Restrict to the largest component of the existing adjacency graph. All
  # category counts and train/test assignments for retained districts stay fixed.
  # The original india object is unchanged: assign the returned object separately.
  .india_require(c("sf", "Matrix"))
  required <- c("Y", "Y_train", "Y_test", "W", "districts", "district_sf", "boundary_sf")
  if (!is.list(india) || !all(required %in% names(india))) {
    stop("Supply an object from prepare_india_mother_tongue().", call. = FALSE)
  }
  .india_check_named_counts(india)
  .india_check_count_order(india$Y_train, india$Y, "india$Y_train")
  .india_check_count_order(india$Y_test, india$Y, "india$Y_test")
  if (any(india$Y_train + india$Y_test != india$Y)) {
    stop("Training and test counts must partition india$Y.", call. = FALSE)
  }
  n <- nrow(india$Y)
  ids <- rownames(india$Y)
  W <- as.matrix(india$W)
  if (!identical(dim(W), c(n, n)) || !is.numeric(W) || is.complex(W) || anyNA(W) ||
      any(!W %in% c(0, 1)) || !isSymmetric(unname(W)) || any(diag(W) != 0)) {
    stop("Expected a symmetric binary adjacency matrix with a zero diagonal.")
  }
  if (length(ids) != n || anyDuplicated(ids) ||
      !identical(as.character(india$districts$district_code), ids) ||
      !identical(as.character(india$district_sf$district_code), ids) ||
      !identical(dim(india$Y_train), dim(india$Y)) ||
      !identical(dim(india$Y_test), dim(india$Y))) {
    stop("The district, geometry and count orders do not match.")
  }
  # Recompute component membership from W, not from potentially stale metadata.
  if (!is.null(rownames(W)) || !is.null(colnames(W))) {
    if (!setequal(rownames(W), ids) || !setequal(colnames(W), ids)) {
      stop("W's district names do not match india$Y.")
    }
    W <- W[ids, ids, drop = FALSE]
  }
  graph <- .india_graph_from_neighbors(lapply(seq_len(n), function(i) which(W[i, ] > 0)), ids)
  sizes <- tabulate(graph$component)
  keep <- which(graph$component == which.max(sizes))
  if (length(keep) < 2L) {
    stop("No connected group has at least two districts. A different, justified graph is needed.")
  }
  result <- india
  result$excluded_districts <- india$districts[-keep, , drop = FALSE]
  result$analysis_selection <- list(
    rule = "Largest connected component of the supplied adjacency graph",
    original_n_districts = n, retained_n_districts = length(keep),
    retained_rows = keep, excluded_rows = setdiff(seq_len(n), keep),
    excluded_district_codes = ids[-keep]
  )
  # Re-running the fix on an already restricted object must preserve the original
  # exclusion record; the full map boundary is also retained across repeated calls.
  if (length(keep) == n && !is.null(india$analysis_selection)) {
    result$analysis_selection <- india$analysis_selection
    result$excluded_districts <- india$excluded_districts
  }
  for (name in c("Y", "Y_train", "Y_test")) {
    result[[name]] <- india[[name]][keep, , drop = FALSE]
  }
  result$n_i <- rowSums(result$Y)
  result$n_i_train <- rowSums(result$Y_train)
  result$n_i_test <- rowSums(result$Y_test)
  result$districts <- india$districts[keep, , drop = FALSE]
  result$district_sf <- india$district_sf[keep, , drop = FALSE]
  W <- W[keep, keep, drop = FALSE]
  subgraph <- .india_graph_from_neighbors(
    lapply(seq_along(keep), function(i) which(W[i, ] > 0)), ids[keep])
  for (name in names(subgraph)) result[[name]] <- subgraph[[name]]
  result$train <- list(Y = result$Y_train, n_i = result$n_i_train)
  result$test <- list(Y = result$Y_test, n_i = result$n_i_test)
  result$component_totals <- rowsum(result$Y_train, result$component)
  result$component_summary <- data.frame(
    component = 1L, n_districts = length(keep),
    training_population = sum(result$Y_train),
    missing_categories = paste(colnames(result$Y)[colSums(result$Y_train) == 0], collapse = ", ")
  )
  result$isolated_districts <- character()
  stopifnot(all(result$degree > 0), max(result$component) == 1L,
            all(result$Y_train + result$Y_test == result$Y))
  if (any(result$component_totals == 0)) warning(
    "A category has zero training counts in the retained component. ",
    "Inspect component_totals before an unconstrained logit fit.", call. = FALSE)
  message("Retained ", length(keep), " of ", n, " districts in one connected component; ",
          n - length(keep), " excluded. Use this returned object for both fits, maps and scores.")
  result
}

.india_probabilities <- function(P, india, label, category_order = NULL) {
  .india_check_named_counts(india)
  P <- as.matrix(P)
  if (!is.numeric(P) || is.complex(P) || !identical(dim(P), dim(india$Y)) ||
      anyNA(P) || any(!is.finite(P) | P < 0 | P > 1) ||
      any(abs(rowSums(P) - 1) > 1e-6)) {
    stop(label, " must be a districts x categories probability matrix with rows summing to 1.")
  }
  expected <- colnames(india$Y)
  # An explicit order describes the meaning of the existing columns, before
  # reordering. It does not merely rename an already reordered matrix to Y order.
  if (!is.null(category_order)) {
    category_order <- as.character(category_order)
    if (length(category_order) != ncol(P) || anyNA(category_order) ||
        anyDuplicated(category_order) || !setequal(category_order, expected)) {
      stop(label, ": category_order must name every india$Y category exactly once, ",
           "in this fitted matrix's current column order. Expected categories: ",
           paste(expected, collapse = ", "), call. = FALSE)
    }
    colnames(P) <- category_order
  }
  # Named outputs are aligned explicitly. Unnamed outputs must already follow Y.
  if (!is.null(rownames(P))) {
    if (anyNA(rownames(P)) || anyDuplicated(rownames(P)) ||
        !setequal(rownames(P), rownames(india$Y))) {
      stop(label, " row names must be india$Y district codes; or supply an unnamed matrix in Y order.")
    }
    P <- P[rownames(india$Y), , drop = FALSE]
  }
  if (!is.null(colnames(P))) {
    if (anyNA(colnames(P)) || anyDuplicated(colnames(P)) ||
        !setequal(colnames(P), expected)) {
      stop(label, " category names do not match india$Y.\n",
           "Fitted column names: ", paste(colnames(P), collapse = ", "), "\n",
           "Expected categories: ", paste(expected, collapse = ", "), "\n",
           "Specify the fitted columns' actual category names using category_order ",
           "(mcar_category_order for the MCAR matrix in plot_india_mother_tongue). ",
           "Use colnames(india$Y) only if the output retained that column order.",
           call. = FALSE)
    }
    P <- P[, expected, drop = FALSE]
  }
  P
}

#' Compare observed and fitted probabilities for a single mother-tongue category.
#' @param india Prepared object matching both model fits.
#' @param P_hat,P_hat_MCAR District-by-category fitted probability matrices.
#' @param j Category name or column index in india$Y.
#' @param observed Whether observations use test counts or full Census counts.
#' @param category_order,mcar_category_order Actual existing column order for
#'   each fit, when its labels are absent or generic. Named outputs are aligned.
#' @return A ggplot with three panels and a common probability scale [0, 1].
#'   Areas outside the analysis subset or with no observed counts are grey.
plot_india_mother_tongue <- function(india, P_hat, P_hat_MCAR, j = 1,
                                     observed = c("test", "full"),
                                     category_order = NULL, mcar_category_order = NULL) {
  .india_require(c("sf", "ggplot2"))
  observed <- match.arg(observed)
  P_hat <- .india_probabilities(P_hat, india, "P_hat", category_order)
  P_hat_MCAR <- .india_probabilities(P_hat_MCAR, india, "P_hat_MCAR", mcar_category_order)
  if (!inherits(india$district_sf, "sf") || !inherits(india$boundary_sf, "sf") ||
      !identical(as.character(india$district_sf$district_code), rownames(india$Y))) {
    stop("Map geometries must follow the district order of india$Y.", call. = FALSE)
  }
  if (is.character(j)) j <- match(j, colnames(india$Y))
  if (!is.numeric(j) || is.complex(j) || length(j) != 1L || !is.finite(j) ||
      j < 1 || j > ncol(india$Y) || j != floor(j)) {
    stop("j must be a category name or column number in india$Y.")
  }
  Y_observed <- if (observed == "test") india$Y_test else india$Y
  .india_check_count_order(Y_observed, india$Y, "Observed counts")
  n_observed <- rowSums(Y_observed)
  p_observed <- ifelse(n_observed > 0, Y_observed[, j] / n_observed, NA_real_)
  values <- list(p_observed, P_hat[, j], P_hat_MCAR[, j])
  panels <- c(if (observed == "test") "Observed (test)" else "Observed (full Census)",
              "Graph-fused multinomial", "MCAR")
  map_data <- do.call(rbind, lapply(seq_along(values), function(k) {
    layer <- india$district_sf
    layer$probability <- values[[k]]
    layer$probability[n_observed == 0] <- NA_real_
    layer$panel <- factor(panels[k], levels = panels)
    layer
  }))
  # District polygons themselves supply the fill: no square grid or bounding-box
  # clipping. All panels use the same geometries, extent, category and colour scale.
  ggplot2::ggplot(map_data) +
    ggplot2::geom_sf(data = india$boundary_sf, fill = "grey90", color = NA) +
    ggplot2::geom_sf(ggplot2::aes(fill = probability), color = NA) +
    ggplot2::geom_sf(data = india$boundary_sf, fill = NA, color = "grey35", linewidth = 0.2) +
    ggplot2::facet_wrap(~panel, nrow = 1) +
    ggplot2::scale_fill_viridis_c(option = "C", limits = c(0, 1),
                                  na.value = "grey90", name = "Probability") +
    ggplot2::coord_sf(datum = NA) +
    ggplot2::labs(title = paste("India, Census 2011:", colnames(india$Y)[j]),
                  subtitle = "District mother-tongue proportions and fitted probabilities",
                  caption = paste("2011 boundaries: DataMeet (CC BY 2.5 India).",
                                   "Grey: excluded areas or no observed records.")) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   strip.text = ggplot2::element_text(face = "bold"))
}

#' Score fitted probabilities against the within-district test counts.
#' @param india Prepared object with Y_test in the same district/category order.
#' @param P_hat Fitted district-by-category probability matrix.
#' @param category_order Actual existing column order, if labels need declaration.
#' @return A one-row data frame containing test_population, mean individual
#'   log_loss, and multicategory brier score. Both scores weight districts by
#'   their test populations. A zero prediction for an observed event gives
#'   infinite log loss; brier is not MSE between district proportions.
score_india_mother_tongue <- function(india, P_hat, category_order = NULL) {
  P <- .india_probabilities(P_hat, india, "P_hat", category_order)
  Y <- india$Y_test
  .india_check_count_order(Y, india$Y, "india$Y_test")
  n <- rowSums(Y)
  if (sum(n) == 0) stop("The test set is empty.")
  positive <- Y > 0
  # A predicted zero for an observed event correctly gives infinite log loss.
  log_loss <- -sum(Y[positive] * log(P[positive])) / sum(n)
  # Average individual multicategory Brier score, not MSE of district proportions.
  brier <- sum(n * rowSums(P^2) - 2 * rowSums(Y * P) + n) / sum(n)
  data.frame(test_population = sum(n), log_loss = log_loss, brier = brier)
}

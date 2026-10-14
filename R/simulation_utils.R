# Truth generation and multinomial simulation for the GFM experiments.

.gfm_check_grid <- function(g, D, nrow_grid, ncol_grid) {
  .gfm_check_scalar(D, "D", lower = 2, integer = TRUE)
  .gfm_check_scalar(nrow_grid, "nrow_grid", lower = 1, integer = TRUE)
  .gfm_check_scalar(ncol_grid, "ncol_grid", lower = 1, integer = TRUE)
  if (!inherits(g, "igraph") ||
      igraph::vcount(g) != nrow_grid * ncol_grid) {
    stop("g must have nrow_grid * ncol_grid vertices.", call. = FALSE)
  }
  row <- igraph::V(g)$row
  col <- igraph::V(g)$col
  if (length(row) != igraph::vcount(g) || length(col) != igraph::vcount(g) ||
      any(!is.finite(row)) || any(!is.finite(col)) ||
      any(row < 1 | row > nrow_grid | row != floor(row)) ||
      any(col < 1 | col > ncol_grid | col != floor(col))) {
    stop("g must have valid integer row and col vertex attributes.", call. = FALSE)
  }
  invisible(g)
}

.gfm_check_X_parameters <- function(D, p_high, p_low, half_width) {
  .gfm_check_scalar(D, "D", lower = 2, integer = TRUE)
  .gfm_check_scalar(p_high, "p_high", lower = 0, upper = 1)
  .gfm_check_scalar(p_low, "p_low", lower = 0, upper = 1)
  if (p_high %in% c(0, 1) || p_low %in% c(0, 1)) {
    stop("p_high and p_low must lie strictly between zero and one.", call. = FALSE)
  }
  .gfm_check_scalar(half_width, "half_width", lower = 0)
}

# Category one has probability p1; all remaining categories share the rest.
.gfm_truth_from_p1 <- function(p1, D) {
  q <- D - 1L
  P_true <- cbind(p1, matrix((1 - p1) / q, length(p1), q))
  beta_true <- matrix(0, length(p1), q)
  beta_true[, 1L] <- stats::qlogis(p1) + log(q)
  list(beta_true = beta_true, P_true = P_true)
}

#' Generate a quadrant signal with a spatial trend and Gaussian perturbations.
#' @param g Grid graph with row and col vertex attributes.
#' @param D Number of categories, at least two.
#' @param nrow_grid,ncol_grid Grid dimensions matching g.
#' @return A list containing beta_true and quadrant. The original base
#'   logit patterns are truncated or recycled when D differs from four.
generate_true_beta <- function(g, D, nrow_grid, ncol_grid) {
  .gfm_check_grid(g, D, nrow_grid, ncol_grid)
  n <- igraph::vcount(g)
  q <- D - 1L
  row <- igraph::V(g)$row
  col <- igraph::V(g)$col
  quadrant <- ifelse(row <= nrow_grid / 2,
                     ifelse(col <= ncol_grid / 2, 1L, 2L),
                     ifelse(col <= ncol_grid / 2, 3L, 4L))
  base_levels <- list(c(1.5, -1.0, 0.5), c(-1.0, 1.5, 0.0),
                      c(0.0, 0.0, 1.5), c(1.0, 1.0, -1.5))
  beta_true <- matrix(0, n, q)
  for (region in seq_len(4L)) {
    index <- which(quadrant == region)
    if (length(index) > 0L) {
      beta_true[index, ] <- matrix(rep(base_levels[[region]], length.out = q),
                                   length(index), q, byrow = TRUE)
    }
  }
  standardize <- function(x) {
    spread <- stats::sd(x)
    if (length(x) < 2L || !is.finite(spread) || spread == 0) return(rep(0, length(x)))
    (x - mean(x)) / spread
  }
  trend <- outer(standardize(row) + standardize(col), rep(1, q)) * 0.15
  beta_true <- beta_true + trend + matrix(stats::rnorm(n * q, sd = 0.05), n, q)
  list(beta_true = beta_true, quadrant = quadrant)
}

#' Generate an X-shaped probability signal on a rectangular grid.
#' @param g Grid graph with row and col vertex attributes.
#' @param D Number of categories, at least two.
#' @param nrow_grid,ncol_grid Grid dimensions matching g.
#' @param p_high,p_low Category-one probabilities outside and inside the X.
#' @param half_width Nonnegative perpendicular band half-width in [0, 1] units.
#' @return beta_true, P_true, and the logical vertex indicator on_X.
#'   A grid axis of length one is placed at its normalized midpoint, 0.5.
generate_true_beta_X <- function(g, D, nrow_grid, ncol_grid,
                                 p_high = 0.70, p_low = 0.08, half_width = 0.06) {
  .gfm_check_grid(g, D, nrow_grid, ncol_grid)
  normalize_axis <- function(x, size) {
    if (size == 1L) rep(0.5, length(x)) else (x - 1) / (size - 1)
  }
  coords <- cbind(normalize_axis(igraph::V(g)$row, nrow_grid),
                  normalize_axis(igraph::V(g)$col, ncol_grid))
  generate_true_beta_X_coords(coords, D, p_high, p_low, half_width)
}

#' Generate the same X-shaped signal at supplied spatial locations.
#' @param coords Finite matrix with at least two columns; the first two
#'   columns are used as supplied, typically on the unit square.
#' @param D Number of categories, at least two.
#' @param p_high,p_low Category-one probabilities outside and inside the X.
#' @param half_width Nonnegative perpendicular band half-width.
#' @return beta_true, P_true, and on_X. Coordinates are not rescaled.
generate_true_beta_X_coords <- function(coords, D, p_high = 0.70, p_low = 0.08,
                                        half_width = 0.06) {
  .gfm_check_matrix(coords, "coords", min_cols = 2L)
  .gfm_check_X_parameters(D, p_high, p_low, half_width)
  x <- coords[, 1L]
  y <- coords[, 2L]
  distance_main <- abs(x - y) / sqrt(2)
  distance_anti <- abs(x + y - 1) / sqrt(2)
  on_X <- pmin(distance_main, distance_anti) < half_width
  truth <- .gfm_truth_from_p1(ifelse(on_X, p_low, p_high), D)
  c(truth, list(on_X = on_X))
}

#' Generate a block signal with node-level perturbations on the logit scale.
#' @param g Graph with positive integer block vertex labels.
#' @param D Number of categories, at least two.
#' @param p1_by_block Category-one center probabilities indexed by block label.
#' @param within_block_sd Nonnegative SD of the Gaussian logit perturbation.
#' @return beta_true, P_true, and block. Other categories have equal probability.
generate_true_beta_blocks <- function(g, D, p1_by_block, within_block_sd = 0.35) {
  .gfm_check_scalar(D, "D", lower = 2, integer = TRUE)
  .gfm_check_scalar(within_block_sd, "within_block_sd", lower = 0)
  if (!inherits(g, "igraph") || igraph::vcount(g) < 1L) {
    stop("g must be a nonempty igraph object.", call. = FALSE)
  }
  block <- igraph::V(g)$block
  if (!is.numeric(block) || length(block) != igraph::vcount(g) ||
      any(!is.finite(block)) || any(block < 1 | block != floor(block))) {
    stop("g must have positive integer block vertex labels.", call. = FALSE)
  }
  if (!is.numeric(p1_by_block) || length(p1_by_block) < max(block) ||
      any(!is.finite(p1_by_block)) || any(p1_by_block <= 0 | p1_by_block >= 1)) {
    stop("p1_by_block must supply a probability in (0, 1) for every block label.",
         call. = FALSE)
  }
  n <- igraph::vcount(g)
  q <- D - 1L
  logit_node <- stats::qlogis(p1_by_block[block]) +
    stats::rnorm(n, mean = 0, sd = within_block_sd)
  beta_true <- matrix(0, n, q)
  beta_true[, 1L] <- logit_node + log(q)
  list(beta_true = beta_true, P_true = softmax_probs(beta_true), block = block)
}

#' Simulate independent multinomial count vectors at each node.
#' @param P Finite n by D probability matrix with nonnegative entries and
#'   row sums within 1e-8 of one; D >= 2.
#' @param n_i Length-n vector of nonnegative integer trial counts.
#' @return An n by D count matrix with row sums n_i.
simulate_data <- function(P, n_i) {
  .gfm_check_matrix(P, "P", min_cols = 2L)
  if (any(P < 0) || any(abs(rowSums(P) - 1) > 1e-8)) {
    stop("P must have nonnegative entries and rows summing to one.", call. = FALSE)
  }
  if (!is.numeric(n_i) || length(n_i) != nrow(P) || any(!is.finite(n_i)) ||
      any(n_i < 0 | n_i != floor(n_i) | n_i > .Machine$integer.max)) {
    stop("n_i must contain one nonnegative integer per node, within the ",
         "range supported by rmultinom().", call. = FALSE)
  }
  Y <- matrix(0, nrow(P), ncol(P), dimnames = dimnames(P))
  for (i in seq_len(nrow(P))) {
    Y[i, ] <- as.numeric(stats::rmultinom(1L, size = n_i[i], prob = P[i, ]))
  }
  Y
}

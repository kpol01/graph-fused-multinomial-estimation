# Graph construction and sparse incidence-matrix utilities.

#' Build an undirected rectangular grid with rook adjacency.
#' @param nrow_grid,ncol_grid Positive integer grid dimensions.
#' @return An igraph object with vertex attributes row and col.
build_grid_graph <- function(nrow_grid, ncol_grid) {
  .gfm_check_scalar(nrow_grid, "nrow_grid", lower = 1, integer = TRUE)
  .gfm_check_scalar(ncol_grid, "ncol_grid", lower = 1, integer = TRUE)
  g <- igraph::make_lattice(dimvector = c(nrow_grid, ncol_grid), nei = 1,
                            directed = FALSE, circular = FALSE)
  igraph::V(g)$row <- rep(seq_len(nrow_grid), times = ncol_grid)
  igraph::V(g)$col <- rep(seq_len(ncol_grid), each = nrow_grid)
  g
}

#' Build a rectangular grid with rook and diagonal (king) adjacency.
#' @param nrow_grid,ncol_grid Positive integer grid dimensions.
#' @return An igraph object; a one-dimensional grid has only rook edges.
build_king_grid_graph <- function(nrow_grid, ncol_grid) {
  g <- build_grid_graph(nrow_grid, ncol_grid)
  if (nrow_grid == 1L || ncol_grid == 1L) return(g)
  vertex_id <- function(row, col) (col - 1L) * nrow_grid + row
  cells <- expand.grid(row = seq_len(nrow_grid - 1L),
                       col = seq_len(ncol_grid - 1L))
  from <- c(vertex_id(cells$row, cells$col),
            vertex_id(cells$row + 1L, cells$col))
  to <- c(vertex_id(cells$row + 1L, cells$col + 1L),
          vertex_id(cells$row, cells$col + 1L))
  igraph::add_edges(g, as.vector(rbind(from, to)))
}

#' Construct an oriented, sparse incidence matrix.
#' @param g An igraph object. Edge weights are not used.
#' @return An |E| by |V| sparse matrix; each stored edge has one row,
#'   with +1 at its first endpoint and -1 at its second endpoint.
#'   The orientation does not affect an absolute-difference penalty.
build_incidence_matrix <- function(g) {
  if (!inherits(g, "igraph")) stop("g must be an igraph object.", call. = FALSE)
  edges <- igraph::as_edgelist(g, names = FALSE)
  n_edges <- nrow(edges)
  Matrix::sparseMatrix(
    i = rep(seq_len(n_edges), 2L),
    j = c(edges[, 1L], edges[, 2L]),
    x = c(rep(1, n_edges), rep(-1, n_edges)),
    dims = c(n_edges, igraph::vcount(g))
  )
}

# A guaranteed bound: ||A||_2^2 <= ||A||_1 * ||A||_infinity.
# For an unweighted incidence matrix this is 2 * maximum vertex degree.
.gfm_operator_bound <- function(D_G) {
  if (nrow(D_G) == 0L || ncol(D_G) == 0L) return(0)
  bound <- as.numeric(Matrix::norm(D_G, "1") * Matrix::norm(D_G, "I"))
  if (!is.finite(bound)) {
    stop("The operator-norm bound exceeds the numeric range.", call. = FALSE)
  }
  bound
}

#' Estimate the squared operator norm by power iteration.
#' @param D_G Numeric or Matrix incidence matrix.
#' @param n_iter Positive integer iteration limit.
#' @param tol Positive relative convergence tolerance.
#' @return A numerical estimate of ||D_G||_2^2. This is a diagnostic,
#'   not a guaranteed upper bound for selecting solver step sizes.
#'   Initialization is reproducible and preserves the caller's RNG state.
op_norm_sq_DG <- function(D_G, n_iter = 200, tol = 1e-10) {
  .gfm_check_operator(D_G)
  .gfm_check_scalar(n_iter, "n_iter", lower = 1, integer = TRUE)
  .gfm_check_scalar(tol, "tol", lower = .Machine$double.xmin)
  if (.gfm_operator_bound(D_G) == 0) return(0)
  v <- .gfm_with_seed(1L, rnorm(ncol(D_G)))
  v <- v / sqrt(sum(v^2))
  previous <- 0
  for (iteration in seq_len(n_iter)) {
    w <- as.numeric(Matrix::crossprod(D_G, D_G %*% v))
    estimate <- sqrt(sum(w^2))
    if (!is.finite(estimate)) {
      stop("Power iteration exceeded the numeric range.", call. = FALSE)
    }
    if (estimate == 0) return(0)
    v <- w / estimate
    if (abs(estimate - previous) <= tol * max(1, estimate)) break
    previous <- estimate
  }
  sum(v * as.numeric(Matrix::crossprod(D_G, D_G %*% v)))
}

#' Build the undirected symmetrized union k-nearest-neighbor graph.
#' @param coords Finite n by d coordinate matrix, with n >= 2 and d >= 1.
#' @param k Integer in {1, ..., n - 1}.
#' @return A simple igraph object with all n vertices and attributes x and y.
#'   An edge is included when either endpoint selects the other. Queries
#'   explicitly remove self matches, including for coincident locations.
#'   Distance ties follow FNN's deterministic query order.
build_knn_graph <- function(coords, k) {
  .gfm_check_matrix(coords, "coords", min_rows = 2L)
  n <- nrow(coords)
  .gfm_check_scalar(k, "k", lower = 1, upper = n - 1L, integer = TRUE)
  if (k == n - 1L) {
    g <- igraph::make_full_graph(n, directed = FALSE, loops = FALSE)
  } else {
    neighbors <- FNN::get.knnx(coords, coords, k = k + 1L)$nn.index
    selected <- matrix(0L, n, k)
    for (i in seq_len(n)) {
      candidates <- neighbors[i, ]
      selected[i, ] <- head(candidates[candidates != i], k)
    }
    from <- rep(seq_len(n), each = k)
    to <- as.vector(t(selected))
    edges <- unique(cbind(pmin(from, to), pmax(from, to)))
    g <- igraph::make_empty_graph(n, directed = FALSE)
    g <- igraph::add_edges(g, as.vector(t(edges)))
  }
  igraph::V(g)$x <- coords[, 1L]
  igraph::V(g)$y <- if (ncol(coords) >= 2L) coords[, 2L] else rep(0, n)
  g
}

#' Sample an SBM and repair isolated vertices for CAR/MCAR comparisons.
#' @param block_sizes Positive integer vector with total size at least two.
#' @param p_in,p_out Within-block and between-block probabilities in [0, 1].
#' @return A simple undirected igraph object with a block vertex attribute
#'   and positive degree at every vertex. Repairs prefer the same block;
#'   singleton blocks use another vertex. The result is a repaired SBM.
build_sbm_graph <- function(block_sizes, p_in = 0.03, p_out = 0.001) {
  if (!is.numeric(block_sizes) || length(block_sizes) == 0L ||
      any(!is.finite(block_sizes)) || any(block_sizes < 1) ||
      any(block_sizes != floor(block_sizes)) || sum(block_sizes) < 2) {
    stop("block_sizes must be positive integers with total size >= 2.",
         call. = FALSE)
  }
  .gfm_check_scalar(p_in, "p_in", lower = 0, upper = 1)
  .gfm_check_scalar(p_out, "p_out", lower = 0, upper = 1)
  n_blocks <- length(block_sizes)
  preference <- matrix(p_out, n_blocks, n_blocks)
  diag(preference) <- p_in
  g <- igraph::sample_sbm(sum(block_sizes), pref.matrix = preference,
                          block.sizes = block_sizes, directed = FALSE,
                          loops = FALSE)
  block <- rep(seq_len(n_blocks), times = block_sizes)
  igraph::V(g)$block <- block
  isolated <- which(igraph::degree(g) == 0)
  if (length(isolated) > 0L) {
    partners <- vapply(isolated, function(vertex) {
      candidates <- setdiff(which(block == block[vertex]), vertex)
      if (length(candidates) == 0L) {
        candidates <- setdiff(seq_len(igraph::vcount(g)), vertex)
      }
      # sample(x, 1) treats a length-one numeric x as the sequence 1:x.
      candidates[sample.int(length(candidates), 1L)]
    }, integer(1L))
    g <- igraph::add_edges(g, as.vector(rbind(isolated, partners)))
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  }
  stopifnot(all(igraph::degree(g) > 0))
  g
}

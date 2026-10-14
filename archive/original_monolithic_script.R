## ============================================================
##  Graph-Fused Multinomial Logistic Regression
##  Primal-Dual (Chambolle-Pock type) Algorithm - Simulation Study
## ============================================================
##
## Model:
##   Y_i ~ Multinomial(n_i, p_i),  i = 1,...,n
##   p_ij = exp(beta_ij) / (1 + sum_k exp(beta_ik)),  j = 1,...,D-1
##   p_iD = 1 / (1 + sum_k exp(beta_ik))
##
## Penalty: total-variation (graph fusion) on beta across edges of G,
##          separately for each of the D-1 logit coordinates.
##
## ============================================================

## ---- 0. Packages -------------------------------------------------------

pkgs <- c("igraph", "Matrix", "ggplot2", "reshape2", "FNN")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(1)

## ============================================================
## 1. Build the graph G  (2D grid, e.g. spatial sensor network)
## ============================================================

build_grid_graph <- function(nrow_grid, ncol_grid) {
  g <- igraph::make_lattice(dimvector = c(nrow_grid, ncol_grid), nei = 1)
  V(g)$row <- rep(1:nrow_grid, times = ncol_grid)
  V(g)$col <- rep(1:ncol_grid, each  = nrow_grid)
  g
}

## Same grid, but also connected diagonally ("king" adjacency: 8 neighbors
## instead of 4). Useful when the true spatial pattern has diagonal edges
## (e.g. an X shape): a fusion penalty on a rook (4-neighbor) graph can only
## propagate smoothness along rows/columns, so it approximates a diagonal
## line with a staircase; adding diagonal edges lets the penalty follow the
## diagonal directly.
build_king_grid_graph <- function(nrow_grid, ncol_grid) {
  g <- build_grid_graph(nrow_grid, ncol_grid)
  vid <- function(r, c) (c - 1) * nrow_grid + r   # matches make_lattice's column-major vertex order
  rc <- expand.grid(r = 1:(nrow_grid - 1), c = 1:(ncol_grid - 1))
  from1 <- vid(rc$r, rc$c);        to1 <- vid(rc$r + 1, rc$c + 1)   # "\" diagonal
  from2 <- vid(rc$r + 1, rc$c);    to2 <- vid(rc$r, rc$c + 1)       # "/" diagonal
  edge_vec <- as.vector(rbind(c(from1, from2), c(to1, to2)))
  igraph::add_edges(g, edge_vec)
}

## Signed edge-vertex incidence matrix D_G  (|E| x n), sparse.
## Row e = (i, i'):  D_G[e, i] = 1, D_G[e, i'] = -1
build_incidence_matrix <- function(g) {
  el <- igraph::as_edgelist(g, names = FALSE)
  n  <- igraph::vcount(g)
  nE <- nrow(el)
  Matrix::sparseMatrix(
    i = rep(1:nE, 2),
    j = c(el[, 1], el[, 2]),
    x = c(rep(1, nE), rep(-1, nE)),
    dims = c(nE, n)
  )
}

## Operator norm^2 of D_G via power iteration on the graph Laplacian
## L = D_G' D_G  (avoids full eigendecomposition for large graphs)
op_norm_sq_DG <- function(D_G, n_iter = 200, tol = 1e-10) {
  n <- ncol(D_G)
  v <- rnorm(n)
  v <- v / sqrt(sum(v^2))
  lambda_old <- 0
  for (k in 1:n_iter) {
    w <- Matrix::crossprod(D_G, D_G %*% v)   # L %*% v
    lambda_new <- sqrt(sum(as.numeric(w)^2))
    v <- as.numeric(w) / lambda_new
    if (abs(lambda_new - lambda_old) < tol) break
    lambda_old <- lambda_new
  }
  lambda_new  # this is the top eigenvalue of L = ||D_G||_op^2
}

## ============================================================
## 2. Generate ground truth beta with a blocky (piecewise-constant)
##    spatial structure on the grid, for D categories
## ============================================================

generate_true_beta <- function(g, D, nrow_grid, ncol_grid) {
  n <- igraph::vcount(g)
  row_v <- V(g)$row
  col_v <- V(g)$col
  
  # Split the grid into 4 quadrants; each quadrant gets its own
  # constant logit vector (length D-1) -> piecewise constant truth,
  # ideal setting for a total-variation / graph-fusion penalty.
  quadrant <- ifelse(row_v <= nrow_grid / 2,
                     ifelse(col_v <= ncol_grid / 2, 1, 2),
                     ifelse(col_v <= ncol_grid / 2, 3, 4))
  
  base_levels <- list(
    c( 1.5, -1.0, 0.5),
    c(-1.0,  1.5, 0.0),
    c( 0.0,  0.0, 1.5),
    c( 1.0,  1.0, -1.5)
  )
  # if D-1 != 3, adjust / recycle
  Dm1 <- D - 1
  base_levels <- lapply(base_levels, function(v) {
    if (length(v) >= Dm1) v[1:Dm1] else rep(v, length.out = Dm1)
  })
  
  beta_true <- matrix(0, n, Dm1)
  for (q in 1:4) beta_true[quadrant == q, ] <- matrix(base_levels[[q]], sum(quadrant == q), Dm1, byrow = TRUE)
  
  # small smooth spatial trend + noise so it's not perfectly piecewise-constant
  trend <- outer(scale(row_v)[, 1] + scale(col_v)[, 1], rep(1, Dm1)) * 0.15
  beta_true <- beta_true + trend + matrix(rnorm(n * Dm1, sd = 0.05), n, Dm1)
  
  list(beta_true = beta_true, quadrant = quadrant)
}

softmax_probs <- function(beta) {
  # beta: n x (D-1) matrix of logits (category D is reference)
  expb <- exp(beta)
  denom <- 1 + rowSums(expb)
  P <- expb / denom
  pD <- 1 / denom
  cbind(P, pD)
}

## Alternative ground truth: category 1's rate is LOW along the two
## diagonal arms of an X and HIGH in the background (matches the uploaded
## figure "True metabolite 1 rates"). The remaining probability mass is
## split evenly across the other D-1 categories, which are left spatially
## flat so the X pattern in category 1 is the only structure being tested.
generate_true_beta_X <- function(g, D, nrow_grid, ncol_grid,
                                 p_high = 0.70, p_low = 0.08, half_width = 0.06) {
  n <- igraph::vcount(g)
  row_v <- V(g)$row
  col_v <- V(g)$col
  
  r <- (row_v - 1) / (nrow_grid - 1)          # normalize both axes to [0, 1]
  c <- (col_v - 1) / (ncol_grid - 1)
  
  dist_main <- abs(r - c) / sqrt(2)           # perpendicular distance to the "\" diagonal
  dist_anti <- abs(r + c - 1) / sqrt(2)       # perpendicular distance to the "/" diagonal
  on_X <- pmin(dist_main, dist_anti) < half_width
  
  p1 <- ifelse(on_X, p_low, p_high)           # metabolite-1 rate
  p_rest <- (1 - p1) / (D - 1)
  P_true <- cbind(p1, matrix(p_rest, n, D - 1))
  
  Dm1 <- D - 1
  beta_true <- log(P_true[, 1:Dm1, drop = FALSE] / P_true[, D])
  
  list(beta_true = beta_true, P_true = P_true, on_X = on_X)
}

## ============================================================
## 3. Simulate multinomial count data on the graph
## ============================================================

simulate_data <- function(P, n_i) {
  n <- nrow(P)
  D <- ncol(P)
  Y <- matrix(0, n, D)
  for (i in 1:n) Y[i, ] <- as.numeric(rmultinom(1, size = n_i[i], prob = P[i, ]))
  Y
}

## ============================================================
## 4. Primal-Dual Algorithm  (Algorithm 1 in the paper)
## ============================================================

proj_box <- function(x, lambda) pmin(pmax(x, -lambda), lambda)

fit_graph_fused_multinomial <- function(Y, D_G, lambda,
                                        tau = NULL, sigma = NULL,
                                        theta = 1, a = 0.5,
                                        T_max = 2000, tol = 1e-6,
                                        beta_init = NULL,
                                        verbose = TRUE) {
  n <- nrow(Y)
  D <- ncol(Y)
  Dm1 <- D - 1
  n_i <- rowSums(Y)
  N <- n_i                       # store as vector; multiply rowwise
  Y_minus_D <- Y[, 1:Dm1, drop = FALSE]
  
  ## --- step-size selection satisfying 1/tau - sigma*||D_G||_op^2 > L_f/2 ---
  op2 <- op_norm_sq_DG(D_G)
  Lf  <- max(n_i)                # conservative Lipschitz bound for multinomial logit gradient
  if (is.null(sigma)) sigma <- 1 / (2 * op2)
  if (is.null(tau)) {
    rhs <- Lf / 2 + sigma * op2
    tau <- 0.9 / rhs              # 0.9 safety margin below the feasibility boundary
  }
  stopifnot(1 / tau - sigma * op2 > Lf / 2)
  
  ## --- initialization ---
  a_const <- a
  if (is.null(beta_init)) {
    beta <- log((Y[, 1:Dm1, drop = FALSE] + a_const) / (Y[, D] + a_const))
  } else {
    beta <- beta_init            # warm start, e.g. from a neighboring lambda / CV fold
  }
  U <- Matrix::Matrix(0, nrow(D_G), Dm1, sparse = FALSE)
  beta_bar <- beta
  
  resid_hist <- numeric(T_max)
  t_stop <- T_max
  
  for (t in 1:T_max) {
    beta_old <- beta
    U_old <- U
    
    ## Dual update
    U <- proj_box(as.matrix(U_old + sigma * (D_G %*% beta_bar)), lambda)
    
    ## Multinomial probabilities at current primal iterate
    expb <- exp(beta_old)
    denom <- 1 + rowSums(expb)
    P <- expb / denom
    
    ## Gradient of the (negative) log-likelihood part
    Gmat <- N * P - Y_minus_D          # elementwise row-scaling: N is a vector -> N*P broadcasts by row
    
    ## Primal update
    beta <- beta_old - tau * (Gmat + as.matrix(Matrix::crossprod(D_G, U)))
    
    ## Extrapolation
    beta_bar <- beta + theta * (beta - beta_old)
    
    ## Relative residual for stopping
    r_beta <- norm(beta - beta_old, "F") / max(1, norm(beta_old, "F"))
    r_U    <- norm(U - U_old, "F") / max(1, norm(U_old, "F"))
    resid_hist[t] <- max(r_beta, r_U)
    
    if (verbose && (t %% 100 == 0)) {
      cat(sprintf("iter %5d   residual = %.3e\n", t, resid_hist[t]))
    }
    
    if (resid_hist[t] < tol) { t_stop <- t; break }
  }
  
  P_hat <- softmax_probs(beta)
  
  list(beta_hat = beta, P_hat = P_hat, U = U,
       n_iter = t_stop, resid_hist = resid_hist[1:t_stop],
       tau = tau, sigma = sigma, Lf = Lf, op2 = op2)
}

## ============================================================
## 5. Cross-validation for lambda (within-node count splitting)
## ============================================================
##
## Rationale: the fusion penalty requires the FULL graph G to be present in
## every fold (deleting a node would delete its edges and break the penalty,
## and there is no notion of "predicting" a held-out node since beta_i is a
## free parameter at each node, not a function of covariates). Instead we
## split each node's total count n_i into K independent replicate count
## vectors and hold out one at a time, keeping the graph intact throughout.
##
## Key fact used: if Y_i ~ Multinomial(n_i, p_i), and each of the n_i trials
## is assigned i.i.d. and uniformly to one of K folds, the resulting fold
## count vectors Y_i^(1),...,Y_i^(K) are themselves INDEPENDENT multinomial
## draws from the same p_i. This can be simulated cell-by-cell:
##   (Y_ij^(1),...,Y_ij^(K)) ~ Multinomial(Y_ij, (1/K,...,1/K))
## independently for each (i,j); summing across folds reproduces Y_ij.

## Split an n x D count matrix into K fold-count arrays via the cellwise
## multinomial-thinning scheme described above.
split_counts_Kfold <- function(Y, K, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(Y); D <- ncol(Y)
  folds <- array(0, dim = c(n, D, K))
  for (i in 1:n) {
    for (j in 1:D) {
      yij <- Y[i, j]
      if (yij > 0) {
        folds[i, j, ] <- as.numeric(rmultinom(1, size = yij, prob = rep(1 / K, K)))
      }
    }
  }
  folds  # folds[, , k] is the n x D count matrix held out in fold k
}

## K-fold CV over a grid of lambda values.
## Returns, for each lambda: the average held-out multinomial deviance
## (per held-out trial, so it's comparable across folds/lambdas), its
## across-fold standard error, and the resulting lambda.min / lambda.1se.
cv_graph_fused_multinomial <- function(Y, D_G, lambda_grid, K = 5, seed = 1,
                                       T_max = 3000, tol = 1e-6,
                                       warm_start = TRUE, verbose = TRUE) {
  n <- nrow(Y); D <- ncol(Y)
  folds <- split_counts_Kfold(Y, K, seed = seed)
  
  cv_mat  <- matrix(NA_real_, K, length(lambda_grid))   # held-out NLL per fold/lambda
  n_test  <- numeric(K)                                  # held-out trial count per fold
  
  # sort lambda increasingly so warm starts move smoothly along the path
  ord <- order(lambda_grid)
  lambda_sorted <- lambda_grid[ord]
  
  for (k in 1:K) {
    Y_test  <- folds[, , k]
    Y_train <- Y - Y_test
    n_test[k] <- sum(Y_test)
    
    beta_warm <- NULL
    for (li in seq_along(lambda_sorted)) {
      lam <- lambda_sorted[li]
      fit <- fit_graph_fused_multinomial(Y_train, D_G, lambda = lam,
                                         T_max = T_max, tol = tol,
                                         beta_init = beta_warm, verbose = FALSE)
      if (warm_start) beta_warm <- fit$beta_hat
      
      P_hat <- pmax(fit$P_hat, 1e-10)                 # guard against log(0)
      nll <- -sum(Y_test * log(P_hat))
      cv_mat[k, match(lam, lambda_grid)] <- nll
      
      if (verbose) {
        cat(sprintf("fold %d/%d | lambda = %7.4f | held-out NLL = %.2f\n",
                    k, K, lam, nll))
      }
    }
  }
  
  cv_fold_rate <- sweep(cv_mat, 1, n_test, "/")          # per-trial deviance, per fold
  cv_mean <- colMeans(cv_fold_rate)
  cv_se   <- apply(cv_fold_rate, 2, sd) / sqrt(K)
  
  min_idx    <- which.min(cv_mean)
  lambda_min <- lambda_grid[min_idx]
  thresh     <- cv_mean[min_idx] + cv_se[min_idx]
  lambda_1se <- max(lambda_grid[cv_mean <= thresh])      # largest lambda within 1 SE
  
  list(lambda_grid = lambda_grid, cv_mean = cv_mean, cv_se = cv_se,
       lambda_min = lambda_min, lambda_1se = lambda_1se,
       cv_mat = cv_mat, n_test = n_test)
}

## ============================================================
## 6. Run the simulation: true "metabolite 1" rate forms an X
## ============================================================
##
## This mirrors the uploaded figure: metabolite 1 has a LOW rate along
## the two diagonal arms of an X and a HIGH rate everywhere else. Since
## the arms are diagonal, this is a good stress test for the ROOK
## (4-neighbor, axis-aligned only) grid graph used above: a graph-fusion
## penalty on a rook graph can only propagate smoothness along rows and
## columns, so a sharp diagonal edge tends to come out as a "staircase"
## rather than a clean diagonal line. We also fit the same data on a
## KING graph (rook graph + diagonal edges) to show how adding diagonal
## connectivity mitigates the staircasing. Which graph to use is a real
## modeling choice based on the geometry you expect in your application,
## not something the algorithm decides on its own.

## --- grid size -------------------------------------------------------
## The uploaded figure looks like roughly a 90x90 grid; we default to
## 50x50 here for speed. Set nrow_grid/ncol_grid <- 90 to match it
## exactly (slower: ~8100 nodes, ~16000+ edges).
nrow_grid <- 50
ncol_grid <- 50
D <- 4                                     # number of categories (e.g. metabolites)

g_rook <- build_grid_graph(nrow_grid, ncol_grid)
g_king <- build_king_grid_graph(nrow_grid, ncol_grid)
n <- igraph::vcount(g_rook)
D_G_rook <- build_incidence_matrix(g_rook)
D_G_king <- build_incidence_matrix(g_king)

## --- true X-shaped spatial pattern for metabolite (category) 1 -------
truth <- generate_true_beta_X(g_rook, D, nrow_grid, ncol_grid,
                              p_high = 0.70, p_low = 0.08, half_width = 0.06)
beta_true <- truth$beta_true
P_true <- truth$P_true

## --- sample sizes per node and observed counts ------------------------
## Larger totals than a blocky-quadrant truth would need: recovering a
## thin diagonal feature requires more signal-to-noise per node than
## recovering large flat regions does.
n_i <- sample(150:300, n, replace = TRUE)
Y <- simulate_data(P_true, n_i)

## --- oracle sweep over lambda (rook graph), for the MSE-vs-lambda plot
lambda_grid <- c(0, 0.05, 0.1, 0.25, 0.5, 1, 2, 4)

mse_P_grid <- numeric(length(lambda_grid))
beta_warm <- NULL
for (k in seq_along(lambda_grid)) {
  fit_k <- fit_graph_fused_multinomial(Y, D_G_rook, lambda = lambda_grid[k],
                                       T_max = 4000, tol = 1e-6,
                                       beta_init = beta_warm, verbose = FALSE)
  beta_warm <- fit_k$beta_hat
  mse_P_grid[k] <- mean((fit_k$P_hat - P_true)^2)
  cat(sprintf("lambda = %-5.2f | iters = %5d | MSE(P) = %.5f\n",
              lambda_grid[k], fit_k$n_iter, mse_P_grid[k]))
}
best_k <- which.min(mse_P_grid)
cat(sprintf("\nOracle-best lambda (rook graph): %.2f\n", lambda_grid[best_k]))

## --- select lambda via cross-validation (rook graph) ------------------
cv_out <- cv_graph_fused_multinomial(Y, D_G_rook, lambda_grid, K = 5, seed = 42,
                                     T_max = 4000, tol = 1e-6, verbose = FALSE)
cat("\n=== Cross-validation results (rook graph) ===\n")
print(data.frame(lambda = cv_out$lambda_grid,
                 cv_mean = round(cv_out$cv_mean, 5),
                 cv_se   = round(cv_out$cv_se, 5)))
cat(sprintf("lambda.min = %.3f | lambda.1se = %.3f\n",
            cv_out$lambda_min, cv_out$lambda_1se))

## --- final fits at the CV-selected lambda, on BOTH graphs -------------
fit_rook <- fit_graph_fused_multinomial(Y, D_G_rook, lambda = cv_out$lambda_min,
                                        T_max = 5000, tol = 1e-6, verbose = FALSE)
fit_king <- fit_graph_fused_multinomial(Y, D_G_king, lambda = cv_out$lambda_min,
                                        T_max = 5000, tol = 1e-6, verbose = FALSE)

mse_P_rook <- mean((fit_rook$P_hat - P_true)^2)
mse_P_king <- mean((fit_king$P_hat - P_true)^2)
cat(sprintf("\nMSE(P) at lambda.min, rook graph:               %.5f\n", mse_P_rook))
cat(sprintf("MSE(P) at lambda.min, king graph (+ diagonals):  %.5f\n", mse_P_king))

## ============================================================
## 7. Diagnostics and plots
## ============================================================

## (a) Convergence of the residual (rook graph fit)
df_conv <- data.frame(iter = seq_along(fit_rook$resid_hist),
                      residual = fit_rook$resid_hist)
p_conv <- ggplot(df_conv, aes(iter, residual)) +
  geom_line(color = "steelblue") +
  scale_y_log10() +
  labs(title = "Primal-Dual convergence (rook graph)",
       x = "Iteration", y = "Relative residual (log scale)") +
  theme_minimal()
print(p_conv)

## (b) MSE(P) vs lambda (oracle sweep, rook graph)
df_lambda <- data.frame(lambda = lambda_grid, mse_P = mse_P_grid)
p_lambda <- ggplot(df_lambda, aes(lambda, mse_P)) +
  geom_line(color = "darkred") + geom_point() +
  labs(title = "Prediction error vs. regularization strength",
       x = expression(lambda), y = "MSE(P_hat, P_true)") +
  theme_minimal()
print(p_lambda)

## (c) Spatial recovery: True X vs Rook-graph estimate vs King-graph estimate
row_v <- V(g_rook)$row
col_v <- V(g_rook)$col
df_map <- data.frame(
  row = row_v, col = col_v,
  true = P_true[, 1],
  rook = fit_rook$P_hat[, 1],
  king = fit_king$P_hat[, 1]
)
df_map_long <- reshape2::melt(df_map, id.vars = c("row", "col"))
df_map_long$variable <- factor(df_map_long$variable, levels = c("true", "rook", "king"))

p_map <- ggplot(df_map_long, aes(col, row, fill = value)) +
  geom_tile() +
  facet_wrap(~variable, labeller = as_labeller(c(
    true = "True p_1 (metabolite 1)",
    rook = "Estimated (rook graph)",
    king = "Estimated (king graph)"
  ))) +
  scale_fill_viridis_c(name = "p_1") +
  scale_y_reverse() +                      # row 1 at the top, matching the uploaded figure's axes
  labs(title = "Recovery of an X-shaped spatial pattern", x = NULL, y = NULL) +
  theme_minimal() +
  coord_fixed()
print(p_map)

## (d) CV error curve with +/-1 SE bands and the two selection rules marked
df_cv <- data.frame(lambda = cv_out$lambda_grid,
                    mean = cv_out$cv_mean,
                    lo = cv_out$cv_mean - cv_out$cv_se,
                    hi = cv_out$cv_mean + cv_out$cv_se)
p_cv <- ggplot(df_cv, aes(lambda, mean)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") +
  geom_line(color = "black") + geom_point() +
  geom_vline(xintercept = cv_out$lambda_min, color = "darkgreen", linetype = "dashed") +
  geom_vline(xintercept = cv_out$lambda_1se, color = "blue", linetype = "dotted") +
  labs(title = "Cross-validated held-out deviance vs. lambda",
       subtitle = "dashed = lambda.min, dotted = lambda.1se",
       x = expression(lambda), y = "Mean held-out deviance per trial") +
  theme_minimal()
print(p_cv)

## ============================================================
## 8. Summary
## ============================================================

cat("\n=== Simulation summary ===\n")
cat(sprintf("Grid: %d x %d, n = %d nodes\n", nrow_grid, ncol_grid, n))
cat(sprintf("Rook graph: |E| = %d edges (4-neighbor, no diagonals)\n", nrow(D_G_rook)))
cat(sprintf("King graph: |E| = %d edges (rook + diagonals)\n", nrow(D_G_king)))
cat(sprintf("Categories D = %d, average n_i = %.1f\n", D, mean(n_i)))
cat(sprintf("Oracle-best lambda (rook) = %.2f | CV lambda.min = %.3f | CV lambda.1se = %.3f\n",
            lambda_grid[best_k], cv_out$lambda_min, cv_out$lambda_1se))
cat(sprintf("MSE(P), rook graph: %.5f | MSE(P), king graph: %.5f\n",
            mse_P_rook, mse_P_king))

## ============================================================
## 9. Variant: build the graph from a k-nearest-neighbor graph
## ============================================================
##
## Real spatial/network data are often not observed on a regular grid --
## sampling locations may be irregularly scattered (e.g., unevenly spaced
## field sites, sensor placements, or subjects' GPS coordinates). This
## section repeats the X-pattern simulation with such a graph: n locations
## are scattered uniformly at random in the unit square, and G is built by
## connecting each node to its k nearest neighbors (Euclidean distance).
## Everything downstream of graph construction -- the incidence matrix, the
## fitting algorithm, and the cross-validation routine -- is unchanged, since
## none of it assumes a grid; only graph construction, truth generation, and
## plotting need to be adapted to irregular coordinates.

## --- k-NN graph builder ------------------------------------------------
## Connects each point to its k nearest neighbors and symmetrizes by taking
## the UNION of the two directed neighbor relations (i.e., an edge (i, i')
## is included if i' is among i's k nearest neighbors OR i is among i''s).
## This tends to give a better-connected graph than requiring MUTUAL nearest
## neighbors, at the cost of some nodes having degree slightly above k.
build_knn_graph <- function(coords, k) {
  n <- nrow(coords)
  nn <- FNN::get.knn(coords, k = k)
  from <- rep(1:n, each = k)
  to   <- as.vector(t(nn$nn.index))
  edges <- unique(t(apply(cbind(from, to), 1, sort)))   # undirected, deduplicated
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  if (igraph::vcount(g) < n) g <- igraph::add_vertices(g, n - igraph::vcount(g))
  g <- igraph::simplify(g)                              # drop any self-loops / repeats
  V(g)$x <- coords[, 1]
  V(g)$y <- coords[, 2]
  g
}

## --- true X-shaped pattern from arbitrary (x, y) coordinates in [0,1]^2 -
## Same construction as generate_true_beta_X(), just parameterized directly
## by continuous coordinates instead of grid row/col indices.
generate_true_beta_X_coords <- function(coords, D, p_high = 0.70, p_low = 0.08,
                                        half_width = 0.06) {
  n <- nrow(coords)
  x <- coords[, 1]; y <- coords[, 2]
  dist_main <- abs(x - y) / sqrt(2)            # perpendicular distance to the "\" diagonal
  dist_anti <- abs(x + y - 1) / sqrt(2)        # perpendicular distance to the "/" diagonal
  on_X <- pmin(dist_main, dist_anti) < half_width
  
  p1 <- ifelse(on_X, p_low, p_high)
  p_rest <- (1 - p1) / (D - 1)
  P_true <- cbind(p1, matrix(p_rest, n, D - 1))
  
  Dm1 <- D - 1
  beta_true <- log(P_true[, 1:Dm1, drop = FALSE] / P_true[, D])
  list(beta_true = beta_true, P_true = P_true, on_X = on_X)
}

## --- random locations and two candidate k-NN graphs (small vs large k) --
set.seed(7)
n_pts <- 2000
coords <- cbind(x = runif(n_pts), y = runif(n_pts))

k_small <- 5
k_large <- 15
g_knn_small <- build_knn_graph(coords, k = k_small)
g_knn_large <- build_knn_graph(coords, k = k_large)

cat(sprintf("kNN graph (k=%d): |E| = %d, connected = %s\n",
            k_small, igraph::ecount(g_knn_small), igraph::is_connected(g_knn_small)))
cat(sprintf("kNN graph (k=%d): |E| = %d, connected = %s\n",
            k_large, igraph::ecount(g_knn_large), igraph::is_connected(g_knn_large)))

D_G_knn_small <- build_incidence_matrix(g_knn_small)
D_G_knn_large <- build_incidence_matrix(g_knn_large)

## --- truth and data (shared across both graphs; only the fitting graph differs)
truth_knn <- generate_true_beta_X_coords(coords, D, p_high = 0.70, p_low = 0.08, half_width = 0.06)
P_true_knn <- truth_knn$P_true

n_i_knn <- sample(150:300, n_pts, replace = TRUE)
Y_knn <- simulate_data(P_true_knn, n_i_knn)

## --- select lambda via CV (small-k graph), then fit both graphs at it ---
lambda_grid_knn <- c(0, 0.05, 0.1, 0.25, 0.5, 1, 2, 4)
cv_knn <- cv_graph_fused_multinomial(Y_knn, D_G_knn_small, lambda_grid_knn, K = 5, seed = 11,
                                     T_max = 4000, tol = 1e-6, verbose = FALSE)
cat(sprintf("\nkNN CV (k=%d): lambda.min = %.3f | lambda.1se = %.3f\n",
            k_small, cv_knn$lambda_min, cv_knn$lambda_1se))

fit_knn_small <- fit_graph_fused_multinomial(Y_knn, D_G_knn_small, lambda = cv_knn$lambda_min,
                                             T_max = 5000, tol = 1e-6, verbose = FALSE)
fit_knn_large <- fit_graph_fused_multinomial(Y_knn, D_G_knn_large, lambda = cv_knn$lambda_min,
                                             T_max = 5000, tol = 1e-6, verbose = FALSE)

mse_knn_small <- mean((fit_knn_small$P_hat - P_true_knn)^2)
mse_knn_large <- mean((fit_knn_large$P_hat - P_true_knn)^2)
cat(sprintf("MSE(P), kNN graph k=%d: %.5f\n", k_small, mse_knn_small))
cat(sprintf("MSE(P), kNN graph k=%d: %.5f\n", k_large, mse_knn_large))

## ============================================================
## 10. kNN-graph diagnostics and plots
## ============================================================

## (a) The graph itself, to see the irregular connectivity structure
el_small <- igraph::as_edgelist(g_knn_small, names = FALSE)
df_edges <- data.frame(x = coords[el_small[, 1], 1], y = coords[el_small[, 1], 2],
                       xend = coords[el_small[, 2], 1], yend = coords[el_small[, 2], 2])
df_nodes <- data.frame(x = coords[, 1], y = coords[, 2])

p_graph <- ggplot() +
  geom_segment(data = df_edges, aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey70", linewidth = 0.15) +
  geom_point(data = df_nodes, aes(x, y), size = 0.4, color = "steelblue") +
  coord_fixed() + theme_minimal() +
  labs(title = sprintf("k-NN graph (k = %d) on %d random locations", k_small, n_pts),
       x = NULL, y = NULL)
print(p_graph)

## (b) True vs. recovered X pattern (small-k vs large-k graph)
df_map_knn <- rbind(
  data.frame(x = coords[, 1], y = coords[, 2], value = P_true_knn[, 1], panel = "true"),
  data.frame(x = coords[, 1], y = coords[, 2], value = fit_knn_small$P_hat[, 1],
             panel = sprintf("kNN estimate (k=%d)", k_small)),
  data.frame(x = coords[, 1], y = coords[, 2], value = fit_knn_large$P_hat[, 1],
             panel = sprintf("kNN estimate (k=%d)", k_large))
)
df_map_knn$panel <- factor(df_map_knn$panel,
                           levels = c("true", sprintf("kNN estimate (k=%d)", k_small),
                                      sprintf("kNN estimate (k=%d)", k_large)))

p_map_knn <- ggplot(df_map_knn, aes(x, y, color = value)) +
  geom_point(size = 0.6) +
  facet_wrap(~panel) +
  scale_color_viridis_c(name = "p_1") +
  coord_fixed() + theme_minimal() +
  labs(title = "Recovery of an X-shaped pattern on a k-NN graph", x = NULL, y = NULL)
print(p_map_knn)

## ============================================================
## 11. Summary (kNN variant)
## ============================================================

cat("\n=== kNN graph simulation summary ===\n")
cat(sprintf("n = %d randomly scattered points in [0,1]^2\n", n_pts))
cat(sprintf("Small-k graph: k=%d, |E|=%d\n", k_small, igraph::ecount(g_knn_small)))
cat(sprintf("Large-k graph: k=%d, |E|=%d\n", k_large, igraph::ecount(g_knn_large)))
cat(sprintf("CV-selected lambda.min = %.3f\n", cv_knn$lambda_min))
cat(sprintf("MSE(P): k=%d -> %.5f | k=%d -> %.5f\n",
            k_small, mse_knn_small, k_large, mse_knn_large))

## ============================================================
## 12. Effect of sample size
## ============================================================

mse_knn <- c()
n_can <- seq(50, 1000, 100)
D <- 4
for(n in n_can)
{
  coords <- cbind(x = runif(n), y = runif(n))
  k <- 5 
  
  g_knn <- build_knn_graph(coords, k)
  D_G_knn <- build_incidence_matrix(g_knn)
  
  truth_knn <- generate_true_beta_X_coords(coords, D, p_high = 0.70, p_low = 0.08, half_width = 0.06)
  P_true_knn <- truth_knn$P_true
  
  n_i_knn <- sample(150:300, n, replace = TRUE)
  Y_knn <- simulate_data(P_true_knn, n_i_knn)
  
  fit_knn <- fit_graph_fused_multinomial(Y_knn, D_G_knn, lambda = 4,
                                               T_max = 5000, tol = 1e-6, verbose = FALSE)
  
  mse <- mean((fit_knn$P_hat - P_true_knn)^2)/mean(P_true_knn^2)
  mse_knn <- c(mse_knn, mse)
  cat(sprintf("MSE(P), n = %d, kNN graph k=%d: %.5f\n", n, k, mse))
}

mse_n_df <- data.frame(n_can, mse_knn)
ggplot(data = mse_n_df, aes(n_can, mse_knn)) + geom_point() + geom_line() +
  labs(title = "MSE for varying values of n, k = 5", x = "Sample size", y = "MSE") 

mse_lat <- c()
n_lat <- seq(5, 75, 10)

D <- 4

for(n in n_lat)
{
  
  g_lat <- build_grid_graph(n, n)
  D_G_lat <- build_incidence_matrix(g_lat)
  
  truth_lat <- generate_true_beta_X(g_lat, D, n, n,
                                p_high = 0.70, p_low = 0.08, half_width = 0.06)
  beta_true <- truth_lat$beta_true
  P_true <- truth_lat$P_true
  
  n_i <- sample(150:300, n^2, replace = TRUE)
  Y <- simulate_data(P_true, n_i)
  
  fit_lat <- fit_graph_fused_multinomial(Y, D_G_lat, lambda = 4,
                                         T_max = 5000, tol = 1e-6, verbose = FALSE)
  
  mse <- mean((fit_lat$P_hat - P_true)^2)/mean(P_true^2)
  mse_lat <- c(mse_lat, mse)
  cat(sprintf("MSE(P), n = %d: %.5f\n", n, mse))
}

mse_n_df <- data.frame(n_lat^2, mse_lat)
ggplot(data = mse_n_df, aes(n_lat^2, mse_lat)) + geom_point() + geom_line() +
  labs(title = "MSE for varying values of n, lattice graph", x = "Sample size", y = "MSE") 


# validity of cv
## ============================================================
## 2. Oracle-vs-CV comparison across sample sizes
## ============================================================

D <- 4                                            # number of categories
lambda_grid <- c(0.5, 1, 2.5, 4, 6) # candidate lambda values
n_targets <- c(100, 500, 1000, 2000)              # target number of grid nodes

results <- data.frame(
  n_target    = n_targets,
  grid_dim    = NA_character_,
  n_actual    = NA_integer_,
  n_edges     = NA_integer_,
  lambda_oracle = NA_real_,
  mse_oracle    = NA_real_,
  lambda_cv     = NA_real_,
  mse_cv        = NA_real_,
  mse_diff      = NA_real_,     # mse_cv - mse_oracle
  pct_increase  = NA_real_      # 100 * mse_diff / mse_oracle
)

for (row in seq_along(n_targets)) {
  
  ## --- build a (roughly) square grid with ~ n_targets[row] nodes ---
  side <- round(sqrt(n_targets[row]))
  g <- build_grid_graph(side, side)
  n_actual <- igraph::vcount(g)
  D_G <- build_incidence_matrix(g)
  
  cat(sprintf("\n=== n_target = %d  (grid %dx%d, n = %d, |E| = %d) ===\n",
              n_targets[row], side, side, n_actual, nrow(D_G)))
  
  ## --- true X-shaped pattern and simulated data ---
  truth <- generate_true_beta_X(g, D, side, side, p_high = 0.70, p_low = 0.08, half_width = 0.06)
  P_true <- truth$P_true
  n_i <- sample(150:300, n_actual, replace = TRUE)
  Y <- simulate_data(P_true, n_i)
  
  ## --- oracle: sweep lambda_grid, pick the one minimizing MSE(P) vs truth ---
  mse_grid <- numeric(length(lambda_grid))
  beta_warm <- NULL
  for (k in seq_along(lambda_grid)) {
    fit_k <- fit_graph_fused_multinomial(Y, D_G, lambda = lambda_grid[k],
                                         T_max = 4000, tol = 1e-6,
                                         beta_init = beta_warm, verbose = FALSE)
    beta_warm <- fit_k$beta_hat
    mse_grid[k] <- mean((fit_k$P_hat - P_true)^2)
  }
  best_idx <- which.min(mse_grid)
  lambda_oracle <- lambda_grid[best_idx]
  mse_oracle <- mse_grid[best_idx]
  
  ## --- CV: select lambda.min via 5-fold count-splitting CV, then fit at it ---
  cv_out <- cv_graph_fused_multinomial(Y, D_G, lambda_grid, K = 5, seed = 42,
                                       T_max = 4000, tol = 1e-6, verbose = FALSE)
  fit_cv <- fit_graph_fused_multinomial(Y, D_G, lambda = cv_out$lambda_min,
                                        T_max = 4000, tol = 1e-6, verbose = FALSE)
  lambda_cv <- cv_out$lambda_min
  mse_cv <- mean((fit_cv$P_hat - P_true)^2)
  
  cat(sprintf("  oracle: lambda = %.2f | MSE = %.5f\n", lambda_oracle, mse_oracle))
  cat(sprintf("  CV:     lambda = %.2f | MSE = %.5f\n", lambda_cv, mse_cv))
  
  results$grid_dim[row]      <- sprintf("%dx%d", side, side)
  results$n_actual[row]      <- n_actual
  results$n_edges[row]       <- nrow(D_G)
  results$lambda_oracle[row] <- lambda_oracle
  results$mse_oracle[row]    <- mse_oracle
  results$lambda_cv[row]     <- lambda_cv
  results$mse_cv[row]        <- mse_cv
  results$mse_diff[row]      <- mse_cv - mse_oracle
  results$pct_increase[row]  <- 100 * (mse_cv - mse_oracle) / mse_oracle
}

## ============================================================
## 3. Report the table
## ============================================================

results_print <- results
results_print$mse_oracle   <- round(results_print$mse_oracle, 5)
results_print$mse_cv       <- round(results_print$mse_cv, 5)
results_print$mse_diff     <- round(results_print$mse_diff, 5)
results_print$pct_increase <- round(results_print$pct_increase, 2)

cat("\n\n=== Oracle vs. CV lambda: MSE comparison (grid graph, X-shaped truth) ===\n\n")
print(results_print, row.names = FALSE)

## Save the table for use in the paper (CSV + a LaTeX table if xtable is available)
## Files are written to the current working directory -- run getwd() to check
## where that is, or set one explicitly with setwd("your/path/here") first.
write.csv(results_print, "oracle_vs_cv_table.csv", row.names = FALSE)

if (requireNamespace("xtable", quietly = TRUE)) {
  tex_tab <- xtable::xtable(
    results_print[, c("n_actual", "n_edges", "lambda_oracle", "mse_oracle",
                      "lambda_cv", "mse_cv", "mse_diff", "pct_increase")],
    caption = "Oracle vs. cross-validated lambda: MSE(P) comparison o
    
    
    n the grid graph across sample sizes.",
    label = "tab:oracle_vs_cv",
    digits = c(0, 0, 0, 2, 5, 2, 5, 5, 2)
  )
  sink("oracle_vs_cv_table.tex")
  print(tex_tab, include.rownames = FALSE)
  sink()
  cat(sprintf("\nLaTeX table written to %s/oracle_vs_cv_table.tex\n", getwd()))
} else {
  cat("\n(Install the 'xtable' package to also generate a ready-to-use LaTeX table.)\n")
}

# ============================================================
## 12. Variant: a general (non-spatial) graph -- stochastic block model
## ============================================================
##
## The grid and k-NN graphs above are both spatially embedded: two nodes
## are connected because they are geometrically close, and the true
## pattern was defined directly as a function of (x, y) coordinates. Many
## graphs relevant to compositional data have no such embedding at all --
## e.g. co-occurrence networks, correlation networks, contact networks,
## or social networks, where connectivity reflects shared community or
## function rather than physical proximity, and there is no natural 2-D
## layout to define a "spatial pattern" on. This section repeats the
## exercise on such a graph, generated from a stochastic block model
## (SBM): nodes are assigned to one of several communities, and edges
## occur far more often within a community than between communities. The
## true metabolite-1 rate is then constant within each community and
## differs across communities -- the natural notion of "smoothness" for
## this kind of graph is shared community membership, not geometric
## adjacency, and nothing in build_incidence_matrix(),fit_graph_fused_
## multinomial(), or cv_graph_fused_multinomial() needs to change to
## accommodate that: they operate on the edge list alone.

## --- stochastic block model graph ---------------------------------------
build_sbm_graph <- function(block_sizes, p_in = 0.03, p_out = 0.001) {
  n_blocks <- length(block_sizes)
  pref_matrix <- matrix(p_out, n_blocks, n_blocks)
  diag(pref_matrix) <- p_in
  g <- igraph::sample_sbm(sum(block_sizes), pref.matrix = pref_matrix, block.sizes = block_sizes)
  V(g)$block <- rep(seq_along(block_sizes), times = block_sizes)
  g
}

## --- true pattern: constant within each community, varies across them --
# --- true pattern: SIMILAR (not identical) within each community,
##     centered on a distinct level per community -------------------------
## Noise is added on the logit scale, centered at each node's community
## mean, so that p1 stays automatically bounded in (0, 1) with no clipping,
## and larger within_block_sd gives more within-community heterogeneity
## relative to the between-community differences in p1_by_block.
generate_true_beta_blocks <- function(g, D, p1_by_block, within_block_sd = 0.35) {
  n <- igraph::vcount(g)
  block <- V(g)$block
  p1_center <- p1_by_block[block]
  logit_center <- log(p1_center / (1 - p1_center))
  logit_node <- logit_center + rnorm(n, mean = 0, sd = within_block_sd)
  p1 <- plogis(logit_node)                      # back to (0, 1), no clipping needed
  
  p_rest <- (1 - p1) / (D - 1)                  # now varies by node, not just by block
  P_true <- cbind(p1, matrix(rep(p_rest, D - 1), n, D - 1))
  
  Dm1 <- D - 1
  beta_true <- log(P_true[, 1:Dm1, drop = FALSE] / P_true[, D])
  list(beta_true = beta_true, P_true = P_true, block = block)
}

## --- build a 5-community graph with a distinct rate per community -------
set.seed(21)
block_sizes <- c(300, 250, 400, 350, 300)          # n = 1600, 5 communities
g_gen <- build_sbm_graph(block_sizes, p_in = 0.03, p_out = 0.001)
n_gen <- igraph::vcount(g_gen)
D_G_gen <- build_incidence_matrix(g_gen)

cat(sprintf("General (SBM) graph: n = %d, |E| = %d, connected = %s\n",
            n_gen, igraph::ecount(g_gen), igraph::is_connected(g_gen)))
## If FALSE, raise p_out slightly (e.g. to 0.002) and rerun -- a
## disconnected graph just means the fusion penalty cannot pool
## information across components, which is still valid but changes
## what the method can achieve.

p1_by_block <- c(0.75, 0.10, 0.55, 0.20, 0.40)     # distinct metabolite-1 rate per community
truth_gen <- generate_true_beta_blocks(g_gen, D, p1_by_block)
P_true_gen <- truth_gen$P_true

n_i_gen <- sample(150:300, n_gen, replace = TRUE)
Y_gen <- simulate_data(P_true_gen, n_i_gen)

## --- select lambda via CV, then fit --------------------------------------
lambda_grid_gen <- c(0, 0.05, 0.1, 0.25, 0.5, 1, 2, 4)
cv_gen <- cv_graph_fused_multinomial(Y_gen, D_G_gen, lambda_grid_gen, K = 5, seed = 5,
                                     T_max = 4000, tol = 1e-6, verbose = FALSE)
cat(sprintf("General graph CV: lambda.min = %.3f | lambda.1se = %.3f\n",
            cv_gen$lambda_min, cv_gen$lambda_1se))

fit_gen <- fit_graph_fused_multinomial(Y_gen, D_G_gen, lambda = cv_gen$lambda_min,
                                       T_max = 5000, tol = 1e-6, verbose = FALSE)
mse_gen <- mean((fit_gen$P_hat - P_true_gen)^2)

## also fit unregularized (lambda = 0), to quantify what the fusion penalty buys you
fit_gen_unreg <- fit_graph_fused_multinomial(Y_gen, D_G_gen, lambda = 0,
                                             T_max = 5000, tol = 1e-6, verbose = FALSE)
mse_gen_unreg <- mean((fit_gen_unreg$P_hat - P_true_gen)^2)

cat(sprintf("MSE(P) at lambda = 0 (unregularized):     %.5f\n", mse_gen_unreg))
cat(sprintf("MSE(P) at lambda.min (graph-fused):       %.5f\n", mse_gen))

## ============================================================
## 13. General-graph diagnostics and plots
## ============================================================
##
## There are no real coordinates for this graph, so we compute a
## force-directed layout purely for visualization -- it plays no role in
## fitting the model and is not used anywhere above this point.
set.seed(99)
layout_xy <- igraph::layout_with_fr(g_gen)
colnames(layout_xy) <- c("x", "y")

## (a) community structure
df_comm <- data.frame(x = layout_xy[, 1], y = layout_xy[, 2], block = factor(V(g_gen)$block))
p_comm <- ggplot(df_comm, aes(x, y, color = block)) +
  geom_point(size = 0.8) +
  coord_fixed() + theme_minimal() +
  labs(title = "Community structure of the general (SBM) graph",
       x = NULL, y = NULL, color = "community")
print(p_comm)

## (b) true vs. unregularized vs. graph-fused estimate of p_1
df_map_gen <- rbind(
  data.frame(x = layout_xy[, 1], y = layout_xy[, 2], value = P_true_gen[, 1], panel = "true"),
  data.frame(x = layout_xy[, 1], y = layout_xy[, 2], value = fit_gen_unreg$P_hat[, 1],
             panel = "unregularized (lambda=0)"),
  data.frame(x = layout_xy[, 1], y = layout_xy[, 2], value = fit_gen$P_hat[, 1],
             panel = "graph-fused (lambda.min)")
)
df_map_gen$panel <- factor(df_map_gen$panel,
                           levels = c("true", "unregularized (lambda=0)", "graph-fused (lambda.min)"))

p_map_gen <- ggplot(df_map_gen, aes(x, y, color = value)) +
  geom_point(size = 0.8) +
  facet_wrap(~panel) +
  scale_color_viridis_c(name = "p_1") +
  coord_fixed() + theme_minimal() +
  labs(title = "Recovery of a community-structured pattern on a general graph",
       x = NULL, y = NULL)
print(p_map_gen)

## ============================================================
## 14. Summary (general graph variant)
## ============================================================

cat("\n=== General (SBM) graph simulation summary ===\n")
cat(sprintf("n = %d, communities = %d, |E| = %d, connected = %s\n",
            n_gen, length(block_sizes), igraph::ecount(g_gen), igraph::is_connected(g_gen)))
cat(sprintf("CV-selected lambda.min = %.3f | lambda.1se = %.3f\n",
            cv_gen$lambda_min, cv_gen$lambda_1se))
cat(sprintf("MSE(P): unregularized = %.5f | graph-fused = %.5f (%.1f%% reduction)\n",
            mse_gen_unreg, mse_gen, 100 * (mse_gen_unreg - mse_gen) / mse_gen_unreg))


mse_sbm <- c()
n_can <- seq(50, 2000, 200)
D <- 4
n<-100
for(n in n_can)
{
  k <- 0
  block_sizes <- c(n/5, 2*n/5, n/5, n/5)         
  g_gen <- build_sbm_graph(block_sizes, p_in = 0.03, p_out = 0.001)
  n_gen <- igraph::vcount(g_gen)
  D_G_gen <- build_incidence_matrix(g_gen)
  
  cat(sprintf("General (SBM) graph: n = %d, |E| = %d, connected = %s\n",
              n_gen, igraph::ecount(g_gen), igraph::is_connected(g_gen)))
  ## If FALSE, raise p_out slightly (e.g. to 0.002) and rerun -- a
  ## disconnected graph just means the fusion penalty cannot pool
  ## information across components, which is still valid but changes
  ## what the method can achieve.
  k
  p1_by_block <- c(0.75, 0.30, 0.55, 0.20)  # distinct metabolite-1 rate per community
  truth_gen <- generate_true_beta_blocks(g_gen, D, p1_by_block)
  P_true_gen <- truth_gen$P_true
  
  n_i_gen <- sample(150:300, n, replace = TRUE)
  
  Y_gen <- simulate_data(P_true_gen, n_i_gen)
 
  
  fit_sbm <- fit_graph_fused_multinomial(Y_gen, D_G_gen, lambda = 0.5,
                                         T_max = 5000, tol = 1e-6, verbose = FALSE)
  
  mse <- mean((fit_sbm$P_hat - P_true_gen)^2)/mean(P_true_gen^2)
  mse_sbm <- c(mse_sbm, mse)
  cat(sprintf("MSE(P), n = %d: %.5f\n", n, mse))
}

mse_n_df <- data.frame(n_can, mse_knn)
ggplot(data = mse_n_df, aes(n_can, mse_knn)) + geom_point() + geom_line() +
  labs(title = "MSE for varying values of n, k = 5", x = "Sample size", y = "MSE") 
# Numerical and integration checks. Run from the project root:
#   Rscript tests/validate.R
# No testing packages are required.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
source(file.path(project_dir, "R", "load_all.R"))

expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}
check <- function(label, expr) {
  force(expr)
  cat("PASS:", label, "\n")
}
objective <- function(Y, D_G, beta, lambda) {
  -sum(Y * log_softmax_probs(beta)) + lambda * sum(abs(D_G %*% beta))
}

check("stable probabilities for extreme logits", {
  beta <- rbind(c(1000, -1000), c(-1000, -1000), c(0, 0))
  P <- softmax_probs(beta)
  log_P <- log_softmax_probs(beta)
  stopifnot(all(is.finite(P)), all(is.finite(log_P)),
            max(abs(rowSums(P) - 1)) < 1e-12,
            max(abs(P - exp(log_P))) < 1e-12,
            P[1, 1] == 1, P[2, 3] == 1)
})

check("rook and king graph topology, including one-dimensional grids", {
  for (nr in 1:4) for (nc in 1:4) {
    rook <- build_grid_graph(nr, nc)
    king <- build_king_grid_graph(nr, nc)
    rook_edges <- (nr - 1) * nc + (nc - 1) * nr
    stopifnot(igraph::ecount(rook) == rook_edges,
              igraph::ecount(king) == rook_edges + 2 * (nr - 1) * (nc - 1),
              igraph::is_simple(rook), igraph::is_simple(king))
  }
})

check("incidence operator and guaranteed norm bound", {
  g <- build_grid_graph(3, 4)
  D_G <- build_incidence_matrix(g)
  dense <- as.matrix(D_G)
  stopifnot(all(rowSums(dense) == 0), all(rowSums(abs(dense)) == 2))
  exact <- max(eigen(crossprod(dense), symmetric = TRUE, only.values = TRUE)$values)
  bound <- .gfm_operator_bound(D_G)
  set.seed(25)
  old_seed <- .Random.seed
  estimated <- op_norm_sq_DG(D_G)
  stopifnot(bound >= exact - 1e-12, abs(estimated - exact) < 1e-6,
            identical(.Random.seed, old_seed))
  empty <- build_incidence_matrix(igraph::make_empty_graph(3, directed = FALSE))
  stopifnot(identical(dim(empty), c(0L, 3L)), op_norm_sq_DG(empty) == 0)
})

check("kNN self exclusion, complete graphs, and component-size bound", {
  coordinates <- rbind(c(0, 0), c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(2, 1))
  for (k in 1:5) {
    g <- build_knn_graph(coordinates, k)
    stopifnot(igraph::vcount(g) == 6, igraph::is_simple(g),
              all(igraph::degree(g) >= k),
              all(igraph::components(g)$csize >= k + 1))
    if (k == 5) stopifnot(igraph::ecount(g) == 15)
  }
  coincident <- build_knn_graph(matrix(0, 8, 2), 2)
  stopifnot(igraph::is_simple(coincident), all(igraph::degree(coincident) >= 2))
  stopifnot(igraph::vcount(build_knn_graph(matrix(1:6, ncol = 1), 2)) == 6)
})

check("SBM isolation repair respects two-vertex blocks", {
  for (seed in 1:30) {
    set.seed(seed)
    g <- build_sbm_graph(c(2, 2), p_in = 0, p_out = 0)
    edges <- igraph::as_edgelist(g, names = FALSE)
    block <- igraph::V(g)$block
    stopifnot(all(igraph::degree(g) > 0), igraph::is_simple(g),
              all(block[edges[, 1]] == block[edges[, 2]]))
  }
  g <- build_sbm_graph(c(1, 1), p_in = 0, p_out = 0)
  stopifnot(igraph::ecount(g) == 1)
  expect_error(build_sbm_graph(1))
})

check("truth generators and simulated multinomial totals", {
  for (D in c(2, 4, 10)) for (nr in c(1, 3)) {
    g <- build_grid_graph(nr, 4)
    truth <- generate_true_beta_X(g, D, nr, 4)
    smooth <- generate_true_beta(g, D, nr, 4)
    stopifnot(all(is.finite(smooth$beta_true)),
              max(abs(softmax_probs(truth$beta_true) - truth$P_true)) < 1e-12)
    n_i <- seq_len(igraph::vcount(g)) * 10
    Y <- simulate_data(truth$P_true, n_i)
    stopifnot(all(rowSums(Y) == n_i), all(Y >= 0))
  }
  g <- build_sbm_graph(c(3, 3))
  block_truth <- generate_true_beta_blocks(g, 5, c(0.8, 0.2), within_block_sd = 0)
  stopifnot(max(abs(block_truth$P_true[, 1] -
                     c(rep(0.8, 3), rep(0.2, 3)))) < 1e-12)
  single <- generate_true_beta_X(build_grid_graph(1, 1), 2, 1, 1)
  stopifnot(all(is.finite(single$beta_true)))
})

check("two-node binary solver agrees with analytic penalized solutions", {
  Y <- rbind(c(18, 2), c(3, 17))
  D_G <- build_incidence_matrix(build_grid_graph(2, 1))
  for (lambda in c(0, 2.5, 10)) {
    fit <- fit_graph_fused_multinomial(Y, D_G, lambda, T_max = 20000,
                                       tol = 1e-10, verbose = FALSE)
    expected <- if (lambda < 7.5) c((18 - lambda) / 20, (3 + lambda) / 20) else
      rep(21 / 40, 2)
    stopifnot(fit$converged, max(abs(fit$P_hat[, 1] - expected)) < 1e-7,
              max(abs(rowSums(fit$P_hat) - 1)) < 1e-12,
              is.finite(objective(Y, D_G, fit$beta_hat, lambda)),
              1 / fit$tau - fit$sigma * fit$op2 > fit$Lf / 2)
  }
})

check("edgeless and zero-count fits", {
  D_G <- build_incidence_matrix(igraph::make_empty_graph(2, directed = FALSE))
  Y <- rbind(c(30, 20, 10), c(10, 30, 20))
  fit <- fit_graph_fused_multinomial(Y, D_G, 5, T_max = 10000,
                                     tol = 1e-10, verbose = FALSE)
  stopifnot(fit$converged, max(abs(fit$P_hat - Y / rowSums(Y))) < 1e-7)
  zero <- fit_graph_fused_multinomial(matrix(0, 2, 3), D_G, 1, verbose = FALSE)
  stopifnot(zero$converged, max(abs(zero$P_hat - 1 / 3)) < 1e-12,
            is.finite(zero$tau),
            is.finite(zero$sigma))
  mixed <- fit_graph_fused_multinomial(rbind(c(0, 0), c(10, 10)),
                                      build_incidence_matrix(build_grid_graph(2, 1)),
                                      1, verbose = FALSE)
  stopifnot(all(is.finite(mixed$P_hat)))
})

check("seeded count splitting preserves counts and the caller's RNG", {
  Y <- rbind(c(40, 30), c(20, 10))
  set.seed(7)
  old_seed <- .Random.seed
  folds <- split_counts_Kfold(Y, 5, seed = 20)
  stopifnot(identical(.Random.seed, old_seed),
            identical(folds, split_counts_Kfold(Y, 5, seed = 20)),
            all(apply(folds, c(1, 2), sum) == Y))
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(split_counts_Kfold(Y, 5, seed = 20))
  stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

check("CV scores match independent closed-form multinomial MLE scores", {
  Y <- matrix(c(60, 40), nrow = 1)
  D_G <- build_incidence_matrix(igraph::make_empty_graph(1, directed = FALSE))
  cv <- cv_graph_fused_multinomial(Y, D_G, c(2, 0, 2, 1), K = 5, seed = 6,
                                  T_max = 10000, tol = 1e-10,
                                  warm_start = FALSE, verbose = FALSE)
  folds <- split_counts_Kfold(Y, 5, seed = 6)
  exact_scores <- vapply(1:5, function(k) {
    test <- matrix(folds[, , k], 1, 2)
    train <- Y - test
    -sum(test * log(train / sum(train)))
  }, numeric(1))
  stopifnot(all(cv$converged), max(abs(cv$cv_mat - exact_scores)) < 1e-7,
            max(abs(cv$cv_mean - mean(exact_scores / cv$n_test))) < 1e-8,
            max(abs(cv$cv_se - sd(exact_scores / cv$n_test) / sqrt(5))) < 1e-8,
            max(abs(cv$cv_mat[, 1] - cv$cv_mat[, 3])) < 1e-12,
            cv$lambda_1se == 2)
})

check("CV ignores empty validation folds and rejects unusable splits", {
  Y <- matrix(c(10, 10), nrow = 1)
  D_G <- build_incidence_matrix(igraph::make_empty_graph(1, directed = FALSE))
  cv <- cv_graph_fused_multinomial(Y, D_G, c(0, 1), K = 50, seed = 4,
                                  T_max = 10000, tol = 1e-9, verbose = FALSE)
  stopifnot(any(cv$n_test == 0), all(is.finite(cv$cv_mean)),
            all(is.finite(cv$cv_se)), all(is.na(cv$cv_fold_rate[cv$n_test == 0, ])))
  expect_error(cv_graph_fused_multinomial(matrix(c(1, 0), 1), D_G,
                                         c(0, 1), verbose = FALSE))
})

check("clear invalid-input errors and convergence diagnostics", {
  Y <- matrix(c(10, 10), nrow = 1)
  D_G <- build_incidence_matrix(igraph::make_empty_graph(1, directed = FALSE))
  expect_error(softmax_probs(matrix(Inf, 1, 1)))
  expect_error(build_knn_graph(matrix(0, 3, 2), 3))
  expect_error(simulate_data(matrix(c(0.2, 0.2), 1), 10))
  expect_error(fit_graph_fused_multinomial(Y, D_G, -1, verbose = FALSE))
  expect_error(fit_graph_fused_multinomial(Y, D_G, 1, tau = 100, verbose = FALSE))
  expect_error(split_counts_Kfold(Y, 1))
  fit <- fit_graph_fused_multinomial(matrix(c(9, 1), 1), D_G, 0,
                                     T_max = 1, tol = 1e-12, verbose = FALSE)
  stopifnot(!fit$converged, fit$n_iter == 1, fit$residual > 1e-12)
})

check("portable loading into a specified environment without attaching packages", {
  previous <- getwd()
  search_before <- search()
  setwd(tempdir())
  isolated_env <- new.env(parent = globalenv())
  source(file.path(project_dir, "R", "load_all.R"), local = isolated_env)
  setwd(previous)
  stopifnot(exists("fit_graph_fused_multinomial", envir = isolated_env,
                   inherits = FALSE), identical(search(), search_before),
            !exists("required_packages", envir = isolated_env, inherits = FALSE),
            !exists("script_dir", envir = isolated_env, inherits = FALSE))
})

cat("\nAll 13 validation groups passed.\n")
cat(R.version.string, "\n")
for (package in c("Matrix", "igraph", "FNN")) {
  cat(package, as.character(utils::packageVersion(package)), "\n")
}

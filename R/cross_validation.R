# K-fold cross-validation by splitting trials within each graph node.
# Every fold retains the full graph. Only the count matrices change.

#' Randomly partition each node's counts into K validation folds.
#' Conditional on Y[i, j], counts are assigned with equal fold probability.
#' This is equivalent to independently assigning the observed trials.
#' @param Y Nonnegative integer n by D count matrix.
#' @param K Integer number of folds, at least two.
#' @param seed Optional local random seed; NULL uses the caller's RNG.
#' @return An n by D by K array whose sum over folds is exactly Y.
split_counts_Kfold <- function(Y, K, seed = NULL) {
  .gfm_check_counts(Y)
  .gfm_check_scalar(K, "K", lower = 2, upper = .Machine$integer.max,
                    integer = TRUE)
  if (any(Y > .Machine$integer.max)) {
    stop("Individual counts exceed the range supported by rmultinom().",
         call. = FALSE)
  }
  .gfm_with_seed(seed, {
    folds <- array(0, dim = c(nrow(Y), ncol(Y), K),
                   dimnames = list(rownames(Y), colnames(Y), NULL))
    fold_probabilities <- rep(1 / K, K)
    for (i in seq_len(nrow(Y))) {
      for (j in seq_len(ncol(Y))) {
        if (Y[i, j] > 0) {
          folds[i, j, ] <- as.numeric(stats::rmultinom(
            1L, size = Y[i, j], prob = fold_probabilities
          ))
        }
      }
    }
    folds
  })
}

#' Select lambda using held-out negative log-likelihood per trial.
#' @param Y Nonnegative integer n by D count matrix.
#' @param D_G Incidence matrix with one column per node.
#' @param lambda_grid Nonempty vector of finite nonnegative candidates.
#' @param K Integer number of folds, at least two.
#' @param seed Local seed for count splitting.
#' @param T_max,tol Iteration limit and tolerance passed to the estimator.
#' @param warm_start Whether to reuse logits along each increasing lambda path.
#' @param verbose Whether to print fold scores.
#' @return Original CV fields, plus fold rates, active-fold indices,
#'   iteration counts, and convergence flags. cv_mean is the unweighted
#'   mean over nonempty folds of NLL divided by held-out trial count.
#'   cv_se is the fold-rate SD divided by sqrt(number of active folds).
#'   The one-SE choice is a heuristic because training folds overlap.
cv_graph_fused_multinomial <- function(Y, D_G, lambda_grid, K = 5, seed = 1,
                                       T_max = 3000, tol = 1e-6,
                                       warm_start = TRUE, verbose = TRUE) {
  .gfm_check_counts(Y)
  .gfm_check_operator(D_G)
  if (ncol(D_G) != nrow(Y)) {
    stop("D_G must have one column for each row of Y.", call. = FALSE)
  }
  if (!is.numeric(lambda_grid) || length(lambda_grid) == 0L ||
      any(!is.finite(lambda_grid)) || any(lambda_grid < 0)) {
    stop("lambda_grid must contain finite nonnegative numbers.", call. = FALSE)
  }
  .gfm_check_scalar(T_max, "T_max", lower = 1, integer = TRUE)
  .gfm_check_scalar(tol, "tol", lower = .Machine$double.xmin)
  .gfm_check_flag(warm_start, "warm_start")
  .gfm_check_flag(verbose, "verbose")
  folds <- split_counts_Kfold(Y, K, seed = seed)
  n_test <- vapply(seq_len(K), function(k) sum(folds[, , k]), numeric(1L))
  active_folds <- which(n_test > 0)
  if (length(active_folds) < 2L) {
    stop("CV requires at least two nonempty validation folds; ",
         "use more trials, fewer folds, or a different split seed.", call. = FALSE)
  }
  n_lambda <- length(lambda_grid)
  cv_mat <- matrix(0, K, n_lambda)
  n_iter <- matrix(NA_integer_, K, n_lambda)
  converged <- matrix(NA, K, n_lambda)
  order_lambda <- order(lambda_grid)

  for (k in active_folds) {
    Y_test <- matrix(folds[, , k], nrow(Y), ncol(Y), dimnames = dimnames(Y))
    Y_train <- Y - Y_test
    positive_counts <- Y_test > 0
    beta_warm <- NULL
    for (index in order_lambda) {
      lambda <- lambda_grid[index]
      fit <- fit_graph_fused_multinomial(
        Y_train, D_G, lambda = lambda, T_max = T_max, tol = tol,
        beta_init = beta_warm, verbose = FALSE
      )
      if (warm_start) beta_warm <- fit$beta_hat
      log_P <- log_softmax_probs(fit$beta_hat)
      # Omit zero cells explicitly: 0 * (-Inf) is undefined in floating point.
      nll <- -sum(Y_test[positive_counts] * log_P[positive_counts])
      cv_mat[k, index] <- nll
      n_iter[k, index] <- fit$n_iter
      converged[k, index] <- fit$converged
      if (verbose) {
        cat(sprintf("fold %d/%d | lambda = %7.4f | held-out NLL = %.2f\n",
                    k, K, lambda, nll))
      }
    }
  }
  cv_fold_rate <- matrix(NA_real_, K, n_lambda)
  cv_fold_rate[active_folds, ] <- sweep(cv_mat[active_folds, , drop = FALSE],
                                       1L, n_test[active_folds], "/")
  active_rates <- cv_fold_rate[active_folds, , drop = FALSE]
  cv_mean <- colMeans(active_rates)
  cv_se <- apply(active_rates, 2L, stats::sd) / sqrt(length(active_folds))
  min_index <- which.min(cv_mean)
  threshold <- cv_mean[min_index] + cv_se[min_index]
  lambda_min <- lambda_grid[min_index]
  lambda_1se <- max(lambda_grid[cv_mean <= threshold])

  failed <- sum(!converged[active_folds, , drop = FALSE])
  if (failed > 0L) {
    warning(failed, " CV fit(s) reached T_max before meeting tol. ",
            "Inspect converged and n_iter before relying on the selected lambda.",
            call. = FALSE)
  }
  list(lambda_grid = lambda_grid, cv_mean = cv_mean, cv_se = cv_se,
       lambda_min = lambda_min, lambda_1se = lambda_1se,
       cv_mat = cv_mat, n_test = n_test, cv_fold_rate = cv_fold_rate,
       active_folds = active_folds, n_iter = n_iter, converged = converged)
}

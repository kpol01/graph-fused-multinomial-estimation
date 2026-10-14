# Probability alignment, multinomial scores, and an optional CARBayes wrapper.

#' Compute raw multinomial proportions, preserving zero-total rows as NA.
multinomial_mle <- function(Y) {
  .gfm_check_counts(Y)
  totals <- rowSums(Y)
  P <- matrix(NA_real_, nrow(Y), ncol(Y), dimnames = dimnames(Y))
  active <- totals > 0
  P[active, ] <- sweep(Y[active, , drop = FALSE], 1L, totals[active], "/")
  P
}

gfm_align_probabilities <- function(P, Y, category_order = NULL) {
  .gfm_check_counts(Y)
  P <- as.matrix(P)
  .gfm_check_matrix(P, "P")
  if (!identical(dim(P), dim(Y)) || any(P < 0 | P > 1) ||
      any(abs(rowSums(P) - 1) > 1e-6)) {
    stop("P must match the count matrix and have probability rows summing to one.",
         call. = FALSE)
  }
  if (!is.null(category_order)) {
    if (is.null(colnames(Y)) || length(category_order) != ncol(P) ||
        anyNA(category_order) || anyDuplicated(category_order) ||
        !setequal(category_order, colnames(Y))) {
      stop("category_order must name the actual fitted columns exactly once.",
           call. = FALSE)
    }
    colnames(P) <- category_order
  }
  for (axis in c("rows", "columns")) {
    fitted <- if (axis == "rows") rownames(P) else colnames(P)
    expected <- if (axis == "rows") rownames(Y) else colnames(Y)
    if (!is.null(fitted) && !is.null(expected)) {
      if (anyNA(fitted) || anyDuplicated(fitted) || !setequal(fitted, expected)) {
        stop("Probability ", axis, " do not match the count matrix names.", call. = FALSE)
      }
      P <- if (axis == "rows") P[expected, , drop = FALSE] else P[, expected, drop = FALSE]
    }
  }
  dimnames(P) <- dimnames(Y)
  P
}

#' Score probabilities on held-out counts.
#' @return Equal-node/category MSE against nonempty test proportions, and
#'   population-weighted individual log loss and multicategory Brier score.
score_multinomial <- function(Y_test, P, category_order = NULL) {
  P <- gfm_align_probabilities(P, Y_test, category_order)
  totals <- rowSums(Y_test)
  active <- totals > 0
  if (!any(active)) stop("The test set contains no observations.", call. = FALSE)
  observed <- multinomial_mle(Y_test)
  positive <- Y_test > 0
  data.frame(
    tested_nodes = sum(active), test_population = sum(totals),
    mse = mean((observed[active, , drop = FALSE] - P[active, , drop = FALSE])^2),
    log_loss = -sum(Y_test[positive] * log(P[positive])) / sum(totals),
    brier = sum(totals * rowSums(P^2) - 2 * rowSums(Y_test * P) + totals) / sum(totals)
  )
}

relative_probability_error <- function(P, P_true) {
  .gfm_check_matrix(P_true, "P_true")
  .gfm_check_matrix(P, "P")
  if (!identical(dim(P), dim(P_true)) || mean(P_true^2) <= 0) {
    stop("P and P_true must have matching dimensions and a positive truth norm.",
         call. = FALSE)
  }
  sqrt(mean((P - P_true)^2) / mean(P_true^2))
}

# CARBayes returns fitted response counts in response-column order.
# Normalize row totals and explicitly declare that order when generic labels
# are returned. Named district rows remain checked rather than overwritten.
gfm_mcar_probabilities <- function(fit, Y) {
  counts <- as.matrix(fit$fitted.values)
  .gfm_check_matrix(counts, "MCAR fitted values")
  if (!identical(dim(counts), dim(Y)) || any(counts < 0) || any(rowSums(counts) <= 0)) {
    stop("MCAR fitted values must be nonnegative counts matching Y, with positive totals.",
         call. = FALSE)
  }
  P <- sweep(counts, 1L, rowSums(counts), "/")
  if (!is.null(colnames(Y)) &&
      (is.null(colnames(P)) || !setequal(colnames(P), colnames(Y)))) {
    colnames(P) <- colnames(Y)
  }
  gfm_align_probabilities(P, Y)
}

#' Fit the original intercept-only multinomial MCAR comparison model.
#' @param Y Training count matrix, with positive totals at every node.
#' @param W Nonnegative symmetric adjacency matrix with no isolated node.
#' @param burnin,n_sample,thin MCMC settings; n_sample includes burn-in.
#' @param seed Local seed for this fit.
#' @return The CARBayes fit, aligned probabilities, and measured elapsed seconds.
fit_mcar_multinomial <- function(Y, W, burnin = 5000L, n_sample = 20000L,
                                 thin = 10L, seed = 123L) {
  .gfm_check_counts(Y)
  W <- as.matrix(W)
  .gfm_check_matrix(W, "W")
  if (!identical(dim(W), c(nrow(Y), nrow(Y))) || any(W < 0) ||
      !isTRUE(isSymmetric(unname(W))) || any(diag(W) != 0) ||
      any(rowSums(W) <= 0) || any(rowSums(Y) <= 0)) {
    stop("MCAR requires an aligned symmetric adjacency with positive neighbor ",
         "and training totals at every node. Choose a common analysis domain explicitly.",
         call. = FALSE)
  }
  if (!is.null(rownames(Y)) && !is.null(rownames(W))) {
    if (!identical(rownames(W), rownames(Y)) || !identical(colnames(W), rownames(Y))) {
      stop("W must follow the training count row order.", call. = FALSE)
    }
  }
  .gfm_check_scalar(burnin, "burnin", lower = 0, integer = TRUE)
  .gfm_check_scalar(n_sample, "n_sample", lower = burnin + 1, integer = TRUE)
  .gfm_check_scalar(thin, "thin", lower = 1, integer = TRUE)
  gfm_require_packages("CARBayes")
  timing <- system.time({
    fit <- .gfm_with_seed(seed, CARBayes::MVS.CARleroux(
      formula = Y ~ 1, family = "multinomial", trials = rowSums(Y), W = W,
      burnin = burnin, n.sample = n_sample, thin = thin, verbose = FALSE
    ))
  })
  list(fit = fit, P_hat = gfm_mcar_probabilities(fit, Y),
       elapsed_seconds = unname(timing["elapsed"]))
}

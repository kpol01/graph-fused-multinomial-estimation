# Primal-dual solver for graph-fused multinomial logistic regression.

#' Fit the graph-fused multinomial estimator.
#' The objective is the multinomial negative log-likelihood (omitting
#' constants) plus lambda * sum(abs(D_G %*% beta)). The last category is
#' the reference category. The likelihood is not divided by trial count.
#' @param Y Nonnegative integer n by D count matrix, D >= 2.
#'   Zero-total rows are allowed, including rows produced by count splitting.
#' @param D_G Numeric or Matrix incidence matrix with n columns.
#' @param lambda Nonnegative regularization parameter.
#' @param tau,sigma Positive primal and dual step sizes, or NULL for defaults.
#' @param theta Extrapolation parameter in [0, 1]; default is one.
#' @param a Positive pseudocount used only for initialization.
#' @param T_max Positive integer iteration limit.
#' @param tol Positive tolerance for relative changes in both iterates.
#' @param beta_init Optional finite n by (D - 1) warm-start matrix.
#' @param verbose Whether to print progress.
#' @return Existing fit fields, plus converged and residual. op2 is a
#'   guaranteed upper bound for the squared operator norm, and Lf is a
#'   conservative gradient Lipschitz bound. Small iterate changes are
#'   a stopping diagnostic, not a certificate of objective optimality.
fit_graph_fused_multinomial <- function(Y, D_G, lambda,
                                        tau = NULL, sigma = NULL,
                                        theta = 1, a = 0.5,
                                        T_max = 2000, tol = 1e-6,
                                        beta_init = NULL,
                                        verbose = TRUE) {
  .gfm_check_counts(Y)
  .gfm_check_operator(D_G)
  .gfm_check_scalar(lambda, "lambda", lower = 0)
  .gfm_check_scalar(theta, "theta", lower = 0, upper = 1)
  .gfm_check_scalar(a, "a", lower = .Machine$double.xmin)
  .gfm_check_scalar(T_max, "T_max", lower = 1, integer = TRUE)
  .gfm_check_scalar(tol, "tol", lower = .Machine$double.xmin)
  .gfm_check_flag(verbose, "verbose")
  n <- nrow(Y)
  D <- ncol(Y)
  q <- D - 1L
  if (ncol(D_G) != n) {
    stop("D_G must have one column for each row of Y.", call. = FALSE)
  }
  n_i <- rowSums(Y)
  Y_nonreference <- Y[, seq_len(q), drop = FALSE]

  # Use a guaranteed norm bound rather than a power-iteration estimate.
  # Preserve the original conservative likelihood bound and safety margin.
  op2 <- .gfm_operator_bound(D_G)
  Lf <- max(n_i)
  if (is.null(sigma)) sigma <- if (op2 > 0) 1 / (2 * op2) else 1
  .gfm_check_scalar(sigma, "sigma", lower = .Machine$double.xmin)
  boundary <- Lf / 2 + sigma * op2
  if (is.null(tau)) tau <- if (boundary > 0) 0.9 / boundary else 1
  .gfm_check_scalar(tau, "tau", lower = .Machine$double.xmin)
  if (!(1 / tau - sigma * op2 > Lf / 2)) {
    stop("Step sizes must satisfy 1/tau - sigma*op2 > Lf/2, ",
         "where op2 is the guaranteed norm bound.", call. = FALSE)
  }

  if (is.null(beta_init)) {
    beta <- sweep(log(Y_nonreference + a), 1L, log(Y[, D] + a), "-")
  } else {
    .gfm_check_matrix(beta_init, "beta_init")
    if (!identical(dim(beta_init), c(n, q))) {
      stop("beta_init must have dimensions n by (D - 1).", call. = FALSE)
    }
    beta <- beta_init
  }
  U <- matrix(0, nrow(D_G), q)
  beta_bar <- beta
  resid_hist <- numeric(T_max)
  converged <- FALSE

  for (iteration in seq_len(T_max)) {
    beta_old <- beta
    U_old <- U
    U <- pmin(pmax(as.matrix(U_old + sigma * (D_G %*% beta_bar)),
                   -lambda), lambda)
    P_nonreference <- exp(.gfm_log_softmax(beta_old)[, seq_len(q),
                                                  drop = FALSE])
    gradient <- sweep(P_nonreference, 1L, n_i, "*") - Y_nonreference
    beta <- beta_old - tau * (gradient + as.matrix(Matrix::crossprod(D_G, U)))
    if (any(!is.finite(beta))) {
      stop("Nonfinite primal iterate at iteration ", iteration,
           "; check the counts, initialization, and step sizes.", call. = FALSE)
    }
    beta_bar <- beta + theta * (beta - beta_old)

    primal_change <- norm(beta - beta_old, "F") / max(1, norm(beta_old, "F"))
    dual_change <- norm(U - U_old, "F") / max(1, norm(U_old, "F"))
    resid_hist[iteration] <- max(primal_change, dual_change)
    if (verbose && iteration %% 100L == 0L) {
      cat(sprintf("iter %5d   residual = %.3e\n", iteration,
                  resid_hist[iteration]))
    }
    if (resid_hist[iteration] < tol) {
      converged <- TRUE
      break
    }
  }
  residual <- resid_hist[iteration]
  if (verbose && !converged) {
    message("Iteration limit reached; final residual = ",
            format(residual, scientific = TRUE))
  }
  P_hat <- softmax_probs(beta)
  dimnames(P_hat) <- dimnames(Y)
  list(beta_hat = beta, P_hat = P_hat, U = U,
       n_iter = iteration, resid_hist = resid_hist[seq_len(iteration)],
       tau = tau, sigma = sigma, Lf = Lf, op2 = op2,
       converged = converged, residual = residual)
}

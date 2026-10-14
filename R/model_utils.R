# Numerical utilities for the graph-fused multinomial model.
# The final category is the reference category, with logit fixed at zero.

# Internal argument checks shared by the project scripts.
.gfm_check_scalar <- function(x, name, lower = -Inf, upper = Inf,
                              integer = FALSE) {
  valid <- is.numeric(x) && !is.complex(x) && length(x) == 1L && is.finite(x) &&
    x >= lower && x <= upper && (!integer || x == floor(x))
  if (!valid) {
    stop(name, " must be a finite ", if (integer) "integer" else "number",
         " in [", lower, ", ", upper, "].", call. = FALSE)
  }
  invisible(x)
}

.gfm_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(x)
}

.gfm_check_matrix <- function(x, name, min_rows = 1L, min_cols = 1L) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) || nrow(x) < min_rows ||
      ncol(x) < min_cols || any(!is.finite(x))) {
    stop(name, " must be a finite numeric matrix with at least ", min_rows,
         " row(s) and ", min_cols, " column(s).", call. = FALSE)
  }
  invisible(x)
}

.gfm_check_counts <- function(Y) {
  .gfm_check_matrix(Y, "Y", min_cols = 2L)
  if (any(Y < 0) || any(Y != floor(Y))) {
    stop("Y must contain nonnegative integer counts.", call. = FALSE)
  }
  invisible(Y)
}

.gfm_check_operator <- function(D_G) {
  if (!(is.matrix(D_G) && is.numeric(D_G)) &&
      !inherits(D_G, "Matrix")) {
    stop("D_G must be a numeric matrix or a Matrix object.", call. = FALSE)
  }
  values <- if (inherits(D_G, "sparseMatrix")) {
    if ("x" %in% methods::slotNames(D_G)) D_G@x else numeric(0)
  } else {
    as.matrix(D_G)
  }
  if (is.complex(values) || any(!is.finite(values))) {
    stop("D_G must have finite entries.", call. = FALSE)
  }
  invisible(D_G)
}

# Explicit seeds are local to the operation; seed = NULL uses the caller's RNG.
.gfm_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  .gfm_check_scalar(seed, "seed", lower = 0, upper = .Machine$integer.max,
                    integer = TRUE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

# Stable log probabilities, including the zero reference logit.
# Internal callers have already checked the matrix dimensions and entries.
.gfm_log_softmax <- function(beta) {
  row_max <- pmax(0, apply(beta, 1L, max))
  shifted <- sweep(beta, 1L, row_max, "-")
  log_denom <- log(rowSums(exp(shifted)) + exp(-row_max))
  cbind(sweep(shifted, 1L, log_denom, "-"), -row_max - log_denom)
}

#' Convert reference-category logits to multinomial probabilities.
#' @param beta Finite n by (D - 1) numeric matrix.
#' @return An n by D matrix whose rows sum to one.
softmax_probs <- function(beta) {
  .gfm_check_matrix(beta, "beta")
  exp(.gfm_log_softmax(beta))
}

#' Convert reference-category logits to log probabilities.
#' @param beta Finite n by (D - 1) numeric matrix.
#' @return An n by D matrix, without clipping small probabilities.
log_softmax_probs <- function(beta) {
  .gfm_check_matrix(beta, "beta")
  .gfm_log_softmax(beta)
}

#' Project each entry onto the interval [-lambda, lambda].
#' @param x Finite numeric vector or matrix.
#' @param lambda Nonnegative scalar box radius.
#' @return The projected object, preserving its dimensions.
proj_box <- function(x, lambda) {
  .gfm_check_scalar(lambda, "lambda", lower = 0)
  if (!is.numeric(x) || is.complex(x) || any(!is.finite(x))) {
    stop("x must contain finite numeric values.", call. = FALSE)
  }
  pmin(pmax(x, -lambda), lambda)
}

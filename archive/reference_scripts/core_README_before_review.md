# Graph-Fused Multinomial: cleaned R scripts

These scripts retain the existing function names, argument order, default
settings, and original output fields. They use the same reference-category
parameterization, graph-fused objective, primal-dual update order, and
within-node trial-splitting cross-validation scheme. Corrections that can affect
numerical results or random-number streams are listed below.

## Setup and use

Install the three computational dependencies once:

```r
install.packages(c("igraph", "Matrix", "FNN"))
```

Keep the six scripts together in the `R` directory and load them with:

```r
source("R/load_all.R")
```

An absolute path also works from any working directory. The loader checks
dependencies and sources functions into the calling environment. It does not
install or attach packages. Add `library(ggplot2)`, `library(reshape2)`, or other
attachments explicitly in analysis scripts that use their unqualified functions.
Plotting packages and CARBayes are not dependencies of these six scripts.

For manual sourcing, load `model_utils.R` first, then `graph_utils.R`, followed
by the remaining function files. Shared checks are private `.gfm_*` helpers.
Using `load_all.R` is recommended.

```r
source("R/load_all.R")
set.seed(1)

g <- build_grid_graph(5, 5)
D_G <- build_incidence_matrix(g)
truth <- generate_true_beta_X(g, D = 4, nrow_grid = 5, ncol_grid = 5)
Y <- simulate_data(truth$P_true, n_i = rep(100, igraph::vcount(g)))

fit <- fit_graph_fused_multinomial(
  Y, D_G, lambda = 2.5, T_max = 8000, verbose = FALSE
)
fit$converged
fit$P_hat

cv <- cv_graph_fused_multinomial(
  Y, D_G, lambda_grid = c(0.5, 1, 2.5, 5),
  K = 5, seed = 1, T_max = 8000, verbose = FALSE
)
cv$lambda_min
cv$lambda_1se

# Refit on the full counts after selecting lambda.
final_fit <- fit_graph_fused_multinomial(
  Y, D_G, lambda = cv$lambda_min, T_max = 8000, verbose = FALSE
)
```

The objective uses the total negative log-likelihood, not an average over
trials. Lambda therefore has the same scale as in the original scripts. CV
returns a selected candidate; it does not rescale lambda for the final refit.

## Files and public functions

| File | Responsibility |
| --- | --- |
| `R/model_utils.R` | Stable `softmax_probs()`, new `log_softmax_probs()`, `proj_box()`, and shared private checks |
| `R/graph_utils.R` | Rook/king grids, union kNN, repaired SBM, sparse incidence matrices, and norm diagnostics |
| `R/simulation_utils.R` | Quadrant, X-shaped, and block truths; independent multinomial draws |
| `R/fit_graph_fused_multinomial.R` | Primal-dual estimator and convergence diagnostics |
| `R/cross_validation.R` | Count splitting, held-out scoring, minimum-CV and one-SE selection |
| `R/load_all.R` | Dependency checking and location-independent sourcing |

Each public function has comments describing its arguments and return values.
Source files use consistent indentation, LF line endings, explicit namespaces,
dimension-preserving subsetting, and informative input errors.

## Corrections and compatibility notes

1. **Stable probabilities and scores.** Softmax calculations subtract the
   largest row logit, including the zero reference logit. CV scores use stable
   log probabilities instead of clipping probabilities at `1e-10`. Scores can
   therefore change when fitted probabilities were previously clipped.
2. **Guaranteed step-size check.** The solver uses
   `op2 = norm(D_G, "1") * norm(D_G, "I")`, an upper bound for the squared
   operator norm. For an ordinary incidence matrix this is twice the maximum
   vertex degree. The original power-iteration estimate can underestimate the
   norm. The new bound can make default steps more conservative and can reject
   user-supplied steps that passed the former estimated check. The `op2` fit
   field now reports this bound. `op_norm_sq_DG()` remains an estimate for
   diagnostics; its initialization preserves the caller's random state.
3. **Zero-edge and zero-count cases.** Edgeless graphs and zero-total count rows
   have finite default step sizes. A completely zero-count dataset returns the
   initialized uniform probabilities; no unique probability estimate is
   implied when there is no likelihood information.
4. **CV indexing and empty folds.** Unsorted and repeated lambda candidates
   retain their original column positions. Empty validation folds have raw
   NLL zero and rate `NA`; they are omitted from the CV mean and SE. At least
   two nonempty validation folds are required. Node-level zero counts never
   require removing graph vertices.
5. **Explicit CV score definition.** The original equal-fold convention is
   retained: each active fold's NLL is divided by its total validation count,
   then rates are averaged without trial-count weighting. `cv_se` is their SD
   divided by the square root of the active-fold count. This is negative
   log-likelihood, rather than deviance. The one-SE rule is a selection
   heuristic, not a confidence interval, because training folds overlap.
6. **Random-number behavior.** A supplied splitting seed is local to that
   split and does not reset or advance the caller's random stream. `seed = NULL`
   uses and advances that stream. The solver no longer consumes random numbers
   to estimate its operator norm. Future simulated datasets in an old combined
   simulation/fitting loop can consequently differ even with the same initial
   seed; fitted objectives and the truth distributions are unchanged.
7. **Graph repairs.** One-dimensional king grids no longer create invalid
   diagonal indices. kNN queries explicitly exclude the node itself, including
   at duplicate coordinates. SBM repair avoids R's length-one `sample()` trap
   and rejects the impossible one-vertex no-isolation case. A repaired SBM is
   not an unmodified draw from the original SBM distribution.
8. **Simulation edge cases.** Constant grid axes use a zero standardized trend
   and a normalized midpoint of `0.5`. The block generator constructs logits
   directly to avoid rounding a perturbed category probability to one before
   taking its log odds. Ordinary nondegenerate truth settings retain their
   definitions, up to floating-point rounding.
9. **Diagnostics.** Fits add `converged` and `residual`. CV adds `cv_fold_rate`,
   `active_folds`, `n_iter`, and `converged`, and warns once if any fitted path
   reaches its iteration limit. The stopping rule remains relative changes in
   primal and dual iterates; it is not a primal-dual gap certificate. Cases
   without a finite logit minimizer can still reach the iteration limit.

## Validation

Run the included checks with:

```bash
Rscript tests/validate.R
```

The checks cover extreme-logit stability; graph edge counts and incidence
properties; duplicate-location kNN; SBM isolation repair; truth generators and
count totals; analytic two-node binary penalized solutions; edgeless and
zero-count fits; fold reconstruction and RNG preservation; closed-form
multinomial CV scores; empty folds; invalid inputs; and portable loading.

All 13 validation groups passed on 2026-10-09 with R 4.3.3, Matrix 1.6.5,
igraph 1.6.0, and FNN 1.1.4. A separate comparison confirmed that all 16
existing public function signatures and defaults were retained. Six ordinary
grid fits, using D in {2, 4, 8} and lambda in {1, 2.5}, agreed with the original
implementation to a maximum probability difference of 8.916e-9 after both
solvers met a 1e-9 residual tolerance. A three-fold binary CV comparison gave
identical mean scores at displayed precision and the same selected lambdas.
These comparisons used the same observed count matrices and matching split
seeds; they do not imply identical future random streams in an old simulation
loop.

These are numerical and integration checks, not a proof of solver convergence
or a rerun of the full GFM/MCAR Monte Carlo study.

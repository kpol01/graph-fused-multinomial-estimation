# Run this script from the repository root (for example, from the RStudio project).
source(file.path("R", "load_all.R"))

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



## Save the table for use in the paper.
write.csv(results_print, file.path("results", "oracle_vs_cv_table.csv"), row.names = FALSE)

if (requireNamespace("xtable", quietly = TRUE)) {
  tex_tab <- xtable::xtable(
    results_print[, c("n_actual", "n_edges", "lambda_oracle", "mse_oracle",
                      "lambda_cv", "mse_cv", "mse_diff", "pct_increase")],
    caption = "Oracle vs. cross-validated lambda: MSE(P) comparison on the grid graph across sample sizes.",
    label = "tab:oracle_vs_cv",
    digits = c(0, 0, 0, 2, 5, 2, 5, 5, 2)
  )
  tex_path <- file.path("results", "oracle_vs_cv_table.tex")
  sink(tex_path)
  print(tex_tab, include.rownames = FALSE)
  sink()
  cat(sprintf("\nLaTeX table written to %s\n", tex_path))
} else {
  cat("\n(Install the 'xtable' package to also generate a ready-to-use LaTeX table.)\n")
}

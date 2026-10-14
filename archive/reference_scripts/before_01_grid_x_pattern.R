# Run this script from the repository root (for example, from the RStudio project).
source(file.path("R", "load_all.R"))

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

# Optional saved outputs
ggsave(file.path("figures", "grid_convergence.png"), p_conv, width = 7, height = 5, dpi = 300)
ggsave(file.path("figures", "grid_mse_vs_lambda.png"), p_lambda, width = 7, height = 5, dpi = 300)
ggsave(file.path("figures", "grid_x_recovery.png"), p_map, width = 10, height = 4, dpi = 300)
ggsave(file.path("figures", "grid_cv_curve.png"), p_cv, width = 7, height = 5, dpi = 300)

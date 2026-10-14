# Run this script from the repository root (for example, from the RStudio project).
source(file.path("R", "load_all.R"))

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

## --- random locations and two candidate k-NN graphs (small vs large k) --
D <- 4
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

# Optional saved outputs
ggsave(file.path("figures", "knn_graph.png"), p_graph, width = 7, height = 6, dpi = 300)
ggsave(file.path("figures", "knn_x_recovery.png"), p_map_knn, width = 10, height = 4, dpi = 300)

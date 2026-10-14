# Run this script from the repository root (for example, from the RStudio project).
source(file.path("R", "load_all.R"))

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


D <- 4

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

# Optional saved outputs
ggsave(file.path("figures", "sbm_community_structure.png"), p_comm, width = 7, height = 6, dpi = 300)
ggsave(file.path("figures", "sbm_recovery.png"), p_map_gen, width = 10, height = 4, dpi = 300)

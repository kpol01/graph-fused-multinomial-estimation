# Run this script from the repository root (for example, from the RStudio project).
setwd("C:/Users/KAUSTAV/Box/Graph Fused Multinomial")
source(file.path("R", "load_all.R"))

library(parallel)
library(doParallel)
library(foreach)
library(tidyr)

n_pt <- 20
D <- 4
mse <- c()
t <- c()

side  <- n_pt                            # NOTE: fixed below, see caveat
g_bc  <- build_grid_graph(side, side)
n_bc  <- igraph::vcount(g_bc)
D_G_bc <- build_incidence_matrix(g_bc)
W_bc  <- as.matrix(igraph::as_adjacency_matrix(g_bc))

truth_bc <- generate_true_beta_X(g_bc, D, side, side,
                                 p_high = 0.70, p_low = 0.08, half_width = 0.06)
P_true_bc <- truth_bc$P_true

n_i_bc <- sample(150:300, n_bc, replace = TRUE)
Y_bc <- simulate_data(P_true_bc, n_i_bc)

t_ours <- system.time({
  fit_ours_bc <- fit_graph_fused_multinomial(Y_bc, D_G_bc, lambda = 2.5,
                                             T_max = 5000, tol = 1e-6, verbose = FALSE)
})


t_car <- system.time({
  fit_car_bc <- CARBayes::MVS.CARleroux(
    formula  = Y_bc ~ 1,
    family   = "multinomial",
    trials   = n_i_bc,
    W        = W_bc,
    burnin   = 5000,
    n.sample = 20000,
    thin     = 10,
    verbose  = FALSE
  )
})


library(ggplot2)
library(dplyr)

## ============================================================
## 1. Extract true and fitted probabilities
## ============================================================

P_true <- as.matrix(P_true_bc)

## Raw multinomial proportions = MLE
P_mle <- sweep(
  as.matrix(Y_bc),
  MARGIN = 1,
  STATS = n_i_bc,
  FUN = "/"
)

## GFM estimated probabilities
P_gfm <- as.matrix(fit_ours_bc$P_hat)

## MCAR fitted values
## Convert fitted counts to probabilities
P_mcar <- sweep(
  as.matrix(fit_car_bc$fitted.values),
  MARGIN = 1,
  STATS = n_i_bc,
  FUN = "/"
)

## ============================================================
## 2. Check dimensions
## ============================================================

print(dim(P_true))
print(dim(P_mle))
print(dim(P_gfm))
print(dim(P_mcar))

stopifnot(
  all(dim(P_true) == dim(P_mle)),
  all(dim(P_true) == dim(P_gfm)),
  all(dim(P_true) == dim(P_mcar))
)

n_plot <- nrow(P_true)
D_plot <- ncol(P_true)

## ============================================================
## 3. Compute MSE values
## ============================================================

RE_mle  <- sqrt(mean((P_mle  - P_true)^2)/mean(P_true^2))
RE_gfm  <- sqrt(mean((P_gfm  - P_true)^2)/mean(P_true^2))
RE_mcar <- sqrt(mean((P_mcar - P_true)^2)/mean(P_true^2))

cat("MLE RE  =", RE_mle,  "\n")
cat("GFM RE  =", RE_gfm,  "\n")
cat("MCAR RE =", RE_mcar, "\n")

## ============================================================
## 4. Category names
## ============================================================

true_names <- colnames(P_true)

if (!is.null(true_names) &&
    length(true_names) == D_plot &&
    !anyDuplicated(true_names)) {
  
  category_names <- true_names
  
} else {
  
  category_names <- paste0("Category ", seq_len(D_plot))
}

print(category_names)

## ============================================================
## 5. Column labels with MSE
## ============================================================

true_label <- "True"

mle_label <- paste0(
  "MLE\n(Error = ",
  formatC(RE_mle*100, format = "f", digits = 2),
  "%)"
)

gfm_label <- paste0(
  "GFM\n(Error = ",
  formatC(RE_gfm*100, format = "f", digits = 2),
  "%)"
)

mcar_label <- paste0(
  "MCAR\n(Error = ",
  formatC(RE_mcar*100, format = "f", digits = 2),
  "%)"
)

method_levels <- c(
  true_label,
  mle_label,
  gfm_label,
  mcar_label
)

## ============================================================
## 6. Grid coordinates
## ============================================================

grid_coords <- data.frame(
  vertex = seq_len(n_plot),
  x = (seq_len(n_plot) - 1) %% side + 1,
  y = (seq_len(n_plot) - 1) %/% side + 1
)

## ============================================================
## 7. Convert matrices to long format
## ============================================================

make_plot_df <- function(P, method_name) {
  
  P <- as.matrix(P)
  
  stopifnot(
    nrow(P) == n_plot,
    ncol(P) == D_plot
  )
  
  out <- data.frame(
    vertex = rep(seq_len(n_plot), times = D_plot),
    Category = rep(category_names, each = n_plot),
    Probability = as.vector(P),
    Method = method_name,
    stringsAsFactors = FALSE
  )
  
  left_join(
    out,
    grid_coords,
    by = "vertex"
  )
}

## ============================================================
## 8. Build plotting data
## ============================================================

plot_true <- make_plot_df(
  P_true,
  true_label
)

plot_mle <- make_plot_df(
  P_mle,
  mle_label
)

plot_gfm <- make_plot_df(
  P_gfm,
  gfm_label
)

plot_mcar <- make_plot_df(
  P_mcar,
  mcar_label
)

plot_df <- bind_rows(
  plot_true,
  plot_mle,
  plot_gfm,
  plot_mcar
)

## ============================================================
## 9. Set facet ordering
## ============================================================

plot_df$Method <- factor(
  plot_df$Method,
  levels = method_levels
)

plot_df$Category <- factor(
  plot_df$Category,
  levels = category_names
)

## ============================================================
## 10. Plot
## ============================================================

p <- ggplot(
  plot_df,
  aes(
    x = x,
    y = y,
    fill = Probability
  )
) +
  
  geom_tile() +
  
  ## Rows = categories
  ## Columns = True, MLE, GFM, MCAR
  facet_grid(
    Category ~ Method,
    switch = "y"
  ) +
  
  ## Common scale across every panel
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1),
    name = "Probability / Proportion",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  
  coord_equal() +
  
  scale_y_reverse() +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(base_size = 11) +
  
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    
    panel.grid = element_blank(),
    
    ## Column headings
    strip.text.x = element_text(
      face = "bold",
      size = 10,
      lineheight = 1
    ),
    
    ## Category names on left
    strip.text.y.left = element_text(
      size = 9,
      angle = 0
    ),
    
    strip.placement = "outside",
    strip.background = element_blank(),
    
    panel.spacing = grid::unit(0.15, "cm"),
    
    legend.position = "bottom"
  ) +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(10, "cm"),
      barheight = unit(0.4, "cm")
    )
  )

p

ggsave(
  "figures/comparison/per_sig_plot.pdf",
  plot = p,
  width = 7.2,
  height = 6.4 ,
  units = "in"
)


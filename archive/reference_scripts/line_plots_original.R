setwd("C:/Users/KAUSTAV/Box/Graph Fused Multinomial")
source(file.path("R", "load_all.R"))

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

mse_grid_n <- read.csv("C:/Users/KAUSTAV/Box/Graph Fused Multinomial/results/mse_grid_n.csv")
mse_knn_n <- read.csv("C:/Users/KAUSTAV/Box/Graph Fused Multinomial/results/mse_knn_n.csv")
mse_sbm_n <- read.csv("C:/Users/KAUSTAV/Box/Graph Fused Multinomial/results/mse_sbm_n.csv")

mse_grid_n[,1] <- c(49, 100, 196, 289, 484)

mse_grid_n[,-1] <- sqrt(mse_grid_n[,-1]/0.109168)
mse_knn_n[,-1] <- sqrt(mse_knn_n[,-1]/0.1120613)
mse_sbm_n[,-1] <- sqrt(mse_sbm_n[,-1]/0.08016409)

mse_grid_n[,c(3,5,7)] <- mse_grid_n[,c(3,5,7)]/40
mse_knn_n[,c(3,5,7)] <- mse_knn_n[,c(3,5,7)]/40
mse_sbm_n[,c(3,5,7)] <- mse_sbm_n[,c(3,5,7)]/40
#------------------------------------------------------------
# Function to convert each MSE data frame to long format
#------------------------------------------------------------

prepare_mse <- function(dat) {
  
  # MSE estimates
  mse_dat <- dat %>%
    select(n, gfm, mcar, mmle) %>%
    pivot_longer(
      cols = -n,
      names_to = "method",
      values_to = "mse"
    )
  
  # Standard deviations / sigma
  sig_dat <- dat %>%
    select(n, gfm_sig, mcar_sig, mmle_sig) %>%
    rename(
      gfm  = gfm_sig,
      mcar = mcar_sig,
      mmle = mmle_sig
    ) %>%
    pivot_longer(
      cols = -n,
      names_to = "method",
      values_to = "sigma"
    )
  
  left_join(mse_dat, sig_dat, by = c("n", "method")) %>%
    mutate(
      method = factor(
        method,
        levels = c("gfm", "mcar", "mmle"),
        labels = c("GFM", "MCAR", "MLE")
      )
    )
}

add_x_shift <- function(dat, shift = 9) {
  dat %>%
    mutate(
      x_plot = n + case_when(
        method == "GFM"  ~ -shift,
        method == "MCAR" ~  0,
        method == "MLE" ~  shift
      )
    )
}




#------------------------------------------------------------
# Prepare the three data sets
#------------------------------------------------------------

grid_dat <- prepare_mse(mse_grid_n)
knn_dat  <- prepare_mse(mse_knn_n)
sbm_dat  <- prepare_mse(mse_sbm_n)


grid_dat <- add_x_shift(grid_dat)
knn_dat  <- add_x_shift(knn_dat)
sbm_dat  <- add_x_shift(sbm_dat)

#------------------------------------------------------------
# Plotting function
#------------------------------------------------------------

#------------------------------------------------------------
# Common y-axis limits for all three plots
#------------------------------------------------------------

y_max <- max(
  grid_dat$mse + grid_dat$sigma,
  knn_dat$mse  + knn_dat$sigma,
  sbm_dat$mse  + sbm_dat$sigma,
  na.rm = TRUE
)

# Add a little space above the highest error bar
y_lim <- c(0, 1.05 * y_max)

#------------------------------------------------------------
# Plotting function
#------------------------------------------------------------

make_mse_plot <- function(dat, title, show_legend = FALSE) {
  
  cap_width <- 6
  
  p <- ggplot(
    dat,
    aes(
      x = x_plot,
      y = mse,
      color = method,
      shape = method,
      group = method
    )
  ) +
    
    # Vertical error bars
    geom_linerange(
      aes(
        ymin = pmax(mse - sigma, 0),
        ymax = mse + sigma
      ),
      linewidth = 0.4
    ) +
    
    # Bottom horizontal caps
    geom_segment(
      aes(
        x = x_plot - cap_width,
        xend = x_plot + cap_width,
        y = pmax(mse - sigma, 0),
        yend = pmax(mse - sigma, 0)
      ),
      linewidth = 0.4
    ) +
    
    # Top horizontal caps
    geom_segment(
      aes(
        x = x_plot - cap_width,
        xend = x_plot + cap_width,
        y = mse + sigma,
        yend = mse + sigma
      ),
      linewidth = 0.4
    ) +
    
    # Lines connecting each method
    geom_line(
      linewidth = 0.6
    ) +
    
    geom_point(
      size = 2.8
    ) +
    
    scale_color_manual(
      values = c(
        "GFM"  = "#0072B2",
        "MCAR" = "#D55E00",
        "MLE" = "#009E73"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        "GFM"  = 16,
        "MCAR" = 17,
        "MLE" = 15
      )
    ) +
    
    scale_y_continuous(
      limits = y_lim,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    # Show the original n values on the axis
    scale_x_continuous(
      breaks = sort(unique(dat$n)),
      labels = sort(unique(dat$n))
    ) +
    
    labs(
      title = title,
      x = "Sample size(n)",
      y = "Relative Error",
      color = NULL,
      shape = NULL
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 13
      ),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        linewidth = 0.25,
        color = "grey85"
      )
    )
  
  if (show_legend) {
    p <- p +
      theme(
        legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "center"
      ) +
      guides(
        color = guide_legend(ncol = 1),
        shape = guide_legend(ncol = 1)
      )
  } else {
    p <- p +
      theme(legend.position = "none")
  }
  
  return(p)
}

#------------------------------------------------------------
# Make plots
#------------------------------------------------------------

p_grid <- make_mse_plot(

    grid_dat,
  "Grid graph",
  show_legend = FALSE
)

p_knn <- make_mse_plot(
  knn_dat,
  "kNN graph",
  show_legend = FALSE
)

p_sbm <- make_mse_plot(
  sbm_dat,
  "SBM graph",
  show_legend = TRUE
)



design <- "
AABB
#CC#
"

p_final <- p_grid + p_knn + p_sbm +
  plot_layout(design = design)

p_final

ggsave(
  "figures/Comparison/mse_vs_n.pdf",
  p_final,
  width = 17,
  height = 13.5,
  units = "cm"
)


# ===========================================================================================

t_grid_n <- read.csv("C:/Users/KAUSTAV/Box/Graph Fused Multinomial/results/t_grid_n.csv")
t_knn_n <- read.csv("C:/Users/KAUSTAV/Box/Graph Fused Multinomial/results/t_knn_n.csv")
t_sbm_n <- read.csv("C:/Users/KAUSTAV/Box/Graph Fused Multinomial/results/t_sbm_n.csv")


t_grid_n[,1] <- c(49, 100, 196, 289, 484)
t_grid_n[,3] <- t_grid_n[,3]/t_grid_n[,2]
t_grid_n[,5] <- t_grid_n[,5]/t_grid_n[,4]
t_grid_n[,7] <- t_grid_n[,7]/t_grid_n[,6]
t_grid_n[,2] <- log(t_grid_n[,2])
t_grid_n[,4] <- log(t_grid_n[,4])
t_grid_n[,6] <- log(t_grid_n[,6])




t_knn_n[,3] <- t_knn_n[,3]/t_knn_n[,2]
t_knn_n[,5] <- t_knn_n[,5]/t_knn_n[,4]
t_knn_n[,7] <- t_knn_n[,7]/t_knn_n[,6]
t_knn_n[,2] <- log(t_knn_n[,2])
t_knn_n[,4] <- log(t_knn_n[,4])
t_knn_n[,6] <- log(t_knn_n[,6])



t_sbm_n[,3] <- t_sbm_n[,3]/t_sbm_n[,2]
t_sbm_n[,5] <- t_sbm_n[,5]/t_sbm_n[,4]
t_sbm_n[,7] <- t_sbm_n[,7]/t_sbm_n[,6]
t_sbm_n[,2] <- log(t_sbm_n[,2])
t_sbm_n[,4] <- log(t_sbm_n[,4])
t_sbm_n[,6] <- log(t_sbm_n[,6])







#------------------------------------------------------------
# Function to convert each MSE data frame to long format
#------------------------------------------------------------

prepare_mse <- function(dat) {
  
  # MSE estimates
  mse_dat <- dat %>%
    select(n, gfm, mcar, mmle) %>%
    pivot_longer(
      cols = -n,
      names_to = "method",
      values_to = "mse"
    )
  
  # Standard deviations / sigma
  sig_dat <- dat %>%
    select(n, gfm_sig, mcar_sig, mmle_sig) %>%
    rename(
      gfm  = gfm_sig,
      mcar = mcar_sig,
      mmle = mmle_sig
    ) %>%
    pivot_longer(
      cols = -n,
      names_to = "method",
      values_to = "sigma"
    )
  
  left_join(mse_dat, sig_dat, by = c("n", "method")) %>%
    mutate(
      method = factor(
        method,
        levels = c("gfm", "mcar", "mmle"),
        labels = c("GFM", "MCAR", "MMLE")
      )
    )
}


#------------------------------------------------------------
# Prepare the three data sets
#------------------------------------------------------------

grid_dat <- prepare_mse(t_grid_n)
knn_dat  <- prepare_mse(t_knn_n)
sbm_dat  <- prepare_mse(t_sbm_n)


#------------------------------------------------------------
# Plotting function
#------------------------------------------------------------

#------------------------------------------------------------
# Common y-axis limits for all three plots
#------------------------------------------------------------

y_max <- max(
  grid_dat$mse + grid_dat$sigma,
  knn_dat$mse  + knn_dat$sigma,
  sbm_dat$mse  + sbm_dat$sigma,
  na.rm = TRUE
)

# Add a little space above the highest error bar
y_lim <- c(-2, 1.05 * y_max)


#------------------------------------------------------------
# Plotting function
#------------------------------------------------------------

make_mse_plot <- function(dat, title, show_legend = FALSE) {
  
  p <- ggplot(
    dat,
    aes(
      x = n,
      y = mse,
      color = method,
      shape = method,
      group = method
    )
  ) +
    
    # ± 1 sigma error bars
    geom_errorbar(
      aes(
        ymin = mse - 6*sigma,
        ymax = mse + 6*sigma
      ),
      width = 15,
      linewidth = 0.4
    ) +
    
    geom_line(
      linewidth = 0.6
    ) +
    
    geom_point(
      size = 2.8
    ) +
    
    scale_color_manual(
      values = c(
        "GFM"  = "#0072B2",
        "MCAR" = "#D55E00",
        "MMLE" = "#009E73"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        "GFM"  = 16,   # circle
        "MCAR" = 17,   # triangle
        "MMLE" = 15    # square
      )
    ) +
    
    # Same y-axis limits for every panel
    scale_y_continuous(
      limits = y_lim,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    labs(
      title = title,
      x = expression(n),
      y = "MSE",
      color = NULL,
      shape = NULL
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 13
      ),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        linewidth = 0.25,
        color = "grey85"
      )
    )
  
  # Legend only for SBM plot
  if (show_legend) {
    
    p <- p +
      theme(
        legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "center"
      ) +
      guides(
        color = guide_legend(ncol = 1),
        shape = guide_legend(ncol = 1)
      )
    
  } else {
    
    p <- p +
      theme(
        legend.position = "none"
      )
  }
  
  return(p)
}


#------------------------------------------------------------
# Make plots
#------------------------------------------------------------

p_grid <- make_mse_plot(
  grid_dat,
  "Grid graph",
  show_legend = FALSE
)

p_knn <- make_mse_plot(
  knn_dat,
  "kNN graph",
  show_legend = FALSE
)

p_sbm <- make_mse_plot(
  sbm_dat,
  "SBM graph",
  show_legend = TRUE
)



design <- "
AABB
#CC#
"

p_final <-  p_knn + p_grid + p_sbm +
  plot_layout(design = design)

p_final

ggsave(
  "figures/Comparison/mse_vs_n.pdf",
  p_final,
  width = 17,
  height = 13.5,
  units = "cm"
)


#===============================================================================


mse_grid_D <- read.csv("C:/Users/KAUSTAV/Box/Graph Fused Multinomial/results/mse_grid_n_D.csv")
mse_knn_D <- read.csv("C:/Users/KAUSTAV/Box/Graph Fused Multinomial/results/mse_knn_n_D.csv")
mse_sbm_D <- read.csv("C:/Users/KAUSTAV/Box/Graph Fused Multinomial/results/mse_sbm_n_D.csv")

mse_grid_D[,1] <- c(5, 75, 150, 300, 500)
mse_knn_D[,1] <- c(5, 75, 150, 300, 500)
mse_sbm_D[,1] <- c(5, 75, 150, 300, 500)

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

#------------------------------------------------------------
# Data preparation
#------------------------------------------------------------

prepare_mse <- function(dat) {
  
  # MSE estimates
  mse_dat <- dat %>%
    select(n, gfm, mcar, mmle) %>%
    pivot_longer(
      cols = -n,
      names_to = "method",
      values_to = "mse"
    )
  
  # Standard deviations / sigma
  sig_dat <- dat %>%
    select(n, gfm_sig, mcar_sig, mmle_sig) %>%
    rename(
      gfm  = gfm_sig,
      mcar = mcar_sig,
      mmle = mmle_sig
    ) %>%
    pivot_longer(
      cols = -n,
      names_to = "method",
      values_to = "sigma"
    )
  
  left_join(mse_dat, sig_dat, by = c("n", "method")) %>%
    mutate(
      method = factor(
        method,
        levels = c("gfm", "mcar", "mmle"),
        labels = c("GFM", "MCAR", "MLE")
      )
    )
}

add_x_shift <- function(dat, shift = 9) {
  dat %>%
    mutate(
      x_plot = n + case_when(
        method == "GFM"  ~ -shift,
        method == "MCAR" ~  0,
        method == "MLE"  ~  shift
      )
    )
}

grid_dat <- add_x_shift(prepare_mse(mse_grid_D))
knn_dat  <- add_x_shift(prepare_mse(mse_knn_D))
sbm_dat  <- add_x_shift(prepare_mse(mse_sbm_D))

#------------------------------------------------------------
# Common y-axis limits for all three plots
#------------------------------------------------------------

y_max <- max(
  grid_dat$mse + grid_dat$sigma,
  knn_dat$mse  + knn_dat$sigma,
  sbm_dat$mse  + sbm_dat$sigma,
  na.rm = TRUE
)

y_lim <- c(0, 1.05 * y_max)

#------------------------------------------------------------
# Plotting function (identical legend setup for every plot)
#------------------------------------------------------------

make_mse_plot <- function(dat, title) {
  
  cap_width <- 6
  
  ggplot(
    dat,
    aes(
      x = x_plot,
      y = mse,
      color = method,
      shape = method,
      group = method
    )
  ) +
    
    # Vertical error bars
    geom_linerange(
      aes(
        ymin = pmax(mse - 6*sigma, 0),
        ymax = mse + 6*sigma
      ),
      linewidth = 0.4
    ) +
    
    # Bottom horizontal caps
    geom_segment(
      aes(
        x = x_plot - cap_width,
        xend = x_plot + cap_width,
        y = pmax(mse -  6*sigma, 0),
        yend = pmax(mse - 6*sigma, 0)
      ),
      linewidth = 0.4
    ) +
    
    # Top horizontal caps
    geom_segment(
      aes(
        x = x_plot - cap_width,
        xend = x_plot + cap_width,
        y = mse + 6*sigma,
        yend = mse + 6*sigma
      ),
      linewidth = 0.4
    ) +
    
    # Lines connecting each method
    geom_line(linewidth = 0.4) +
    geom_point(size = 1.8) +
    
    scale_color_manual(
      values = c(
        "GFM"  = "#0072B2",
        "MCAR" = "#D55E00",
        "MLE"  = "#009E73"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        "GFM"  = 16,
        "MCAR" = 17,
        "MLE"  = 15
      )
    ) +
    
    scale_y_continuous(
      limits = y_lim,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    # Show the original D values on the axis
    scale_x_continuous(
      breaks = sort(unique(dat$n)),
      labels = sort(unique(dat$n))
    ) +
    
    labs(
      title = title,
      x = "Sample Size(n)",
      y = "Relative Error",
      color = NULL,
      shape = NULL
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 13
      ),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        linewidth = 0.25,
        color = "grey85"
      ),
      legend.position  = "bottom",
      legend.direction = "horizontal"
    )
}

#------------------------------------------------------------
# Make plots
#------------------------------------------------------------

p_grid <- make_mse_plot(grid_dat, "Grid graph")
p_knn  <- make_mse_plot(knn_dat,  "kNN graph")
p_sbm  <- make_mse_plot(sbm_dat,  "SBM graph")

#------------------------------------------------------------
# Layout: two plots on top, one centered below,
# single horizontal legend at the bottom
#------------------------------------------------------------

top    <- p_sbm + p_grid
bottom <- plot_spacer() + p_knn + plot_spacer() +
  plot_layout(widths = c(1, 2, 1))

p_final <- (top / bottom) +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      legend.position  = "bottom",
      legend.direction = "horizontal"
    )
  )

p_final

ggsave(
  "figures/Comparison/mse_vs_n_D.pdf",
  p_final,
  width = 17,
  height = 15,
  units = "cm"
)

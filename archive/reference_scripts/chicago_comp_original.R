# Install once
install.packages(c("sf", "Matrix", "igraph", "httr2"))
setwd("C:/Users/KAUSTAV/Box/Graph Fused Multinomial")
source(file.path("R", "load_all.R"))
# Select the downloaded script
source("C:/Users/KAUSTAV/Box/Graph Fused Multinomial/Real Data Analysis/prepare_chicago_crimes.R")

chicago <- prepare_chicago_crimes(
  year = 2024,
  cellsize_m = 1000,
  min_count = 1,
  output_dir = "chicago_crime_data"
)

# Inputs for your existing code
Y_train <- chicago$Y_train
Y_test  <- chicago$Y_test
D_G     <- chicago$D_G
n      <- chicago$n
D      <- chicago$D
g      <- chicago$g
W      <- as.matrix(chicago$W)
coords <- chicago$coords
n_i_train <- chicago$n_i_train
# Inspect the prepared data
chicago$prep_info

# fit GFM
t_ours <- system.time({
  fit_chicago <- fit_graph_fused_multinomial(Y_train, D_G, lambda = 1,
                                              T_max = 5000, tol = 1e-6, verbose = FALSE)
})
fit_chicago
P_hat <- fit_chicago$P_hat


fit_chicago_prop <- fit_graph_fused_multinomial(Y_test, D_G, lambda = 0,
                                                T_max = 5000, tol = 1e-6, verbose = FALSE)
P_hat_prop <- fit_chicago_prop$P_hat


# fit MCAR
t_car <- system.time({
  fit_car_gen<- CARBayes::MVS.CARleroux(
    formula  = Y_train ~ 1,
    family   = "multinomial",
    trials   = n_i_train,
    W        = W,
    burnin   = 5000,
    n.sample = 20000,
    thin     = 10,
    verbose  = FALSE
  )
})

fitted_counts_car <- as.matrix(fit_car_gen$fitted.values)
stopifnot(all(dim(fitted_counts_car) == dim(Y_train)))
P_hat_car <- fitted_counts_car / rowSums(fitted_counts_car)

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ggh4x)

## ============================================================
## Observed vs Graph-Fused Multinomial vs MCAR
## All crime categories:
##   rows    = crime categories
##   columns = Observed, Graph-Fused Multinomial, MCAR
## ============================================================

if (is.null(chicago$Y_test)) {
  stop("chicago$Y_test is missing. First prepare or split the Chicago counts.")
}

Y_test <- as.matrix(chicago$Y_test)
P_ours <- as.matrix(P_hat)
P_mcar <- as.matrix(P_hat_car)

## ------------------------------------------------------------
## Checks
## ------------------------------------------------------------

stopifnot(
  nrow(Y_test) == nrow(chicago$grid_sf),
  ncol(Y_test) == ncol(chicago$Y),
  identical(dim(P_ours), dim(Y_test)),
  identical(dim(P_mcar), dim(Y_test)),
  is.numeric(Y_test),
  all(is.finite(Y_test)),
  all(Y_test >= 0),
  all(Y_test == floor(Y_test))
)

stopifnot(
  is.numeric(P_ours),
  is.numeric(P_mcar),
  all(is.finite(P_ours)),
  all(is.finite(P_mcar)),
  all(P_ours >= 0 & P_ours <= 1),
  all(P_mcar >= 0 & P_mcar <= 1),
  all(abs(rowSums(P_ours) - 1) < 1e-6),
  all(abs(rowSums(P_mcar) - 1) < 1e-6)
)

## ------------------------------------------------------------
## Observed test-set proportions
## ------------------------------------------------------------

n_test <- rowSums(Y_test)
has_test <- n_test > 0

if (!any(has_test)) {
  stop("There are no test crimes to compare with the fits.")
}

P_obs <- matrix(
  NA_real_,
  nrow = nrow(Y_test),
  ncol = ncol(Y_test),
  dimnames = dimnames(Y_test)
)

P_obs[has_test, ] <-
  Y_test[has_test, , drop = FALSE] /
  n_test[has_test]

## Use identical geographical support for all methods
P_ours_plot <- P_ours
P_mcar_plot <- P_mcar

P_ours_plot[!has_test, ] <- NA_real_
P_mcar_plot[!has_test, ] <- NA_real_


## ------------------------------------------------------------
## Overall MSE
##
## Average squared difference between observed test proportion
## and fitted probability over all retained locations and all
## crime categories.
## ------------------------------------------------------------

mse_ours <- mean(
  (P_obs[has_test, , drop = FALSE] -
     P_ours[has_test, , drop = FALSE])^2,
  na.rm = TRUE
)

mse_mcar <- mean(
  (P_obs[has_test, , drop = FALSE] -
     P_mcar[has_test, , drop = FALSE])^2,
  na.rm = TRUE
)

cat("Graph-Fused Multinomial MSE:", mse_ours, "\n")
cat("MCAR MSE:", mse_mcar, "\n")


## ------------------------------------------------------------
## Crime-category names
## ------------------------------------------------------------

category_names <- colnames(chicago$Y)

if (is.null(category_names)) {
  category_names <- paste0("Category ", seq_len(ncol(Y_test)))
}


## ------------------------------------------------------------
## Construct long spatial data
## ------------------------------------------------------------

map_list <- vector("list", length(category_names) * 3)

counter <- 1L

for (j in seq_along(category_names)) {
  
  ## Observed
  map_list[[counter]] <- st_sf(
    Category = category_names[j],
    Panel = "Observed",
    probability = P_obs[, j],
    geometry = st_geometry(chicago$grid_sf)
  )
  counter <- counter + 1L
  
  ## Graph-fused multinomial
  map_list[[counter]] <- st_sf(
    Category = category_names[j],
    Panel = "GFM",
    probability = P_ours_plot[, j],
    geometry = st_geometry(chicago$grid_sf)
  )
  counter <- counter + 1L
  
  ## MCAR
  map_list[[counter]] <- st_sf(
    Category = category_names[j],
    Panel = "MCAR",
    probability = P_mcar_plot[, j],
    geometry = st_geometry(chicago$grid_sf)
  )
  counter <- counter + 1L
}

map_long <- do.call(rbind, map_list)
mse_ours <- 0.01075
mse_mcar <- 0.01122

## ------------------------------------------------------------
## Add MSE to column headings
## ------------------------------------------------------------

ours_label <- paste0(
  "GFM\n(MSE = ",
  formatC(mse_ours, format = "f", digits = 5),
  ")"
)

mcar_label <- paste0(
  "MCAR\n(MSE = ",
  formatC(mse_mcar, format = "f", digits = 5),
  ")"
)

map_long$Panel <- factor(
  map_long$Panel,
  levels = c(
    "Observed",
    "GFM",
    "MCAR"
  ),
  labels = c(
    "Observed",
    ours_label,
    mcar_label
  )
)

## Preserve the desired category order
map_long$Category <- factor(
  map_long$Category,
  levels = category_names
)


## ------------------------------------------------------------
## Chicago boundary
## ------------------------------------------------------------

city_geom <- st_geometry(chicago$city_boundary) %>%
  st_transform(26916) %>%
  st_make_valid() %>%
  st_union()

if (any(st_geometry_type(city_geom) == "GEOMETRYCOLLECTION")) {
  city_geom <- st_union(
    st_collection_extract(city_geom, "POLYGON")
  )
}

city_plot <- st_sf(geometry = city_geom)


## ------------------------------------------------------------
## Clip maps to Chicago boundary
## ------------------------------------------------------------

map_long <- st_transform(map_long, st_crs(city_plot))
map_long <- st_make_valid(map_long)

st_agr(map_long) <- "constant"

map_long <- st_intersection(
  map_long,
  st_geometry(city_plot)
)

## Remove empty or non-polygon geometries
map_long <- map_long[
  !st_is_empty(map_long),
  ,
  drop = FALSE
]

if (any(st_geometry_type(map_long) == "GEOMETRYCOLLECTION")) {
  map_long <- st_collection_extract(
    map_long,
    "POLYGON",
    warn = FALSE
  )
}

map_long <- map_long[
  !st_is_empty(map_long) &
    st_dimension(map_long) == 2L,
  ,
  drop = FALSE
]

if (!nrow(map_long)) {
  stop("No grid polygons overlap the Chicago boundary.")
}


## ------------------------------------------------------------
## Plot
## ------------------------------------------------------------

p <- ggplot() +
  
  ## Grey city background
  geom_sf(
    data = city_plot,
    fill = "grey90",
    color = NA
  ) +
  
  ## Probabilities
  geom_sf(
    data = map_long,
    aes(fill = probability),
    color = NA
  ) +
  
  ## Chicago boundary
  geom_sf(
    data = city_plot,
    fill = NA,
    color = "grey30",
    linewidth = 0.25
  ) +
  
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1),
    na.value = "grey90",
    name = "Probability /\nproportion"
  ) +
  
  ## Rows = categories, columns = methods
  facet_grid(
    Category ~ Panel
  ) +
  
  coord_sf(
    crs = st_crs(city_plot),
    datum = NA
  ) +
  
  theme_void(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5
    ),
    
    ## Column headings
    strip.text.x = element_text(
      size = 11,
      lineheight = 1.2
    ),
    
    ## Crime names on rows
    strip.text.y = element_text(
      size = 11,
    ),
    
    strip.background = element_blank(),
    
    legend.position = "bottom",
    
    panel.spacing.x = unit(1, "cm"),
    panel.spacing.y = unit(0.15, "cm"),
    plot.margin = margin(
      5, 2, 5, 2
    )
  ) +
  
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(10, "cm"),
      barheight = unit(0.4, "cm")
    )
  ) +
  force_panelsizes(rows = unit(3.1, "cm"), cols = unit(3, "cm"))

print(p)

ggsave(
  "chicago_crime_plot.pdf",
  plot = p,
  width = 6.69,
  height = 9.0,     # leaves ~1.1in for caption + spacing
  units = "in"
)
# Optional export
# ggsave("chicago_observed_ours_mcar.png", p,
#        width = 16, height = 8, dpi = 300, bg = "white")

mse_t <- as.data.frame(matrix(c(0.01075022, 0.0112219, 10.5, 1618.85), nrow = 2, byrow = "TRUE"))
row.names(mse_t) <- c("MSE", "Time (secs)")
colnames(mse_t) <- c("Graph Fused multinomial", "MCAR")
mse_t

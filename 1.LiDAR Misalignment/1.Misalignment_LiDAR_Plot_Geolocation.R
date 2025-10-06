#install.packages("lidR")
#install.packages("terra")
#install.packages("ggplot2")
#install.packages("readr")
#install.packages("randomForest")
#install.packages('RANN')
#install.packages("sp")
#install.packages("raster")
#install.packages("dbscan")
#install.packages("ggnewscale")
library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(sf)
library(lidR)
library(RANN)
library(sp)
library(raster)
library(gridExtra)
library(dbscan) 
library(Morpho)

#-------------------------------------------------------------------------------
# LiDAR data
las_folder <- "2023.5.18_3DEP LiDAR"

# Forest inventory data
tree_file <- "Howland_Data/Howland_Tree_Data.csv"
plot_file <- "Howland_Data/Howland_Plot_Data.csv"

# Define the plot ID
plot_id <- 11

# Define the plot radius (m)
plot_radius <- 8.92

# Define the buffer distance (m)
buffer_dist <- 6

# Define the fraction of big trees
f_bigtree   <- 0.8

# Define the canopy size (m)
canopy_dist <- 3

# Define the stem distance (m)
tree_dist   <- 3

# Column names for GNSS in plot file
east  <- "Easting"
north <- "Northing"

#-------------------------------------------------------------------------------
# Read & merge LiDAR tiles (LAS/LAZ). Drop NULL tiles safely.
tile_paths <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE, ignore.case = TRUE)
las_list   <- lapply(tile_paths, lidR::readLAS)
las_list   <- Filter(Negate(is.null), las_list)
stopifnot(length(las_list) > 0)
las_m <- do.call(rbind, las_list)

#-------------------------------------------------------------------------------
# Tree and plot data
tree_data <- read.csv(tree_file, stringsAsFactors = FALSE) %>%
  mutate(Plot     = as.integer(Plot),
         Azimuth  = as.numeric(Azimuth),
         Distance = as.numeric(Distance),
         DBH      = as.numeric(DBH),
         Species  = as.factor(Species)) %>%
  filter(Plot == plot_id)

plot_data <- read.csv(plot_file, stringsAsFactors = FALSE) %>%
  filter(Plot == plot_id)

stopifnot(nrow(plot_data) == 1)

# Convert polar (Azimuth, Distance) to local Cartesian (x,y)
tree_data <- tree_data %>%
  mutate(Azimuth_rad = Azimuth * pi/180,
         x = Distance * sin(Azimuth_rad),
         y = Distance * cos(Azimuth_rad))

# Local plot circle for the map
circle_data <- data.frame(
  x = plot_radius * cos(seq(0, 2*pi, length.out = 200)),
  y = plot_radius * sin(seq(0, 2*pi, length.out = 200))
)

# Map trees in local coordinates
ggplot(tree_data, aes(x = x, y = y, size = DBH, fill = Species)) +
  geom_point(shape = 21, color = "black") +
  geom_path(data = circle_data, aes(x = x, y = y), color = "black", linetype = "dashed", inherit.aes = FALSE) +
  scale_size_continuous(name = "DBH (cm)") +
  scale_fill_manual(values = scales::hue_pal()(length(unique(tree_data$Species))),
                    name = "Species",
                    guide = guide_legend(override.aes = list(size = 6),
                                         keyheight = unit(2, 'lines'),
                                         keywidth  = unit(2, 'lines'))) +
  labs(title = paste("Tree Distribution in Plot", plot_id), x = "X (m)", y = "Y (m)") +
  coord_fixed() +
  theme_minimal() +
  theme(legend.position = "right",
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14))

#-------------------------------------------------------------------------------
# LiDAR processing around the plot center
e0 <- plot_data[[east]][1]
n0 <- plot_data[[north]][1]

# Build an sf point for the plot center with the same CRS as LAS (set EPSG as needed)
# If your LAS has a WKT, convert to EPSG or use st_crs(lidR::crs(las_m))
las_crs <- tryCatch(sf::st_crs(lidR::crs(las_m)), error = function(e) NA)
if (is.na(las_crs)) {
  # Fallback: set to UTM19N if that's your data (change if needed)
  las_crs <- sf::st_crs(26919)
}
center_sf <- st_sf(geometry = st_sfc(st_point(c(e0, n0)), crs = las_crs))

# Terrain & normalization
dtm_m <- rasterize_terrain(las_m, res = 1, algorithm = knnidw())
nlas_m <- normalize_height(las_m, dtm_m)

# Clip LAS to a buffer around the plot
r_buf  <- plot_radius + buffer_dist
buf_las <- filter_poi(nlas_m, (X - e0)^2 + (Y - n0)^2 <= r_buf^2)

# Quick visual
#lidR::plot(buf_las, bg = "white", size = 2)  

# Tree tops & CHM
ttops <- locate_trees(buf_las, lmf(ws = 2.5))
thr   <- c(0, 2, 5, 10, 15)
edg   <- c(0, 1.5)
chm_c <- rasterize_canopy(buf_las, res = 0.5, pitfree(thr, edg))

# Individual tree segmentation
trees <- segment_trees(buf_las, dalponte2016(chm_c, ttops, max_cr = canopy_dist,
                                             th_tree = 10, th_seed = 0.2, th_cr = 0.5))

# Crowns and centroids
crowns     <- delineate_crowns(trees, attribute = "treeID")
crowns_sf  <- st_as_sf(crowns)
centroids <- sf::st_centroid(crowns_sf, of_largest_polygon = TRUE)
coords <- sf::st_coordinates(centroids)
centroids$X <- coords[, 1]
centroids$Y <- coords[, 2]

# Extents for plotting
bb    <- st_bbox(crowns_sf)
x_min <- bb[['xmin']]; x_max <- bb[['xmax']]
y_min <- bb[['ymin']]; y_max <- bb[['ymax']]

ggplot() +
  geom_sf(data = crowns_sf, fill = NA, color = "blue", linewidth = 0.5) +
  geom_sf(data = centroids, shape = 21, fill = "green", color = "black", size = 3, alpha = 0.8) +
  coord_sf(xlim = c(x_min, x_max), ylim = c(y_min, y_max), datum = NA, expand = FALSE) +
  theme_minimal() +
  labs(title = paste("Plot", plot_id, "– canopy crowns and stem locations"),
       x = "Easting (m)", y = "Northing (m)")

#-------------------------------------------------------------------------------
# Inventory big trees in absolute coordinates
db_thresh <- stats::quantile(tree_data$DBH, 1 - f_bigtree, na.rm = TRUE)

big_trees <- tree_data %>%
  filter(DBH >= db_thresh) %>%
  mutate(Xinv = e0 + x, Yinv = n0 + y)

# Further filter big trees so any two trees are at least 3 m apart
inv_abs <- big_trees %>%
  transmute(X_abs = Xinv, Y_abs = Yinv, DBH = DBH, Species = Species)

set.seed(42)                     
ord <- sample(nrow(inv_abs))
                   
keep <- integer(0)
for (i in ord) {
  if (length(keep) == 0) {
    keep <- i
  } else {
    d <- sqrt((inv_abs$X_abs[i] - inv_abs$X_abs[keep])^2 +
                (inv_abs$Y_abs[i] - inv_abs$Y_abs[keep])^2)
    if (all(d >= tree_dist)) keep <- c(keep, i)
  }
}

filter_trees <- inv_abs[keep, ]
cat("Number of big trees after 3 m thinning:", nrow(filter_trees), "\n")

# Plot circle as sf polygon
circle_orig <- st_buffer(center_sf$geometry, dist = plot_radius)
circle_sf   <- st_sf(type = "original", geometry = circle_orig)

ggplot() +
  geom_sf(data = crowns_sf, fill = NA, color = "green") +
  geom_sf(data = circle_sf, fill = NA, color = "black", linetype = "solid", linewidth = 1) +
  geom_sf(data = centroids, shape = 21, fill = "skyblue", color = "blue", size = 3, alpha = 0.8) +
  geom_point(data = filter_trees, aes(x = X_abs, y = Y_abs), shape = 21, fill = NA, color = "black", size = 3, stroke = 0.8) +
  guides(size = "none", fill = "none") +
  coord_sf(xlim = c(x_min, x_max), ylim = c(y_min, y_max), datum = NA, expand = FALSE) +
  theme_minimal() +
  labs(title = paste("Plot", plot_id, "– Crowns, Centroids, Inventory Trees & Plot Boundary"),
       x = "Easting (m)", y = "Northing (m)")

#-------------------------------------------------------------------------------
# Correct plot location 
# Build inventory points in absolute coordinates
inv_spatial <- filter_trees %>%transmute(X_abs = X_abs, Y_abs = Y_abs, DBH = DBH)

# Build pure numeric matrices
pred_mat  <- as.matrix(st_coordinates(centroids))
field_mat <- as.matrix(inv_spatial[, c("X_abs", "Y_abs"), drop = FALSE])

# Find the nearest LiDAR centroid for each field point
nn  <- RANN::nn2(data = pred_mat, query = field_mat, k = 1)
idx <- nn$nn.idx[, 1]
d   <- nn$nn.dists[, 1]

# Keep only those pairs within some distance threshold
thresh <- 5   # e.g., 5 m maximum
keep_pairs <- which(d < thresh)

# If not enough matches, fall back to using the original field points without transform
if (length(keep_pairs) < 3) {
  warning("Fewer than 3 matched pairs within threshold; skipping transform and using original inventory points.")
  adj_field_mat <- field_mat
  new_center    <- c(e0, n0)
} else {
  field_match <- field_mat[keep_pairs, , drop = FALSE]
  pred_match  <- pred_mat[idx[keep_pairs], , drop = FALSE]
  
  # Do the Procrustes fit on those matched subsets
  trafo_raw <- Morpho::computeTransform(x = field_match, y = pred_match, type = "rigid", reflection = FALSE)
  R_raw <- trafo_raw[1:2, 1:2]
  t_raw <- trafo_raw[1:2, 3]
  
  # Extract raw rotation and clamp
  angle_raw <- atan2(R_raw[2, 1], R_raw[1, 1]) * 180 / pi
  angle_lim <- max(min(angle_raw, 20), -20)
  
  # Rebuild a pure-rotation matrix from the clamped angle
  theta <- angle_lim * pi / 180
  R_lim <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), ncol = 2, byrow = TRUE)
  
  # Recompute translation so the matched centroids still align
  centF <- colMeans(field_match)
  centP <- colMeans(pred_match)
  t_lim <- centP - (R_lim %*% centF)
  
  # Build the limited 3×3 transform and apply to the full stem cloud
  trafo_lim <- diag(3)
  trafo_lim[1:2, 1:2] <- R_lim
  trafo_lim[1:2, 3]   <- t_lim
  
  adj_field_mat <- Morpho::applyTransform(field_mat, trafo_lim)
  
  # Compute new plot centre
  new_center <- as.numeric(R_lim %*% c(e0, n0) + t_lim)
}

# Compute circles
orig_pt    <- st_sfc(st_point(c(e0, n0)), crs = st_crs(crowns_sf))
calib_pt   <- st_sfc(st_point(new_center), crs = st_crs(crowns_sf))
circle_orig  <- st_buffer(orig_pt,  dist = plot_radius)
circle_calib <- st_buffer(calib_pt, dist = plot_radius)

# Combine into one sf
circles_sf <- rbind(
  st_sf(type = "original",   geometry = circle_orig),
  st_sf(type = "calibrated", geometry = circle_calib)
)

# Compute dx and dy
dx <- new_center[1] - e0
dy <- new_center[2] - n0

# Determine movement direction
dir_x <- if (dx >= 0) "East" else "West"
dir_y <- if (dy >= 0) "North" else "South"

# Determine rotation direction (if transform happened, angle_lim exists; otherwise set to 0)
if (exists("angle_lim")) {
  angle <- angle_lim
} else {
  angle <- 0
}
rot_dir <- if (angle >= 0) "counterclockwise" else "clockwise"

# Print summary
cat("Original center:     ", round(e0, 1), ",", round(n0, 1), "\n",
    "Calibrated center:   ", round(new_center[1], 1), ",", round(new_center[2], 1), "\n",
    "Moved:               ", round(abs(dx), 1), "m", dir_x, ",", round(abs(dy), 1), "m", dir_y, "\n",
    "Rotation:            ", round(abs(angle), 1), "degrees", rot_dir, "\n")

# Plot adjusted stems
adj_df <- data.frame(X = adj_field_mat[, 1], Y = adj_field_mat[, 2], DBH = inv_spatial$DBH)

# Compute bounding box to include crowns and calibrated circle
bb1  <- st_bbox(crowns_sf)
bb2  <- st_bbox(circle_calib)
xlim <- range(c(bb1["xmin"], bb1["xmax"], bb2["xmin"], bb2["xmax"]))
ylim <- range(c(bb1["ymin"], bb1["ymax"], bb2["ymin"], bb2["ymax"]))

ggplot() +
  geom_sf(data = crowns_sf, fill  = NA, color = "green", size = 0.5) +
  geom_sf(data = filter(circles_sf, type == "calibrated"),
          fill = NA, color = "black", linetype = "dashed", size = 1) +
  geom_sf(data = centroids, shape = 21, fill = "skyblue", color = "blue", size = 3, alpha = 0.8) +
  geom_point(data = adj_df, aes(x = X, y = Y), shape = 21, fill = NA, color = "black", size = 3, stroke = 0.8) +
  guides(size = "none", fill = "none") +
  coord_sf(xlim = xlim, ylim = ylim, datum = NA, expand = FALSE) +
  theme_minimal() +
  labs(title = paste("Plot", plot_id, "— Calibrated Boundary & Trees"),
       x = "Easting (m)", y = "Northing (m)")

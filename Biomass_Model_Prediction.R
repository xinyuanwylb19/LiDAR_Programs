#install.packages("lidR")
#install.packages("terra")
#install.packages("ggplot2")
#install.packages("readr")
#install.packages("randomForest")
#install.packages("viridis")

library(lidR)
library(sf)
library(terra)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(randomForest)
library(viridis)
library(purrr)
#-------------------------------------------------------------------------------
# Read forest inventory data for model calibration
inv_data <- read_csv("Howland_Data/Howland_Plot_Calibration.csv")

# Read the validation data
val_data <- read_csv("Howland_Data/Howland_Plot_Validation.csv")

# Define the location data
east  <- "Easting_3"
north <- "Northing_3"

# Radius for the inventory plot (meters) 8.92m or 12.61m
plot_radius <- 8.92

# Read and merge LiDAR tiles to a single LAS object
las_folder <- "2023.5.18_3DEP LiDAR"
tile_paths <- list.files(las_folder, pattern = "\\.las$", full.names = TRUE)
las_list   <- lapply(tile_paths, readLAS)
las_m      <- do.call(rbind, las_list)

# t-test (times) and 
n_test    <- 100
n_samples <- 15

# Save the model file
savemodel <- "Howland_Results/rf_model_3.rds"
savefile   <- "Howland_Results/Howland_Biomass_2023_3.tif"
val_paired <- "Howland_Results/Validation_Paired_3.csv"
val_result <- "Howland_Results/Validation_Results_3.csv"
#-------------------------------------------------------------------------------
# Calibrate the random forest model
# Build DTM and normalize heights
dtm_m   <- rasterize_terrain(las_m, res = 1, algorithm = knnidw())
nlas_m  <- normalize_height(las_m, dtm_m)

# Convert to sf points and reproject into LiDAR CRS, adjust to match LAS EPSG
inv_points <- st_as_sf(inv_data, coords = c(east,north), crs = 32619, remove = FALSE) %>%
  st_transform(st_crs(nlas_m))

# Unified metrics
compute_metrics <- function(z, th = 2) {
  z <- z[is.finite(z)]
  if (length(z) == 0) return(data.frame(
    mean_height=NA_real_, max_height=NA_real_, p50=NA_real_,
    p75=NA_real_, p90=NA_real_, canopy_cover=NA_real_
  ))
  cf <- sum(z > th) / length(z)
  data.frame(
    mean_height  = mean(z),
    max_height   = max(z),
    p50          = quantile(z, 0.50),
    p75          = quantile(z, 0.75),
    p90          = quantile(z, 0.90),
    canopy_cover = cf
  )
}

# Plot-level metrics wrapper (LAS -> uses Z and th = 2 m)
metrics <- function(las) compute_metrics(las$Z, th = 2)

# Loop over each plot
metrics_list <- lapply(seq_len(nrow(inv_points)), function(i) {
  pt   <- inv_points[i, ]
  coords <- st_coordinates(pt)
  x0 <- coords[1]; y0 <- coords[2]
  las_clip <- filter_poi(nlas_m, (X - x0)^2 + (Y - y0)^2 <= plot_radius^2)
  if (is.null(las_clip) || las_clip@header@PHB$`Number of point records` == 0) {
    return(NULL)
  }
  df <- metrics(las_clip)
  df$Plot <- pt$Plot
  df
})
lidar_metrics <- bind_rows(metrics_list)
feats <- setdiff(names(lidar_metrics), "Plot")

# Merge LiDAR metrics to the plot biomass data, fit a Random Forest model
bio_rf_data <- inv_data %>%
  dplyr::select(Plot, Biomass) %>%
  dplyr::inner_join(lidar_metrics, by = "Plot")

set.seed(42)
frm <- reformulate(termlabels = feats, response = "Biomass")
rf_model <- randomForest(frm, data = bio_rf_data, importance = TRUE, ntree = 500)
print(rf_model)

# Save model
saveRDS(rf_model, savemodel)
#-------------------------------------------------------------------------------
# Apply the random forest model and generate the biomass map (15×15 m grid)
metrics_grid <- grid_metrics(nlas_m, ~compute_metrics(Z, th = 3), res = 15)
metrics_df   <- as.data.frame(metrics_grid, xy = TRUE) %>%
  drop_na() %>%
  mutate(across(where(is.list), unlist))
metrics_df$biomass <- predict(rf_model, newdata = metrics_df[, feats, drop = FALSE])

# Rasterize and save the biomass map
biomass_raster <- rast(
  metrics_df[, c("x","y","biomass")],
  type = "xyz", crs = st_crs(nlas_m)$wkt
)
plot(biomass_raster, main = "Forest Aboveground Biomass (Mg/ha)")
writeRaster(biomass_raster, savefile, overwrite = TRUE)
#-------------------------------------------------------------------------------
# Validate the biomass map using the t-test
bio_r <- rast(savefile)
names(bio_r) <- 'AGB'

east  <- "Easting"
north <- "Northing"

# Points from validation data -> sf in EPSG:32619, then to raster CRS
val_pts <- st_as_sf(val_data, coords = c(east, north), crs = 32619, remove = FALSE) |>
  st_transform(st_crs(bio_r))

# Extract predicted biomass at plot locations
xy  <- st_coordinates(val_pts)
pred <- terra::extract(bio_r, val_pts, ID = FALSE)$AGB

# Build paired table and save
paired <- val_data |>
  mutate(Predicted = pred) |>
  select(Plot, Biomass, Predicted, Easting, Northing)
paired_cc <- paired |>
  filter(is.finite(Biomass), is.finite(Predicted))
write_csv(paired_cc, val_paired)

# Repeated t-tests on random samples (paired)
res <- map_dfr(1:n_test, function(i) {
  samp <- dplyr::slice_sample(paired_cc, n = n_samples, replace = FALSE)
  tt   <- t.test(samp$Predicted, samp$Biomass, paired = TRUE)
  tibble(
    iter       = i,
    n          = nrow(samp),
    t          = unname(tt$statistic),
    df         = unname(tt$parameter),
    p_value    = tt$p.value,
    mean_obs   = mean(samp$Biomass),
    mean_pred  = mean(samp$Predicted),
    mean_diff  = mean(samp$Predicted - samp$Biomass),
    sd_diff    = sd(samp$Predicted - samp$Biomass)
  )
})
#print(res)
write_csv(res, val_result)

#-------------------------------------------------------------------------------
# Compare biomass maps
# Define the standard map and boundary
standard_file <- "Howland_Results/Howland_Biomass_2023_i.tif"
boundary_file <- "Howland_Data/Howland_Boundary/Howland Boundary.shp"

# Load all biomass maps from the results folder
tif_paths <- list.files("Howland_Results", pattern = "Howland_Biomass_2023.*\\.tif$",
                        full.names = TRUE)
tif_names <- tools::file_path_sans_ext(basename(tif_paths))

# Load standard raster and boundary, reproject boundary to raster CRS
std_r    <- rast(standard_file)
boundary <- st_read(boundary_file, quiet = TRUE) %>%
  st_transform(crs(std_r))

# Crop and mask each raster to the boundary region
rast_list <- lapply(tif_paths, function(f) {
  r <- rast(f)
  mask(crop(r, boundary), vect(boundary))
})
names(rast_list) <- tif_names

# Crop and mask the standard map
std_masked <- mask(crop(std_r, boundary), vect(boundary))

#-------------------------------------------------------------------------------
# Average AGB per map and percentage difference from the standard map
mean_standard <- global(std_masked, "mean", na.rm = TRUE)$mean

agb_summary <- data.frame(
  Map      = tif_names,
  Mean_AGB = sapply(rast_list, function(r) global(r, "mean", na.rm = TRUE)$mean),
  Std_AGB  = sapply(rast_list, function(r) global(r, "sd",   na.rm = TRUE)$sd)
) %>%
  mutate(Pct_Diff = (Mean_AGB - mean_standard) / mean_standard * 100)

print(agb_summary)

#-------------------------------------------------------------------------------
# Per-cell percentage difference from the standard map
std_name <- tools::file_path_sans_ext(basename(standard_file))

cell_diff_summary <- lapply(tif_names, function(nm) {
  # Skip the standard map itself
  if (nm == std_name) return(NULL)

  # Resample to the standard grid if resolutions differ
  r_aligned <- resample(rast_list[[nm]], std_masked, method = "bilinear")

  # Per-cell absolute percentage difference: |map - standard| / standard * 100
  pct_diff <- abs(r_aligned - std_masked) / std_masked * 100
  vals     <- values(pct_diff, na.rm = TRUE)

  data.frame(
    Map      = nm,
    Mean_Pct = round(mean(vals), 2),
    Std_Pct  = round(sd(vals),   2)
  )
})
cell_diff_result <- bind_rows(cell_diff_summary)
print(cell_diff_result)

#-------------------------------------------------------------------------------
# Plot observed and predicted AGB
# Read validation data with all predicted columns
val_plot <- read_csv("Howland_Data/Howland_Plot_Validation.csv")

# Reshape to long format for faceted plot
val_long <- val_plot %>%
  pivot_longer(
    cols      = starts_with("Predicted AGB"),
    names_to  = "Model",
    values_to = "Predicted"
  ) %>%
  rename(Observed = `Observed AGB`)

# Define panel order and labels
model_levels <- c("Predicted AGB (~0m)", "Predicted AGB (1m)",  "Predicted AGB (3m)",
                  "Predicted AGB (5m)",  "Predicted AGB (10m)", "Predicted AGB (15m)")
panel_labels <- setNames(
  c("~m", "1 m", "3 m", "5 m", "10 m", "15 m"),
  model_levels
)

val_long$Model <- factor(val_long$Model, levels = model_levels)

# Common axis range so the 1:1 line is meaningful across all panels
agb_range <- range(c(val_long$Observed, val_long$Predicted), na.rm = TRUE)

# Six-panel scatter plot (2 rows x 3 columns)
ggplot(val_long, aes(x = Observed, y = Predicted)) +
  geom_point(shape = 16, color = "black", size = 3) +
  geom_abline(intercept = 0, slope = 1, color = "blue", linewidth = 0.5) +
  facet_wrap(~ Model, nrow = 2, ncol = 3, labeller = as_labeller(panel_labels)) +
  scale_x_continuous(limits = agb_range) +
  scale_y_continuous(limits = agb_range) +
  coord_fixed() +
  labs(title = NULL, x = "Observed AGB (MgC/ha)", y = "Predicted AGB (MgC/ha)") +
  theme_bw() +
  theme(
    text         = element_text(family = "Times New Roman", size = 18, color = "black"),
    strip.text   = element_text(family = "Times New Roman", size = 18, color = "black"),
    axis.text    = element_text(family = "Times New Roman", size = 18, color = "black"),
    axis.title   = element_text(family = "Times New Roman", size = 18, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )    
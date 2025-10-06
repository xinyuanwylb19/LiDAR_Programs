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
library(randomForest)
library(viridis)
library(purrr)
#-------------------------------------------------------------------------------
# Read forest inventory data for model calibration
inv_data <- read_csv("Howland_Data/Howland_Plot_Calibration.csv")

# Read the validation data
val_data <- read_csv("Howland_Data/Howland_Plot_Validation.csv")

# Define the location data
east  <- "Easting_15"
north <- "Northing_15"

# Radius for the inventory plot (meters) 8.92m or 12.61m
plot_radius <- 8.92

# Read and merge LiDAR tiles to a single LAS object
las_folder <- "2023.5.18_3DEP LiDAR"
tile_paths <- list.files(las_folder, pattern = "\\.las$", full.names = TRUE)
las_list   <- lapply(tile_paths, readLAS)
las_m      <- do.call(rbind, las_list)

# t-test (times) and 
n_test    <- 100
n_samples <- 13

# Save the model file
savemodel <- "Howland_Results/rf_model_15m.rds"
savefile   <- "Howland_Results/Howland_Biomass_2023_15m.tif"
val_paired <- "Howland_Results/Validation_Paired_15m.csv"
val_result <- "Howland_Results/Validation_Results_15m.csv"
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
# Apply the random forest model and generate the biomass map
# Apply over a 10×10 m grid and predict biomass 
metrics_grid <- grid_metrics(nlas_m, ~compute_metrics(Z, th = 3), res = 10)
metrics_df   <- as.data.frame(metrics_grid, xy = TRUE) %>%
  drop_na() %>%
  mutate(across(where(is.list), unlist))
metrics_df$biomass <- predict(rf_model, newdata = metrics_df[, feats, drop = FALSE])

# Rasterize and save the biomass map
biomass_raster <- rast(
  metrics_df[, c("x","y","biomass")],
  type = "xyz", crs = st_crs(nlas_m)$wkt
)
plot(biomass_raster, main = "Forest Biomass (kgC/m²)")
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
# Compare biomass maps and summarize within Howland Forest
res_folder <- "Howland_Results"
std_path   <- file.path(res_folder, "Howland_Biomass_2023_im.tif")

# Candidate maps (excluding the standard)
other_files <- c("Howland_Biomass_2023_no.tif",
                 "Howland_Biomass_2023_1m.tif",
                 "Howland_Biomass_2023_3m.tif",
                 "Howland_Biomass_2023_5m.tif",
                 "Howland_Biomass_2023_10m.tif",
                 "Howland_Biomass_2023_15m.tif")

# Read standard raster
std_rst <- raster::raster(std_path)

# Read boundary
bnd <- sf::st_read("Howland_Data/Howland_Boundary/Howland Boundary.shp", quiet = TRUE)

for (fname in other_files) {
  # Read target raster
  tgt_path <- file.path(res_folder, fname)
  tgt_rst  <- raster::raster(tgt_path)
  
  # Resample standard to target grid (continuous variable => bilinear)
  std_on_tgt <- raster::resample(std_rst, tgt_rst, method = "bilinear")
  
  # Difference rasters
  diff_rst     <- tgt_rst - std_on_tgt
  abs_diff_rst <- abs(diff_rst)
  
  # Per-pixel % difference (abs)
  pct_diff_rst <- raster::overlay(diff_rst, std_on_tgt, fun = function(d, s) {
    out <- rep(NA_real_, length(d))
    ok  <- is.finite(s) & s != 0 & is.finite(d)
    out[ok] <- 100 * abs(d[ok]) / s[ok]
    out
  })
  
  # Save rasters
  out_name_diff <- sub("\\.tif$", "_minus_non.tif", fname)
  out_name_abs  <- sub("\\.tif$", "_abs_diff.tif", fname)
  out_name_pct  <- sub("\\.tif$", "_pct_diff.tif", fname)
  
  raster::writeRaster(diff_rst, filename = file.path(res_folder, out_name_diff),
                      overwrite = TRUE, options = c("COMPRESS=LZW"), datatype = "FLT4S")
  raster::writeRaster(abs_diff_rst, filename = file.path(res_folder, out_name_abs),
                      overwrite = TRUE, options = c("COMPRESS=LZW"), datatype = "FLT4S")
  raster::writeRaster(pct_diff_rst, filename = file.path(res_folder, out_name_pct),
                      overwrite = TRUE, options = c("COMPRESS=LZW"), datatype = "FLT4S")
  
  # Transform boundary
  bnd_tgt <- sf::st_transform(bnd, sf::st_crs(raster::crs(diff_rst)))
  bnd_sp  <- as(bnd_tgt, "Spatial")
  
  # Mask/crop within Howland
  diff_crop     <- raster::mask(raster::crop(diff_rst,     bnd_sp), bnd_sp)
  abs_diff_crop <- raster::mask(raster::crop(abs_diff_rst, bnd_sp), bnd_sp)
  pct_diff_crop <- raster::mask(raster::crop(pct_diff_rst, bnd_sp), bnd_sp)
  
  tgt_crop  <- raster::mask(raster::crop(tgt_rst,  bnd_sp), bnd_sp)
  std_crop  <- raster::mask(raster::crop(std_on_tgt, bnd_sp), bnd_sp)
  
  # Values
  vals_diff_roi     <- raster::getValues(diff_crop)
  vals_abs_diff_roi <- raster::getValues(abs_diff_crop)
  vals_pct_diff_roi <- raster::getValues(pct_diff_crop)
  vals_tgt_roi      <- raster::getValues(tgt_crop)
  vals_std_roi      <- raster::getValues(std_crop)
  
  # Stats
  mean_diff_roi     <- mean(vals_diff_roi,     na.rm = TRUE)
  sd_diff_roi       <- stats::sd(vals_diff_roi, na.rm = TRUE)
  
  mean_abs_diff_roi <- mean(vals_abs_diff_roi, na.rm = TRUE)
  sd_abs_diff_roi   <- stats::sd(vals_abs_diff_roi, na.rm = TRUE)
  
  mean_pct_diff_roi <- mean(vals_pct_diff_roi, na.rm = TRUE)
  sd_pct_diff_roi   <- stats::sd(vals_pct_diff_roi, na.rm = TRUE)
  
  mean_tgt_roi <- mean(vals_tgt_roi, na.rm = TRUE)
  sd_tgt_roi   <- stats::sd(vals_tgt_roi, na.rm = TRUE)
  
  mean_std_roi <- mean(vals_std_roi, na.rm = TRUE)
  sd_std_roi   <- stats::sd(vals_std_roi, na.rm = TRUE)
  
  pct_diff_roi <- if (is.finite(mean_std_roi) && mean_std_roi != 0) {
    100 * (mean_tgt_roi - mean_std_roi) / mean_std_roi
  } else {
    NA_real_
  }
  
  # Print summary
  cat(
    "\nFile: ", fname, "\n",
    "Mean diff (tgt - std):     ", round(mean_diff_roi,     4), "\n",
    "Std  diff:                 ", round(sd_diff_roi,       4), "\n",
    "Mean abs diff:             ", round(mean_abs_diff_roi, 4), "\n",
    "Std  abs diff:             ", round(sd_abs_diff_roi,   4), "\n",
    "Mean % diff (per-pixel):   ", round(mean_pct_diff_roi, 2), "%\n",
    "Std  % diff (per-pixel):   ", round(sd_pct_diff_roi,   2), "%\n",
    "Mean target:               ", round(mean_tgt_roi,      4), "\n",
    "Std  target:               ", round(sd_tgt_roi,        4), "\n",
    "Mean standard:             ", round(mean_std_roi,      4), "\n",
    "Std  standard:             ", round(sd_std_roi,        4), "\n",
    "Overall % diff (mean tgt vs mean std): ",
    round(pct_diff_roi, 2), "%\n",
    sep = ""
  )
}

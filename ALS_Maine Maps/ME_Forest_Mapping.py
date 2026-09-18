# Created on Sep 15, 2026, Author: Xinyuan Wei

import glob
import json
import os
import time
import numpy as np
import pandas as pd
import rasterio
import requests
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import ListedColormap
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import from_bounds
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from rasterio.warp import reproject
from rasterio.warp import transform
from rasterio.windows import Window

# -----------------------------------------------------------------------------
# Parameters.
# -----------------------------------------------------------------------------
# Sections to run.
data_prepare = False      # True or False
reconstruction = False
delineation = False
validation = False
figure = True

# Folders.
Root = r"C:\Users\xinyuan.wei\Desktop\LiDAR-5.Maine Biomass"
Work_dir = os.path.join(Root, "1.Data", "2.Model and Programs")

Lcms_dir = os.path.join(Work_dir, "LCMS")
Cover_dir = os.path.join(Work_dir, "Cover")
Layer_dir = os.path.join(Work_dir, "Layers")
Stats_dir = os.path.join(Work_dir, "Stats")
Flight_agb_dir = os.path.join(Root, "1.Data", "3.3DEP_AGB")
Flight_vol_dir = os.path.join(Root, "1.Data", "4.3DEP_Volume")
Agb_dir = os.path.join(Root, "1.Data", "5.Annual AGB")
Volume_dir = os.path.join(Root, "1.Data", "6.Annual Volume")
Figure_dir = os.path.join(Root, "2.Final Paper", "1.Figures")

Loss_file = "Loss Calibration.xlsx"
Mask_table = "Forest Mask.xlsx"

State_file = "State Comparison.xlsx"
State_figure = "Figure 5_Validation.png"
Mask_path = os.path.join(Layer_dir, "Forest Mask.tif")

# Years to map, and the years with lidar data.
Years = list(range(2012, 2025))
Anchor_years = [2012, 2013, 2015, 2016, 2017, 2018, 2020, 2021, 2022, 2023, 2024]

# Pairs of years where the same ground was flown twice, newest pair first.
Pairs = [(2016, 2024), (2015, 2023), (2013, 2022), (2012, 2021)]

# Output grid, the same one the biomass and volume maps use.
Cell = 15.0
Maine_bounds = (335985.0, 4767990.0, 660015.0, 5257995.0)
Maine_crs = "EPSG:6348"

# Ground covered by one 15 m cell, in hectares and in square kilometres.
Cell_ha = 0.0225
Cell_km2 = 0.000225

# The LCMS data.
Change_zip = "https://data.fs.usda.gov/geodata/LCMS/LCMS_CONUS_v%s_Change_Annual_%d.zip"
Change_tif = "LCMS_CONUS_v%s_Change_%d.tif"
Cover_zip = "https://data.fs.usda.gov/geodata/LCMS/LCMS_CONUS_v%s_Land_Cover_Annual_%d.zip"
Cover_tif = "LCMS_CONUS_v%s_Land_Cover_%d.tif"
Version = "2025-11"

# Change classes, as the Forest Service names them.
Change_names = {1: "slow loss", 2: "fast loss drought", 3: "fast loss mining", 4: "desiccation",
                5: "inundation", 6: "prescribed fire", 7: "wildfire",
                8: "mechanical transformation", 9: "tree removal", 10: "defoliation",
                11: "southern pine beetle", 12: "insect disease or drought", 13: "other loss",
                14: "vegetation growth", 15: "stable", 16: "non processing mask"}

# The classes that mean the stand lost something.
Loss_classes = [1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13]

# The classes that take wood off the ground.
Wood_loss_classes = [1, 2, 3, 6, 7, 8, 9, 13]

# Land cover classes.
Cover_names = {1: "trees", 2: "tall shrubs and trees", 3: "shrubs and trees",
               4: "grass and trees", 5: "barren and trees", 6: "tall shrubs", 7: "shrubs",
               8: "grass and shrubs", 9: "barren and shrubs", 10: "grass forb herb",
               11: "barren and grass", 12: "barren or impervious", 13: "snow or ice",
               14: "water"}

# Ground the survey did not process, written as no data.
Skip_class = 15

# Classes with trees.
Forest_classes = [1, 2, 3, 4, 5]

# A cell the forest survey never calls forest is still kept.
Min_carbon = 10.0

Acre_to_ha = 0.404686

# Smallest standing biomass worth following, in Mg per ha.
Min_start = 20.0

# Bin width and range for biomass change, in Mg per ha.
Bin_width = 5.0
Bin_limit = 150.0

# Loss sizes reported separately, in Mg per ha.
Loss_steps = [10.0, 20.0, 40.0]

# Standing biomass classes used for the growth table, in Mg per ha.
Stock_edges = [20, 40, 60, 80, 100, 130, 1000]

# Loss that counts as a cut the survey missed, in Mg per ha.
Quiet_loss = 10.0

# Yearly net change as a fraction of standing biomass, by biomass class in Mg per ha.
Growth_edges = [0, 20, 40, 60, 80, 100, 130, 10000]
Growth_rates = [0.029, 0.029, 0.015, 0.012, 0.009, 0.005, 0.004]

# Share of standing stock taken out by one flagged year, by change class.
Loss_fraction = {1: 0.25, 2: 0.25, 3: 0.55, 6: 0.43, 7: 0.25, 8: 0.55, 9: 0.43, 13: 0.25}

# Insect damage, a smaller loss applied a few years later.
Slow_classes = [10, 11, 12]
Slow_fraction = 0.09
Slow_delay = 3

# Cells below this are treated as open ground, in Mg per ha.
Min_stock = 1.0

# Cap on the reconstructed past, in Mg per ha.
Max_stock = 250.0

# Volume to biomass ratio, used only for that cap.
Volume_ratio = 3.4

# The survey estimate comes from EVALIDator, the estimator the Forest Service publishes.
Evalidator_url = "https://apps.fs.usda.gov/fiadb-api/fullreport"

# Evaluation group, the state code followed by the year its five annual panels end on.
Eval_group = "23%d"

# Estimate attributes, forest area, aboveground carbon, and total stem wood and bark volume.
# Carbon and volume both start at one inch, so the same trees stand behind the two numbers.
Area_attribute = 2
Carbon_attribute = 53000
Volume_attribute = 11069

# Row grouping for the county estimates, and no grouping for a state total.
County_rows = "County code and name"
State_rows = "None"

# EVALIDator publishes acres, short tons and cubic feet.
Ston_to_Mg = 0.90718474
Cf_to_m3 = 0.028316846592

# The year the county totals are compared on.
County_year = 2023

# Map products, with folder, file pattern, figure name, colour bar label, range and tick step.
Products = [("AGB", Agb_dir, "%s_AGB_Maine.tif", "Figure 3_AGB_maps.png",
             "Aboveground biomass (Mg C ha$^{-1}$)", (0, 150), 30),
            ("Volume", Volume_dir, "%s_Vol_Maine.tif", "Figure 4_Volume_maps.png",
             "Volume (m$^{3}$ ha$^{-1}$)", (0, 600), 200)]

# A coarse map cell is drawn when at least this share of its 15 m cells are forest.
Mask_share = 0.5

# Map window in degrees and the pixel grid the maps are drawn on.
Map_extent = (-71.2, -66.8, 42.9, 47.55)
Map_size = (1600, 2400)

# Width of the aerial photo in pixels, its height follows the map window.
Imagery_width = 2000

# County boundaries and the aerial photo, fetched once into Layer_dir. Delete them to refetch.
County_url = "https://tigerweb.geo.census.gov/arcgis/rest/services/Generalized_ACS2023/State_County/MapServer/11/query"
Imagery_url = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/export"
County_file = "Maine Counties.geojson"
Imagery_file = "Maine Imagery.jpg"

# The annual maps on the map grid, cached in Layer_dir.
Stack_file = "Annual Map Grids.npz"

# 15 m cells averaged into one coarse cell before drawing.
Coarsen = 10

# Colour map, cold to warm. Only the lower part is used, so the top is a dark yellow.
Colormap = "Spectral_r"
Colormap_top = 0.62

# Panels per figure as rows and columns. The last slot holds the colour bar.
Layout = (3, 5)

# Longitude and latitude ticks, in degrees.
Lon_ticks = [-70, -68]
Lat_ticks = [44, 46]

# Panel width in inches. The height follows the map window.
Panel_width = 1.35

# The coverage figure, one panel per flight year and one for all of them.
Coverage_figure = "Figure 1_ALS_coverage.png"
Coverage_stack = "Flight Coverage Grids.npz"
Coverage_layout = (3, 4)
Coverage_colour = "#87cefa"
Coverage_label = "All flights"

Panel_labels = ["(a)", "(b)", "(c)", "(d)"]

Font_size = 10

# Text size in the validation figure, and the colours of its two series.
Validation_font = 14
Map_colour = "green"
Fia_colour = "blue"

Figure_dpi = 600

# Values are stored as whole numbers after multiplying by this.
Scale = 10.0

# No data values, one per map type.
Flight_nodata = -9999.0
Annual_nodata = -32768
Lcms_nodata = 255

# Rows held at a time, one setting per section.
Download_strip = 4096
Loss_strip = 512
Annual_strip = 64
Mask_strip = 128

State_strip = 512

# -----------------------------------------------------------------------------
# The growth rule, shared by Sections 2 and 4.
# -----------------------------------------------------------------------------
# Yearly growth as a fraction of the biomass standing at the start of the year.
def growth_fraction(stock):
    band = np.digitize(stock, Growth_edges) - 1
    return np.take(np.array(Growth_rates, np.float32), np.clip(band, 0, len(Growth_rates) - 1))

# -----------------------------------------------------------------------------
# Section 1. Data preprocessing.
# -----------------------------------------------------------------------------
if data_prepare:

    # Output folders and the grid both maps are written on.
    os.makedirs(Lcms_dir, exist_ok=True)
    os.makedirs(Cover_dir, exist_ok=True)
    mx0, my0, mx1, my1 = Maine_bounds
    width = int(round((mx1 - mx0) / Cell))
    height = int(round((my1 - my0) / Cell))
    grid = from_origin(mx0, my1, Cell, Cell)
    print("output grid %d x %d at %g m, %s" % (width, height, Cell, Maine_crs))

    # One change map per year, on the Maine grid.
    profile = {"driver": "GTiff", "height": height, "width": width, "count": 1, "dtype": "uint8",
               "crs": Maine_crs, "nodata": Lcms_nodata, "transform": grid,
               "compress": "deflate", "tiled": True, "blockxsize": 256, "blockysize": 256,
               "BIGTIFF": "IF_SAFER"}

    for year in Years:
        out_path = os.path.join(Lcms_dir, "%d_Change_Maine.tif" % year)
        if os.path.exists(out_path):
            print("%d already done" % year)
            continue
        start = time.time()
        source = "/vsizip//vsicurl/" + (Change_zip % (Version, year)) + "/" + (Change_tif % (Version, year))
        counts = np.zeros(256, np.int64)

        # Written under a working name and renamed when whole, so a break leaves no half file.
        part_path = out_path + ".part"
        try:
            with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_MULTIRANGE="YES",
                              VSI_CACHE="TRUE", VSI_CACHE_SIZE="200000000", GDAL_HTTP_MAX_RETRY="3",
                              GDAL_HTTP_RETRY_DELAY="5"):
                with rasterio.open(source) as src:
                    # A class, so the nearest cell is taken rather than an average.
                    with WarpedVRT(src, crs=Maine_crs, transform=grid, width=width, height=height,
                                   resampling=Resampling.nearest, src_nodata=src.nodata,
                                   nodata=Lcms_nodata) as vrt:
                        with rasterio.open(part_path, "w", **profile) as dst:
                            for row in range(0, height, Download_strip):
                                rows = min(Download_strip, height - row)
                                window = Window(0, row, width, rows)
                                block = vrt.read(1, window=window)
                                dst.write(block, 1, window=window)
                                counts += np.bincount(block.ravel(), minlength=256)
            os.replace(part_path, out_path)
        except Exception as trouble:
            if os.path.exists(part_path):
                os.remove(part_path)
            print("%d failed, run the script again to retry it: %s" % (year, trouble))
            continue
        named = [(counts[c], c) for c in Loss_classes if counts[c] > 0]
        total = counts.sum() - counts[Lcms_nodata]
        print("%d written in %.0f s, %s" % (year, time.time() - start, out_path))
        for n, c in sorted(named, reverse=True):
            print("    %-28s %9d cells, %.3f%% of the state" %
                  (Change_names[c], n, 100.0 * n / max(total, 1)))

    # One land cover map per year, on the Maine grid.
    profile = {"driver": "GTiff", "height": height, "width": width, "count": 1, "dtype": "uint8",
               "crs": Maine_crs, "nodata": Lcms_nodata, "transform": grid,
               "compress": "deflate", "tiled": True, "blockxsize": 256, "blockysize": 256,
               "BIGTIFF": "IF_SAFER"}

    for year in Years:
        out_path = os.path.join(Cover_dir, "%d_Cover_Maine.tif" % year)
        if os.path.exists(out_path):
            print("%d already done" % year)
            continue
        start = time.time()
        source = "/vsizip//vsicurl/" + (Cover_zip % (Version, year)) + "/" + (Cover_tif % (Version, year))
        counts = np.zeros(256, np.int64)

        # Written under a working name and renamed when whole, so a break leaves no half file.
        part_path = out_path + ".part"
        try:
            with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_MULTIRANGE="YES",
                              VSI_CACHE="TRUE", VSI_CACHE_SIZE="200000000", GDAL_HTTP_MAX_RETRY="3",
                              GDAL_HTTP_RETRY_DELAY="5"):
                with rasterio.open(source) as src:
                    # A class, so the nearest cell is taken rather than an average.
                    with WarpedVRT(src, crs=Maine_crs, transform=grid, width=width, height=height,
                                   resampling=Resampling.nearest, src_nodata=src.nodata, nodata=Lcms_nodata) as vrt:
                        with rasterio.open(part_path, "w", **profile) as dst:
                            for row in range(0, height, Download_strip):
                                rows = min(Download_strip, height - row)
                                window = Window(0, row, width, rows)
                                block = vrt.read(1, window=window)
                                block[block == Skip_class] = Lcms_nodata
                                dst.write(block, 1, window=window)
                                counts += np.bincount(block.ravel(), minlength=256)
            os.replace(part_path, out_path)
        except Exception as trouble:
            if os.path.exists(part_path):
                os.remove(part_path)
            print("%d failed, run the script again to retry it: %s" % (year, trouble))
            continue
        total = counts.sum() - counts[Lcms_nodata]
        print("%d written in %.0f s, %s" % (year, time.time() - start, out_path))
        for code in sorted(Cover_names):
            if counts[code] > 0:
                print("    %-24s %11d cells, %5.2f%% of the state" %
                      (Cover_names[code], counts[code], 100.0 * counts[code] / max(total, 1)))

# -----------------------------------------------------------------------------
# Section 2. Reconstruction of annual maps.
# -----------------------------------------------------------------------------
if reconstruction:

    # Totals collected while the rasters are read.
    edges = np.arange(-Bin_limit, Bin_limit + Bin_width, Bin_width)
    n_bins = len(edges) + 1
    centres = np.concatenate([[edges[0] - Bin_width / 2], edges[:-1] + Bin_width / 2,
                              [edges[-1] + Bin_width / 2]])

    # Median from a histogram, which avoids holding every cell in memory.
    def histogram_median(counts):
        total = counts.sum()
        if total == 0:
            return np.nan
        return float(centres[np.searchsorted(np.cumsum(counts), total / 2.0)])

    # Counts go into one store, keyed by an integer, so a single pass serves every table.
    def collect(store, keys, change, start, index):
        for code in np.unique(keys):
            pick = keys == code
            entry = store.setdefault(int(code), [0, 0.0, 0.0, np.zeros(n_bins, np.int64)])
            entry[0] += int(pick.sum())
            entry[1] += float(change[pick].sum())
            entry[2] += float(start[pick].sum())
            entry[3] += np.bincount(index[pick], minlength=n_bins)

    rows = []
    events = []
    growth = []
    detection = []
    for first, second in Pairs:
        span = second - first
        agb_paths = [os.path.join(Flight_agb_dir, "%d_AGB_Maine.tif" % y) for y in (first, second)]
        years = list(range(first + 1, second + 1))
        change_paths = [os.path.join(Lcms_dir, "%d_Change_Maine.tif" % y) for y in years]
        if not all(os.path.exists(p) for p in agb_paths + change_paths):
            print("%d to %d: a map is missing, skipped" % (first, second))
            continue
        with rasterio.open(agb_paths[0]) as src:
            height, width = src.height, src.width

        # One store per table, all filled in the same pass over the rasters.
        by_since = {}
        by_events = {}
        by_stock = {}
        flagged_loss = np.zeros((len(Loss_steps), 2), np.int64)
        flagged_total = np.zeros(2, np.int64)
        for row in range(0, height, Loss_strip):
            take = min(Loss_strip, height - row)
            window = Window(0, row, width, take)
            with rasterio.open(agb_paths[0]) as src:
                before = src.read(1, window=window)
            with rasterio.open(agb_paths[1]) as src:
                after = src.read(1, window=window)
            middle = (before + after) / 2.0
            keep = ((before != Flight_nodata) & (after != Flight_nodata) &
                    ((before >= Min_start) | (middle >= Min_start)))
            if not keep.any():
                continue

            # The first loss sets the class and date, and the count says how often it was entered.
            agent_map = np.zeros(before.shape, np.uint8)
            year_map = np.zeros(before.shape, np.int16)
            entries = np.zeros(before.shape, np.int16)
            wood = np.zeros(before.shape, bool)
            for year, path in zip(years, change_paths):
                with rasterio.open(path) as src:
                    year_class = src.read(1, window=window)
                is_loss = np.isin(year_class, Loss_classes)
                fresh = is_loss & (agent_map == 0)
                agent_map[fresh] = year_class[fresh]
                year_map[fresh] = year
                is_wood = np.isin(year_class, Wood_loss_classes)
                entries += is_wood
                wood |= is_wood
            change = (after - before)[keep]
            agent = agent_map[keep]
            since = np.where(year_map[keep] > 0, second - year_map[keep], 0).astype(np.int16)
            start = before[keep]
            mid = middle[keep]
            index = np.digitize(change, edges)

            # The loss tables keep the original choice of cells, so they stay comparable.
            on_first = start >= Min_start
            collect(by_since, (agent.astype(np.int32) * 100 + since)[on_first],
                    change[on_first], start[on_first], index[on_first])
            collect(by_events, (agent.astype(np.int32) * 100 + np.minimum(entries[keep], 9))[on_first],
                    change[on_first], start[on_first], index[on_first])

            # Binned on the midpoint of the two flights, which carries no tie to their difference.
            quiet = (agent == 0) & (change > -Quiet_loss)
            for basis, level in ((0, start), (1, mid)):
                band = np.digitize(level, Stock_edges) - 1
                ok = (agent == 0) & (band >= 0) & (band < len(Stock_edges) - 1)
                collect(by_stock, np.where(ok, (basis * 100 + band) * 10 + quiet, -1),
                        change, level, index)
            marked = wood[keep][on_first].astype(np.int64)
            flagged_total += np.bincount(marked, minlength=2)
            for i, step in enumerate(Loss_steps):
                big = change[on_first] <= -step
                flagged_loss[i] += np.bincount(marked[big], minlength=2)

        for code, (n, change_sum, start_sum, counts) in by_since.items():
            key, since = divmod(code, 100)
            rows.append({"first year": first, "second year": second, "class": key,
                         "agent": Change_names.get(key, "no loss seen"), "years since loss": since,
                         "cells": n, "km2": n * Cell_km2, "start AGB": start_sum / n,
                         "mean change": change_sum / n, "mean per year": change_sum / n / span,
                         "median change": histogram_median(counts)})
        for code, (n, change_sum, start_sum, counts) in by_events.items():
            key, count = divmod(code, 100)
            events.append({"first year": first, "second year": second,
                           "agent": Change_names.get(key, "no loss seen"), "loss years": count,
                           "cells": n, "km2": n * Cell_km2, "start AGB": start_sum / n,
                           "mean change": change_sum / n,
                           "mean per loss year": change_sum / n / max(count, 1)})
        for code, (n, change_sum, level_sum, counts) in by_stock.items():
            if code < 0:
                continue
            key, quiet = divmod(code, 10)
            basis, band = divmod(key, 100)
            level = level_sum / n
            growth.append({"first year": first, "second year": second,
                           "binned by": "first flight" if basis == 0 else "midpoint",
                           "stock from": Stock_edges[band], "stock to": Stock_edges[band + 1],
                           "no quiet loss": bool(quiet), "cells": n, "km2": n * Cell_km2,
                           "mean stock": level, "mean per year": change_sum / n / span,
                           "rate per year": change_sum / n / span / max(level, 1.0)})
        for i, step in enumerate(Loss_steps):
            seen, missed = flagged_loss[i][1], flagged_loss[i][0]
            detection.append({"first year": first, "second year": second, "loss at least": step,
                              "cells": int(seen + missed), "flagged by the survey": int(seen),
                              "detection rate": seen / max(seen + missed, 1),
                              "cells flagged in all": int(flagged_total[1]),
                              "of those losing this much": seen / max(flagged_total[1], 1)})
        print("%d to %d done, %d cells followed" % (first, second, sum(v[0] for v in by_since.values())))

    # Save the tables.
    detail = pd.DataFrame(rows)
    per_event = pd.DataFrame(events)
    by_start = pd.DataFrame(growth)
    found = pd.DataFrame(detection)
    detail["change_total"] = detail["mean change"] * detail.cells
    detail["start_total"] = detail["start AGB"] * detail.cells
    by_agent = detail.groupby(["first year", "second year", "agent"], as_index=False).agg(
        cells=("cells", "sum"), km2=("km2", "sum"),
        change_total=("change_total", "sum"), start_total=("start_total", "sum"))
    by_agent["start AGB"] = by_agent.start_total / by_agent.cells
    by_agent["mean change"] = by_agent.change_total / by_agent.cells
    by_agent["mean per year"] = by_agent["mean change"] / (by_agent["second year"] - by_agent["first year"])
    by_agent = by_agent.drop(columns=["change_total", "start_total"])
    big = detail[detail.cells > 5000]
    by_since = big.pivot_table(index=["agent", "years since loss"], columns="first year",
                               values="mean change")
    print("\nbiomass change by change class and pair, Mg ha-1 per year:")
    print(by_agent.pivot_table(index="agent", columns="first year", values="mean per year").round(2).to_string())
    print("\nbiomass change against time since the loss, Mg ha-1 in total:")
    print(by_since.round(1).to_string())
    print("\nloss against the number of years the survey flagged, Mg ha-1 per loss year:")
    wood = per_event[(per_event["loss years"] > 0) & (per_event.cells > 5000)]
    print(wood.pivot_table(index=["agent", "loss years"], columns="first year",
                           values="mean per loss year").round(1).to_string())
    quiet = by_start[by_start["no quiet loss"]]
    print("\ngrowth on cells that also held their biomass, as a fraction of stock per year:")
    print(quiet.pivot_table(index=["stock from", "stock to"], columns=["binned by", "first year"],
                            values="rate per year").round(4).to_string())

    # Quiet losses are kept in, because the rate is applied to every cell the survey did not flag.
    pooled = by_start.copy()
    pooled["change_total"] = pooled["mean per year"] * pooled.cells
    pooled["stock_total"] = pooled["mean stock"] * pooled.cells
    pooled = pooled.groupby(["first year", "second year", "binned by", "stock from", "stock to"],
                            as_index=False).agg(cells=("cells", "sum"), km2=("km2", "sum"),
                                                change_total=("change_total", "sum"),
                                                stock_total=("stock_total", "sum"))
    pooled["mean stock"] = pooled.stock_total / pooled.cells
    pooled["mean per year"] = pooled.change_total / pooled.cells
    pooled["rate per year"] = pooled["mean per year"] / pooled["mean stock"]
    pooled = pooled.drop(columns=["change_total", "stock_total"])
    print("\nnet change on every cell the survey never flagged, as a fraction of stock per year:")
    print(pooled.pivot_table(index=["stock from", "stock to"], columns=["binned by", "first year"],
                             values="rate per year").round(4).to_string())
    print("\nthe same net change, Mg ha-1 per year:")
    print(pooled.pivot_table(index=["stock from", "stock to"], columns=["binned by", "first year"],
                             values="mean per year").round(2).to_string())
    print("\ncells behind each class, millions:")
    print((pooled.pivot_table(index=["stock from", "stock to"], columns=["binned by", "first year"],
                              values="cells") / 1e6).round(2).to_string())
    print("\nhow much of the real loss the survey flagged:")
    print(found.round(3).to_string(index=False))
    os.makedirs(Stats_dir, exist_ok=True)
    out_path = os.path.join(Stats_dir, Loss_file)
    with pd.ExcelWriter(out_path) as writer:
        by_agent.to_excel(writer, sheet_name="By Agent", index=False)
        detail.to_excel(writer, sheet_name="By Loss Age", index=False)
        per_event.to_excel(writer, sheet_name="By Loss Years", index=False)
        by_start.to_excel(writer, sheet_name="Growth By Stock", index=False)
        pooled.to_excel(writer, sheet_name="Net By Stock", index=False)
        found.to_excel(writer, sheet_name="Detection", index=False)
    print("\nsaved %s" % out_path)

    # Open every map once.
    for folder in (Agb_dir, Volume_dir):
        os.makedirs(folder, exist_ok=True)
    agb_in = {y: rasterio.open(os.path.join(Flight_agb_dir, "%d_AGB_Maine.tif" % y)) for y in Anchor_years}
    vol_in = {y: rasterio.open(os.path.join(Flight_vol_dir, "%d_Vol_Maine.tif" % y)) for y in Anchor_years}
    change_in = {y: rasterio.open(os.path.join(Lcms_dir, "%d_Change_Maine.tif" % y)) for y in Years}
    sample = agb_in[Anchor_years[0]]
    height, width = sample.height, sample.width
    profile = {"driver": "GTiff", "height": height, "width": width, "count": 1, "dtype": "int16",
               "crs": sample.crs, "nodata": Annual_nodata, "transform": sample.transform,
               "compress": "deflate", "predictor": 2, "tiled": True,
               "blockxsize": 256, "blockysize": 256, "BIGTIFF": "IF_SAFER"}
    agb_out = {y: rasterio.open(os.path.join(Agb_dir, "%d_AGB_Maine.tif" % y), "w", **profile) for y in Years}
    vol_out = {y: rasterio.open(os.path.join(Volume_dir, "%d_Vol_Maine.tif" % y), "w", **profile) for y in Years}
    print("grid %d x %d, %d years, %d flights" % (width, height, len(Years), len(Anchor_years)))

    # Walk the state in strips and fill in every year.
    n_years = len(Years)
    started = time.time()
    filled = 0
    for row in range(0, height, Annual_strip):
        take = min(Annual_strip, height - row)
        window = Window(0, row, width, take)
        shape = (n_years, take * width)

        # Flights, held as not a number away from a flight year.
        agb_seen = np.full(shape, np.nan, np.float32)
        vol_seen = np.full(shape, np.nan, np.float32)
        for i, year in enumerate(Years):
            if year not in agb_in:
                continue
            a = agb_in[year].read(1, window=window).ravel()
            v = vol_in[year].read(1, window=window).ravel()
            good = (a != Flight_nodata) & (v != Flight_nodata)
            agb_seen[i][good] = a[good]
            vol_seen[i][good] = v[good]
        if not np.isfinite(agb_seen).any():
            for year in Years:
                blank = np.full((take, width), Annual_nodata, np.int16)
                agb_out[year].write(blank, 1, window=window)
                vol_out[year].write(blank, 1, window=window)
            continue

        # The share taken out in each year, from the change class the survey gives.
        cut = np.zeros(shape, np.float32)
        for i, year in enumerate(Years):
            code = change_in[year].read(1, window=window).ravel()
            for key, share in Loss_fraction.items():
                cut[i][code == key] = share
            slow = np.isin(code, Slow_classes)
            if slow.any() and i + Slow_delay < n_years:
                cut[i + Slow_delay][slow] = np.maximum(cut[i + Slow_delay][slow], Slow_fraction)

        # Forward from the first flight, then backward from the last.
        forward_agb = np.full(shape, np.nan, np.float32)
        forward_vol = np.full(shape, np.nan, np.float32)
        state_a = np.full(take * width, np.nan, np.float32)
        state_v = np.full(take * width, np.nan, np.float32)
        for i in range(n_years):
            live = np.isfinite(state_a)
            if live.any():
                step = growth_fraction(state_a[live])
                state_a[live] = state_a[live] * (1.0 + step) * (1.0 - cut[i][live])
                state_v[live] = state_v[live] * (1.0 + step) * (1.0 - cut[i][live])
            seen = np.isfinite(agb_seen[i])
            state_a[seen] = agb_seen[i][seen]
            state_v[seen] = vol_seen[i][seen]
            forward_agb[i] = state_a
            forward_vol[i] = state_v
        back_agb = np.full(shape, np.nan, np.float32)
        back_vol = np.full(shape, np.nan, np.float32)
        state_a = np.full(take * width, np.nan, np.float32)
        state_v = np.full(take * width, np.nan, np.float32)
        for i in range(n_years - 1, -1, -1):
            seen = np.isfinite(agb_seen[i])
            state_a[seen] = agb_seen[i][seen]
            state_v[seen] = vol_seen[i][seen]
            back_agb[i] = state_a
            back_vol[i] = state_v
            live = np.isfinite(state_a)
            if live.any() and i > 0:
                # Undo the year just passed, so the state moves back to the year before it.
                opened = np.maximum(1.0 - cut[i][live], 0.2)
                back_a = state_a[live] / opened
                back_v = state_v[live] / opened
                step = growth_fraction(back_a)
                state_a[live] = np.minimum(back_a / (1.0 + step), Max_stock)
                state_v[live] = np.minimum(back_v / (1.0 + step), Max_stock * Volume_ratio)

        # Distance in years to the nearest flight on each side, used to weight the two passes.
        gap_forward = np.full(shape, n_years, np.int8)
        gap_back = np.full(shape, n_years, np.int8)
        run = np.full(take * width, n_years, np.int8)
        for i in range(n_years):
            run = np.where(np.isfinite(agb_seen[i]), 0, np.minimum(run + 1, n_years)).astype(np.int8)
            gap_forward[i] = run
        run = np.full(take * width, n_years, np.int8)
        for i in range(n_years - 1, -1, -1):
            run = np.where(np.isfinite(agb_seen[i]), 0, np.minimum(run + 1, n_years)).astype(np.int8)
            gap_back[i] = run

        # The two passes are blended year by year, so only one year sits in memory at a time.
        for i, year in enumerate(Years):
            weight_f = np.where(np.isfinite(forward_agb[i]), 1.0 / (gap_forward[i] + 1.0), 0.0)
            weight_b = np.where(np.isfinite(back_agb[i]), 1.0 / (gap_back[i] + 1.0), 0.0)
            total = weight_f + weight_b
            for ahead, behind, store in ((forward_agb[i], back_agb[i], agb_out[year]),
                                         (forward_vol[i], back_vol[i], vol_out[year])):
                mixed = np.nan_to_num(ahead) * weight_f + np.nan_to_num(behind) * weight_b
                values = np.where(total > 0, mixed / np.maximum(total, 1e-9), np.nan)
                values = np.where(values < Min_stock, 0.0, values)
                out = np.where(np.isfinite(values), np.rint(values * Scale), Annual_nodata)
                store.write(np.clip(out, Annual_nodata, 32767).astype(np.int16).reshape(take, width), 1, window=window)
            if i == 0:
                filled += int((total > 0).sum())
        if (row // Annual_strip) % 50 == 0:
            done = (row + take) / height
            used = time.time() - started
            print("  %.0f%% of the state, used %dm, left about %dm" %
                  (100 * done, used // 60, (used / max(done, 1e-6) - used) // 60), flush=True)

    for store in list(agb_in.values()) + list(vol_in.values()) + list(change_in.values()):
        store.close()
    for store in list(agb_out.values()) + list(vol_out.values()):
        store.close()
    print("finished in %dm, %s and %s" % ((time.time() - started) // 60, Agb_dir, Volume_dir))

# -----------------------------------------------------------------------------
# Section 3. Delineation of forest area.
# -----------------------------------------------------------------------------
if delineation:

    # Open every map once.
    cover_in = {}
    carbon_in = {}
    volume_in = {}
    for year in Years:
        cover_path = os.path.join(Cover_dir, "%d_Cover_Maine.tif" % year)
        carbon_path = os.path.join(Agb_dir, "%d_AGB_Maine.tif" % year)
        volume_path = os.path.join(Volume_dir, "%d_Vol_Maine.tif" % year)
        if not all(os.path.exists(p) for p in (cover_path, carbon_path, volume_path)):
            raise SystemExit("%d is missing a cover, carbon or volume map" % year)
        cover_in[year] = rasterio.open(cover_path)
        carbon_in[year] = rasterio.open(carbon_path)
        volume_in[year] = rasterio.open(volume_path)
    sample = carbon_in[Years[0]]
    height, width = sample.height, sample.width
    profile = {"driver": "GTiff", "height": height, "width": width, "count": 1, "dtype": "uint8",
               "crs": sample.crs, "nodata": 0, "transform": sample.transform,
               "compress": "deflate", "tiled": True, "blockxsize": 256, "blockysize": 256,
               "BIGTIFF": "IF_SAFER"}
    print("grid %d x %d, %d years" % (width, height, len(Years)))

    # Walk the state in strips and decide each cell.
    os.makedirs(Stats_dir, exist_ok=True)
    os.makedirs(Layer_dir, exist_ok=True)
    totals = {year: np.zeros(4) for year in Years}
    dropped_class = np.zeros(256, np.int64)
    kept_cells = 0
    all_cells = 0
    started = time.time()
    with rasterio.open(Mask_path, "w", **profile) as mask_out:
        for row in range(0, height, Mask_strip):
            take = min(Mask_strip, height - row)
            window = Window(0, row, width, take)

            # The survey is read year by year, so only the running answers are held.
            ever_forest = np.zeros((take, width), bool)
            seen_cover = np.zeros((take, width), bool)
            last_code = np.zeros((take, width), np.uint8)
            for year in Years:
                code = cover_in[year].read(1, window=window)
                good = code != Lcms_nodata
                seen_cover |= good
                ever_forest |= np.isin(code, Forest_classes)
                last_code = np.where(good, code, last_code)

            # The maps are held, because the yearly totals are taken from them after the decision.
            carbon = np.empty((len(Years), take, width), np.int16)
            volume = np.empty((len(Years), take, width), np.int16)
            for i, year in enumerate(Years):
                carbon[i] = carbon_in[year].read(1, window=window)
                volume[i] = volume_in[year].read(1, window=window)
            mapped = (carbon != Annual_nodata).any(axis=0)
            # Gaps hold the lowest possible number, so the highest year is always a real one.
            top_carbon = carbon.max(axis=0) / Scale

            # Drop a cell only when it was never forest and the maps found nothing there.
            drop = ~ever_forest & seen_cover & (top_carbon < Min_carbon)
            keep = mapped & ~drop
            mask_out.write(keep.astype(np.uint8), 1, window=window)
            all_cells += int(mapped.sum())
            kept_cells += int(keep.sum())
            dropped_class += np.bincount(last_code[mapped & drop], minlength=256)
            for i, year in enumerate(Years):
                good = carbon[i] != Annual_nodata
                run = totals[year]
                run += [int((good & keep).sum()),
                        float(carbon[i][good & keep].sum()) / Scale,
                        float(volume[i][good & keep].sum()) / Scale,
                        float(carbon[i][good].sum()) / Scale]
            if (row // Mask_strip) % 100 == 0:
                done = (row + take) / height
                used = time.time() - started
                print("  %.0f%% of the state, used %dm, left about %dm" %
                      (100 * done, used // 60, (used / max(done, 1e-6) - used) // 60), flush=True)
    for store in list(cover_in.values()) + list(carbon_in.values()) + list(volume_in.values()):
        store.close()

    # What the mask changed.
    print("\nmapped cells %d over %.2f million ha" % (all_cells, all_cells * Cell_ha / 1e6))
    print("kept as forest %d over %.2f million ha, %.1f%% of the mapped area"
          % (kept_cells, kept_cells * Cell_ha / 1e6, 100.0 * kept_cells / max(all_cells, 1)))
    print("\nwhat the survey calls the ground that was removed, by its latest class:")
    removed = dropped_class.sum()
    for code in np.argsort(dropped_class)[::-1]:
        if dropped_class[code] == 0:
            break
        print("    %-20s %11d cells, %5.1f%% of what was removed, %.0f thousand ha"
              % (Cover_names.get(int(code), str(int(code))), dropped_class[code],
                 100.0 * dropped_class[code] / max(removed, 1), dropped_class[code] * Cell_ha / 1e3))
    series = []
    for year in Years:
        n, carbon_sum, volume_sum, carbon_all = totals[year]
        series.append({"year": year, "forest cells": int(n), "forest Mha": n * Cell_ha / 1e6,
                       "carbon mean": carbon_sum / max(n, 1),
                       "carbon total Mg": carbon_sum * Cell_ha,
                       "volume mean": volume_sum / max(n, 1),
                       "volume total m3": volume_sum * Cell_ha,
                       "carbon total before mask Mg": carbon_all * Cell_ha})
    yearly = pd.DataFrame(series)
    yearly["share of carbon kept"] = yearly["carbon total Mg"] / yearly["carbon total before mask Mg"]
    print("\nstatewide series on forest land:")
    print(yearly.round(2).to_string(index=False))
    out_path = os.path.join(Stats_dir, Mask_table)
    with pd.ExcelWriter(out_path) as writer:
        yearly.to_excel(writer, sheet_name="Yearly Series", index=False)
    print("\nsaved %s and %s" % (Mask_path, out_path))

# -----------------------------------------------------------------------------
# Section 4. Map validation.
# -----------------------------------------------------------------------------
if validation:
    # The survey estimate comes from EVALIDator, so the maps are held against the number the
    # Forest Service publishes. The estimates are saved with the comparison they feed.
    os.makedirs(Stats_dir, exist_ok=True)

    # One estimate, for a year and an attribute, as a state total or split by county.
    def evalidator(year, attribute, rows):
        reply = requests.get(Evalidator_url,
                             params={"wc": Eval_group % year, "snum": attribute,
                                     "rselected": rows, "cselected": "None",
                                     "outputFormat": "NJSON"}, timeout=180)
        reply.raise_for_status()

        # An evaluation that is not published yet comes back as a page rather than as data.
        found = reply.json()["estimates"] \
            if "json" in reply.headers.get("Content-Type", "") else []
        print("  %d %d %s, %d rows" % (year, attribute, rows, len(found)), flush=True)
        return found

    # The survey estimate for each year.
    survey = []
    for year in Years:
        area = evalidator(year, Area_attribute, State_rows)
        carbon = evalidator(year, Carbon_attribute, State_rows)
        volume = evalidator(year, Volume_attribute, State_rows)
        if not (area and carbon and volume):
            print("%d has no published evaluation, left out" % year)
            continue
        hectares = area[0]["ESTIMATE"] * Acre_to_ha
        stock = carbon[0]["ESTIMATE"] * Ston_to_Mg
        wood = volume[0]["ESTIMATE"] * Cf_to_m3
        survey.append({"year": year, "plots": area[0]["PLOT_COUNT"],
                       "FIA forest Mha": hectares / 1e6, "FIA AGB Mg": stock,
                       "FIA AGB per ha": stock / hectares, "FIA volume m3": wood,
                       "FIA volume per ha": wood / hectares,
                       "SE% area": area[0]["SE_PERCENT"], "SE% AGB": carbon[0]["SE_PERCENT"],
                       "SE% volume": volume[0]["SE_PERCENT"]})
    survey = pd.DataFrame(survey)
    print("\nthe survey estimate for Maine, from EVALIDator:")
    print(survey.round(3).to_string(index=False))

    # The same totals from the maps, on forest land.
    mask_in = rasterio.open(Mask_path) if os.path.exists(Mask_path) else None
    if mask_in is None:
        print("\nno forest mask found, the map totals cover every mapped cell")
    maps = []
    for year in Years:
        agb_path = os.path.join(Agb_dir, "%d_AGB_Maine.tif" % year)
        vol_path = os.path.join(Volume_dir, "%d_Vol_Maine.tif" % year)
        if not (os.path.exists(agb_path) and os.path.exists(vol_path)):
            continue
        cells = 0
        agb_run = 0.0
        vol_run = 0.0
        with rasterio.open(agb_path) as agb_src, rasterio.open(vol_path) as vol_src:
            for row in range(0, agb_src.height, State_strip):
                take = min(State_strip, agb_src.height - row)
                window = Window(0, row, agb_src.width, take)
                agb = agb_src.read(1, window=window)
                vol = vol_src.read(1, window=window)
                good = agb != Annual_nodata
                if mask_in is not None:
                    good &= mask_in.read(1, window=window) == 1
                cells += int(good.sum())
                agb_run += float(agb[good].sum()) / Scale
                vol_run += float(vol[good].sum()) / Scale
        maps.append({"year": year, "map forest Mha": cells * Cell_ha / 1e6,
                     "map AGB Mg": agb_run * Cell_ha, "map AGB per ha": agb_run / max(cells, 1),
                     "map volume m3": vol_run * Cell_ha,
                     "map volume per ha": vol_run / max(cells, 1)})
        print("%d read, %.2f Mha, %.1f million Mg C, %.2f billion m3"
              % (year, cells * Cell_ha / 1e6, agb_run * Cell_ha / 1e6, vol_run * Cell_ha / 1e9),
              flush=True)
    maps = pd.DataFrame(maps)

    # The difference between them.
    both = maps.merge(survey, on="year", how="outer").sort_values("year")
    both["AGB ratio"] = both["map AGB Mg"] / both["FIA AGB Mg"]
    both["volume ratio"] = both["map volume m3"] / both["FIA volume m3"]
    both["area ratio"] = both["map forest Mha"] / both["FIA forest Mha"]
    both["AGB density ratio"] = both["map AGB per ha"] / both["FIA AGB per ha"]
    both["volume density ratio"] = both["map volume per ha"] / both["FIA volume per ha"]
    show = ["year", "map forest Mha", "FIA forest Mha", "area ratio", "map AGB per ha",
            "FIA AGB per ha", "AGB density ratio", "map volume per ha", "FIA volume per ha",
            "volume density ratio"]
    print("\nthe maps against the survey, statewide:")
    print(both[show].round(3).to_string(index=False))
    paired = both.dropna(subset=["AGB ratio"])
    for name in ("area ratio", "AGB density ratio", "volume density ratio",
                 "AGB ratio", "volume ratio"):
        print("%-21s mean %.3f, from %.3f to %.3f"
              % (name, paired[name].mean(), paired[name].min(), paired[name].max()))

    # County boundaries on the map grid, so the map totals split the way the survey does.
    os.makedirs(Layer_dir, exist_ok=True)
    county_path = os.path.join(Layer_dir, County_file)
    if not os.path.exists(county_path):
        reply = requests.get(County_url, params={"where": "STATE='23'", "outFields": "NAME",
                                                 "returnGeometry": "true", "outSR": "4326",
                                                 "f": "geojson"}, timeout=120)
        reply.raise_for_status()
        with open(county_path, "w") as f:
            f.write(reply.text)
    with open(county_path) as f:
        borders = json.load(f)["features"]
    shapes = []
    names = {}
    for number, border in enumerate(borders, start=1):
        names[number] = border["properties"]["NAME"].replace(" County", "").strip()
        geometry = border["geometry"]
        polygons = geometry["coordinates"] if geometry["type"] == "MultiPolygon" \
            else [geometry["coordinates"]]
        moved = []
        for polygon in polygons:
            rings = []
            for ring in polygon:
                xs, ys = transform("EPSG:4326", Maine_crs, [p[0] for p in ring],
                                   [p[1] for p in ring])
                rings.append(list(zip(xs, ys)))
            moved.append(rings)
        shapes.append(({"type": "MultiPolygon", "coordinates": moved}, number))
    print("\n%d counties on the map grid" % len(shapes))

    # The map totals of the county year, split by county.
    agb_path = os.path.join(Agb_dir, "%d_AGB_Maine.tif" % County_year)
    vol_path = os.path.join(Volume_dir, "%d_Vol_Maine.tif" % County_year)
    county_cells = np.zeros(len(shapes) + 1)
    county_agb = np.zeros(len(shapes) + 1)
    county_vol = np.zeros(len(shapes) + 1)
    with rasterio.open(agb_path) as agb_src, rasterio.open(vol_path) as vol_src:
        for row in range(0, agb_src.height, State_strip):
            take = min(State_strip, agb_src.height - row)
            window = Window(0, row, agb_src.width, take)
            grid = from_origin(agb_src.transform.c, agb_src.transform.f - row * Cell, Cell, Cell)
            number = rasterize(shapes, out_shape=(take, agb_src.width), transform=grid,
                               fill=0, dtype="uint8")
            agb = agb_src.read(1, window=window)
            vol = vol_src.read(1, window=window)
            good = (agb != Annual_nodata) & (number > 0)
            if mask_in is not None:
                good &= mask_in.read(1, window=window) == 1
            county_cells += np.bincount(number[good], minlength=len(shapes) + 1)
            county_agb += np.bincount(number[good], weights=agb[good] / Scale,
                                      minlength=len(shapes) + 1)
            county_vol += np.bincount(number[good], weights=vol[good] / Scale,
                                      minlength=len(shapes) + 1)
    if mask_in is not None:
        mask_in.close()

    # The survey totals of the same year, by county. The label carries the name last.
    rows = {}
    for attribute, column, factor in [(Carbon_attribute, "FIA AGB Mg", Ston_to_Mg),
                                      (Volume_attribute, "FIA volume m3", Cf_to_m3),
                                      (Area_attribute, "FIA forest Mha", Acre_to_ha / 1e6)]:
        for entry in evalidator(County_year, attribute, County_rows):
            name = entry["GRP1"].split(" ME ")[-1].strip()
            rows.setdefault(name, {"county": name})[column] = entry["ESTIMATE"] * factor
            if attribute == Carbon_attribute:
                rows[name]["plots"] = entry["PLOT_COUNT"]
                rows[name]["SE% AGB"] = entry["SE_PERCENT"]
    by_county = pd.DataFrame(list(rows.values()))
    drawn = pd.DataFrame({"county": [names[n] for n in range(1, len(shapes) + 1)],
                          "map forest Mha": county_cells[1:] * Cell_ha / 1e6,
                          "map AGB Mg": county_agb[1:] * Cell_ha,
                          "map volume m3": county_vol[1:] * Cell_ha})
    counties = drawn.merge(by_county, on="county", how="inner").sort_values("county")
    counties["AGB ratio"] = counties["map AGB Mg"] / counties["FIA AGB Mg"]
    counties["volume ratio"] = counties["map volume m3"] / counties["FIA volume m3"]
    print("\nthe maps against the survey by county, %d:" % County_year)
    print(counties.round(3).to_string(index=False))
    for name in ("AGB ratio", "volume ratio"):
        print("county %s from %.3f to %.3f"
              % (name, counties[name].min(), counties[name].max()))

    # Save the tables.
    out_path = os.path.join(Stats_dir, State_file)
    with pd.ExcelWriter(out_path) as writer:
        both.to_excel(writer, sheet_name="Comparison", index=False)
        counties.to_excel(writer, sheet_name="Counties", index=False)
        survey.to_excel(writer, sheet_name="FIA", index=False)
        maps.to_excel(writer, sheet_name="Maps", index=False)
    print("\nsaved %s" % out_path)

# -----------------------------------------------------------------------------
# Section 5. Figure generation.
# -----------------------------------------------------------------------------
if figure:

    # Map layers, county boundaries and aerial photo.
    os.makedirs(Layer_dir, exist_ok=True)
    os.makedirs(Figure_dir, exist_ok=True)
    west, east, south, north = Map_extent
    width, height = Map_size
    grid = from_bounds(west, south, east, north, width, height)

    county_path = os.path.join(Layer_dir, County_file)
    if not os.path.exists(county_path):
        reply = requests.get(County_url, params={"where": "STATE='23'", "outFields": "NAME", "returnGeometry": "true",
                                                 "outSR": "4326", "f": "geojson"}, timeout=120)
        reply.raise_for_status()
        with open(county_path, "w") as f:
            f.write(reply.text)
    with open(county_path) as f:
        counties = json.load(f)["features"]
    print("%d counties" % len(counties))

    # Boundary rings for drawing, and a state mask so the maps stop at the state line.
    rings = []
    for county in counties:
        geometry = county["geometry"]
        polygons = geometry["coordinates"] if geometry["type"] == "MultiPolygon" else [geometry["coordinates"]]
        for polygon in polygons:
            for ring in polygon:
                rings.append(np.array(ring))
    inside = rasterize([county["geometry"] for county in counties], out_shape=(height, width), transform=grid) == 1

    # The rendered extent is kept in a world file beside the image.
    imagery_path = os.path.join(Layer_dir, Imagery_file)
    if not os.path.exists(imagery_path):
        imagery_height = int(round(Imagery_width * (north - south) / (east - west)))
        window = {"bbox": "%s,%s,%s,%s" % (west, south, east, north), "bboxSR": "4326", "imageSR": "4326",
                  "size": "%d,%d" % (Imagery_width, imagery_height), "format": "jpg"}
        info = requests.get(Imagery_url, params=dict(window, f="json"), timeout=300).json()
        reply = requests.get(Imagery_url, params=dict(window, f="image"), timeout=300)
        reply.raise_for_status()
        with open(imagery_path, "wb") as f:
            f.write(reply.content)
        box = info["extent"]
        x_size = (box["xmax"] - box["xmin"]) / info["width"]
        y_size = (box["ymax"] - box["ymin"]) / info["height"]
        with open(os.path.splitext(imagery_path)[0] + ".jgw", "w") as f:
            f.write("%.10f\n0\n0\n%.10f\n%.10f\n%.10f\n" %
                    (x_size, -y_size, box["xmin"] + x_size / 2, box["ymax"] - y_size / 2))
    with rasterio.open(imagery_path) as src:
        imagery = src.read().transpose(1, 2, 0)
        imagery_extent = (src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top)
    print("aerial photo %d x %d, window %.3f to %.3f east, %.3f to %.3f north" %
          (imagery.shape[1], imagery.shape[0], imagery_extent[0], imagery_extent[1], imagery_extent[2], imagery_extent[3]))

    # Read the annual maps onto the map grid.
    # Each map is averaged to coarse cells, then onto the map grid, and cached until a map changes.
    paths = {}
    for tag, folder, pattern, _, _, _, _ in Products:
        # A four digit year only, so a per project map is never read as a year.
        for path in sorted(glob.glob(os.path.join(folder, pattern % "[0-9][0-9][0-9][0-9]"))):
            paths[tag + "_" + os.path.basename(path).split("_")[0]] = path
    if os.path.exists(Mask_path):
        paths["forest"] = Mask_path
    else:
        print("no forest mask found, the maps are drawn over every mapped cell")
    stack_path = os.path.join(Layer_dir, Stack_file)
    fresh = os.path.exists(stack_path) and os.path.getmtime(stack_path) > max(os.path.getmtime(p) for p in paths.values())
    if fresh:
        stack = np.load(stack_path)
        fresh = list(stack["extent"]) == list(Map_extent) and list(stack["size"]) == list(Map_size) and \
            all(key in stack for key in paths) and len(stack.files) == len(paths) + 2
    if not fresh:
        stack = {"extent": Map_extent, "size": Map_size}
        for key, path in paths.items():
            start = time.time()
            with rasterio.Env(GDAL_NUM_THREADS="ALL_CPUS"):
                with rasterio.open(path) as src:
                    coarse = src.read(1, out_shape=(src.height // Coarsen, src.width // Coarsen),
                                      resampling=Resampling.average, masked=True)
                    coarse_grid = src.transform * src.transform.scale(src.width / coarse.shape[1],
                                                                      src.height / coarse.shape[0])
                    # Whole numbers, so gaps become not a number only after a cast.
                    coarse = coarse.astype("float32").filled(np.nan)
                    values = np.full((height, width), np.nan, dtype="float32")
                    reproject(coarse, values, src_transform=coarse_grid, src_crs=src.crs,
                              src_nodata=np.nan, dst_transform=grid, dst_crs="EPSG:4326", dst_nodata=np.nan,
                              resampling=Resampling.average)
            # Averaging the mask gives the forest share of a map cell.
            stack[key] = values if key == "forest" else values / Scale
            print("%s read in %.0f s" % (os.path.basename(path), time.time() - start))
        np.savez_compressed(stack_path, **stack)

    wooded = np.ones((height, width), bool) if "forest" not in paths else \
        np.nan_to_num(stack["forest"]) >= Mask_share
    maps = {}
    for key in paths:
        if key == "forest":
            continue
        tag, year = key.split("_")
        values = np.ma.masked_where(~inside | ~wooded | np.isnan(stack[key]), stack[key])
        maps.setdefault(tag, {})[year] = values
        print("%s %s: %.0f%% of the state drawn, mean %.0f, 99th percentile %.0f" %
              (tag, year, 100.0 * values.count() / inside.sum(), values.mean(),
               np.percentile(values.compressed(), 99)))

    # Figure settings and the map panel.
    rcParams["font.family"] = "Times New Roman"
    rcParams["font.size"] = Font_size
    rcParams["text.color"] = "black"
    rcParams["axes.labelcolor"] = "black"
    rcParams["axes.edgecolor"] = "black"
    rcParams["xtick.color"] = "black"
    rcParams["ytick.color"] = "black"
    rcParams["mathtext.fontset"] = "custom"
    rcParams["mathtext.rm"] = "Times New Roman"
    rcParams["mathtext.it"] = "Times New Roman"
    rcParams["mathtext.bf"] = "Times New Roman"
    rcParams["mathtext.cal"] = "Times New Roman"
    rcParams["mathtext.default"] = "regular"

    # A degree of latitude is longer on the page than a degree of longitude at this latitude.
    stretch = 1.0 / np.cos(np.radians((south + north) / 2.0))

    # Draw one map panel and return the map image.
    def draw_panel(ax, values, low, high, colormap, lon_ticks, lat_ticks, text):
        ax.imshow(imagery, extent=imagery_extent, zorder=1)
        image = ax.imshow(values, extent=(west, east, south, north), cmap=colormap, vmin=low, vmax=high, zorder=2)
        for ring in rings:
            ax.plot(ring[:, 0], ring[:, 1], color="black", linewidth=0.3, zorder=3)
        ax.text(0.96, 0.04, text, transform=ax.transAxes, ha="right", va="bottom", zorder=4,
                color="white")
        ax.set_xlim(west, east)
        ax.set_ylim(south, north)
        ax.set_aspect(stretch)
        ax.set_xticks(lon_ticks)
        ax.set_yticks(lat_ticks)
        ax.set_xticklabels(["%d°W" % -lon for lon in lon_ticks])
        ax.set_yticklabels(["%d°N" % lat for lat in lat_ticks])
        return image

    # Figure 1, where the lidar flights reach.

    # Margins around the block of panels, in inches, shared by every figure below.
    left, right, bottom, top, gap = 0.55, 0.1, 0.45, 0.1, 0.12

    # Each flight year is carried onto the map grid, so a cell holds the share it reached.
    coverage_path = os.path.join(Layer_dir, Coverage_stack)
    flight_paths = {"%d" % year: os.path.join(Flight_agb_dir, "%d_AGB_Maine.tif" % year)
                    for year in Anchor_years}
    fresh = os.path.exists(coverage_path) and \
        os.path.getmtime(coverage_path) > max(os.path.getmtime(p) for p in flight_paths.values())
    if fresh:
        flights = np.load(coverage_path)
        fresh = list(flights["extent"]) == list(Map_extent) and list(flights["size"]) == list(Map_size) \
            and all(key in flights for key in flight_paths) and len(flights.files) == len(flight_paths) + 2
    if not fresh:
        flights = {"extent": Map_extent, "size": Map_size}
        for key, path in flight_paths.items():
            start = time.time()
            with rasterio.Env(GDAL_NUM_THREADS="ALL_CPUS"):
                with rasterio.open(path) as src:
                    coarse = src.read(1, out_shape=(src.height // Coarsen, src.width // Coarsen),
                                      resampling=Resampling.average, masked=True)
                    coarse_grid = src.transform * src.transform.scale(src.width / coarse.shape[1],
                                                                      src.height / coarse.shape[0])
                    covered = (~np.ma.getmaskarray(coarse)).astype("float32")
                    values = np.full((height, width), np.nan, dtype="float32")
                    reproject(covered, values, src_transform=coarse_grid, src_crs=src.crs,
                              src_nodata=np.nan, dst_transform=grid, dst_crs="EPSG:4326",
                              dst_nodata=np.nan, resampling=Resampling.average)
            flights[key] = values
            print("%s read in %.0f s" % (os.path.basename(path), time.time() - start))
        np.savez_compressed(coverage_path, **flights)

    # One panel per flight year, and a last panel holding every flight together.
    reached = {key: np.nan_to_num(flights[key]) >= Mask_share for key in flight_paths}
    every = np.zeros((height, width), bool)
    filled = np.ones((height, width), np.float32)
    panels = []
    for key in sorted(flight_paths):
        every |= reached[key]
        panels.append((key, np.ma.masked_where(~(inside & reached[key]), filled)))
        print("%s reaches %.1f%% of the state" %
              (key, 100.0 * (inside & reached[key]).sum() / inside.sum()))
    panels.append((Coverage_label, np.ma.masked_where(~(inside & every), filled)))
    print("every flight together reaches %.1f%% of the state"
          % (100.0 * (inside & every).sum() / inside.sum()))

    # The same grid of panels, with one fill colour because a panel only says reached or not.
    rows, cols = Coverage_layout
    panel_height = Panel_width * (north - south) / (east - west) * stretch
    figure_width = left + cols * Panel_width + (cols - 1) * gap + right
    figure_height = bottom + rows * panel_height + (rows - 1) * gap + top
    fig, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))
    fig.subplots_adjust(left=left / figure_width, right=1 - right / figure_width,
                        bottom=bottom / figure_height, top=1 - top / figure_height,
                        wspace=gap / Panel_width, hspace=gap / panel_height)
    for slot, ax in enumerate(axes.flat):
        if slot >= len(panels):
            ax.set_visible(False)
            continue
        label, values = panels[slot]
        draw_panel(ax, values, 0.0, 1.0, ListedColormap([Coverage_colour]), Lon_ticks, Lat_ticks,
                   label)
        # Tick labels only on the left column and on the lowest panel of each column.
        row, col = divmod(slot, cols)
        ax.tick_params(labelleft=col == 0, labelbottom=row == rows - 1 or slot + cols >= len(panels))
    fig.savefig(os.path.join(Figure_dir, Coverage_figure), dpi=Figure_dpi)
    plt.close(fig)
    print("%s saved, %d panels" % (Coverage_figure, len(panels)))

    # Figures 3 and 4, one panel per year.
    rows, cols = Layout
    panel_height = Panel_width * (north - south) / (east - west) * stretch
    figure_width = left + cols * Panel_width + (cols - 1) * gap + right
    figure_height = bottom + rows * panel_height + (rows - 1) * gap + top

    ramp = ListedColormap(plt.get_cmap(Colormap)(np.linspace(0.0, Colormap_top, 256)))
    for tag, _, _, figure_name, label, (low, high), step in Products:
        years = sorted(maps[tag])
        if len(years) >= rows * cols:
            raise SystemExit("%d years need a larger Layout, the colour bar takes one slot" % len(years))
        fig, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))
        fig.subplots_adjust(left=left / figure_width, right=1 - right / figure_width,
                            bottom=bottom / figure_height, top=1 - top / figure_height,
                            wspace=gap / Panel_width, hspace=gap / panel_height)
        for slot, ax in enumerate(axes.flat):
            if slot >= len(years):
                ax.set_visible(False)
                continue
            year = years[slot]
            image = draw_panel(ax, maps[tag][year], low, high, ramp, Lon_ticks, Lat_ticks, year)
            # Tick labels only on the left column and on the lowest panel of each column.
            row, col = divmod(slot, cols)
            ax.tick_params(labelleft=col == 0, labelbottom=row == rows - 1 or slot + cols >= len(years))
        # The colour bar lies flat across the slots the years leave empty in the last row.
        spare = axes.flat[len(years)].get_position()
        last = axes.flat[-1].get_position()
        bar_ax = fig.add_axes([spare.x0 + 0.05 * (last.x1 - spare.x0),
                               spare.y0 + 0.46 * spare.height,
                               0.90 * (last.x1 - spare.x0), 0.09 * spare.height])
        bar = fig.colorbar(image, cax=bar_ax, orientation="horizontal",
                           ticks=np.arange(low, high + step, step))
        bar.set_label(label)
        fig.savefig(os.path.join(Figure_dir, figure_name), dpi=Figure_dpi)
        plt.close(fig)
        print("%s saved, %d years" % (figure_name, len(years)))

    # Figure 5, the maps against the forest survey.
    rcParams["font.size"] = Validation_font
    both = pd.read_excel(os.path.join(Stats_dir, State_file), sheet_name="Comparison")
    by_county = pd.read_excel(os.path.join(Stats_dir, State_file), sheet_name="Counties")
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 8.2))

    # Map column, survey column, axis label, and the range and tick step of the axis.
    yearly = [("map AGB per ha", "FIA AGB per ha", "Aboveground biomass (Mg C ha$^{-1}$)",
               20, 60, 10),
              ("map volume per ha", "FIA volume per ha", "Stem volume (m$^{3}$ ha$^{-1}$)",
               150, 200, 10)]
    for ax, panel_label, (map_col, fia_col, label, low, high, step) in zip(
            axes[0], Panel_labels, yearly):
        good = both[fia_col].notna()
        ax.plot(both.year[good], both[fia_col][good], color=Fia_colour, linewidth=1.2,
                marker="o", markersize=4, label="FIA-derived")
        ax.plot(both.year, both[map_col], color=Map_colour, linewidth=1.2, linestyle="--",
                marker="s", markersize=4, label="Map-derived")
        ax.set_xlabel("Year")
        ax.set_ylabel(label)
        ax.set_xlim(Years[0] - 0.5, Years[-1] + 0.5)
        ax.set_xticks(range(Years[0], Years[-1] + 1, 3))
        ax.set_ylim(low, high)
        ax.set_yticks(np.arange(low, high + step, step))
        ax.set_box_aspect(1)
        ax.text(0.96, 0.04, panel_label, transform=ax.transAxes, ha="right", va="bottom")

    # Both panels carry the same two series, so one legend serves them.
    axes[0][1].legend(frameon=False, loc="upper left")

    # The county totals of one year, on axes that share a range so the line runs at 45 degrees.
    totals = [("map AGB Mg", "FIA AGB Mg", 1e6, "Aboveground biomass (million Mg C)", 0, 100, 20),
              ("map volume m3", "FIA volume m3", 1e6, "Stem volume (million m$^{3}$)", 0, 400, 100)]
    for ax, panel_label, (map_col, fia_col, scale, label, low, high, step) in zip(
            axes[1], Panel_labels[2:], totals):
        ax.plot([low, high], [low, high], color="black", linewidth=0.8)
        ax.scatter(by_county[fia_col] / scale, by_county[map_col] / scale, s=22,
                   facecolor="white", edgecolor="black", linewidth=0.8, zorder=3)
        ax.set_xlim(low, high)
        ax.set_ylim(low, high)
        ax.set_xticks(np.arange(low, high + step, step))
        ax.set_yticks(np.arange(low, high + step, step))
        ax.set_box_aspect(1)
        ax.set_xlabel("FIA - %s" % label)
        ax.set_ylabel("Map - %s" % label)
        ax.text(0.96, 0.04, panel_label, transform=ax.transAxes, ha="right", va="bottom")
    fig.tight_layout(h_pad=0.6)
    os.makedirs(Figure_dir, exist_ok=True)
    fig.savefig(os.path.join(Figure_dir, State_figure), dpi=Figure_dpi)
    plt.close(fig)
    print("%s saved" % State_figure)

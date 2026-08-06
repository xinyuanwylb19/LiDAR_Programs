"""
Data_Prepare.py
Generate sub-plots from original 500 m² circular inventory plots, assign trees,
compute plot-level AGB (oven-dry Mg ha-1), extract normalised ALS metrics, and save CSV.
Inventory data and plot centres are read from Demeritt_Inventory_Data.xlsx.

Figures (Pseudo_Plot_Illustration.py is merged in here):

    python Data_Prepare.py --figs1               # Figure S1, pseudo-plot layout

Also builds a 1 m canopy height model (CHM) GeoTIFF over the inventory plots
from the Demeritt_ALS tiles:

    python Data_Prepare.py --chm                 # extent from the plot table
    python Data_Prepare.py --chm X0 Y0 X1 Y1     # explicit extent

Author: Xinyuan Wei
"""

import os
import gc
import struct
import sys
import numpy as np
import pandas as pd

try:
    import laspy                 # used by the sub-plot pipeline further down;
except ImportError:              # the CHM section needs only numpy, so a
    laspy = None                 # missing laspy must not block --chm

base_dir = os.path.dirname(os.path.abspath(__file__))
ALS_dir = os.path.join(base_dir, "Demeritt_ALS")
save_file = os.path.join(base_dir, "Prepared_Dataset.csv")
results_dir = os.path.join(base_dir, "Results")
workbook = os.path.join(base_dir, "Demeritt_Inventory_Data.xlsx")


seed = 42
plot_radius = 12.62             # 500 m² circle
height_thresh = 2.0             # remove points < 2 m after normalisation
subplots_per_parent = 5         # evenly distributed pseudo-plots per parent plot
plot_size = [100, 150, 200, 225, 250, 300]  # m²

# carbon-mass columns in inventory (kg per tree, already converted to C).
# Four allometric methods; "_Mean" below is the mean of all four, matching the
# Biomass (m) column of the workbook.
allometry = {"Y": "Biomass (Y)", "J": "Biomass (J)",
             "C": "Biomass (C)", "W": "Biomass (W)"}

# Method-W carbon fraction by inventory species code, from the
# species-specific values of Westfall et al. (2024) (Allometry_W).
W_CARBON_FRACTION = {
    "ABBA": 0.505515,     "ACRU": 0.485733,     "ACSA": 0.485464,     "BEAL": 0.487783,
    "BEPA": 0.519087,     "BEPO": 0.477250,     "FAGR": 0.478167,     "FRNI": 0.480000,
    "PIAB": 0.478307,     "PIMA": 0.479700,     "PIRE": 0.532800,     "PIRU": 0.480050,
    "PIST": 0.507067,     "POBA": 0.482150,     "POGR": 0.480400,     "POTR": 0.479237,
    "QURU": 0.478312,     "THOC": 0.500743,     "TSCA": 0.479700,     "UNK": 0.505515,
}


# ------------------------------ CHM settings ---------------------------------
chm_res = 1.0                   # output cell size (m)
chm_buffer = 250.0              # margin around the outermost plot centres (m)
chm_max_height = 50.0           # cap; anything taller is treated as an outlier
chm_coarse = 10                 # block size (cells) for the coarse ground fill
chm_epsg = 32619                # WGS 84 / UTM 19N, read from the LAS geo keys
chm_file = os.path.join(results_dir, "Demeritt_CHM_1m.tif")
GROUND_CLASS = 2
NOISE_CLASSES = (7, 18)         # low / high noise


#==============================================================================
# PSEUDO-PLOT GEOMETRY
#
# Deterministic placement: pseudo-plots sit at equal angular intervals on a ring
# centred on the parent-plot centre, the ring radius chosen so each one stays
# fully inside the 500 m² parent. Shared by the dataset build and by Figure S1,
# so the illustration cannot drift from what is actually generated.
#==============================================================================
def even_offsets_circle(parent_r, sub_r, n):
    """Ring offsets keeping a circular pseudo-plot inside the parent circle."""
    max_r = parent_r - sub_r
    if max_r <= 0:
        return [(0.0, 0.0)] * n
    angles = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return [(max_r * np.cos(a), max_r * np.sin(a)) for a in angles]


def even_offsets_square(parent_r, half_side, n):
    """Ring offsets keeping all four square corners inside the parent circle."""
    max_r = parent_r - half_side * np.sqrt(2)
    if max_r <= 0:
        return [(0.0, 0.0)] * n
    angles = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return [(max_r * np.cos(a), max_r * np.sin(a)) for a in angles]


#==============================================================================
# CANOPY HEIGHT MODEL
#
# The tiles total ~22 GB, so points are streamed in chunks and reduced straight
# onto the output grid; nothing larger than one chunk is ever held in memory.
# LAS records are read directly with numpy (fixed-length records, classification
# at byte 15) and the GeoTIFF is written by hand, so this section needs only
# numpy and runs wherever the rest of the script does.
#==============================================================================
def las_header(path):
    """Public header block: bounds, scale/offset, point count and record size."""
    with open(path, "rb") as fh:
        b = fh.read(375)
    if b[:4] != b"LASF":
        raise ValueError(f"not a LAS file: {path}")
    h = {"version": (b[24], b[25])}
    h["header_size"], h["offset"] = struct.unpack("<HI", b[94:100])
    h["fmt"] = b[104] & 0b00111111
    h["reclen"], = struct.unpack("<H", b[105:107])
    n32, = struct.unpack("<I", b[107:111])
    (h["sx"], h["sy"], h["sz"],
     h["ox"], h["oy"], h["oz"]) = struct.unpack("<6d", b[131:179])
    (h["xmax"], h["xmin"], h["ymax"],
     h["ymin"], h["zmax"], h["zmin"]) = struct.unpack("<6d", b[179:227])
    if h["version"] >= (1, 4):
        n64, = struct.unpack("<Q", b[247:255])
        h["count"] = n64 or n32
    else:
        h["count"] = n32
    if h["fmt"] > 5:                     # classification moves in formats 6+
        raise ValueError(f"point format {h['fmt']} not supported")
    return h


def las_chunks(path, chunk=2_000_000):
    """Stream (x, y, z, classification, withheld) in georeferenced units."""
    h = las_header(path)
    dt = np.dtype({"names": ["X", "Y", "Z", "c"],
                   "formats": ["<i4", "<i4", "<i4", "u1"],
                   "offsets": [0, 4, 8, 15], "itemsize": h["reclen"]})
    left = h["count"]
    with open(path, "rb") as fh:
        fh.seek(h["offset"])
        while left > 0:
            n = min(chunk, left)
            raw = fh.read(n * h["reclen"])
            if len(raw) < n * h["reclen"]:
                n = len(raw) // h["reclen"]
                if n == 0:
                    break
            a = np.frombuffer(raw, dtype=dt, count=n)
            yield (a["X"] * h["sx"] + h["ox"],
                   a["Y"] * h["sy"] + h["oy"],
                   a["Z"] * h["sz"] + h["oz"],
                   a["c"] & 0b00011111,             # classification bits
                   (a["c"] & 0b10000000) != 0)      # withheld flag
            left -= n


def _scatter(target, idx, val, how):
    """Per-cell max or min into a flat grid.

    np.maximum.at would be the obvious call but is far too slow at this point
    count. Sorting by (cell, value) puts the winner of each cell at a known end
    of its run, so one deduplicated fancy-index assignment does the job."""
    if idx.size == 0:
        return
    order = np.lexsort((val, idx))
    i, v = idx[order], val[order]
    keep = np.empty(i.size, bool)
    if how == "max":                     # ascending value -> last of each run
        keep[-1] = True
        keep[:-1] = i[1:] != i[:-1]
        ii, vv = i[keep], v[keep]
        np.maximum(target[ii], vv, out=vv)
    else:                                # ascending value -> first of each run
        keep[0] = True
        keep[1:] = i[1:] != i[:-1]
        ii, vv = i[keep], v[keep]
        np.minimum(target[ii], vv, out=vv)
    target[ii] = vv


def _fill_nearest(a, max_iter=200):
    """Grow valid values outward into NaN cells using the 4-neighbour mean."""
    out = a.copy()
    for _ in range(max_iter):
        bad = np.isnan(out)
        if not bad.any():
            break
        s = np.zeros_like(out)
        n = np.zeros(out.shape, np.int16)
        for sh, ax in ((1, 0), (-1, 0), (1, 1), (-1, 1)):
            r = np.roll(out, sh, axis=ax)
            if ax == 0:
                r[0 if sh > 0 else -1, :] = np.nan
            else:
                r[:, 0 if sh > 0 else -1] = np.nan
            ok = ~np.isnan(r)
            s[ok] += r[ok]
            n[ok] += 1
        fill = bad & (n > 0)
        if not fill.any():
            break
        out[fill] = s[fill] / n[fill]
    return out


def _block_min(a, f):
    """NaN-aware minimum over f x f blocks, padding to a whole multiple.

    Blocks with no ground return at all come back as NaN, which is expected
    here and filled in afterwards, so the all-NaN warning is silenced."""
    import warnings
    r, c = a.shape
    p = np.pad(a, ((0, (-r) % f), (0, (-c) % f)), constant_values=np.nan)
    p = p.reshape(p.shape[0] // f, f, p.shape[1] // f, f)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmin(p, axis=(1, 3))


def _upsample(coarse, shape, f):
    """Bilinear expansion of a coarse grid back onto the fine grid."""
    nr, nc = shape
    rr = (np.arange(nr) + 0.5) / f - 0.5
    cc = (np.arange(nc) + 0.5) / f - 0.5
    r0 = np.clip(np.floor(rr).astype(np.int32), 0, coarse.shape[0] - 1)
    c0 = np.clip(np.floor(cc).astype(np.int32), 0, coarse.shape[1] - 1)
    r1 = np.clip(r0 + 1, 0, coarse.shape[0] - 1)
    c1 = np.clip(c0 + 1, 0, coarse.shape[1] - 1)
    wr = (rr - r0).astype(np.float32)[:, None]
    wc = (cc - c0).astype(np.float32)[None, :]
    top = coarse[np.ix_(r0, c0)] * (1 - wc) + coarse[np.ix_(r0, c1)] * wc
    bot = coarse[np.ix_(r1, c0)] * (1 - wc) + coarse[np.ix_(r1, c1)] * wc
    return top * (1 - wr) + bot * wr


def build_chm(las_files, extent, res=chm_res, max_height=chm_max_height,
              coarse=chm_coarse, log=print):
    """DSM (max return) minus DTM (min ground return) on a common grid.

    extent = (xmin, ymin, xmax, ymax). Returns (chm, (xmin, ymax))."""
    xmin, ymin, xmax, ymax = extent
    xmin, ymin = np.floor(xmin / res) * res, np.floor(ymin / res) * res
    xmax, ymax = np.ceil(xmax / res) * res, np.ceil(ymax / res) * res
    ncols = int(round((xmax - xmin) / res))
    nrows = int(round((ymax - ymin) / res))
    log(f"grid {ncols} x {nrows} cells at {res} m "
        f"({xmin:.0f}-{xmax:.0f} E, {ymin:.0f}-{ymax:.0f} N)")

    dsm = np.full(nrows * ncols, -np.inf, np.float32)
    dtm = np.full(nrows * ncols, np.inf, np.float32)

    used = total = 0
    for k, f in enumerate(las_files, 1):
        h = las_header(f)
        if (h["xmax"] < xmin or h["xmin"] > xmax
                or h["ymax"] < ymin or h["ymin"] > ymax):
            continue                                   # tile misses the extent
        used += 1
        npts = 0
        for x, y, z, cls, withheld in las_chunks(f):
            m = ((x >= xmin) & (x < xmax) & (y >= ymin) & (y < ymax)
                 & ~withheld
                 & (cls != NOISE_CLASSES[0]) & (cls != NOISE_CLASSES[1]))
            if not m.any():
                continue
            x, y, z, cls = x[m], y[m], z[m], cls[m]
            col = np.clip(((x - xmin) / res).astype(np.int64), 0, ncols - 1)
            row = np.clip(((ymax - y) / res).astype(np.int64), 0, nrows - 1)
            idx = row * ncols + col
            _scatter(dsm, idx, z.astype(np.float32), "max")
            g = cls == GROUND_CLASS
            if g.any():
                _scatter(dtm, idx[g], z[g].astype(np.float32), "min")
            npts += int(m.sum())
        total += npts
        log(f"  [{k}/{len(las_files)}] {os.path.basename(f)}: {npts:,} points")

    dsm = dsm.reshape(nrows, ncols)
    dtm = dtm.reshape(nrows, ncols)
    dsm[np.isinf(dsm)] = np.nan
    dtm[np.isinf(dtm)] = np.nan
    log(f"{used} tiles intersected, {total:,} points; "
        f"surface cells {np.isfinite(dsm).mean():.1%}, "
        f"ground cells {np.isfinite(dtm).mean():.1%}")

    # Ground returns are sparse under closed canopy, so gaps in the terrain are
    # filled from a coarse ground surface rather than left as holes.
    cg = _fill_nearest(_block_min(dtm, coarse))
    dtm = np.where(np.isnan(dtm),
                   _upsample(cg, dtm.shape, coarse)[:nrows, :ncols], dtm)

    chm = dsm - dtm
    chm[chm < 0] = 0.0                                 # ground noise below DTM
    hi = np.isfinite(chm) & (chm > max_height)
    if hi.any():
        log(f"clipped {int(hi.sum()):,} cells above {max_height} m")
        chm[hi] = max_height
    log(f"CHM {np.nanmin(chm):.2f}-{np.nanmax(chm):.2f} m, "
        f"mean {np.nanmean(chm):.2f} m")
    return chm, (xmin, ymax)


def write_geotiff(path, arr, origin, res=chm_res, epsg=chm_epsg, nodata=-9999.0):
    """Single-band float32 GeoTIFF, written directly so no GDAL is required."""
    arr = np.ascontiguousarray(np.where(np.isfinite(arr), arr, nodata), np.float32)
    nrows, ncols = arr.shape
    x0, y0 = origin
    geokeys = np.array([1, 1, 0, 3,            # directory version, 3 keys
                        1024, 0, 1, 1,         # GTModelType     = projected
                        1025, 0, 1, 1,         # GTRasterType    = PixelIsArea
                        3072, 0, 1, int(epsg)  # ProjectedCSType = EPSG
                        ], "<u2").tobytes()
    scale = struct.pack("<3d", res, res, 0.0)
    tie = struct.pack("<6d", 0.0, 0.0, 0.0, x0, y0, 0.0)
    nod = f"{nodata:g}".encode() + b"\0"

    tags = [(256, 4, 1, ncols), (257, 4, 1, nrows), (258, 3, 1, 32),
            (259, 3, 1, 1), (262, 3, 1, 1), (273, 4, 1, None),
            (277, 3, 1, 1), (278, 4, 1, nrows), (279, 4, 1, arr.nbytes),
            (339, 3, 1, 3), (33550, 12, 3, scale), (33922, 12, 6, tie),
            (34735, 3, len(geokeys) // 2, geokeys), (42113, 2, len(nod), nod)]

    ifd_off = 8
    blob_off = ifd_off + 2 + 12 * len(tags) + 4
    blobs, entries, cur = [], [], blob_off
    for tag, typ, cnt, val in tags:
        if isinstance(val, bytes) and len(val) > 4:
            entries.append((tag, typ, cnt, struct.pack("<I", cur)))
            blobs.append(val)
            cur += len(val) + (len(val) & 1)
        elif isinstance(val, bytes):
            entries.append((tag, typ, cnt, val.ljust(4, b"\0")))
        else:
            entries.append((tag, typ, cnt, val))
    data_off = cur + (-cur % 4)
    entries = [(t, ty, c, struct.pack("<I", data_off) if (t == 273 and v is None) else v)
               for t, ty, c, v in entries]

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(b"II" + struct.pack("<HI", 42, ifd_off))
        fh.write(struct.pack("<H", len(entries)))
        for tag, typ, cnt, val in sorted(entries):
            packed = (val if isinstance(val, bytes) else
                      (struct.pack("<H", val) + b"\0\0" if typ == 3
                       else struct.pack("<I", val)))
            fh.write(struct.pack("<HHI", tag, typ, cnt) + packed)
        fh.write(struct.pack("<I", 0))
        for b in blobs:
            fh.write(b + (b"\0" if len(b) & 1 else b""))
        fh.write(b"\0" * (data_off - fh.tell()))
        fh.write(arr.tobytes())
    return path


def plot_locations():
    """Plot centres from the inventory workbook (falls back to the CSV)."""
    if os.path.exists(workbook):
        raw = pd.read_excel(workbook, sheet_name="Tree", header=None)
        hdr = [c for c in range(raw.shape[1]) if raw.iloc[0, c] == "Plot"]
        col = hdr[-1]                       # the plot/coordinate block, not col A
        df = raw.iloc[2:, [col, col + 1, col + 2]].copy()
        df.columns = ["Plot", "Northing", "Easting"]
        df = df.apply(pd.to_numeric, errors="coerce").dropna()
        return df[df.Plot > 0].astype({"Plot": int}).reset_index(drop=True)
    csv = os.path.join(base_dir, "Demeritt_Plot_Location.csv")
    return pd.read_csv(csv).rename(columns={"Plot": "Plot"})


def make_chm(extent=None, out_path=chm_file, log=print):
    """Build and write the plot-area CHM. Returns (path, chm, origin, plots)."""
    plots = plot_locations()
    if extent is None:
        extent = (plots.Easting.min() - chm_buffer, plots.Northing.min() - chm_buffer,
                  plots.Easting.max() + chm_buffer, plots.Northing.max() + chm_buffer)
    log(f"{len(plots)} plots ({plots.Plot.min()}-{plots.Plot.max()}), "
        f"buffer {chm_buffer:.0f} m")
    files = sorted(os.path.join(ALS_dir, f) for f in os.listdir(ALS_dir)
                   if f.lower().endswith(".las"))
    chm, origin = build_chm(files, extent, log=log)
    write_geotiff(out_path, chm, origin)
    log(f"wrote {out_path} ({os.path.getsize(out_path) / 1e6:.1f} MB)")
    return out_path, chm, origin, plots


#==============================================================================
# FIGURE STYLE
#
# Times New Roman where it exists (Windows); the fallbacks are metric-compatible
# clones, so a figure rendered on Linux and re-rendered on Windows lays out
# identically -- only the glyphs change. All text is black.
#==============================================================================
FIG_SERIF = ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"]
FIG_FONTSIZE = 11               # one size for every piece of text in a figure
INK = "#000000"                 # every label, tick and annotation
PLOT_FACE = "#87CEFA"           # light sky blue -- plot symbols
BG_NODATA = "#eceff1"           # cells with no ALS return


def _fig_style():
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif", "font.serif": FIG_SERIF,
        "font.size": FIG_FONTSIZE, "axes.labelsize": FIG_FONTSIZE,
        "axes.titlesize": FIG_FONTSIZE, "xtick.labelsize": FIG_FONTSIZE,
        "ytick.labelsize": FIG_FONTSIZE, "legend.fontsize": FIG_FONTSIZE,
        "figure.dpi": 500, "mathtext.fontset": "stix",
        "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK,
        "axes.edgecolor": INK,
    })


#==============================================================================
# FIGURE 1 -- inventory plots over the canopy height model
#==============================================================================
def read_geotiff(path):
    """Read back a GeoTIFF written by write_geotiff. Returns (array, extent)."""
    b = open(path, "rb").read()
    if b[:2] != b"II":
        raise ValueError("expected a little-endian TIFF")
    ifd = struct.unpack("<I", b[4:8])[0]
    tags = {}
    for k in range(struct.unpack("<H", b[ifd:ifd + 2])[0]):
        o = ifd + 2 + 12 * k
        tag, typ, cnt = struct.unpack("<HHI", b[o:o + 8])
        raw = b[o + 8:o + 12]
        size = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 12: 8}.get(typ, 1) * cnt
        if size > 4:
            off = struct.unpack("<I", raw)[0]
            raw = b[off:off + size]
        tags[tag] = raw
    w = struct.unpack("<I", tags[256][:4])[0]
    h = struct.unpack("<I", tags[257][:4])[0]
    off = struct.unpack("<I", tags[273][:4])[0]
    nod = float(tags[42113].split(b"\0")[0]) if 42113 in tags else np.nan
    sx, sy, _ = struct.unpack("<3d", tags[33550])
    tie = struct.unpack("<6d", tags[33922])
    a = np.frombuffer(b, "<f4", count=w * h, offset=off).reshape(h, w)
    return np.where(a == nod, np.nan, a), (tie[3], tie[3] + w * sx,
                                           tie[4] - h * sy, tie[4])


def make_figure1(tif=chm_file, out=None, log=print):
    """Plot centres over the CHM; written to ../2.Final Paper/1.Figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D

    if out is None:
        out = os.path.abspath(os.path.join(base_dir, "..", "2.Final Paper",
                                           "1.Figure", "Fig1_CHM_plots.png"))
    _fig_style()
    MUTED = INK                       # all text, ticks and frame in black
    FACE, EDGE, BG = PLOT_FACE, INK, BG_NODATA

    chm, (xmin, xmax, ymin, ymax) = read_geotiff(tif)
    p = plot_locations()
    W, H = xmax - xmin, ymax - ymin
    vmax = float(np.nanpercentile(chm, 99.5))

    # Keep only plots the ALS actually covers, decided from the raster itself
    # rather than a hard-coded list, so it stays right if the extent changes.
    res = W / chm.shape[1]
    keep = []
    for _, r in p.iterrows():
        c = int((r.Easting - xmin) / res)
        rw = int((ymax - r.Northing) / res)
        rad = int(round(plot_radius / res))
        win = chm[max(0, rw - rad):rw + rad + 1, max(0, c - rad):c + rad + 1]
        keep.append(win.size > 0 and bool(np.isfinite(win).any()))
    dropped = p.loc[~np.array(keep), "Plot"].tolist()
    p = p[np.array(keep)]
    if dropped:
        log(f"  {len(dropped)} plots outside ALS coverage, not plotted: {dropped}")

    # Cool-to-warm. viridis is monotonic in lightness and colour-vision safe,
    # unlike a blue-white-red ramp (which would also imply a meaningful midpoint
    # that canopy height does not have).
    cmap = ListedColormap(plt.get_cmap("viridis")(np.linspace(0.0, 1.0, 256)))
    cmap.set_bad(BG)

    fig = plt.figure(figsize=(5.6, 7.0))
    ax = fig.add_axes([0.115, 0.075, 0.75, 0.885])
    ax.set_facecolor(BG)
    im = ax.imshow(np.ma.masked_invalid(chm), cmap=cmap, vmin=0, vmax=vmax,
                   extent=[xmin, xmax, ymin, ymax], origin="upper",
                   interpolation="nearest", rasterized=True)
    ax.scatter(p.Easting, p.Northing, s=26, marker="o", facecolor=FACE,
               edgecolor=EDGE, linewidth=0.8, zorder=5)
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    ax.set_xlabel("Easting (m)", color=INK)
    ax.set_ylabel("Northing (m)", color=INK)
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.set_xticks(np.arange(np.ceil(xmin / 1000) * 1000, xmax, 1000))
    ax.set_yticks(np.arange(np.ceil(ymin / 1000) * 1000, ymax, 1000))
    ax.tick_params(colors=MUTED, labelsize=FIG_FONTSIZE,
                   direction="out", length=3)
    for lab in ax.get_yticklabels():                 # read along the axis
        lab.set_rotation(90)
        lab.set_va("center")
    for s in ax.spines.values():
        s.set_color(MUTED); s.set_linewidth(0.6)

    ax.legend(handles=[Line2D([], [], marker="o", linestyle="none",
                              markersize=6.5, markerfacecolor=FACE,
                              markeredgecolor=EDGE, markeredgewidth=0.8,
                              label=f"Inventory plot (n = {len(p)})")],
              loc="lower right", frameon=True, fontsize=FIG_FONTSIZE,
              handletextpad=0.5,
              borderpad=0.5).get_frame().set(facecolor="white", edgecolor=MUTED,
                                             linewidth=0.5, alpha=0.9)

    cb_h = 0.42                                       # short bar, vertically centred
    cax = fig.add_axes([0.885, 0.075 + (0.885 - cb_h) / 2, 0.026, cb_h])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Canopy height (m)", color=INK, fontsize=FIG_FONTSIZE)
    cb.ax.tick_params(colors=MUTED, labelsize=FIG_FONTSIZE, length=3)
    cb.outline.set(edgecolor=MUTED, linewidth=0.6)

    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=500, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log(f"wrote {out}")
    return out


#==============================================================================
# FIGURE S1 -- how the pseudo-plots are placed inside a parent plot
#
# Merged in from Pseudo_Plot_Illustration.py. Pure geometry, no ALS data: it
# calls the same even_offsets_* functions the dataset build uses, so the panel
# count follows subplots_per_parent automatically.
#==============================================================================
def make_figure_s1(out=None, n_sub=None, log=print):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Rectangle
    from matplotlib.lines import Line2D

    n_sub = subplots_per_parent if n_sub is None else n_sub
    if out is None:
        out = os.path.abspath(os.path.join(
            base_dir, "..", "2.Final Paper", "1.Figure",
            "FigS1_pseudo_plot_generation.png"))
    _fig_style()
    circle_color, square_color, parent_color = INK, "#4C9A2A", INK

    fig, axes = plt.subplots(2, 3, figsize=(9.0, 6.4))
    labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
    lim = plot_radius + 1.0

    for i, area in enumerate(plot_size):
        ax = axes[i // 3, i % 3]
        # parent boundary, heavier so it stays distinct from the pseudo-plots
        ax.add_patch(Circle((0, 0), plot_radius, fill=False,
                            edgecolor=parent_color, lw=2.2, zorder=5))

        sub_r = np.sqrt(area / np.pi)
        cc = np.array(even_offsets_circle(plot_radius, sub_r, n_sub))
        for dx, dy in cc:
            ax.add_patch(Circle((dx, dy), sub_r, facecolor="none",
                                edgecolor=circle_color, lw=0.8, zorder=2))
        ax.scatter(cc[:, 0], cc[:, 1], s=7, color=circle_color, zorder=6)

        hs = np.sqrt(area) / 2.0
        sc = np.array(even_offsets_square(plot_radius, hs, n_sub))
        for dx, dy in sc:
            ax.add_patch(Rectangle((dx - hs, dy - hs), 2 * hs, 2 * hs,
                                   facecolor="none", edgecolor=square_color,
                                   lw=0.9, ls=(0, (4, 2)), zorder=3))
        ax.scatter(sc[:, 0], sc[:, 1], s=7, marker="s", color=square_color,
                   zorder=6)

        ax.scatter([0], [0], s=22, color="black", marker="+", lw=0.9, zorder=7)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.5)
        ax.text(0.03, 0.975, f"{labels[i]} {area} m$^2$", transform=ax.transAxes,
                va="top", ha="left", fontsize=FIG_FONTSIZE)

    fig.legend(handles=[
        Line2D([0], [0], color=parent_color, lw=2.2,
               label="Parent plot (500 m$^2$)"),
        Line2D([0], [0], color=circle_color, lw=0.9,
               label=f"Circular pseudo-plot (n = {n_sub})"),
        Line2D([0], [0], color=square_color, lw=1.2, ls=(0, (4, 2)),
               label=f"Square pseudo-plot (n = {n_sub})")],
        loc="lower center", ncol=3, frameon=False, fontsize=FIG_FONTSIZE,
        bbox_to_anchor=(0.5, 0.005))

    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.subplots_adjust(wspace=0.06, hspace=0.10)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=500, bbox_inches="tight")
    plt.close(fig)
    log(f"wrote {out}")
    return out


if {"--chm", "--fig", "--figs1"} & set(sys.argv):
    if "--chm" in sys.argv:
        rest = [a for a in sys.argv[sys.argv.index("--chm") + 1:]
                if not a.startswith("-")]
        ext = tuple(float(v) for v in rest[:4]) if len(rest) >= 4 else None
        make_chm(extent=ext)
    if "--chm" in sys.argv or "--fig" in sys.argv:
        make_figure1()
    if "--figs1" in sys.argv:
        make_figure_s1()
    sys.exit(0)

#==============================================================================
# PSEUDO-PLOT DATASET
#
# Inventory and plot centres come from the workbook; the standalone
# Demeritt_Inventory.csv / Demeritt_Plot_Location.csv are no longer used.
#==============================================================================
print("Loading inventory ...")
raw = pd.read_excel(workbook, sheet_name="Tree", header=None)
head = [raw.iloc[0, c] for c in range(raw.shape[1])]
tcols = {name: c for c, name in enumerate(head[:16]) if isinstance(name, str)}
inv = raw.iloc[2:, [tcols["Plot"], tcols["Azimuth"], tcols["Distance"],
                    tcols["Species"], tcols["DBH"]]
                   + [head.index(v) for v in allometry.values()]].copy()
inv.columns = ["Plot", "Azimuth", "Distance", "Species", "DBH"] + list(allometry)
for c in ["Plot", "Azimuth", "Distance", "DBH"] + list(allometry):
    inv[c] = pd.to_numeric(inv[c], errors="coerce")
inv = inv.dropna(subset=["Plot", "Azimuth", "Distance", "DBH"]).copy()
inv["Plot"] = inv["Plot"].astype(int)
az = np.deg2rad(inv["Azimuth"].values.astype(float))
dist = inv["Distance"].values.astype(float)
inv["dx"] = dist * np.sin(az)   # local X offset from plot centre
inv["dy"] = dist * np.cos(az)   # local Y offset from plot centre

#------------------------------------------------------------------------------
# Carbon -> oven-dry biomass.  The workbook stores carbon: Jenkins, Young and
# Chojnacky are the published biomass equations divided by 2, so biomass is
# twice the stored value; method W applies the species-specific carbon fractions
# of Westfall et al. (2024), so its biomass is the stored value divided by that
# fraction.  Reporting is in oven-dry biomass, so the conversion happens here,
# once, per tree.
CARBON_FRACTION_JYC = 0.5
for key in allometry:
    if key != "W":
        inv[key] = inv[key] / CARBON_FRACTION_JYC
if "W" in allometry:
    cf = inv["Species"].map(lambda s: W_CARBON_FRACTION.get(s, np.nan))
    missing = sorted(set(inv.loc[cf.isna(), "Species"].dropna()))
    if missing:
        raise KeyError(f"no method-W carbon fraction for species {missing}")
    inv["W"] = inv["W"] / cf
print(f"  converted to oven-dry biomass (J/Y/C x{1/CARBON_FRACTION_JYC:.0f}, "
      f"W / species carbon fraction {cf.min():.4f}-{cf.max():.4f})")

locs = plot_locations().rename(columns={"Northing": "y", "Easting": "x"})
print(f"  {len(inv)} trees, {locs['Plot'].nunique()} plot locations, "
      f"{len(allometry)} allometric methods")

#------------------------------------------------------------------------------
# Load ALS per parent plot.  Tiles are read once and the clipped, height
# normalised cloud for each plot is cached, so re-running with different
# pseudo-plot settings does not re-read the whole ~22 GB acquisition.
os.makedirs(results_dir, exist_ok=True)
cache_file = os.path.join(results_dir, "plot_clouds.npz")
plot_coords = {int(r["Plot"]): (r["x"], r["y"]) for _, r in locs.iterrows()}
buffer_m = plot_radius + 1.0

if os.path.exists(cache_file):
    print("Loading cached plot clouds ...")
    z = np.load(cache_file)
    plot_clouds = {int(k): z[k] for k in z.files}
else:
    print("Loading ALS per parent plot ...")
    las_files = sorted(os.path.join(ALS_dir, f)
                       for f in os.listdir(ALS_dir) if f.lower().endswith(".las"))
    raw_pts = {p: [] for p in plot_coords}
    raw_cls = {p: [] for p in plot_coords}
    for las_path in las_files:
        h = las_header(las_path)
        relevant = [p for p, (cx, cy) in plot_coords.items()
                    if cx + buffer_m >= h["xmin"] and cx - buffer_m <= h["xmax"]
                    and cy + buffer_m >= h["ymin"] and cy - buffer_m <= h["ymax"]]
        if not relevant:
            continue
        print(f"  {os.path.basename(las_path)} ...")
        for xs, ys, zs, cls, wh in las_chunks(las_path):
            for pid in relevant:
                cx, cy = plot_coords[pid]
                mask = ((xs >= cx - buffer_m) & (xs <= cx + buffer_m) &
                        (ys >= cy - buffer_m) & (ys <= cy + buffer_m))
                if not mask.any():
                    continue
                dx = xs[mask] - cx
                dy = ys[mask] - cy
                circ = dx**2 + dy**2 <= buffer_m**2
                if circ.any():
                    raw_pts[pid].append(
                        np.column_stack([dx[circ], dy[circ], zs[mask][circ]]))
                    raw_cls[pid].append(cls[mask][circ])

    # Normalise heights against the per-plot ground elevation
    print("  Normalising heights ...")
    plot_clouds = {}
    for pid in plot_coords:
        if not raw_pts[pid]:
            continue
        pts = np.vstack(raw_pts[pid])
        cls_arr = np.concatenate(raw_cls[pid])
        raw_pts[pid] = []
        raw_cls[pid] = []
        gnd = cls_arr == GROUND_CLASS
        if gnd.sum() >= 3:
            ground_z = float(np.median(pts[gnd, 2]))
        else:
            ground_z = float(np.percentile(pts[:, 2], 1))
        valid = (cls_arr != NOISE_CLASSES[0]) & (cls_arr != NOISE_CLASSES[1])
        vp = pts[valid]
        h = vp[:, 2] - ground_z
        pos = (h >= -0.05) & (h <= 35.0)
        plot_clouds[pid] = np.column_stack([vp[pos, 0], vp[pos, 1], h[pos]]
                                           ).astype(np.float32)
    del raw_pts, raw_cls
    gc.collect()
    np.savez_compressed(cache_file, **{str(k): v for k, v in plot_clouds.items()})
    print(f"  cached -> {cache_file}")

all_pids = sorted(plot_clouds.keys())
missing = [p for p in sorted(plot_coords) if p not in plot_clouds]
print(f"  {len(all_pids)} plots with ALS returns"
      + (f"; NO coverage for plots {missing}" if missing else ""))

#------------------------------------------------------------------------------
# Generate pseudo-plots
print("Generating pseudo-plots ...")

records = []
for shape in ["circle", "square"]:
    for area in plot_size:
        prefix = "c" if shape == "circle" else "s"
        if shape == "circle":
            radius = np.sqrt(area / np.pi)
            offsets = even_offsets_circle(plot_radius, radius, subplots_per_parent)
        else:
            half_side = np.sqrt(area) / 2.0
            offsets = even_offsets_square(plot_radius, half_side, subplots_per_parent)
        for idx in range(len(all_pids) * subplots_per_parent):
            pid = all_pids[idx // subplots_per_parent]
            dx, dy = offsets[idx % subplots_per_parent]
            records.append({
                "Plot_ID": f"{prefix}_{area}_{idx + 1:03d}",
                "Shape": shape, "Plot_size_m2": area,
                "parent_plot": pid, "offset_dx": dx, "offset_dy": dy,
                "Easting": plot_coords[pid][0] + dx,
                "Northing": plot_coords[pid][1] + dy,
            })
subplots = pd.DataFrame(records)
print(f"  {len(subplots)} pseudo-plots generated "
      f"({len(all_pids)} parents x {subplots_per_parent} each "
      f"x {len(plot_size)} sizes x 2 shapes)")

#------------------------------------------------------------------------------
# Assign trees and compute AGB (oven-dry Mg ha-1)
print("Computing plot-level AGB ...")
agb_cols = [f"AGB(Mg/ha)_{k}" for k in allometry] + ["AGB(Mg/ha)_Mean"]
KG_M2_TO_MG_HA = 10.0            # 1 kg m-2 = 10 Mg ha-1
by_plot = {p: g for p, g in inv.groupby("Plot")}
agb_results = []
for _, sp in subplots.iterrows():
    trees = by_plot.get(sp["parent_plot"])
    if trees is None or trees.empty:
        agb_results.append({k: 0.0 for k in agb_cols})
        continue
    tx = trees["dx"].values - sp["offset_dx"]
    ty = trees["dy"].values - sp["offset_dy"]
    area = sp["Plot_size_m2"]
    if sp["Shape"] == "circle":
        r = np.sqrt(area / np.pi)
        inside = tx**2 + ty**2 <= r**2
    else:
        hs = np.sqrt(area) / 2.0
        inside = (np.abs(tx) <= hs) & (np.abs(ty) <= hs)
    sel = trees[inside]
    row, vals = {}, []
    for key in allometry:
        # tree values are oven-dry biomass (kg); per unit area, then to Mg ha-1
        agb = (sel[key].sum() / area) * KG_M2_TO_MG_HA if not sel.empty else 0.0
        row[f"AGB(Mg/ha)_{key}"] = agb
        vals.append(agb)
    row["AGB(Mg/ha)_Mean"] = float(np.mean(vals))
    agb_results.append(row)
subplots = pd.concat([subplots.reset_index(drop=True),
                      pd.DataFrame(agb_results)], axis=1)

#------------------------------------------------------------------------------
# Extract ALS metrics per pseudo-plot
print("Extracting ALS metrics ...")
metric_names = ["Point_density", "Vegetation_density", "Max_height",
                "Mean_height", "P25", "P50", "P75", "P95"]
metrics_list = []
for _, sp in subplots.iterrows():
    cloud = plot_clouds.get(sp["parent_plot"])
    if cloud is None:
        metrics_list.append({m: np.nan for m in metric_names})
        continue
    cx, cy = sp["offset_dx"], sp["offset_dy"]
    area = sp["Plot_size_m2"]
    if sp["Shape"] == "circle":
        r = np.sqrt(area / np.pi)
        mask = (cloud[:, 0] - cx)**2 + (cloud[:, 1] - cy)**2 <= r**2
    else:
        hs = np.sqrt(area) / 2.0
        mask = (np.abs(cloud[:, 0] - cx) <= hs) & (np.abs(cloud[:, 1] - cy) <= hs)
    pts = cloud[mask]
    if len(pts) < 5:
        metrics_list.append({m: np.nan for m in metric_names})
        continue
    h = pts[:, 2]
    total = len(h)
    # Vegetation density: proportion of returns above HEIGHT_THRESH
    # relative to all returns in the plot (consistent across all plots)
    veg = h[h >= height_thresh]
    veg_density = len(veg) / total if total > 0 else 0.0
    if len(veg) < 3:
        metrics_list.append({m: np.nan for m in metric_names})
        continue
    metrics_list.append({
        "Point_density": total / area,
        "Vegetation_density": veg_density,
        "Max_height": np.max(veg),
        "Mean_height": np.mean(veg),
        "P25": np.percentile(veg, 25),
        "P50": np.percentile(veg, 50),
        "P75": np.percentile(veg, 75),
        "P95": np.percentile(veg, 95),
    })
subplots = pd.concat([subplots.reset_index(drop=True),
                      pd.DataFrame(metrics_list)], axis=1)

#------------------------------------------------------------------------------
# Save final dataset
out_cols = (["Plot_ID", "Shape", "Plot_size_m2", "parent_plot",
             "Northing", "Easting"] + agb_cols + metric_names)
final = subplots[out_cols].dropna()
final.to_csv(save_file, index=False)
print(f"\nSaved {len(final)} rows to {save_file}")
print(f"  Shapes: {final['Shape'].value_counts().to_dict()}")
print(f"  Sizes:  {sorted(final['Plot_size_m2'].unique())}")
print(f"  Per shape-size: "
      f"{sorted(final.groupby(['Shape','Plot_size_m2']).size().unique())}")
print("Done.")

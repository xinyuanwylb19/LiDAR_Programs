"""
RF_Model_Analysis.py
Three-objective RF biomass model analysis:
  1. Allometric method comparison (300 m²; 120 AGB-evenly distributed
     calibration plots + 30 random independent validation plots).
  2. Plot size & shape effect (120 AGB-evenly distributed calibration plots
     per shape × size; common validation = 30 random square 225 m² plots
     and 30 random square 100 m² plots, excluding any calibration overlap).
  3. Sample size & variance effect (plots 20–160 step 20; 150 candidate
     calibration sets per n_cal sorted into 5 variance quintile groups of 30,
     plot size 300 m²).

Author: Xinyuan Wei
"""

import os
# Silence all warnings: must be done BEFORE importing sklearn so child
# processes spawned by joblib (n_jobs > 1) inherit the setting.
os.environ.setdefault("PYTHONWARNINGS", "ignore")
import warnings
warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import ttest_rel
warnings.filterwarnings("ignore")

#------------------------------------------------------------------------------
# Paths & parameters
data_dir = os.path.dirname(os.path.abspath(__file__))
res_dir = os.path.join(data_dir, "Results")
# Figures now live in the manuscript folder: ../2.Final Paper/1.Figure
fig_dir = os.path.abspath(os.path.join(data_dir, "..", "2.Final Paper", "1.Figure"))
os.makedirs(res_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

seed = 42
n_trees = 100
min_leaf = 3
n_iter = 30                         # iterations per variance group (Objective 3)
predictors = ["Point_density", "Vegetation_density", "Max_height",
              "Mean_height", "P25", "P50", "P75", "P95"]
allo_keys = ["Y", "J", "C", "W"]
# Four allometric methods, keyed as in the inventory workbook (J/Y/C/W)
allo_display = {"Y": "Young et al. 1980", "J": "Jenkins et al. 2003",
                "C": "Chojnacky et al. 2014", "W": "Westfall et al. 2024"}

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"],
    "font.size": 9, "axes.labelsize": 10,
    "axes.titlesize": 11, "figure.dpi": 500,
    "mathtext.fontset": "stix",
})

#------------------------------------------------------------------------------
# Load prepared dataset
print("Loading prepared dataset ...")
df = pd.read_csv(os.path.join(data_dir, "Prepared_Dataset.csv"))
print(f"  {len(df)} plots, shapes={df['Shape'].unique().tolist()}, "
      f"sizes={sorted(df['Plot_size_m2'].unique())}, "
      f"parents={df['parent_plot'].nunique()}")

#------------------------------------------------------------------------------
# Helper functions
def fit_rf(x, y, rs=seed):
    rf = RandomForestRegressor(n_estimators=n_trees, min_samples_leaf=min_leaf,
                               max_features="sqrt", random_state=rs, n_jobs=-1)
    rf.fit(x, y); return rf

def cv_predict(x, y, n_folds=5, rs=seed):
    """Return OOF predictions via k-fold CV (per-fold seeding)."""
    kf = KFold(n_splits=min(n_folds, len(y)), shuffle=True, random_state=rs)
    pred = np.full(len(y), np.nan)
    for k, (tr, te) in enumerate(kf.split(x)):
        rf = fit_rf(x[tr], y[tr], rs=rs + k); pred[te] = rf.predict(x[te])
    return pred

def cal_metrics(y, pred):
    # Negative R² preserved (the actual sklearn r2_score)
    return r2_score(y, pred), np.sqrt(mean_squared_error(y, pred))

def val_metrics(y, pred):
    rmse = np.sqrt(mean_squared_error(y, pred))
    bias = float(np.mean(pred - y))
    _, p = ttest_rel(pred, y)
    return rmse, bias, float(p)

def stratified_split(data, y_col, n_cal=60, rs=seed):
    """Stratified split: 3 quantile bins, ~n_cal/3 from each, no input mutation."""
    rng = np.random.default_rng(rs)
    bins = pd.qcut(data[y_col], q=3, labels=False, duplicates="drop")
    n_bins = bins.nunique()
    per_bin = n_cal // n_bins
    cal_idx = []
    for b in range(n_bins):
        qidx = np.array(data.index[bins == b].to_numpy(), copy=True)
        rng.shuffle(qidx)
        cal_idx.extend(qidx[:per_bin].tolist())
    if len(cal_idx) < n_cal:
        remaining = np.array(data.index.difference(cal_idx).to_numpy(), copy=True)
        rng.shuffle(remaining)
        cal_idx.extend(remaining[:n_cal - len(cal_idx)].tolist())
    val_idx = data.index.difference(cal_idx)
    return data.loc[cal_idx], data.loc[val_idx]

def even_agb_split(data, y_col, n_cal=120):
    """Pick n_cal calibration plots whose AGB values are evenly distributed
    from the minimum to the maximum (selected at equally-spaced positions in
    AGB-sorted order). Remaining plots are returned as the residual pool."""
    sorted_df = data.sort_values(y_col).reset_index()  # 'index' = original idx
    N = len(sorted_df)
    if n_cal >= N: return data.copy(), data.iloc[0:0].copy()
    pick_pos = np.unique(np.linspace(0, N - 1, n_cal).round().astype(int))
    if len(pick_pos) < n_cal:
        remaining_pos = np.setdiff1d(np.arange(N), pick_pos)
        pick_pos = np.concatenate([pick_pos, remaining_pos[:n_cal - len(pick_pos)]])
    cal_orig_idx = sorted_df.iloc[pick_pos]["index"].values
    pool_orig_idx = data.index.difference(cal_orig_idx)
    return data.loc[cal_orig_idx].copy(), data.loc[pool_orig_idx].copy()

def random_subset(data, n, rs=seed):
    """Random subset of size n from `data`, no replacement."""
    rng = np.random.default_rng(rs)
    if n >= len(data): return data.copy()
    idx = rng.choice(data.index.to_numpy(), size=n, replace=False)
    return data.loc[idx].copy()

def make_scatter_grid(data_dict, keys_row, keys_col, col_labels,
                      annotation_fn, fig_title, out_name, row_labels,
                      xlabel="Observed AGB (Mg ha$^{-1}$)",
                      ylabel="Predicted AGB (Mg ha$^{-1}$)",
                      font_scale=1.0, lim=560, tick=None):
    step = tick if tick else lim / 4
    _hidden = 0
    nrows = len(keys_row); ncols = len(keys_col)
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(2.1*ncols + 0.7, 2.1*nrows + 0.65),
                              sharex=True, sharey=True, squeeze=False)
    shape_marker = {"circle": "o", "square": "s"}
    for r, rk in enumerate(keys_row):
        for c, ck in enumerate(keys_col):
            ax = axes[r, c]
            obs, pred = data_dict[(rk, ck)]
            _hidden += int(np.sum((np.asarray(obs) > lim) | (np.asarray(pred) > lim)))
            ax.scatter(obs, pred, s=20, alpha=0.3, color="blue",
                       edgecolor="black", linewidth=0.3,
                       marker=shape_marker.get(rk, "o"))
            ax.plot([0, lim], [0, lim], "r--", lw=1)
            txt = annotation_fn(obs, pred)
            ax.text(0.05, 0.95, txt, transform=ax.transAxes, fontsize=8*font_scale,
                    va="top", bbox=dict(boxstyle="round,pad=0.3",
                                        fc="none", ec="none"))
            ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
            ax.set_xticks(np.arange(0, lim + 0.1, step))
            ax.set_yticks(np.arange(0, lim + 0.1, step))
            if r == 0: ax.set_title(col_labels[c], fontsize=10*font_scale)
    for r in range(nrows):
        for c in range(ncols):
            if c > 0: axes[r, c].tick_params(labelleft=False)
            if r < nrows - 1: axes[r, c].tick_params(labelbottom=False)
    # Panels sit shoulder to shoulder, so the trailing "16" of one panel and the
    # leading "0" of the next print as "160". Blank the leading label on panels 2+.
    for c in range(1, ncols):
        lab = [f"{v:g}" for v in np.arange(0, lim + 0.1, step)]
        lab[0] = ""
        axes[nrows - 1, c].set_xticklabels(lab)
    for r in range(nrows): axes[r, 0].set_ylabel(ylabel)
    for c in range(ncols): axes[nrows-1, c].set_xlabel(xlabel)
    fig.tight_layout(rect=[0.07, 0, 1, 1])
    fig.subplots_adjust(wspace=0.05, hspace=0.08)
    for r, rl in enumerate(row_labels):
        axes[r, 0].annotate(rl, xy=(-0.34, 0.5), xycoords="axes fraction",
                            fontsize=10*font_scale, rotation=90, ha="center", va="center")
    if _hidden:
        print(f"    NOTE {out_name}: {_hidden} point(s) fall outside the 0-{lim:g} axis")
    fig.savefig(os.path.join(fig_dir, out_name), dpi=500, bbox_inches="tight")
    plt.close(fig)

def ann_cal(obs, pred):
    r2, rmse = cal_metrics(obs, pred)
    return f"R² = {r2:.3f}\nRMSE = {rmse:.2f}"

def ann_val(obs, pred):
    vrmse, vbias, vp = val_metrics(obs, pred)
    ptxt = "p-value < 0.001" if vp < 0.001 else f"p-value = {vp:.3f}"
    return f"RMSE = {vrmse:.2f}\nbias = {vbias:.2f}\n{ptxt}"

def ann_val_sub(obs, pred, sub):
    """Validation annotation with a subscript label (e.g. '225' or '100')."""
    vrmse, vbias, vp = val_metrics(obs, pred)
    ptxt = (f"p-value$_{{{sub}}}$ < 0.001" if vp < 0.001
            else f"p-value$_{{{sub}}}$ = {vp:.3f}")
    return (f"RMSE$_{{{sub}}}$ = {vrmse:.2f}\n"
            f"bias$_{{{sub}}}$ = {vbias:.2f}\n{ptxt}")

def make_scatter_grid_two(data_dict, keys_row, keys_col, col_labels,
                          tl_sub, br_sub, fig_title, out_name, row_labels,
                          tl_color="#1f77b4", br_color="#d95f02",
                          xlabel="Observed AGB (Mg ha$^{-1}$)",
                          ylabel="Predicted AGB (Mg ha$^{-1}$)",
                          font_scale=1.0, lim=560, tick=None):
    """data_dict[(rk, ck)] = ((obs_tl, pred_tl), (obs_br, pred_br))."""
    step = tick if tick else lim / 4
    _hidden = 0
    nrows = len(keys_row); ncols = len(keys_col)
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(2.1*ncols + 0.7, 2.1*nrows + 0.65),
                              sharex=True, sharey=True, squeeze=False)
    shape_marker = {"circle": "o", "square": "s"}
    for r, rk in enumerate(keys_row):
        for c, ck in enumerate(keys_col):
            ax = axes[r, c]
            (obs_tl, pred_tl), (obs_br, pred_br) = data_dict[(rk, ck)]
            _hidden += int(np.sum((np.asarray(obs_tl) > lim) | (np.asarray(pred_tl) > lim)))
            _hidden += int(np.sum((np.asarray(obs_br) > lim) | (np.asarray(pred_br) > lim)))
            ax.scatter(obs_tl, pred_tl, s=20, alpha=0.35, color=tl_color,
                       edgecolor="black", linewidth=0.3,
                       marker=shape_marker.get(rk, "o"), label=f"{tl_sub} m²")
            ax.scatter(obs_br, pred_br, s=20, alpha=0.35, color=br_color,
                       edgecolor="black", linewidth=0.3,
                       marker=shape_marker.get(rk, "o"), label=f"{br_sub} m²")
            ax.plot([0, lim], [0, lim], "r--", lw=1)
            ax.text(0.05, 0.95, ann_val_sub(obs_tl, pred_tl, tl_sub),
                    transform=ax.transAxes, fontsize=7.5*font_scale,
                    va="top", ha="left", color=tl_color,
                    bbox=dict(boxstyle="round,pad=0.2", fc="none", ec="none"))
            ax.text(0.95, 0.05, ann_val_sub(obs_br, pred_br, br_sub),
                    transform=ax.transAxes, fontsize=7.5*font_scale,
                    va="bottom", ha="right", color=br_color,
                    bbox=dict(boxstyle="round,pad=0.2", fc="none", ec="none"))
            ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
            ax.set_xticks(np.arange(0, lim + 0.1, step))
            ax.set_yticks(np.arange(0, lim + 0.1, step))
            if r == 0: ax.set_title(col_labels[c], fontsize=10*font_scale)
    for r in range(nrows):
        for c in range(ncols):
            if c > 0: axes[r, c].tick_params(labelleft=False)
            if r < nrows - 1: axes[r, c].tick_params(labelbottom=False)
    # Panels sit shoulder to shoulder, so the trailing "16" of one panel and the
    # leading "0" of the next print as "160". Blank the leading label on panels 2+.
    for c in range(1, ncols):
        lab = [f"{v:g}" for v in np.arange(0, lim + 0.1, step)]
        lab[0] = ""
        axes[nrows - 1, c].set_xticklabels(lab)
    for r in range(nrows): axes[r, 0].set_ylabel(ylabel)
    for c in range(ncols): axes[nrows-1, c].set_xlabel(xlabel)
    fig.tight_layout(rect=[0.07, 0, 1, 1])
    fig.subplots_adjust(wspace=0.05, hspace=0.08)
    for r, rl in enumerate(row_labels):
        axes[r, 0].annotate(rl, xy=(-0.34, 0.5), xycoords="axes fraction",
                            fontsize=10*font_scale, rotation=90,
                            ha="center", va="center")
    if _hidden:
        print(f"    NOTE {out_name}: {_hidden} point(s) fall outside the 0-{lim:g} axis")
    fig.savefig(os.path.join(fig_dir, out_name), dpi=500, bbox_inches="tight")
    plt.close(fig)

#------------------------------------------------------------------------------
# OBJECTIVE 1: Allometric method comparison (300 m²; 120 cal + 30 val)
#------------------------------------------------------------------------------
print("\n" + "-" * 70)
print("OBJECTIVE 1 — Allometric method comparison")
print("-" * 70)

n_cal_o1 = 120                                # AGB-evenly distributed plots
n_val_o1 = 30                                 # random independent validation
o1_cal, o1_val = [], []
o1_cal_data, o1_val_data = {}, {}

for shape in ["circle", "square"]:
    d300 = df[(df["Shape"] == shape) & (df["Plot_size_m2"] == 300)].copy()
    for allo in allo_keys:
        y_col = f"AGB(Mg/ha)_{allo}"
        cal, pool = even_agb_split(d300, y_col, n_cal=n_cal_o1)
        val = random_subset(pool, n=n_val_o1, rs=seed)
        x_cal = cal[predictors].values; y_cal = cal[y_col].values
        x_val = val[predictors].values; y_val = val[y_col].values
        # Fit RF once with OOB scoring; OOB predictions give the calibration
        # R² / RMSE while the same fitted model predicts the held-out val plots.
        rf = RandomForestRegressor(n_estimators=n_trees, min_samples_leaf=min_leaf,
                                   max_features="sqrt", random_state=seed,
                                   n_jobs=-1, oob_score=True)
        rf.fit(x_cal, y_cal)
        oob_pred = rf.oob_prediction_
        r2, rmse = cal_metrics(y_cal, oob_pred)
        o1_cal.append({"shape": shape, "allometric": allo_display[allo],
                       "R2": r2, "RMSE": rmse})
        o1_cal_data[(shape, allo)] = (y_cal, oob_pred)
        pred_val = rf.predict(x_val)
        vrmse, vbias, vpval = val_metrics(y_val, pred_val)
        o1_val.append({"shape": shape, "allometric": allo_display[allo],
                       "RMSE": vrmse, "bias": vbias, "p_value": vpval})
        o1_val_data[(shape, allo)] = (y_val, pred_val)
        print(f"  {shape}/{allo_display[allo]}: cal R²={r2:.3f} RMSE={rmse:.2f} | "
              f"val RMSE={vrmse:.2f} bias={vbias:.2f} p={vpval:.3f}")

pd.DataFrame(o1_cal).to_csv(os.path.join(res_dir, "Objective1_Calibration_Results.csv"), index=False)
pd.DataFrame(o1_val).to_csv(os.path.join(res_dir, "Objective1_Validation_Results.csv"), index=False)

# Figure 2 & 3
make_scatter_grid(o1_cal_data, ["circle", "square"], allo_keys,
    [allo_display[a] for a in allo_keys], ann_cal,
    "Figure 2: Experiment 1 — Calibration", "Fig2_O1_calibration.png",
    ["Circle Plots", "Square Plots"], lim=300, tick=100)
print("  Fig2_O1_calibration.png")

make_scatter_grid(o1_val_data, ["circle", "square"], allo_keys,
    [allo_display[a] for a in allo_keys], ann_val,
    "Figure 3: Experiment 1 — Validation", "Fig3_O1_validation.png",
    ["Circle Plots", "Square Plots"], lim=300, tick=100)
print("  Fig3_O1_validation.png")

#------------------------------------------------------------------------------
# OBJECTIVE 2: Plot size & shape effect (120 AGB-even cal per shape × size;
#              independent validation = 30 random square 225 m² and 30 random
#              square 100 m² sub-plots, excluding any calibration overlap).
#------------------------------------------------------------------------------
print("\n" + "-" * 70)
print("OBJECTIVE 2 — Plot size and shape effect")
print("-" * 70)

cal_sizes = [100, 150, 200, 225, 250, 300]
y_col2 = "AGB(Mg/ha)_Mean"
n_cal_o2 = 120                                # AGB-evenly distributed
n_val_o2 = 30                                 # random per validation size

# Calibrate every (shape, size) on its own 120 AGB-evenly distributed plots.
# Track which sub-plot indices were used as calibration so we can exclude
# them when building the validation sets for square 225 and square 100.
cal_used_idx = {}                                   # (shape, area) -> Index
rfs = {}                                            # (shape, area) -> fitted RF

o2_cal = []
o2_cal_data = {}

for shape in ["circle", "square"]:
    for area in cal_sizes:
        sub = df[(df["Shape"] == shape) & (df["Plot_size_m2"] == area)].copy()
        cal, _ = even_agb_split(sub, y_col2, n_cal=n_cal_o2)
        cal_used_idx[(shape, area)] = cal.index
        x = cal[predictors].values; y = cal[y_col2].values
        rf = RandomForestRegressor(n_estimators=n_trees, min_samples_leaf=min_leaf,
                                   max_features="sqrt", random_state=seed,
                                   n_jobs=-1, oob_score=True)
        rf.fit(x, y)
        rfs[(shape, area)] = rf
        oob_pred = rf.oob_prediction_
        r2, rmse = cal_metrics(y, oob_pred)
        o2_cal.append({"shape": shape, "area_m2": area,
                       "R2": r2, "RMSE": rmse, "n_train": len(y)})
        o2_cal_data[(shape, area)] = (y, oob_pred)
        print(f"  {shape}/{area}m²: cal R²={r2:.3f} RMSE={rmse:.2f} (n={len(y)})")

# Build the two independent validation sets (square 225 m² and square 100 m²).
# Each set is 30 randomly drawn plots, excluding plots used as calibration for
# the model trained at the same shape × size.
val_sets = {}                                       # size_m2 -> (x, y)
val_dfs  = {}
for area in (225, 100):
    pool = df[(df["Shape"] == "square") & (df["Plot_size_m2"] == area)].copy()
    excl = cal_used_idx.get(("square", area), pd.Index([]))
    pool = pool.loc[~pool.index.isin(excl)].copy()
    vs = random_subset(pool, n=n_val_o2, rs=seed + area)
    val_dfs[area] = vs
    val_sets[area] = (vs[predictors].values, vs[y_col2].values)
print(f"  Validation sets: 225 m² n={len(val_dfs[225])}, "
      f"100 m² n={len(val_dfs[100])}")

# Apply each (shape, size) model to both validation sets.
o2_app = []
o2_app_data = {}
for shape in ["circle", "square"]:
    for area in cal_sizes:
        rf = rfs[(shape, area)]
        x225, y225 = val_sets[225]; pred225 = rf.predict(x225)
        x100, y100 = val_sets[100]; pred100 = rf.predict(x100)
        r225, b225, p225 = val_metrics(y225, pred225)
        r100, b100, p100 = val_metrics(y100, pred100)
        o2_app.append({"shape": shape, "area_m2": area,
                       "RMSE_225": r225, "bias_225": b225, "p_value_225": p225,
                       "RMSE_100": r100, "bias_100": b100, "p_value_100": p100,
                       "n_val_225": len(y225), "n_val_100": len(y100)})
        o2_app_data[(shape, area)] = ((y225, pred225), (y100, pred100))
        print(f"    {shape}/{area}m² | 225-val RMSE={r225:.2f} bias={b225:.2f} "
              f"p={p225:.3f}  100-val RMSE={r100:.2f} bias={b100:.2f} p={p100:.3f}")

pd.DataFrame(o2_cal).to_csv(os.path.join(res_dir, "Objective2_Calibration_Results.csv"), index=False)
pd.DataFrame(o2_app).to_csv(os.path.join(res_dir, "Objective2_Application_Results.csv"), index=False)

# Figure 4 & 5
make_scatter_grid(o2_cal_data, ["circle", "square"], cal_sizes,
    [f"{a} m²" for a in cal_sizes], ann_cal,
    "Figure 4: Experiment 2 — Calibration (OOB predictions)",
    "Fig4_O2_calibration.png", ["Circle Plots", "Square Plots"],
    font_scale=1.2, lim=400, tick=100)
print("  Fig4_O2_calibration.png")

make_scatter_grid_two(o2_app_data, ["circle", "square"], cal_sizes,
    [f"{a} m²" for a in cal_sizes], tl_sub="225", br_sub="100",
    fig_title="Figure 5: Experiment 2 — Application to 225 m² and 100 m² square plots",
    out_name="Fig5_O2_application.png",
    row_labels=["Circle Plots", "Square Plots"],
    tl_color="#1f77b4", br_color="#d95f02", font_scale=1.2,
    lim=400, tick=100)
print("  Fig5_O2_application.png")

#------------------------------------------------------------------------------
# OBJECTIVE 3: Sample size & biomass variance effect.
# At each sample size n_cal: draw 150 random calibration sets of size n_cal,
# rank them by AGB variance, split into 5 quintile groups (low → high) of 30
# sets each. For each (n_cal, variance) cell, the 30 sets in that group serve
# as the 30 iterations of the analysis. Plot size 300 m²
#------------------------------------------------------------------------------
print("\n" + "-" * 70)
print("OBJECTIVE 3 — Sample size and biomass variance")
print("-" * 70)

sample_sizes = [20, 40, 60, 80, 100, 120, 140, 160]   # capped at 160
var_levels = ["low", "low-medium", "medium", "medium-high", "high"]
n_candidates = 150                                   # candidate cal sets per n_cal
sets_per_level = n_candidates // len(var_levels)     # 150 / 5 = 30
y_col3 = "AGB(Mg/ha)_Mean"
o3_cal_rows, o3_val_rows = [], []
o3_raw = []                                          # per-iteration data for boxplots

# Single rng across all combinations (no per-loop reseeding)
rng = np.random.default_rng(seed)

for shape in ["circle", "square"]:
    pool = df[(df["Shape"] == shape) & (df["Plot_size_m2"] == 300)].copy()
    pool = pool.reset_index(drop=True)
    n_total = len(pool)
    x_pool = pool[predictors].values; y_pool = pool[y_col3].values

    for n_cal in sample_sizes:
        if n_total - n_cal < 5:
            print(f"  {shape}/n={n_cal}: skipped (too few validation plots)")
            continue
        # Draw 150 random candidate calibration sets, rank by AGB variance,
        # and split into 5 quintile groups of 30 sets each (no overlap).
        candidates = []
        for _ in range(n_candidates):
            idx = rng.choice(n_total, size=n_cal, replace=False)
            candidates.append((idx, float(np.var(y_pool[idx]))))
        candidates.sort(key=lambda t: t[1])
        groups = {vl: [c[0] for c in candidates[i*sets_per_level:(i+1)*sets_per_level]]
                  for i, vl in enumerate(var_levels)}

        for var_level in var_levels:
            pool_subs = groups[var_level]
            actual_iter = len(pool_subs)
            if actual_iter < 5:
                print(f"  {shape}/n={n_cal}/{var_level}: skipped")
                continue
            cal_r2s, cal_rmses = [], []
            val_rmses, val_biases, val_ps = [], [], []
            for it in range(actual_iter):
                cal_idx = pool_subs[it]
                val_idx = np.setdiff1d(np.arange(n_total), cal_idx,
                                       assume_unique=True)
                x_c, y_c = x_pool[cal_idx], y_pool[cal_idx]
                x_v, y_v = x_pool[val_idx], y_pool[val_idx]
                # OOB predictions on the calibration set give cal R²/RMSE;
                # the same fitted RF predicts the held-out val plots.
                rf = RandomForestRegressor(n_estimators=n_trees,
                                           min_samples_leaf=min_leaf,
                                           max_features="sqrt",
                                           random_state=seed + it,
                                           n_jobs=-1, oob_score=True)
                rf.fit(x_c, y_c)
                oob_pred = rf.oob_prediction_
                r2, rmse = cal_metrics(y_c, oob_pred)
                cal_r2s.append(r2); cal_rmses.append(rmse)
                pred_v = rf.predict(x_v)
                vrmse, vbias, vp = val_metrics(y_v, pred_v)
                val_rmses.append(vrmse); val_biases.append(vbias); val_ps.append(vp)
                o3_raw.append({"shape": shape, "n_cal": n_cal, "variance": var_level,
                               "cal_R2": r2, "cal_RMSE": rmse,
                               "val_RMSE": vrmse, "val_bias": vbias})

            o3_cal_rows.append({
                "shape": shape, "n_cal": n_cal, "variance": var_level,
                "R2_mean": np.mean(cal_r2s), "R2_std": np.std(cal_r2s),
                "RMSE_mean": np.mean(cal_rmses), "RMSE_std": np.std(cal_rmses),
                "n_iter": actual_iter})
            o3_val_rows.append({
                "shape": shape, "n_cal": n_cal, "variance": var_level,
                "RMSE_mean": np.mean(val_rmses), "RMSE_std": np.std(val_rmses),
                "bias_mean": np.mean(val_biases), "bias_std": np.std(val_biases),
                "p_value_mean": np.mean(val_ps), "p_value_std": np.std(val_ps),
                "n_iter": actual_iter})
            print(f"  {shape}/n={n_cal}/{var_level}: "
                  f"cal R²={np.mean(cal_r2s):.3f}±{np.std(cal_r2s):.3f} "
                  f"RMSE={np.mean(cal_rmses):.2f} | "
                  f"val RMSE={np.mean(val_rmses):.2f} bias={np.mean(val_biases):.2f}")

pd.DataFrame(o3_cal_rows).to_csv(os.path.join(res_dir, "Objective3_Calibration_Results.csv"), index=False)
pd.DataFrame(o3_val_rows).to_csv(os.path.join(res_dir, "Objective3_Validation_Results.csv"), index=False)
o3_df = pd.DataFrame(o3_raw)

# Figures 5 and 6
var_colors = {"low": "#2166AC", "low-medium": "#87CEFA",
              "medium": "#4CAF50", "medium-high": "#FFD700", "high": "#E31A1C"}
var_labels_display = {"low":"Low","low-medium":"Low-medium","medium":"Medium","medium-high":"Medium-high","high":"High"}

box_lw = 0.5                                         # boxplot border / median / whisker / cap
spine_lw = 0.5                                       # axes spine / tick width
box_width = 0.13
offsets = np.array([-0.32, -0.16, 0.0, 0.16, 0.32])
legend_elements = [Patch(facecolor=var_colors[vl], alpha=0.7,
                         label=var_labels_display[vl]) for vl in var_levels]

def styled_boxplot(ax, data, positions, color):
    bp = ax.boxplot(data, positions=positions, widths=box_width,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", lw=box_lw),
                    boxprops=dict(lw=box_lw),
                    whiskerprops=dict(lw=box_lw),
                    capprops=dict(lw=box_lw))
    for patch in bp["boxes"]:
        patch.set_facecolor(color); patch.set_alpha(0.7)

def style_axes(ax):
    for s in ax.spines.values(): s.set_linewidth(spine_lw)
    ax.tick_params(width=spine_lw)

# Figure 6: Calibration boxplots (R² and RMSE) — smaller text, polished layout
print("  Generating Fig6 ...")
sb_var_colors = {"low":"#2166AC","low-medium":"#87CEFA","medium":"#4CAF50",
                  "medium-high":"#FFD700","high":"#E31A1C"}
sb_var_labels = {"low":"Low","low-medium":"Low-medium","medium":"Medium",
                  "medium-high":"Medium-high","high":"High"}
sb_box_lw = 0.5; sb_spine_lw = 0.5; sb_box_width = 0.13
sb_offsets = np.array([-0.32,-0.16,0.0,0.16,0.32])
sb_legend_elements = [Patch(facecolor=sb_var_colors[vl], alpha=0.7,
                             label=sb_var_labels[vl]) for vl in var_levels]

def _styled_boxplot(ax, data, positions, color):
    bp = ax.boxplot(data, positions=positions, widths=sb_box_width,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", lw=sb_box_lw),
                    boxprops=dict(lw=sb_box_lw),
                    whiskerprops=dict(lw=sb_box_lw),
                    capprops=dict(lw=sb_box_lw))
    for patch in bp["boxes"]:
        patch.set_facecolor(color); patch.set_alpha(0.7)
def _style_axes(ax):
    for sp in ax.spines.values(): sp.set_linewidth(sb_spine_lw)
    ax.tick_params(width=sb_spine_lw)

# Font sizes (kept constant; panels shrunk so text appears relatively larger)
fs_tick = 8; fs_label = 9; fs_title = 10; fs_suptitle = 10; fs_legend = 8
with plt.rc_context({"font.size": fs_tick, "axes.labelsize": fs_label,
                      "axes.titlesize": fs_title, "xtick.labelsize": fs_tick,
                      "ytick.labelsize": fs_tick, "legend.fontsize": fs_legend}):
    fig, axes = plt.subplots(2, 2, figsize=(5.6, 3.6))
    sparse_xticklabels = [str(s) if i % 2 == 0 else "" for i, s in enumerate(sample_sizes)]
    for col_idx, shape in enumerate(["circle","square"]):
        sub = o3_df[o3_df["shape"]==shape]
        for vi, vl in enumerate(var_levels):
            vs = sub[sub["variance"]==vl]
            r2_data = [vs[vs["n_cal"]==n]["cal_R2"].values for n in sample_sizes]
            rmse_data = [vs[vs["n_cal"]==n]["cal_RMSE"].values for n in sample_sizes]
            positions = np.arange(len(sample_sizes)) + sb_offsets[vi]
            _styled_boxplot(axes[0,col_idx], r2_data, positions, sb_var_colors[vl])
            _styled_boxplot(axes[1,col_idx], rmse_data, positions, sb_var_colors[vl])
        for r in [0,1]:
            ax = axes[r,col_idx]
            ax.set_xticks(range(len(sample_sizes)))
            ax.set_xticklabels(sparse_xticklabels)
            _style_axes(ax)
        axes[0,col_idx].tick_params(labelbottom=False)
        axes[1,col_idx].set_xlabel("Number of Plots")
        axes[0,col_idx].set_ylim(-0.5, 1); axes[0,col_idx].set_yticks(np.arange(-0.5, 1.01, 0.5))
        axes[1,col_idx].set_ylim(0, 80); axes[1,col_idx].set_yticks([0, 20, 40, 60, 80])
        if col_idx == 0:
            axes[0,col_idx].set_ylabel("R²")
            axes[1,col_idx].set_ylabel("RMSE (Mg ha$^{-1}$)")
        else:
            axes[0,col_idx].set_ylabel(""); axes[0,col_idx].set_yticklabels([])
            axes[1,col_idx].set_ylabel(""); axes[1,col_idx].set_yticklabels([])
        axes[0,col_idx].set_title(f"{shape.title()} Plots", fontsize=fs_title)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.96, bottom=0.22,
                         wspace=0.10, hspace=0.18)
    bb0 = axes[1,0].get_position(); bb1 = axes[1,1].get_position()
    cx = (bb0.x0 + bb1.x1) / 2
    fig.legend(handles=sb_legend_elements, ncol=5, loc="lower center",
                bbox_to_anchor=(cx, 0.01), frameon=False,
                fontsize=fs_legend, handlelength=1.4,
                handletextpad=0.4, columnspacing=1.0)
    fig.savefig(os.path.join(fig_dir, "Fig6_O3_calibration.png"),
                 dpi=500, bbox_inches="tight")
    plt.close(fig); print("  Fig6_O3_calibration.png")

# Figure 7: Validation boxplots (RMSE and bias) — same polishing
print("  Generating Fig7 ...")
with plt.rc_context({"font.size": fs_tick, "axes.labelsize": fs_label,
                      "axes.titlesize": fs_title, "xtick.labelsize": fs_tick,
                      "ytick.labelsize": fs_tick, "legend.fontsize": fs_legend}):
    fig, axes = plt.subplots(2, 2, figsize=(5.6, 3.6))
    sparse_xticklabels = [str(s) if i % 2 == 0 else "" for i, s in enumerate(sample_sizes)]
    for col_idx, shape in enumerate(["circle","square"]):
        sub = o3_df[o3_df["shape"]==shape]
        for vi, vl in enumerate(var_levels):
            vs = sub[sub["variance"]==vl]
            rmse_data = [vs[vs["n_cal"]==n]["val_RMSE"].values for n in sample_sizes]
            bias_data = [vs[vs["n_cal"]==n]["val_bias"].values for n in sample_sizes]
            positions = np.arange(len(sample_sizes)) + sb_offsets[vi]
            _styled_boxplot(axes[0,col_idx], rmse_data, positions, sb_var_colors[vl])
            _styled_boxplot(axes[1,col_idx], bias_data, positions, sb_var_colors[vl])
        for r in [0,1]:
            ax = axes[r,col_idx]
            ax.set_xticks(range(len(sample_sizes)))
            ax.set_xticklabels(sparse_xticklabels)
            _style_axes(ax)
        axes[0,col_idx].tick_params(labelbottom=False)
        axes[1,col_idx].set_xlabel("Number of Plots")
        axes[0,col_idx].set_ylim(0, 80); axes[0,col_idx].set_yticks([0, 20, 40, 60, 80])
        axes[1,col_idx].set_ylim(-40,40); axes[1,col_idx].set_yticks([-40,-20,0,20,40])
        axes[1,col_idx].axhline(0, color="gray", ls="--", lw=sb_spine_lw)
        if col_idx == 0:
            axes[0,col_idx].set_ylabel("RMSE (Mg ha$^{-1}$)")
            axes[1,col_idx].set_ylabel("bias (Mg ha$^{-1}$)")
        else:
            axes[0,col_idx].set_ylabel(""); axes[0,col_idx].set_yticklabels([])
            axes[1,col_idx].set_ylabel(""); axes[1,col_idx].set_yticklabels([])
        axes[0,col_idx].set_title(f"{shape.title()} Plots", fontsize=fs_title)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.96, bottom=0.22,
                         wspace=0.10, hspace=0.18)
    bb0 = axes[1,0].get_position(); bb1 = axes[1,1].get_position()
    cx = (bb0.x0 + bb1.x1) / 2
    fig.legend(handles=sb_legend_elements, ncol=5, loc="lower center",
                bbox_to_anchor=(cx, 0.01), frameon=False,
                fontsize=fs_legend, handlelength=1.4,
                handletextpad=0.4, columnspacing=1.0)
    fig.savefig(os.path.join(fig_dir, "Fig7_O3_validation.png"),
                 dpi=500, bbox_inches="tight")
    plt.close(fig); print("  Fig7_O3_validation.png")

#------------------------------------------------------------------------------
# FIGURE 11 — Hierarchical structure and pseudoreplication diagnostics.
# The dataset is built from 30 parent inventory plots (12.6 m radius). Within
# each parent, 8 sub-plot configurations were extracted at each (shape × size)
# combination, giving 240 sub-plots per (shape × size) cell but only 30 unique
# parents. A reviewer questioned whether sub-plots inside the same parent are
# spatially autocorrelated, which would inflate the effective sample size and
# constitute pseudoreplication ("n = 30, not 240"). Figure 11 quantifies this
# hierarchical structure and evaluates its impact with four complementary
# diagnostics:
#   (a-b) Intraclass correlation coefficient (ICC) for AGB and the mean LiDAR
#         feature as a function of sub-plot size. ICC = between-parent / total
#         variance: low values (≈0) indicate that within-parent sub-plots
#         sample largely independent forest patches, high values (→1) indicate
#         that the 8 sub-plots within a parent are near-replicates.
#   (c-d) For the 300 m² pool (the basis of Objective 3), calibration R² and
#         RMSE under three evaluation schemes:
#             • Random 5-fold CV  – the standard estimator used in Obj 3.
#             • Grouped 5-fold CV (leave-parent-out) – an unbiased estimator
#               that breaks within-parent autocorrelation by holding out all
#               sub-plots of 6 parents per fold.
#             • One sub-plot per parent (n = 30) – a strictly-independent
#               design that matches the reviewer's effective sample size.
#------------------------------------------------------------------------------
print("\n" + "-" * 70)
print("FIGURE 11 — Spatial autocorrelation / pseudoreplication diagnostics")
print("-" * 70)

from sklearn.model_selection import GroupKFold

y_colS1   = "AGB(Mg/ha)_Mean"
n_repsS1a = 30          # repetitions for random-vs-grouped CV
n_bootS1c = 100         # bootstrap iterations for one-per-parent control
icc_sizes = [100, 150, 200, 225, 250, 300]

# 300 m² pools per shape (used for panels c & d)
poolsS1 = {sh: df[(df["Shape"] == sh) & (df["Plot_size_m2"] == 300)]
                .reset_index(drop=True) for sh in ["circle", "square"]}
print(f"  300 m² pool sizes: circle n={len(poolsS1['circle'])}, "
      f"square n={len(poolsS1['square'])}, "
      f"unique parents per shape: "
      f"{poolsS1['circle']['parent_plot'].nunique()}/"
      f"{poolsS1['square']['parent_plot'].nunique()}")

#---- Panels (a-b): ICC of AGB and LiDAR features vs sub-plot size ----------
def _icc(g, col):
    """One-way ANOVA ICC: between-parent / total variance."""
    n_groups  = g["parent_plot"].nunique()
    grp_means = g.groupby("parent_plot")[col].mean()
    grp_sizes = g.groupby("parent_plot")[col].count()
    n_per     = grp_sizes.mean()
    mu        = g[col].mean()
    msb = ((grp_means - mu) ** 2 * grp_sizes).sum() / (n_groups - 1)
    msw = g.groupby("parent_plot")[col].apply(
              lambda x: ((x - x.mean()) ** 2).sum()).sum() / (len(g) - n_groups)
    return (msb - msw) / (msb + (n_per - 1) * msw)

icc_agb = {sh: [_icc(df[(df.Shape == sh) & (df.Plot_size_m2 == sz)],
                      "AGB(Mg/ha)_Mean") for sz in icc_sizes]
            for sh in ["circle", "square"]}
icc_lid = {sh: [np.mean([_icc(df[(df.Shape == sh) & (df.Plot_size_m2 == sz)], f)
                          for f in predictors]) for sz in icc_sizes]
            for sh in ["circle", "square"]}
for sh in ["circle", "square"]:
    print(f"  {sh:6s} ICC_AGB by size: "
          + ", ".join([f"{sz}m²={v:.2f}" for sz, v in zip(icc_sizes, icc_agb[sh])]))
    print(f"  {sh:6s} ICC_LiDAR (mean): "
          + ", ".join([f"{sz}m²={v:.2f}" for sz, v in zip(icc_sizes, icc_lid[sh])]))

#---- Panels (c-d): Random vs Grouped 5-fold CV at 300 m² ------------------
randomR2  = {sh: [] for sh in ["circle", "square"]}
groupedR2 = {sh: [] for sh in ["circle", "square"]}
randomRM  = {sh: [] for sh in ["circle", "square"]}
groupedRM = {sh: [] for sh in ["circle", "square"]}

for sh in ["circle", "square"]:
    pool = poolsS1[sh]
    x_pool = pool[predictors].values
    y_pool = pool[y_colS1].values
    g_pool = pool["parent_plot"].values
    n_groups = len(np.unique(g_pool))
    n_splits = min(5, n_groups)
    for rep in range(n_repsS1a):
        # Standard random 5-fold CV
        pred_r = cv_predict(x_pool, y_pool, n_folds=5, rs=seed + rep)
        r2r, rmser = cal_metrics(y_pool, pred_r)
        randomR2[sh].append(r2r); randomRM[sh].append(rmser)
        # Grouped 5-fold CV by parent_plot (deterministic given groups, but we
        # shuffle the group order via a permutation so we get a distribution).
        rng_rep = np.random.default_rng(seed + rep)
        perm = rng_rep.permutation(n_groups)
        g_map = {og: perm[i] for i, og in enumerate(np.unique(g_pool))}
        g_shuf = np.array([g_map[g] for g in g_pool])
        gkf = GroupKFold(n_splits=n_splits)
        pred_g = np.full(len(y_pool), np.nan)
        for tr, te in gkf.split(x_pool, y_pool, g_shuf):
            rf = fit_rf(x_pool[tr], y_pool[tr], rs=seed + rep)
            pred_g[te] = rf.predict(x_pool[te])
        r2g, rmseg = cal_metrics(y_pool, pred_g)
        groupedR2[sh].append(r2g); groupedRM[sh].append(rmseg)
    print(f"  {sh}: random  R²={np.mean(randomR2[sh]):.3f}±{np.std(randomR2[sh]):.3f}, "
          f"RMSE={np.mean(randomRM[sh]):.2f}")
    print(f"  {sh}: grouped R²={np.mean(groupedR2[sh]):.3f}±{np.std(groupedR2[sh]):.3f}, "
          f"RMSE={np.mean(groupedRM[sh]):.2f}")

#---- Panels (c-d): One sub-plot per parent (strictly independent) ----------
indepR2 = {sh: [] for sh in ["circle", "square"]}
indepRM = {sh: [] for sh in ["circle", "square"]}
for sh in ["circle", "square"]:
    pool = poolsS1[sh]
    parents = pool["parent_plot"].unique()
    for b in range(n_bootS1c):
        rng_b = np.random.default_rng(seed + b)
        sel_idx = []
        for p in parents:
            cand = pool.index[pool["parent_plot"] == p].to_numpy()
            sel_idx.append(int(rng_b.choice(cand, size=1)[0]))
        sel = pool.loc[sel_idx]
        x = sel[predictors].values
        y = sel[y_colS1].values
        if len(y) < 6: continue
        rf = RandomForestRegressor(n_estimators=n_trees, min_samples_leaf=min_leaf,
                                    max_features="sqrt", random_state=seed + b,
                                    n_jobs=-1, oob_score=True)
        rf.fit(x, y)
        oob = rf.oob_prediction_
        r2, rmse = cal_metrics(y, oob)
        indepR2[sh].append(r2); indepRM[sh].append(rmse)
    print(f"  {sh} one-per-parent (n={len(parents)}): "
          f"R²={np.mean(indepR2[sh]):.3f}±{np.std(indepR2[sh]):.3f}, "
          f"RMSE={np.mean(indepRM[sh]):.2f}")

# Save numerical summary
figS1_summary = pd.DataFrame([
    {"diagnostic": "random_5fold_CV", "shape": sh,
     "R2_mean": float(np.mean(randomR2[sh])), "R2_std": float(np.std(randomR2[sh])),
     "RMSE_mean": float(np.mean(randomRM[sh])), "RMSE_std": float(np.std(randomRM[sh])),
     "n_iter": n_repsS1a}
    for sh in ["circle", "square"]
] + [
    {"diagnostic": "grouped_5fold_CV_by_parent", "shape": sh,
     "R2_mean": float(np.mean(groupedR2[sh])), "R2_std": float(np.std(groupedR2[sh])),
     "RMSE_mean": float(np.mean(groupedRM[sh])), "RMSE_std": float(np.std(groupedRM[sh])),
     "n_iter": n_repsS1a}
    for sh in ["circle", "square"]
] + [
    {"diagnostic": "one_subplot_per_parent", "shape": sh,
     "R2_mean": float(np.mean(indepR2[sh])), "R2_std": float(np.std(indepR2[sh])),
     "RMSE_mean": float(np.mean(indepRM[sh])), "RMSE_std": float(np.std(indepRM[sh])),
     "n_iter": n_bootS1c}
    for sh in ["circle", "square"]
])
figS1_summary.to_csv(os.path.join(res_dir, "FigureS1_autocorrelation_summary.csv"),
                     index=False)

# ICC table
icc_tbl = pd.DataFrame({
    "Plot_size_m2":     icc_sizes,
    "ICC_AGB_circle":   icc_agb["circle"],
    "ICC_AGB_square":   icc_agb["square"],
    "ICC_LiDAR_circle": icc_lid["circle"],
    "ICC_LiDAR_square": icc_lid["square"]})
icc_tbl.to_csv(os.path.join(res_dir, "FigureS1_ICC_by_size.csv"), index=False)

#---- Render Figure 11 (2x2) ------------------------------------------------
print("  Generating Fig11 ...")
def _add_panel_label(ax, label, xy=(0.03, 0.95)):
    ax.text(xy[0], xy[1], label, transform=ax.transAxes,
             fontsize=13, va="top", ha="left")

# Figure 11 is the widest panel grid in the paper, so its text was the smallest
# once placed at column width. Sizes raised by about a third.
with plt.rc_context({"font.size": 12, "axes.labelsize": 12.5,
                      "xtick.labelsize": 11.5, "ytick.labelsize": 11.5,
                      "legend.fontsize": 10.5}):
    fig, axes = plt.subplots(2, 2, figsize=(7.6, 5.6))

    # Panel (a): ICC of AGB vs sub-plot size
    ax = axes[0, 0]
    ax.plot(icc_sizes, icc_agb["circle"], "o-", color="#1f77b4", lw=1.2, ms=5, label="Circle")
    ax.plot(icc_sizes, icc_agb["square"], "s-", color="#d95f02", lw=1.2, ms=5, label="Square")
    ax.set_xlabel("Sub-plot size (m$^2$)")
    ax.set_ylabel("ICC of AGB across parent plots")
    ax.set_ylim(0, 1.05); ax.set_xticks(icc_sizes)
    ax.grid(alpha=0.3, lw=0.4)
    ax.legend(loc="lower right", frameon=False)
    _add_panel_label(ax, "(a)")

    # Panel (b): ICC of mean LiDAR feature vs sub-plot size
    ax = axes[0, 1]
    ax.plot(icc_sizes, icc_lid["circle"], "o-", color="#1f77b4", lw=1.2, ms=5, label="Circle")
    ax.plot(icc_sizes, icc_lid["square"], "s-", color="#d95f02", lw=1.2, ms=5, label="Square")
    ax.set_xlabel("Sub-plot size (m$^2$)")
    ax.set_ylabel("Mean ICC of ALS Metrics")
    ax.set_ylim(0.55, 1.05); ax.set_yticks(np.arange(0.6, 1.01, 0.1))
    ax.set_xticks(icc_sizes)
    ax.grid(alpha=0.3, lw=0.4)
    ax.legend(loc="lower right", frameon=False)
    _add_panel_label(ax, "(b)")

    # Panel (c): Calibration R² at 300 m²
    ax = axes[1, 0]
    pos = [0, 1, 2, 3.6, 4.6, 5.6]
    data = [randomR2["circle"], groupedR2["circle"], indepR2["circle"],
            randomR2["square"], groupedR2["square"], indepR2["square"]]
    colors = ["#1f77b4", "#7fbf7b", "#d95f02"] * 2
    bp = ax.boxplot(data, positions=pos, widths=0.55, patch_artist=True,
                     showfliers=False,
                     medianprops=dict(color="black", lw=0.7),
                     boxprops=dict(lw=0.5), whiskerprops=dict(lw=0.5),
                     capprops=dict(lw=0.5))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.7)
    ax.set_xticks([1, 4.6]); ax.set_xticklabels(["Circle", "Square"])
    ax.set_ylabel("Calibration R$^2$ (300 m$^2$)")
    # headroom so the legend sits clear of every box and whisker
    ax.set_ylim(-0.1, 1.45); ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.axhline(0, color="gray", ls="--", lw=0.4)
    leg_c = [Patch(facecolor="#1f77b4", alpha=0.7, label="Random 5-fold CV"),
             Patch(facecolor="#7fbf7b", alpha=0.7, label="Grouped 5-fold CV"),
             Patch(facecolor="#d95f02", alpha=0.7, label="One sub-plot per parent")]
    ax.legend(handles=leg_c, loc="upper center", bbox_to_anchor=(0.5, 1.0),
              frameon=False, fontsize=10)
    _add_panel_label(ax, "(c)")

    # Panel (d): Calibration RMSE at 300 m² (no legend; shares with c)
    ax = axes[1, 1]
    data = [randomRM["circle"], groupedRM["circle"], indepRM["circle"],
            randomRM["square"], groupedRM["square"], indepRM["square"]]
    bp = ax.boxplot(data, positions=pos, widths=0.55, patch_artist=True,
                     showfliers=False,
                     medianprops=dict(color="black", lw=0.7),
                     boxprops=dict(lw=0.5), whiskerprops=dict(lw=0.5),
                     capprops=dict(lw=0.5))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.7)
    ax.set_xticks([1, 4.6]); ax.set_xticklabels(["Circle", "Square"])
    ax.set_ylabel("Calibration RMSE (Mg ha$^{-1}$, 300 m$^2$)")
    ax.set_ylim(0, 70)
    _add_panel_label(ax, "(d)")

    for ax in axes.ravel():
        for sp in ax.spines.values(): sp.set_linewidth(0.5)
        ax.tick_params(width=0.5)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.27, wspace=0.27)
    fig.savefig(os.path.join(fig_dir, "Fig11_O3_autocorrelation_test.png"),
                 dpi=500, bbox_inches="tight")
    plt.close(fig); print("  Fig11_O3_autocorrelation_test.png")

#------------------------------------------------------------------------------
print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print(f"  Results: {res_dir}")
print(f"  Figures: {fig_dir}")
print("=" * 70)

"""
ML_Comparison.py
Objective 4: Compare six machine-learning models (Random Forest, XGBoost,
Support Vector Regression, Cubist, K-Nearest Neighbour Regression and
Gaussian Process Regression) for ALS-based AGB calibration as a function
of plot shape, calibration sample size, and biomass variance, using the
300 m² pseudo-plots from Prepared_Dataset.csv.

Data: Prepared_Dataset.csv filtered to Plot_size_m2 == 300 m² (200 pseudo-plots
      per shape).  For each shape, 30 random pseudo-plots are drawn aside as
      the independent validation set and the remaining 170 form the
      calibration candidate pool.

Analysis 1 — Sample-size effect.  At each n_cal in {40, 60, 80, 100, 120,
             140, 160}, n_cal calibration plots are drawn at random (without
             replacement) from the shape's calibration pool.  Six models are
             fitted per (shape, n_cal).  Drives Figures 8 (calibration) and
             9 (validation).

Analysis 2 — Variance sensitivity.  At each n_cal, 150 random candidate
             calibration sets are drawn, ranked by AGB variance, and
             split into five quintile groups of 30 sets each.  Six models
             are fitted on every set; calibration (R², RMSE) and
             validation (RMSE, bias, paired-t p-value) are recorded.

Performance synthesis — Composite Performance Index (CPI), Variance
             Sensitivity Index (VSI), Variance Robustness Index (VRI)
             aggregated across the variance levels.  Drives Figure 10
             (VRI line plot) and Figure S2 (CPI heatmaps).

Author: Xinyuan Wei
"""

import os
# Silence all warnings: must be done BEFORE importing sklearn / xgboost so
# child processes spawned by joblib (n_jobs > 1) inherit the setting.
os.environ.setdefault("PYTHONWARNINGS", "ignore")
import warnings
warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
import contextlib, io as _io
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.compose import TransformedTargetRegressor
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, WhiteKernel
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import r2_score, mean_squared_error
from xgboost import XGBRegressor
from cubist import Cubist
import cubist._make_names_string as _cmns
import cubist._make_data_string as _cmds
import re as _re_mod
def _escapes_safe(x):
    """NaN / float-safe replacement for cubist._escapes — coerce to str first.
    Older pandas versions leave NaN as floats inside Series.astype(str), which
    breaks the original implementation's '.replace()' call.  Coercing here
    fixes the issue without touching Cubist's installed source."""
    chars = ["\\", '"', "\'"]
    x = [str(c) for c in x]
    for i in chars:
        x = [c.replace(i, "\\" + i) for c in x]
    return [_re_mod.escape(c) for c in x]
_cmns._escapes = _escapes_safe
_cmds._escapes = _escapes_safe                            # fix the cached import too
from sklearn.base import BaseEstimator, RegressorMixin

class SafeCubist(BaseEstimator, RegressorMixin):
    """Wrapper around cubist.Cubist that always converts input to a pandas
    DataFrame before fit/predict.  Avoids 'float object has no attribute
    replace' errors that occur when the underlying package receives a plain
    NumPy array (a known incompatibility with newer pandas/numpy versions)."""
    def __init__(self, n_committees=1, n_rules=100, random_state=None):
        self.n_committees = n_committees
        self.n_rules = n_rules
        self.random_state = random_state
    @staticmethod
    def _to_df(X):
        if isinstance(X, pd.DataFrame):
            return X
        return pd.DataFrame(X, columns=[f"x{i}" for i in range(X.shape[1])])
    def fit(self, X, y):
        self._inner = Cubist(n_committees=self.n_committees,
                              n_rules=self.n_rules,
                              random_state=self.random_state)
        ys = y if isinstance(y, pd.Series) else pd.Series(np.asarray(y), name="y")
        # Cubist prints diagnostic messages directly to stdout/stderr,
        # bypassing the Python warnings system; redirect them to dev/null.
        with contextlib.redirect_stdout(_io.StringIO()), \
             contextlib.redirect_stderr(_io.StringIO()):
            self._inner.fit(self._to_df(X), ys)
        return self
    def predict(self, X):
        with contextlib.redirect_stdout(_io.StringIO()), \
             contextlib.redirect_stderr(_io.StringIO()):
            return np.asarray(self._inner.predict(self._to_df(X)))

from scipy.stats import ttest_rel
# Suppress XGBoost C-level chatter (covers a few warnings that bypass Python)
try:
    import xgboost as _xgb; _xgb.set_config(verbosity=0)
except Exception:
    pass
warnings.filterwarnings("ignore")

#------------------------------------------------------------------------------
# Paths & parameters
data_dir = os.path.dirname(os.path.abspath(__file__))
lidar_dir = os.path.join(data_dir, "Demeritt_ALS")
res_dir = os.path.join(data_dir, "Results")
# Figures now live in the manuscript folder: ../2.Final Paper/1.Figure
fig_dir = os.path.abspath(os.path.join(data_dir, "..", "2.Final Paper", "1.Figure"))
ml_data_file = os.path.join(data_dir, "Prepared_Dataset.csv")
os.makedirs(res_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

seed = 42                                    # fixed so the published run is reproducible
ml_plot_size = 300                           # m² (filter applied to the loaded dataset)
n_val = 30                                   # validation plots per shape
sample_sizes = [40, 60, 80, 100, 120, 140, 160]
var_levels = ["low", "low-medium", "medium", "medium-high", "high"]
n_iter = 30
n_candidates_var = 150                       # candidate cal sets per (shape, n_cal)
sets_per_level = n_candidates_var // len(var_levels)   # = 30
predictors = ["Point_density", "Vegetation_density", "Max_height",
              "Mean_height", "P25", "P50", "P75", "P95"]
y_col = "AGB(Mg/ha)_Mean"       # mean of methods J, Y, C and W

model_names = ["RF", "XGBoost", "SVR", "Cubist", "KNN", "GPR"]
model_colors = {"RF":"#000000", "XGBoost":"#E41A1C", "SVR":"#FFD700",
                "Cubist":"#4CAF50", "KNN":"#87CEFA", "GPR":"#0000FF"}

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"],
    "font.size": 9, "axes.labelsize": 10,
    "axes.titlesize": 11, "figure.dpi": 500,
    "mathtext.fontset": "stix",
})

#------------------------------------------------------------------------------
# Load prepared dataset and split into calibration pool / validation set
print("Loading Prepared_Dataset.csv ...")
df = pd.read_csv(ml_data_file)
df = df[df["Plot_size_m2"] == ml_plot_size].copy()
print(f"  {len(df)} rows at {ml_plot_size} m²; shapes={df['Shape'].unique().tolist()}")

val_pool = {}; cal_pool = {}
rng_split = np.random.default_rng(seed)
for shape in ["circle", "square"]:
    sub = df[df["Shape"] == shape].reset_index(drop=True)
    idx = np.arange(len(sub)); rng_split.shuffle(idx)
    val_idx = sorted(idx[:n_val]); cal_idx = sorted(idx[n_val:])
    val_pool[shape] = sub.iloc[val_idx].reset_index(drop=True)
    cal_pool[shape] = sub.iloc[cal_idx].reset_index(drop=True)
    print(f"  {shape}: validation={len(val_pool[shape])}, calibration pool={len(cal_pool[shape])}")

def random_select(data, n_cal, rs):
    """Pick n_cal calibration plots at random (without replacement) from
    `data`.  Returns positional indices into the (already reset_index'd)
    calibration pool DataFrame."""
    rng = np.random.default_rng(rs)
    if n_cal >= len(data): return np.arange(len(data))
    return rng.choice(len(data), size=n_cal, replace=False)

#------------------------------------------------------------------------------
# Model factory & metric helpers
def make_models(rs):
    # XGBoost and GPR are tamed (n_estimators=50, max_depth=3 for XGBoost;
    # WhiteKernel noise term for GPR) so they no longer interpolate training data
    return {
        "RF":      RandomForestRegressor(n_estimators=100, min_samples_leaf=3,
                                          max_features="sqrt", random_state=rs, n_jobs=1),
        "XGBoost": XGBRegressor(n_estimators=50, max_depth=3, learning_rate=0.1,
                                 random_state=rs, n_jobs=1, verbosity=0),
        # C and epsilon are absolute on the response scale, so the response is
        # standardised inside each fit; without this the same model would give a
        # different answer purely because AGB is reported in Mg/ha rather than
        # kg m-2. Predictions are back-transformed automatically.
        "SVR":     TransformedTargetRegressor(
                       regressor=Pipeline([("sc", StandardScaler()),
                                           ("m",  SVR(kernel="rbf", C=1.0,
                                                      epsilon=0.1))]),
                       transformer=StandardScaler()),
        "Cubist":  SafeCubist(n_committees=1, n_rules=100, random_state=rs),
        "KNN":     Pipeline([("sc", StandardScaler()),
                             ("m",  KNeighborsRegressor(n_neighbors=5))]),
        "GPR":     Pipeline([("sc", StandardScaler()),
                             ("m",  GaussianProcessRegressor(
                                 kernel=ConstantKernel(1.0) * RBF(length_scale=1.0)
                                       + WhiteKernel(noise_level=0.5),
                                 random_state=rs, normalize_y=True))]),
    }

def cal_metrics(y, pred):
    # Negative R² preserved (no clipping)
    return r2_score(y, pred), np.sqrt(mean_squared_error(y, pred))

def val_metrics(y, pred):
    rmse = np.sqrt(mean_squared_error(y, pred))
    bias = float(np.mean(pred - y))
    _, p = ttest_rel(pred, y)
    return rmse, bias, float(p)

#------------------------------------------------------------------------------
# ANALYSIS 1: Sample size effect (random calibration sets, 30 iterations each)
# 6 models  x  7 sizes  x  2 shapes  x  30 iterations  =  2,520 fitted models

print("\n" + "-" * 70)
print("ANALYSIS 1 — Sample size effect (random calibration, 30 iterations)")
print("-" * 70)

n_iter_ss = 30                                # iterations per (shape, n_cal)
ss_cal_rows = []; ss_val_rows = []

for shape in ["circle", "square"]:
    cdf = cal_pool[shape]; vdf = val_pool[shape]
    n_total = len(cdf)
    x_pool = cdf[predictors].values; y_pool = cdf[y_col].values
    x_val  = vdf[predictors].values; y_val  = vdf[y_col].values

    for n_cal in sample_sizes:
        for it in range(n_iter_ss):
            # Distinct seed per (shape, n_cal, iter) so every iteration draws
            # an independent random calibration set.  Pass None when seed is
            # None so the RNG uses system entropy.
            iter_rs = (None if seed is None
                       else seed + abs(hash((shape, n_cal, it))) % 1_000_000)
            sel_idx = random_select(cdf, n_cal, rs=iter_rs)
            x_c = x_pool[sel_idx]; y_c = y_pool[sel_idx]
            sample_var = float(np.var(y_c))

            cv = KFold(n_splits=5, shuffle=True, random_state=iter_rs)
            models = make_models(rs=iter_rs)
            for name in model_names:
                m = models[name]
                pc = cross_val_predict(m, x_c, y_c, cv=cv)
                r2, rmse = cal_metrics(y_c, pc)
                ss_cal_rows.append({"shape":shape, "n_cal":n_cal, "model":name,
                                     "iter":it, "sample_var":sample_var,
                                     "R2":r2, "RMSE":rmse})
                m.fit(x_c, y_c); pv = m.predict(x_val)
                vrmse, vbias, vp = val_metrics(y_val, pv)
                ss_val_rows.append({"shape":shape, "n_cal":n_cal, "model":name,
                                     "iter":it, "RMSE":vrmse, "bias":vbias,
                                     "abs_bias":abs(vbias), "p_value":vp})
        print(f"  {shape}/n={n_cal}: {n_iter_ss} iter × 6 models = "
              f"{n_iter_ss * 6} fits")

ss_cal = pd.DataFrame(ss_cal_rows)
ss_val = pd.DataFrame(ss_val_rows)
ss_cal.to_csv(os.path.join(res_dir, "Objective4_SampleSize_Calibration.csv"), index=False)
ss_val.to_csv(os.path.join(res_dir, "Objective4_SampleSize_Validation.csv"), index=False)
print(f"  Saved {len(ss_cal)} sample-size calibration rows, {len(ss_val)} validation rows")

#------------------------------------------------------------------------------
# Figure helpers
spine_lw = 0.5
def style_axes(ax):
    for s in ax.spines.values(): s.set_linewidth(spine_lw)
    ax.tick_params(width=spine_lw)

def heatmap(ax, mat, row_labels, col_labels, cmap, vmin, vmax,
            divider_after=None, annot_fmt="{:.2f}", annot_size=6,
            text_color_threshold=None):
    """Draw a numeric heatmap with annotations and an internal divider line."""
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_xticklabels(col_labels, fontsize=8)
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_xticks(np.arange(mat.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(mat.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=0.4)
    ax.tick_params(which="minor", length=0)
    style_axes(ax)
    if divider_after is not None:
        ax.axhline(divider_after + 0.5, color="black", lw=0.8)
    norm = (mat - vmin) / max(vmax - vmin, 1e-12)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            if text_color_threshold is None:
                tc = "white" if norm[i, j] > 0.55 else "black"
            else:
                tc = text_color_threshold(v)
            ax.text(j, i, annot_fmt.format(v), ha="center", va="center",
                    color=tc, fontsize=annot_size)
    return im

def build_matrix(df, value_col, model_order, shapes_order, n_cals):
    """Stack rows = (shape × model) order: all circles first then squares."""
    rows = []; row_labels = []
    for shape in shapes_order:
        for m in model_order:
            r = []
            for n in n_cals:
                v = df[(df["shape"]==shape) & (df["model"]==m) & (df["n_cal"]==n)][value_col]
                r.append(v.values[0] if len(v) else np.nan)
            rows.append(r)
            row_labels.append(f"{m} ({shape[0].upper()})")
    return np.array(rows), row_labels

#------------------------------------------------------------------------------
# Figures 7 & 8 — colour matrices
from matplotlib.colors import LinearSegmentedColormap
cmap_sby = LinearSegmentedColormap.from_list("skyblue_yellow",
                                              ["#87CEFA", "#FFD700"])
cmap_sby_div = LinearSegmentedColormap.from_list("skyblue_white_yellow",
                                                  ["#87CEFA", "#FFFFFF", "#FFD700"])

def shape_matrix(df, value_col, models, shape, n_cals):
    """Per (model, n_cal) cell, return the MEAN value across iterations."""
    rows = []
    for m in models:
        r = []
        for n in n_cals:
            v = df[(df["shape"]==shape) & (df["model"]==m) & (df["n_cal"]==n)][value_col]
            r.append(float(v.mean()) if len(v) else np.nan)
        rows.append(r)
    return np.array(rows)

def heatmap_panel(ax, mat, row_labels, col_labels, cmap, vmin, vmax,
                   show_yticks=True, annot_fmt="{:.2f}", annot_size=8.0):
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(mat.shape[1])); ax.set_xticklabels(col_labels, fontsize=8)
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels(row_labels if show_yticks else [""]*mat.shape[0], fontsize=8)
    ax.set_xticks(np.arange(mat.shape[1]+1)-0.5, minor=True)
    ax.set_yticks(np.arange(mat.shape[0]+1)-0.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=0.4); ax.tick_params(which="minor", length=0)
    for s in ax.spines.values(): s.set_linewidth(0.5)
    ax.tick_params(width=0.5)
    norm = (mat - vmin) / max(vmax - vmin, 1e-12)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            tc = "white" if norm[i, j] > 0.55 else "black"
            ax.text(j, i, annot_fmt.format(mat[i, j]), ha="center", va="center",
                    color=tc, fontsize=annot_size)
    return im

#------------------------------------------------------------------------------
# Figure 8: Calibration R² and RMSE — 2 rows (R², RMSE) × 2 columns (Circle, Square).
# Each panel is a line chart with one line per ML model.  The line shows the
# mean across the 30 iterations and the shaded band is ±1 SD.
print("Generating Fig8 (line plots: mean ± SD across iterations) ...")

shape_label = {"circle": "Circle", "square": "Square"}
fig7_marker = {"circle": "o", "square": "s"}

# Smaller panels with larger relative text: shrink figsize and bump fonts via rc_context.
fig7_fs = {"font.size": 10, "axes.labelsize": 11, "axes.titlesize": 12,
           "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9}
sparse_fig7_xticklabels = [str(s) if i % 2 == 0 else "" for i, s in enumerate(sample_sizes)]

with plt.rc_context(fig7_fs):
    fig, axes = plt.subplots(2, 2, figsize=(6.4, 4.4),
                              gridspec_kw={"hspace": 0.12, "wspace": 0.08})

    for row_idx, metric in enumerate(["R2", "RMSE"]):
        for col_idx, shape in enumerate(["circle", "square"]):
            ax = axes[row_idx, col_idx]
            sub_all = ss_cal[ss_cal["shape"] == shape]
            for mn in model_names:
                sub = sub_all[sub_all["model"] == mn]
                grp = sub.groupby("n_cal")[metric]
                means = grp.mean().reindex(sample_sizes).values
                stds  = grp.std().reindex(sample_sizes).values
                ax.errorbar(sample_sizes, means, yerr=stds,
                            color=model_colors[mn], lw=1.2,
                            marker=fig7_marker[shape], markersize=5.0,
                            markeredgecolor="black", markeredgewidth=0.4,
                            capsize=2.0, capthick=0.8, elinewidth=0.8,
                            label=mn)
            ax.set_xticks(sample_sizes)
            ax.set_xticklabels(sparse_fig7_xticklabels)
            for s in ax.spines.values(): s.set_linewidth(0.5)
            ax.tick_params(width=0.5)
            if row_idx == 0:
                ax.set_title(f"{shape_label[shape]} Plots")
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("Number of Plots")
            if col_idx == 0:
                ax.set_ylabel("R²" if metric == "R2" else "RMSE (Mg ha$^{-1}$)")
            else:
                ax.tick_params(labelleft=False)

    # Y-limits / ticks per row
    for col_idx in range(2):
        axes[0, col_idx].set_ylim(0.0, 1.0)
        axes[0, col_idx].set_yticks(np.arange(0.0, 1.01, 0.2))
        axes[1, col_idx].set_ylim(0.0, 80.0)
        axes[1, col_idx].set_yticks([0, 20, 40, 60, 80])

    # Legends on both bottom-row panels (Circle and Square), upper-right corner
    for col_idx in (0, 1):
        axes[1, col_idx].legend(loc="upper right", frameon=False, ncol=2,
                                 handletextpad=0.4, columnspacing=0.8)

    fig.subplots_adjust(left=0.10, right=0.98, top=0.97, bottom=0.10)
    fig.savefig(os.path.join(fig_dir, "Fig8_O4_SS_calibration.png"),
                dpi=500, bbox_inches="tight")
    plt.close(fig); print("  Fig8_O4_SS_calibration.png")

#------------------------------------------------------------------------------
# Figure 9: Validation RMSE / bias / p-value (2 rows = shape, 3 cols = metric)
print("Generating Fig9 (split circle/square, bottom colorbars) ...")
col_labels = [str(n) for n in sample_sizes]
rmse_vmin, rmse_vmax = 10, 70
bias_vmin, bias_vmax = -16, 16

fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.2),
                          gridspec_kw={"hspace":0.08, "wspace":0.04})
im_handles = {"rmse":[], "bias":[], "p":[]}
for ri, shape in enumerate(["circle","square"]):
    mat_rmse = shape_matrix(ss_val, "RMSE",     model_names, shape, sample_sizes)
    mat_bias = shape_matrix(ss_val, "bias",     model_names, shape, sample_sizes)
    mat_p    = shape_matrix(ss_val, "p_value",  model_names, shape, sample_sizes)
    im_handles["rmse"].append(heatmap_panel(axes[ri,0], mat_rmse, model_names, col_labels,
                                              "YlGn_r", rmse_vmin, rmse_vmax,
                                              show_yticks=True))
    im_handles["bias"].append(heatmap_panel(axes[ri,1], mat_bias, model_names, col_labels,
                                              "RdBu_r", bias_vmin, bias_vmax,
                                              show_yticks=False))
    im_handles["p"].append(heatmap_panel(axes[ri,2], mat_p, model_names, col_labels,
                                          "YlGn", 0.0, 1.0, show_yticks=False))
    if ri == 0:
        axes[ri, 0].set_title("RMSE (Mg ha$^{-1}$)", fontsize=11)
        axes[ri, 1].set_title("bias (Mg ha$^{-1}$)",  fontsize=11)
        axes[ri, 2].set_title("Paired t-test p-value",     fontsize=11)
    if ri == 1:
        axes[ri, 0].set_xlabel("Number of Plots")
        axes[ri, 1].set_xlabel("Number of Plots")
        axes[ri, 2].set_xlabel("Number of Plots")
    axes[ri, 0].set_ylabel(shape_label[shape], fontsize=11, labelpad=8)
    axes[ri, 1].set_ylabel(""); axes[ri, 2].set_ylabel("")
    if ri == 0:                                    # hide x-tick labels on upper row
        for c_ in range(3): axes[ri, c_].tick_params(labelbottom=False)
fig.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.18)
bb0 = axes[1, 0].get_position()
bb1 = axes[1, 1].get_position()
bb2 = axes[1, 2].get_position()
# Inset each colorbar slightly so the end tick of one does not touch the start
# tick of the next ("3.5" and "0.8" otherwise print as "3.50.8").
_cb_pad = 0.018
cax1 = fig.add_axes([bb0.x0 + _cb_pad, 0.07, bb0.width - 2 * _cb_pad, 0.022])
cax2 = fig.add_axes([bb1.x0 + _cb_pad, 0.07, bb1.width - 2 * _cb_pad, 0.022])
cax3 = fig.add_axes([bb2.x0 + _cb_pad, 0.07, bb2.width - 2 * _cb_pad, 0.022])
cb1 = fig.colorbar(im_handles["rmse"][0], cax=cax1, orientation="horizontal", label="RMSE (Mg ha$^{-1}$)")
cb2 = fig.colorbar(im_handles["bias"][0], cax=cax2, orientation="horizontal", label="bias (Mg ha$^{-1}$)")
cb3 = fig.colorbar(im_handles["p"][0],    cax=cax3, orientation="horizontal", label="p-value")
cb1.set_ticks(np.arange(rmse_vmin, rmse_vmax + 1e-9, 10))
cb2.set_ticks(np.arange(bias_vmin, bias_vmax + 1e-9, 8))
for cb in (cb1, cb2, cb3):
    cb.outline.set_linewidth(0.5); cb.ax.tick_params(width=0.5)
fig.savefig(os.path.join(fig_dir, "Fig9_O4_SS_validation.png"),
            dpi=500, bbox_inches="tight")
plt.close(fig); print("  Fig9_O4_SS_validation.png")

#------------------------------------------------------------------------------
# ANALYSIS 2: Variance sensitivity (12,600 fitted models)

vs_cal_path = os.path.join(res_dir, "Objective4_Calibration_Results.csv")
vs_val_path = os.path.join(res_dir, "Objective4_Validation_Results.csv")
def _nonempty(p): return os.path.exists(p) and os.path.getsize(p) > 0
if _nonempty(vs_cal_path) and _nonempty(vs_val_path):
    vs_cal = pd.read_csv(vs_cal_path)
    vs_val = pd.read_csv(vs_val_path)
    print(f"Loaded existing variance-sensitivity results: "
          f"cal={len(vs_cal)} rows, val={len(vs_val)} rows")
else:
    print("\n" + "-" * 70)
    print("ANALYSIS 2 — Variance sensitivity")
    print("-" * 70)
    cal_rows = []; val_rows = []
    rng = np.random.default_rng(seed)
    for shape in ["circle", "square"]:
        cdf = cal_pool[shape]; vdf = val_pool[shape]
        n_total = len(cdf)
        x_pool = cdf[predictors].values; y_pool = cdf[y_col].values
        x_val  = vdf[predictors].values; y_val  = vdf[y_col].values
        for n_cal in sample_sizes:
            # Draw n_candidates_var random calibration sets, rank them by AGB
            # variance, and split into 5 quintile groups of sets_per_level each.
            candidates = []
            for _ in range(n_candidates_var):
                idx = rng.choice(n_total, size=n_cal, replace=False)
                candidates.append((idx, float(np.var(y_pool[idx]))))
            candidates.sort(key=lambda t: t[1])
            groups = {vl: [c[0] for c in
                            candidates[i*sets_per_level:(i+1)*sets_per_level]]
                       for i, vl in enumerate(var_levels)}
            for vl in var_levels:
                pool_subs = groups[vl]
                actual_iter = len(pool_subs)
                if actual_iter < 5: continue
                for it in range(actual_iter):
                    cal_idx = pool_subs[it]
                    x_c, y_c = x_pool[cal_idx], y_pool[cal_idx]
                    rs_it = None if seed is None else seed + it
                    cv = KFold(n_splits=5, shuffle=True, random_state=rs_it)
                    models = make_models(rs=rs_it)
                    for name in model_names:
                        m = models[name]
                        pc = cross_val_predict(m, x_c, y_c, cv=cv)
                        r2, rmse = cal_metrics(y_c, pc)
                        cal_rows.append({"shape":shape, "n_cal":n_cal, "variance":vl,
                                          "model":name, "iter":it, "R2":r2, "RMSE":rmse})
                        m.fit(x_c, y_c); pv = m.predict(x_val)
                        vrmse, vbias, vp = val_metrics(y_val, pv)
                        val_rows.append({"shape":shape, "n_cal":n_cal, "variance":vl,
                                          "model":name, "iter":it,
                                          "RMSE":vrmse, "bias":vbias, "p_value":vp})
                print(f"  {shape}/n={n_cal}/{vl}: {actual_iter} iter done")
    vs_cal = pd.DataFrame(cal_rows); vs_val = pd.DataFrame(val_rows)
    vs_cal.to_csv(vs_cal_path, index=False)
    vs_val.to_csv(vs_val_path, index=False)

#------------------------------------------------------------------------------
# STEP 3: Composite Performance Index (CPI), Variance Sensitivity Index (VSI),
print("Computing CPI, VSI, VRI ...")

# Merge cal and val rows on (shape, n_cal, variance, model, iter)
key = ["shape","n_cal","variance","model","iter"]
merged = vs_cal.merge(vs_val.rename(columns={"RMSE":"RMSE_val"}),
                       on=key, suffixes=("_cal", "_val"))
merged = merged.rename(columns={"RMSE":"RMSE_cal", "R2":"R2_cal"})
merged["abs_bias"] = merged["bias"].abs()

# Min-max normalisation across all rows (one normalisation table for the whole experiment)
def norm_higher_better(s):
    lo, hi = s.min(), s.max()
    if hi - lo < 1e-12: return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - lo) / (hi - lo)
def norm_lower_better(s):
    lo, hi = s.min(), s.max()
    if hi - lo < 1e-12: return pd.Series(np.ones(len(s)), index=s.index)
    return 1.0 - (s - lo) / (hi - lo)

merged["S_R2"]       = norm_higher_better(merged["R2_cal"])
merged["S_RMSEcal"]  = norm_lower_better(merged["RMSE_cal"])
merged["S_RMSEval"]  = norm_lower_better(merged["RMSE_val"])
merged["S_bias"]     = norm_lower_better(merged["abs_bias"])
merged["S_p"]        = (merged["p_value"] / 0.05).clip(upper=1.0)

w_R2, w_RMSEcal, w_RMSEval, w_bias, w_p = 0.15, 0.10, 0.35, 0.25, 0.15
merged["CPI"] = (w_R2*merged["S_R2"] + w_RMSEcal*merged["S_RMSEcal"]
                 + w_RMSEval*merged["S_RMSEval"] + w_bias*merged["S_bias"]
                 + w_p*merged["S_p"])

# Mean CPI per (shape, n_cal, variance, model) (aggregating over 30 iterations)
cpi_mean = (merged.groupby(["shape","n_cal","variance","model"])["CPI"]
                  .mean().reset_index().rename(columns={"CPI":"CPI_mean"}))

# VSI per (shape, n_cal, model): SD of 5 variance-level means / Mean of those 5
def vsi_row(g):
    vals = g.set_index("variance").reindex(var_levels)["CPI_mean"].values
    sd = float(np.std(vals, ddof=0)); mu = float(np.mean(vals))
    return pd.Series({"VSI": sd / mu if mu > 1e-12 else np.nan,
                       "CPI_avg": mu})
vsi = (cpi_mean.groupby(["shape","n_cal","model"])
                .apply(vsi_row).reset_index())

# VRI: 1 - min-max-normalised VSI (over all (shape, n_cal, model) rows)
vsi_min, vsi_max = vsi["VSI"].min(), vsi["VSI"].max()
vsi["VRI"] = 1.0 - (vsi["VSI"] - vsi_min) / max(vsi_max - vsi_min, 1e-12)

merged.to_csv(os.path.join(res_dir, "Objective4_CPI_per_iteration.csv"), index=False)
cpi_mean.to_csv(os.path.join(res_dir, "Objective4_CPI_mean.csv"), index=False)
vsi.to_csv(os.path.join(res_dir, "Objective4_VSI_VRI.csv"), index=False)
print(f"  Wrote CPI per-iteration ({len(merged)} rows), CPI mean ({len(cpi_mean)} rows), "
      f"VSI/VRI ({len(vsi)} rows)")

#------------------------------------------------------------------------------
# Figure 10: VRI line plots only (no CPI heatmaps).
print("Generating Fig10 (VRI line plots only) ...")
fig9_colors = dict(model_colors)             # same palette as Figure 8
fig9_marker = {"circle":"o", "square":"s"}

fig, axes = plt.subplots(1, 2, figsize=(8, 4.2), sharey=True,
                          gridspec_kw={"wspace":0.05})
for ax, shape in zip(axes, ["circle", "square"]):
    mk = fig9_marker[shape]
    for mn in model_names:
        sub = vsi[(vsi["shape"]==shape) & (vsi["model"]==mn)].sort_values("n_cal")
        ax.plot(sub["n_cal"], sub["VRI"], "-", color=fig9_colors[mn], lw=1.2,
                marker=mk, markersize=5.5, markeredgecolor="black",
                markeredgewidth=0.4, label=mn)
    ax.set_xlabel("Number of Plots")
    ax.set_xticks(sample_sizes)
    ax.set_ylim(0, 1.02); ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_title(f"{shape.title()} Plots", fontsize=11)
    for s in ax.spines.values(): s.set_linewidth(0.5)
    ax.tick_params(width=0.5)
axes[0].set_ylabel("Variance Robustness Index (VRI)")
for ax in axes: ax.set_box_aspect(1)
axes[0].legend(fontsize=8, loc="lower right", frameon=False, ncol=2)
axes[1].legend(fontsize=8, loc="lower right", frameon=False, ncol=2)
fig.tight_layout(rect=[0, 0, 1, 1], w_pad=0.3)
fig.savefig(os.path.join(fig_dir, "Fig10_O4_indices.png"),
            dpi=500, bbox_inches="tight")
plt.close(fig); print("  Fig10_O4_indices.png")

#------------------------------------------------------------------------------
# Figure S2 — heatmap of mean CPI for every (shape, sample size, model, variance).
# 2 rows (shape) x 7 cols (sample size); each cell is a 6-models x 5-variance heatmap.
print("Generating FigS2 (CPI heatmaps) ...")
import matplotlib.gridspec as gridspec
sb_var_lbl_short = ["L","LM","M","MH","H"]
fig = plt.figure(figsize=(14, 7.0))
gs  = gridspec.GridSpec(2, 7, left=0.05, right=0.94, top=0.92, bottom=0.06,
                         hspace=0.05, wspace=0.08)
shape_label_s1 = {"circle":"Circle", "square":"Square"}
last_im = None
for ri, shape in enumerate(["circle","square"]):
    for ci, n_cal in enumerate(sample_sizes):
        ax = fig.add_subplot(gs[ri, ci])
        sub = cpi_mean[(cpi_mean["shape"]==shape) & (cpi_mean["n_cal"]==n_cal)]
        mat = np.full((len(model_names), len(var_levels)), np.nan)
        for mi, mn in enumerate(model_names):
            row = sub[sub["model"]==mn].set_index("variance")
            for vi, vl in enumerate(var_levels):
                if vl in row.index: mat[mi, vi] = row.loc[vl, "CPI_mean"]
        last_im = ax.imshow(mat, aspect="auto", cmap="viridis", vmin=0, vmax=1)
        ax.set_xticks(range(len(var_levels)))
        ax.set_xticklabels(sb_var_lbl_short if ri == 1 else [""]*len(var_levels), fontsize=8)
        ax.set_yticks(range(len(model_names)))
        ax.set_yticklabels(model_names if ci == 0 else [""]*len(model_names), fontsize=8)
        for sp in ax.spines.values(): sp.set_linewidth(0.5)
        ax.tick_params(width=0.5, length=2)
        if ri == 0: ax.set_title(f"n = {n_cal}", fontsize=9)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if not np.isnan(mat[i, j]):
                    tc = "white" if mat[i, j] < 0.55 else "black"
                    ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center",
                             color=tc, fontsize=7)
        if ci == 0:
            ax.set_ylabel(shape_label_s1[shape], fontsize=10, labelpad=8)
cax = fig.add_axes([0.95, 0.20, 0.012, 0.55])
cb  = fig.colorbar(last_im, cax=cax)
cb.set_label("CPI", fontsize=9); cb.outline.set_linewidth(0.5)
cb.ax.tick_params(width=0.5, labelsize=8)
fig.savefig(os.path.join(fig_dir, "FigS2_O4_CPI_heatmap.png"),
            dpi=500, bbox_inches="tight")
plt.close(fig); print("  FigS2_O4_CPI_heatmap.png")

#------------------------------------------------------------------------------
print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print(f"  Results: {res_dir}")
print(f"  Figures: {fig_dir}")
print("=" * 70)

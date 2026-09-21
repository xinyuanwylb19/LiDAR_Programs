# ML_Comparison.py
# Author: Xinyuan Wei.
# Sep. 02 2026

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import xgboost
import cubist._make_names_string as cubist_names
import cubist._make_data_string as cubist_data
from cubist import Cubist
from scipy.stats import ttest_rel
from sklearn.compose import TransformedTargetRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, WhiteKernel
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from xgboost import XGBRegressor

matplotlib.use("Agg")
os.environ.setdefault("PYTHONWARNINGS", "ignore")   # Inherited by the joblib workers.
warnings.filterwarnings("ignore")
xgboost.set_config(verbosity=0)

# Cubist hands float NaN to its text-escaping helper when it predicts, which
# fails on some versions, so values are turned into text first.
def cubist_escapes(x, escapes=cubist_names._escapes):
    return escapes([str(c) for c in x])

cubist_names._escapes = cubist_escapes
cubist_data._escapes = cubist_escapes

#------------------------------------------------------------------------------
# Paths and parameters.
data_dir = os.path.dirname(os.path.abspath(__file__))
res_dir = os.path.join(data_dir, "Results")
fig_dir = os.path.abspath(os.path.join(data_dir, "..", "2.Final Paper", "1.Figure"))
data_file = os.path.join(data_dir, "Prepared_Dataset.csv")
os.makedirs(res_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

seed = 42
plot_size = 300                              # Pseudo-plot size used (m2).
plot_shapes = ["circle", "square"]
n_val = 30                                   # Validation plots per shape.
n_iter = 30                                  # Random calibration sets per shape and sample size.
sample_sizes = [40, 60, 80, 100, 120, 140, 160]
var_levels = ["low", "low-medium", "medium", "medium-high", "high"]
n_candidates = 150                           # Candidate calibration sets per shape and sample size.
sets_per_level = n_candidates // len(var_levels)
predictors = ["Point_density", "Vegetation_density", "Max_height",
              "Mean_height", "P25", "P50", "P75", "P95"]
y_col = "AGB(Mg/ha)_Mean"                    # Mean of the methods J, Y, C and W.
cpi_weights = {"R2": 0.15, "RMSE_cal": 0.10, "RMSE_val": 0.35, "bias": 0.25, "p": 0.15}

model_names = ["RF", "XGBoost", "SVR", "Cubist", "KNN", "GPR"]
model_colors = {"RF": "#000000", "XGBoost": "#E41A1C", "SVR": "#FFD700",
                "Cubist": "#4CAF50", "KNN": "#87CEFA", "GPR": "#0000FF"}
shape_marker = {"circle": "o", "square": "s"}
shape_title = {"circle": "Circle", "square": "Square"}
spine_lw = 0.5

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"],
    "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 11,
    "figure.dpi": 500, "mathtext.fontset": "stix",
})

#------------------------------------------------------------------------------
# Load the prepared dataset and split each shape into a validation set and a
# calibration pool.
print("Loading Prepared_Dataset.csv ...")
df = pd.read_csv(data_file)
df = df[df["Plot_size_m2"] == plot_size]
print(f"  {len(df)} pseudo-plots at {plot_size} m2")

val_pool, cal_pool = {}, {}
rng_split = np.random.default_rng(seed)
for shape in plot_shapes:
    sub = df[df["Shape"] == shape].reset_index(drop=True)
    idx = np.arange(len(sub))
    rng_split.shuffle(idx)
    val_pool[shape] = sub.iloc[sorted(idx[:n_val])].reset_index(drop=True)
    cal_pool[shape] = sub.iloc[sorted(idx[n_val:])].reset_index(drop=True)
    print(f"  {shape}: validation {len(val_pool[shape])}, calibration pool {len(cal_pool[shape])}")

#------------------------------------------------------------------------------
# Models and metrics.
def make_models(rs):
    # The six algorithms with hyperparameters fixed a priori (Table S1). XGBoost
    # is kept shallow and GPR carries a noise term so neither interpolates the
    # calibration data. SVR standardises the response inside each fit because
    # C and epsilon are on the response scale.
    return {
        "RF": RandomForestRegressor(n_estimators=100, min_samples_leaf=3,
                                    max_features="sqrt", random_state=rs, n_jobs=1),
        "XGBoost": XGBRegressor(n_estimators=50, max_depth=3, learning_rate=0.1,
                                random_state=rs, n_jobs=1, verbosity=0),
        "SVR": TransformedTargetRegressor(
            regressor=Pipeline([("sc", StandardScaler()),
                                ("m", SVR(kernel="rbf", C=1.0, epsilon=0.1))]),
            transformer=StandardScaler()),
        "Cubist": Cubist(n_committees=1, n_rules=100, random_state=rs),
        "KNN": Pipeline([("sc", StandardScaler()),
                         ("m", KNeighborsRegressor(n_neighbors=5))]),
        "GPR": Pipeline([("sc", StandardScaler()),
                         ("m", GaussianProcessRegressor(
                             kernel=ConstantKernel(1.0) * RBF(length_scale=1.0)
                             + WhiteKernel(noise_level=0.5),
                             random_state=rs, normalize_y=True))]),
    }

def cal_metrics(y, pred):
    # Calibration R2 (negative values kept) and RMSE.
    return r2_score(y, pred), np.sqrt(mean_squared_error(y, pred))

def val_metrics(y, pred):
    # Validation RMSE, mean bias and paired t-test p-value.
    rmse = np.sqrt(mean_squared_error(y, pred))
    bias = float(np.mean(pred - y))
    _, p = ttest_rel(pred, y)
    return rmse, bias, float(p)

def evaluate_models(x_c, y_c, x_v, y_v, rs):
    # Fit the six models on one calibration set. Calibration metrics come from
    # 5-fold cross-validation, validation metrics from the held-out plots.
    # Returns (model, R2, RMSE, val RMSE, val bias, val p-value) per model.
    cv = KFold(n_splits=5, shuffle=True, random_state=rs)
    out = []
    for name, m in make_models(rs).items():
        r2, rmse = cal_metrics(y_c, cross_val_predict(m, x_c, y_c, cv=cv))
        m.fit(x_c, y_c)
        vrmse, vbias, vp = val_metrics(y_v, m.predict(x_v))
        out.append((name, r2, rmse, vrmse, vbias, vp))
    return out

def pools(shape):
    # Predictor and response arrays of the calibration pool and validation set.
    return (cal_pool[shape][predictors].values, cal_pool[shape][y_col].values,
            val_pool[shape][predictors].values, val_pool[shape][y_col].values)

#------------------------------------------------------------------------------
# Analysis 1: sample size effect, random calibration sets.
print("\n" + "-" * 70)
print("Analysis 1 - Sample size effect")
print("-" * 70)

rng = np.random.default_rng(seed)
ss_cal, ss_val = [], []
for shape in plot_shapes:
    x_pool, y_pool, x_val, y_val = pools(shape)
    for n in sample_sizes:
        for it in range(n_iter):
            sel = rng.choice(len(y_pool), size=n, replace=False)
            sample_var = float(np.var(y_pool[sel]))
            for name, r2, rmse, vrmse, vbias, vp in evaluate_models(
                    x_pool[sel], y_pool[sel], x_val, y_val, seed + it):
                ss_cal.append({"shape": shape, "n_cal": n, "model": name, "iter": it,
                               "sample_var": sample_var, "R2": r2, "RMSE": rmse})
                ss_val.append({"shape": shape, "n_cal": n, "model": name, "iter": it,
                               "RMSE": vrmse, "bias": vbias, "abs_bias": abs(vbias),
                               "p_value": vp})
        print(f"  {shape}/n {n}: {n_iter} iterations x {len(model_names)} models")

ss_cal = pd.DataFrame(ss_cal)
ss_val = pd.DataFrame(ss_val)
ss_cal.to_csv(os.path.join(res_dir, "Objective4_SampleSize_Calibration.csv"), index=False)
ss_val.to_csv(os.path.join(res_dir, "Objective4_SampleSize_Validation.csv"), index=False)

#------------------------------------------------------------------------------
# Figure helpers.
def style_axes(ax):
    for s in ax.spines.values():
        s.set_linewidth(spine_lw)
    ax.tick_params(width=spine_lw)

def mean_matrix(data, value, shape):
    # Models (rows) by sample sizes (columns), mean of value over the iterations.
    return np.array([[data[(data["shape"] == shape) & (data["model"] == m)
                           & (data["n_cal"] == n)][value].mean()
                      for n in sample_sizes] for m in model_names])

def heatmap_panel(ax, mat, row_labels, col_labels, cmap, vmin, vmax):
    # Annotated colour matrix with white cell borders.
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels(col_labels, fontsize=8)
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_xticks(np.arange(mat.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(mat.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=0.4)
    ax.tick_params(which="minor", length=0)
    style_axes(ax)
    norm = (mat - vmin) / max(vmax - vmin, 1e-12)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                    color="white" if norm[i, j] > 0.55 else "black", fontsize=8.0)
    return im

#------------------------------------------------------------------------------
# Figure 8: calibration R2 and RMSE against sample size, mean and SD over the
# iterations, one line per model.
print("Generating Fig8 ...")
sparse_labels = [str(s) if i % 2 == 0 else "" for i, s in enumerate(sample_sizes)]
with plt.rc_context({"font.size": 10, "axes.labelsize": 11, "axes.titlesize": 12,
                     "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9}):
    fig, axes = plt.subplots(2, 2, figsize=(6.4, 4.4),
                             gridspec_kw={"hspace": 0.12, "wspace": 0.08})
    for row, metric in enumerate(["R2", "RMSE"]):
        for col, shape in enumerate(plot_shapes):
            ax = axes[row, col]
            sub = ss_cal[ss_cal["shape"] == shape]
            for mn in model_names:
                grp = sub[sub["model"] == mn].groupby("n_cal")[metric]
                ax.errorbar(sample_sizes, grp.mean().reindex(sample_sizes).values,
                            yerr=grp.std().reindex(sample_sizes).values,
                            color=model_colors[mn], lw=1.2, marker=shape_marker[shape],
                            markersize=5.0, markeredgecolor="black", markeredgewidth=0.4,
                            capsize=2.0, capthick=0.8, elinewidth=0.8, label=mn)
            ax.set_xticks(sample_sizes)
            ax.set_xticklabels(sparse_labels)
            style_axes(ax)
            if row == 0:
                ax.set_title(f"{shape_title[shape]} Plots")
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("Number of Plots")
            if col == 0:
                ax.set_ylabel("R²" if metric == "R2" else "RMSE (Mg ha$^{-1}$)")
            else:
                ax.tick_params(labelleft=False)
    for col in range(2):
        axes[0, col].set_ylim(0.0, 1.0)
        axes[0, col].set_yticks(np.arange(0.0, 1.01, 0.2))
        axes[1, col].set_ylim(0.0, 80.0)
        axes[1, col].set_yticks([0, 20, 40, 60, 80])
        axes[1, col].legend(loc="upper right", frameon=False, ncol=2,
                            handletextpad=0.4, columnspacing=0.8)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.97, bottom=0.10)
    fig.savefig(os.path.join(fig_dir, "Fig8_O4_SS_calibration.png"), dpi=500, bbox_inches="tight")
    plt.close(fig)
print("  Fig8_O4_SS_calibration.png")

#------------------------------------------------------------------------------
# Figure 9: validation RMSE, bias and p-value matrices, one row per shape.
print("Generating Fig9 ...")
col_labels = [str(n) for n in sample_sizes]
panels = [("RMSE", "RMSE (Mg ha$^{-1}$)", "YlGn_r", 10, 70),
          ("bias", "bias (Mg ha$^{-1}$)", "RdBu_r", -16, 16),
          ("p_value", "Paired t-test p-value", "YlGn", 0, 1)]
fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.2), gridspec_kw={"hspace": 0.08, "wspace": 0.04})
images = []
for ri, shape in enumerate(plot_shapes):
    for ci, (value, title, cmap, vmin, vmax) in enumerate(panels):
        ax = axes[ri, ci]
        im = heatmap_panel(ax, mean_matrix(ss_val, value, shape),
                           model_names if ci == 0 else [""] * len(model_names),
                           col_labels, cmap, vmin, vmax)
        if ri == 0:
            images.append(im)
            ax.set_title(title, fontsize=11)
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("Number of Plots")
    axes[ri, 0].set_ylabel(shape_title[shape], fontsize=11, labelpad=8)
fig.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.18)
# The colour bars are inset a little so neighbouring end ticks do not touch.
pad = 0.018
for ci, (value, title, cmap, vmin, vmax) in enumerate(panels):
    bb = axes[1, ci].get_position()
    cax = fig.add_axes([bb.x0 + pad, 0.07, bb.width - 2 * pad, 0.022])
    label = "p-value" if value == "p_value" else title
    cb = fig.colorbar(images[ci], cax=cax, orientation="horizontal", label=label)
    if value == "RMSE":
        cb.set_ticks(np.arange(vmin, vmax + 1e-9, 10))
    if value == "bias":
        cb.set_ticks(np.arange(vmin, vmax + 1e-9, 8))
    cb.outline.set_linewidth(spine_lw)
    cb.ax.tick_params(width=spine_lw)
fig.savefig(os.path.join(fig_dir, "Fig9_O4_SS_validation.png"), dpi=500, bbox_inches="tight")
plt.close(fig)
print("  Fig9_O4_SS_validation.png")

#------------------------------------------------------------------------------
# Analysis 2: variance sensitivity. Reused from the result files when present.
vs_cal_path = os.path.join(res_dir, "Objective4_Calibration_Results.csv")
vs_val_path = os.path.join(res_dir, "Objective4_Validation_Results.csv")
if all(os.path.exists(p) and os.path.getsize(p) > 0 for p in (vs_cal_path, vs_val_path)):
    vs_cal = pd.read_csv(vs_cal_path)
    vs_val = pd.read_csv(vs_val_path)
    print(f"Loaded existing variance-sensitivity results: "
          f"{len(vs_cal)} calibration rows, {len(vs_val)} validation rows")
else:
    print("\n" + "-" * 70)
    print("Analysis 2 - Variance sensitivity")
    print("-" * 70)
    rng = np.random.default_rng(seed)
    vs_cal, vs_val = [], []
    for shape in plot_shapes:
        x_pool, y_pool, x_val, y_val = pools(shape)
        for n in sample_sizes:
            # Candidate sets ranked by AGB variance and split into five groups.
            candidates = []
            for _ in range(n_candidates):
                idx = rng.choice(len(y_pool), size=n, replace=False)
                candidates.append((idx, float(np.var(y_pool[idx]))))
            candidates.sort(key=lambda t: t[1])
            groups = {vl: [c[0] for c in candidates[i * sets_per_level:(i + 1) * sets_per_level]]
                      for i, vl in enumerate(var_levels)}
            for vl in var_levels:
                for it, sel in enumerate(groups[vl]):
                    for name, r2, rmse, vrmse, vbias, vp in evaluate_models(
                            x_pool[sel], y_pool[sel], x_val, y_val, seed + it):
                        vs_cal.append({"shape": shape, "n_cal": n, "variance": vl,
                                       "model": name, "iter": it, "R2": r2, "RMSE": rmse})
                        vs_val.append({"shape": shape, "n_cal": n, "variance": vl,
                                       "model": name, "iter": it, "RMSE": vrmse,
                                       "bias": vbias, "p_value": vp})
                print(f"  {shape}/n {n}/{vl}: {len(groups[vl])} iterations done")
    vs_cal = pd.DataFrame(vs_cal)
    vs_val = pd.DataFrame(vs_val)
    vs_cal.to_csv(vs_cal_path, index=False)
    vs_val.to_csv(vs_val_path, index=False)

#------------------------------------------------------------------------------
# CPI, VSI and VRI. Every metric is scaled to 0 to 1 over all rows (larger is
# better), the p-value scores 1 once it is above 0.05, and the CPI is their
# weighted sum. VSI is the coefficient of variation of the mean CPI across the
# five variance levels and VRI its min-max complement, so higher is more robust.
print("Computing CPI, VSI and VRI ...")
key = ["shape", "n_cal", "variance", "model", "iter"]
merged = vs_cal.rename(columns={"R2": "R2_cal", "RMSE": "RMSE_cal"}).merge(
    vs_val.rename(columns={"RMSE": "RMSE_val"}), on=key)
merged["abs_bias"] = merged["bias"].abs()

def scale_up(s):
    # 0 to 1, larger is better.
    lo, hi = s.min(), s.max()
    return (s - lo) / (hi - lo) if hi - lo >= 1e-12 else pd.Series(0.0, index=s.index)

def scale_down(s):
    # 0 to 1, smaller is better.
    lo, hi = s.min(), s.max()
    return 1.0 - (s - lo) / (hi - lo) if hi - lo >= 1e-12 else pd.Series(1.0, index=s.index)

merged["S_R2"] = scale_up(merged["R2_cal"])
merged["S_RMSEcal"] = scale_down(merged["RMSE_cal"])
merged["S_RMSEval"] = scale_down(merged["RMSE_val"])
merged["S_bias"] = scale_down(merged["abs_bias"])
merged["S_p"] = (merged["p_value"] / 0.05).clip(upper=1.0)
merged["CPI"] = (cpi_weights["R2"] * merged["S_R2"]
                 + cpi_weights["RMSE_cal"] * merged["S_RMSEcal"]
                 + cpi_weights["RMSE_val"] * merged["S_RMSEval"]
                 + cpi_weights["bias"] * merged["S_bias"]
                 + cpi_weights["p"] * merged["S_p"])

cpi_mean = (merged.groupby(["shape", "n_cal", "variance", "model"])["CPI"].mean()
            .reset_index().rename(columns={"CPI": "CPI_mean"}))
by_level = cpi_mean.pivot(index=["shape", "n_cal", "model"], columns="variance",
                          values="CPI_mean")[var_levels]
vsi = by_level.index.to_frame(index=False)
by_level = by_level.to_numpy()
vsi["VSI"] = by_level.std(axis=1) / by_level.mean(axis=1)
vsi["CPI_avg"] = by_level.mean(axis=1)
vsi["VRI"] = 1.0 - (vsi["VSI"] - vsi["VSI"].min()) / max(vsi["VSI"].max() - vsi["VSI"].min(), 1e-12)

merged.to_csv(os.path.join(res_dir, "Objective4_CPI_per_iteration.csv"), index=False)
cpi_mean.to_csv(os.path.join(res_dir, "Objective4_CPI_mean.csv"), index=False)
vsi.to_csv(os.path.join(res_dir, "Objective4_VSI_VRI.csv"), index=False)
print(f"  wrote CPI per iteration ({len(merged)} rows), CPI mean ({len(cpi_mean)} rows), "
      f"VSI and VRI ({len(vsi)} rows)")

#------------------------------------------------------------------------------
# Figure 10: VRI against sample size, one line per model.
print("Generating Fig10 ...")
fig, axes = plt.subplots(1, 2, figsize=(8, 4.2), sharey=True, gridspec_kw={"wspace": 0.05})
for ax, shape in zip(axes, plot_shapes):
    for mn in model_names:
        sub = vsi[(vsi["shape"] == shape) & (vsi["model"] == mn)].sort_values("n_cal")
        ax.plot(sub["n_cal"], sub["VRI"], "-", color=model_colors[mn], lw=1.2,
                marker=shape_marker[shape], markersize=5.5, markeredgecolor="black",
                markeredgewidth=0.4, label=mn)
    ax.set_xlabel("Number of Plots")
    ax.set_xticks(sample_sizes)
    ax.set_ylim(0, 1.02)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_title(f"{shape_title[shape]} Plots", fontsize=11)
    ax.set_box_aspect(1)
    ax.legend(fontsize=8, loc="lower right", frameon=False, ncol=2)
    style_axes(ax)
axes[0].set_ylabel("Variance Robustness Index (VRI)")
fig.tight_layout(rect=[0, 0, 1, 1], w_pad=0.3)
fig.savefig(os.path.join(fig_dir, "Fig10_O4_indices.png"), dpi=500, bbox_inches="tight")
plt.close(fig)
print("  Fig10_O4_indices.png")

#------------------------------------------------------------------------------
# Figure S2: mean CPI of every model and variance level, one panel per shape
# and sample size.
print("Generating FigS2 ...")
level_labels = ["L", "LM", "M", "MH", "H"]
fig = plt.figure(figsize=(14, 7.0))
gs = gridspec.GridSpec(2, 7, left=0.05, right=0.94, top=0.92, bottom=0.06,
                       hspace=0.05, wspace=0.08)
for ri, shape in enumerate(plot_shapes):
    for ci, n in enumerate(sample_sizes):
        ax = fig.add_subplot(gs[ri, ci])
        sub = cpi_mean[(cpi_mean["shape"] == shape) & (cpi_mean["n_cal"] == n)]
        mat = np.full((len(model_names), len(var_levels)), np.nan)
        for mi, mn in enumerate(model_names):
            row = sub[sub["model"] == mn].set_index("variance")
            for vi, vl in enumerate(var_levels):
                if vl in row.index:
                    mat[mi, vi] = row.loc[vl, "CPI_mean"]
        im = ax.imshow(mat, aspect="auto", cmap="viridis", vmin=0, vmax=1)
        ax.set_xticks(range(len(var_levels)))
        ax.set_xticklabels(level_labels if ri == 1 else [""] * len(var_levels), fontsize=8)
        ax.set_yticks(range(len(model_names)))
        ax.set_yticklabels(model_names if ci == 0 else [""] * len(model_names), fontsize=8)
        style_axes(ax)
        ax.tick_params(length=2)
        if ri == 0:
            ax.set_title(f"n = {n}", fontsize=9)
        if ci == 0:
            ax.set_ylabel(shape_title[shape], fontsize=10, labelpad=8)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if not np.isnan(mat[i, j]):
                    ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                            color="white" if mat[i, j] < 0.55 else "black", fontsize=7)
cax = fig.add_axes([0.95, 0.20, 0.012, 0.55])
cb = fig.colorbar(im, cax=cax)
cb.set_label("CPI", fontsize=9)
cb.outline.set_linewidth(spine_lw)
cb.ax.tick_params(width=spine_lw, labelsize=8)
fig.savefig(os.path.join(fig_dir, "FigS2_O4_CPI_heatmap.png"), dpi=500, bbox_inches="tight")
plt.close(fig)
print("  FigS2_O4_CPI_heatmap.png")

print("\n" + "-" * 70)
print(f"Analysis complete. Results: {res_dir}; figures: {fig_dir}")
print("-" * 70)

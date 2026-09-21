# RF_Model_Analysis.py
# Author: Xinyuan Wei.
# Sep. 02 2026

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import ttest_rel
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.model_selection import KFold, GroupKFold

matplotlib.use("Agg")
os.environ.setdefault("PYTHONWARNINGS", "ignore")   # Inherited by the joblib workers.
warnings.filterwarnings("ignore")

#------------------------------------------------------------------------------
# Paths and parameters.
data_dir = os.path.dirname(os.path.abspath(__file__))
res_dir = os.path.join(data_dir, "Results")
fig_dir = os.path.abspath(os.path.join(data_dir, "..", "2.Final Paper", "1.Figure"))
os.makedirs(res_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

seed = 42
n_trees = 100
min_leaf = 3
plot_shapes = ["circle", "square"]
predictors = ["Point_density", "Vegetation_density", "Max_height",
              "Mean_height", "P25", "P50", "P75", "P95"]
allo_keys = ["Y", "J", "C", "W"]             # As keyed in the inventory workbook.
allo_display = {"Y": "Young et al. 1980", "J": "Jenkins et al. 2003",
                "C": "Chojnacky et al. 2014", "W": "Westfall et al. 2024"}
y_mean = "AGB(Mg/ha)_Mean"                   # Response of Objectives 2 and 3 and Figure 11.

n_cal = 120                                  # Calibration plots, Objectives 1 and 2.
n_val = 30                                   # Independent validation plots.
plot_sizes = [100, 150, 200, 225, 250, 300]
grid_sizes = [225, 100]                      # Application grids of Objective 2.
sample_sizes = [20, 40, 60, 80, 100, 120, 140, 160]
var_levels = ["low", "low-medium", "medium", "medium-high", "high"]
n_candidates = 150                           # Candidate calibration sets per sample size.
sets_per_level = n_candidates // len(var_levels)
n_reps_cv = 30                               # Repeats of random and grouped CV, Figure 11.
n_boot = 100                                 # One-pseudo-plot-per-parent draws, Figure 11.

shape_marker = {"circle": "o", "square": "s"}
shape_title = {s: f"{s.title()} Plots" for s in plot_shapes}
series_colors = ["#1f77b4", "#d95f02"]
var_colors = {"low": "#2166AC", "low-medium": "#87CEFA", "medium": "#4CAF50",
              "medium-high": "#FFD700", "high": "#E31A1C"}
var_legend = [Patch(facecolor=var_colors[vl], alpha=0.7,
                    label=vl[0].upper() + vl[1:]) for vl in var_levels]
box_lw = 0.5
spine_lw = 0.5
box_width = 0.13
box_offsets = np.array([-0.32, -0.16, 0.0, 0.16, 0.32])

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"],
    "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 11,
    "figure.dpi": 500, "mathtext.fontset": "stix",
})

#------------------------------------------------------------------------------
# Load the prepared dataset.
print("Loading prepared dataset ...")
df = pd.read_csv(os.path.join(data_dir, "Prepared_Dataset.csv"))
n_parents = df["parent_plot"].nunique()
sizes_found = sorted(df["Plot_size_m2"].unique())
print(f"  {len(df)} pseudo-plots, {n_parents} parent plots, sizes {sizes_found}")

#------------------------------------------------------------------------------
# Helpers.
def fit_rf(x, y, rs, oob=False):
    # Random Forest with fixed settings; oob=True keeps the out-of-bag predictions.
    rf = RandomForestRegressor(n_estimators=n_trees, min_samples_leaf=min_leaf,
                               max_features="sqrt", random_state=rs, n_jobs=-1,
                               oob_score=oob)
    rf.fit(x, y)
    return rf

def cv_predict(x, y, rs):
    # Out-of-fold predictions from random 5-fold cross-validation.
    kf = KFold(n_splits=5, shuffle=True, random_state=rs)
    pred = np.full(len(y), np.nan)
    for k, (tr, te) in enumerate(kf.split(x)):
        pred[te] = fit_rf(x[tr], y[tr], rs + k).predict(x[te])
    return pred

def cal_metrics(y, pred):
    # Calibration R2 (negative values kept) and RMSE.
    return r2_score(y, pred), np.sqrt(mean_squared_error(y, pred))

def val_metrics(y, pred):
    # Validation RMSE, mean bias and paired t-test p-value.
    rmse = np.sqrt(mean_squared_error(y, pred))
    bias = float(np.mean(pred - y))
    _, p = ttest_rel(pred, y)
    return rmse, bias, float(p)

def even_agb_split(data, y_col, n):
    # n calibration plots at equally spaced positions of the AGB-sorted data, and the rest.
    order = data.sort_values(y_col).index
    pos = np.unique(np.linspace(0, len(order) - 1, n).round().astype(int))
    if len(pos) < n:
        pos = np.concatenate([pos, np.setdiff1d(np.arange(len(order)), pos)[:n - len(pos)]])
    cal_idx = order[pos]
    return data.loc[cal_idx].copy(), data.drop(cal_idx).copy()

def random_subset(data, n, rs):
    # n rows drawn at random without replacement.
    rng = np.random.default_rng(rs)
    return data.loc[rng.choice(data.index.to_numpy(), size=n, replace=False)].copy()

def ann_cal(obs, pred, tag=None):
    # Panel annotation of a calibration scatter.
    r2, rmse = cal_metrics(obs, pred)
    return f"R² = {r2:.3f}\nRMSE = {rmse:.2f}"

def ann_val(obs, pred, tag=None):
    # Panel annotation of a validation scatter; tag is written as a subscript.
    rmse, bias, p = val_metrics(obs, pred)
    sub = f"$_{{{tag}}}$" if tag else ""
    ptxt = f"p-value{sub} < 0.001" if p < 0.001 else f"p-value{sub} = {p:.3f}"
    return f"RMSE{sub} = {rmse:.2f}\nbias{sub} = {bias:.2f}\n{ptxt}"

def scatter_grid(panels, row_keys, col_keys, col_labels, annotate, out_name,
                 lim, tick, font_scale=1.0):
    # Grid of observed against predicted panels. panels[(row, col)] is a list
    # of (obs, pred, tag) series: one series is drawn in blue, two series in
    # blue and orange with their annotations in opposite corners.
    nrows, ncols = len(row_keys), len(col_keys)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(2.1 * ncols + 0.7, 2.1 * nrows + 0.65),
                             sharex=True, sharey=True, squeeze=False)
    hidden = 0
    for r, rk in enumerate(row_keys):
        for c, ck in enumerate(col_keys):
            ax = axes[r, c]
            series = panels[(rk, ck)]
            for k, (obs, pred, tag) in enumerate(series):
                color = "blue" if len(series) == 1 else series_colors[k]
                hidden += int(np.sum((np.asarray(obs) > lim) | (np.asarray(pred) > lim)))
                ax.scatter(obs, pred, s=20, alpha=0.3, color=color,
                           edgecolor="black", linewidth=0.3, marker=shape_marker[rk])
                if k == 0:
                    ax.text(0.05, 0.95, annotate(obs, pred, tag), transform=ax.transAxes,
                            fontsize=8 * font_scale, va="top", ha="left",
                            color="black" if len(series) == 1 else color)
                else:
                    ax.text(0.95, 0.05, annotate(obs, pred, tag), transform=ax.transAxes,
                            fontsize=8 * font_scale, va="bottom", ha="right", color=color)
            ax.plot([0, lim], [0, lim], "r--", lw=1)
            ax.set_xlim(0, lim)
            ax.set_ylim(0, lim)
            ax.set_aspect("equal")
            ax.set_xticks(np.arange(0, lim + 0.1, tick))
            ax.set_yticks(np.arange(0, lim + 0.1, tick))
            if r == 0:
                ax.set_title(col_labels[c], fontsize=10 * font_scale)
            if c > 0:
                ax.tick_params(labelleft=False)
            if r < nrows - 1:
                ax.tick_params(labelbottom=False)
    # Panels touch, so the leading "0" of every panel after the first is
    # blanked; otherwise it runs into the last label of the panel before it.
    labels = [f"{v:g}" for v in np.arange(0, lim + 0.1, tick)]
    labels[0] = ""
    for c in range(1, ncols):
        axes[nrows - 1, c].set_xticklabels(labels)
    for r in range(nrows):
        axes[r, 0].set_ylabel("Predicted AGB (Mg ha$^{-1}$)")
    for c in range(ncols):
        axes[nrows - 1, c].set_xlabel("Observed AGB (Mg ha$^{-1}$)")
    fig.tight_layout(rect=[0.07, 0, 1, 1])
    fig.subplots_adjust(wspace=0.05, hspace=0.08)
    for r, rk in enumerate(row_keys):
        axes[r, 0].annotate(shape_title[rk], xy=(-0.34, 0.5), xycoords="axes fraction",
                            fontsize=10 * font_scale, rotation=90, ha="center", va="center")
    if hidden:
        print(f"    Note {out_name}: {hidden} points fall outside the 0 to {lim:g} axis")
    fig.savefig(os.path.join(fig_dir, out_name), dpi=500, bbox_inches="tight")
    plt.close(fig)
    print(f"  {out_name}")

def styled_boxplot(ax, data, positions, color):
    # Filled boxes without fliers, thin black lines.
    bp = ax.boxplot(data, positions=positions, widths=box_width, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=box_lw),
                    boxprops=dict(lw=box_lw), whiskerprops=dict(lw=box_lw),
                    capprops=dict(lw=box_lw))
    for patch in bp["boxes"]:
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

def style_axes(ax):
    for s in ax.spines.values():
        s.set_linewidth(spine_lw)
    ax.tick_params(width=spine_lw)

def variance_boxplots(raw, metrics, ylabels, ylims, yticks, out_name, zero_line=False):
    # Two metric rows by two shape columns; at every sample size one box per
    # variance group, summarising the replicate iterations in raw.
    with plt.rc_context({"font.size": 8, "axes.labelsize": 9, "axes.titlesize": 10,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8}):
        fig, axes = plt.subplots(2, 2, figsize=(5.6, 3.6))
        xlabels = [str(s) if i % 2 == 0 else "" for i, s in enumerate(sample_sizes)]
        for col, shape in enumerate(plot_shapes):
            sub = raw[raw["shape"] == shape]
            for vi, vl in enumerate(var_levels):
                vs = sub[sub["variance"] == vl]
                positions = np.arange(len(sample_sizes)) + box_offsets[vi]
                for row, metric in enumerate(metrics):
                    data = [vs[vs["n_cal"] == n][metric].values for n in sample_sizes]
                    styled_boxplot(axes[row, col], data, positions, var_colors[vl])
            for row in (0, 1):
                ax = axes[row, col]
                ax.set_xticks(range(len(sample_sizes)))
                ax.set_xticklabels(xlabels)
                ax.set_ylim(*ylims[row])
                ax.set_yticks(yticks[row])
                style_axes(ax)
                if col == 0:
                    ax.set_ylabel(ylabels[row])
                else:
                    ax.set_yticklabels([])
            axes[0, col].tick_params(labelbottom=False)
            axes[0, col].set_title(shape_title[shape])
            axes[1, col].set_xlabel("Number of Plots")
            if zero_line:
                axes[1, col].axhline(0, color="gray", ls="--", lw=spine_lw)
        fig.subplots_adjust(left=0.12, right=0.98, top=0.96, bottom=0.22,
                            wspace=0.10, hspace=0.18)
        bb0 = axes[1, 0].get_position()
        bb1 = axes[1, 1].get_position()
        fig.legend(handles=var_legend, ncol=5, loc="lower center",
                   bbox_to_anchor=((bb0.x0 + bb1.x1) / 2, 0.01), frameon=False,
                   handlelength=1.4, handletextpad=0.4, columnspacing=1.0)
        fig.savefig(os.path.join(fig_dir, out_name), dpi=500, bbox_inches="tight")
        plt.close(fig)
    print(f"  {out_name}")

#------------------------------------------------------------------------------
# Objective 1: allometric method comparison at 300 m2.
print("\n" + "-" * 70)
print("Objective 1 - Allometric method comparison")
print("-" * 70)

o1_cal, o1_val = [], []
o1_cal_panels, o1_val_panels = {}, {}
for shape in plot_shapes:
    d300 = df[(df["Shape"] == shape) & (df["Plot_size_m2"] == 300)]
    for allo in allo_keys:
        y_col = f"AGB(Mg/ha)_{allo}"
        cal, pool = even_agb_split(d300, y_col, n_cal)
        val = random_subset(pool, n_val, seed)
        y_cal, y_val = cal[y_col].values, val[y_col].values
        # OOB predictions give the calibration metrics; the same model predicts the validation plots.
        rf = fit_rf(cal[predictors].values, y_cal, seed, oob=True)
        r2, rmse = cal_metrics(y_cal, rf.oob_prediction_)
        pred = rf.predict(val[predictors].values)
        vrmse, vbias, vp = val_metrics(y_val, pred)
        o1_cal.append({"shape": shape, "allometric": allo_display[allo],
                       "R2": r2, "RMSE": rmse})
        o1_val.append({"shape": shape, "allometric": allo_display[allo],
                       "RMSE": vrmse, "bias": vbias, "p_value": vp})
        o1_cal_panels[(shape, allo)] = [(y_cal, rf.oob_prediction_, None)]
        o1_val_panels[(shape, allo)] = [(y_val, pred, None)]
        print(f"  {shape}/{allo_display[allo]}: cal R2 {r2:.3f}, RMSE {rmse:.2f} | "
              f"val RMSE {vrmse:.2f}, bias {vbias:.2f}, p {vp:.3f}")

pd.DataFrame(o1_cal).to_csv(os.path.join(res_dir, "Objective1_Calibration_Results.csv"), index=False)
pd.DataFrame(o1_val).to_csv(os.path.join(res_dir, "Objective1_Validation_Results.csv"), index=False)

allo_labels = [allo_display[a] for a in allo_keys]
scatter_grid(o1_cal_panels, plot_shapes, allo_keys, allo_labels, ann_cal,
             "Fig2_O1_calibration.png", lim=300, tick=100)
scatter_grid(o1_val_panels, plot_shapes, allo_keys, allo_labels, ann_val,
             "Fig3_O1_validation.png", lim=300, tick=100)

#------------------------------------------------------------------------------
# Objective 2: plot size and shape effect. One model per shape and size, all
# applied to the same square 225 m2 and 100 m2 validation sets.
print("\n" + "-" * 70)
print("Objective 2 - Plot size and shape effect")
print("-" * 70)

rfs, cal_index = {}, {}
o2_cal, o2_cal_panels = [], {}
for shape in plot_shapes:
    for area in plot_sizes:
        sub = df[(df["Shape"] == shape) & (df["Plot_size_m2"] == area)]
        cal, _ = even_agb_split(sub, y_mean, n_cal)
        cal_index[(shape, area)] = cal.index
        y = cal[y_mean].values
        rf = fit_rf(cal[predictors].values, y, seed, oob=True)
        rfs[(shape, area)] = rf
        r2, rmse = cal_metrics(y, rf.oob_prediction_)
        o2_cal.append({"shape": shape, "area_m2": area, "R2": r2, "RMSE": rmse,
                       "n_train": len(y)})
        o2_cal_panels[(shape, area)] = [(y, rf.oob_prediction_, None)]
        print(f"  {shape}/{area} m2: cal R2 {r2:.3f}, RMSE {rmse:.2f} (n {len(y)})")

# Validation sets: 30 random square plots per grid size, excluding the plots
# used to calibrate the square model of that size.
val_sets = {}
for area in grid_sizes:
    pool = df[(df["Shape"] == "square") & (df["Plot_size_m2"] == area)]
    pool = pool.loc[~pool.index.isin(cal_index[("square", area)])]
    vs = random_subset(pool, n_val, seed + area)
    val_sets[area] = (vs[predictors].values, vs[y_mean].values)

o2_app, o2_app_panels = [], {}
for shape in plot_shapes:
    for area in plot_sizes:
        row = {"shape": shape, "area_m2": area}
        panel, note = [], []
        for grid in grid_sizes:
            x_v, y_v = val_sets[grid]
            pred = rfs[(shape, area)].predict(x_v)
            rmse, bias, p = val_metrics(y_v, pred)
            row.update({f"RMSE_{grid}": rmse, f"bias_{grid}": bias, f"p_value_{grid}": p})
            panel.append((y_v, pred, str(grid)))
            note.append(f"{grid} m2 grid RMSE {rmse:.2f}, bias {bias:.2f}")
        row.update({f"n_val_{grid}": len(val_sets[grid][1]) for grid in grid_sizes})
        o2_app.append(row)
        o2_app_panels[(shape, area)] = panel
        print(f"    {shape}/{area} m2 | " + " | ".join(note))

pd.DataFrame(o2_cal).to_csv(os.path.join(res_dir, "Objective2_Calibration_Results.csv"), index=False)
pd.DataFrame(o2_app).to_csv(os.path.join(res_dir, "Objective2_Application_Results.csv"), index=False)

size_labels = [f"{a} m²" for a in plot_sizes]
scatter_grid(o2_cal_panels, plot_shapes, plot_sizes, size_labels, ann_cal,
             "Fig4_O2_calibration.png", lim=400, tick=100, font_scale=1.2)
scatter_grid(o2_app_panels, plot_shapes, plot_sizes, size_labels, ann_val,
             "Fig5_O2_application.png", lim=400, tick=100, font_scale=1.2)

#------------------------------------------------------------------------------
# Objective 3: sample size and AGB variance effect at 300 m2. For every sample
# size, 150 random calibration sets are ranked by AGB variance and split into
# five groups of 30; each set is one iteration, validated on the plots left out.
print("\n" + "-" * 70)
print("Objective 3 - Sample size and AGB variance")
print("-" * 70)

rng = np.random.default_rng(seed)            # One stream across all combinations.
o3_cal, o3_val, o3_rows = [], [], []
for shape in plot_shapes:
    pool = df[(df["Shape"] == shape) & (df["Plot_size_m2"] == 300)].reset_index(drop=True)
    n_total = len(pool)
    x_pool, y_pool = pool[predictors].values, pool[y_mean].values

    for n_sub in sample_sizes:
        candidates = []
        for _ in range(n_candidates):
            idx = rng.choice(n_total, size=n_sub, replace=False)
            candidates.append((idx, float(np.var(y_pool[idx]))))
        candidates.sort(key=lambda t: t[1])
        groups = {vl: [c[0] for c in candidates[i * sets_per_level:(i + 1) * sets_per_level]]
                  for i, vl in enumerate(var_levels)}

        for vl in var_levels:
            cal_r2, cal_rmse, val_rmse, val_bias, val_p = [], [], [], [], []
            for it, cal_idx in enumerate(groups[vl]):
                val_idx = np.setdiff1d(np.arange(n_total), cal_idx, assume_unique=True)
                rf = fit_rf(x_pool[cal_idx], y_pool[cal_idx], seed + it, oob=True)
                r2, rmse = cal_metrics(y_pool[cal_idx], rf.oob_prediction_)
                vrmse, vbias, vp = val_metrics(y_pool[val_idx], rf.predict(x_pool[val_idx]))
                cal_r2.append(r2)
                cal_rmse.append(rmse)
                val_rmse.append(vrmse)
                val_bias.append(vbias)
                val_p.append(vp)
                o3_rows.append({"shape": shape, "n_cal": n_sub, "variance": vl,
                                "cal_R2": r2, "cal_RMSE": rmse,
                                "val_RMSE": vrmse, "val_bias": vbias})
            o3_cal.append({"shape": shape, "n_cal": n_sub, "variance": vl,
                           "R2_mean": np.mean(cal_r2), "R2_std": np.std(cal_r2),
                           "RMSE_mean": np.mean(cal_rmse), "RMSE_std": np.std(cal_rmse),
                           "n_iter": len(cal_r2)})
            o3_val.append({"shape": shape, "n_cal": n_sub, "variance": vl,
                           "RMSE_mean": np.mean(val_rmse), "RMSE_std": np.std(val_rmse),
                           "bias_mean": np.mean(val_bias), "bias_std": np.std(val_bias),
                           "p_value_mean": np.mean(val_p), "p_value_std": np.std(val_p),
                           "n_iter": len(val_rmse)})
            print(f"  {shape}/n {n_sub}/{vl}: cal R2 {np.mean(cal_r2):.3f} "
                  f"(sd {np.std(cal_r2):.3f}), RMSE {np.mean(cal_rmse):.2f} | "
                  f"val RMSE {np.mean(val_rmse):.2f}, bias {np.mean(val_bias):.2f}")

pd.DataFrame(o3_cal).to_csv(os.path.join(res_dir, "Objective3_Calibration_Results.csv"), index=False)
pd.DataFrame(o3_val).to_csv(os.path.join(res_dir, "Objective3_Validation_Results.csv"), index=False)
o3 = pd.DataFrame(o3_rows)

variance_boxplots(o3, ["cal_R2", "cal_RMSE"], ["R²", "RMSE (Mg ha$^{-1}$)"],
                  [(-0.5, 1), (0, 80)],
                  [np.arange(-0.5, 1.01, 0.5), [0, 20, 40, 60, 80]],
                  "Fig6_O3_calibration.png")
variance_boxplots(o3, ["val_RMSE", "val_bias"], ["RMSE (Mg ha$^{-1}$)", "bias (Mg ha$^{-1}$)"],
                  [(0, 80), (-40, 40)],
                  [[0, 20, 40, 60, 80], [-40, -20, 0, 20, 40]],
                  "Fig7_O3_validation.png", zero_line=True)

#------------------------------------------------------------------------------
# Figure 11: hierarchical structure and pseudoreplication diagnostics.
# (a, b) Intraclass correlation (ICC) of AGB and of the ALS metrics across
#        parent plots by pseudo-plot size: near 0 means the pseudo-plots of a
#        parent sample different forest, near 1 means they are near-replicates.
# (c, d) Calibration R2 and RMSE of the 300 m2 pool under random 5-fold CV,
#        grouped 5-fold CV (all pseudo-plots of a parent held out together)
#        and a strictly independent design with one pseudo-plot per parent.
print("\n" + "-" * 70)
print("Figure 11 - Spatial autocorrelation and pseudoreplication diagnostics")
print("-" * 70)

def icc(g, col):
    # One-way ANOVA ICC: between-parent variance over total variance.
    n_groups = g["parent_plot"].nunique()
    grp_means = g.groupby("parent_plot")[col].mean()
    grp_sizes = g.groupby("parent_plot")[col].count()
    n_per = grp_sizes.mean()
    mu = g[col].mean()
    msb = ((grp_means - mu) ** 2 * grp_sizes).sum() / (n_groups - 1)
    msw = g.groupby("parent_plot")[col].apply(
        lambda x: ((x - x.mean()) ** 2).sum()).sum() / (len(g) - n_groups)
    return (msb - msw) / (msb + (n_per - 1) * msw)

icc_agb, icc_als = {}, {}
for shape in plot_shapes:
    icc_agb[shape] = [icc(df[(df["Shape"] == shape) & (df["Plot_size_m2"] == sz)], y_mean)
                      for sz in plot_sizes]
    icc_als[shape] = [np.mean([icc(df[(df["Shape"] == shape) & (df["Plot_size_m2"] == sz)], f)
                               for f in predictors]) for sz in plot_sizes]
    print(f"  {shape} ICC of AGB by size: "
          + ", ".join(f"{sz} m2 {v:.2f}" for sz, v in zip(plot_sizes, icc_agb[shape])))
    print(f"  {shape} mean ICC of ALS metrics: "
          + ", ".join(f"{sz} m2 {v:.2f}" for sz, v in zip(plot_sizes, icc_als[shape])))

# Random and grouped 5-fold CV on the 300 m2 pool, repeated with shuffled folds.
cv_r2 = {k: {s: [] for s in plot_shapes} for k in ["random", "grouped", "independent"]}
cv_rmse = {k: {s: [] for s in plot_shapes} for k in ["random", "grouped", "independent"]}
for shape in plot_shapes:
    pool = df[(df["Shape"] == shape) & (df["Plot_size_m2"] == 300)].reset_index(drop=True)
    x_pool, y_pool = pool[predictors].values, pool[y_mean].values
    g_pool = pool["parent_plot"].values
    parents = np.unique(g_pool)
    for rep in range(n_reps_cv):
        r2, rmse = cal_metrics(y_pool, cv_predict(x_pool, y_pool, seed + rep))
        cv_r2["random"][shape].append(r2)
        cv_rmse["random"][shape].append(rmse)
        # The parents are relabelled by a random permutation so the grouped folds differ per repeat.
        perm = np.random.default_rng(seed + rep).permutation(len(parents))
        g_shuf = perm[np.searchsorted(parents, g_pool)]
        pred = np.full(len(y_pool), np.nan)
        for tr, te in GroupKFold(n_splits=5).split(x_pool, y_pool, g_shuf):
            pred[te] = fit_rf(x_pool[tr], y_pool[tr], seed + rep).predict(x_pool[te])
        r2, rmse = cal_metrics(y_pool, pred)
        cv_r2["grouped"][shape].append(r2)
        cv_rmse["grouped"][shape].append(rmse)
    # One random pseudo-plot per parent plot, OOB metrics, repeated draws.
    for b in range(n_boot):
        rng_b = np.random.default_rng(seed + b)
        sel = [int(rng_b.choice(pool.index[pool["parent_plot"] == p].to_numpy(), size=1)[0])
               for p in pool["parent_plot"].unique()]
        x, y = pool.loc[sel, predictors].values, pool.loc[sel, y_mean].values
        rf = fit_rf(x, y, seed + b, oob=True)
        r2, rmse = cal_metrics(y, rf.oob_prediction_)
        cv_r2["independent"][shape].append(r2)
        cv_rmse["independent"][shape].append(rmse)
    for k in ["random", "grouped", "independent"]:
        print(f"  {shape} {k}: R2 {np.mean(cv_r2[k][shape]):.3f} "
              f"(sd {np.std(cv_r2[k][shape]):.3f}), RMSE {np.mean(cv_rmse[k][shape]):.2f}")

diag_names = {"random": "random_5fold_CV", "grouped": "grouped_5fold_CV_by_parent",
              "independent": "one_subplot_per_parent"}
summary = [{"diagnostic": diag_names[k], "shape": shape,
            "R2_mean": float(np.mean(cv_r2[k][shape])), "R2_std": float(np.std(cv_r2[k][shape])),
            "RMSE_mean": float(np.mean(cv_rmse[k][shape])), "RMSE_std": float(np.std(cv_rmse[k][shape])),
            "n_iter": len(cv_r2[k][shape])}
           for k in ["random", "grouped", "independent"] for shape in plot_shapes]
pd.DataFrame(summary).to_csv(os.path.join(res_dir, "FigureS1_autocorrelation_summary.csv"), index=False)
pd.DataFrame({"Plot_size_m2": plot_sizes,
              "ICC_AGB_circle": icc_agb["circle"], "ICC_AGB_square": icc_agb["square"],
              "ICC_LiDAR_circle": icc_als["circle"], "ICC_LiDAR_square": icc_als["square"]}
             ).to_csv(os.path.join(res_dir, "FigureS1_ICC_by_size.csv"), index=False)

# Figure 11 is the widest panel grid of the paper, so its text is set larger.
with plt.rc_context({"font.size": 12, "axes.labelsize": 12.5, "xtick.labelsize": 11.5,
                     "ytick.labelsize": 11.5, "legend.fontsize": 10.5}):
    fig, axes = plt.subplots(2, 2, figsize=(7.6, 5.6))
    for ax, values, ylabel, ylim, yticks, label in [
            (axes[0, 0], icc_agb, "ICC of AGB across parent plots", (0, 1.05), None, "(a)"),
            (axes[0, 1], icc_als, "Mean ICC of ALS Metrics", (0.55, 1.05), np.arange(0.6, 1.01, 0.1), "(b)")]:
        for shape, color in zip(plot_shapes, series_colors):
            ax.plot(plot_sizes, values[shape], "-", marker=shape_marker[shape], color=color,
                    lw=1.2, ms=5, label=shape.title())
        ax.set_xlabel("Sub-plot size (m$^2$)")
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        if yticks is not None:
            ax.set_yticks(yticks)
        ax.set_xticks(plot_sizes)
        ax.grid(alpha=0.3, lw=0.4)
        ax.legend(loc="lower right", frameon=False)
        ax.text(0.03, 0.95, label, transform=ax.transAxes, fontsize=13, va="top", ha="left")

    positions = [0, 1, 2, 3.6, 4.6, 5.6]
    colors = ["#1f77b4", "#7fbf7b", "#d95f02"] * 2
    for ax, values, ylabel, label in [(axes[1, 0], cv_r2, "Calibration R$^2$ (300 m$^2$)", "(c)"),
                                      (axes[1, 1], cv_rmse, "Calibration RMSE (Mg ha$^{-1}$, 300 m$^2$)", "(d)")]:
        data = [values[k][shape] for shape in plot_shapes for k in ["random", "grouped", "independent"]]
        bp = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True,
                        showfliers=False, medianprops=dict(color="black", lw=0.7),
                        boxprops=dict(lw=0.5), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        ax.set_xticks([1, 4.6])
        ax.set_xticklabels(["Circle", "Square"])
        ax.set_ylabel(ylabel)
        ax.text(0.03, 0.95, label, transform=ax.transAxes, fontsize=13, va="top", ha="left")
    axes[1, 0].set_ylim(-0.1, 1.45)          # Headroom for the legend.
    axes[1, 0].set_yticks(np.arange(0, 1.01, 0.2))
    axes[1, 0].axhline(0, color="gray", ls="--", lw=0.4)
    axes[1, 0].legend(handles=[Patch(facecolor="#1f77b4", alpha=0.7, label="Random 5-fold CV"),
                               Patch(facecolor="#7fbf7b", alpha=0.7, label="Grouped 5-fold CV"),
                               Patch(facecolor="#d95f02", alpha=0.7, label="One sub-plot per parent")],
                      loc="upper center", bbox_to_anchor=(0.5, 1.0), frameon=False, fontsize=10)
    axes[1, 1].set_ylim(0, 70)
    for ax in axes.ravel():
        style_axes(ax)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.27, wspace=0.27)
    fig.savefig(os.path.join(fig_dir, "Fig11_O3_autocorrelation_test.png"), dpi=500, bbox_inches="tight")
    plt.close(fig)
    print("  Fig11_O3_autocorrelation_test.png")

print("\n" + "-" * 70)
print(f"Analysis complete. Results: {res_dir}; figures: {fig_dir}")
print("-" * 70)

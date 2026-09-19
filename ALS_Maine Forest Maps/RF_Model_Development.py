# Created on Sep 7, 2026, Author: Xinyuan Wei

import os
import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GroupKFold
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

# ----------------------------------------------------------
# Parameters.
# ----------------------------------------------------------

Data_dir = r"C:\Users\xinyuan.wei\Desktop\LiDAR-5.Maine Biomass\1.Data\1.DigiSylva Data"
Out_dir = r"C:\Users\xinyuan.wei\Desktop\LiDAR-5.Maine Biomass\1.Data\2.Model and Programs\RF_Models"
Stats_dir = r"C:\Users\xinyuan.wei\Desktop\LiDAR-5.Maine Biomass\1.Data\2.Model and Programs\Stats"
Out_file = "Model Comparison.xlsx"

# 3DEP campaigns near the field years, so the response is measured rather than back-cast.
Campaigns = ["2021_3DEP", "2022_3DEP", "2023_3DEP", "2024_3DEP"]

# Predictors and the plot table column each one comes from, with the campaign filled in.
# The plot table stops at p90, so p90 stands in for P95.
Predictor_columns = {"p25": "p25_%s", "p50": "p50_%s", "p75": "p75_%s", "p90": "p90_%s",
                     "max_height": "max_height_%s", "mean_height": "mean_height_%s",
                     "density": "point_density_m2_%s"}
Predictors = list(Predictor_columns)

# The single height metric the allometric form is built on.
Height_metric = "p50"

# Models to build, with the response column, the factor that reaches the reported unit, and that unit.
# The field carbon is stored as kg m-2, so a factor of 10 gives Mg C ha-1. Volume is already m3 ha-1.
Models = [("ME_RF_AGB", "Carbon", "AGB_%s_mean", 10.0, "Mg C ha-1"),
          ("ME_RF_Vol", "Volume", "Volume_%s_W", 1.0, "m3 ha-1")]

# Candidates compared on the same folds. The kind says whether the response is fitted on the
# log scale, and the columns say how much of the metric set each one is allowed to see.
Candidates = [("Power law on height", "log", [Height_metric]),
              ("Log log on all metrics", "log", Predictors),
              ("Multiple linear regression", "plain", Predictors),
              ("k nearest neighbours", "plain", Predictors),
              ("Support vector regression", "plain", Predictors),
              ("Gradient boosting", "plain", Predictors),
              ("Random forest", "plain", Predictors)]

# The candidate carried through to the maps.
Chosen = "Random forest"

# Conifer species in the field data. Everything else is counted as hardwood.
Conifers = ["ABBA", "PIAB", "PIBA", "PIRE", "PIRU", "PIST", "THOC", "TSCA"]

# Share of plot carbon held by conifers, used to sort each plot into a forest type.
Softwood_share = 0.75
Hardwood_share = 0.25

# Sites grouped so a regional model can be fitted. WM is western Maine and SD is Schoodic.
Western_sites = ["WM"]
Central_sites = ["DD", "HL", "PN"]

# Folds for the plot level cross validation.
Folds = 5

# A log scale needs a positive number, so metrics and responses are floored here.
Min_positive = 0.01

# A group needs this many rows before its score is worth reporting.
Min_group = 30

Seed = 42

# ----------------------------------------------------------
# Load the plot table.
# ----------------------------------------------------------

# The plot table is wide, so the site column is joined on rather than inserted into it.
plot = pd.read_csv(os.path.join(Data_dir, "DigiSylva_Plot_Data.csv"))
plot = pd.concat([plot, plot["Plot"].str.split("_").str[0].rename("Site")], axis=1)

# Plots with no trees are kept. They are real openings and they anchor the low end.
print("plots with no trees kept in: %s" % ", ".join(plot.loc[plot["N_Trees"] == 0, "Plot"]))

# ----------------------------------------------------------
# Forest type of each plot, from the field trees.
# ----------------------------------------------------------

tree = pd.read_csv(os.path.join(Data_dir, "DigiSylva_Tree_Data.csv"))
live = tree[tree["Snag"].isna() & tree["Biomass_mean"].notna()]
# The tree table names a tree, so the plot it stands on is the name with the tree number taken off.
stand = live["Plot"].str.rsplit("_", n=1).str[0]
conifer = live["Species"].isin(Conifers) * live["Biomass_mean"]
share = conifer.groupby(stand).sum() / live["Biomass_mean"].groupby(stand).sum()
kind = pd.cut(share, [-0.01, Hardwood_share, Softwood_share, 1.01],
              labels=["Hardwood", "Mixedwood", "Softwood"])
plot = plot.merge(kind.rename("Forest type"), left_on="Plot", right_index=True, how="left")
plot["Forest type"] = plot["Forest type"].astype(object).fillna("Hardwood")
print("plots by forest type: %s" % plot["Forest type"].value_counts().to_dict())
print("plots by site: %s" % plot["Site"].value_counts().to_dict())

# ----------------------------------------------------------
# Training table, one row per plot and campaign.
# ----------------------------------------------------------

frames = []
for camp in Campaigns:
    year = camp.split("_")[0]
    keep = plot["n_points_" + camp].notna()
    # Built as one dictionary and added in a single step, so the frame is not rewritten column by column.
    columns = {"Campaign": camp}
    for _, target, pattern, scale, _ in Models:
        columns[target] = plot.loc[keep, pattern % year].values * scale
    for predictor, pattern in Predictor_columns.items():
        columns[predictor] = plot.loc[keep, pattern % camp].values
    frames.append(plot.loc[keep, ["Plot", "Site", "Forest type"]].assign(**columns))

data = pd.concat(frames, ignore_index=True)
print("rows %d, plots %d, sites %s" % (len(data), data["Plot"].nunique(), sorted(data["Site"].unique())))

# ----------------------------------------------------------
# Candidate models on the same folds.
# ----------------------------------------------------------

# The estimator behind a candidate, built fresh for every fold so nothing leaks between them.
def make_model(label):
    if label == "k nearest neighbours":
        return make_pipeline(StandardScaler(), KNeighborsRegressor(n_neighbors=5, weights="distance"))
    if label == "Support vector regression":
        return make_pipeline(StandardScaler(), SVR(C=100.0, epsilon=2.0))
    if label == "Gradient boosting":
        return GradientBoostingRegressor(n_estimators=500, learning_rate=0.05, max_depth=3,
                                         random_state=Seed)
    if label == "Random forest":
        return RandomForestRegressor(n_estimators=500, min_samples_leaf=3, max_features=0.5,
                                     random_state=Seed, n_jobs=-1)
    return LinearRegression()

# Cross validated predictions folded by plot, so no plot sits in both halves of a fold.
# A log fit is trained on the stocked plots only and carried back with the smearing factor,
# because the mean of a log fit is not the log of the mean.
def run_cv(label, kind, columns, target):
    x = np.maximum(data[columns].values, Min_positive)
    y = data[target].values
    g = data["Plot"].values
    splitter = GroupKFold(n_splits=Folds)
    pred = np.zeros(len(y))
    for train_i, test_i in splitter.split(x, y, g):
        model = make_model(label)
        if kind == "log":
            stocked = train_i[y[train_i] > Min_positive]
            model.fit(np.log(x[stocked]), np.log(y[stocked]))
            smear = np.exp(np.log(y[stocked]) - model.predict(np.log(x[stocked]))).mean()
            pred[test_i] = np.exp(model.predict(np.log(x[test_i]))) * smear
        else:
            model.fit(x[train_i], y[train_i])
            pred[test_i] = model.predict(x[test_i])
    return np.maximum(pred, 0.0)

# A set of scores for one vector of predictions against one vector of field values.
def score(y, pred):
    rmse = np.sqrt(((y - pred) ** 2).mean())
    return {"n": len(y), "field mean": y.mean(), "R2": 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum(),
            "RMSE": rmse, "rRMSE %": 100 * rmse / y.mean(),
            "MAE": np.abs(y - pred).mean(), "bias": (pred - y).mean()}

rows = []
held = {}
for name, target, _, _, unit in Models:
    for label, kind, columns in Candidates:
        pred = run_cv(label, kind, columns, target)
        run = {"response": target, "unit": unit, "model": label, "predictors": len(columns)}
        run.update(score(data[target].values, pred))
        rows.append(run)
        held[(target, label)] = pred
        print("%-10s %-28s R2 %.3f  rRMSE %5.1f%%  bias %6.2f" %
              (target, label, run["R2"], run["rRMSE %"], run["bias"]), flush=True)

compare = pd.DataFrame(rows)
show = ["response", "model", "predictors", "n", "R2", "RMSE", "rRMSE %", "MAE", "bias"]
print("\ncandidate models, cross validated by plot, in the reported unit of each response:")
print(compare[show].round(3).to_string(index=False))

# ----------------------------------------------------------
# How the candidates hold up by forest type and by site.
# ----------------------------------------------------------

groups = []
for name, target, _, _, unit in Models:
    for label, _, _ in Candidates:
        pred = held[(target, label)]
        for field in ("Forest type", "Site"):
            for level, pick in data.groupby(field).groups.items():
                take = data.index.isin(pick)
                if take.sum() < Min_group:
                    continue
                run = {"response": target, "model": label, "split": field, "group": str(level)}
                run.update(score(data[target].values[take], pred[take]))
                groups.append(run)
by_group = pd.DataFrame(groups)
for field in ("Forest type", "Site"):
    print("\nrelative error by %s, carbon, %% of the field mean:" % field.lower())
    part = by_group[(by_group.split == field) & (by_group.response == "Carbon")]
    print(part.pivot_table(index="model", columns="group", values="rRMSE %").round(1).to_string())

# The chosen model has to work everywhere, not only on average, so it is judged on its worst group.
carbon = by_group[by_group.response == "Carbon"].groupby("model")
steady = pd.DataFrame({"worst group rRMSE %": carbon["rRMSE %"].max(),
                       "spread of bias": carbon["bias"].max() - carbon["bias"].min()})
steady = steady.sort_values("worst group rRMSE %")
print("\nhow far each candidate slips on its weakest group, carbon:")
print(steady.round(2).to_string())
print("\ncarried through to the maps: %s" % Chosen)

# ----------------------------------------------------------
# The chosen model fitted on parts of the plot set.
# ----------------------------------------------------------

# Each part is a rule on the training table. The last one is the model used in the study.
parts = [("Softwood", data["Forest type"] == "Softwood"),
         ("Hardwood", data["Forest type"] == "Hardwood"),
         ("Mixedwood", data["Forest type"] == "Mixedwood"),
         ("Western Maine", data["Site"].isin(Western_sites)),
         ("Central Maine", data["Site"].isin(Central_sites)),
         ("All plots", data["Site"].notna())]

# Cross validated predictions inside one part, folded by plot so no plot sits in both halves.
def run_part(take, target):
    x = np.maximum(data.loc[take, Predictors].values, Min_positive)
    y = data.loc[take, target].values
    g = data.loc[take, "Plot"].values
    splitter = GroupKFold(n_splits=min(Folds, len(np.unique(g))))
    pred = np.zeros(len(y))
    for train_i, test_i in splitter.split(x, y, g):
        model = make_model(Chosen)
        model.fit(x[train_i], y[train_i])
        pred[test_i] = model.predict(x[test_i])
    return y, np.maximum(pred, 0.0)

fitted = []
for label, take in parts:
    run = {"part": label, "plots": int(data.loc[take, "Plot"].nunique())}
    for _, target, _, _, _ in Models:
        y, pred = run_part(take, target)
        for key, value in score(y, pred).items():
            run["%s %s" % (target, key)] = value
        # The model fitted on every plot, scored on the same rows, so the two are comparable.
        run["%s pooled rRMSE %%" % target] = score(y, held[(target, Chosen)][take.values])["rRMSE %"]
    fitted.append(run)
    print("%-14s %4d plots  AGB R2 %.3f  volume R2 %.3f" %
          (label, run["plots"], run["Carbon R2"], run["Volume R2"]), flush=True)

stratified = pd.DataFrame(fitted)
show_parts = ["part", "plots", "Carbon R2", "Carbon rRMSE %", "Carbon bias",
              "Volume R2", "Volume rRMSE %", "Volume bias",
              "Carbon pooled rRMSE %", "Volume pooled rRMSE %"]
print("\nthe chosen model fitted on parts of the plot set, cross validated by plot:")
print(stratified[show_parts].round(3).to_string(index=False))

# ----------------------------------------------------------
# Fit and save the chosen models.
# ----------------------------------------------------------

for name, target, _, _, unit in Models:
    forest = RandomForestRegressor(n_estimators=1000, min_samples_leaf=3, max_features=0.5,
                                   random_state=Seed, n_jobs=-1)
    forest.fit(data[Predictors], data[target])
    bundle = {"model": forest,
              "predictors": Predictors,
              "response": target,
              "unit": unit,
              "n_rows": len(data),
              "n_plots": int(data["Plot"].nunique()),
              "campaigns": Campaigns,
              "plot_area_m2": 250.0,
              "sklearn": sklearn.__version__}
    joblib.dump(bundle, os.path.join(Out_dir, name + ".joblib"))
    print("\n%s importance" % name)
    print(pd.Series(forest.feature_importances_, index=Predictors).sort_values(ascending=False).round(3).to_string())

# ----------------------------------------------------------
# Save the tables.
# ----------------------------------------------------------

os.makedirs(Stats_dir, exist_ok=True)
out_path = os.path.join(Stats_dir, Out_file)
with pd.ExcelWriter(out_path) as writer:
    compare[show].to_excel(writer, sheet_name="Candidates", index=False)
    by_group.to_excel(writer, sheet_name="By Group", index=False)
    steady.to_excel(writer, sheet_name="Weakest Group")
    stratified[show_parts].to_excel(writer, sheet_name="Parts", index=False)
    data.to_excel(writer, sheet_name="Training Table", index=False)
print("\nsaved %s" % out_path)

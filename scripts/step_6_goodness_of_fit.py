#!/usr/bin/env python
"""
step_6_goodness_of_fit.py

Generate goodness-of-fit plots for the three modeling steps:
  1) Hill dose–response fits
  2) Derepression ODE fits (REV)
  3) Repression ODE fits (KD)

Requires the output parameter files from previous steps.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

import step_1c_fit_hill_curves as step1c
import step_2_simulate_derepression as step2
import step_3_simulate_repression as step3

# ---------------------------------------------------------------------
# Paths and small helpers
# ---------------------------------------------------------------------
PARAM_PATH = Path("parameters")
PLOTS_ROOT = Path("plots")
OUT_PATH = PLOTS_ROOT / "goodness_of_fit"


def ensure_outdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def save_current(fig, out_path: Path, filename: str) -> None:
    out_file = out_path / filename
    fig.tight_layout()
    fig.savefig(out_file, bbox_inches="tight", dpi=300)
    plt.close(fig)


def r_squared(y_true, y_pred):
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    y_true = y_true[mask]
    y_pred = y_pred[mask]
    if len(y_true) < 2:
        return np.nan
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot == 0:
        return np.nan
    return 1.0 - ss_res / ss_tot


def mae(y_true, y_pred):
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    y_true = y_true[mask]
    y_pred = y_pred[mask]
    if len(y_true) == 0:
        return np.nan
    return float(np.mean(np.abs(y_true - y_pred)))


# ---------------------------------------------------------------------
# 1. Hill fits: observed vs predicted fc
# ---------------------------------------------------------------------
def build_hill_dataset():
    """
    Rebuild the day-4 dose–response dataset.

    Returns a DataFrame with at least columns:
      plasmid, guide, dTAG, norm.bfp, fc

    """
    # NFC background
    mBFP_neg, mmCherry_neg = step1c.compute_nfc_background("fcs_files/NFC")

    # Replicates
    rep1 = step1c.load_replicate(
        "fcs_files/dose_response_data/R1", mBFP_neg, mmCherry_neg
    ).assign(rep=1, day=4)
    rep2 = step1c.load_replicate(
        "fcs_files/dose_response_data/R2", mBFP_neg, mmCherry_neg
    ).assign(rep=2, day=4)
    rep3 = step1c.load_replicate(
        "fcs_files/dose_response_data/R3", mBFP_neg, mmCherry_neg
    ).assign(rep=3, day=4)

    d4 = pd.concat([rep1, rep2, rep3], ignore_index=True)

    CH_BFP = step1c.CH_BFP  # "BV421-A"
    CH_mCh = step1c.CH_mCh  # "PE-A"

    # 1) mean over NTC rows, grouped only by plasmid
    meanNTC_pl = (
        d4.query("guide == 'N'")
        .groupby(["plasmid"])[CH_mCh]
        .mean()
        .rename("meanNTC")
        .reset_index()
    )

    # 2) expand to exact (plasmid, dTAG) pairs in NTC and join
    ntc_pairs = (
        d4.query("guide == 'N'")[["plasmid", "dTAG"]]
        .drop_duplicates()
        .merge(meanNTC_pl, on="plasmid", how="left")
    )
    d4 = d4.merge(ntc_pairs, on=["plasmid", "dTAG"], how="left")

    # 3) fold-change vs NTC
    d4["fc"] = d4[CH_mCh] / d4["meanNTC"]

    # We now apply the same min-max scaling used during fitting.
    d4["norm.bfp"] = (
        d4.groupby(["plasmid", "guide"], group_keys=False)[CH_BFP]
        .apply(lambda v: (v - v.min()) / (v.max() - v.min())
        if v.max() > v.min() else 0.0)
    )

    # Keep only rows actually used for Hill fits: guide == 'G'
    d4g = d4.query("guide == 'G'").copy()

    return d4g


def gof_hill(out_path: Path):
    """
    Observed vs predicted fold-change for Hill fits, across plasmids.

    Args
    -----
    out_path : Path
        Directory to save output plots.

    Returns
    -------
    None
    """
    hill_par = pd.read_csv(PARAM_PATH / "Hill_parameters.csv")
    d4g = build_hill_dataset()

    name_map = {
        "430": "SP430",
        "411": "SP411",
        "427": "SP427",
        "428": "SP428",
        "430ABA": "SP430A"
    }
    # Apply map; drop rows that don't match known targets
    d4g["plasmid"] = d4g["plasmid"].map(name_map)
    d4g = d4g.dropna(subset=["plasmid"])

    # Now the merge will find matches
    merged = d4g.merge(hill_par, on="plasmid", how="inner")

    # Predicted fc from K,n
    merged["fc_pred"] = step1c.hill_func(
        merged["norm.bfp"].to_numpy(float),
        merged["K"].to_numpy(float),
        merged["n"].to_numpy(float),
    )

    obs = merged["fc"].to_numpy(float)
    pred = merged["fc_pred"].to_numpy(float)
    mask = np.isfinite(obs) & np.isfinite(pred)
    obs, pred = obs[mask], pred[mask]
    combined = np.concatenate([obs, pred])
    combined = combined[np.isfinite(combined)]

    if combined.size == 0:
        print("[WARN] No valid data for Hill GOF plot (merged is empty).")
        return

    # Global R²
    r2 = r_squared(obs, pred)

    # Scatter plot (all plasmids together)
    fig, ax = plt.subplots()
    ax.scatter(obs, pred, alpha=0.5, s=10)
    max_val = float(combined.max()) * 1.05
    ax.plot([0, max_val], [0, max_val], linestyle="--", linewidth=1)
    ax.set_xlabel("Observed fc (mCherry / NTC)")
    ax.set_ylabel("Predicted fc (Hill model)")
    ax.set_title(f"Hill fits: observed vs predicted fc (R² = {r2:.2f})")
    save_current(fig, out_path, "gof_hill_fc_obs_vs_pred.pdf")

    # Optional: per-plasmid R² barplot
    rows = []
    for pl, sub in merged.groupby("plasmid"):
        o = sub["fc"].to_numpy(float)
        p = sub["fc_pred"].to_numpy(float)
        m = np.isfinite(o) & np.isfinite(p)
        if m.sum() > 2:
            rows.append({"plasmid": pl, "R2": r_squared(o[m], p[m])})

    if rows:
        r2_df = pd.DataFrame(rows).sort_values("R2")
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.barh(r2_df["plasmid"], r2_df["R2"])
        ax.set_xlabel("R²")
        ax.set_title("Hill fits: R² per plasmid")
        save_current(fig, out_path, "gof_hill_r2_per_plasmid.pdf")


# ---------------------------------------------------------------------
# 2. Derepression ODE (REV) GOF
# ---------------------------------------------------------------------
def build_rev_dataset():
    """
    Build the REV dataset for derepression fits.
    """
    mBFP_neg, mmCherry_neg = step2.compute_nfc_background("fcs_files/NFC")
    df_tc = step2.load_timecourse_with_bg(mBFP_neg, mmCherry_neg)
    rev = step2.compute_rev_transforms(df_tc)
    return rev


def gof_derepression(out_path: Path):
    """
    Observed vs predicted mCherry trajectories for derepression ODE fits.
    """
    rev = build_rev_dataset()
    pars = step2.load_parameters()  # same logic as in step_2_simulate_derepression
    delays_de = pd.read_csv(PARAM_PATH / "delays_derepression.csv")

    name_map = {
        "430": "SP430",
        "411": "SP411",
        "427": "SP427",
        "428": "SP428",
        "430ABA": "SP430A"
    }
    # Apply map; drop rows that don't match known targets
    rev["plasmid"] = rev["plasmid"].map(name_map)
    rev = rev.dropna(subset=["plasmid"])

    all_obs = []
    all_pred = []
    rows_stats = []

    for _, row in delays_de.iterrows():
        pl = row["plasmid"]
        best_delay = float(row["d_rev"])

        sub = rev.query("plasmid == @pl").copy()
        if sub.empty:
            continue

        # Mean over replicates at each time
        mean_fc = (
            sub.groupby("time", as_index=False)[["fc.cherry", "norm.bfp"]].mean()
        )
        # Initial conditions from time 0
        t0_rows = mean_fc.query("time == 0")
        if t0_rows.empty:
            continue
        R0 = float(t0_rows["norm.bfp"].iloc[0])
        Y0 = float(t0_rows["fc.cherry"].iloc[0])

        # Parameter row (K, n, t_down, alpha, etc.)
        par_row = pars.query("plasmid == @pl")
        if par_row.empty:
            continue
        par_row = par_row.iloc[0]

        # Simulate base ODE (no delay applied yet)
        sim = step2.simulate_ode(
            R0=R0,
            Y0=Y0,
            pars=par_row,
            tmax=150.0,
            step=0.05,
            delay=0.0,
        )
        # Interpolator for Y(t)
        pchipY = PchipInterpolator(sim["time"], sim["Y"], extrapolate=False)

        # Align with best delay: experiment time t >= d
        m = mean_fc["time"] >= best_delay
        t_data = mean_fc.loc[m, "time"].to_numpy(float)
        y_obs = mean_fc.loc[m, "fc.cherry"].to_numpy(float)
        t_query = t_data - best_delay
        y_pred = pchipY(t_query)

        mask = np.isfinite(y_obs) & np.isfinite(y_pred)
        y_obs = y_obs[mask]
        y_pred = y_pred[mask]

        if len(y_obs) == 0:
            continue

        all_obs.append(y_obs)
        all_pred.append(y_pred)

        r2_pl = r_squared(y_obs, y_pred)
        mae_pl = mae(y_obs, y_pred)
        rows_stats.append({"plasmid": pl, "R2": r2_pl, "MAE": mae_pl})

        # Per-plasmid scatter
        fig, ax = plt.subplots()
        ax.scatter(y_obs, y_pred, s=20, alpha=0.8)
        vmax = float(np.nanmax([y_obs.max(), y_pred.max()])) * 1.05
        ax.plot([0, vmax], [0, vmax], linestyle="--", linewidth=1)
        ax.set_xlabel("Observed fc.cherry")
        ax.set_ylabel("Predicted fc.cherry (ODE REV)")
        ax.set_title(f"Derepression GOF – {pl} (R²={r2_pl:.2f}, MAE={mae_pl:.3g})")
        save_current(fig, out_path, f"gof_rev_obs_vs_pred_{pl}.pdf")

        # Residual vs time
        residuals = y_obs - y_pred
        fig, ax = plt.subplots()
        ax.axhline(0.0, linestyle="--", linewidth=1)
        ax.scatter(t_data[mask], residuals, s=20)
        ax.set_xlabel("Time (h)")
        ax.set_ylabel("Residual (obs - pred)")
        ax.set_title(f"Derepression residuals – {pl}")
        save_current(fig, out_path, f"gof_rev_residuals_{pl}.pdf")

    # Combined scatter
    if all_obs:
        obs_all = np.concatenate(all_obs)
        pred_all = np.concatenate(all_pred)
        r2_all = r_squared(obs_all, pred_all)

        fig, ax = plt.subplots()
        ax.scatter(obs_all, pred_all, s=12, alpha=0.6)
        vmax = float(np.nanmax([obs_all.max(), pred_all.max()])) * 1.05
        ax.plot([0, vmax], [0, vmax], linestyle="--", linewidth=1)
        ax.set_xlabel("Observed fc.cherry")
        ax.set_ylabel("Predicted fc.cherry (ODE REV)")
        ax.set_title(f"Derepression GOF – all plasmids (R²={r2_all:.2f})")
        save_current(fig, out_path, "gof_rev_obs_vs_pred_all.pdf")

    # MAE / R² barplots per plasmid
    if rows_stats:
        stats_df = pd.DataFrame(rows_stats).sort_values("MAE")

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.bar(stats_df["plasmid"], stats_df["MAE"])
        ax.set_ylabel("MAE")
        ax.set_title("Derepression ODE – MAE per plasmid")
        save_current(fig, out_path, "gof_rev_mae_per_plasmid.pdf")

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.bar(stats_df["plasmid"], stats_df["R2"])
        ax.set_ylabel("R²")
        ax.set_title("Derepression ODE – R² per plasmid")
        save_current(fig, out_path, "gof_rev_r2_per_plasmid.pdf")


# ---------------------------------------------------------------------
# 3. Repression ODE (KD) GOF
# ---------------------------------------------------------------------
def build_kd_dataset():
    """
    Build the KD dataset for repression fits.
    """
    kd = step3.build_kd_dataset(nfc_dir="fcs_files/NFC",
                                tc_dir="fcs_files/time-course_data")
    return kd


def gof_repression(out_path: Path):
    """
    Observed vs predicted mCherry trajectories for repression ODE fits.
    """
    kd = build_kd_dataset()
    pars = step3.load_parameters()
    delays_rep = pd.read_csv(PARAM_PATH / "delays_repression.csv")

    name_map = {
        "430": "SP430",
        "411": "SP411",
        "427": "SP427",
        "428": "SP428",
        "430ABA": "SP430A"
    }
    # Apply map; drop rows that don't match known targets
    kd["plasmid"] = kd["plasmid"].map(name_map)
    kd = kd.dropna(subset=["plasmid"])

    all_obs = []
    all_pred = []
    rows_stats = []

    for _, row in delays_rep.iterrows():
        pl = row["plasmid"]
        best_delay = float(row["d_rev"])  # same column name as derepression

        sub = kd.query("plasmid == @pl").copy()
        if sub.empty:
            continue

        mean_fc = (
            sub.groupby("time", as_index=False)[["fc.cherry", "norm.bfp"]].mean()
        )

        t0_rows = mean_fc.query("time == 0")
        if t0_rows.empty:
            continue
        R0 = float(t0_rows["norm.bfp"].iloc[0])
        Y0 = float(t0_rows["fc.cherry"].iloc[0])

        par_row = pars.query("plasmid == @pl")
        if par_row.empty:
            continue
        par_row = par_row.iloc[0]

        sim = step3.simulate_ode(
            R0=R0,
            Y0=Y0,
            pars=par_row,
            tmax=150.0,
            step=0.05,
            delay=0.0,
        )

        pchipY = PchipInterpolator(sim["time"], sim["Y"], extrapolate=False)

        m = mean_fc["time"] >= best_delay
        t_data = mean_fc.loc[m, "time"].to_numpy(float)
        y_obs = mean_fc.loc[m, "fc.cherry"].to_numpy(float)
        t_query = t_data - best_delay
        y_pred = pchipY(t_query)

        mask = np.isfinite(y_obs) & np.isfinite(y_pred)
        y_obs = y_obs[mask]
        y_pred = y_pred[mask]

        if len(y_obs) == 0:
            continue

        all_obs.append(y_obs)
        all_pred.append(y_pred)

        r2_pl = r_squared(y_obs, y_pred)
        mae_pl = mae(y_obs, y_pred)
        rows_stats.append({"plasmid": pl, "R2": r2_pl, "MAE": mae_pl})

        fig, ax = plt.subplots()
        ax.scatter(y_obs, y_pred, s=20, alpha=0.8)
        vmax = float(np.nanmax([y_obs.max(), y_pred.max()])) * 1.05
        ax.plot([0, vmax], [0, vmax], linestyle="--", linewidth=1)
        ax.set_xlabel("Observed fc.cherry")
        ax.set_ylabel("Predicted fc.cherry (ODE KD)")
        ax.set_title(f"Repression GOF – {pl} (R²={r2_pl:.2f}, MAE={mae_pl:.3g})")
        save_current(fig, out_path, f"gof_kd_obs_vs_pred_{pl}.pdf")

        residuals = y_obs - y_pred
        fig, ax = plt.subplots()
        ax.axhline(0.0, linestyle="--", linewidth=1)
        ax.scatter(t_data[mask], residuals, s=20)
        ax.set_xlabel("Time (h)")
        ax.set_ylabel("Residual (obs - pred)")
        ax.set_title(f"Repression residuals – {pl}")
        save_current(fig, out_path, f"gof_kd_residuals_{pl}.pdf")

    if all_obs:
        obs_all = np.concatenate(all_obs)
        pred_all = np.concatenate(all_pred)
        r2_all = r_squared(obs_all, pred_all)

        fig, ax = plt.subplots()
        ax.scatter(obs_all, pred_all, s=12, alpha=0.6)
        vmax = float(np.nanmax([obs_all.max(), pred_all.max()])) * 1.05
        ax.plot([0, vmax], [0, vmax], linestyle="--", linewidth=1)
        ax.set_xlabel("Observed fc.cherry")
        ax.set_ylabel("Predicted fc.cherry (ODE KD)")
        ax.set_title(f"Repression GOF – all plasmids (R²={r2_all:.2f})")
        save_current(fig, out_path, "gof_kd_obs_vs_pred_all.pdf")

    if rows_stats:
        stats_df = pd.DataFrame(rows_stats).sort_values("MAE")

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.bar(stats_df["plasmid"], stats_df["MAE"])
        ax.set_ylabel("MAE")
        ax.set_title("Repression ODE – MAE per plasmid")
        save_current(fig, out_path, "gof_kd_mae_per_plasmid.pdf")

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.bar(stats_df["plasmid"], stats_df["R2"])
        ax.set_ylabel("R²")
        ax.set_title("Repression ODE – R² per plasmid")
        save_current(fig, out_path, "gof_kd_r2_per_plasmid.pdf")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    ensure_outdir(OUT_PATH)

    # 1) Hill dose–response fits
    gof_hill(OUT_PATH)

    # 2) Derepression ODE fits (REV)
    gof_derepression(OUT_PATH)

    # 3) Repression ODE fits (KD)
    gof_repression(OUT_PATH)


if __name__ == "__main__":
    main()

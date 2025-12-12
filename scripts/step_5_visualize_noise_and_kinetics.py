#!/usr/bin/env python
"""
step_5_visualize_noise_and_kinetics.py
------------------------------------

Generates plots:
  1. Noise vs mean plots (BFP / mCherry) – empirical time-series
  2. Hill K vs n
  3. Up vs down half-times
  4. Repression vs derepression delays
  5. Hierarchical noise summary (per plasmid / channel / experiment)
  6. Noise vs kinetic parameters (design-guided view)
  7. Time series summaries (mean ± SE) per plasmid
  8. Goodness-of-fit between empirical and hierarchical noise

Requires input CSV files in "parameters/" directory:
    - Hill_parameters.csv
    - half_times_upregulation.csv
    - half_times_downregulation.csv
    - delays_repression.csv
    - delays_derepression.csv
    - single_cell_noise_timeseries.csv
    - single_cell_noise_hierarchical.csv
"""

import os
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------
PARAM_PATH = Path("parameters")
PLOTS_ROOT = Path("plots")
OUT_PATH = PLOTS_ROOT / "noise_kinetics"


def ensure_outdir(path: Path) -> None:
    """
    Ensure output directory exists.

    Args:
        path (Path): Path to the output directory.
    """

    path.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------
def save_current(fig, out_path: Path, filename: str) -> None:
    out_file = out_path / filename
    fig.tight_layout()
    fig.savefig(out_file, bbox_inches="tight", dpi=300)
    plt.close(fig)


def labeled_scatter(ax, x, y, labels, xlabel, ylabel, title):
    x = np.asarray(x)
    y = np.asarray(y)
    labels = np.asarray(labels)

    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    labels = labels[mask]

    ax.scatter(x, y)
    for xi, yi, lab in zip(x, y, labels):
        ax.text(xi, yi, str(lab), fontsize=8)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


# ---------------------------------------------------------------------
# Load all parameter / noise tables
# ---------------------------------------------------------------------
def load_tables() -> Dict[str, pd.DataFrame]:
    """
    Load all required CSV tables into a dictionary of DataFrames.

    1. Hill parameters
    2. Half-times upregulation
    3. Half-times downregulation
    4. Alpha values
    5. Delays derepression
    6. Delays repression
    7. Single-cell noise time series
    8. Single-cell noise hierarchical

    Returns:
        dict: Dictionary of DataFrames keyed by table name.
    """
    hill = pd.read_csv(PARAM_PATH / "Hill_parameters.csv")  # plasmid, K, n
    up = pd.read_csv(PARAM_PATH / "half_times_upregulation.csv")  # plasmid, halftime, se
    down = pd.read_csv(PARAM_PATH / "half_times_downregulation.csv")  # plasmid, halftime, se
    alpha = pd.read_csv(PARAM_PATH / "alphamcherry.csv")  # alpha
    delays_de = pd.read_csv(PARAM_PATH / "delays_derepression.csv")  # plasmid, d_rev
    delays_rep = pd.read_csv(PARAM_PATH / "delays_repression.csv")  # plasmid, d_rev
    design = pd.read_csv(PARAM_PATH / "design_space_scan_repression.csv")  # not used here
    noise_ts = pd.read_csv(PARAM_PATH / "single_cell_noise_timeseries.csv")
    noise_hier = pd.read_csv(PARAM_PATH / "single_cell_noise_hierarchical.csv")
    return {
        "hill": hill,
        "up": up,
        "down": down,
        "alpha": alpha,
        "delays_de": delays_de,
        "delays_rep": delays_rep,
        "design": design,
        "noise_ts": noise_ts,
        "noise_hier_raw": noise_hier,
    }


# ---------------------------------------------------------------------
# Convert hierarchical noise (wide) → long with 'channel'
# ---------------------------------------------------------------------
def make_noise_hier_long(noise_hier_raw: pd.DataFrame) -> pd.DataFrame:
    """
    Build a long-format table with explicit 'channel' from the wide
    single_cell_noise_hierarchical.csv.

    Args:
        noise_hier_raw (pd.DataFrame): Raw hierarchical noise DataFrame.

    Returns:
        pd.DataFrame: Long-format hierarchical noise DataFrame.
    """
    records = []
    for _, row in noise_hier_raw.iterrows():
        for channel, prefix in [("BFP", "BFP"), ("mCherry", "mCherry")]:
            records.append(
                {
                    "plasmid": row["plasmid"],
                    "exp": row["exp"],
                    "n_groups": row["n_groups"],
                    "channel": channel,
                    "cv2": row[f"cv2_{prefix}"],
                    "cv2_se": row[f"cv2_{prefix}_se"],
                    "cv2_ci_low": row[f"cv2_{prefix}_ci_low"],
                    "cv2_ci_high": row[f"cv2_{prefix}_ci_high"],
                    "mean": row[f"mean_{prefix}"],
                    "mean_se": row[f"mean_{prefix}_se"],
                    "mean_ci_low": row[f"mean_{prefix}_ci_low"],
                    "mean_ci_high": row[f"mean_{prefix}_ci_high"],
                }
            )
    df_long = pd.DataFrame(records)
    return df_long


# ---------------------------------------------------------------------
# 1. Noise vs mean plots (BFP / mCherry) – empirical time-series
# ---------------------------------------------------------------------
def plot_noise_vs_mean(noise_ts: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots()
    ax.scatter(noise_ts["mean_BFP"], noise_ts["cv2_BFP"], alpha=0.5)
    ax.set_xlabel("Mean BFP")
    ax.set_ylabel("CV² (BFP)")
    ax.set_title("Single-cell noise vs mean (BFP)")
    save_current(fig, out_path, "noise_vs_mean_BFP.pdf")

    fig, ax = plt.subplots()
    ax.scatter(noise_ts["mean_mCherry"], noise_ts["cv2_mCherry"], alpha=0.5)
    ax.set_xlabel("Mean mCherry")
    ax.set_ylabel("CV² (mCherry)")
    ax.set_title("Single-cell noise vs mean (mCherry)")
    save_current(fig, out_path, "noise_vs_mean_mCherry.pdf")


# ---------------------------------------------------------------------
# 2. Hill K vs n
# ---------------------------------------------------------------------
def plot_hill_params(hill: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots()
    labeled_scatter(
        ax,
        hill["K"],
        hill["n"],
        hill["plasmid"],
        xlabel="K (midpoint)",
        ylabel="n (Hill coefficient)",
        title="Hill parameters per plasmid",
    )
    save_current(fig, out_path, "hill_K_vs_n.pdf")


# ---------------------------------------------------------------------
# 3. Up vs down half-times
# ---------------------------------------------------------------------
def plot_half_times(up: pd.DataFrame, down: pd.DataFrame, out_path: Path) -> None:
    up_ = up.rename(columns={"halftime": "halftime_up"})
    down_ = down.rename(columns={"halftime": "halftime_down"})
    merged = pd.merge(
        up_[["plasmid", "halftime_up"]],
        down_[["plasmid", "halftime_down"]],
        on="plasmid",
    )

    fig, ax = plt.subplots()
    labeled_scatter(
        ax,
        merged["halftime_up"],
        merged["halftime_down"],
        merged["plasmid"],
        xlabel="Upregulation half-time (h)",
        ylabel="Downregulation half-time (h)",
        title="Kinetic comparison: up vs down",
    )
    save_current(fig, out_path, "half_times_up_vs_down.pdf")


# ---------------------------------------------------------------------
# 4. Repression vs derepression delays
# ---------------------------------------------------------------------
def plot_delays(delays_de: pd.DataFrame,
                delays_rep: pd.DataFrame,
                out_path: Path) -> None:
    de = delays_de.rename(columns={"d_rev": "d_derep"})
    rep = delays_rep.rename(columns={"d_rev": "d_rep"})
    merged = pd.merge(de, rep, on="plasmid")

    fig, ax = plt.subplots()
    labeled_scatter(
        ax,
        merged["d_derep"],
        merged["d_rep"],
        merged["plasmid"],
        xlabel="Derepression delay (h)",
        ylabel="Repression delay (h)",
        title="Repression vs derepression delays",
    )
    save_current(fig, out_path, "delays_repression_vs_derepression.pdf")


# ---------------------------------------------------------------------
# 5. Hierarchical noise summary (per plasmid / channel / experiment)
# ---------------------------------------------------------------------
def plot_noise_hier_summary(noise_hier_long: pd.DataFrame, out_path: Path) -> None:
    ordered_plasmids = sorted(noise_hier_long["plasmid"].unique())
    ordered_channels = ["BFP", "mCherry"]
    ordered_exps = sorted(noise_hier_long["exp"].unique())

    # ---- CV² per plasmid / channel / experiment ----
    plt.figure(figsize=(8, 5))
    x_positions = np.arange(len(ordered_plasmids))

    for i, channel in enumerate(ordered_channels):
        subset_ch = noise_hier_long[noise_hier_long["channel"] == channel]
        for j, exp in enumerate(ordered_exps):
            sub = subset_ch[subset_ch["exp"] == exp]
            sub = sub.set_index("plasmid").reindex(ordered_plasmids)
            means = sub["cv2"].values
            low = sub["cv2_ci_low"].values
            high = sub["cv2_ci_high"].values
            errs = [means - low, high - means]

            offset = (i + j * 0.2 - 0.3) * 0.25
            plt.errorbar(
                x_positions + offset,
                means,
                yerr=errs,
                fmt="o",
                label=f"{channel} {exp}",
            )

    plt.xticks(x_positions, ordered_plasmids)
    plt.ylabel("CV² (hierarchical estimate)")
    plt.title("Hierarchical noise per plasmid / channel / condition")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path / "noise_hier_cv2_summary.png", dpi=300)
    plt.close()

    # ---- Mean vs CV² scatter (hierarchical) ----
    plt.figure(figsize=(7, 5))
    for channel in ordered_channels:
        for exp in ordered_exps:
            s = noise_hier_long[
                (noise_hier_long["channel"] == channel) &
                (noise_hier_long["exp"] == exp)
                ]
            plt.scatter(
                s["mean"],
                s["cv2"],
                label=f"{channel} {exp}",
                alpha=0.7,
            )

    plt.xlabel("Mean expression (hierarchical)")
    plt.ylabel("CV² (hierarchical)")
    plt.title("Noise vs mean (hierarchical, BFP / mCherry)")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path / "noise_hier_mean_vs_cv2.png", dpi=300)
    plt.close()

    # ---- Per-plasmid CV²: BFP vs mCherry ----
    plt.figure(figsize=(6, 6))
    bfp = noise_hier_long[noise_hier_long["channel"] == "BFP"]
    mch = noise_hier_long[noise_hier_long["channel"] == "mCherry"]

    merged = (
        bfp[["plasmid", "exp", "cv2"]].rename(columns={"cv2": "cv2_BFP"})
        .merge(
            mch[["plasmid", "exp", "cv2"]].rename(columns={"cv2": "cv2_mCherry"}),
            on=["plasmid", "exp"],
        )
    )

    plt.scatter(merged["cv2_BFP"], merged["cv2_mCherry"])
    for _, r in merged.iterrows():
        plt.text(
            r["cv2_BFP"],
            r["cv2_mCherry"],
            f"{r['plasmid']}_{r['exp']}",
            fontsize=8,
        )

    max_val = float(
        np.nanmax([merged["cv2_BFP"].max(), merged["cv2_mCherry"].max()])
    ) * 1.1
    plt.plot([0, max_val], [0, max_val], linewidth=1, linestyle="--")
    plt.xlim(0, max_val)
    plt.ylim(0, max_val)
    plt.xlabel("CV² BFP (hierarchical)")
    plt.ylabel("CV² mCherry (hierarchical)")
    plt.title("Per-plasmid noise: BFP vs mCherry")
    plt.tight_layout()
    plt.savefig(out_path / "noise_hier_bfp_vs_mcherry.png", dpi=300)
    plt.close()


# ---------------------------------------------------------------------
# 6. Noise vs kinetic parameters
# ---------------------------------------------------------------------
def plot_noise_vs_kinetics(noise_hier_long: pd.DataFrame,
                           hill: pd.DataFrame,
                           up: pd.DataFrame,
                           down: pd.DataFrame,
                           delays_rep: pd.DataFrame,
                           out_path: Path) -> None:
    # Focus on mCherry channel for noise–kinetic links
    # (as it's the regulated reporter)
    hier_mc = noise_hier_long[noise_hier_long["channel"] == "mCherry"].copy()
    hier_mc = hier_mc.rename(columns={"cv2": "cv2_mCherry"})

    up_ = up.rename(columns={"halftime": "halftime_up"})
    down_ = down.rename(columns={"halftime": "halftime_down"})
    rep = delays_rep.rename(columns={"d_rev": "d_rep"})

    merged = (
        hier_mc
        .merge(hill, on="plasmid", how="left")
        .merge(up_[["plasmid", "halftime_up"]], on="plasmid", how="left")
        .merge(down_[["plasmid", "halftime_down"]], on="plasmid", how="left")
        .merge(rep[["plasmid", "d_rep"]], on="plasmid", how="left")
    )

    # (a) cv2 vs K
    fig, ax = plt.subplots()
    labeled_scatter(
        ax,
        merged["K"],
        merged["cv2_mCherry"],
        merged["plasmid"],
        xlabel="K (Hill midpoint)",
        ylabel="CV² mCherry (hierarchical)",
        title="Noise vs repression strength (K)",
    )
    save_current(fig, out_path, "noise_vs_K_mCherry.pdf")

    # (b) cv2 vs Hill coefficient n
    fig, ax = plt.subplots()
    labeled_scatter(
        ax,
        merged["n"],
        merged["cv2_mCherry"],
        merged["plasmid"],
        xlabel="n (Hill coefficient)",
        ylabel="CV² mCherry (hierarchical)",
        title="Noise vs cooperativity (n)",
    )
    save_current(fig, out_path, "noise_vs_n_mCherry.pdf")

    # (c) cv2 vs up half-time
    fig, ax = plt.subplots()
    labeled_scatter(
        ax,
        merged["halftime_up"],
        merged["cv2_mCherry"],
        merged["plasmid"],
        xlabel="Upregulation half-time (h)",
        ylabel="CV² mCherry (hierarchical)",
        title="Noise vs upregulation kinetics",
    )
    save_current(fig, out_path, "noise_vs_half_time_up_mCherry.pdf")

    # (d) cv2 vs repression delay
    fig, ax = plt.subplots()
    labeled_scatter(
        ax,
        merged["d_rep"],
        merged["cv2_mCherry"],
        merged["plasmid"],
        xlabel="Repression delay (h)",
        ylabel="CV² mCherry (hierarchical)",
        title="Noise vs repression delay",
    )
    save_current(fig, out_path, "noise_vs_repression_delay_mCherry.pdf")


# ---------------------------------------------------------------------
# 7. Time series summaries (mean ± SE) per plasmid
# ---------------------------------------------------------------------
def compute_ts_summary(noise_ts: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate single-cell time series to per-(plasmid, time) means and SEs.
    """

    def agg(group):
        mean_BFP = group["mean_BFP"].mean()
        mean_mC = group["mean_mCherry"].mean()
        se_BFP = np.sqrt((group["var_BFP"] / group["n_cells"]).mean())
        se_mC = np.sqrt((group["var_mCherry"] / group["n_cells"]).mean())
        return pd.Series(
            {
                "mean_BFP": mean_BFP,
                "se_BFP": se_BFP,
                "mean_mCherry": mean_mC,
                "se_mCherry": se_mC,
            }
        )

    summary = (
        noise_ts
        .groupby(["plasmid", "time"], as_index=False)
        .apply(agg, include_groups=False)
    )
    return summary


def plot_timeseries(summary: pd.DataFrame, out_path: Path) -> None:
    for plasmid, sub in summary.groupby("plasmid"):
        sub = sub.sort_values("time")

        fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
        ax1, ax2 = axes

        # BFP
        ax1.errorbar(
            sub["time"],
            sub["mean_BFP"],
            yerr=sub["se_BFP"],
            fmt="o-",
            capsize=3,
        )
        ax1.set_xlabel("Time (h)")
        ax1.set_ylabel("Mean BFP")
        ax1.set_title(f"{plasmid} – BFP dynamics")

        # mCherry
        ax2.errorbar(
            sub["time"],
            sub["mean_mCherry"],
            yerr=sub["se_mCherry"],
            fmt="o-",
            capsize=3,
        )
        ax2.set_xlabel("Time (h)")
        ax2.set_ylabel("Mean mCherry")
        ax2.set_title(f"{plasmid} – mCherry dynamics")

        fig.suptitle(f"Time-series mean ± SE – {plasmid}")
        save_current(fig, out_path, f"timeseries_{plasmid}.pdf")


# ---------------------------------------------------------------------
# 8. Goodness-of-fit between empirical and hierarchical noise
# ---------------------------------------------------------------------
def r_squared(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot == 0:
        return np.nan
    return 1.0 - ss_res / ss_tot


def plot_noise_gof(noise_ts: pd.DataFrame,
                   noise_hier_long: pd.DataFrame,
                   out_path: Path) -> None:
    # Empirical time-averaged mean and CV² for mCherry (across time & reps)
    emp = (
        noise_ts
        .groupby("plasmid", as_index=False)
        .agg(
            mean_mCherry_emp=("mean_mCherry", "mean"),
            cv2_mCherry_emp=("cv2_mCherry", "mean"),
        )
    )

    # Hierarchical estimates (mCherry channel, pooled over exp)
    hier_mc = noise_hier_long[noise_hier_long["channel"] == "mCherry"].copy()
    hier_mc = (
        hier_mc
        .groupby("plasmid", as_index=False)
        .agg(
            mean_mCherry_hier=("mean", "mean"),
            cv2_mCherry_hier=("cv2", "mean"),
        )
    )

    merged = emp.merge(hier_mc, on="plasmid", how="inner")

    # Mean comparison
    fig, ax = plt.subplots()
    ax.scatter(merged["mean_mCherry_emp"], merged["mean_mCherry_hier"])
    max_val = float(
        np.nanmax(
            [merged["mean_mCherry_emp"].max(), merged["mean_mCherry_hier"].max()]
        )
    )
    ax.plot([0, max_val], [0, max_val], linestyle="--")
    ax.set_xlabel("Empirical mean mCherry (time-avg.)")
    ax.set_ylabel("Hierarchical mean mCherry")
    rsq_mean = r_squared(
        merged["mean_mCherry_emp"].values,
        merged["mean_mCherry_hier"].values,
    )
    ax.set_title(f"Mean mCherry: empirical vs hierarchical (R² = {rsq_mean:.2f})")
    save_current(fig, out_path, "gof_mean_mCherry_emp_vs_hier.pdf")

    # CV² comparison
    fig, ax = plt.subplots()
    ax.scatter(merged["cv2_mCherry_emp"], merged["cv2_mCherry_hier"])
    max_val = float(
        np.nanmax(
            [merged["cv2_mCherry_emp"].max(), merged["cv2_mCherry_hier"].max()]
        )
    )
    ax.plot([0, max_val], [0, max_val], linestyle="--")
    ax.set_xlabel("Empirical CV² mCherry (time-avg.)")
    ax.set_ylabel("Hierarchical CV² mCherry")
    rsq_cv2 = r_squared(
        merged["cv2_mCherry_emp"].values,
        merged["cv2_mCherry_hier"].values,
    )
    ax.set_title(f"CV² mCherry: empirical vs hierarchical (R² = {rsq_cv2:.2f})")
    save_current(fig, out_path, "gof_cv2_mCherry_emp_vs_hier.pdf")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    ensure_outdir(OUT_PATH)
    tables = load_tables()

    hill = tables["hill"]
    up = tables["up"]
    down = tables["down"]
    delays_de = tables["delays_de"]
    delays_rep = tables["delays_rep"]
    noise_ts = tables["noise_ts"]
    noise_hier_raw = tables["noise_hier_raw"]

    # Build long-format hierarchical noise table once
    noise_hier_long = make_noise_hier_long(noise_hier_raw)

    # 1–4: global parameter relationships
    plot_noise_vs_mean(noise_ts, OUT_PATH)
    plot_hill_params(hill, OUT_PATH)
    plot_half_times(up, down, OUT_PATH)
    plot_delays(delays_de, delays_rep, OUT_PATH)

    # 5–6: hierarchical noise & noise–kinetic links
    plot_noise_hier_summary(noise_hier_long, OUT_PATH)
    plot_noise_vs_kinetics(noise_hier_long, hill, up, down, delays_rep, OUT_PATH)

    # 7: time-series summaries (mean ± SE)
    ts_summary = compute_ts_summary(noise_ts)
    plot_timeseries(ts_summary, OUT_PATH)

    # 8: goodness-of-fit empirical vs hierarchical noise
    plot_noise_gof(noise_ts, noise_hier_long, OUT_PATH)


if __name__ == "__main__":
    main()

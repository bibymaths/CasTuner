#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_7_design_selection_and_map.py
----------------------------------
This script analyzes the simulated design space of genetic constructs
to identify top-performing designs based on dynamic range and response time. It ranks designs using a composite score that favors high dynamic range
and low response time (t50). The top candidates are visualized in a performance map alongside real experimental constructs for context.
The final selected designs are saved for further experimental validation.
----------------------------------
"""

import os
import math
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pathlib import Path

warnings.filterwarnings("ignore")

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
PARAM_PATH = Path("parameters")
OUT_PATH = Path("plots")
os.makedirs(OUT_PATH, exist_ok=True)


# -----------------------------------------------------------------------------
# Simulation Logic
# -----------------------------------------------------------------------------
def ode_rhs(t, y, t_up, K, n, alpha):
    """
    ODE right-hand side for the genetic construct model.

    Args
        t (float): Current time.
        y (list): Current state [R, Y].
        t_up (float): Half-time for upregulation.
        K (float): Hill constant.
        n (float): Hill coefficient.
        alpha (float): Degradation rate of Y.

    Returns
        list: Derivatives [dR/dt, dY/dt].
    """
    R, Y = y
    t_up = max(float(t_up), 1e-6)
    beta = math.log(2.0) / t_up

    R_pos = max(R, 0.0)
    Kn = float(K) ** n
    Rn = R_pos ** n

    dR = beta - R_pos * (math.log(2.0) / t_up)
    dY = (Kn / (Kn + Rn)) - alpha * Y
    return [dR, dY]


def simulate_construct(t_up, K, n, alpha, t_end=72.0, dt=0.05):
    """
    Simulate the genetic construct dynamics over time.

    Args:
        t_up (float): Half-time for upregulation.
        K (float): Hill constant.
        n (float): Hill coefficient.
        alpha (float): Degradation rate of Y.
        t_end (float): End time for simulation.
        dt (float): Time step for evaluation.
    Returns:
        pd.DataFrame: Time series of Y over time.
    """
    t_eval = np.arange(0.0, t_end + dt / 2, dt)
    y0 = [0.0, 1.0 / max(float(alpha), 1e-12)]

    sol = solve_ivp(
        fun=lambda t, y: ode_rhs(t, y, t_up, K, n, alpha),
        t_span=(0, t_end),
        y0=y0,
        t_eval=t_eval,
        # method="LSODA"
    )
    if not sol.success:
        return None

    Y = sol.y[1] * alpha
    return pd.DataFrame({"time": t_eval, "Y": Y})


def compute_metrics(ts):
    """
    Compute performance metrics from the time series data.

    Args:
        ts (pd.DataFrame): Time series data with columns 'time' and 'Y'.
    Returns:
        dict: Computed metrics: dynamic_range, t50, overshoot.
    """
    if ts is None or len(ts) < 5:
        return {k: np.nan for k in ["dynamic_range", "t50", "overshoot"]}

    y = ts["Y"].to_numpy()
    t = ts["time"].to_numpy()

    # 1. Dynamic Range
    dyn_range = y.max() - y.min()

    # 2. t50 (Time to 50% change)
    y0 = y[0]
    y_ss = y[-1]
    total_change = y_ss - y0
    if abs(total_change) < 1e-6:
        t50 = np.nan
    else:
        target = y0 + 0.5 * total_change
        idx = np.where(np.sign(y - target) != np.sign(y[0] - target))[0]
        if len(idx) > 0:
            i = idx[0]
            t1, t2 = t[i - 1], t[i]
            y1, y2 = y[i - 1], y[i]
            if y2 != y1:
                t50 = t1 + (target - y1) * (t2 - t1) / (y2 - y1)
            else:
                t50 = t1
        else:
            t50 = np.nan

    # 3. Overshoot
    tail_mean = y[int(0.8 * len(y)):].mean()
    overshoot = y.max() - tail_mean

    return {
        "dynamic_range": dyn_range,
        "t50": t50,
        "overshoot": overshoot
    }


# -----------------------------------------------------------------------------
# Data Loading
# -----------------------------------------------------------------------------
def load_real_constructs():
    try:
        hill = pd.read_csv(PARAM_PATH / "Hill_parameters.csv")
        up = pd.read_csv(PARAM_PATH / "half_times_upregulation.csv")
        alpha_df = pd.read_csv(PARAM_PATH / "alphamcherry.csv")

        if "halftime" in up.columns: up = up.rename(columns={"halftime": "t_up"})
        if "K" in hill.columns and "k" not in hill.columns: hill = hill.rename(columns={"K": "k"})

        df = hill.merge(up, on="plasmid", how="inner")

        if "plasmid" in alpha_df.columns:
            df = df.merge(alpha_df[["plasmid", "alpha"]], on="plasmid", how="left")
        else:
            df["alpha"] = float(alpha_df["alpha"].iloc[0])

        return df
    except Exception as e:
        print(f"[WARN] Could not load real constructs: {e}")
        return pd.DataFrame()


# -----------------------------------------------------------------------------
# Main Analysis
# -----------------------------------------------------------------------------
def main():
    print("--- Step 7: Design Selection & Mapping (Robust) ---")

    # 1. Load Simulated Design Space
    sim_file = PARAM_PATH / "design_space_scan_repression.csv"
    if not sim_file.exists():
        print(f"[ERROR] {sim_file} not found. Run step_4b first.")
        return

    sim_df = pd.read_csv(sim_file)
    # Drop rows where simulation failed (NaN metrics)
    sim_df = sim_df.dropna(subset=["dynamic_range", "t50"])

    print(f"Analyzing {len(sim_df)} valid simulated designs.")
    if len(sim_df) == 0:
        print("[ERROR] No valid simulations found. Check step_4b output.")
        return

    # 2. Ranking Strategy: Performance Score
    # We want HIGH Dynamic Range and LOW t50.
    # Score = Range / (t50 + 1.0)
    # (Adding 1.0 prevents division by zero/small numbers and smooths the score)
    sim_df["score"] = sim_df["dynamic_range"] / (sim_df["t50"] + 1.0)

    # Sort by Score Descending
    candidates = sim_df.sort_values("score", ascending=False)
    top_10 = candidates.head(10)

    # 3. Simulate Real Constructs for Context
    real_df = load_real_constructs()
    real_metrics = []

    if not real_df.empty:
        for _, row in real_df.iterrows():
            ts = simulate_construct(
                t_up=row["t_up"],
                K=row["k"],
                n=row["n"],
                alpha=row["alpha"]
            )
            m = compute_metrics(ts)
            m["plasmid"] = row["plasmid"]
            real_metrics.append(m)
        real_results = pd.DataFrame(real_metrics)
    else:
        real_results = pd.DataFrame()

    # 4. Visualization: Performance Map
    plt.figure(figsize=(9, 7))

    # Background: All designs (Grey)
    plt.scatter(sim_df["t50"], sim_df["dynamic_range"],
                c="lightgrey", alpha=0.5, s=15, label="Simulated Space")

    # Highlight: Top 10 Candidates (Blue)
    plt.scatter(top_10["t50"], top_10["dynamic_range"],
                c="#4DBBD5", s=50, marker="D", label="Top 10 Candidates")

    # Overlay: Real Constructs (Red)
    if not real_results.empty:
        r = real_results.dropna(subset=["t50", "dynamic_range"])
        plt.scatter(r["t50"], r["dynamic_range"],
                    c="#E64B35", s=100, edgecolors="black", label="Observed Constructs")

        for _, row in r.iterrows():
            plt.text(row["t50"] + 0.5, row["dynamic_range"],
                     str(row["plasmid"]), fontsize=9, fontweight="bold", color="#E64B35")

    plt.xlabel("Response Time ($t_{50}$, hours) [Lower is Faster]")
    plt.ylabel("Dynamic Range (Norm. Units) [Higher is Better]")
    plt.title("Design Space Map: Selecting High-Performance Tuners")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend(loc="upper right")
    plt.tight_layout()

    plot_file = OUT_PATH / "design_space_map.pdf"
    plt.savefig(plot_file, dpi=300)
    print(f"[PLOT] Saved map to {plot_file}")

    # 5. Save Output
    csv_out = PARAM_PATH / "candidate_selection_top10.csv"
    # Select readable columns
    cols = ["design_id", "K", "n", "t_up", "dynamic_range", "t50", "overshoot", "score"]
    top_10[cols].to_csv(csv_out, index=False)

    print("\n--- Top 5 Recommended Designs (Ranked by Score) ---")
    print(top_10[cols].head(5).to_string(index=False))
    print(f"\n[OUTPUT] Saved top 10 candidates to {csv_out}")


if __name__ == "__main__":
    main()

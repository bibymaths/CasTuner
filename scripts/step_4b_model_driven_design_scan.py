#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step 4b – Model-driven design scan for analog gene tuners.

This script performs the following steps:
1) Loads measured parameter ranges from CSV files.
2) Samples a design space over key parameters (K, n, t_up, alpha, delay_kd).
3) Simulates repression dynamics for each sampled design.
4) Computes design metrics: dynamic range, t10-90, t50, overshoot.
5) Saves the results to a CSV file.
"""

import os, math, warnings
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

warnings.filterwarnings("ignore")

PARAM_PATH = "parameters"
os.makedirs(PARAM_PATH, exist_ok=True)


# -----------------------------------------------------------------------------
# Load existing parameter tables
# -----------------------------------------------------------------------------
def _read_norm(path: str) -> pd.DataFrame:
    """
    Read a CSV file and raise FileNotFoundError if not found.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    return df


def load_measured_parameter_ranges() -> Dict[str, Tuple[float, float]]:
    """
    Load measured parameter tables and derive min/max ranges for each.

    Returns
    -------
    dict
        Keys: 'K', 'n', 't_up', 't_down', 'alpha', 'delay_rev', 'delay_kd'
        Values: (min, max) tuples.
    """
    # t_down (Rev half-times)
    t_down = _read_norm(os.path.join(PARAM_PATH, "half_times_downregulation.csv"))
    if "halftime" in t_down.columns and "t_down" not in t_down.columns:
        t_down = t_down.rename(columns={"halftime": "t_down"})
    # t_up (KD half-times)
    t_up = _read_norm(os.path.join(PARAM_PATH, "half_times_upregulation.csv"))
    if "halftime" in t_up.columns and "t_up" not in t_up.columns:
        t_up = t_up.rename(columns={"halftime": "t_up"})

    # Hill parameters (K, n)
    hills = _read_norm(os.path.join(PARAM_PATH, "Hill_parameters.csv"))
    if "K" in hills.columns and "k" not in hills.columns:
        hills = hills.rename(columns={"K": "k"})
    if "n" not in hills.columns:
        raise KeyError("Hill_parameters.csv must contain a column 'n'.")

    # alpha (mCherry decay)
    alpha_df = _read_norm(os.path.join(PARAM_PATH, "alphamcherry.csv"))
    if "alpha" not in alpha_df.columns:
        for alt in ("alpha_mcherry", "alphamcherry", "a"):
            if alt in alpha_df.columns:
                alpha_df = alpha_df.rename(columns={alt: "alpha"})
                break
    if "alpha" not in alpha_df.columns:
        raise KeyError("alphamcherry.csv must contain 'alpha' or equivalent.")

    alpha_vals = alpha_df["alpha"].to_numpy(float)

    # optional delays
    delay_rev_path = os.path.join(PARAM_PATH, "delays_derepression.csv")
    delay_kd_path = os.path.join(PARAM_PATH, "delays_repression.csv")

    def _delay_range(path: str, col_name: str) -> Tuple[float, float]:
        if not os.path.exists(path):
            return (0.0, 8.0)  # fallback wide default
        df = pd.read_csv(path)
        if col_name not in df.columns:
            num_cols = [c for c in df.columns if df[c].dtype != "O"]
            if not num_cols:
                return (0.0, 8.0)
            col_name = num_cols[0]
        vals = df[col_name].to_numpy(float)
        return float(vals.min()), float(vals.max())

    delay_rev_min, delay_rev_max = _delay_range(delay_rev_path, "d_rev")
    delay_kd_min, delay_kd_max = _delay_range(delay_kd_path, "d_rev")

    ranges = {
        "K": (float(hills["k"].min()), float(hills["k"].max())),
        "n": (float(hills["n"].min()), float(hills["n"].max())),
        "t_up": (float(t_up["t_up"].min()), float(t_up["t_up"].max())),
        "t_down": (float(t_down["t_down"].min()), float(t_down["t_down"].max())),
        "alpha": (float(alpha_vals.min()), float(alpha_vals.max())),
        "delay_rev": (delay_rev_min, delay_rev_max),
        "delay_kd": (delay_kd_min, delay_kd_max),
    }
    return ranges


# -----------------------------------------------------------------------------
# Repression ODE (KD)
# -----------------------------------------------------------------------------
def ode_rhs(t: float, y: np.ndarray,
            t_up: float, K: float, n: float, alpha: float) -> np.ndarray:
    """
    ODE right-hand side for repression dynamics:

    dR = beta - R*(ln2/t_up)
    dY = (K^n/(K^n + R^n)) - alpha*Y

    We set beta = ln(2)/t_up so that R rises with half-time t_up.

    Parameters
    ----------
    t : float
        Current time (not used, as ODE is time-invariant).
    y : np.ndarray
        Current state vector [R, Y].
    t_up, K, n, alpha : float
        Model parameters.

    Returns
    -------
    np.ndarray
        Derivatives [dR, dY].
    """
    R, Y = y
    t_up = max(float(t_up), 1e-6)
    K = max(float(K), 1e-12)
    n = max(float(n), 1e-6)
    alpha = max(float(alpha), 1e-12)

    beta = math.log(2.0) / t_up
    R_pos = max(R, 0.0)
    Kn = K ** n
    Rn = R_pos ** n

    dR = beta - R_pos * (math.log(2.0) / t_up)
    dY = (Kn / (Kn + Rn)) - alpha * Y
    return np.array([dR, dY], dtype=float)


def simulate_repression(t_up: float, K: float, n: float, alpha: float,
                        delay: float,
                        t_end: float = 72.0,
                        dt: float = 0.05) -> pd.DataFrame:
    """
    Simulate repression trajectory Y(t) with a display delay (time shift).

    Parameters
    ----------
    t_up, K, n, alpha : float
        Model parameters.
    delay : float
        Display-only shift on the t-axis (>=0).
    t_end : float
        End time (hours) for simulation.
    dt : float
        Time step for output.

    Returns
    -------
    pd.DataFrame with columns: time, R, Y
    """
    t_eval = np.arange(0.0, t_end + dt / 2, dt)
    # Initial condition: R(0)=0, Y(0)=1/alpha → alpha*Y starts at 1 (matching step_3)
    y0 = np.array([0.0, 1.0 / max(alpha, 1e-12)], dtype=float)

    sol = solve_ivp(
        fun=lambda t, y: ode_rhs(t, y, t_up=t_up, K=K, n=n, alpha=alpha),
        t_span=(t_eval[0], t_eval[-1]),
        y0=y0,
        t_eval=t_eval,
        method="LSODA",
        rtol=1e-7,
        atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")

    R = sol.y[0, :]
    Y = sol.y[1, :] * alpha  # scale back (Y*alpha)
    df = pd.DataFrame({"time": t_eval + delay, "R": R, "Y": Y})
    return df


# -----------------------------------------------------------------------------
# Design metrics
# -----------------------------------------------------------------------------
def compute_design_metrics(ts: pd.DataFrame) -> Dict[str, float]:
    """
    Given a time series (time, Y), compute dynamic range and timing metrics.

    Parameters
    ----------
    ts : pd.DataFrame
        Time series with columns: time, Y

    Returns
    -------
    dict
        Keys: dynamic_range, t10_90, t50, overshoot
    """
    t = ts["time"].to_numpy(float)
    y = ts["Y"].to_numpy(float)

    if len(t) < 3:
        return {
            "dynamic_range": np.nan,
            "t10_90": np.nan,
            "t50": np.nan,
            "overshoot": np.nan,
        }

    y_min = float(y.min())
    y_max = float(y.max())
    dynamic_range = y_max - y_min

    # Define initial and final levels from early and late time windows
    y0 = float(y[0])
    y_ss = float(y[-1])
    total_change = y_ss - y0
    if abs(total_change) < 1e-9:
        t10_90 = np.nan
        t50 = np.nan
    else:
        y10 = y0 + 0.10 * total_change
        y90 = y0 + 0.90 * total_change
        y50 = y0 + 0.50 * total_change

        def _t_at_level(level):
            idx = np.where(np.sign(y - level) != np.sign(y[0] - level))[0]
            if len(idx) == 0:
                return np.nan
            i = idx[0]
            if i == 0:
                return float(t[0])
            # linear interpolation
            t1, t2 = t[i - 1], t[i]
            y1, y2 = y[i - 1], y[i]
            if y2 == y1:
                return float(t2)
            return float(t1 + (level - y1) * (t2 - t1) / (y2 - y1))

        t10 = _t_at_level(y10)
        t90 = _t_at_level(y90)
        t50 = _t_at_level(y50)
        t10_90 = t90 - t10 if (np.isfinite(t10) and np.isfinite(t90)) else np.nan

    # Overshoot relative to late-time mean
    y_tail = y[int(0.8 * len(y)):]
    y_tail_mean = float(y_tail.mean())
    overshoot = float(y.max() - y_tail_mean)

    return {
        "dynamic_range": dynamic_range,
        "t10_90": t10_90,
        "t50": t50,
        "overshoot": overshoot,
    }


# -----------------------------------------------------------------------------
# Design-space sampling
# -----------------------------------------------------------------------------
def sample_design_space(ranges: Dict[str, Tuple[float, float]],
                        n_samples: int = 500,
                        seed: int = 13) -> pd.DataFrame:
    """
    Sample a design space over (K, n, t_up, t_down, delay_kd) using uniform
    sampling within measured min/max ranges, slightly padded.

    Parameters
    ----------
    ranges : dict
        Keys: 'K', 'n', 't_up', 't_down', 'alpha', 'delay_kd'
        Values: (min, max) tuples.
    n_samples : int
        Number of samples to generate.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame with columns: design_id, K, n, t_up, t_down, alpha, delay_kd

    """
    rng = np.random.default_rng(seed)

    def pad(min_val, max_val, frac=0.25):
        span = max_val - min_val
        if span <= 0:
            span = max(abs(max_val), 1.0)
        return max(min_val * (1 - frac), 1e-6), max_val * (1 + frac)

    K_min, K_max = pad(*ranges["K"])
    n_min, n_max = pad(*ranges["n"])
    t_up_min, t_up_max = pad(*ranges["t_up"])
    t_down_min, t_down_max = pad(*ranges["t_down"])
    alpha_min, alpha_max = pad(*ranges["alpha"])
    delay_kd_min, delay_kd_max = pad(*ranges["delay_kd"])

    rows = []
    for i in range(n_samples):
        K = rng.uniform(K_min, K_max)
        n = rng.uniform(n_min, n_max)
        t_up = rng.uniform(t_up_min, t_up_max)
        t_down = rng.uniform(t_down_min, t_down_max)  # kept for completeness
        alpha = rng.uniform(alpha_min, alpha_max)
        delay_kd = rng.uniform(delay_kd_min, delay_kd_max)
        rows.append({
            "design_id": f"D{i:04d}",
            "K": K,
            "n": n,
            "t_up": t_up,
            "t_down": t_down,
            "alpha": alpha,
            "delay_kd": delay_kd,
        })
    return pd.DataFrame.from_records(rows)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    # 1) Load measured ranges
    ranges = load_measured_parameter_ranges()
    print("[ranges]", ranges)

    # 2) Sample design space
    designs = sample_design_space(ranges, n_samples=500, seed=13)

    # 3) Simulate each design and compute metrics
    metrics_records = []
    for _, row in designs.iterrows():
        K = float(row["K"])
        n = float(row["n"])
        t_up = float(row["t_up"])
        alpha = float(row["alpha"])
        delay = float(row["delay_kd"])

        try:
            ts = simulate_repression(t_up=t_up, K=K, n=n, alpha=alpha, delay=delay,
                                     t_end=72.0, dt=0.05)
            metrics = compute_design_metrics(ts)
        except Exception as e:
            # Mark design as failed
            metrics = {
                "dynamic_range": np.nan,
                "t10_90": np.nan,
                "t50": np.nan,
                "overshoot": np.nan,
            }

        rec = {"design_id": row["design_id"],
               "K": K, "n": n, "t_up": t_up, "t_down": float(row["t_down"]),
               "alpha": alpha, "delay_kd": delay}
        rec.update(metrics)
        metrics_records.append(rec)

    design_df = pd.DataFrame.from_records(metrics_records)
    out_path = os.path.join(PARAM_PATH, "design_space_scan_repression.csv")
    design_df.to_csv(out_path, index=False)
    print(f"[write] design_space_scan_repression.csv shape={design_df.shape}")


if __name__ == "__main__":
    main()

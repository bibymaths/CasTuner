#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step 4a – Single-cell noise and hierarchical summaries for CasTuner constructs.

This script processes single-cell time-course flow cytometry data to compute
noise metrics (mean, variance, CV²) for BFP and mCherry expression. It applies
gating identical to prior steps, subtracts NFC background, and summarizes
noise metrics per (plasmid, experiment, replicate, time) group. Finally, it
computes hierarchical summaries per construct using a simple normal–normal
partial pooling approach.

Outputs:
  - parameters/single_cell_noise_timeseries.csv : Per-(plasmid, exp, rep, time) noise metrics.
  - parameters/single_cell_noise_hierarchical.csv : Hierarchical summaries per construct.
"""

import os, re, glob, warnings
from typing import Tuple

import numpy as np
import pandas as pd
from FlowCytometryTools import FCMeasurement

warnings.filterwarnings("ignore")

# -----------------------------------------------------------------------------
# Paths and config
# -----------------------------------------------------------------------------
OUT_PATH = "plots"
PARAM_PATH = "parameters"

# FCS files stay at the project root
FCS_TC_DIR = os.path.join("fcs_files", "time-course_data")
FCS_NFC_DIR = os.path.join("fcs_files", "NFC")

os.makedirs(OUT_PATH, exist_ok=True)
os.makedirs(PARAM_PATH, exist_ok=True)

# Channels
CH_FSC_A = "FSC-A"
CH_SSC_A = "SSC-A"
CH_FSC_H = "FSC-H"
CH_BFP = "BV421-A"
CH_mCh = "PE-A"

# Gating bounds (same as Steps 1–3)
BOUND_MIN = {CH_FSC_A: 0.4e5, CH_SSC_A: 0.20e5}
BOUND_MAX = {CH_FSC_A: 2.0e5, CH_SSC_A: 1.3e5}
SINGLET_RATIO_LOW, SINGLET_RATIO_HIGH = 0.85, 1.15

_SPLIT = re.compile(r"_")

# -----------------------------------------------------------------------------
# Gating + NFC background
# -----------------------------------------------------------------------------
def apply_boundary_gate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply boundary gate on FSC-A and SSC-A.
    If no events pass, return raw data.

    Parameters
    ----------
    df : pd.DataFrame
        Flow cytometry events.

    Returns
    -------
    pd.DataFrame
        Gated events.
    """
    m = (
            (df[CH_FSC_A] >= BOUND_MIN[CH_FSC_A]) & (df[CH_FSC_A] <= BOUND_MAX[CH_FSC_A]) &
            (df[CH_SSC_A] >= BOUND_MIN[CH_SSC_A]) & (df[CH_SSC_A] <= BOUND_MAX[CH_SSC_A])
    )
    out = df.loc[m]
    if out.empty:
        # fall back to raw
        return df
    return out


def apply_singlet_gate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply singlet gate based on FSC-H / FSC-A ratio.
    If no events pass, return raw data.

    Parameters
    ----------
    df : pd.DataFrame
        Flow cytometry events.

    Returns
    -------
    pd.DataFrame
        Gated events.
    """
    ratio = df[CH_FSC_H] / df[CH_FSC_A].replace(0, np.nan)
    m = (ratio >= SINGLET_RATIO_LOW) & (ratio <= SINGLET_RATIO_HIGH)
    out = df.loc[m]
    if out.empty:
        return df
    return out


def compute_nfc_background(nfc_dir: str) -> Tuple[float, float]:
    """
    Compute NFC background medians for BFP and mCherry.
    Parameters
    ----------
    nfc_dir : str
        Directory containing NFC .fcs files.
    Returns
    -------
    Tuple[float, float]
        (mBFP_neg, mmCherry_neg)
    """
    files = sorted(glob.glob(os.path.join(nfc_dir, "*.fcs")))
    if not files:
        raise FileNotFoundError(f"No NFC .fcs files in {nfc_dir}")
    rows = []
    for fpath in files:
        dat = FCMeasurement(ID=os.path.basename(fpath), datafile=fpath).data
        gated = apply_singlet_gate(apply_boundary_gate(dat))
        med = gated.median(numeric_only=True)
        rows.append(med)
    df = pd.DataFrame(rows).reset_index(drop=True)
    mBFP_neg = float(df.loc[:2, CH_BFP].mean())
    mmCherry_neg = float(df.loc[:2, CH_mCh].mean())
    return mBFP_neg, mmCherry_neg


def parse_timecourse_name(name: str):
    """
    Extract (plasmid, exp, rep, time) from filename stem.

    Expected pattern:
      <prefix>_<prefix>_<plasmid>_<exp>_<rep>_<time>[_...]
    """
    parts = _SPLIT.split(name)
    plasmid = parts[2] if len(parts) > 2 else ""
    exp = parts[3] if len(parts) > 3 else ""
    rep = parts[4] if len(parts) > 4 else ""
    time_s = parts[5] if len(parts) > 5 else ""
    try:
        t = float(time_s)
    except Exception:
        m = re.search(r"(\d+(?:\.\d+)?)", time_s)
        t = float(m.group(1)) if m else np.nan
    return plasmid, exp, rep, t


# -----------------------------------------------------------------------------
# Single-cell time-course load
# -----------------------------------------------------------------------------
def load_single_cell_timecourse(exp_filter=None) -> pd.DataFrame:
    """
    Load single-cell events from time-course FCS, gate, subtract NFC, and attach
    metadata (plasmid, exp, rep, time).

    Parameters
    ----------
    exp_filter : {None, "Rev", "KD"}, optional
        If provided, restrict to that experiment type.

    Returns
    -------
    pd.DataFrame
        Columns:
          plasmid, exp, rep, time, BFP, mCherry
    """
    mBFP_neg, mmCherry_neg = compute_nfc_background(FCS_NFC_DIR)
    files = sorted(glob.glob(os.path.join(FCS_TC_DIR, "*.fcs")))
    if not files:
        raise FileNotFoundError(f"No .fcs files in {FCS_TC_DIR}")

    rows = []
    for fpath in files:
        stem = os.path.splitext(os.path.basename(fpath))[0]
        plasmid, exp, rep, t = parse_timecourse_name(stem)
        if exp_filter is not None and exp != exp_filter:
            continue

        dat = FCMeasurement(ID=stem, datafile=fpath).data
        gated = apply_singlet_gate(apply_boundary_gate(dat))

        if CH_BFP not in gated.columns or CH_mCh not in gated.columns:
            continue  # skip malformed files

        # Background subtraction
        bfp = gated[CH_BFP].to_numpy() - mBFP_neg
        mch = gated[CH_mCh].to_numpy() - mmCherry_neg

        # Drop negative values (optional but common for log/noise metrics)
        mask = (bfp > 0) & (mch > 0)
        bfp = bfp[mask]
        mch = mch[mask]

        n = len(bfp)
        if n == 0:
            continue

        df_local = pd.DataFrame({
            "plasmid": plasmid,
            "exp": exp,
            "rep": rep,
            "time": float(t),
            "BFP": bfp,
            "mCherry": mch,
        })
        rows.append(df_local)

    if not rows:
        raise RuntimeError("No gated events left after filtering – check FCS layout.")
    out = pd.concat(rows, ignore_index=True)
    return out


# -----------------------------------------------------------------------------
# Per-group noise metrics
# -----------------------------------------------------------------------------
def summarize_noise_per_group(events: pd.DataFrame) -> pd.DataFrame:
    """
    Compute noise metrics per (plasmid, exp, rep, time).

    Parameters
    ----------
    events : pd.DataFrame
        Single-cell events with columns:
          plasmid, exp, rep, time, BFP, mCherry

    Returns
    -------
    pd.DataFrame with columns:
      plasmid, exp, rep, time,
      mean_BFP, var_BFP, cv2_BFP,
      mean_mCherry, var_mCherry, cv2_mCherry,
      n_cells
    """

    def _stats(x: pd.Series):
        """
        Compute mean, variance, CV² for a series.

        Parameters
        ----------
        x : pd.Series
            Numeric values.
        Returns
        -------
        Tuple[float, float, float]
            (mean, variance, CV²)
        """
        m = float(x.mean())
        v = float(x.var(ddof=1)) if len(x) > 1 else 0.0
        cv2 = (v / (m ** 2)) if m > 0 else np.nan
        return m, v, cv2

    records = []
    for (pl, ex, rep, t), sub in events.groupby(["plasmid", "exp", "rep", "time"], dropna=False):
        mB, vB, cv2B = _stats(sub["BFP"])
        mC, vC, cv2C = _stats(sub["mCherry"])
        records.append({
            "plasmid": pl,
            "exp": ex,
            "rep": rep,
            "time": t,
            "mean_BFP": mB,
            "var_BFP": vB,
            "cv2_BFP": cv2B,
            "mean_mCherry": mC,
            "var_mCherry": vC,
            "cv2_mCherry": cv2C,
            "n_cells": int(len(sub)),
        })
    return pd.DataFrame.from_records(records)


# -----------------------------------------------------------------------------
# Hierarchical summaries – simple normal–normal partial pooling
# -----------------------------------------------------------------------------
def hierarchical_summary(noise_df: pd.DataFrame) -> pd.DataFrame:
    """
    For each plasmid (and exp), compute pooled noise metrics and uncertainty.

    Parameters
    ----------
    noise_df : pd.DataFrame
        Per-(plasmid, exp, rep, time) noise metrics with columns:
          plasmid, exp, rep, time,
          mean_BFP, var_BFP, cv2_BFP,
          mean_mCherry, var_mCherry, cv2_mCherry,
          n_cells

    Returns
    -------
    pd.DataFrame with columns:
      plasmid, exp, n_groups,
      mean_BFP, mean_BFP_se, mean_BFP_ci_low, mean_BFP_ci_high,
      mean_mCherry, mean_mCherry_se, mean_mCherry_ci_low, mean_mCherry_ci_high,
      cv2_BFP, cv2_BFP_se, cv2_BFP_ci_low, cv2_BFP_ci_high,
      cv2_mCherry, cv2_mCherry_se, cv2_mCherry_ci_low, cv2_mCherry_ci_high
    """
    metrics = ["cv2_BFP", "cv2_mCherry", "mean_BFP", "mean_mCherry"]
    rows = []

    for (pl, ex), sub in noise_df.groupby(["plasmid", "exp"]):
        n_groups = len(sub)
        if n_groups == 0:
            continue

        rec = {"plasmid": pl, "exp": ex, "n_groups": n_groups}
        for m in metrics:
            vals = sub[m].dropna().to_numpy()
            if len(vals) == 0:
                rec[m] = np.nan
                rec[m + "_se"] = np.nan
                rec[m + "_ci_low"] = np.nan
                rec[m + "_ci_high"] = np.nan
                continue
            mbar = float(vals.mean())
            s_between = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
            se = s_between / np.sqrt(max(len(vals), 1))
            ci_low = mbar - 1.96 * se
            ci_high = mbar + 1.96 * se
            rec[m] = mbar
            rec[m + "_se"] = se
            rec[m + "_ci_low"] = ci_low
            rec[m + "_ci_high"] = ci_high
        rows.append(rec)

    return pd.DataFrame.from_records(rows)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    # 1) Load Rev + KD single-cell events
    events = load_single_cell_timecourse(exp_filter=None)

    # 2) Per-group noise metrics
    noise_ts = summarize_noise_per_group(events)
    noise_ts.to_csv(
        os.path.join(PARAM_PATH, "single_cell_noise_timeseries.csv"),
        index=False
    )
    print(f"[write] single_cell_noise_timeseries.csv  shape={noise_ts.shape}")

    # 3) Hierarchical summaries per construct
    hier = hierarchical_summary(noise_ts)
    hier.to_csv(
        os.path.join(PARAM_PATH, "single_cell_noise_hierarchical.csv"),
        index=False
    )
    print(f"[write] single_cell_noise_hierarchical.csv shape={hier.shape}")


if __name__ == "__main__":
    main()

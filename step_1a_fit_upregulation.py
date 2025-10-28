#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
KD (up-regulation) time-course: fit exponential half-times on min-max scaled BFP.

Requires (Py3.9-compatible):
  - numpy<2.0, pandas, scipy
  - FlowCytometryTools, flowio, flowutils
  - plotnine

Folder structure:
  fcs_files/
    NFC/*.fcs
    time-course_data/*.fcs

Outputs:
  plots/KD_KRAB-Split-dCas9_fitting.pdf
  plots/KD_dCas9_fitting.pdf
  plots/KD_KRAB-dCas9_fitting.pdf
  plots/KD_HDAC4-dCas9_fitting.pdf
  plots/KD_CasRx_fitting.pdf
  parameters/half_times_upregulation.csv
"""

import os
import re
import glob
import math
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from FlowCytometryTools import FCMeasurement
from plotnine import (
    ggplot, aes, geom_point, geom_line, labs,
    scale_y_continuous, coord_cartesian, theme_classic
)
warnings.filterwarnings("ignore")

# ----------------------------
# Config / constants
# ----------------------------
OUT_PATH = "plots"
PARAM_PATH = "parameters"
os.makedirs(OUT_PATH, exist_ok=True)
os.makedirs(PARAM_PATH, exist_ok=True)

# Channels
CH_FSC_A = "FSC-A"
CH_SSC_A = "SSC-A"
CH_FSC_H = "FSC-H"
CH_BFP   = "BV421-A"  # tagBFP
CH_mCh   = "PE-A"     # mCherry

# Boundary gate (matches R numeric limits)
BOUND_MIN = {CH_FSC_A: 0.4e5, CH_SSC_A: 0.20e5}
BOUND_MAX = {CH_FSC_A: 2.0e5, CH_SSC_A: 1.3e5}

# Singlet gate (FSC-H ~ FSC-A)
SINGLET_RATIO_LOW, SINGLET_RATIO_HIGH = 0.85, 1.15

# Plot sizing (inches)
PLOT_W = 1.5 * 1.618
PLOT_H = 1.5

# Aesthetics
POINT_COLOR = "#4DBBD5FF"
THEME = theme_classic()

# ----------------------------
# Gating helpers
# ----------------------------
def apply_boundary_gate(df: pd.DataFrame) -> pd.DataFrame:
    m = (
        (df[CH_FSC_A] >= BOUND_MIN[CH_FSC_A]) & (df[CH_FSC_A] <= BOUND_MAX[CH_FSC_A]) &
        (df[CH_SSC_A] >= BOUND_MIN[CH_SSC_A]) & (df[CH_SSC_A] <= BOUND_MAX[CH_SSC_A])
    )
    return df.loc[m]

def apply_singlet_gate(df: pd.DataFrame) -> pd.DataFrame:
    ratio = df[CH_FSC_H] / (df[CH_FSC_A].replace(0, np.nan))
    m = (ratio >= SINGLET_RATIO_LOW) & (ratio <= SINGLET_RATIO_HIGH)
    return df.loc[m]

def median_channels_for_file(fpath: str, mBFP_neg: float = 0.0, mmCherry_neg: float = 0.0) -> pd.Series:
    sample = FCMeasurement(ID=os.path.basename(fpath), datafile=fpath).data
    for ch in (CH_FSC_A, CH_SSC_A, CH_FSC_H, CH_BFP, CH_mCh):
        if ch not in sample.columns:
            raise ValueError(f"Channel '{ch}' not found in: {fpath}")

    gated = apply_boundary_gate(sample)
    if len(gated) == 0:
        gated = sample
    singlets = apply_singlet_gate(gated)
    if len(singlets) == 0:
        singlets = gated

    # Apply subtraction BEFORE median, if values are non-zero
    if mBFP_neg != 0.0 or mmCherry_neg != 0.0:
        singlets_sub = singlets.copy()
        singlets_sub[CH_BFP] = singlets_sub[CH_BFP] - mBFP_neg
        singlets_sub[CH_mCh] = singlets_sub[CH_mCh] - mmCherry_neg
        s = singlets_sub.median(numeric_only=True)
    else:
        # Original behavior for loading NFC files
        s = singlets.median(numeric_only=True)

    s["__filename"] = os.path.splitext(os.path.basename(fpath))[0]
    return s

def load_flowset_medians(folder: str, mBFP_neg: float = 0.0, mmCherry_neg: float = 0.0) -> pd.DataFrame:
    files = sorted(glob.glob(os.path.join(folder, "*.fcs")))
    if not files:
        raise FileNotFoundError(f"No FCS files found in: {folder}")
    # Pass background values (or defaults of 0) to the helper
    return pd.DataFrame([median_channels_for_file(f, mBFP_neg, mmCherry_neg) for f in files])

# ----------------------------
# NFC background (mean of first 3 medians)
# ----------------------------
def compute_nfc_background(nfc_dir: str):
    df = load_flowset_medians(nfc_dir)
    mBFP_neg = float(df.iloc[:3][CH_BFP].mean())
    mmCherry_neg = float(df.iloc[:3][CH_mCh].mean())
    return mBFP_neg, mmCherry_neg

# ----------------------------
# Parse filename: tokens
# R: separate(rowname, c(NA,NA,"plasmid","exp","rep","time",NA,NA), "_")
# -> plasmid=parts[2], exp=parts[3], rep=parts[4], time=parts[5]
# ----------------------------
FILENAME_SPLIT_RE = re.compile(r"_")

def parse_timecourse_name(name: str):
    parts = FILENAME_SPLIT_RE.split(name)
    plasmid = parts[2] if len(parts) > 2 else ""
    exp     = parts[3] if len(parts) > 3 else ""
    rep     = parts[4] if len(parts) > 4 else ""
    time_s  = parts[5] if len(parts) > 5 else ""
    try:
        time = float(time_s)
    except Exception:
        m = re.search(r"(\d+(?:\.\d+)?)", time_s)
        time = float(m.group(1)) if m else np.nan
    return plasmid, exp, rep, time

# ----------------------------
# Build KD dataset
# ----------------------------
def load_kd_timecourse(mBFP_neg: float, mmCherry_neg: float) -> pd.DataFrame:
    med = load_flowset_medians("fcs_files/time-course_data").copy()

    med = load_flowset_medians(
        "fcs_files/time-course_data",
        mBFP_neg,
        mmCherry_neg
    ).copy()

    parsed = med["__filename"].apply(parse_timecourse_name).tolist()
    parsed_df = pd.DataFrame(parsed, columns=["plasmid", "exp", "rep", "time"])
    df = pd.concat([med[[CH_BFP, CH_mCh]], parsed_df], axis=1)

    # keep KD experiments
    df = df[df["exp"] == "KD"].copy()
    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df = df.dropna(subset=["time"])
    # plasmid order (factor levels in R)
    cat = pd.CategoricalDtype(["SP430", "SP428", "SP430ABA", "SP427", "SP411"], ordered=True)
    df["plasmid"] = df["plasmid"].astype(cat)
    return df

# ----------------------------
# Normalization (KD):
# mean.final = mean(BFP | time>10)
# mean.init  = mean(BFP | time==0)
# norm.bfp = (BFP - mean.init) / (mean.final - mean.init)
# ----------------------------
def add_minmax_norm_kd(df: pd.DataFrame) -> pd.DataFrame:
    mean_final = (
        df.query("time > 10")
          .groupby("plasmid")[CH_BFP]
          .mean().rename("mean.final").reset_index()
    )
    mean_init = (
        df.query("time == 0")
          .groupby("plasmid")[CH_BFP]
          .mean().rename("mean.init").reset_index()
    )
    out = df.merge(mean_final, on="plasmid", how="left").merge(mean_init, on="plasmid", how="left")
    denom = (out["mean.final"] - out["mean.init"])
    out["norm.bfp"] = (out[CH_BFP] - out["mean.init"]) / denom.replace(0, np.nan)
    return out

# ----------------------------
# Model & fitting (KD):
# y = 1 - exp(-t * ln(2) / t_half)
# ----------------------------
def exp_rise_to_one(t, t_half):
    return 1.0 - np.exp(-t * (math.log(2.0) / t_half))

def fit_half_time_rise(t, y, start=0.8):
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    m = np.isfinite(t) & np.isfinite(y)
    t, y = t[m], y[m]
    if len(t) < 3:
        raise RuntimeError("Not enough points for growth fit.")
    popt, pcov = curve_fit(exp_rise_to_one, t, y, p0=[float(start)], maxfev=20000)
    se = float(np.sqrt(np.diag(pcov))[0]) if pcov is not None else np.nan
    return float(popt[0]), se

# ----------------------------
# Plot helper (match R look)
# ----------------------------
def save_kd_plot(df, t_half, out_pdf):
    tmin, tmax = float(df["time"].min()), float(df["time"].max())
    curve = pd.DataFrame({"time": np.linspace(tmin, tmax, 200)})
    curve["fit"] = exp_rise_to_one(curve["time"], t_half)

    p = (
        ggplot(df, aes("time", "norm.bfp"))
        + geom_point(size=0.4, alpha=0.7, color=POINT_COLOR)
        + geom_line(data=curve, mapping=aes(x="time", y="fit"), color="black")
        + coord_cartesian(ylim=(-0.15, 1.2))
        + scale_y_continuous(breaks=[0, 0.25, 0.5, 0.75, 1.0])
        + labs(x="Time (hours)", y="tagBFP (% of final)")
        + THEME
    )
    p.save(os.path.join(OUT_PATH, out_pdf), width=PLOT_W, height=PLOT_H, units="in")

# ----------------------------
# Main
# ----------------------------
def main():
    # 1) NFC background
    mBFP_neg, mmCherry_neg = compute_nfc_background("fcs_files/NFC")

    # 2) Load KD time-course with background-subtracted medians
    kd = load_kd_timecourse(mBFP_neg, mmCherry_neg)

    # 3) Min–max normalize per plasmid (KD definition)
    kd = add_minmax_norm_kd(kd)

    # 4) Fit per plasmid & save plots
    targets = [
        ("SP430ABA", "KD_KRAB-Split-dCas9_fitting.pdf", "SP430A"),
        ("SP430",    "KD_dCas9_fitting.pdf",            "SP430"),
        ("SP428",    "KD_KRAB-dCas9_fitting.pdf",       "SP428"),
        ("SP427",    "KD_HDAC4-dCas9_fitting.pdf",      "SP427"),
        ("SP411",    "KD_CasRx_fitting.pdf",            "SP411"),
    ]

    rows = []
    for plasmid, pdfname, label in targets:
        dfp = kd[kd["plasmid"] == plasmid].copy()
        if dfp.empty:
            print(f"[WARN] No KD data for {plasmid}; skipping.")
            continue

        t_half, se = fit_half_time_rise(dfp["time"], dfp["norm.bfp"], start=0.8)
        print(f"{label}: t1/2 = {t_half:.6g}  (SE={se:.3g})")

        rows.append({"plasmid": label, "halftime": t_half, "se": se})
        save_kd_plot(dfp, t_half, pdfname)

    # 5) Save half-times
    if rows:
        ht = pd.DataFrame(rows, columns=["plasmid", "halftime", "se"])
        ht.to_csv(os.path.join(PARAM_PATH, "half_times_upregulation.csv"), index=False)
        print("\nSaved:", os.path.join(PARAM_PATH, "half_times_upregulation.csv"))
        print(ht)
    else:
        print("[WARN] No half-times estimated.")

if __name__ == "__main__":
    main()

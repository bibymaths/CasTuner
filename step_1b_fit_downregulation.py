#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Reverse (Rev) time-course: fit exponential half-times on min-max scaled BFP.

Requires (Py3.9-compatible):
  - numpy<2.0, pandas, scipy
  - FlowCytometryTools, flowio, flowutils
  - plotnine

Folder structure:
  fcs_files/
    NFC/*.fcs
    time-course_data/*.fcs

Outputs:
  plots/REV_SP430ABA_fitting.pdf
  plots/REV_dCas9_fitting.pdf
  plots/REV_KRAB-dCas9_fitting.pdf
  plots/REV_HDAC4-dCas9_fitting.pdf
  plots/REV_CasRx_fitting.pdf
  parameters/half_times_downregulation.csv
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

# Plot sizing (inches), like set_panel_size in R
PLOT_W = 1.5 * 1.618
PLOT_H = 1.5

# Aesthetics
POINT_COLOR = "#4DBBD5FF"  # NPG-like
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

def median_channels_for_file(fpath: str) -> pd.Series:
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

    s = singlets.median(numeric_only=True)
    s["__filename"] = os.path.splitext(os.path.basename(fpath))[0]
    return s

def load_flowset_medians(folder: str) -> pd.DataFrame:
    files = sorted(glob.glob(os.path.join(folder, "*.fcs")))
    if not files:
        raise FileNotFoundError(f"No FCS files found in: {folder}")
    return pd.DataFrame([median_channels_for_file(f) for f in files])

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
# Rowname split in R: NA, NA, plasmid, exp, rep, time, NA, NA
# -> plasmid = parts[2], exp = parts[3], rep = parts[4], time = parts[5]
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
        # fallback: extract trailing digits from token
        m = re.search(r"(\d+(?:\.\d+)?)", time_s)
        time = float(m.group(1)) if m else np.nan
    return plasmid, exp, rep, time

# ----------------------------
# Build Rev dataset
# ----------------------------
def load_rev_timecourse(mBFP_neg: float, mmCherry_neg: float) -> pd.DataFrame:
    med = load_flowset_medians("fcs_files/time-course_data").copy()

    # background-subtract medians (equivalent to per-event linear transform)
    med[CH_BFP] = med[CH_BFP] - mBFP_neg
    med[CH_mCh] = med[CH_mCh] - mmCherry_neg

    parsed = med["__filename"].apply(parse_timecourse_name).tolist()
    parsed_df = pd.DataFrame(parsed, columns=["plasmid", "exp", "rep", "time"])
    df = pd.concat([med[[CH_BFP, CH_mCh]], parsed_df], axis=1)

    # keep Rev experiments
    df = df[df["exp"] == "Rev"].copy()
    # ensure numeric
    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df = df.dropna(subset=["time"])
    return df

# ----------------------------
# Normalization: min-max scaled BFP per plasmid
# mean.final = mean(BFP) for time>10; mean.init = mean(BFP) for time==0
# norm.bfp = (BFP - mean.final) / (mean.init - mean.final)
# ----------------------------
def add_minmax_norm(df: pd.DataFrame) -> pd.DataFrame:
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
    denom = (out["mean.init"] - out["mean.final"])
    out["norm.bfp"] = (out[CH_BFP] - out["mean.final"]) / denom.replace(0, np.nan)
    return out

# ----------------------------
# Model & fitting: y = exp(-t * ln(2) / t_half)
# ----------------------------
def exp_decay(t, t_half):
    return np.exp(-t * (math.log(2.0) / t_half))

def fit_half_time(t, y, start=0.1):
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    m = np.isfinite(t) & np.isfinite(y)
    t, y = t[m], y[m]
    if len(t) < 3:
        raise RuntimeError("Not enough points for decay fit.")
    popt, pcov = curve_fit(exp_decay, t, y, p0=[float(start)], maxfev=20000)
    # standard error for t_half
    se = float(np.sqrt(np.diag(pcov))[0]) if pcov is not None else np.nan
    return float(popt[0]), se

# ----------------------------
# Plot helper (match R look)
# ----------------------------
def save_rev_plot(df, t_half, out_pdf):
    tmin, tmax = float(df["time"].min()), float(df["time"].max())
    curve = pd.DataFrame({"time": np.linspace(tmin, tmax, 200)})
    curve["fit"] = exp_decay(curve["time"], t_half)

    p = (
        ggplot(df, aes("time", "norm.bfp"))
        + geom_point(size=0.4, alpha=0.7, color=POINT_COLOR)
        + geom_line(data=curve, mapping=aes(x="time", y="fit"), color="black")
        + coord_cartesian(ylim=(-0.15, 1.4))
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

    # 2) Load Rev time-course with background-subtracted medians
    rev = load_rev_timecourse(mBFP_neg, mmCherry_neg)

    # 3) Min–max normalize per plasmid
    rev = add_minmax_norm(rev)

    # 4) Fit per required plasmid & save plots
    targets = [
        ("SP430ABA", "REV_SP430ABA_fitting.pdf", "SP430A"),
        ("SP430",    "REV_dCas9_fitting.pdf",    "SP430"),
        ("SP428",    "REV_KRAB-dCas9_fitting.pdf","SP428"),
        ("SP427",    "REV_HDAC4-dCas9_fitting.pdf","SP427"),
        ("SP411",    "REV_CasRx_fitting.pdf",    "SP411"),
    ]

    rows = []
    for plasmid, pdfname, label in targets:
        dfp = rev[rev["plasmid"] == plasmid].copy()
        if dfp.empty:
            print(f"[WARN] No Rev data for {plasmid}; skipping.")
            continue

        t_half, se = fit_half_time(dfp["time"], dfp["norm.bfp"], start=0.1)
        print(f"{label}: t1/2 = {t_half:.6g}  (SE={se:.3g})")

        rows.append({"plasmid": label, "halftime": t_half, "se": se})
        save_rev_plot(dfp, t_half, pdfname)

    # 5) Save half-times
    if rows:
        ht = pd.DataFrame(rows, columns=["plasmid", "halftime", "se"])
        ht.to_csv(os.path.join(PARAM_PATH, "half_times_downregulation.csv"), index=False)
        print("\nSaved:", os.path.join(PARAM_PATH, "half_times_downregulation.csv"))
        print(ht)
    else:
        print("[WARN] No half-times estimated.")

if __name__ == "__main__":
    main()

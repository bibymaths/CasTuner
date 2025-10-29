#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------------
#
# This script was independently developed by Abhinav Mishra as a Python
# reimplementation of the R-based CasTuner analysis pipeline, originally used in:
#   “CasTuner: a degron and CRISPR/Cas-based toolkit for analog tuning of
#    endogenous gene expression” — Gemma Noviello, Rutger A. F. Gjaltema,
#    and Edda G. Schulz.
#
# Purpose:
#   This code is part of a prospective research extension designed to demonstrate
#   technical capability and methodological thinking relevant to a potential PhD
#   position in the Schulz group. I am not currently affiliated with the group;
#   this work is intended to showcase independent reproducibility, modeling rigor,
#   and the ability to extend CasTuner’s first analytical stage in Python.
#
# Original R workflow (for reference):
#   Step 1a  Fit tagBFP upregulation dynamics after dTAG-13 withdrawal
#            (file: "Fitting tagBFP upregulation.R")
#   Step 1b  Fit tagBFP downregulation dynamics after dTAG-13 addition
#   Step 1c  Fit steady-state dose–response (Hill curves)
#   Step 2   ODE simulations for mCherry derepression ("ODE_REV.R")
#   Step 3   ODE simulations for mCherry repression ("ODE_KD.R")
#
# This Python implementation reproduces Step 1b:
#   • Loads .fcs files, applies the same gating logic and background subtraction
#   • Performs per-plasmid min–max normalization for downregulation
#   • Fits an exponential decay model using scipy.optimize.curve_fit
#   • Generates matching plots and parameter tables for validation
#
# All thresholds, formulas, and transformations are preserved from the original
# R code to maintain reproducibility and comparability.
# -----------------------------------------------------------------------------

"""
Reverse (Rev) time-course analysis.

This module loads flow cytometry FCS files, applies boundary and singlet gates,
computes per-file medians, subtracts NFC-derived background, constructs a Rev-only
time-course dataset, performs per-plasmid min–max normalization on BFP, fits an
exponential decay model to estimate half-times (t1/2), generates diagnostic plots,
and writes the estimated parameters to CSV.

Outputs:
  - plots/REV_SP430ABA_fitting.pdf
  - plots/REV_dCas9_fitting.pdf
  - plots/REV_KRAB-dCas9_fitting.pdf
  - plots/REV_HDAC4-dCas9_fitting.pdf
  - plots/REV_CasRx_fitting.pdf
  - parameters/half_times_downregulation.csv
"""


import os, re, glob, math, warnings
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
OUT_PATH = "plots"                                 # directory for plot outputs
PARAM_PATH = "parameters"                          # directory for parameter CSVs
os.makedirs(OUT_PATH, exist_ok=True)               # ensure plot dir exists
os.makedirs(PARAM_PATH, exist_ok=True)             # ensure parameter dir exists

# Channels
CH_FSC_A = "FSC-A"                                 # forward scatter area
CH_SSC_A = "SSC-A"                                 # side scatter area
CH_FSC_H = "FSC-H"                                 # forward scatter height
CH_BFP   = "BV421-A"                               # tagBFP channel name
CH_mCh   = "PE-A"                                  # mCherry channel name

# Boundary gate (matches R numeric limits)
BOUND_MIN = {CH_FSC_A: 0.4e5, CH_SSC_A: 0.20e5}    # lower thresholds for FSC/SSC
BOUND_MAX = {CH_FSC_A: 2.0e5, CH_SSC_A: 1.3e5}     # upper thresholds for FSC/SSC

# Singlet gate (FSC-H ~ FSC-A)
SINGLET_RATIO_LOW, SINGLET_RATIO_HIGH = 0.85, 1.15 # acceptable FSC-H/FSC-A ratio

# Plot sizing (inches), like set_panel_size in R
PLOT_W = 1.5 * 1.618                               # plot width (golden-ish)
PLOT_H = 1.5                                       # plot height

# Aesthetics
POINT_COLOR = "#4DBBD5FF"                          # point color (NPG-like)
THEME = theme_classic()                            # clean ggplot-like theme

# ----------------------------
# Gating helpers
# ----------------------------
def apply_boundary_gate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply coarse boundary gate on FSC-A and SSC-A.

    Parameters
    ----------
    df : pd.DataFrame
        Events table with FSC-A and SSC-A.

    Returns
    -------
    pd.DataFrame
        Subset passing the boundary gate.
    """
    m = (                                            # build mask within bounds
        (df[CH_FSC_A] >= BOUND_MIN[CH_FSC_A]) & (df[CH_FSC_A] <= BOUND_MAX[CH_FSC_A]) &
        (df[CH_SSC_A] >= BOUND_MIN[CH_SSC_A]) & (df[CH_SSC_A] <= BOUND_MAX[CH_SSC_A])
    )
    out = df.loc[m]                                  # subset rows by mask
    print(f"[gate] boundary: {len(out):,}/{len(df):,} kept")  # debug summary
    return out                                       # return gated events


def apply_singlet_gate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep singlets using FSC-H / FSC-A ratio.

    Parameters
    ----------
    df : pd.DataFrame
        Events table with FSC-H and FSC-A.

    Returns
    -------
    pd.DataFrame
        Subset passing the singlet gate.
    """
    ratio = df[CH_FSC_H] / (df[CH_FSC_A].replace(0, np.nan))  # compute ratio, avoid /0
    m = (ratio >= SINGLET_RATIO_LOW) & (ratio <= SINGLET_RATIO_HIGH)  # mask within window
    out = df.loc[m]                                            # subset
    print(f"[gate] singlet : {len(out):,}/{len(df):,} kept")   # debug summary
    return out                                                 # return singlets


def median_channels_for_file(fpath: str) -> pd.Series:
    """
    Compute per-file channel medians after gating.

    Parameters
    ----------
    fpath : str
        Path to .fcs file.

    Returns
    -------
    pd.Series
        Median values (numeric-only) with '__filename' field.
    """
    print(f"[read] {os.path.basename(fpath)}")                 # which file
    sample = FCMeasurement(ID=os.path.basename(fpath),         # load FCS events
                           datafile=fpath).data
    print(f"[read] events={len(sample):,}, cols={len(sample.columns)}")  # size preview
    for ch in (CH_FSC_A, CH_SSC_A, CH_FSC_H, CH_BFP, CH_mCh):  # validate channels
        if ch not in sample.columns:
            raise ValueError(f"Channel '{ch}' not found in: {fpath}")
    gated = apply_boundary_gate(sample)                        # boundary gate
    if len(gated) == 0:                                        # fallback if empty
        print("[gate] boundary empty → using raw sample")
        gated = sample
    singlets = apply_singlet_gate(gated)                       # singlet gate
    if len(singlets) == 0:                                     # fallback if empty
        print("[gate] singlet empty  → using boundary-gated")
        singlets = gated
    s = singlets.median(numeric_only=True)                     # per-channel medians
    print(f"[median] BFP~{s.get(CH_BFP, np.nan):.3g}, mCh~{s.get(CH_mCh, np.nan):.3g}")  # preview
    s["__filename"] = os.path.splitext(os.path.basename(fpath))[0]  # store basename
    return s                                                   # return medians series


def load_flowset_medians(folder: str) -> pd.DataFrame:
    """
    Load all FCS files in a folder and compute gated medians.

    Parameters
    ----------
    folder : str
        Directory with .fcs files.

    Returns
    -------
    pd.DataFrame
        One row per file with medians and '__filename'.
    """
    files = sorted(glob.glob(os.path.join(folder, "*.fcs")))   # discover files
    print(f"[scan] {folder} → {len(files)} files")             # debug
    if not files:                                              # guard: no inputs
        raise FileNotFoundError(f"No FCS files found in: {folder}")
    out = pd.DataFrame([median_channels_for_file(f) for f in files])  # compute all medians
    print(f"[table] medians shape={out.shape}")                # debug
    return out                                                 # return table

# ----------------------------
# NFC background (mean of first 3 medians)
# ----------------------------
def compute_nfc_background(nfc_dir: str):
    """
    Estimate background medians from NFC files (first 3 medians).

    Parameters
    ----------
    nfc_dir : str
        Directory with NFC .fcs files.

    Returns
    -------
    (float, float)
        (mBFP_neg, mmCherry_neg) background medians.
    """
    df = load_flowset_medians(nfc_dir)                          # load NFC medians
    print(f"[NFC] files={len(df)}; using first 3 for background")  # debug
    mBFP_neg = float(df.iloc[:3][CH_BFP].mean())                # mean BFP background
    mmCherry_neg = float(df.iloc[:3][CH_mCh].mean())            # mean mCherry background
    print(f"[NFC] mBFP_neg={mBFP_neg:.6g}, mmCherry_neg={mmCherry_neg:.6g}")  # debug
    return mBFP_neg, mmCherry_neg                               # return tuple

# ----------------------------
# Parse filename: tokens
# ----------------------------
FILENAME_SPLIT_RE = re.compile(r"_")                            # split on underscores

def parse_timecourse_name(name: str):
    """
    Parse plasmid/exp/rep/time tokens from filename stem.

    Parameters
    ----------
    name : str
        Basename without extension.

    Returns
    -------
    (str, str, str, float)
        (plasmid, exp, rep, time)
    """
    parts = FILENAME_SPLIT_RE.split(name)                       # token list
    plasmid = parts[2] if len(parts) > 2 else ""                # plasmid token
    exp     = parts[3] if len(parts) > 3 else ""                # experiment token
    rep     = parts[4] if len(parts) > 4 else ""                # replicate token
    time_s  = parts[5] if len(parts) > 5 else ""                # time token (string)
    try:
        time = float(time_s)                                    # direct float cast
    except Exception:
        m = re.search(r"(\d+(?:\.\d+)?)", time_s)               # fallback: first number
        time = float(m.group(1)) if m else np.nan               # NaN if not found
    return plasmid, exp, rep, time                              # return tuple

# ----------------------------
# Build Rev dataset
# ----------------------------
def load_rev_timecourse(mBFP_neg: float, mmCherry_neg: float) -> pd.DataFrame:
    """
    Build Rev-only dataset with background-subtracted medians.

    Parameters
    ----------
    mBFP_neg : float
        BFP background (NFC mean of first 3 medians).
    mmCherry_neg : float
        mCherry background (NFC mean of first 3 medians).

    Returns
    -------
    pd.DataFrame
        Columns: [BV421-A, PE-A, plasmid, exp, rep, time]
    """
    med = load_flowset_medians("fcs_files/time-course_data").copy()  # load medians
    print(f"[REV] medians (pre-subtraction) shape={med.shape}")      # debug
    med[CH_BFP] = med[CH_BFP] - mBFP_neg                             # subtract BFP bg
    med[CH_mCh] = med[CH_mCh] - mmCherry_neg                         # subtract mCh bg
    print(f"[REV] bg-sub preview:\n{med[[CH_BFP, CH_mCh]].head().to_string(index=False)}")  # preview
    parsed = med["__filename"].apply(parse_timecourse_name).tolist() # parse tokens
    parsed_df = pd.DataFrame(parsed, columns=["plasmid", "exp", "rep", "time"])  # to frame
    df = pd.concat([med[[CH_BFP, CH_mCh]], parsed_df], axis=1)       # combine data+meta
    print(f"[REV] combined shape={df.shape}")                        # debug
    df = df[df["exp"] == "Rev"].copy()                               # keep Rev experiments
    print(f"[REV] after exp=='Rev': n={len(df)}")                    # debug
    df["time"] = pd.to_numeric(df["time"], errors="coerce")          # numeric time
    before = len(df)                                                 # count before dropna
    df = df.dropna(subset=["time"])                                  # drop invalid times
    print(f"[REV] dropped {before - len(df)} rows with non-numeric time")  # debug
    return df                                                        # return Rev dataset

# ----------------------------
# Normalization: min-max scaled BFP per plasmid
# ----------------------------
def add_minmax_norm(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add per-plasmid min–max normalization for downregulation:

      mean.final = mean(BFP | time > 10)
      mean.init  = mean(BFP | time == 0)
      norm.bfp   = (BFP - mean.final) / (mean.init - mean.final)

    Parameters
    ----------
    df : pd.DataFrame
        Rev dataset with columns [plasmid, time, BV421-A].

    Returns
    -------
    pd.DataFrame
        Input with added columns mean.final, mean.init, norm.bfp.
    """
    mean_final = (                                              # per-plasmid final mean
        df.query("time > 10").groupby("plasmid")[CH_BFP]
          .mean().rename("mean.final").reset_index()
    )
    mean_init = (                                               # per-plasmid init mean
        df.query("time == 0").groupby("plasmid")[CH_BFP]
          .mean().rename("mean.init").reset_index()
    )
    print(f"[norm] mean.final rows={len(mean_final)}, mean.init rows={len(mean_init)}")  # debug
    out = df.merge(mean_final, on="plasmid", how="left")         # attach final mean
    out = out.merge(mean_init, on="plasmid", how="left")         # attach init mean
    denom = (out["mean.init"] - out["mean.final"]).replace(0, np.nan)  # protect /0
    out["norm.bfp"] = (out[CH_BFP] - out["mean.final"]) / denom # compute normalized BFP
    print("[norm] preview:\n" +                                  # print small preview
          out[["plasmid", "time", CH_BFP, "mean.init", "mean.final", "norm.bfp"]]
             .head().to_string(index=False))
    return out                                                   # return normalized table

# ----------------------------
# Model & fitting: y = exp(-t * ln(2) / t_half)
# ----------------------------
def exp_decay(t, t_half):
    """
    Exponential decay model:
      y(t) = exp(-t * ln(2) / t_half)

    Parameters
    ----------
    t : array-like
        Time points.
    t_half : float
        Half-time parameter.

    Returns
    -------
    np.ndarray
        Model values y(t).
    """
    return np.exp(-t * (math.log(2.0) / t_half))                 # vectorized decay


def fit_half_time(t, y, start=0.1):
    """
    Fit exponential decay to (t, y) to estimate half-time.

    Parameters
    ----------
    t : array-like
        Time points.
    y : array-like
        Normalized response (ideally in [0, 1]).
    start : float, optional
        Initial guess for t_half, by default 0.1.

    Returns
    -------
    (float, float)
        (t_half estimate, standard error). SE is NaN if covariance missing.

    Raises
    ------
    RuntimeError
        If fewer than 3 finite points are available.
    """
    t = np.asarray(t, float)                                     # ensure float array
    y = np.asarray(y, float)                                     # ensure float array
    m = np.isfinite(t) & np.isfinite(y)                          # finite mask
    t, y = t[m], y[m]                                            # filtered pairs
    print(f"[fit] points used: {len(t)}")                        # debug
    if len(t) < 3:                                               # guard against underfit
        raise RuntimeError("Not enough points for decay fit.")
    popt, pcov = curve_fit(exp_decay, t, y, p0=[float(start)], maxfev=20000)  # NLS fit
    se = float(np.sqrt(np.diag(pcov))[0]) if pcov is not None else np.nan     # SE of t_half
    print(f"[fit] t_half={popt[0]:.6g}, SE={se:.3g}")            # debug
    return float(popt[0]), se                                    # return estimate & SE

# ----------------------------
# Plot helper (match R look)
# ----------------------------
def save_rev_plot(df: pd.DataFrame, t_half: float, out_pdf: str):
    """
    Save a Rev plot (norm.bfp vs time) with fitted decay curve.

    Parameters
    ----------
    df : pd.DataFrame
        Data for a single plasmid with 'time' and 'norm.bfp'.
    t_half : float
        Estimated half-time to draw the curve.
    out_pdf : str
        Output PDF filename (saved in OUT_PATH).
    """
    tmin, tmax = float(df["time"].min()), float(df["time"].max())  # plotting range
    curve = pd.DataFrame({"time": np.linspace(tmin, tmax, 200)})   # dense time grid
    curve["fit"] = exp_decay(curve["time"], t_half)                # model predictions
    p = (                                                          # build plot
        ggplot(df, aes("time", "norm.bfp"))
        + geom_point(size=0.4, alpha=0.7, color=POINT_COLOR)
        + geom_line(data=curve, mapping=aes(x="time", y="fit"), color="black")
        + coord_cartesian(ylim=(-0.15, 1.4))
        + scale_y_continuous(breaks=[0, 0.25, 0.5, 0.75, 1.0])
        + labs(x="Time (hours)", y="tagBFP (% of final)")
        + THEME
    )
    out_path = os.path.join(OUT_PATH, out_pdf)                    # full output path
    p.save(out_path, width=PLOT_W, height=PLOT_H, units="in")    # write PDF
    print(f"[plot] saved: {out_path}")                            # debug

# ----------------------------
# Main
# ----------------------------
def main():
    # 1) NFC background
    mBFP_neg, mmCherry_neg = compute_nfc_background("fcs_files/NFC")   # get backgrounds
    # 2) Load Rev time-course with background-subtracted medians
    rev = load_rev_timecourse(mBFP_neg, mmCherry_neg)                  # build dataset
    print(f"[REV] dataset rows={len(rev)}, cols={rev.shape[1]}")       # debug shape
    # 3) Min–max normalize per plasmid
    rev = add_minmax_norm(rev)                                         # add norm.bfp
    # 4) Fit per required plasmid & save plots
    targets = [                                                         # (plasmid, pdf, label)
        ("SP430ABA", "REV_SP430ABA_fitting.pdf",  "SP430A"),
        ("SP430",    "REV_dCas9_fitting.pdf",     "SP430"),
        ("SP428",    "REV_KRAB-dCas9_fitting.pdf","SP428"),
        ("SP427",    "REV_HDAC4-dCas9_fitting.pdf","SP427"),
        ("SP411",    "REV_CasRx_fitting.pdf",     "SP411"),
    ]
    rows = []                                                           # collect results
    for plasmid, pdfname, label in targets:                             # iterate targets
        dfp = rev[rev["plasmid"] == plasmid].copy()                     # subset per plasmid
        print(f"[target] {plasmid}: n={len(dfp)}")                      # debug n
        if dfp.empty:                                                   # skip if no data
            print(f"[WARN] No Rev data for {plasmid}; skipping.")
            continue
        t_half, se = fit_half_time(dfp["time"], dfp["norm.bfp"], start=0.1)  # fit decay
        print(f"{label}: t1/2 = {t_half:.6g}  (SE={se:.3g})")           # concise result
        rows.append({"plasmid": label, "halftime": t_half, "se": se})   # add to table
        save_rev_plot(dfp, t_half, pdfname)                             # save plot
    # 5) Save half-times
    if rows:                                                            # if we have results
        ht = pd.DataFrame(rows, columns=["plasmid", "halftime", "se"])  # build DataFrame
        out_csv = os.path.join(PARAM_PATH, "half_times_downregulation.csv")  # output path
        ht.to_csv(out_csv, index=False)                                 # write csv
        print("\nSaved:", out_csv)                                      # echo path
        print(ht.to_string(index=False))                                # show table
    else:                                                               # nothing fitted
        print("[WARN] No half-times estimated.")


if __name__ == "__main__":
    main()

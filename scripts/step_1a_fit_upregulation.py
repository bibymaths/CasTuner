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
# This Python implementation reproduces Step 1a:
#   • Loads .fcs files, applies the same gating logic and normalization
#   • Fits an exponential rise-to-one model using scipy.optimize.curve_fit
#   • Generates matching plots and parameter tables for validation
#
# All thresholds, formulas, and transformations are preserved from the original
# R code to maintain reproducibility and comparability.
# -----------------------------------------------------------------------------

"""
KD (up-regulation) time-course analysis.

This module loads flow cytometry FCS files, applies boundary and singlet gates,
computes per-file medians for relevant channels, constructs a KD-only time-course
dataset, performs per-plasmid min–max normalization on BFP, fits an exponential
rise-to-one model to estimate half-times (t1/2), generates diagnostic plots, and
writes the estimated parameters to CSV.

Outputs:
  - plots/KD_KRAB-Split-dCas9_fitting.pdf
  - plots/KD_dCas9_fitting.pdf
  - plots/KD_KRAB-dCas9_fitting.pdf
  - plots/KD_HDAC4-dCas9_fitting.pdf
  - plots/KD_CasRx_fitting.pdf
  - parameters/half_times_upregulation.csv
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
OUT_PATH = "plots"  # directory to save plots
PARAM_PATH = "parameters"  # directory to save parameter tables
os.makedirs(OUT_PATH, exist_ok=True)  # ensure plots directory exists
os.makedirs(PARAM_PATH, exist_ok=True)  # ensure parameters directory exists

# Channel names used in FCS files
CH_FSC_A = "FSC-A"      # forward scatter area
CH_SSC_A = "SSC-A"      # side scatter area
CH_FSC_H = "FSC-H"      # forward scatter height
CH_BFP   = "BV421-A"    # tagBFP channel
CH_mCh   = "PE-A"       # mCherry channel

# Boundary gate thresholds for FSC-A and SSC-A
BOUND_MIN = {CH_FSC_A: 0.4e5, CH_SSC_A: 0.20e5}  # lower bounds
BOUND_MAX = {CH_FSC_A: 2.0e5, CH_SSC_A: 1.3e5}   # upper bounds

# Singlet gate based on FSC-H/FSC-A ratio
SINGLET_RATIO_LOW, SINGLET_RATIO_HIGH = 0.85, 1.15  # acceptable ratio range

# Plot sizing (inches)
PLOT_W = 1.5 * 1.618  # golden-ish width
PLOT_H = 1.5          # height

# Aesthetics
POINT_COLOR = "#4DBBD5FF"  # point color for scatter
THEME = theme_classic()    # clean theme

# ----------------------------
# Gating helpers
# ----------------------------
def apply_boundary_gate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep events within coarse FSC/SSC bounds.

    Parameters
    ----------
    df : pd.DataFrame
        Raw events table with at least FSC-A and SSC-A columns.

    Returns
    -------
    pd.DataFrame
        Subset of events passing the boundary gate.
    """
    # Build boolean mask for events within both FSC-A and SSC-A ranges
    m = (
        (df[CH_FSC_A] >= BOUND_MIN[CH_FSC_A]) & (df[CH_FSC_A] <= BOUND_MAX[CH_FSC_A]) &
        (df[CH_SSC_A] >= BOUND_MIN[CH_SSC_A]) & (df[CH_SSC_A] <= BOUND_MAX[CH_SSC_A])
    )
    out = df.loc[m]  # apply mask
    print(f"[gate] boundary: kept {len(out):,}/{len(df):,} events")  # debug
    return out  # gated events


def apply_singlet_gate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep singlets using FSC-H/FSC-A ratio.

    Parameters
    ----------
    df : pd.DataFrame
        Events table with FSC-H and FSC-A.

    Returns
    -------
    pd.DataFrame
        Subset of events passing the singlet gate.
    """
    # Avoid division by zero by turning zeros in FSC-A into NaN temporarily
    ratio = df[CH_FSC_H] / (df[CH_FSC_A].replace(0, np.nan))  # compute ratio
    # Keep events whose ratio lies in the pre-defined window
    m = (ratio >= SINGLET_RATIO_LOW) & (ratio <= SINGLET_RATIO_HIGH)
    out = df.loc[m]  # apply mask
    print(f"[gate] singlet:  kept {len(out):,}/{len(df):,} events")  # debug
    return out  # singlets


def median_channels_for_file(fpath: str, mBFP_neg: float = 0.0, mmCherry_neg: float = 0.0) -> pd.Series:
    """
    Compute per-channel medians for one FCS file, with gating and optional background subtraction.

    Parameters
    ----------
    fpath : str
        Path to the .fcs file.
    mBFP_neg : float, optional
        Background (negative control) median for BFP to subtract before median, by default 0.0.
    mmCherry_neg : float, optional
        Background (negative control) median for mCherry to subtract before median, by default 0.0.

    Returns
    -------
    pd.Series
        Median values (numeric-only) per channel with an extra '__filename' field (basename without extension).

    Raises
    ------
    ValueError
        If any required channel is missing in the FCS data.
    """
    print(f"[read] {os.path.basename(fpath)}")  # debug: which file is being read
    # Load events table from FCS
    sample = FCMeasurement(ID=os.path.basename(fpath), datafile=fpath).data  # raw events
    print(f"[read] events: {len(sample):,}, columns: {list(sample.columns)[:6]}...")  # preview

    # Validate required channels exist
    for ch in (CH_FSC_A, CH_SSC_A, CH_FSC_H, CH_BFP, CH_mCh):
        if ch not in sample.columns:
            raise ValueError(f"Channel '{ch}' not found in: {fpath}")

    # Apply boundary gate; if empty, fallback to original
    gated = apply_boundary_gate(sample)
    if len(gated) == 0:
        print("[gate] boundary produced 0 events; using raw sample")  # debug
        gated = sample

    # Apply singlet gate; if empty, fallback to boundary-gated
    singlets = apply_singlet_gate(gated)
    if len(singlets) == 0:
        print("[gate] singlet produced 0 events; using boundary-gated")  # debug
        singlets = gated

    # Background subtraction (if provided) BEFORE computing median
    if mBFP_neg != 0.0 or mmCherry_neg != 0.0:
        singlets_sub = singlets.copy()  # copy to avoid mutating upstream data
        singlets_sub[CH_BFP] = singlets_sub[CH_BFP] - mBFP_neg  # subtract BFP background
        singlets_sub[CH_mCh] = singlets_sub[CH_mCh] - mmCherry_neg  # subtract mCherry background
        s = singlets_sub.median(numeric_only=True)  # compute medians
        print(f"[median] (bg-sub) BFP~{s.get(CH_BFP, np.nan):.3g}, mCh~{s.get(CH_mCh, np.nan):.3g}")  # debug
    else:
        s = singlets.median(numeric_only=True)  # medians without subtraction
        print(f"[median] (raw)    BFP~{s.get(CH_BFP, np.nan):.3g}, mCh~{s.get(CH_mCh, np.nan):.3g}")  # debug

    # Store filename (sans extension) for later parsing
    s["__filename"] = os.path.splitext(os.path.basename(fpath))[0]
    return s  # series of medians for this file


def load_flowset_medians(folder: str, mBFP_neg: float = 0.0, mmCherry_neg: float = 0.0) -> pd.DataFrame:
    """
    Load all .fcs files in a folder and compute gated channel medians per file.

    Parameters
    ----------
    folder : str
        Directory containing .fcs files.
    mBFP_neg : float, optional
        Background median to subtract from BFP before per-file median, by default 0.0.
    mmCherry_neg : float, optional
        Background median to subtract from mCherry before per-file median, by default 0.0.

    Returns
    -------
    pd.DataFrame
        One row per file with medians and '__filename'.

    Raises
    ------
    FileNotFoundError
        If no .fcs files are found in `folder`.
    """
    files = sorted(glob.glob(os.path.join(folder, "*.fcs")))  # discover files
    print(f"[scan] folder={folder} files={len(files)}")  # debug: number of files
    if not files:
        raise FileNotFoundError(f"No FCS files found in: {folder}")
    # Compute medians for each file and combine into DataFrame
    out = pd.DataFrame([median_channels_for_file(f, mBFP_neg, mmCherry_neg) for f in files])
    print(f"[table] medians shape={out.shape}")  # debug
    return out  # table of per-file medians

# ----------------------------
# NFC background (mean of first 3 medians)
# ----------------------------
def compute_nfc_background(nfc_dir: str):
    """
    Estimate background medians from NFC (negative control) files.

    Uses the mean of the first 3 NFC medians for BFP and mCherry.

    Parameters
    ----------
    nfc_dir : str
        Directory with NFC .fcs files.

    Returns
    -------
    (float, float)
        Tuple of (mBFP_neg, mmCherry_neg).
    """
    df = load_flowset_medians(nfc_dir)  # load medians for NFC files
    print(f"[NFC] using first 3 files for background, total NFC files={len(df)}")  # debug
    mBFP_neg = float(df.iloc[:3][CH_BFP].mean())  # mean BFP of first 3
    mmCherry_neg = float(df.iloc[:3][CH_mCh].mean())  # mean mCherry of first 3
    print(f"[NFC] mBFP_neg={mBFP_neg:.6g}, mmCherry_neg={mmCherry_neg:.6g}")  # debug
    return mBFP_neg, mmCherry_neg  # background levels

# ----------------------------
# Parse filename: tokens
# ----------------------------
FILENAME_SPLIT_RE = re.compile(r"_")  # split on underscore


def parse_timecourse_name(name: str):
    """
    Parse tokens from a filename (already stripped of extension) to extract:
    plasmid, experiment label, replicate, and numeric time.

    Expected pattern (R analogy):
      separate(rowname, c(NA,NA,"plasmid","exp","rep","time",NA,NA), "_")

    Parameters
    ----------
    name : str
        Basename of file without extension.

    Returns
    -------
    (str, str, str, float)
        plasmid, exp, rep, time (np.nan if not parsed).
    """
    parts = FILENAME_SPLIT_RE.split(name)  # split by underscores
    plasmid = parts[2] if len(parts) > 2 else ""  # third token
    exp     = parts[3] if len(parts) > 3 else ""  # fourth token
    rep     = parts[4] if len(parts) > 4 else ""  # fifth token
    time_s  = parts[5] if len(parts) > 5 else ""  # sixth token (string)
    try:
        time = float(time_s)  # direct float cast if clean
    except Exception:
        # Fallback: extract first numeric group like "12" or "12.5"
        m = re.search(r"(\d+(?:\.\d+)?)", time_s)
        time = float(m.group(1)) if m else np.nan
    return plasmid, exp, rep, time  # parsed tokens

# ----------------------------
# Build KD dataset
# ----------------------------
def load_kd_timecourse(mBFP_neg: float, mmCherry_neg: float) -> pd.DataFrame:
    """
    Build the KD-only time-course dataset with background-subtracted medians.

    Parameters
    ----------
    mBFP_neg : float
        NFC-derived BFP background to subtract.
    mmCherry_neg : float
        NFC-derived mCherry background to subtract.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns [BV421-A, PE-A, plasmid, exp, rep, time].
    """
    # Load medians with background subtraction from the time-course directory
    med = load_flowset_medians(
        "fcs_files/time-course_data",
        mBFP_neg,
        mmCherry_neg
    ).copy()
    print(f"[KD] raw medians shape={med.shape}")  # debug

    # Parse filename tokens into a structured DataFrame
    parsed = med["__filename"].apply(parse_timecourse_name).tolist()
    parsed_df = pd.DataFrame(parsed, columns=["plasmid", "exp", "rep", "time"])
    print(f"[KD] parsed tokens shape={parsed_df.shape}")  # debug

    # Combine measured channels with parsed metadata
    df = pd.concat([med[[CH_BFP, CH_mCh]], parsed_df], axis=1)
    print(f"[KD] combined table shape={df.shape}")  # debug
    print(df.head().to_string(index=False))  # preview

    # Keep only KD experiments
    df = df[df["exp"] == "KD"].copy()
    print(f"[KD] after exp=='KD': {len(df)} rows")  # debug

    # Ensure time is numeric and drop rows with missing/invalid times
    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    before = len(df)
    df = df.dropna(subset=["time"])
    print(f"[KD] dropped {before - len(df)} rows with non-numeric time")  # debug

    # Set ordered categorical for plasmid to match R factor levels
    cat = pd.CategoricalDtype(["SP430", "SP428", "SP430ABA", "SP427", "SP411"], ordered=True)
    df["plasmid"] = df["plasmid"].astype(cat)
    print(f"[KD] plasmids present: {df['plasmid'].dropna().unique().tolist()}")  # debug
    return df  # KD-only dataset

# ----------------------------
# Normalization (KD)
# ----------------------------
def add_minmax_norm_kd(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add per-plasmid min–max normalization of BFP:
      mean.final = mean(BFP | time > 10)
      mean.init  = mean(BFP | time == 0)
      norm.bfp   = (BFP - mean.init) / (mean.final - mean.init)

    Parameters
    ----------
    df : pd.DataFrame
        KD dataset with columns [plasmid, time, BV421-A].

    Returns
    -------
    pd.DataFrame
        Input df with extra columns: mean.final, mean.init, norm.bfp.
    """
    # Compute per-plasmid final mean (time > 10)
    mean_final = (
        df.query("time > 10")
          .groupby("plasmid")[CH_BFP]
          .mean().rename("mean.final").reset_index()
    )
    print(f"[norm] mean.final rows: {len(mean_final)}")  # debug

    # Compute per-plasmid initial mean (time == 0)
    mean_init = (
        df.query("time == 0")
          .groupby("plasmid")[CH_BFP]
          .mean().rename("mean.init").reset_index()
    )
    print(f"[norm] mean.init rows:  {len(mean_init)}")  # debug

    # Merge means back to each row
    out = df.merge(mean_final, on="plasmid", how="left").merge(mean_init, on="plasmid", how="left")
    # Denominator for min–max scaling; protect against division by zero
    denom = (out["mean.final"] - out["mean.init"]).replace(0, np.nan)
    # Compute normalized BFP
    out["norm.bfp"] = (out[CH_BFP] - out["mean.init"]) / denom

    # Quick sanity preview
    print("[norm] preview:")
    print(out[["plasmid", "time", CH_BFP, "mean.init", "mean.final", "norm.bfp"]]
          .head().to_string(index=False))
    return out  # normalized table

# ----------------------------
# Model & fitting (KD)
# ----------------------------
def exp_rise_to_one(t, t_half):
    """
    Exponential rise-to-one model:
      y(t) = 1 - exp(-t * ln(2) / t_half)

    Parameters
    ----------
    t : array-like
        Time points (hours).
    t_half : float
        Half-time parameter.

    Returns
    -------
    np.ndarray
        Model values y(t).
    """
    return 1.0 - np.exp(-t * (math.log(2.0) / t_half))  # vectorized formula


def fit_half_time_rise(t, y, start=0.8):
    """
    Fit the exponential rise-to-one model to (t, y) and estimate half-time.

    Parameters
    ----------
    t : array-like
        Time points.
    y : array-like
        Normalized response in [0, 1] (ideally).
    start : float, optional
        Initial guess for t_half, by default 0.8.

    Returns
    -------
    (float, float)
        Estimated t_half and its standard error (SE). SE is NaN if covariance unavailable.

    Raises
    ------
    RuntimeError
        If there are fewer than 3 finite points to fit.
    """
    t = np.asarray(t, float)  # ensure float array
    y = np.asarray(y, float)  # ensure float array
    m = np.isfinite(t) & np.isfinite(y)  # finite mask
    t, y = t[m], y[m]  # keep finite pairs
    print(f"[fit] points used: {len(t)}")  # debug
    if len(t) < 3:
        raise RuntimeError("Not enough points for growth fit.")
    # Nonlinear least squares fit
    popt, pcov = curve_fit(exp_rise_to_one, t, y, p0=[float(start)], maxfev=20000)
    # Standard error of t_half from covariance matrix
    se = float(np.sqrt(np.diag(pcov))[0]) if pcov is not None else np.nan
    print(f"[fit] t_half={popt[0]:.6g}, SE={se:.3g}")  # debug
    return float(popt[0]), se  # parameter and SE

# ----------------------------
# Plot helper (match R look)
# ----------------------------
def save_kd_plot(df: pd.DataFrame, t_half: float, out_pdf: str):
    """
    Save a KD plot with points (norm.bfp vs time) and fitted curve to a PDF.

    Parameters
    ----------
    df : pd.DataFrame
        Subset for one plasmid with columns time, norm.bfp.
    t_half : float
        Fitted half-time to draw the curve.
    out_pdf : str
        Output PDF filename (basename only; saved under OUT_PATH).
    """
    # Define prediction grid for a smooth curve
    tmin, tmax = float(df["time"].min()), float(df["time"].max())
    curve = pd.DataFrame({"time": np.linspace(tmin, tmax, 200)})
    curve["fit"] = exp_rise_to_one(curve["time"], t_half)  # model predictions

    # Build plot
    p = (
        ggplot(df, aes("time", "norm.bfp"))
        + geom_point(size=0.4, alpha=0.7, color=POINT_COLOR)
        + geom_line(data=curve, mapping=aes(x="time", y="fit"), color="black")
        + coord_cartesian(ylim=(-0.15, 1.2))
        + scale_y_continuous(breaks=[0, 0.25, 0.5, 0.75, 1.0])
        + labs(x="Time (hours)", y="tagBFP (% of final)")
        + THEME
    )
    # Save to disk
    out_path = os.path.join(OUT_PATH, out_pdf)
    p.save(out_path, width=PLOT_W, height=PLOT_H, units="in")
    print(f"[plot] saved: {out_path}")  # debug

# ----------------------------
# Main
# ----------------------------
def main():
    """
    Orchestrate KD up-regulation analysis end-to-end:
      1) Compute NFC background medians (BFP, mCherry)
      2) Load KD time-course medians with background subtraction
      3) Add per-plasmid min–max normalization on BFP
      4) Fit exponential rise model per target plasmid and save plots
      5) Save half-times to CSV
    """
    # 1) NFC background
    print("\n== Step 1: NFC background ==")  # banner
    mBFP_neg, mmCherry_neg = compute_nfc_background("fcs_files/NFC")

    # 2) Load KD time-course with background-subtracted medians
    print("\n== Step 2: Load KD time-course ==")  # banner
    kd = load_kd_timecourse(mBFP_neg, mmCherry_neg)
    print(f"[KD] rows={len(kd)}, cols={kd.shape[1]}")  # debug

    # 3) Min–max normalize per plasmid (KD definition)
    print("\n== Step 3: Normalize BFP (min–max per plasmid) ==")  # banner
    kd = add_minmax_norm_kd(kd)

    # 4) Fit per plasmid & save plots
    print("\n== Step 4: Fit half-times and plot ==")  # banner
    targets = [
        ("SP430ABA", "KD_KRAB-Split-dCas9_fitting.pdf", "SP430A"),
        ("SP430",    "KD_dCas9_fitting.pdf",            "SP430"),
        ("SP428",    "KD_KRAB-dCas9_fitting.pdf",       "SP428"),
        ("SP427",    "KD_HDAC4-dCas9_fitting.pdf",      "SP427"),
        ("SP411",    "KD_CasRx_fitting.pdf",            "SP411"),
    ]

    rows = []  # collect results for CSV
    for plasmid, pdfname, label in targets:
        dfp = kd[kd["plasmid"] == plasmid].copy()  # subset per plasmid
        print(f"[target] {plasmid}: n={len(dfp)}")  # debug
        if dfp.empty:
            print(f"[WARN] No KD data for {plasmid}; skipping.")  # warn and continue
            continue

        # Fit exponential rise model; start guess 0.8 hours (adjust if needed)
        t_half, se = fit_half_time_rise(dfp["time"], dfp["norm.bfp"], start=0.8)
        print(f"{label}: t1/2 = {t_half:.6g}  (SE={se:.3g})")  # concise result line

        # Record row for CSV output
        rows.append({"plasmid": label, "halftime": t_half, "se": se})

        # Save diagnostic plot for this plasmid
        save_kd_plot(dfp, t_half, pdfname)

    # 5) Save half-times
    print("\n== Step 5: Save parameters ==")  # banner
    if rows:
        ht = pd.DataFrame(rows, columns=["plasmid", "halftime", "se"])  # build result table
        out_csv = os.path.join(PARAM_PATH, "half_times_upregulation.csv")  # destination path
        ht.to_csv(out_csv, index=False)  # write to disk
        print("Saved:", out_csv)  # confirm
        print(ht.to_string(index=False))  # echo table
    else:
        print("[WARN] No half-times estimated.")  # nothing to write

# Entrypoint
if __name__ == "__main__":
    main()  # run the pipeline

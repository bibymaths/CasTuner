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
#   this work showcases independent reproducibility, modeling rigor, and the
#   ability to extend CasTuner’s analytical stages in Python.
#
# Original R workflow (for reference):
#   Step 1a  Fit tagBFP upregulation dynamics after dTAG-13 withdrawal
#   Step 1b  Fit tagBFP downregulation dynamics after dTAG-13 addition
#   Step 1c  Fit steady-state dose–response (Hill curves)  ← (this file)
#   Step 2   ODE simulations for mCherry derepression ("ODE_REV.R")
#   Step 3   ODE simulations for mCherry repression ("ODE_KD.R")
#
# This Python implementation reproduces Step 1c:
#   • Loads .fcs files, applies boundary + singlet gating and NFC background subtraction
#   • Constructs day-4 steady-state dose–response tables (dTAG-13 vs reporter)
#   • Normalizes repressor (BFP) and reporter (mCherry) signals
#   • Fits Hill functions (K, n) via scipy.optimize.curve_fit and exports parameters
#   • Generates diagnostic plots for each plasmid/guide target
#
# All thresholds, formulas, and transformations are preserved from the original
# R code to maintain reproducibility and comparability.
# -----------------------------------------------------------------------------

"""
Fit Hill functions to dTAG-13 dose–response curves of the reporter–repressor.

This module:
  1) Loads FCS files for dose–response replicates (R1–R3) and NFC.
  2) Applies boundary and singlet gates; subtracts NFC-derived backgrounds.
  3) Maps plate indices to dTAG-13 concentrations; constructs day-4 steady state.
  4) Computes fold-change relative to NTC and normalizes BFP (repressor) levels.
  5) Fits Hill curves (fc vs normalized BFP) to estimate K (EC50-like) and Hill n.
  6) Saves per-plasmid parameters (CSV) and diagnostic plots (PDFs).

Outputs:
  plots/Hill_dCas9.pdf
  plots/Hill_CasRx.pdf
  plots/Hill-HDAC4-dCas9.pdf
  plots/Hill-KRAB-dCas9.pdf
  plots/Hill-KRAB-Split-dCas9.pdf
  parameters/Hill_parameters.csv

"""


import os, re, glob, warnings
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
# Config
# ----------------------------
OUT_PATH = "plots"                               # directory for plot outputs
PARAM_PATH = "parameters"                        # directory for parameter CSVs
os.makedirs(OUT_PATH, exist_ok=True)             # ensure plot dir exists
os.makedirs(PARAM_PATH, exist_ok=True)           # ensure parameter dir exists

CH_FSC_A = "FSC-A"                               # forward scatter area
CH_SSC_A = "SSC-A"                               # side scatter area
CH_FSC_H = "FSC-H"                               # forward scatter height
CH_BFP   = "BV421-A"                             # tagBFP channel
CH_mCh   = "PE-A"                                # mCherry channel

BOUND_MIN = {CH_FSC_A: 0.4e5, CH_SSC_A: 0.20e5}  # lower gating bounds
BOUND_MAX = {CH_FSC_A: 2.0e5, CH_SSC_A: 1.3e5}   # upper gating bounds
SINGLET_RATIO_LOW, SINGLET_RATIO_HIGH = 0.85, 1.15  # singlet ratio window

# plate-index → dTAG concentration (nM)
DTAG_MAP = {
    "1": 0, "2": 0.5, "3": 1, "4": 2, "5": 3, "6": 5,
    "7": 8, "8": 10, "9": 25, "10": 50, "11": 100, "12": 500
}

PLOT_W = 1.5 * 1.618                              # plot width (inches)
PLOT_H = 1.5                                      # plot height (inches)

def theme_castuner_like():
    """Return a clean plot theme resembling the original R plots."""
    return theme_classic()

# ----------------------------
# Gating helpers
# ----------------------------
def apply_boundary_gate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply coarse FSC/SSC boundary gate.

    Parameters
    ----------
    df : pd.DataFrame
        Event table with FSC-A and SSC-A.

    Returns
    -------
    pd.DataFrame
        Subset of events within the rectangle gate.
    """
    m = (  # build mask within numeric limits
        (df[CH_FSC_A] >= BOUND_MIN[CH_FSC_A]) & (df[CH_FSC_A] <= BOUND_MAX[CH_FSC_A]) &
        (df[CH_SSC_A] >= BOUND_MIN[CH_SSC_A]) & (df[CH_SSC_A] <= BOUND_MAX[CH_SSC_A])
    )
    out = df.loc[m]                                # filter events
    print(f"[gate] boundary: kept {len(out):,}/{len(df):,}")  # debug
    return out                                     # return gated events

def apply_singlet_gate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep singlets using FSC-H/FSC-A ratio window.

    Parameters
    ----------
    df : pd.DataFrame
        Event table with FSC-H and FSC-A.

    Returns
    -------
    pd.DataFrame
        Events passing the singlet gate.
    """
    ratio = df[CH_FSC_H] / (df[CH_FSC_A].replace(0, np.nan))   # robust ratio
    m = (ratio >= SINGLET_RATIO_LOW) & (ratio <= SINGLET_RATIO_HIGH)  # in-window
    out = df.loc[m]                                            # filter
    print(f"[gate] singlet : kept {len(out):,}/{len(df):,}")   # debug
    return out                                                 # singlets only

def median_channels_for_file(fpath: str) -> pd.Series:
    """
    Compute per-file medians after boundary + singlet gating.

    Parameters
    ----------
    fpath : str
        Path to an .fcs file.

    Returns
    -------
    pd.Series
        Numeric medians per channel with '__filename' set to basename (sans extension).

    Raises
    ------
    ValueError
        If required channels are missing.
    """
    print(f"[read] {os.path.basename(fpath)}")                 # which file
    sample = FCMeasurement(ID=os.path.basename(fpath), datafile=fpath).data  # read FCS
    print(f"[read] events={len(sample):,}, cols={len(sample.columns)}")      # size
    for ch in (CH_FSC_A, CH_SSC_A, CH_FSC_H, CH_BFP, CH_mCh):  # channel check
        if ch not in sample.columns:
            raise ValueError(f"Channel '{ch}' not found in: {fpath}")
    gated = apply_boundary_gate(sample)                        # boundary gate
    if len(gated) == 0:                                        # fallback
        print("[gate] boundary empty → using raw sample")
        gated = sample
    singlets = apply_singlet_gate(gated)                       # singlet gate
    if len(singlets) == 0:                                     # fallback
        print("[gate] singlet empty  → using boundary-gated")
        singlets = gated
    s = singlets.median(numeric_only=True)                     # medians
    print(f"[median] BFP~{s.get(CH_BFP, np.nan):.3g} | mCh~{s.get(CH_mCh, np.nan):.3g}")  # preview
    s["__filename"] = os.path.splitext(os.path.basename(fpath))[0]  # filename stem
    return s                                                   # return series

def load_flowset_medians(folder: str) -> pd.DataFrame:
    """
    Load all .fcs files under a folder and compute gated medians.

    Parameters
    ----------
    folder : str
        Directory path containing .fcs files.

    Returns
    -------
    pd.DataFrame
        One row per file with medians and '__filename'.

    Raises
    ------
    FileNotFoundError
        If no .fcs files were found.
    """
    files = sorted(glob.glob(os.path.join(folder, "*.fcs")))   # discover FCS files
    print(f"[scan] {folder} → {len(files)} files")             # debug
    if not files:                                              # guard
        raise FileNotFoundError(f"No FCS files found in: {folder}")
    out = pd.DataFrame([median_channels_for_file(f) for f in files])  # compute medians
    print(f"[table] medians shape={out.shape}")                # debug
    return out                                                 # return table

# ----------------------------
# Background subtraction (NFC)
# ----------------------------
def compute_nfc_background(nfc_dir: str):
    """
    Compute NFC background medians (mean of first 3 files).

    Parameters
    ----------
    nfc_dir : str
        Directory of NFC .fcs files.

    Returns
    -------
    (float, float)
        Tuple (mBFP_neg, mmCherry_neg).
    """
    df = load_flowset_medians(nfc_dir)                         # NFC medians
    print(f"[NFC] files={len(df)}; using first 3 for background")  # debug
    mBFP_neg = float(df.iloc[:3][CH_BFP].mean())               # BFP background
    mmCherry_neg = float(df.iloc[:3][CH_mCh].mean())           # mCherry background
    print(f"[NFC] mBFP_neg={mBFP_neg:.6g}, mmCherry_neg={mmCherry_neg:.6g}")  # debug
    return mBFP_neg, mmCherry_neg                              # return tuple

# ----------------------------
# Filename parsing
# ----------------------------
FILENAME_SPLIT_RE = re.compile(r"_")                           # underscore splitter

def parse_filename_tokens(name: str):
    """
    Extract plasmid, guide, and plate-index token (for dTAG concentration).

    Parameters
    ----------
    name : str
        Basename of file without extension.

    Returns
    -------
    (str, str, str)
        (plasmid, guide, dtag_idx_token)
    """
    parts = FILENAME_SPLIT_RE.split(name)                      # tokens
    plasmid = parts[0] if len(parts) > 0 else ""               # plasmid code
    guide   = parts[1] if len(parts) > 1 else ""               # guide label 'G'/'N'
    dtag_token = parts[2] if len(parts) > 2 else ""            # token like 'dTAG11'
    m = re.search(r"([A-Za-z]+)(\d+)$", dtag_token)            # split alpha+digits
    dtag_idx = m.group(2) if m else dtag_token                 # keep index digits
    return plasmid, guide, dtag_idx                            # return tokens

# ----------------------------
# Load replicate
# ----------------------------
def load_replicate(folder: str, mBFP_neg: float, mmCherry_neg: float) -> pd.DataFrame:
    """
    Load one replicate folder: medians → background subtraction → token parsing.

    Parameters
    ----------
    folder : str
        Path to replicate folder with .fcs files.
    mBFP_neg : float
        NFC background for BFP.
    mmCherry_neg : float
        NFC background for mCherry.

    Returns
    -------
    pd.DataFrame
        Columns: [BV421-A, PE-A, plasmid, guide, dTAG]
    """
    raw = load_flowset_medians(folder).copy()                  # per-file medians
    print(f"[rep] {folder} shape={raw.shape}")                 # debug
    raw[CH_BFP] = raw[CH_BFP] - mBFP_neg                       # subtract BFP background
    raw[CH_mCh] = raw[CH_mCh] - mmCherry_neg                   # subtract mCh background
    print(f"[rep] bg-sub head:\n{raw[[CH_BFP, CH_mCh]].head().to_string(index=False)}")  # preview
    parsed = raw["__filename"].apply(parse_filename_tokens).tolist()  # parse tokens
    parsed_df = pd.DataFrame(parsed, columns=["plasmid", "guide", "dTAG_token"])         # to frame
    df = pd.concat([raw[[CH_BFP, CH_mCh]], parsed_df], axis=1) # combine
    df["dTAG"] = df["dTAG_token"].map(DTAG_MAP).astype(float)  # map to concentration
    df = df.drop(columns=["dTAG_token"])                       # drop helper
    print(f"[rep] combined shape={df.shape}")                   # debug
    return df                                                  # return replicate table

# ----------------------------
# Hill model
# ----------------------------
def hill_func(R, K, n):
    """
    Hill repression function:
      y = K^n / (K^n + R^n)

    Parameters
    ----------
    R : array-like
        Repressor proxy (normalized BFP).
    K : float
        Half-maximal constant.
    n : float
        Hill coefficient.

    Returns
    -------
    np.ndarray
        Predicted normalized reporter level.
    """
    return (K**n) / (K**n + np.power(R, n))                    # vectorized equation

def fit_hill(R_vals, y_vals, start=(0.1, 1.0)):
    """
    Fit Hill function to data using nonlinear least squares.

    Parameters
    ----------
    R_vals : array-like
        Normalized repressor levels (predictor).
    y_vals : array-like
        Fold-change (response).
    start : (float, float), optional
        Initial guesses (K, n), by default (0.1, 1.0).

    Returns
    -------
    (np.ndarray, np.ndarray | None)
        (popt, pcov) where popt=[K, n]; pcov is covariance matrix or None.

    Raises
    ------
    RuntimeError
        If fewer than 3 finite data points are available.
    """
    R_vals = np.asarray(R_vals, dtype=float)                   # to float array
    y_vals = np.asarray(y_vals, dtype=float)                   # to float array
    m = np.isfinite(R_vals) & np.isfinite(y_vals)              # finite mask
    R_vals, y_vals = R_vals[m], y_vals[m]                      # filter
    print(f"[fit] hill points used: {len(R_vals)}")            # debug
    if len(R_vals) < 3:                                        # guard
        raise RuntimeError("Not enough points for Hill fit.")
    popt, pcov = curve_fit(hill_func, R_vals, y_vals,          # LM fit (like nlsLM)
                           p0=np.array(start, float), maxfev=20000)
    print(f"[fit] K={popt[0]:.6g}, n={popt[1]:.6g}")           # debug
    return popt, pcov                                          # params + covariance

def save_hill_plot(df: pd.DataFrame, fit_params, out_pdf: str):
    """
    Save Hill fit plot (fc vs normalized BFP) to PDF.

    Parameters
    ----------
    df : pd.DataFrame
        Data with columns 'norm.bfp' and 'fc'.
    fit_params : (float, float)
        Tuple (K, n) of fitted parameters.
    out_pdf : str
        Output PDF filename (saved under OUT_PATH).
    """
    K, n = fit_params                                          # unpack params
    tmin, tmax = float(df["norm.bfp"].min()), float(df["norm.bfp"].max())  # x-range
    curve = pd.DataFrame({"norm.bfp": np.linspace(tmin, tmax, 200)})       # grid
    curve["fc_fit"] = hill_func(curve["norm.bfp"], K, n)       # predictions
    p = (                                                      # build plot
        ggplot(df, aes("norm.bfp", "fc"))
        + geom_point(size=0.8, alpha=0.4)
        + geom_line(data=curve, mapping=aes(x="norm.bfp", y="fc_fit"))
        + coord_cartesian(ylim=(0, 1))
        + labs(x="Normalized repressor level", y="Normalized reporter level")
        + scale_y_continuous()
        + theme_castuner_like()
    )
    out_path = os.path.join(OUT_PATH, out_pdf)                 # destination
    p.save(out_path, width=PLOT_W, height=PLOT_H, units="in")  # write PDF
    print(f"[plot] saved: {out_path}")                         # debug

# ----------------------------
# Main
# ----------------------------
def main():
    # NFC background
    mBFP_neg, mmCherry_neg = compute_nfc_background("fcs_files/NFC")

    # Replicates with background subtraction & parsing
    rep1 = load_replicate("fcs_files/dose_response_data/R1", mBFP_neg, mmCherry_neg).assign(rep=1, day=4)
    rep2 = load_replicate("fcs_files/dose_response_data/R2", mBFP_neg, mmCherry_neg).assign(rep=2, day=4)
    rep3 = load_replicate("fcs_files/dose_response_data/R3", mBFP_neg, mmCherry_neg).assign(rep=3, day=4)
    d4 = pd.concat([rep1, rep2, rep3], ignore_index=True)
    print(f"[data] pooled day4 shape={d4.shape}")

    # 1) mean over NTC rows, grouped only by plasmid
    meanNTC_pl = (
        d4.query("guide == 'N'")
        .groupby(["plasmid"])[CH_mCh]
        .mean()
        .rename("meanNTC")
        .reset_index()
    )
    print(f"[norm] meanNTC_pl rows={len(meanNTC_pl)}")

    # 2) expand to the exact (plasmid, dTAG) pairs present in NTC
    ntc_pairs = (
        d4.query("guide == 'N'")[["plasmid", "dTAG"]]
        .drop_duplicates()
        .merge(meanNTC_pl, on="plasmid", how="left")
    )
    print(f"[norm] ntc_pairs rows={len(ntc_pairs)}")

    # 3) join back by (plasmid, dTAG) exactly like the R left_join
    d4 = d4.merge(ntc_pairs, on=["plasmid", "dTAG"], how="left")
    print(f"[norm] after join shape={d4.shape}")

    # 4) compute fold-change
    d4["fc"] = d4[CH_mCh] / d4["meanNTC"]
    print("[norm] fc preview:\n" + d4[["plasmid", "guide", "dTAG", "fc"]].head().to_string(index=False))

    # max.bfp at dTAG==0 per (plasmid, guide)
    max_bfp = (
        d4.query("dTAG == 0")
        .groupby(["plasmid", "guide"])[CH_BFP]
        .mean().rename("max.bfp").reset_index()
    )
    d4 = d4.merge(max_bfp, on=["plasmid", "guide"], how="left")
    d4["fc.bfp"] = d4[CH_BFP] / d4["max.bfp"]
    print("[norm] max_bfp merge preview:\n" +
          d4[["plasmid", "guide", "dTAG", CH_BFP, "max.bfp", "fc.bfp"]].head().to_string(index=False))

    # norm.bfp: min–max per (plasmid, guide)
    d4["norm.bfp"] = (
        d4.groupby(["plasmid", "guide"], group_keys=False)[CH_BFP]
        .apply(lambda v: (v - v.min()) / (v.max() - v.min()) if v.max() > v.min() else 0.0)
    )
    print("[norm] norm.bfp preview:\n" +
          d4[["plasmid", "guide", "dTAG", CH_BFP, "norm.bfp"]].head().to_string(index=False))

    # keep targeting guides (guide=='G')
    d4g = d4.query("guide == 'G'").copy()
    print(f"[data] guide=='G' rows={len(d4g)}")

    # Fit & plot per plasmid
    records = []
    targets = [
        ("430",    "Hill_dCas9.pdf",            "SP430"),
        ("411",    "Hill_CasRx.pdf",            "SP411"),
        ("427",    "Hill-HDAC4-dCas9.pdf",      "SP427"),
        ("428",    "Hill-KRAB-dCas9.pdf",       "SP428"),
        ("430ABA", "Hill-KRAB-Split-dCas9.pdf", "SP430A"),
    ]

    for plasmid, pdfname, label in targets:
        subset = d4g[d4g["plasmid"] == plasmid].copy()
        print(f"[target] {plasmid}: n={len(subset)}")
        if subset.empty:
            print(f"[WARN] No data for plasmid {plasmid}, skipping.")
            continue

        R_vals = subset["norm.bfp"]
        y_vals = subset["fc"]

        try:
            (K, n), _ = fit_hill(R_vals, y_vals, start=(0.1, 1.0))
            print(f"{label}: K={K:.6g}  n={n:.6g}")
            records.append({"plasmid": label, "K": float(K), "n": float(n)})
            save_hill_plot(subset, (K, n), pdfname)
        except Exception as e:
            print(f"[ERROR] Fit failed for {label}: {e}")

    # Save parameters
    if records:
        params_df = pd.DataFrame(records, columns=["plasmid", "K", "n"])
        out_csv = os.path.join(PARAM_PATH, "Hill_parameters.csv")
        params_df.to_csv(out_csv, index=False)
        print("\nSaved:", out_csv)
        print(params_df.to_string(index=False))
    else:
        print("[WARN] No parameters estimated; no CSV written.")

if __name__ == "__main__":
    main()

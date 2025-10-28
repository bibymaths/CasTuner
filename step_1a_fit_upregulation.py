#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fit Hill functions to dTAG-13 dose–response curves of the reporter–repressor.

Python dependencies:
  - FlowCytometryTools
  - flowio, flowutils
  - numpy, pandas, scipy
  - plotnine

Directory layout expected (same as R script):
  fcs_files/
    NFC/
      *.fcs
    dose_response_data/
      R1/*.fcs
      R2/*.fcs
      R3/*.fcs

Outputs:
  plots/Hill_dCas9.pdf
  plots/Hill_CasRx.pdf
  plots/Hill-HDAC4-dCas9.pdf
  plots/Hill-KRAB-dCas9.pdf
  plots/Hill-KRAB-Split-dCas9.pdf
  parameters/Hill_parameters.csv
"""

import os
import re
import glob
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from FlowCytometryTools import FCMeasurement
from plotnine import (
    ggplot, aes, geom_point, geom_line, labs,
    scale_y_continuous, coord_cartesian, theme_classic
)

# ----------------------------
# Config & small utilities
# ----------------------------
OUT_PATH = "plots"
PARAM_PATH = "parameters"
os.makedirs(OUT_PATH, exist_ok=True)
os.makedirs(PARAM_PATH, exist_ok=True)

# Channel names (as in the R code)
CH_FSC_A = "FSC-A"
CH_SSC_A = "SSC-A"
CH_FSC_H = "FSC-H"
CH_BFP   = "BV421-A"  # reporter (BFP)
CH_mCh   = "PE-A"     # reporter (mCherry)

# Boundary gate ranges (same numeric limits as R script)
BOUND_MIN = {CH_FSC_A: 0.4e5, CH_SSC_A: 0.20e5}
BOUND_MAX = {CH_FSC_A: 2.0e5, CH_SSC_A: 1.3e5}

# Singlet gate mimic (FSC-H ~ FSC-A)
SINGLET_RATIO_LOW  = 0.85
SINGLET_RATIO_HIGH = 1.15

# dTAG mapping from filename token to numeric concentration
DTAG_MAP = {
    "1": 0, "2": 0.5, "3": 1, "4": 2, "5": 3, "6": 5,
    "7": 8, "8": 10, "9": 25, "10": 50, "11": 100, "12": 500
}

# Plot sizes (inches) ~ R set_panel_size: width = 1.5*1.618, height = 1.5
PLOT_W = 1.5 * 1.618
PLOT_H = 1.5

# Theme mimic
def theme_castuner_like():
    return theme_classic()

# ----------------------------
# Gating helpers
# ----------------------------
def apply_boundary_gate(df: pd.DataFrame) -> pd.DataFrame:
    """Rectangular gate on FSC-A / SSC-A."""
    m = (
        (df[CH_FSC_A] >= BOUND_MIN[CH_FSC_A]) & (df[CH_FSC_A] <= BOUND_MAX[CH_FSC_A]) &
        (df[CH_SSC_A] >= BOUND_MIN[CH_SSC_A]) & (df[CH_SSC_A] <= BOUND_MAX[CH_SSC_A])
    )
    return df.loc[m]

def apply_singlet_gate(df: pd.DataFrame) -> pd.DataFrame:
    """Singlet gate approximation via FSC-H ≈ FSC-A."""
    ratio = df[CH_FSC_H] / (df[CH_FSC_A].replace(0, np.nan))
    m = (ratio >= SINGLET_RATIO_LOW) & (ratio <= SINGLET_RATIO_HIGH)
    return df.loc[m]

def median_channels_for_file(fpath: str) -> pd.Series:
    """Load one FCS file and return medians of all channels after gating."""
    sample = FCMeasurement(ID=os.path.basename(fpath), datafile=fpath).data
    # Drop missing columns gracefully
    for ch in (CH_FSC_A, CH_SSC_A, CH_FSC_H, CH_BFP, CH_mCh):
        if ch not in sample.columns:
            raise ValueError(f"Channel '{ch}' not found in: {fpath}")

    gated = apply_boundary_gate(sample)
    if len(gated) == 0:
        gated = sample  # fallback if gate is too strict

    singlets = apply_singlet_gate(gated)
    if len(singlets) == 0:
        singlets = gated  # fallback

    return singlets.median(numeric_only=True)

def load_flowset_medians(folder: str) -> pd.DataFrame:
    """Load all FCS files in a folder and compute medians after gating."""
    files = sorted(glob.glob(os.path.join(folder, "*.fcs")))
    if not files:
        raise FileNotFoundError(f"No FCS files found in: {folder}")
    rows = []
    for f in files:
        med = median_channels_for_file(f)
        med["__filename"] = os.path.splitext(os.path.basename(f))[0]
        rows.append(med)
    return pd.DataFrame(rows)

# ----------------------------
# Background subtraction (NFC)
# ----------------------------
def compute_nfc_background(nfc_dir: str):
    df = load_flowset_medians(nfc_dir)
    # R code used mean of first 3 rows for both channels
    mBFP_neg = df.iloc[:3][CH_BFP].mean()
    mmCherry_neg = df.iloc[:3][CH_mCh].mean()
    return float(mBFP_neg), float(mmCherry_neg)

# ----------------------------
# Dose-response replicate processing
# ----------------------------
FILENAME_SPLIT_RE = re.compile(r"_")

def parse_filename_tokens(name: str):
    """
    R split logic:
      rowname -> separate into plasmid, guide, dTAG, NA
      then split dTAG token into letters + number (we keep the numeric part).
    Expected pattern: <plasmid>_<guide>_<dTAG><index>_...
    Example: 430_G_dTAG3_...
    """
    parts = FILENAME_SPLIT_RE.split(name)
    # defensive defaults
    plasmid = parts[0] if len(parts) > 0 else ""
    guide = parts[1] if len(parts) > 1 else ""
    dtag_token = parts[2] if len(parts) > 2 else ""

    # extract trailing digits from the dTAG token
    m = re.search(r"([A-Za-z]+)(\d+)$", dtag_token)
    dtag_idx = m.group(2) if m else dtag_token
    return plasmid, guide, dtag_idx

def load_replicate(folder: str, mBFP_neg: float, mmCherry_neg: float) -> pd.DataFrame:
    """
    Load replicate folder, subtract background, compute medians,
    and parse filename into columns: plasmid, guide, dTAG (numeric).
    """
    raw = load_flowset_medians(folder)

    # Background subtraction (linear transform)
    raw[CH_BFP] = raw[CH_BFP] - mBFP_neg
    raw[CH_mCh] = raw[CH_mCh] - mmCherry_neg

    # Keep only BFP and mCherry columns (positions 7:8 in R, but safer by names)
    df = raw[[CH_BFP, CH_mCh, "__filename"]].copy()

    # Parse filename pieces
    parsed = df["__filename"].apply(parse_filename_tokens).tolist()
    parsed_df = pd.DataFrame(parsed, columns=["plasmid", "guide", "dTAG_token"])
    df = pd.concat([df, parsed_df], axis=1)

    # Map dTAG token to numeric concentration
    df["dTAG"] = df["dTAG_token"].map(DTAG_MAP).astype(float)
    df = df.drop(columns=["__filename", "dTAG_token"])
    return df.rename(columns={CH_BFP: f"{CH_BFP}", CH_mCh: f"{CH_mCh}"})

# ----------------------------
# Hill model fitting
# ----------------------------
def hill_func(R, K, n):
    # y = K^n / (K^n + R^n)  (R: normalized repressor level)
    return (K**n) / (K**n + np.power(R, n))

def fit_hill(R_vals, y_vals, start=(0.1, 1.0)):
    # Ensure arrays
    R_vals = np.asarray(R_vals, dtype=float)
    y_vals = np.asarray(y_vals, dtype=float)
    # clean NaNs/Infs
    m = np.isfinite(R_vals) & np.isfinite(y_vals)
    R_vals = R_vals[m]
    y_vals = y_vals[m]
    if len(R_vals) < 3:
        raise RuntimeError("Not enough points for Hill fit.")

    # positive bounds for K and n
    popt, pcov = curve_fit(hill_func, R_vals, y_vals,
                           p0=[0.1, 1.0], maxfev=20000)
    return popt, pcov

def save_hill_plot(df, fit_params, out_pdf):
    K, n = fit_params
    # Build smooth curve across observed R
    tmin, tmax = float(df["norm.bfp"].min()), float(df["norm.bfp"].max())
    curve = pd.DataFrame({"norm.bfp": np.linspace(tmin, tmax, 200)})
    curve["fc_fit"] = hill_func(curve["norm.bfp"], K, n)

    p = (
        ggplot(df, aes("norm.bfp", "fc"))
        + geom_point(size=0.8, alpha=0.4)
        + geom_line(data=curve, mapping=aes(x="norm.bfp", y="fc_fit"))
        + coord_cartesian(ylim=(0, 1))
        + labs(x="Normalized repressor level", y="Normalized reporter level")
        + scale_y_continuous()
        + theme_castuner_like()
    )
    p.save(os.path.join(OUT_PATH, out_pdf), width=PLOT_W, height=PLOT_H, units="in")

# ----------------------------
# Main flow (mirrors R)
# ----------------------------
def main():
    # 1) NFC background
    mBFP_neg, mmCherry_neg = compute_nfc_background("fcs_files/NFC")

    # 2) Replicates R1..R3 with background subtraction & parsing
    rep1 = load_replicate("fcs_files/dose_response_data/R1", mBFP_neg, mmCherry_neg)
    rep2 = load_replicate("fcs_files/dose_response_data/R2", mBFP_neg, mmCherry_neg)
    rep3 = load_replicate("fcs_files/dose_response_data/R3", mBFP_neg, mmCherry_neg)

    rep1 = rep1.assign(rep=1, day=4)
    rep2 = rep2.assign(rep=2, day=4)
    rep3 = rep3.assign(rep=3, day=4)

    d4 = pd.concat([rep1, rep2, rep3], ignore_index=True)

    # 3) Factor/ordering (optional in Python)
    # levels = ["430", "428", "430ABA", "427", "411"]

    # meanNTC from NTC (guide == 'N')
    meanNTC = (
        d4.query("guide == 'N'")
        .groupby(['plasmid', 'dTAG'])[CH_mCh]
        .mean()
        .rename('meanNTC')
        .reset_index()
    )

    # join by (plasmid, dTAG) and compute fc
    d4 = d4.merge(meanNTC, on=['plasmid', 'dTAG'], how='left')
    d4['fc'] = d4[CH_mCh] / d4['meanNTC']

    # 6) fc.bfp and norm.bfp
    # max.bfp at dTAG==0 per (plasmid, guide)
    max_bfp = (d4.query("dTAG == 0")
               .groupby(['plasmid', 'guide'])[CH_BFP]
               .mean().rename('max.bfp').reset_index())
    d4 = d4.merge(max_bfp, on=['plasmid', 'guide'], how='left')
    d4['fc.bfp'] = d4[CH_BFP] / d4['max.bfp']

    # norm.bfp: min–max per (plasmid, guide)
    d4['norm.bfp'] = (
        d4.groupby(['plasmid', 'guide'], group_keys=False)[CH_BFP]
        .apply(lambda v: (v - v.min()) / (v.max() - v.min()) if v.max() > v.min() else 0.0)
    )

    # 7) Filter targeting guides (guide == 'G')
    d4g = d4.query("guide == 'G'").copy()

    # 8) Fit per plasmid and plot
    records = []
    targets = [
        ("430",    "Hill_dCas9.pdf",             "SP430"),
        ("411",    "Hill_CasRx.pdf",             "SP411"),
        ("427",    "Hill-HDAC4-dCas9.pdf",       "SP427"),
        ("428",    "Hill-KRAB-dCas9.pdf",        "SP428"),
        ("430ABA", "Hill-KRAB-Split-dCas9.pdf",  "SP430A"),
    ]

    for plasmid, pdfname, label in targets:
        subset = d4g[d4g["plasmid"] == plasmid].copy()
        if subset.empty:
            print(f"[WARN] No data for plasmid {plasmid}, skipping.")
            continue

        R_vals = subset["norm.bfp"]
        y_vals = subset["fc"]

        try:
            (K, n), pcov = fit_hill(R_vals, y_vals, start=(0.1, 1.0))
            print(f"Summary {label}: K={K:.4g}, n={n:.4g}")
            records.append({"plasmid": label, "K": K, "n": n})

            save_hill_plot(subset, (K, n), pdfname)
        except Exception as e:
            print(f"[ERROR] Fit failed for {label}: {e}")

    # 9) Save parameters
    if records:
        params_df = pd.DataFrame.from_records(records, columns=["plasmid", "K", "n"])
        params_df.to_csv(os.path.join(PARAM_PATH, "Hill_parameters.csv"), index=False)
        print(params_df)
    else:
        print("[WARN] No parameters estimated; no CSV written.")

if __name__ == "__main__":
    main()

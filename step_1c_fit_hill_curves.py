#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fit Hill functions to dTAG-13 dose–response curves of the reporter–repressor.

Requires (versions compatible with Py3.9):
  FlowCytometryTools, flowio, flowutils
  numpy<2.0, pandas, scipy
  plotnine
"""

import os
import re
import glob
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
# Config
# ----------------------------
OUT_PATH = "plots"
PARAM_PATH = "parameters"
os.makedirs(OUT_PATH, exist_ok=True)
os.makedirs(PARAM_PATH, exist_ok=True)

CH_FSC_A = "FSC-A"
CH_SSC_A = "SSC-A"
CH_FSC_H = "FSC-H"
CH_BFP   = "BV421-A"  # tagBFP
CH_mCh   = "PE-A"     # mCherry

BOUND_MIN = {CH_FSC_A: 0.4e5, CH_SSC_A: 0.20e5}
BOUND_MAX = {CH_FSC_A: 2.0e5, CH_SSC_A: 1.3e5}
SINGLET_RATIO_LOW, SINGLET_RATIO_HIGH = 0.85, 1.15

DTAG_MAP = {
    "1": 0, "2": 0.5, "3": 1, "4": 2, "5": 3, "6": 5,
    "7": 8, "8": 10, "9": 25, "10": 50, "11": 100, "12": 500
}

PLOT_W = 1.5 * 1.618
PLOT_H = 1.5

def theme_castuner_like():
    return theme_classic()

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
# Background subtraction (NFC)
# ----------------------------
def compute_nfc_background(nfc_dir: str):
    df = load_flowset_medians(nfc_dir)
    # match R: mean of first 3 rows
    mBFP_neg = float(df.iloc[:3][CH_BFP].mean())
    mmCherry_neg = float(df.iloc[:3][CH_mCh].mean())
    return mBFP_neg, mmCherry_neg

# ----------------------------
# Filename parsing
# ----------------------------
FILENAME_SPLIT_RE = re.compile(r"_")

def parse_filename_tokens(name: str):
    # Expect: <plasmid>_<guide>_<dTAG<index>>_...
    parts = FILENAME_SPLIT_RE.split(name)
    plasmid = parts[0] if len(parts) > 0 else ""
    guide   = parts[1] if len(parts) > 1 else ""
    dtag_token = parts[2] if len(parts) > 2 else ""
    m = re.search(r"([A-Za-z]+)(\d+)$", dtag_token)
    dtag_idx = m.group(2) if m else dtag_token
    return plasmid, guide, dtag_idx

# ----------------------------
# Load replicate
# ----------------------------
def load_replicate(folder: str, mBFP_neg: float, mmCherry_neg: float) -> pd.DataFrame:
    raw = load_flowset_medians(folder).copy()

    # background subtraction (linear)
    raw[CH_BFP] = raw[CH_BFP] - mBFP_neg
    raw[CH_mCh] = raw[CH_mCh] - mmCherry_neg

    parsed = raw["__filename"].apply(parse_filename_tokens).tolist()
    parsed_df = pd.DataFrame(parsed, columns=["plasmid", "guide", "dTAG_token"])

    df = pd.concat([raw[[CH_BFP, CH_mCh]], parsed_df], axis=1)
    df["dTAG"] = df["dTAG_token"].map(DTAG_MAP).astype(float)
    df = df.drop(columns=["dTAG_token"])
    return df

# ----------------------------
# Hill model
# ----------------------------
def hill_func(R, K, n):
    # y = K^n / (K^n + R^n)
    return (K**n) / (K**n + np.power(R, n))

def fit_hill(R_vals, y_vals, start=(0.1, 1.0)):
    R_vals = np.asarray(R_vals, dtype=float)
    y_vals = np.asarray(y_vals, dtype=float)
    m = np.isfinite(R_vals) & np.isfinite(y_vals)
    R_vals, y_vals = R_vals[m], y_vals[m]
    if len(R_vals) < 3:
        raise RuntimeError("Not enough points for Hill fit.")
    # Levenberg–Marquardt (unbounded), like nlsLM
    popt, pcov = curve_fit(hill_func, R_vals, y_vals, p0=np.array(start, float), maxfev=20000)
    return popt, pcov

def save_hill_plot(df, fit_params, out_pdf):
    K, n = fit_params
    tmin, tmax = float(df["norm.bfp"].min()), float(df["norm.bfp"].max())
    curve = pd.DataFrame({"norm.bfp": np.linspace(tmin, tmax, 200)})
    curve["fc_fit"] = hill_func(curve["norm.bfp"], K, n)

    p = (
        ggplot(df, aes("norm.bfp", "fc"))
        + geom_point(size=0.8, alpha=0.4)   # gray in R; leave default here
        + geom_line(data=curve, mapping=aes(x="norm.bfp", y="fc_fit"))
        + coord_cartesian(ylim=(0, 1))
        + labs(x="Normalized repressor level", y="Normalized reporter level")
        + scale_y_continuous()
        + theme_castuner_like()
    )
    p.save(os.path.join(OUT_PATH, out_pdf), width=PLOT_W, height=PLOT_H, units="in")

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

    # 1) mean over NTC rows, grouped only by plasmid
    meanNTC_pl = (
        d4.query("guide == 'N'")
        .groupby(["plasmid"])[CH_mCh]
        .mean()
        .rename("meanNTC")
        .reset_index()
    )

    # 2) expand to the exact (plasmid, dTAG) pairs present in NTC
    ntc_pairs = (
        d4.query("guide == 'N'")[["plasmid", "dTAG"]]
        .drop_duplicates()
        .merge(meanNTC_pl, on="plasmid", how="left")
    )

    # 3) join back by (plasmid, dTAG) exactly like the R left_join
    d4 = d4.merge(ntc_pairs, on=["plasmid", "dTAG"], how="left")

    # 4) compute fold-change
    d4["fc"] = d4[CH_mCh] / d4["meanNTC"]

    # max.bfp at dTAG==0 per (plasmid, guide)
    max_bfp = (
        d4.query("dTAG == 0")
        .groupby(["plasmid", "guide"])[CH_BFP]
        .mean().rename("max.bfp").reset_index()
    )
    d4 = d4.merge(max_bfp, on=["plasmid", "guide"], how="left")
    d4["fc.bfp"] = d4[CH_BFP] / d4["max.bfp"]

    # norm.bfp: min–max per (plasmid, guide)
    d4["norm.bfp"] = (
        d4.groupby(["plasmid", "guide"], group_keys=False)[CH_BFP]
        .apply(lambda v: (v - v.min()) / (v.max() - v.min()) if v.max() > v.min() else 0.0)
    )

    # keep targeting guides (guide=='G')
    d4g = d4.query("guide == 'G'").copy()

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
        params_df.to_csv(os.path.join(PARAM_PATH, "Hill_parameters.csv"), index=False)
        print("\nSaved:", os.path.join(PARAM_PATH, "Hill_parameters.csv"))
        print(params_df)
    else:
        print("[WARN] No parameters estimated; no CSV written.")

if __name__ == "__main__":
    main()

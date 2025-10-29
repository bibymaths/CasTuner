#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------------
#
# This script was independently developed by Abhinav Mishra as a Python
# reimplementation/extension of the CasTuner analysis pipeline, originally used in:
#   “CasTuner: a degron and CRISPR/Cas-based toolkit for analog tuning of
#    endogenous gene expression” — Gemma Noviello, Rutger A. F. Gjaltema,
#    and Edda G. Schulz.
#
# Purpose:
#   This code is part of a prospective research extension designed to demonstrate
#   technical capability and methodological thinking relevant to a potential PhD
#   position in the Schulz group. I am not currently affiliated with the group;
#   this work showcases independent reproducibility, modeling rigor, and the
#   ability to extend CasTuner’s ODE-based repression analysis in Python.
#
# Reference R workflow (context):
#   Step 1a  Upregulation fit (BFP) after dTAG-13 withdrawal
#   Step 1b  Downregulation fit (BFP) after dTAG-13 addition
#   Step 1c  Dose–response (Hill curves) at steady state
#   Step 2   ODE simulation for derepression (REV)
#   Step 3   ODE simulation for repression (KD)  ← this script
#
# This Python step (KD / repression) does:
#   • Loads time-course medians from .fcs, applies gates, background subtraction
#   • Builds KD-only dataset, computes mCherry fold-change & BFP normalization
#   • Loads fitted parameters (t_up, Hill K/n, alpha)
#   • Simulates ODE system, scans delays to minimize MAE vs data
#   • Exports delay per plasmid and diagnostic plots
#
# All thresholds, formulas, and transformations follow the original logic to
# maintain reproducibility and comparability.
# -----------------------------------------------------------------------------

"""
Python port of the R ‘KD (repression)’ step.

Pipeline:
1) NFC background from fcs_files/NFC → compute mBFP_neg, mmCherry_neg
   (mean of first 3 medians).
2) Load time-course medians from fcs_files/time-course_data with boundary
   + singlet gating, then subtract NFC background.
3) Parse filename tokens → (plasmid, exp, rep, time); keep KD rows (exp=="KD").
4) Normalize signals (KD):
   • fc.cherry = mCherry / mean_mCherry_at_time0  (per plasmid)
   • norm.bfp  = (BFP - mean.init) / (mean.final - mean.init)
     where mean.init = mean(BFP | t==0), mean.final = mean(BFP | t>10) per plasmid.
5) Load parameters from CSVs:
   • parameters/half_times_upregulation.csv → t_up
   • parameters/Hill_parameters.csv         → K, n
   • parameters/alphamcherry.csv            → alpha (global or per-plasmid)
6) Simulate ODEs (repression model) for each plasmid (SP430A/427/428/411/430):
   • beta = ln(2) / t_up
   • dR/dt = beta − R*(ln2/t_up)
   • dY/dt = K^n / (K^n + R^n) − alpha*Y
   • Outputs: R(t) (tagBFP proxy) and alpha*Y(t) (mCherry).
7) Delay scan: test delays 0..25 h (step 0.5); choose delay minimizing MAE
   between simulated mCherry and mean experimental fc.cherry per time.
8) Plots:
   • plots/MAE_KD_<system>_mcherry.pdf     (MAE vs delay; best delay highlighted)
   • plots/KD_ODE_mCherry_<system>.pdf     (mCherry data; base + delayed sim)
   • plots/KD_ODE_tagBFP_<system>.pdf      (normalized BFP vs simulated R)
9) Outputs:
   • parameters/delays_repression.csv with columns [plasmid, d_rev]
   • intermediate parameters are read from parameters/ as above.
"""

import os, re, glob, math, warnings
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from FlowCytometryTools import FCMeasurement
from plotnine import (
    ggplot, aes, geom_point, geom_line, labs, coord_cartesian, theme_classic
)

warnings.filterwarnings("ignore")

# ----------------------------
# Paths and constants
# ----------------------------
OUT_PATH = "plots"                                           # output directory for plots
PARAM_PATH = "parameters"                                    # output directory for parameters
os.makedirs(OUT_PATH, exist_ok=True)                         # ensure plot dir exists
os.makedirs(PARAM_PATH, exist_ok=True)                       # ensure param dir exists

CH_FSC_A = "FSC-A"                                           # forward scatter area
CH_SSC_A = "SSC-A"                                           # side scatter area
CH_FSC_H = "FSC-H"                                           # forward scatter height
CH_BFP   = "BV421-A"                                         # tagBFP channel
CH_mCh   = "PE-A"                                            # mCherry channel

BOUND_MIN = {CH_FSC_A: 0.4e5, CH_SSC_A: 0.20e5}              # lower bounds for rectangle gate
BOUND_MAX = {CH_FSC_A: 2.0e5, CH_SSC_A: 1.3e5}               # upper bounds for rectangle gate
SINGLET_RATIO_LOW, SINGLET_RATIO_HIGH = 0.85, 1.15           # singlet ratio window (FSC-H/FSC-A)

PLOT_W = 1.5 * 1.618                                         # plot width (inches)
PLOT_H = 1.5                                                 # plot height (inches)
COL_BFP = "#4DBBD5FF"                                        # BFP color
COL_CH  = "#E64B35FF"                                        # mCherry color
THEME = theme_classic()                                      # clean plot theme

# ----------------------------
# Utilities: gating & medians
# ----------------------------
def apply_boundary_gate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply coarse rectangle gate on FSC-A and SSC-A.

    Parameters
    ----------
    df : pd.DataFrame
        Raw event table.

    Returns
    -------
    pd.DataFrame
        Events within the boundary gate.
    """
    m = (  # build boolean mask for the rectangular region
        (df[CH_FSC_A] >= BOUND_MIN[CH_FSC_A]) & (df[CH_FSC_A] <= BOUND_MAX[CH_FSC_A]) &
        (df[CH_SSC_A] >= BOUND_MIN[CH_SSC_A]) & (df[CH_SSC_A] <= BOUND_MAX[CH_SSC_A])
    )
    out = df.loc[m]                                           # filter events
    print(f"[gate] boundary: kept {len(out):,}/{len(df):,}")  # debug summary
    return out                                                # return gated df

def apply_singlet_gate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep singlets based on FSC-H/FSC-A ratio.

    Parameters
    ----------
    df : pd.DataFrame
        Events (boundary-gated recommended).

    Returns
    -------
    pd.DataFrame
        Singlet events only.
    """
    ratio = df[CH_FSC_H] / (df[CH_FSC_A].replace(0, np.nan))  # robust ratio
    m = (ratio >= SINGLET_RATIO_LOW) & (ratio <= SINGLET_RATIO_HIGH)  # window
    out = df.loc[m]                                            # filter
    print(f"[gate] singlet : kept {len(out):,}/{len(df):,}")   # debug summary
    return out                                                 # return singlets

def median_channels_for_file(fpath: str) -> pd.Series:
    """
    Compute per-file medians after gating.

    Parameters
    ----------
    fpath : str
        Path to .fcs file.

    Returns
    -------
    pd.Series
        Medians per channel with '__filename' stem.
    """
    print(f"[read] {os.path.basename(fpath)}")                # which file
    data = FCMeasurement(ID=os.path.basename(fpath),          # load FCS
                         datafile=fpath).data
    print(f"[read] events={len(data):,}, cols={len(data.columns)}")  # size
    for ch in (CH_FSC_A, CH_SSC_A, CH_FSC_H, CH_BFP, CH_mCh): # channel check
        if ch not in data.columns:
            raise ValueError(f"Channel '{ch}' not in {fpath}")
    g = apply_boundary_gate(data)                             # boundary gate
    if len(g) == 0:                                           # fallback if empty
        print("[gate] boundary empty → using raw sample")
        g = data
    s = apply_singlet_gate(g)                                 # singlet gate
    if len(s) == 0:                                           # fallback if empty
        print("[gate] singlet empty  → using boundary-gated")
        s = g
    med = s.median(numeric_only=True)                         # compute medians
    print(f"[median] BFP~{med.get(CH_BFP, np.nan):.3g} | mCh~{med.get(CH_mCh, np.nan):.3g}")  # preview
    med["__filename"] = os.path.splitext(os.path.basename(fpath))[0]  # store stem
    return med                                                # return series

def load_flowset_medians(folder: str) -> pd.DataFrame:
    """
    Load all .fcs files in a folder and compute gated medians.

    Parameters
    ----------
    folder : str
        Directory containing .fcs files.

    Returns
    -------
    pd.DataFrame
        One row per file with medians and '__filename'.
    """
    files = sorted(glob.glob(os.path.join(folder, "*.fcs")))  # discover files
    print(f"[scan] {folder} → {len(files)} files")            # debug
    if not files:                                             # guard
        raise FileNotFoundError(f"No FCS in {folder}")
    out = pd.DataFrame([median_channels_for_file(f) for f in files])  # compute all
    print(f"[table] medians shape={out.shape}")               # debug
    return out                                                # return table

# ----------------------------
# NFC background (mean of first 3 medians)
# ----------------------------
def compute_nfc_background(nfc_dir: str):
    """
    Estimate NFC background medians (mean of first 3 rows).

    Parameters
    ----------
    nfc_dir : str
        Directory of NFC .fcs files.

    Returns
    -------
    (float, float)
        (mBFP_neg, mmCherry_neg) background values.
    """
    df = load_flowset_medians(nfc_dir)                        # NFC medians
    print(f"[NFC] files={len(df)}; using first 3 for background")  # debug
    mBFP_neg = float(df.iloc[:3][CH_BFP].mean())              # BFP background
    mmCherry_neg = float(df.iloc[:3][CH_mCh].mean())          # mCherry background
    print(f"[NFC] mBFP_neg={mBFP_neg:.6g}, mmCherry_neg={mmCherry_neg:.6g}")  # debug
    return mBFP_neg, mmCherry_neg                             # return tuple

# ----------------------------
# Parse filename to tokens: ... _ plasmid _ exp _ rep _ time _ ...
# ----------------------------
_SPLIT = re.compile(r"_")                                     # underscore splitter

def parse_timecourse_name(name: str):
    """
    Extract (plasmid, exp, rep, time) tokens from filename stem.

    Parameters
    ----------
    name : str
        Basename of the file (no extension).

    Returns
    -------
    (str, str, str, float)
        (plasmid, exp, rep, time)
    """
    parts = _SPLIT.split(name)                                # tokenize
    plasmid = parts[2] if len(parts) > 2 else ""              # plasmid token
    exp     = parts[3] if len(parts) > 3 else ""              # experiment token
    rep     = parts[4] if len(parts) > 4 else ""              # replicate token
    time_s  = parts[5] if len(parts) > 5 else ""              # time token (string)
    try:
        time = float(time_s)                                  # try direct parse
    except Exception:
        m = re.search(r"(\d+(?:\.\d+)?)", time_s)             # fallback: first number
        time = float(m.group(1)) if m else np.nan             # NaN if none
    return plasmid, exp, rep, time                            # return tuple

# ----------------------------
# Load time-course KD dataset with background subtraction
# ----------------------------
def load_kd(mBFP_neg: float, mmCherry_neg: float) -> pd.DataFrame:
    """
    Load time-course medians, subtract NFC background, keep KD-only rows.

    Parameters
    ----------
    mBFP_neg : float
        Background for BFP to subtract.
    mmCherry_neg : float
        Background for mCherry to subtract.

    Returns
    -------
    pd.DataFrame
        KD-only dataset with columns [BV421-A, PE-A, plasmid, exp, rep, time].
    """
    med = load_flowset_medians("fcs_files/time-course_data").copy()     # medians table
    print(f"[KD] medians pre-sub shape={med.shape}")                    # debug
    med[CH_BFP] = med[CH_BFP] - mBFP_neg                                # subtract BFP bg
    med[CH_mCh] = med[CH_mCh] - mmCherry_neg                            # subtract mCh bg
    print(f"[KD] bg-sub head:\n{med[[CH_BFP, CH_mCh]].head().to_string(index=False)}")  # preview
    parsed = med["__filename"].apply(parse_timecourse_name).tolist()    # parse tokens
    parsed_df = pd.DataFrame(parsed, columns=["plasmid", "exp", "rep", "time"])         # to frame
    df = pd.concat([med[[CH_BFP, CH_mCh]], parsed_df], axis=1)          # combine data/meta
    df = df[df["exp"] == "KD"].copy()                                   # keep KD rows
    print(f"[KD] KD-only rows: {len(df)}")                              # debug
    df["time"] = pd.to_numeric(df["time"], errors="coerce")             # numeric time
    before = len(df)                                                    # count before drop
    df = df.dropna(subset=["time"])                                     # drop invalid time
    print(f"[KD] dropped {before - len(df)} rows with non-numeric time")# debug
    cat = pd.CategoricalDtype(["SP430", "SP428", "SP430ABA", "SP427", "SP411"], ordered=True)  # order
    df["plasmid"] = df["plasmid"].astype(cat)                           # set category order
    return df                                                           # return KD table

# ----------------------------
# Compute mCherry fold-change & BFP min-max norm (KD)
# ----------------------------
def add_fc_and_norm(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add mCherry fold-change (vs t==0 mean per plasmid) and BFP min–max normalization.

      fc.cherry = mCh / mean_t0(plasmid)
      norm.bfp  = (BFP - mean.init) / (mean.final - mean.init), where
                  mean.init = mean(BFP | t==0), mean.final = mean(BFP | t>10)

    Parameters
    ----------
    df : pd.DataFrame
        KD time-course with [plasmid, time, BV421-A, PE-A].

    Returns
    -------
    pd.DataFrame
        Input with extra columns: mmcherry, fc.cherry, mean.final, mean.init, norm.bfp.
    """
    t0 = (df.query("time == 0")                                 # rows at t==0
            .groupby("plasmid")[CH_mCh]                         # by plasmid
            .mean().rename("mmcherry").reset_index())           # mean mCherry
    out = df.merge(t0, on="plasmid", how="left")                # attach mmcherry
    out["fc.cherry"] = out[CH_mCh] / out["mmcherry"]            # fold-change
    mean_final = (out.query("time > 10")                         # t>10 subset
                  .groupby("plasmid")[CH_BFP]                    # by plasmid
                  .mean().rename("mean.final").reset_index())    # final mean BFP
    mean_init  = (out.query("time == 0")                         # t==0 subset
                  .groupby("plasmid")[CH_BFP]                    # by plasmid
                  .mean().rename("mean.init").reset_index())     # initial mean BFP
    out = out.merge(mean_final, on="plasmid", how="left")        # attach mean.final
    out = out.merge(mean_init,  on="plasmid", how="left")        # attach mean.init
    denom = (out["mean.final"] - out["mean.init"]).replace(0, np.nan)  # protect /0
    out["norm.bfp"] = (out[CH_BFP] - out["mean.init"]) / denom   # min–max BFP
    print("[norm] preview:\n" +
          out[["plasmid", "time", CH_BFP, CH_mCh, "fc.cherry", "mean.init", "mean.final", "norm.bfp"]]
            .head().to_string(index=False))                      # debug preview
    return out                                                   # return enriched table

# ----------------------------
# Read parameters and merge (t_up, Hill K,n, alpha)
# ----------------------------
def load_parameters() -> pd.DataFrame:
    """
    Load and merge parameter CSVs required for KD simulation.

    Expected files under 'parameters/':
      - half_times_downregulation.csv : columns [plasmid, halftime] → t_down (not used here)
      - half_times_upregulation.csv   : columns [plasmid, halftime] → t_up
      - Hill_parameters.csv           : columns [plasmid, K, n]
      - alphamcherry.csv              : either [plasmid, alpha] or a single-row [alpha]

    Returns
    -------
    pd.DataFrame
        Merged table with columns at least: [plasmid, t_up, k, n, alpha].

    Raises
    ------
    KeyError
        If required columns/files are missing or malformed.
    """
    base = "parameters"                                           # root directory

    def _read_norm(path: str) -> pd.DataFrame:
        df = pd.read_csv(path)                                    # read CSV
        df.columns = (df.columns.str.strip().str.lower()
                                .str.replace(r"[\s\-]+", "_", regex=True))  # normalize names
        print(f"[params] loaded {os.path.basename(path)} shape={df.shape}") # debug
        return df

    # Load per-plasmid tables
    t_down = _read_norm(os.path.join(base, "half_times_downregulation.csv"))
    t_up   = _read_norm(os.path.join(base, "half_times_upregulation.csv"))
    hills  = _read_norm(os.path.join(base, "Hill_parameters.csv"))

    # Normalize expected columns
    if "halftime" in t_down.columns: t_down = t_down.rename(columns={"halftime": "t_down"})  # rename
    if "halftime" in t_up.columns:   t_up   = t_up.rename(columns={"halftime": "t_up"})      # rename
    if "k" not in hills.columns and "K" in hills.columns: hills = hills.rename(columns={"K": "k"})  # fix K

    # Validate presence of 'plasmid'
    for name, df in [("half_times_downregulation.csv", t_down),
                     ("half_times_upregulation.csv", t_up),
                     ("Hill_parameters.csv", hills)]:
        if "plasmid" not in df.columns:
            raise KeyError(f"'{name}' must have a 'plasmid' column. Found: {list(df.columns)}")

    # Merge core parameter tables
    params = (t_down.merge(t_up, on="plasmid", how="outer")
                    .merge(hills, on="plasmid", how="outer"))
    print(f"[params] merged core shape={params.shape}")           # debug

    # Load alpha; allow global or per-plasmid
    alpha_path = os.path.join(base, "alphamcherry.csv")           # path
    alpha = _read_norm(alpha_path)                                # load alpha

    if "alpha" not in alpha.columns:                              # tolerate variants
        for alt in ("alpha_mcherry", "alphamcherry", "a"):
            if alt in alpha.columns:
                alpha = alpha.rename(columns={alt: "alpha"})
                break
    if "alpha" not in alpha.columns:
        raise KeyError(f"'alphamcherry.csv' must contain 'alpha'. Found: {list(alpha.columns)}")

    if "plasmid" in alpha.columns:                                # per-plasmid alpha
        params = params.merge(alpha[["plasmid", "alpha"]], on="plasmid", how="left")
    else:                                                         # global alpha
        if len(alpha) != 1:
            raise KeyError("Global 'alphamcherry.csv' must contain exactly 1 row.")
        params["alpha"] = float(alpha["alpha"].iloc[0])           # broadcast

    # Final checks
    required = ["plasmid", "t_up", "k", "n", "alpha"]
    missing = [c for c in required if c not in params.columns]
    if missing:
        raise KeyError(f"Merged parameters missing columns: {missing}. Columns: {list(params.columns)}")

    if params[required].isna().any(axis=1).any():                 # warn on NA rows
        print("[WARN] Some parameter rows contain NA values.")
        print(params.loc[params[required].isna().any(axis=1), ["plasmid"] + required[1:]])

    print(f"[params] final shape={params.shape}")                 # debug
    return params                                                 # return merged params

# ----------------------------
# ODE system and simulation
# dR = beta - R*(ln2/t_up)
# dY = (K^n/(K^n + R^n)) - alpha*Y
# initial: R(0)=0, Y(0)=1/alpha (so Y*alpha starts at 1)
# ----------------------------
def ode_rhs(t: float, y: list, beta: float, t_up: float, K: float, n: float, alpha: float):
    """
    Right-hand side of the repression ODE system.

    Parameters
    ----------
    t : float
        Time (hours).
    y : list[float]
        State vector [R, Y], where R≈repressor proxy, Y≈scaled reporter (pre-scaling).
    beta : float
        Production rate chosen so that R rises with half-time t_up (beta=ln2/t_up).
    t_up : float
        Upregulation half-time for R.
    K : float
        Hill half-maximal constant.
    n : float
        Hill coefficient.
    alpha : float
        Reporter decay rate (per hour).

    Returns
    -------
    list[float]
        Derivatives [dR/dt, dY/dt].
    """
    R, Y = y                                                    # unpack state
    dR = beta - R * (math.log(2.0) / t_up)                      # repressor kinetics
    dY = (K**n) / (K**n + max(R, 0.0)**n) - alpha * Y           # reporter kinetics
    return [dR, dY]                                             # return derivatives

def simulate(beta: float, t_up: float, K: float, n: float, alpha: float,
             t_start: float, t_end: float, dt: float = 0.01):
    """
    Simulate the ODE system with LSODA and return (t, R(t), alpha*Y(t)).

    Parameters
    ----------
    beta, t_up, K, n, alpha : float
        Model parameters.
    t_start, t_end : float
        Time window (hours).
    dt : float, optional
        Output sampling step, by default 0.01 h.

    Returns
    -------
    (np.ndarray, np.ndarray, np.ndarray)
        (time grid, R(t), alpha*Y(t)).

    Raises
    ------
    RuntimeError
        If the IVP solver fails.
    """
    t_eval = np.arange(t_start, t_end + dt/2, dt)               # time grid
    y0 = [0.0, 1.0/alpha]                                       # initial state
    print(f"[sim] t∈[{t_start},{t_end}] dt={dt}, params: beta={beta:.4g}, t_up={t_up:.4g}, K={K:.4g}, n={n:.4g}, alpha={alpha:.4g}")  # debug
    sol = solve_ivp(                                            # solve IVP
        ode_rhs, (t_start, t_end), y0,
        t_eval=t_eval, method="LSODA", max_step=0.5,
        args=(beta, t_up, K, n, alpha)
    )
    if not sol.success:                                         # propagate failure
        raise RuntimeError(sol.message)
    R = sol.y[0]                                                # repressor
    Y = sol.y[1] * alpha                                        # scale reporter by alpha
    print(f"[sim] steps={len(sol.t)} R[0..2]={R[:3]} Y[0..2]={Y[:3]}")  # debug
    return sol.t, R, Y                                          # return time series

# ----------------------------
# Delay scan: compute MAE between simulated Y and experimental mean fc.cherry
# for delays in [0,25] step 0.5
# ----------------------------
def delay_scan(pl_df: pd.DataFrame, pars: pd.Series,
               t_grid=(0, 150, 0.01),
               delays=np.arange(0, 25.0 + 1e-9, 0.5)):
    """
    Scan a range of delays (shift in simulation time) to minimize MAE vs. data.

    Parameters
    ----------
    pl_df : pd.DataFrame
        KD dataset for one plasmid with columns ['time', 'fc.cherry', 'norm.bfp'].
    pars : pd.Series
        Parameter row with fields (t_up, k/K, n, alpha).
    t_grid : (float, float, float), optional
        (t0, t1, dt) for base simulation, by default (0, 150, 0.01).
    delays : array-like, optional
        Delay values (hours) to test, by default 0..25 step 0.5.

    Returns
    -------
    (pd.DataFrame, float, float, (np.ndarray, np.ndarray, np.ndarray))
        (delay→MAE table, best_delay, best_mae, (base_t, base_R, base_Y))
    """
    t0, t1, dt = t_grid                                        # unpack time grid

    def pick(series, *keys):                                   # helper to fetch param with aliases
        for k in keys:
            if k in series and pd.notna(series[k]):
                return float(series[k])
        raise KeyError(f"Missing parameter(s) {keys} in row: {series.to_dict()}")

    # Extract parameters with tolerant keys
    t_up  = pick(pars, "t_up", "T_up", "tup")                  # upregulation half-time
    K     = pick(pars, "K", "k")                               # Hill K
    n     = pick(pars, "n", "N", "hill_n")                     # Hill n
    alpha = pick(pars, "alpha", "Alpha", "alphamcherry")       # decay rate
    beta = math.log(2.0) / t_up                                # beta from t_up
    print(f"[delay] params: t_up={t_up:.4g}, K={K:.4g}, n={n:.4g}, alpha={alpha:.4g}, beta={beta:.4g}")  # debug

    # Build experimental mean curve vs time
    exp_mean = (pl_df.groupby("time")["fc.cherry"]             # group by time
                .mean().reset_index().sort_values("time"))     # mean fc per time
    exp_t = exp_mean["time"].to_numpy(float)                   # times
    exp_y = exp_mean["fc.cherry"].to_numpy(float)              # responses
    print(f"[delay] exp points: {len(exp_t)} t[0..3]={exp_t[:3]} y[0..3]={exp_y[:3]}")  # debug

    # Simulate base (no delay)
    base_t, base_R, base_Y = simulate(beta, t_up, K, n, alpha, t0, t1, dt)  # base series

    rows = []                                                  # collect MAE rows
    for d in delays:                                           # iterate candidate delays
        t_del = base_t + d                                     # shifted time grid
        fY = interp1d(t_del, base_Y, kind="linear",            # interpolator on shifted Y
                      bounds_error=False, fill_value=(base_Y[0], base_Y[-1]))
        y_hat = fY(exp_t)                                      # predicted Y at exp times
        resid = exp_y - y_hat                                  # residuals
        N = np.isfinite(resid).sum()                           # number of finite residuals
        mae = np.inf if N <= 1 else np.nansum(np.abs(resid)) / (N - 1)  # robust MAE
        rows.append((float(d), int(N), float(mae)))            # append row

    out = pd.DataFrame(rows, columns=["t", "N", "MAE"])        # delay→MAE table
    best_row = out.loc[out["MAE"].idxmin()]                    # locate min MAE
    best_d, best_mae = float(best_row["t"]), float(best_row["MAE"])  # unpack
    print(f"[delay] best delay={best_d:.3g} h, MAE={best_mae:.4g}")  # debug
    return out, best_d, best_mae, (base_t, base_R, base_Y)     # return results

# ----------------------------
# Plot helpers
# ----------------------------
def save_mae_plot(df_mae: pd.DataFrame, out_pdf: str, y_max: float):
    """
    Save a plot of MAE vs delay and highlight the minimum.

    Parameters
    ----------
    df_mae : pd.DataFrame
        Table with columns ['t', 'N', 'MAE'].
    out_pdf : str
        Output filename (saved under OUT_PATH).
    y_max : float
        Upper y-limit for MAE axis for readability.
    """
    best_row = df_mae.loc[df_mae["MAE"].idxmin()]              # best row for highlight
    p = (                                                       # build plot
        ggplot(df_mae, aes("t", "MAE"))
        + geom_point(size=0.1, alpha=0.4, color="black")
        + geom_point(data=best_row.to_frame().T, mapping=aes("t", "MAE"),
                     size=0.8, alpha=1.0, color=COL_CH)
        + geom_line(data=df_mae.assign(minMAE=df_mae["MAE"].min()),
                    mapping=aes("t", "minMAE"), linetype="dashed")
        + coord_cartesian(xlim=(0, 25), ylim=(0, y_max))
        + labs(y="MAE", x="Delay (hours)")
        + THEME
    )
    out_path = os.path.join(OUT_PATH, out_pdf)                 # full path
    p.save(out_path, width=PLOT_W, height=PLOT_H, units="in")  # save PDF
    print(f"[plot] saved: {out_path}")                         # debug

def save_fit_plots(pl_df: pd.DataFrame, base_t: np.ndarray, base_R: np.ndarray,
                   base_Y: np.ndarray, delay: float, tag_prefix: str):
    """
    Save diagnostic plots comparing experimental data and simulated curves.

    Parameters
    ----------
    pl_df : pd.DataFrame
        KD data for a single plasmid (columns include 'time', 'fc.cherry', 'norm.bfp').
    base_t : np.ndarray
        Base simulation time vector.
    base_R : np.ndarray
        Base simulated R(t).
    base_Y : np.ndarray
        Base simulated alpha*Y(t).
    delay : float
        Best delay (hours) to overlay as dashed line.
    tag_prefix : str
        Tag for filenames (e.g., "KRAB_dCas9").
    """
    # mCherry fit
    p1 = (
        ggplot(pl_df, aes("time", "fc.cherry"))
        + geom_point(size=0.8, alpha=0.4, color=COL_CH)
        + coord_cartesian(xlim=(0, 150), ylim=(0, 1.3))
        + labs(x="Time (hours)", y="mCherry")
        + THEME
        + geom_line(pd.DataFrame({"time": base_t, "Y": base_Y}), aes("time", "Y"))
        + geom_line(pd.DataFrame({"time": base_t + delay, "Y": base_Y}), aes("time", "Y"), linetype="dashed")
    )
    out1 = os.path.join(OUT_PATH, f"KD_ODE_mCherry_{tag_prefix}.pdf")  # filename
    p1.save(out1, width=PLOT_W, height=PLOT_H, units="in")             # save
    print(f"[plot] saved: {out1}")                                     # debug

    # tagBFP fit (R trajectory vs normalized BFP proxy)
    p2 = (
        ggplot(pl_df, aes("time", "norm.bfp"))
        + geom_point(size=0.8, alpha=0.4, color=COL_BFP)
        + coord_cartesian(xlim=(0, 150), ylim=(0, 1.3))
        + labs(x="Time (hours)", y="tagBFP")
        + THEME
        + geom_line(pd.DataFrame({"time": base_t, "R": base_R}), aes("time", "R"))
    )
    out2 = os.path.join(OUT_PATH, f"KD_ODE_tagBFP_{tag_prefix}.pdf")   # filename
    p2.save(out2, width=PLOT_W, height=PLOT_H, units="in")             # save
    print(f"[plot] saved: {out2}")                                     # debug

# ----------------------------
# Main
# ----------------------------
def main():
    # Background from NFC
    mBFP_neg, mmCherry_neg = compute_nfc_background("fcs_files/NFC")
    # KD dataset
    KD = add_fc_and_norm(load_kd(mBFP_neg, mmCherry_neg))

    # Parameters (t_up, K,n, alpha)
    all_par = load_parameters()

    # Targets and pretty names
    target_map = [
        ("SP430ABA", "SP430A",  "KRAB-Split-dCas9"),
        ("SP427",    "SP427",   "HDAC4-dCas9"),
        ("SP428",    "SP428",   "KRAB-dCas9"),
        ("SP411",    "SP411",   "CasRx"),
        ("SP430",    "SP430",   "dCas9"),
    ]

    delay_rows = []

    for pl_raw, pl_label, tag in target_map:
        pl_df = KD[KD["plasmid"] == pl_raw].copy()
        print(f"[target] {pl_raw}: n={len(pl_df)}")
        if pl_df.empty:
            print(f"[WARN] No KD data for {pl_raw}; skipping.")
            continue

        pars = all_par[all_par["plasmid"] == pl_label]
        if pars.empty:
            print(f"[WARN] Missing parameters for {pl_label}; skipping.")
            continue
        pars = pars.iloc[0]

        # delay scan
        mae_df, best_delay, best_mae, (bt, bR, bY) = delay_scan(pl_df, pars)

        # record chosen delay
        delay_rows.append({"plasmid": pl_label, "d_rev": best_delay})
        print(f"[result] {pl_label}: delay={best_delay:.3g}h, MAE={best_mae:.4g}")

        # save MAE plot with reasonable y-limits tuned per case
        ycap = 0.4 if pl_label in ("SP427",) else 0.3
        save_mae_plot(mae_df, f"MAE_KD_{tag}_mcherry.pdf", ycap)

        # save fit plots (mCherry & BFP)
        save_fit_plots(pl_df, bt, bR, bY, best_delay, tag_prefix=tag.replace("-", "_"))

    # write delays
    if delay_rows:
        out = pd.DataFrame(delay_rows, columns=["plasmid", "d_rev"])
        path = os.path.join(PARAM_PATH, "delays_repression.csv")
        out.to_csv(path, index=False)
        print("Saved:", path)
        print(out.to_string(index=False))
    else:
        print("[WARN] No delays estimated.")

if __name__ == "__main__":
    main()

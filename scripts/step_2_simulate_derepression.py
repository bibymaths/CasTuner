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
# This Python implementation reproduces the REV (derepression) step and ODE sims:
#   • Loads .fcs files, applies gating and NFC background subtraction
#   • Computes fc.cherry and norm.bfp transforms for Rev experiments
#   • Fits mCherry degradation rate (alpha) from SP411 time-course
#   • Loads previously fitted half-times / Hill parameters and runs ODE sims
#   • Scans delays by MAE and saves plots and CSV summaries
#
# All thresholds, formulas, and transformations are preserved from the original
# R code to maintain reproducibility and comparability.
# -----------------------------------------------------------------------------

"""
Python port of the R 'REV (derepression)' step.

Pipeline:
1) NFC background from fcs_files/NFC → mBFP_neg, mmCherry_neg
2) Time-course medians from fcs_files/time-course_data with gating & background subtraction
3) Parse filenames → plasmid, exp, rep, time
4) REV (exp == "Rev"): compute fc.cherry (relative to SP411 at 150h) and norm.bfp (min-max via time 0 vs >10h)
5) Fit CasRx (SP411) mCherry half-life → alpha, write parameters/alphamcherry.csv
6) Load parameter CSVs (half_times_downregulation/upregulation, Hill_parameters, alphamcherry)
7) Simulate ODEs for each plasmid (SP411/430/430A/428/427), plot R (tagBFP proxy) & Y (mCherry)
8) Delay scan (0..25h, step 0.5) by MAE on mCherry vs mean per time; write MAE plots; run delayed sims
9) Save mCherry/tagBFP plots; write parameters/delays_derepression.csv
"""

import os, re, glob, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import odeint
from scipy.interpolate import PchipInterpolator
from FlowCytometryTools import FCMeasurement
from plotnine import (
    ggplot, aes, geom_point, geom_line, labs, coord_cartesian,
    theme_classic
)

# ----------------------------
# Config
# ----------------------------
OUT_PATH = "plots"                                   # directory for plots
PARAM_PATH = "parameters"                            # directory for CSV outputs
os.makedirs(OUT_PATH, exist_ok=True)                 # ensure plot dir exists
os.makedirs(PARAM_PATH, exist_ok=True)               # ensure parameter dir exists

CH_FSC_A = "FSC-A"                                   # forward scatter area
CH_SSC_A = "SSC-A"                                   # side scatter area
CH_FSC_H = "FSC-H"                                   # forward scatter height
CH_BFP   = "BV421-A"                                 # tagBFP channel
CH_mCh   = "PE-A"                                    # mCherry channel

BOUND_MIN = {CH_FSC_A: 0.4e5, CH_SSC_A: 0.20e5}      # lower rectangle gate
BOUND_MAX = {CH_FSC_A: 2.0e5, CH_SSC_A: 1.3e5}       # upper rectangle gate
SINGLET_RATIO_LOW  = 0.85                            # FSC-H/FSC-A lower bound
SINGLET_RATIO_HIGH = 1.15                            # FSC-H/FSC-A upper bound

PLOT_W = 1.5 * 1.618                                 # plot width (inches)
PLOT_H = 1.5                                         # plot height (inches)

# ----------------------------
# Helpers: gating & IO
# ----------------------------
def apply_boundary_gate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply coarse FSC/SSC rectangle gate; if empty, fall back to raw.

    Parameters
    ----------
    df : pd.DataFrame
        Events with FSC/SSC channels.

    Returns
    -------
    pd.DataFrame
        Boundary-gated events (or original if gate yields 0).
    """
    m = (                                              # within rectangle mask
        (df[CH_FSC_A] >= BOUND_MIN[CH_FSC_A]) & (df[CH_FSC_A] <= BOUND_MAX[CH_FSC_A]) &
        (df[CH_SSC_A] >= BOUND_MIN[CH_SSC_A]) & (df[CH_SSC_A] <= BOUND_MAX[CH_SSC_A])
    )
    sub = df.loc[m]                                    # apply mask
    print(f"[gate] boundary: kept {len(sub):,}/{len(df):,}")  # debug
    return sub if not sub.empty else df                # fallback to df if empty

def apply_singlet_gate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep singlets using FSC-H/FSC-A ratio; if empty, fall back to input.

    Parameters
    ----------
    df : pd.DataFrame
        Events with FSC-H and FSC-A.

    Returns
    -------
    pd.DataFrame
        Singlet-gated events (or input if gate yields 0).
    """
    ratio = df[CH_FSC_H] / df[CH_FSC_A].replace(0, np.nan)      # robust ratio
    m = (ratio >= SINGLET_RATIO_LOW) & (ratio <= SINGLET_RATIO_HIGH)  # mask
    sub = df.loc[m]                                             # apply mask
    print(f"[gate] singlet : kept {len(sub):,}/{len(df):,}")    # debug
    return sub if not sub.empty else df                         # fallback

def median_channels_for_file(fpath: str) -> pd.Series:
    """
    Read an .fcs file, apply gates, and compute per-channel medians.

    Parameters
    ----------
    fpath : str
        Path to .fcs file.

    Returns
    -------
    pd.Series
        Numeric medians and '__filename' stem.

    Raises
    ------
    ValueError
        If any required channel is missing.
    """
    print(f"[read] {os.path.basename(fpath)}")                  # which file
    dat = FCMeasurement(ID=os.path.basename(fpath), datafile=fpath).data  # load events
    print(f"[read] events={len(dat):,}, cols={len(dat.columns)}")         # size
    for ch in (CH_FSC_A, CH_SSC_A, CH_FSC_H, CH_BFP, CH_mCh):             # verify channels
        if ch not in dat.columns:
            raise ValueError(f"Missing channel '{ch}' in {fpath}")
    gated = apply_singlet_gate(apply_boundary_gate(dat))        # boundary→singlet
    med = gated.median(numeric_only=True)                       # channel medians
    print(f"[median] BFP~{med.get(CH_BFP, np.nan):.3g} | mCh~{med.get(CH_mCh, np.nan):.3g}")  # preview
    med["__filename"] = os.path.splitext(os.path.basename(fpath))[0]  # filename stem
    return med                                                  # return series

def load_flowset_medians(folder: str) -> pd.DataFrame:
    """
    Load all FCS files in a folder and compute gated medians.

    Parameters
    ----------
    folder : str
        Directory containing .fcs files.

    Returns
    -------
    pd.DataFrame
        One row per file with medians + '__filename'.

    Raises
    ------
    FileNotFoundError
        If folder has no .fcs inputs.
    """
    files = sorted(glob.glob(os.path.join(folder, "*.fcs")))    # enumerate files
    print(f"[scan] {folder} → {len(files)} files")              # debug
    if not files:
        raise FileNotFoundError(f"No FCS in {folder}")          # guard
    rows = [median_channels_for_file(f) for f in files]         # compute medians
    out = pd.DataFrame(rows)                                    # to DataFrame
    print(f"[table] medians shape={out.shape}")                 # debug
    return out                                                  # table

def compute_nfc_background(nfc_dir: str):
    """
    NFC background from first 3 medians (match R logic).

    Parameters
    ----------
    nfc_dir : str
        Path to NFC control folder.

    Returns
    -------
    (float, float)
        (mBFP_neg, mmCherry_neg) background medians.
    """
    df = load_flowset_medians(nfc_dir).reset_index(drop=True)   # load NFC
    mBFP_neg = df.loc[:2, CH_BFP].mean()                        # mean BFP (first 3)
    mmCherry_neg = df.loc[:2, CH_mCh].mean()                    # mean mCh (first 3)
    print(f"[NFC] mBFP_neg={mBFP_neg:.6g}, mmCherry_neg={mmCherry_neg:.6g}")  # debug
    return float(mBFP_neg), float(mmCherry_neg)                 # tuple

# ----------------------------
# Name parsing: *_*_*_*_*_*_*
# ----------------------------
def parse_name(name: str):
    """
    Parse filename stem to (plasmid, exp, rep, time).

    Parameters
    ----------
    name : str
        Basename without extension.

    Returns
    -------
    (str, str, str, float)
        (plasmid, exp, rep, time)
    """
    parts = name.split("_")                                     # split tokens
    plasmid = parts[2] if len(parts) > 2 else ""                # plasmid token
    exp     = parts[3] if len(parts) > 3 else ""                # experiment token
    rep     = parts[4] if len(parts) > 4 else ""                # replicate token
    time_s  = parts[5] if len(parts) > 5 else ""                # time token
    try:
        t = float(time_s)                                       # direct cast
    except Exception:
        m = re.search(r"(\d+(\.\d+)?)", time_s)                 # fallback number
        t = float(m.group(1)) if m else np.nan                  # NaN if none
    return plasmid, exp, rep, t                                 # tuple

# ----------------------------
# Time course load (with background subtraction)
# ----------------------------
def load_timecourse_with_bg(mBFP_neg: float, mmCherry_neg: float) -> pd.DataFrame:
    """
    Load time-course medians, subtract NFC backgrounds, and attach parsed metadata.

    Parameters
    ----------
    mBFP_neg : float
        Background for BFP.
    mmCherry_neg : float
        Background for mCherry.

    Returns
    -------
    pd.DataFrame
        Columns: BV421-A, PE-A, plasmid, exp, rep, time
    """
    raw = load_flowset_medians("fcs_files/time-course_data")    # medians table
    df = raw[[CH_BFP, CH_mCh, "__filename"]].copy()             # select needed columns
    df[CH_BFP] = df[CH_BFP] - mBFP_neg                          # subtract BFP bg
    df[CH_mCh] = df[CH_mCh] - mmCherry_neg                      # subtract mCh bg
    print("[tc] bg-sub head:\n" + df[[CH_BFP, CH_mCh]].head().to_string(index=False))  # debug
    parsed = df["__filename"].apply(parse_name)                 # parse tokens
    meta = pd.DataFrame(parsed.tolist(), columns=["plasmid", "exp", "rep", "time"])    # to DataFrame
    df = pd.concat([df.drop(columns="__filename"), meta], axis=1)                      # combine
    print(f"[tc] combined shape={df.shape}")                     # debug
    return df                                                   # return table

# ----------------------------
# CasRx α (mCherry degradation) from REV SP411
# ----------------------------
def exp_recovery(t, t12, y0, yf):
    """
    Exponential recovery curve: y(t) = yf + (y0 - yf) * exp(-t * ln2 / t12)

    Parameters
    ----------
    t : array-like
        Time values (hours).
    t12 : float
        Half-time.
    y0 : float
        Initial value at t=0.
    yf : float
        Final asymptote.

    Returns
    -------
    np.ndarray
        Modeled y(t).
    """
    return yf + (y0 - yf) * np.exp(-t * (np.log(2)/t12))        # vectorized formula

def fit_alpha_from_casrx(REV: pd.DataFrame) -> float:
    """
    Fit mCherry half-life from SP411 Rev data and convert to alpha=ln2/t12.

    Parameters
    ----------
    REV : pd.DataFrame
        Rev-only table with columns ['plasmid','time','fc.cherry'].

    Returns
    -------
    float
        Alpha (1/h). Also writes alphamcherry.csv (rounded to 3 decimals).

    Raises
    ------
    RuntimeError
        If SP411 records are missing.
    """
    sp = REV[REV["plasmid"] == "SP411"].copy()                  # filter SP411
    if sp.empty:
        raise RuntimeError("No SP411 in REV for alpha fitting") # guard
    sp = sp.sort_values("time")                                 # ensure time order
    t = sp["time"].to_numpy(float)                              # time vector
    y = sp["fc.cherry"].to_numpy(float)                         # response vector
    yf = 1.0                                                    # final level
    y0 = float(sp.loc[sp["time"] == 0, "fc.cherry"].mean())     # initial level
    print(f"[alpha] y0={y0:.6g}, yf={yf:.6g}, N={len(t)}")      # debug
    popt, _ = curve_fit(lambda tt, t12: exp_recovery(tt, t12, y0=y0, yf=yf),
                        t, y, p0=[0.1], maxfev=20000)           # fit t12
    t12_hat = float(popt[0])                                    # estimated half-time
    alpha = float(np.log(2)/t12_hat)                            # convert to alpha
    print(f"[alpha] t1/2={t12_hat:.6g} h → alpha={alpha:.6g} 1/h")  # debug
    pd.DataFrame({"alpha": [round(alpha, 3)]}).to_csv(          # save CSV
        os.path.join(PARAM_PATH, "alphamcherry.csv"), index=False
    )
    return alpha                                                # return alpha

# ----------------------------
# Parameter load/merge
# ----------------------------
def load_parameters():
    """
    Load fitted parameter CSVs and merge into one table; broadcast alpha.

    Returns
    -------
    pd.DataFrame
        Columns: plasmid, t_down, t_up, K, n, alpha

    Raises
    ------
    KeyError
        If required columns are missing in any input CSV.
    """
    t_down = pd.read_csv(os.path.join(PARAM_PATH, "half_times_downregulation.csv"))  # load decay half-times
    if "se" in t_down.columns: t_down = t_down.drop(columns=["se"])                  # drop SE if present
    t_down = t_down.rename(columns={"halftime": "t_down"})                           # rename

    t_up = pd.read_csv(os.path.join(PARAM_PATH, "half_times_upregulation.csv"))      # load rise half-times
    if "se" in t_up.columns: t_up = t_up.drop(columns=["se"])
    t_up = t_up.rename(columns={"halftime": "t_up"})

    hill = pd.read_csv(os.path.join(PARAM_PATH, "Hill_parameters.csv"))              # K, n
    alpha = pd.read_csv(os.path.join(PARAM_PATH, "alphamcherry.csv"))                # alpha scalar

    if "plasmid" not in t_down.columns or "plasmid" not in t_up.columns or "plasmid" not in hill.columns:
        raise KeyError("t_down/t_up/Hill_parameters must have 'plasmid' column")     # schema guard

    df = t_down.merge(t_up, on="plasmid", how="outer").merge(hill, on="plasmid", how="outer")  # merge all
    if "plasmid" not in df.columns:
        raise KeyError("Merged parameter table missing 'plasmid'")                   # sanity guard
    if "alpha" not in alpha.columns:
        raise KeyError("alphamcherry.csv must have column 'alpha' (single value)")   # schema guard

    alpha_val = float(alpha["alpha"].iloc[0])                                        # extract scalar
    df["alpha"] = alpha_val                                                          # broadcast to rows
    print(f"[pars] merged shape={df.shape}; alpha={alpha_val:.6g}")                  # debug
    return df                                                                         # parameters table

# ----------------------------
# REV transforms: fc.cherry & norm.bfp
# ----------------------------
def compute_rev_transforms(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute Rev-only transforms:
      • fc.cherry = PE-A / mean(PE-A of SP411 at t=150)
      • norm.bfp  = (BFP - mean.init) / (mean.final - mean.init)

    Parameters
    ----------
    df : pd.DataFrame
        Full time-course with background-subtracted BFP/mCherry.

    Returns
    -------
    pd.DataFrame
        Rev-only table with columns fc.cherry and norm.bfp added.

    Raises
    ------
    RuntimeError
        If Rev data or SP411@150 reference is missing.
    """
    REV = df[df["exp"] == "Rev"].copy()                         # keep Rev only
    if REV.empty:
        raise RuntimeError("No 'Rev' records in time-course")   # guard

    sp411_150 = REV[(REV["plasmid"] == "SP411") & (REV["time"] == 150)]  # reference rows
    if sp411_150.empty:
        raise RuntimeError("No SP411 at 150h to define mCherry reference")  # guard
    mmcherry = float(sp411_150[CH_mCh].mean())                  # reference mean
    REV["fc.cherry"] = REV[CH_mCh] / mmcherry                   # compute fold change
    print(f"[rev] mCherry ref (SP411@150h)={mmcherry:.6g}")     # debug

    finals = (REV[REV["time"] > 10]                             # mean.final per plasmid
              .groupby("plasmid")[CH_BFP]
              .mean().rename("mean.final").reset_index())
    inits = (REV[REV["time"] == 0]                              # mean.init per plasmid
             .groupby("plasmid")[CH_BFP]
             .mean().rename("mean.init").reset_index())
    REV = REV.merge(finals, on="plasmid", how="left").merge(inits, on="plasmid", how="left")  # attach means
    REV["norm.bfp"] = (REV[CH_BFP] - REV["mean.init"]) / (REV["mean.final"] - REV["mean.init"])  # min–max
    print("[rev] transforms preview:\n" +                       # preview few rows
          REV[["plasmid","time",CH_BFP,"mean.init","mean.final","norm.bfp","fc.cherry"]]
          .head().to_string(index=False))
    return REV                                                  # return transformed table

# ----------------------------
# ODE system and simulators
# ----------------------------
def rhs(y, t, t_down, K, n, alpha):
    """
    Right-hand side of ODE system:
      dR/dt = -R * ln2 / t_down
      dY/dt = (K^n)/(K^n + R^n) - alpha * Y

    Parameters
    ----------
    y : array-like
        State [R, Y/alpha] (internal scaling for numerical stability).
    t : float
        Time (hours).
    t_down : float
        Half-time for R decay (hours).
    K : float
        Hill constant.
    n : float
        Hill coefficient.
    alpha : float
        mCherry degradation rate (1/h).

    Returns
    -------
    list[float]
        Derivatives [dR, d(Y/alpha)].
    """
    R, Y = float(y[0]), float(y[1])                             # unpack state
    t_down = max(float(t_down), 1e-6)                           # avoid zero
    K = max(float(K), 1e-12)                                    # avoid zero
    n = max(float(n), 1e-6)                                     # avoid zero
    alpha = max(float(alpha), 1e-12)                            # avoid zero

    Rpos = max(R, 0.0)                                          # nonnegative R
    Kn = K**n                                                   # precompute
    Rn = Rpos**n                                                # precompute
    dR = -Rpos * (np.log(2) / t_down)                           # R decay
    dY = (Kn)/(Kn + Rn) - Y * alpha                             # Y/alpha dynamics
    return [dR, dY]                                             # derivatives

def simulate_ode(R0: float, Y0: float, pars: pd.Series,
                 t0: float = 0.0, tmax: float = 150.0, step: float = 0.05,
                 delay: float = 0.0) -> pd.DataFrame:
    """
    Integrate ODEs and return time series for R and Y.

    Notes
    -----
    Internally integrates Y/alpha for numerical stability; rescales on output.

    Parameters
    ----------
    R0 : float
        Initial normalized repressor level at t=0.
    Y0 : float
        Initial mCherry fold-change at t=0.
    pars : pd.Series | dict
        Parameters with keys: t_down, K, n, alpha.
    t0 : float, optional
        Start time, by default 0.0.
    tmax : float, optional
        End time, by default 150.0.
    step : float, optional
        Time step for output grid, by default 0.05.
    delay : float, optional
        Display-only shift applied to the returned 'time' column.

    Returns
    -------
    pd.DataFrame
        Columns: time, R, Y
    """
    y0 = [max(float(R0), 0.0), max(float(Y0), 0.0) / float(pars["alpha"])]  # scaled ICs
    t_eval = np.arange(t0, tmax + step/2, step)                   # time grid
    print(f"[ode] R0={R0:.4g}, Y0={Y0:.4g}, pars={dict(pars)}, delay={delay}")  # debug

    sol = odeint(                                                 # integrate
        rhs, y0, t_eval,
        args=(pars["t_down"], pars["K"], pars["n"], pars["alpha"]),
        atol=1e-10, rtol=1e-9,
    )

    R = np.maximum(sol[:, 0], 0.0)                                # nonnegative R
    Y = np.maximum(sol[:, 1] * float(pars["alpha"]), 0.0)         # rescale Y
    df = pd.DataFrame({"time": t_eval + float(delay), "R": R, "Y": Y})  # assemble
    return df                                                     # series

# ----------------------------
# Delay scan MAE vs mean(fc.cherry by time)
# ----------------------------
def delay_scan(pl_df: pd.DataFrame, pars: dict, tmax=150.0, step=0.005):
    """
    Scan delays (0..25h) to minimize MAE between simulated Y and mean observed mCherry.

    Parameters
    ----------
    pl_df : pd.DataFrame
        Rev data for one plasmid with 'time','fc.cherry','norm.bfp'.
    pars : dict
        Parameters dict with keys 't_down','K','n','alpha'.
    tmax : float, optional
        Simulation horizon, by default 150.0.
    step : float, optional
        Simulation time step, by default 0.005.

    Returns
    -------
    (pd.DataFrame, float|None, float|None, pd.DataFrame|None)
        (mae_table, best_delay, best_mae, shifted_sim_at_best)
    """
    mean_fc = (pl_df.groupby("time")["fc.cherry"]                # mean per time
               .mean().rename("m_fc").reset_index().sort_values("time"))
    delays = np.arange(0.0, 25.0 + 1e-9, 0.5)                    # scan grid
    rows, best = [], (None, np.inf, None)                        # init best

    R0 = float(pl_df.loc[pl_df["time"] == 0, "norm.bfp"].mean()) # initial R
    Y0 = float(pl_df.loc[pl_df["time"] == 0, "fc.cherry"].mean())# initial Y
    sim = simulate_ode(R0, Y0, pars, tmax=tmax, step=step)       # base sim

    pchipY = PchipInterpolator(sim["time"].to_numpy(),           # monotone interp
                               sim["Y"].to_numpy(), extrapolate=False)

    for d in delays:                                             # iterate delays
        sel = mean_fc["time"] >= d                               # valid rows
        if not np.any(sel):
            continue
        t_data = mean_fc.loc[sel, "time"].to_numpy()             # observed times
        t_query = t_data - d                                     # shift for model
        in_dom  = (t_query >= sim["time"].iloc[0]) & (t_query <= sim["time"].iloc[-1])  # support
        if not np.any(in_dom):
            continue

        y_pred = np.full_like(t_query, np.nan, dtype=float)      # alloc
        y_pred[in_dom] = pchipY(t_query[in_dom])                 # predict
        df_cmp = pd.DataFrame({"time": t_data,
                               "m_fc": mean_fc.loc[sel, "m_fc"].to_numpy(),
                               "Y": y_pred}).dropna(subset=["Y"])  # compare
        N = len(df_cmp)                                          # count
        if N <= 1:
            continue
        mae = (df_cmp["m_fc"] - df_cmp["Y"]).abs().sum() / (N - 1)  # MAE (R-like)
        rows.append({"t": d, "MAE": mae})                        # record

        if mae < best[1]:                                        # update best
            shifted = sim.copy(); shifted["time"] = shifted["time"] + d
            best = (d, mae, shifted)

    mae_df = pd.DataFrame(rows)                                  # full table
    print(f"[delay] scanned {len(mae_df)} delays; best={best[0]} h, MAE={best[1]:.6g}")  # debug
    return mae_df, best[0], best[1], best[2]                     # results

# ----------------------------
# Plot helpers
# ----------------------------
def p_theme():
    """Return a clean, classic theme."""
    return theme_classic()

def save_scatter_with_lines(df_pts, x, y, lines=None, fname="plot.pdf",
                            xlim=None, ylim=None, y_label="", x_label="Time (hours)"):
    """
    Scatter of observed points with optional model lines; save to PDF.

    Parameters
    ----------
    df_pts : pd.DataFrame
        Points to plot.
    x, y : str
        Column names for x and y.
    lines : list[tuple] | None
        Optional [(df_line, aes_y, linetype), ...].
    fname : str
        Output filename under OUT_PATH.
    xlim, ylim : tuple | None
        Axis limits, if any.
    y_label, x_label : str
        Axis labels.
    """
    p = (
        ggplot(df_pts, aes(x=x, y=y))
        + geom_point(size=0.6, alpha=0.4)
        + labs(x=x_label, y=y_label)
        + p_theme()
    )
    if xlim or ylim:
        p = p + coord_cartesian(xlim=xlim, ylim=ylim)            # set limits
    if lines is not None:
        for (df_line, aes_y, lty) in lines:                      # add each line
            p = p + geom_line(data=df_line, mapping=aes(x="time", y=aes_y), linetype=lty)
    out_path = os.path.join(OUT_PATH, fname)                     # build path
    p.save(out_path, width=PLOT_W, height=PLOT_H, units="in")    # save plot
    print(f"[plot] saved: {out_path}")                           # debug

def save_mae_plot(mae_df, best_delay, fname, ylim=(0, 0.3)):
    """
    Save MAE vs delay scatter with highlighted best point.

    Parameters
    ----------
    mae_df : pd.DataFrame
        Columns ['t','MAE'].
    best_delay : float
        Best delay to highlight.
    fname : str
        Output filename under OUT_PATH.
    ylim : tuple
        Y-axis limits.
    """
    if mae_df.empty:
        print("[plot] MAE table empty; skipping MAE plot.")      # guard
        return
    y_min = mae_df["MAE"].min()                                  # baseline line
    p = (
        ggplot(mae_df, aes(x="t", y="MAE"))
        + geom_point(size=0.1, alpha=0.4)
        + geom_point(data=mae_df[mae_df["t"] == best_delay], mapping=aes(x="t", y="MAE"),
                     size=0.8, alpha=1.0)
        + labs(x="Delay (hours)", y="MAE")
        + coord_cartesian(xlim=(0, 25), ylim=ylim)
        + p_theme()
    )
    p = p + geom_line(data=pd.DataFrame({"t": mae_df["t"], "MAE": y_min}),
                      mapping=aes(x="t", y="MAE"), linetype="dashed")  # min line
    out_path = os.path.join(OUT_PATH, fname)                    # output path
    p.save(out_path, width=PLOT_W, height=PLOT_H, units="in")   # save plot
    print(f"[plot] saved: {out_path}")                          # debug

# ----------------------------
# Main
# ----------------------------
def main():
    # 1) Background from NFC
    mBFP_neg, mmCherry_neg = compute_nfc_background("fcs_files/NFC")

    # 2) Time-course with background subtraction
    tc = load_timecourse_with_bg(mBFP_neg, mmCherry_neg)

    # 3) REV transforms
    REV = compute_rev_transforms(tc)

    # 4) Fit mCherry alpha from CasRx (SP411)
    alpha_val = fit_alpha_from_casrx(REV)
    print(f"[INFO] alpha (mCherry) = {alpha_val:.4f}  (1/h)")

    # 5) Parameters
    all_par = load_parameters()

    # 6) Simulate each plasmid without delay and save plots
    #    Order: SP411 (CasRx), SP430 (dCas9), SP430A (KRAB-Split-dCas9), SP428 (KRAB-dCas9), SP427 (HDAC4-dCas9)
    plasmids = ["SP411", "SP430", "SP430A", "SP428", "SP427"]
    for pl in plasmids:
        pars_row = all_par[all_par["plasmid"] == pl]
        if pars_row.empty:
            print(f"[WARN] No params for {pl}, skipping.")
            continue
        pars = {k: float(pars_row.iloc[0][k]) for k in ["t_down", "K", "n", "alpha"]}

        sub = REV[REV["plasmid"] == ("SP430ABA" if pl == "SP430A" else pl)].copy()
        if sub.empty:
            print(f"[WARN] No REV data for {pl}, skipping.")
            continue

        R0 = float(sub.loc[sub["time"] == 0, "norm.bfp"].mean())
        Y0 = float(sub.loc[sub["time"] == 0, "fc.cherry"].mean())

        sim = simulate_ode(R0, Y0, pars)

        # mCherry
        save_scatter_with_lines(
            df_pts=sub, x="time", y="fc.cherry",
            lines=[(sim, "Y", "solid")],
            fname=f"REV_mCherry_{pl}_Hill.pdf",
            xlim=(0, 150), ylim=(0, 1.3),
            y_label="mCherry"
        )
        # tagBFP proxy (R)
        save_scatter_with_lines(
            df_pts=sub, x="time", y="norm.bfp",
            lines=[(sim, "R", "solid")],
            fname=f"REV_tagBFP_{pl}_Hill.pdf",
            xlim=(0, 150), ylim=(0, 1.3),
            y_label="tagBFP"
        )

    # 7) Delay scan & delayed sims for SP430A/SP428/SP427/SP430 (match R flow)
    delay_rows = []
    scans = [
        ("SP430A", "SP430ABA"),  # parameter label → data label
        ("SP428",  "SP428"),
        ("SP427",  "SP427"),
        ("SP430",  "SP430"),
    ]
    for pl_param, pl_data in scans:
        pars_row = all_par[all_par["plasmid"] == pl_param]
        if pars_row.empty:
            print(f"[WARN] No params for {pl_param}, skipping delay scan.")
            continue
        pars = {k: float(pars_row.iloc[0][k]) for k in ["t_down", "K", "n", "alpha"]}

        sub = REV[REV["plasmid"] == pl_data].copy()
        if sub.empty:
            print(f"[WARN] No REV data for {pl_data}, skipping delay scan.")
            continue

        mae_df, best_delay, best_mae, sim_best = delay_scan(sub, pars)
        print(f"[INFO] Best delay {pl_param}: {best_delay} h (MAE={best_mae:.4f})")
        delay_rows.append({"plasmid": pl_param, "d_rev": best_delay})

        # MAE plot
        save_mae_plot(mae_df, best_delay, fname=f"MAE_REV_{pl_param}_mcherry.pdf")

        # Undelayed sim for comparison
        R0 = float(sub.loc[sub["time"] == 0, "norm.bfp"].mean())
        Y0 = float(sub.loc[sub["time"] == 0, "fc.cherry"].mean())
        sim_no_delay = simulate_ode(R0, Y0, pars)

        # mCherry delayed vs undelayed
        save_scatter_with_lines(
            df_pts=sub, x="time", y="fc.cherry",
            lines=[(sim_best, "Y", "dashed"), (sim_no_delay, "Y", "solid")],
            fname=f"REV_mCherry_{pl_param}_{best_delay}h_delay.pdf",
            xlim=(0, 150), ylim=(0, 1.3),
            y_label="mCherry"
        )
        # tagBFP (R) undelayed (to mirror R figs)
        save_scatter_with_lines(
            df_pts=sub, x="time", y="norm.bfp",
            lines=[(sim_no_delay, "R", "solid")],
            fname=f"REV_tagBFP_{pl_param}_{best_delay}h_delay.pdf",
            xlim=(0, 150), ylim=(0, 1.3),
            y_label=""
        )

    # SP411 had no delay scan in R; its non-delayed plots were saved above.

    # 8) Save delays table
    if delay_rows:
        out_csv = os.path.join(PARAM_PATH, "delays_derepression.csv")
        pd.DataFrame(delay_rows).to_csv(out_csv, index=False)
        print("\nSaved:", out_csv)
        print(pd.DataFrame(delay_rows).to_string(index=False))
    else:
        print("[WARN] No delays computed; no CSV written.")

if __name__ == "__main__":
    main()

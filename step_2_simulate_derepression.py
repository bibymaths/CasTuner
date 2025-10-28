#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

Requirements:
  pandas, numpy, scipy, plotnine
  FlowCytometryTools (plus flowio/flowutils)
"""
import math
import os
import re
import glob
import warnings
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
OUT_PATH = "plots"
PARAM_PATH = "parameters"
os.makedirs(OUT_PATH, exist_ok=True)
os.makedirs(PARAM_PATH, exist_ok=True)

CH_FSC_A = "FSC-A"
CH_SSC_A = "SSC-A"
CH_FSC_H = "FSC-H"
CH_BFP   = "BV421-A"
CH_mCh   = "PE-A"

BOUND_MIN = {CH_FSC_A: 0.4e5, CH_SSC_A: 0.20e5}
BOUND_MAX = {CH_FSC_A: 2.0e5, CH_SSC_A: 1.3e5}
SINGLET_RATIO_LOW  = 0.85
SINGLET_RATIO_HIGH = 1.15

PLOT_W = 1.5 * 1.618
PLOT_H = 1.5

# ----------------------------
# Helpers: gating & IO
# ----------------------------
def apply_boundary_gate(df: pd.DataFrame) -> pd.DataFrame:
    m = (
        (df[CH_FSC_A] >= BOUND_MIN[CH_FSC_A]) & (df[CH_FSC_A] <= BOUND_MAX[CH_FSC_A]) &
        (df[CH_SSC_A] >= BOUND_MIN[CH_SSC_A]) & (df[CH_SSC_A] <= BOUND_MAX[CH_SSC_A])
    )
    sub = df.loc[m]
    return sub if not sub.empty else df

def apply_singlet_gate(df: pd.DataFrame) -> pd.DataFrame:
    ratio = df[CH_FSC_H] / df[CH_FSC_A].replace(0, np.nan)
    m = (ratio >= SINGLET_RATIO_LOW) & (ratio <= SINGLET_RATIO_HIGH)
    sub = df.loc[m]
    return sub if not sub.empty else df

def median_channels_for_file(fpath: str) -> pd.Series:
    dat = FCMeasurement(ID=os.path.basename(fpath), datafile=fpath).data
    for ch in (CH_FSC_A, CH_SSC_A, CH_FSC_H, CH_BFP, CH_mCh):
        if ch not in dat.columns:
            raise ValueError(f"Missing channel '{ch}' in {fpath}")
    gated = apply_singlet_gate(apply_boundary_gate(dat))
    med = gated.median(numeric_only=True)
    med["__filename"] = os.path.splitext(os.path.basename(fpath))[0]
    return med

def load_flowset_medians(folder: str) -> pd.DataFrame:
    files = sorted(glob.glob(os.path.join(folder, "*.fcs")))
    if not files:
        raise FileNotFoundError(f"No FCS in {folder}")
    rows = [median_channels_for_file(f) for f in files]
    return pd.DataFrame(rows)

def compute_nfc_background(nfc_dir: str):
    df = load_flowset_medians(nfc_dir).reset_index(drop=True)
    mBFP_neg = df.loc[:2, CH_BFP].mean()
    mmCherry_neg = df.loc[:2, CH_mCh].mean()
    return float(mBFP_neg), float(mmCherry_neg)

# ----------------------------
# Name parsing: *_*_*_*_*_*_*
# Expect tokens: (..)_ (..)_ (plasmid)_ (exp)_ (rep)_ (time) _ .. _
# ----------------------------
def parse_name(name: str):
    parts = name.split("_")
    plasmid = parts[2] if len(parts) > 2 else ""
    exp     = parts[3] if len(parts) > 3 else ""
    rep     = parts[4] if len(parts) > 4 else ""
    time_s  = parts[5] if len(parts) > 5 else ""
    # time numeric
    try:
        t = float(time_s)
    except:
        # fallback: capture trailing digits
        m = re.search(r"(\d+(\.\d+)?)", time_s)
        t = float(m.group(1)) if m else np.nan
    return plasmid, exp, rep, t

# ----------------------------
# Time course load (with background subtraction)
# ----------------------------
def load_timecourse_with_bg(mBFP_neg: float, mmCherry_neg: float) -> pd.DataFrame:
    raw = load_flowset_medians("fcs_files/time-course_data")
    df = raw[[CH_BFP, CH_mCh, "__filename"]].copy()
    df[CH_BFP] = df[CH_BFP] - mBFP_neg
    df[CH_mCh] = df[CH_mCh] - mmCherry_neg
    parsed = df["__filename"].apply(parse_name)
    meta = pd.DataFrame(parsed.tolist(), columns=["plasmid", "exp", "rep", "time"])
    df = pd.concat([df.drop(columns="__filename"), meta], axis=1)
    return df

# ----------------------------
# CasRx α (mCherry degradation) from REV SP411
# y(t) = yf + (y0 - yf) * exp(-t * (log(2)/t12))
# Here we fit t12, then alpha = log(2)/t12
# ----------------------------
def exp_recovery(t, t12, y0, yf):
    return yf + (y0 - yf) * np.exp(-t * (np.log(2)/t12))

def fit_alpha_from_casrx(REV: pd.DataFrame) -> float:
    sp = REV[REV["plasmid"] == "SP411"].copy()
    if sp.empty:
        raise RuntimeError("No SP411 in REV for alpha fitting")
    # y = fc.cherry, t = time
    sp = sp.sort_values("time")
    t = sp["time"].to_numpy(float)
    y = sp["fc.cherry"].to_numpy(float)
    yf = 1.0
    # y0 = mean at t==0 (as in R)
    y0 = float(sp.loc[sp["time"] == 0, "fc.cherry"].mean())
    # initial guess t12=0.1
    popt, _ = curve_fit(lambda tt, t12: exp_recovery(tt, t12, y0=y0, yf=yf),
                        t, y, p0=[0.1], maxfev=20000)
    t12_hat = float(popt[0])
    alpha = float(np.log(2)/t12_hat)
    # write to CSV like R
    pd.DataFrame({"alpha": [round(alpha, 3)]}).to_csv(
        os.path.join(PARAM_PATH, "alphamcherry.csv"), index=False
    )
    return alpha

# ----------------------------
# Parameter load/merge
# ----------------------------
def load_parameters():
    t_down = pd.read_csv(os.path.join(PARAM_PATH, "half_times_downregulation.csv"))
    if "se" in t_down.columns: t_down = t_down.drop(columns=["se"])
    t_down = t_down.rename(columns={"halftime": "t_down"})

    t_up = pd.read_csv(os.path.join(PARAM_PATH, "half_times_upregulation.csv"))
    if "se" in t_up.columns: t_up = t_up.drop(columns=["se"])
    t_up = t_up.rename(columns={"halftime": "t_up"})

    hill = pd.read_csv(os.path.join(PARAM_PATH, "Hill_parameters.csv"))
    alpha = pd.read_csv(os.path.join(PARAM_PATH, "alphamcherry.csv"))

    # Expected columns:
    # t_down: plasmid, t_down
    # t_up:   plasmid, t_up
    # hill:   plasmid, K, n
    # alpha:  alpha (single scalar)
    if "plasmid" not in t_down.columns or "plasmid" not in t_up.columns or "plasmid" not in hill.columns:
        raise KeyError("t_down/t_up/Hill_parameters must have 'plasmid' column")

    df = t_down.merge(t_up, on="plasmid", how="outer").merge(hill, on="plasmid", how="outer")
    if "plasmid" not in df.columns:
        raise KeyError("Merged parameter table missing 'plasmid'")

    if "alpha" not in alpha.columns:
        raise KeyError("alphamcherry.csv must have column 'alpha' (single value)")
    alpha_val = float(alpha["alpha"].iloc[0])
    df["alpha"] = alpha_val
    return df

# ----------------------------
# REV transforms: fc.cherry & norm.bfp
# ----------------------------
def compute_rev_transforms(df: pd.DataFrame) -> pd.DataFrame:
    REV = df[df["exp"] == "Rev"].copy()
    if REV.empty:
        raise RuntimeError("No 'Rev' records in time-course")

    # mCherry reference = mean PE-A of SP411 at time==150
    sp411_150 = REV[(REV["plasmid"] == "SP411") & (REV["time"] == 150)]
    if sp411_150.empty:
        raise RuntimeError("No SP411 at 150h to define mCherry reference")
    mmcherry = float(sp411_150[CH_mCh].mean())
    REV["fc.cherry"] = REV[CH_mCh] / mmcherry

    # mean-max BFP scaling: mean(final) at time>10, mean(init) at time==0
    finals = (REV[REV["time"] > 10]
              .groupby("plasmid")[CH_BFP]
              .mean().rename("mean.final").reset_index())
    inits = (REV[REV["time"] == 0]
             .groupby("plasmid")[CH_BFP]
             .mean().rename("mean.init").reset_index())
    REV = REV.merge(finals, on="plasmid", how="left").merge(inits, on="plasmid", how="left")
    REV["norm.bfp"] = (REV[CH_BFP] - REV["mean.init"]) / (REV["mean.final"] - REV["mean.init"])

    return REV

# ----------------------------
# ODE system and simulators
# dR/dt = -R * (log(2)/t12)                 (derepression)
# dY/dt = (K^n)/(K^n + R^n) - alpha * Y
# ----------------------------
def rhs(y, t, t_down, K, n, alpha):
    # y = [R, Y]
    R, Y = float(y[0]), float(y[1])
    t_down = max(float(t_down), 1e-6)
    K = max(float(K), 1e-12)
    n = max(float(n), 1e-6)
    alpha = max(float(alpha), 1e-12)

    Rpos = max(R, 0.0)
    Kn = K**n
    Rn = Rpos**n
    denom = Kn + Rn + 1e-15

    dR = -Rpos * (np.log(2) / t_down)
    dY = (Kn)/(Kn + Rn) - Y * alpha
    return [dR, dY]


def simulate_ode(R0: float, Y0: float, pars: pd.Series,
                 t0: float = 0.0, tmax: float = 150.0, step: float = 0.05,
                 delay: float = 0.0) -> pd.DataFrame:
    """
    Return a DataFrame with columns: time, R, Y.
    Delay is applied as a shift to the reported 'time' axis (ICs unchanged).
    """
    y0 = [max(float(R0), 0.0), max(float(Y0), 0.0) / float(pars["alpha"])]
    t_eval = np.arange(t0, tmax + step/2, step)

    sol = odeint(
        rhs,
        y0,
        t_eval,
        args=(pars["t_down"], pars["K"], pars["n"], pars["alpha"]),
        atol=1e-10,
        rtol=1e-9,
        # hmax=0.5,
        # mxstep=10000
    )

    R = np.maximum(sol[:, 0], 0.0)
    Y = np.maximum(sol[:, 1] * float(pars["alpha"]), 0.0)  # rescale back

    df = pd.DataFrame({
        "time": t_eval + float(delay),
        "R": R,
        "Y": Y
    })
    return df

# ----------------------------
# Delay scan MAE vs mean(fc.cherry by time)
# ----------------------------
def delay_scan(pl_df: pd.DataFrame, pars: dict, tmax=150.0, step=0.005):
    mean_fc = (pl_df.groupby("time")["fc.cherry"]
               .mean().rename("m_fc").reset_index().sort_values("time"))
    delays = np.arange(0.0, 25.0 + 1e-9, 0.5)
    rows = []
    best = (None, np.inf, None)

    R0 = float(pl_df.loc[pl_df["time"] == 0, "norm.bfp"].mean())
    Y0 = float(pl_df.loc[pl_df["time"] == 0, "fc.cherry"].mean())

    # t_end = mean_fc["time"].max()
    sim = simulate_ode(R0, Y0, pars, tmax=tmax, step=step)

    # monotone interpolator on Y vs time
    pchipY = PchipInterpolator(sim["time"].to_numpy(), sim["Y"].to_numpy(), extrapolate=False)

    for d in delays:
        sel = mean_fc["time"] >= d
        if not np.any(sel):
            continue
        t_data = mean_fc.loc[sel, "time"].to_numpy()

        # evaluate only within support; mask out those beyond sim range
        t_query = t_data - d
        in_dom  = (t_query >= sim["time"].iloc[0]) & (t_query <= sim["time"].iloc[-1])
        if not np.any(in_dom):
            continue

        y_pred = np.full_like(t_query, np.nan, dtype=float)
        y_pred[in_dom] = pchipY(t_query[in_dom])

        df_cmp = pd.DataFrame({"time": t_data, "m_fc": mean_fc.loc[sel, "m_fc"].to_numpy(), "Y": y_pred})
        df_cmp = df_cmp.dropna(subset=["Y"])
        N = len(df_cmp)
        if N <= 1:
            continue
        mae = (df_cmp["m_fc"] - df_cmp["Y"]).abs().sum() / (N - 1)
        rows.append({"t": d, "MAE": mae})

        if mae < best[1]:
            shifted = sim.copy()
            shifted["time"] = shifted["time"] + d
            best = (d, mae, shifted)

    mae_df = pd.DataFrame(rows)
    return mae_df, best[0], best[1], best[2]

# ----------------------------
# Plot helpers
# ----------------------------
def p_theme():
    return theme_classic()

def save_scatter_with_lines(df_pts, x, y, lines=None, fname="plot.pdf",
                            xlim=None, ylim=None, y_label="", x_label="Time (hours)"):
    p = (
        ggplot(df_pts, aes(x=x, y=y))
        + geom_point(size=0.6, alpha=0.4)
        + labs(x=x_label, y=y_label)
        + p_theme()
    )
    if xlim or ylim:
        p = p + coord_cartesian(xlim=xlim, ylim=ylim)
    if lines is not None:
        # lines: list of (df_line, aes_y, linetype)
        for (df_line, aes_y, lty) in lines:
            p = p + geom_line(data=df_line, mapping=aes(x="time", y=aes_y), linetype=lty)
    p.save(os.path.join(OUT_PATH, fname), width=PLOT_W, height=PLOT_H, units="in")

def save_mae_plot(mae_df, best_delay, fname, ylim=(0, 0.3)):
    if mae_df.empty:
        return
    y_min = mae_df["MAE"].min()
    p = (
        ggplot(mae_df, aes(x="t", y="MAE"))
        + geom_point(size=0.1, alpha=0.4)
        + geom_point(data=mae_df[mae_df["t"] == best_delay], mapping=aes(x="t", y="MAE"), size=0.8, alpha=1.0)
        + labs(x="Delay (hours)", y="MAE")
        + coord_cartesian(xlim=(0, 25), ylim=ylim)
        + p_theme()
    )
    # horizontal line at min MAE:
    p = p + geom_line(data=pd.DataFrame({"t": mae_df["t"], "MAE": y_min}), mapping=aes(x="t", y="MAE"), linetype="dashed")
    p.save(os.path.join(OUT_PATH, fname), width=PLOT_W, height=PLOT_H, units="in")

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

    # 7) Delay scan & delayed sims for SP430A/SP428/SP427/SP430/SP411 (match R flow)
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

    # SP411 had no delay scan in R (it’s just derepression to 0 quickly), but we still saved non-delayed above.

    # 8) Save delays table
    if delay_rows:
        pd.DataFrame(delay_rows).to_csv(os.path.join(PARAM_PATH, "delays_derepression.csv"), index=False)
        print(pd.DataFrame(delay_rows))
    else:
        print("[WARN] No delays computed; no CSV written.")

if __name__ == "__main__":
    main()

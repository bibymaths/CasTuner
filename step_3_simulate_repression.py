#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step: ODE-based repression simulation (KD), delay scan, plots, and delays_repression.csv.

Deps (Py3.9):
  numpy<2.0, pandas, scipy
  FlowCytometryTools, flowio, flowutils
  plotnine
"""

import os
import re
import glob
import math
import warnings

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
OUT_PATH = "plots"
PARAM_PATH = "parameters"
os.makedirs(OUT_PATH, exist_ok=True)
os.makedirs(PARAM_PATH, exist_ok=True)

CH_FSC_A = "FSC-A"
CH_SSC_A = "SSC-A"
CH_FSC_H = "FSC-H"
CH_BFP   = "BV421-A"   # tagBFP
CH_mCh   = "PE-A"      # mCherry

BOUND_MIN = {CH_FSC_A: 0.4e5, CH_SSC_A: 0.20e5}
BOUND_MAX = {CH_FSC_A: 2.0e5, CH_SSC_A: 1.3e5}
SINGLET_RATIO_LOW, SINGLET_RATIO_HIGH = 0.85, 1.15

PLOT_W = 1.5 * 1.618
PLOT_H = 1.5
COL_BFP = "#4DBBD5FF"
COL_CH  = "#E64B35FF"
THEME = theme_classic()

# ----------------------------
# Utilities: gating & medians
# ----------------------------
def apply_boundary_gate(df):
    m = (
        (df[CH_FSC_A] >= BOUND_MIN[CH_FSC_A]) & (df[CH_FSC_A] <= BOUND_MAX[CH_FSC_A]) &
        (df[CH_SSC_A] >= BOUND_MIN[CH_SSC_A]) & (df[CH_SSC_A] <= BOUND_MAX[CH_SSC_A])
    )
    return df.loc[m]

def apply_singlet_gate(df):
    ratio = df[CH_FSC_H] / (df[CH_FSC_A].replace(0, np.nan))
    m = (ratio >= SINGLET_RATIO_LOW) & (ratio <= SINGLET_RATIO_HIGH)
    return df.loc[m]

def median_channels_for_file(fpath):
    data = FCMeasurement(ID=os.path.basename(fpath), datafile=fpath).data
    for ch in (CH_FSC_A, CH_SSC_A, CH_FSC_H, CH_BFP, CH_mCh):
        if ch not in data.columns:
            raise ValueError(f"Channel '{ch}' not in {fpath}")
    g = apply_boundary_gate(data)
    if len(g) == 0: g = data
    s = apply_singlet_gate(g)
    if len(s) == 0: s = g
    med = s.median(numeric_only=True)
    med["__filename"] = os.path.splitext(os.path.basename(fpath))[0]
    return med

def load_flowset_medians(folder):
    files = sorted(glob.glob(os.path.join(folder, "*.fcs")))
    if not files:
        raise FileNotFoundError(f"No FCS in {folder}")
    return pd.DataFrame([median_channels_for_file(f) for f in files])

# ----------------------------
# NFC background (mean of first 3 medians)
# ----------------------------
def compute_nfc_background(nfc_dir):
    df = load_flowset_medians(nfc_dir)
    mBFP_neg = float(df.iloc[:3][CH_BFP].mean())
    mmCherry_neg = float(df.iloc[:3][CH_mCh].mean())
    return mBFP_neg, mmCherry_neg

# ----------------------------
# Parse filename to tokens: ... _ plasmid _ exp _ rep _ time _ ...
# ----------------------------
_SPLIT = re.compile(r"_")
def parse_timecourse_name(name: str):
    parts = _SPLIT.split(name)
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
# Load time-course KD dataset with background subtraction
# ----------------------------
def load_kd(mBFP_neg, mmCherry_neg):
    med = load_flowset_medians("fcs_files/time-course_data").copy()
    med[CH_BFP] = med[CH_BFP] - mBFP_neg
    med[CH_mCh] = med[CH_mCh] - mmCherry_neg
    parsed = med["__filename"].apply(parse_timecourse_name).tolist()
    parsed_df = pd.DataFrame(parsed, columns=["plasmid", "exp", "rep", "time"])
    df = pd.concat([med[[CH_BFP, CH_mCh]], parsed_df], axis=1)
    df = df[df["exp"] == "KD"].copy()
    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df = df.dropna(subset=["time"])
    # plasmid order as in R
    cat = pd.CategoricalDtype(["SP430", "SP428", "SP430ABA", "SP427", "SP411"], ordered=True)
    df["plasmid"] = df["plasmid"].astype(cat)
    return df

# ----------------------------
# Compute mCherry fold-change & BFP min-max norm (KD)
# ----------------------------
def add_fc_and_norm(df: pd.DataFrame) -> pd.DataFrame:
    # mCherry fold-change vs mean at time 0 (per plasmid)
    t0 = (df.query("time == 0")
            .groupby("plasmid")[CH_mCh]
            .mean().rename("mmcherry").reset_index())
    out = df.merge(t0, on="plasmid", how="left")
    out["fc.cherry"] = out[CH_mCh] / out["mmcherry"]

    # BFP min-max: mean.init @ t==0; mean.final @ t>10
    mean_final = (out.query("time > 10").groupby("plasmid")[CH_BFP]
                  .mean().rename("mean.final").reset_index())
    mean_init  = (out.query("time == 0").groupby("plasmid")[CH_BFP]
                  .mean().rename("mean.init").reset_index())
    out = out.merge(mean_final, on="plasmid", how="left").merge(mean_init, on="plasmid", how="left")
    denom = (out["mean.final"] - out["mean.init"]).replace(0, np.nan)
    out["norm.bfp"] = (out[CH_BFP] - out["mean.init"]) / denom
    return out

# ----------------------------
# Read parameters and merge (t_up, Hill K,n, alpha)
# ----------------------------
def load_parameters():
    """
    Load/merge parameter CSVs. If alphamcherry.csv contains only a single
    global 'alpha' value (no 'plasmid' column), broadcast it to all plasmids.
    """
    import os
    import pandas as pd

    base = "parameters"

    def _read_norm(path):
        df = pd.read_csv(path)
        df.columns = (
            df.columns.str.strip()
                      .str.lower()
                      .str.replace(r"[\s\-]+", "_", regex=True)
        )
        return df

    # --- Load the three per-plasmid tables first ---
    t_down = _read_norm(os.path.join(base, "half_times_downregulation.csv"))
    t_up   = _read_norm(os.path.join(base, "half_times_upregulation.csv"))
    hills  = _read_norm(os.path.join(base, "Hill_parameters.csv"))

    # normalize expected columns
    if "halftime" in t_down.columns:
        t_down = t_down.rename(columns={"halftime": "t_down"})
    if "halftime" in t_up.columns:
        t_up = t_up.rename(columns={"halftime": "t_up"})
    if "k" not in hills.columns and "K" in hills.columns:
        hills = hills.rename(columns={"K": "k"})

    # ensure 'plasmid' column present in the three tables
    for name, df in [("half_times_downregulation.csv", t_down),
                     ("half_times_upregulation.csv", t_up),
                     ("Hill_parameters.csv", hills)]:
        if "plasmid" not in df.columns:
            raise KeyError(f"'{name}' must have a 'plasmid' column. Found: {list(df.columns)}")

    # merge those three to get the plasmid universe
    params = (t_down.merge(t_up, on="plasmid", how="outer")
                    .merge(hills, on="plasmid", how="outer"))

    # --- Load alpha, handle global or per-plasmid ---
    alpha_path = os.path.join(base, "alphamcherry.csv")
    alpha = _read_norm(alpha_path)

    if "alpha" not in alpha.columns:
        # tolerate variants
        for alt in ("alpha_mcherry", "alphamcherry", "a"):
            if alt in alpha.columns:
                alpha = alpha.rename(columns={alt: "alpha"})
                break

    if "alpha" not in alpha.columns:
        raise KeyError(f"'alphamcherry.csv' must contain an 'alpha' column. Found: {list(alpha.columns)}")

    if "plasmid" in alpha.columns:
        # per-plasmid alpha — merge directly
        params = params.merge(alpha[["plasmid", "alpha"]], on="plasmid", how="left")
    else:
        # global alpha — must be a single row; broadcast to all plasmids
        if len(alpha) != 1:
            raise KeyError(
                f"'alphamcherry.csv' has no 'plasmid' column and {len(alpha)} rows; "
                f"expected exactly 1 global alpha value."
            )
        global_alpha = float(alpha["alpha"].iloc[0])
        params["alpha"] = global_alpha

    # sanity checks
    required = ["plasmid", "t_up", "k", "n", "alpha"]
    missing = [c for c in required if c not in params.columns]
    if missing:
        raise KeyError(f"Merged parameters missing columns: {missing}. Columns: {list(params.columns)}")

    # warn on NA rows
    na_rows = params[required].isna().any(axis=1)
    if na_rows.any():
        print("[WARN] Some parameter rows contain NA values:\n", params.loc[na_rows, ["plasmid"] + required[1:]])

    return params


# ----------------------------
# ODE system and simulation
# dR = beta - R*(ln2/t_up)
# dY = (K^n/(K^n + R^n)) - alpha*Y
# initial: R(0)=0, Y(0)=1/alpha (so Y*alpha starts at 1)
# ----------------------------
def ode_rhs(t, y, beta, t_up, K, n, alpha):
    R, Y = y
    dR = beta - R * (math.log(2.0) / t_up)
    dY = (K**n) / (K**n + max(R, 0.0)**n) - alpha * Y
    return [dR, dY]

def simulate(beta, t_up, K, n, alpha, t_start, t_end, dt=0.01):
    t_eval = np.arange(t_start, t_end + dt/2, dt)
    y0 = [0.0, 1.0/alpha]
    sol = solve_ivp(
        ode_rhs, (t_start, t_end), y0,
        t_eval=t_eval, method="LSODA", max_step=0.5,
        args=(beta, t_up, K, n, alpha)
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    # scale Y by alpha like in R: out[,3] <- out[,3]*alpha
    R = sol.y[0]
    Y = sol.y[1] * alpha
    return sol.t, R, Y

# ----------------------------
# Delay scan: compute MAE between simulated Y and experimental mean fc.cherry
# for delays in [0,25] step 0.5
# ----------------------------
def delay_scan(pl_df, pars, t_grid=(0, 150, 0.01), delays=np.arange(0, 25.0 + 1e-9, 0.5)):
    t0, t1, dt = t_grid

    def pick(series, *keys):
        for k in keys:
            if k in series and pd.notna(series[k]):
                return float(series[k])
        raise KeyError(f"Missing parameter(s) {keys} in row: {series.to_dict()}")

    # tolerate K/k, n/N (and t_up/t.down case), alpha
    t_up  = pick(pars, "t_up", "T_up", "tup")
    K     = pick(pars, "K", "k")
    n     = pick(pars, "n", "N", "hill_n")
    alpha = pick(pars, "alpha", "Alpha", "alphamcherry")

    beta = math.log(2.0) / t_up

    # experimental mean fc per time
    exp_mean = (pl_df.groupby("time")["fc.cherry"]
                .mean().reset_index().sort_values("time"))
    exp_t = exp_mean["time"].to_numpy(float)
    exp_y = exp_mean["fc.cherry"].to_numpy(float)

    # simulate base (no delay)
    base_t, base_R, base_Y = simulate(beta, t_up, K, n, alpha, t0, t1, dt)

    rows = []
    for d in delays:
        t_del = base_t + d
        fY = interp1d(t_del, base_Y, kind="linear", bounds_error=False,
                      fill_value=(base_Y[0], base_Y[-1]))
        y_hat = fY(exp_t)
        resid = exp_y - y_hat
        N = np.isfinite(resid).sum()
        mae = np.inf if N <= 1 else np.nansum(np.abs(resid)) / (N - 1)
        rows.append((float(d), int(N), float(mae)))

    out = pd.DataFrame(rows, columns=["t", "N", "MAE"])
    best_row = out.loc[out["MAE"].idxmin()]
    best_d, best_mae = float(best_row["t"]), float(best_row["MAE"])
    return out, best_d, best_mae, (base_t, base_R, base_Y)


# ----------------------------
# Plot helpers
# ----------------------------
def save_mae_plot(df_mae, out_pdf, y_max):
    best_row = df_mae.loc[df_mae["MAE"].idxmin()]
    p = (
        ggplot(df_mae, aes("t", "MAE"))
        + geom_point(size=0.1, alpha=0.4, color="black")
        + geom_point(data=best_row.to_frame().T, mapping=aes("t", "MAE"),
                     size=0.8, alpha=1.0, color=COL_CH)
        + geom_line(data=df_mae.assign(minMAE=df_mae["MAE"].min()), mapping=aes("t", "minMAE"), linetype="dashed")
        + coord_cartesian(xlim=(0, 25), ylim=(0, y_max))
        + labs(y="MAE", x="Delay (hours)")
        + THEME
    )
    p.save(os.path.join(OUT_PATH, out_pdf), width=PLOT_W, height=PLOT_H, units="in")

def save_fit_plots(pl_df, base_t, base_R, base_Y, delay, tag_prefix):
    # mCherry
    p1 = (
        ggplot(pl_df, aes("time", "fc.cherry"))
        + geom_point(size=0.8, alpha=0.4, color=COL_CH)
        + coord_cartesian(xlim=(0, 150), ylim=(0, 1.3))
        + labs(x="Time (hours)", y="mCherry")
        + THEME
        + geom_line(pd.DataFrame({"time": base_t, "Y": base_Y}), aes("time", "Y"))
        + geom_line(pd.DataFrame({"time": base_t + delay, "Y": base_Y}), aes("time", "Y"), linetype="dashed")
    )
    p1.save(os.path.join(OUT_PATH, f"KD_ODE_mCherry_{tag_prefix}.pdf"), width=PLOT_W, height=PLOT_H, units="in")

    # tagBFP
    p2 = (
        ggplot(pl_df, aes("time", "norm.bfp"))
        + geom_point(size=0.8, alpha=0.4, color=COL_BFP)
        + coord_cartesian(xlim=(0, 150), ylim=(0, 1.3))
        + labs(x="Time (hours)", y="tagBFP")
        + THEME
        + geom_line(pd.DataFrame({"time": base_t, "R": base_R}), aes("time", "R"))
    )
    p2.save(os.path.join(OUT_PATH, f"KD_ODE_tagBFP_{tag_prefix}.pdf"), width=PLOT_W, height=PLOT_H, units="in")

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

        # record chosen delay (R script comments: SP430A~6h, SP428~3h, others 0h typically)
        delay_rows.append({"plasmid": pl_label, "d_rev": best_delay})

        # save MAE plot with reasonable y-limits tuned per case
        ycap = 0.4 if pl_label in ("SP427",) else 0.3
        save_mae_plot(mae_df, f"MAE_KD_{tag}_mcherry.pdf", ycap)

        # save fit plots (mCherry & BFP)
        save_fit_plots(pl_df, bt, bR, bY, best_delay, tag_prefix=tag.replace("-", "_"))

    # write delays
    if delay_rows:
        out = pd.DataFrame(delay_rows, columns=["plasmid", "d_rev"])
        out.to_csv(os.path.join(PARAM_PATH, "delays_repression.csv"), index=False)
        print("Saved:", os.path.join(PARAM_PATH, "delays_repression.csv"))
        print(out)
    else:
        print("[WARN] No delays estimated.")

if __name__ == "__main__":
    main()

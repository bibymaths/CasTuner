#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
KD (up-regulation) time-course: fit exponential half-times on min-max scaled BFP.

Requires (Py3.9-compatible):
  - numpy<2.0, pandas, scipy
  - FlowCytometryTools, flowio, flowutils
  - plotnine

Folder structure:
  fcs_files/
    NFC/*.fcs
    time-course_data/*.fcs

Outputs:
  plots/KD_KRAB-Split-dCas9_fitting.pdf
  plots/KD_dCas9_fitting.pdf
  plots/KD_KRAB-dCas9_fitting.pdf
  plots/KD_HDAC4-dCas9_fitting.pdf
  plots/KD_CasRx_fitting.pdf
  parameters/half_times_upregulation.csv
"""

import os
import re
import glob
import math
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, differential_evolution
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

# Plot sizing (inches)
PLOT_W = 1.5 * 1.618
PLOT_H = 1.5

# Aesthetics
POINT_COLOR = "#4DBBD5FF"
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

def median_channels_for_file(fpath: str, mBFP_neg: float = 0.0, mmCherry_neg: float = 0.0) -> pd.Series:
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

    # Apply subtraction BEFORE median, if values are non-zero
    if mBFP_neg != 0.0 or mmCherry_neg != 0.0:
        singlets_sub = singlets.copy()
        singlets_sub[CH_BFP] = singlets_sub[CH_BFP] - mBFP_neg
        singlets_sub[CH_mCh] = singlets_sub[CH_mCh] - mmCherry_neg
        s = singlets_sub.median(numeric_only=True)
    else:
        # Original behavior for loading NFC files
        s = singlets.median(numeric_only=True)

    s["__filename"] = os.path.splitext(os.path.basename(fpath))[0]
    return s

def load_flowset_medians(folder: str, mBFP_neg: float = 0.0, mmCherry_neg: float = 0.0) -> pd.DataFrame:
    files = sorted(glob.glob(os.path.join(folder, "*.fcs")))
    if not files:
        raise FileNotFoundError(f"No FCS files found in: {folder}")
    # Pass background values (or defaults of 0) to the helper
    return pd.DataFrame([median_channels_for_file(f, mBFP_neg, mmCherry_neg) for f in files])

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
# R: separate(rowname, c(NA,NA,"plasmid","exp","rep","time",NA,NA), "_")
# -> plasmid=parts[2], exp=parts[3], rep=parts[4], time=parts[5]
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
        m = re.search(r"(\d+(?:\.\d+)?)", time_s)
        time = float(m.group(1)) if m else np.nan
    return plasmid, exp, rep, time

# ----------------------------
# Build KD dataset
# ----------------------------
def load_kd_timecourse(mBFP_neg: float, mmCherry_neg: float) -> pd.DataFrame:
    med = load_flowset_medians("fcs_files/time-course_data").copy()

    med = load_flowset_medians(
        "fcs_files/time-course_data",
        mBFP_neg,
        mmCherry_neg
    ).copy()

    parsed = med["__filename"].apply(parse_timecourse_name).tolist()
    parsed_df = pd.DataFrame(parsed, columns=["plasmid", "exp", "rep", "time"])
    df = pd.concat([med[[CH_BFP, CH_mCh]], parsed_df], axis=1)

    # keep KD experiments
    df = df[df["exp"] == "KD"].copy()
    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df = df.dropna(subset=["time"])
    # plasmid order (factor levels in R)
    cat = pd.CategoricalDtype(["SP430", "SP428", "SP430ABA", "SP427", "SP411"], ordered=True)
    df["plasmid"] = df["plasmid"].astype(cat)
    return df

# ----------------------------
# Normalization (KD):
# mean.final = mean(BFP | time>10)
# mean.init  = mean(BFP | time==0)
# norm.bfp = (BFP - mean.init) / (mean.final - mean.init)
# ----------------------------
def add_minmax_norm_kd(df: pd.DataFrame, final_thr=10.0, use_median=True) -> pd.DataFrame:
    agg = "median" if use_median else "mean"
    grp = df.query(f"time >= {final_thr}").groupby("plasmid")[CH_BFP]
    final_stat = getattr(grp, agg)().rename("final").reset_index()

    init_stat = (
        df.query("time == 0")
          .groupby("plasmid")[CH_BFP]
    )
    init_stat = getattr(init_stat, agg)().rename("init").reset_index()

    out = df.merge(final_stat, on="plasmid", how="left").merge(init_stat, on="plasmid", how="left")
    denom = (out["final"] - out["init"])
    out["norm.bfp"] = (out[CH_BFP] - out["init"]) / denom.replace(0, np.nan)

    # Optional: clamp tiny negatives/overshoots like many R workflows
    out["norm.bfp"] = out["norm.bfp"].clip(lower=0.0, upper=1.2)
    return out

# ----------------------------
# Model & fitting (KD):
# y = 1 - exp(-t * ln(2) / t_half)
# ----------------------------
def exp_rise_to_one(t, tau):
    t = np.asarray(t, float)
    return 1.0 - np.exp(-t / float(tau))

from scipy.optimize import least_squares

def fit_tau_global_then_local(t, y, sigma=None, bounds=(1e-6, 1e3)):
    t = np.asarray(t, float); y = np.asarray(y, float)
    def model(t, tau): return 1.0 - np.exp(-t / tau)
    def sse(x):
        tau = x[0]
        return np.sum((y - model(t, tau))**2)

    res_g = differential_evolution(sse, bounds=[bounds], polish=False, seed=0)
    # polish with least_squares
    return fit_half_time_rise(t, y, sigma=sigma, p0_tau=float(res_g.x[0]),
                                 bounds=bounds, loss="soft_l1", f_scale=0.1)

def fit_half_time_rise(t, y, sigma=None, p0_tau=1.0,
                          bounds=(1e-6, 1e3),
                          loss="soft_l1", f_scale=0.1, max_nfev=100000):
    t = np.asarray(t, float); y = np.asarray(y, float)
    # weights
    if sigma is None:
        w = 1.0
    else:
        s = np.asarray(sigma, float)
        if s.ndim == 0: s = np.full_like(y, float(s))
        w = 1.0 / np.clip(s, 1e-12, np.inf)

    def model(t, tau): return 1.0 - np.exp(-t / tau)
    def residuals(x):
        tau = x[0]
        return (y - model(t, tau)) * w

    res = least_squares(
        residuals, x0=[p0_tau], bounds=bounds,
        loss=loss, f_scale=f_scale, max_nfev=max_nfev, method="trf"
    )
    tau = float(res.x[0])
    # covariance approx: (J^T J)^-1 * s^2 (Gauss–Newton)
    if res.jac is not None and res.jac.size:
        J = res.jac
        # robust loss rescales residuals; use plain RSS for s^2
        rss = float(np.sum((y - model(t, tau))**2))
        dof = max(len(y) - 1, 1)
        s2 = rss / dof
        JTJ = J.T @ J
        try:
            pcov = np.linalg.inv(JTJ) * s2
            se_tau = float(np.sqrt(pcov[0,0]))
        except np.linalg.LinAlgError:
            pcov = None; se_tau = np.nan
    else:
        pcov = None; se_tau = np.nan

    t_half = tau * np.log(2.0)
    se_th  = se_tau * np.log(2.0) if np.isfinite(se_tau) else np.nan
    yhat   = model(t, tau)
    rss    = float(np.sum((y - yhat)**2))
    return t, y, t_half, se_th, pcov, yhat, rss, len(y)

def _clean_xy(t, y):
    t = np.asarray(t, float); y = np.asarray(y, float)
    m = np.isfinite(t) & np.isfinite(y)
    t, y = t[m], y[m]
    if t.size < 3:
        raise ValueError(f"Not enough valid points after cleaning: {t.size}")
    return t, y
# ----------------------------
# Sigma (weighting) schemes
# ----------------------------
def make_sigma_scheme(dfp, scheme: str):
    """
    Returns sigma vector aligned to dfp rows (or None).
    Schemes:
      - 'none'        : unweighted
      - 'prop_y'      : sigma ∝ |y| (heteroscedastic)
      - 'const_frac'  : sigma = c * (max(y)-min(y)) with c=0.05
      - 'sem_by_time' : per-timepoint SEM across replicates (fallback to prop_y if single)
    """
    y = np.asarray(dfp["norm.bfp"], float)
    if scheme == "none":
        return None

    if scheme == "prop_y":
        s = np.clip(np.abs(y), 1e-3, None)
        return s

    if scheme == "const_frac":
        rng = float(np.nanmax(y) - np.nanmin(y)) if y.size else 1.0
        s = np.full_like(y, max(1e-3, 0.05 * rng))
        return s

    if scheme == "sem_by_time":
        # group by time (replicates at same time)
        g = dfp.groupby("time")["norm.bfp"]
        sem_by_t = g.transform(lambda v: np.std(v, ddof=1) / np.sqrt(max(1, len(v))))
        s = sem_by_t.to_numpy(dtype=float)
        # if all SEM are zero/NaN (single replicate), fallback to prop_y
        if not np.isfinite(s).any() or np.nanmax(s) == 0.0:
            s = np.clip(np.abs(y), 1e-3, None)
        s = np.nan_to_num(s, nan=np.nanmedian(s) if np.isfinite(np.nanmedian(s)) else 1.0)
        s[s <= 0] = np.nanmedian(s[s > 0]) if np.any(s > 0) else 1.0
        return s

    raise ValueError(f"Unknown sigma scheme: {scheme}")

# ----------------------------
# Information criteria
# ----------------------------
def aic_bic_aicc(rss, n, k):
    # Gaussian LS: AIC = n ln(RSS/n) + 2k ; BIC = n ln(RSS/n) + k ln n
    if n <= k + 1:
        n_eff = k + 2  # guard
    else:
        n_eff = n
    aic = n_eff * math.log(max(rss / n_eff, 1e-30)) + 2 * k
    bic = n_eff * math.log(max(rss / n_eff, 1e-30)) + k * math.log(n_eff)
    aicc = aic + (2 * k * (k + 1)) / max(n_eff - k - 1, 1e-6)
    return aic, bic, aicc

# ----------------------------
# Residual bootstrap
# ----------------------------
def bootstrap_half_time(t, y, t_half_hat, B=1000, rng=None, sigma=None, **fit_kwargs):
    """
    Residual bootstrap:
      y* = yhat + r_resampled
      refit to get t_half*
    Returns: (t_half_mean, t_half_std, ci_low, ci_high, samples)
    """
    rng = np.random.default_rng(rng)
    yhat = exp_rise_to_one(t, t_half_hat)
    resid = (y - yhat)
    if not np.isfinite(resid).all():
        resid = np.nan_to_num(resid, nan=0.0)

    samples = []
    for _ in range(int(B)):
        r_star = rng.choice(resid, size=resid.size, replace=True)
        y_star = yhat + r_star
        try:
            _, _, t_star, _, _, _, _, _ = fit_half_time_rise(
                t, y_star, sigma=sigma, **fit_kwargs
            )
            samples.append(t_star)
        except Exception:
            continue

    if len(samples) == 0:
        return np.nan, np.nan, np.nan, np.nan, []

    s_arr = np.array(samples, float)
    mean = float(np.mean(s_arr))
    std = float(np.std(s_arr, ddof=1)) if s_arr.size > 1 else 0.0
    ci_low, ci_high = np.percentile(s_arr, [2.5, 97.5]).tolist()
    return mean, std, float(ci_low), float(ci_high), samples

# ----------------------------
# Evaluate multiple sigma schemes and pick best by AICc
# ----------------------------
def evaluate_schemes(dfp, schemes=("none", "sem_by_time", "prop_y", "const_frac"),
                     p0=0.8, bootstrap_B=500, rng=0, **fit_kwargs):
    results = []
    best = None
    best_aicc = np.inf

    for scheme in schemes:
        sigma_vec = make_sigma_scheme(dfp, scheme)
        try:
            t, y, t_half, se, pcov, yhat, rss, n = fit_tau_global_then_local(
                dfp["time"], dfp["norm.bfp"], sigma=sigma_vec
            )
        except Exception as e:
            results.append({"scheme": scheme, "ok": False, "error": str(e)})
            continue

        k = 1  # one parameter: t_half
        aic, bic, aicc = aic_bic_aicc(rss, n, k)

        # Bootstrap (use same sigma weighting)
        bt_mean, bt_std, ci_lo, ci_hi, _ = bootstrap_half_time(
            t, y, t_half, B=bootstrap_B, rng=rng, sigma=sigma_vec, **fit_kwargs
        )

        row = {
            "scheme": scheme, "ok": True,
            "t_half": t_half, "se_cov": se,
            "rss": rss, "n": n, "aic": aic, "bic": bic, "aicc": aicc,
            "bt_mean": bt_mean, "bt_std": bt_std,
            "ci_low": ci_lo, "ci_high": ci_hi
        }
        results.append(row)

        if aicc < best_aicc:
            best_aicc = aicc
            best = row

    if best is None:
        raise RuntimeError("All sigma schemes failed.")
    return best, pd.DataFrame(results)

# ----------------------------
# Plot helper (match R look)
# ----------------------------
def save_kd_plot(df, t_half, out_pdf):
    tmin, tmax = float(df["time"].min()), float(df["time"].max())
    curve = pd.DataFrame({"time": np.linspace(tmin, tmax, 200)})
    curve["fit"] = exp_rise_to_one(curve["time"], t_half)

    p = (
        ggplot(df, aes("time", "norm.bfp"))
        + geom_point(size=0.4, alpha=0.7, color=POINT_COLOR)
        + geom_line(data=curve, mapping=aes(x="time", y="fit"), color="black")
        + coord_cartesian(ylim=(-0.15, 1.2))
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

    # 2) Load KD time-course with background-subtracted medians
    kd = load_kd_timecourse(mBFP_neg, mmCherry_neg)

    # 3) Min–max normalize per plasmid (KD definition)
    kd = add_minmax_norm_kd(kd, final_thr=10.0, use_median=True)

    # DEBUG: show counts
    print("KD rows:", len(kd))
    print(kd.groupby(["plasmid"])["time"].agg(["count","min","max"]).reset_index())

    # 4) Fit per plasmid using multiple sigma schemes, bootstrap, and save plots
    targets = [
        ("SP430ABA", "KD_KRAB-Split-dCas9_fitting.pdf", "SP430A"),
        ("SP430",    "KD_dCas9_fitting.pdf",            "SP430"),
        ("SP428",    "KD_KRAB-dCas9_fitting.pdf",       "SP428"),
        ("SP427",    "KD_HDAC4-dCas9_fitting.pdf",      "SP427"),
        ("SP411",    "KD_CasRx_fitting.pdf",            "SP411"),
    ]

    rng_seed = 0
    rows, diag_rows = [], []
    for plasmid, pdfname, label in targets:
        dfp = kd[kd["plasmid"] == plasmid].copy()
        if dfp.empty:
            print(f"[WARN] No KD data for {plasmid}; skipping.")
            continue

        try:
            best, diag = evaluate_schemes(
                dfp,
                schemes=("none", "sem_by_time", "prop_y", "const_frac"),
                p0=0.8,
                bootstrap_B=500,   # <- quick check; raise to 500+ later
                rng=2,
            )
        except Exception as e:
            print(f"[ERROR] {label}: {e}")
            continue

        t_half = float(best["t_half"])
        se = float(best["bt_std"]) if np.isfinite(best["bt_std"]) and best["bt_std"] > 0 else float(best["se_cov"])
        ci_low = best.get("ci_low", np.nan)
        ci_high = best.get("ci_high", np.nan)
        scheme_used = best["scheme"]

        print(f"{label}: t1/2 = {t_half:.6g}  (SE≈{se:.3g})  [{ci_low:.3g}, {ci_high:.3g}] via {scheme_used}")

        rows.append({
            "plasmid": label,
            "halftime": t_half,
            "se": se,
            "ci_low": ci_low,
            "ci_high": ci_high,
            "scheme": scheme_used
        })

        try:
            save_kd_plot(dfp, t_half, pdfname)
        except Exception as e:
            print(f"[WARN] Plot failed for {label}: {e}")

        diag = diag.copy()
        diag.insert(0, "plasmid", label)
        diag_rows.append(diag)

    # 5) Save outputs
    if rows:
        ht = pd.DataFrame(rows, columns=["plasmid", "halftime", "se", "ci_low", "ci_high", "scheme"])
        ht.to_csv(os.path.join(PARAM_PATH, "half_times_upregulation.csv"), index=False)
        print("\nSaved:", os.path.join(PARAM_PATH, "half_times_upregulation.csv"))
        print(ht)

        if diag_rows:
            diag_all = pd.concat(diag_rows, ignore_index=True)
            diag_all.to_csv(os.path.join(PARAM_PATH, "half_times_upregulation_diagnostics.csv"), index=False)
            print("Saved:", os.path.join(PARAM_PATH, "half_times_upregulation_diagnostics.csv"))
    else:
        print("[WARN] No half-times estimated.")

if __name__ == "__main__":
    main()

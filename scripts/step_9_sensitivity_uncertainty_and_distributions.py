#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_9_sensitivity_uncertainty_and_distributions.py
---------------------------------------------------
Rigorous sensitivity + uncertainty quantification for CasTuner-Python Port pipeline.

Purpose:
    - For each plasmid/model (REV and KD), perform Sobol sensitivity analysis
        on key metrics (dynamic range, t50, t10-90, overshoot, AUC, end level).
    - For each plasmid/model, perform bootstrap parameter sampling to
        generate predictive bands and parameter/metric distributions.

Requires:
    - step_2_simulate_derepression.py
    - step_3_simulate_repression.py

Inputs:
  parameters/half_times_upregulation.csv        (plasmid, halftime, se)  -> t_up
  parameters/half_times_downregulation.csv      (plasmid, halftime, se)  -> t_down
  parameters/Hill_parameters.csv                (plasmid, K, n)
  parameters/alphamcherry.csv                   (alpha) or (plasmid, alpha)
  parameters/delays_derepression.csv            (plasmid, d_rev)
  parameters/delays_repression.csv              (plasmid, d_rev)

Notes:
  - Half-times have SE from Step 1a/1b, so bootstrap uses Normal(mean, se).
    (Clipped to positive.)
  - Hill K and n have no SE in your pipeline; bootstrap uses a conservative
    lognormal perturbation around fitted values (default CV=0.20).
  - Delay bootstrap: uses Normal(mean=best_delay, sd=0.5h) clipped [0, 25].
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse, os, contextlib, math, sys
from SALib.sample import saltelli
from SALib.analyze import sobol
from joblib import Parallel, delayed
import step_2_simulate_derepression as step2
import step_3_simulate_repression as step3


def parse_args():
    """
    Parse command-line arguments for directories and parameters.

    Returns:
        argparse.Namespace: Parsed arguments with attributes for each parameter.

    """
    ap = argparse.ArgumentParser()
    ap.add_argument("--params-dir", default="parameters", help="Parameters output dir (from config.yaml)")
    ap.add_argument("--plots-dir", default="plots", help="Plots output dir (from config.yaml)")
    ap.add_argument("--nfc-dir", default="fcs_files/NFC", help="NFC directory (from config.yaml)")
    ap.add_argument("--timecourse-dir", default="fcs_files/time-course_data", help="Timecourse directory (from config.yaml)")
    ap.add_argument("--sobol-base", type=int, default=512)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=13)
    return ap.parse_args()

@contextlib.contextmanager
def quiet():
    """
    Suppress stdout/stderr (Step 2 prints per ODE call).

    Usage:
        with quiet():
            # code that prints
            ...
    """
    with open(os.devnull, "w") as devnull:
        with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
            yield

# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------
SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))
ARGS = parse_args()

PARAM_PATH = Path(ARGS.params_dir)
PLOTS_PATH = Path(ARGS.plots_dir)

OUT_PARAM = PARAM_PATH / "uq_sensitivity"
OUT_PLOTS = PLOTS_PATH / "uq_sensitivity"
OUT_PARAM.mkdir(parents=True, exist_ok=True)
OUT_PLOTS.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------
def _clip_pos(x, eps=1e-9):
    """
    Clip input to be at least eps (default 1e-9).

    Args:
        x: array-like input
        eps: minimum value

    Returns:
        np.ndarray: clipped array
    """
    return np.maximum(np.asarray(x, float), eps)


def _safe_float(x, default=np.nan):
    """
    Safely convert x to float, return default on failure.

    Args:
        x: input value
        default: value to return on failure

    Returns:
        float or default
    """
    try:
        return float(x)
    except Exception:
        return default


def _norm_cols(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize DataFrame column names: strip, lower, replace spaces/hyphens with underscores.

    Args:
        df (pd.DataFrame): input DataFrame

    Returns:
        pd.DataFrame: DataFrame with normalized column names
    """
    df = df.copy()
    df.columns = (df.columns.astype(str).str.strip().str.lower()
                  .str.replace(r"[\s\-]+", "_", regex=True))
    return df

def _normalize_plasmid_labels(s: pd.Series) -> pd.Series:
    """
    Normalize plasmid labels to standard SP* format.

    Args:
        s (pd.Series): Series of plasmid labels

    Returns:
        pd.Series: Series with normalized plasmid labels
    """
    s = s.astype(str)
    # If already SP*, keep it.
    m = s.str.startswith("SP")
    out = s.copy()

    # Only map “raw” labels when not already normalized
    raw_map = {
        "430": "SP430",
        "411": "SP411",
        "427": "SP427",
        "428": "SP428",
        "430ABA": "SP430A",
        "SP430ABA": "SP430A",  # common alias in your REV data
    }
    out.loc[~m] = out.loc[~m].replace(raw_map)
    # Also handle SP430ABA -> SP430A even if it starts with SP
    out = out.replace({"SP430ABA": "SP430A"})
    return out

def saveplot(path: Path):
    """
    Save current matplotlib figure to path with tight layout and high DPI.

    Args:
        path (Path): Output file path

    Returns:
        None
    """
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def metric_summary(y: np.ndarray, t: np.ndarray):
    """
    Compute metrics used for sensitivity/UQ.

    Metrics:
      dynamic_range = max(y) - min(y)
      t50           = time to reach 50% of total change
      t10_90        = t90 - t10
      overshoot     = max(y) - mean(last 20% of y)
      auc           = integral y(t) dt
      y_end         = y(t_end)


    Args:
        y (np.ndarray): output trajectory values
        t (np.ndarray): time points

    Returns:
        dict: computed metrics
    """
    y = np.asarray(y, float)
    t = np.asarray(t, float)
    if y.size < 5 or t.size != y.size or not np.all(np.isfinite(t)):
        return {k: np.nan for k in ["dynamic_range", "t50", "t10_90", "overshoot", "auc", "y_end"]}

    y0 = y[0]
    yend = y[-1]
    dyn = float(np.nanmax(y) - np.nanmin(y))

    total = yend - y0
    if not np.isfinite(total) or abs(total) < 1e-12:
        t10 = t50 = t90 = np.nan
    else:
        y10 = y0 + 0.10 * total
        y50 = y0 + 0.50 * total
        y90 = y0 + 0.90 * total

        def _t_at(level):
            """
            Find time when y crosses level (linear interpolation).

            Args:
                level (float): target y level

            Returns:
                float: interpolated time of crossing, or np.nan if not found
            """
            s = np.sign(y - level)
            # crossing index (first sign change)
            idx = np.where(s[:-1] * s[1:] <= 0)[0]
            if idx.size == 0:
                return np.nan
            i = int(idx[0])
            # linear interp
            t1, t2 = t[i], t[i + 1]
            y1, y2 = y[i], y[i + 1]
            if not np.isfinite([t1, t2, y1, y2]).all() or y2 == y1:
                return float(t2)
            return float(t1 + (level - y1) * (t2 - t1) / (y2 - y1))

        t10 = _t_at(y10)
        t50 = _t_at(y50)
        t90 = _t_at(y90)

    t10_90 = float(t90 - t10) if (np.isfinite(t10) and np.isfinite(t90)) else np.nan

    tail = y[int(0.8 * len(y)):] if len(y) >= 10 else y
    overshoot = float(np.nanmax(y) - np.nanmean(tail))

    auc = float(np.trapz(y, t))
    return {
        "dynamic_range": dyn,
        "t50": float(t50) if np.isfinite(t50) else np.nan,
        "t10_90": float(t10_90) if np.isfinite(t10_90) else np.nan,
        "overshoot": overshoot,
        "auc": auc,
        "y_end": float(yend),
    }


# ---------------------------------------------------------------------
# Load fitted parameter tables
# ---------------------------------------------------------------------
def load_fitted_tables():
    up = _norm_cols(pd.read_csv(PARAM_PATH / "half_times_upregulation.csv"))  # plasmid, halftime, se
    down = _norm_cols(pd.read_csv(PARAM_PATH / "half_times_downregulation.csv"))  # plasmid, halftime, se
    hill = _norm_cols(pd.read_csv(PARAM_PATH / "Hill_parameters.csv"))  # plasmid, k/K, n
    alpha = _norm_cols(pd.read_csv(PARAM_PATH / "alphamcherry.csv"))  # alpha (global or per plasmid)

    delays_rev = _norm_cols(pd.read_csv(PARAM_PATH / "delays_derepression.csv"))  # plasmid, d_rev
    delays_kd = _norm_cols(pd.read_csv(PARAM_PATH / "delays_repression.csv"))  # plasmid, d_rev

    # Ensure columns exist
    for df, need, name in [
        (up, ["plasmid", "halftime"], "half_times_upregulation.csv"),
        (down, ["plasmid", "halftime"], "half_times_downregulation.csv"),
        (hill, ["plasmid", "n"], "Hill_parameters.csv"),
        (delays_rev, ["plasmid", "d_rev"], "delays_derepression.csv"),
        (delays_kd, ["plasmid", "d_rev"], "delays_repression.csv"),
    ]:
        missing = [c for c in need if c not in df.columns]
        if missing:
            raise KeyError(f"{name} missing columns: {missing}. Found={list(df.columns)}")

    # Map Hill midpoint column to "K" (it can be K or k)
    # after reading hill:
    if "k" in hill.columns:
        hill = hill.rename(columns={"k": "K"})
    elif "K" in hill.columns:  # optional, if you ever skip _norm_cols
        hill = hill.rename(columns={"K": "K"})
    else:
        raise KeyError(f"Hill_parameters.csv must contain K or k. Found={list(hill.columns)}")

    # Alpha: global or per plasmid
    if "alpha" not in alpha.columns:
        raise KeyError(f"alphamcherry.csv missing 'alpha'. Found={list(alpha.columns)}")
    if "plasmid" not in alpha.columns:
        if len(alpha) != 1:
            raise KeyError("alphamcherry.csv is global but has != 1 row.")
        alpha_mode = "global"
        alpha_val = float(alpha["alpha"].iloc[0])
    else:
        alpha_mode = "per_plasmid"
        alpha_val = None

    return up, down, hill, alpha, alpha_mode, alpha_val, delays_rev, delays_kd


# ---------------------------------------------------------------------
# Build datasets
# ---------------------------------------------------------------------
def build_rev_dataset():
    """
    Build REV dataset from NFC and timecourse data using step2 functions.

    Returns:
        pd.DataFrame: REV timecourse dataset with background correction and normalization.
    """
    mBFP_neg, mmCherry_neg = step2.compute_nfc_background(ARGS.nfc_dir)
    df_tc = step2.load_timecourse_with_bg(mBFP_neg, mmCherry_neg)
    rev = step2.compute_rev_transforms(df_tc)
    return rev


def build_kd_dataset():
    """
    Build KD dataset from NFC and timecourse data using step3 functions.

    Returns:
        pd.DataFrame: KD timecourse dataset with background correction and normalization.
    """
    return step3.build_kd_dataset(nfc_dir=ARGS.nfc_dir, tc_dir=ARGS.timecourse_dir)


# ---------------------------------------------------------------------
# Parameter ranges for SALib (using fitted value ± band)
# ---------------------------------------------------------------------
def make_bounds(center, se=None, rel=0.25, lo=None, hi=None, clip=(1e-9, np.inf)):
    """
    Generate parameter bounds for SALib based on fitted center and SE or relative range.

    Args:
        center (float): Fitted center value.
        se (float, optional): Standard error of the fitted value. If provided and positive,
                                bounds are set to center ± 2*se.
        rel (float, optional): Relative range for bounds if se is not provided.
        lo (float, optional): Minimum bound override.
        hi (float, optional): Maximum bound override.
        clip (tuple, optional): Overall clipping range for bounds.

    Returns:
        tuple: (min_bound, max_bound) for the parameter.
    """
    c = float(center)
    if not np.isfinite(c) or c <= 0:
        return (clip[0], max(1.0, clip[0]))

    if se is not None and np.isfinite(se) and se > 0:
        a = c - 2.0 * float(se)
        b = c + 2.0 * float(se)
    else:
        a = c * (1.0 - rel)
        b = c * (1.0 + rel)

    if lo is not None:
        a = max(a, float(lo))
    if hi is not None:
        b = min(b, float(hi))

    a = max(a, clip[0])
    b = max(b, a * 1.0001)
    if np.isfinite(clip[1]):
        b = min(b, clip[1])
    return (float(a), float(b))


# ---------------------------------------------------------------------
# Bootstrap sampling for parameter distributions
# ---------------------------------------------------------------------
def bootstrap_params(rng, center, se=None, kind="normal", cv=0.20, clip_lo=1e-9, clip_hi=np.inf):
    """
    Draw bootstrap samples for a positive parameter.

    Args:
        rng: np.random.Generator instance for reproducibility.
        center (float): Fitted center value.
        se (float, optional): Standard error for normal sampling.
        kind (str): "normal" or "lognormal" sampling.
        cv (float): Coefficient of variation for lognormal sampling.
        clip_lo (float): Minimum value to clip the sample.
        clip_hi (float): Maximum value to clip the sample.

    Returns:
        float: Bootstrap sampled parameter value.
    """
    c = float(center)
    if kind == "normal" and se is not None and np.isfinite(se) and se > 0:
        x = rng.normal(loc=c, scale=float(se))
    elif kind == "lognormal":
        # lognormal parameterization with mean ~ c and given cv
        # sigma^2 = ln(1+cv^2); mu = ln(c) - 0.5*sigma^2
        sigma2 = math.log(1.0 + float(cv) ** 2)
        sigma = math.sqrt(sigma2)
        mu = math.log(max(c, clip_lo)) - 0.5 * sigma2
        x = rng.lognormal(mean=mu, sigma=sigma)
    else:
        # fallback: small normal noise
        x = rng.normal(loc=c, scale=0.05 * abs(c))

    x = float(np.clip(x, clip_lo, clip_hi))
    return x


# ---------------------------------------------------------------------
# Core runner
# ---------------------------------------------------------------------
def run_one_model(
        model_name: str,
        plasmid: str,
        data_df: pd.DataFrame,
        fitted: dict,
        bounds: dict,
        sim_func,
        obs_y_col="fc.cherry",
        time_col="time",
        n_sobol_base=ARGS.sobol_base,
        n_boot=ARGS.n_boot,
        seed=ARGS.seed
):
    """
    Run Sobol sensitivity analysis and bootstrap UQ for one model/plasmid.

    Args:
        model_name (str): Model name ("REV" or "KD").
        plasmid (str): Plasmid identifier.
        data_df (pd.DataFrame): Observed timecourse data for the plasmid.
        fitted (dict): Fitted parameter center values and SEs.
        bounds (dict): Parameter bounds for SALib.
        sim_func (callable): Function sim(pars_draw)->df(time,R,Y) to simulate the model.
        obs_y_col (str): Column name for observed output in data_df.
        time_col (str): Column name for time in data_df.
        n_sobol_base (int): Base sample size for Sobol analysis.
        n_boot (int): Number of bootstrap samples.
        seed (int): Random seed for reproducibility.

    Returns:
        sobol_df (pd.DataFrame): Sobol sensitivity indices per metric.
        boot_df (pd.DataFrame): Bootstrap sampled parameters and metrics.
    """
    rng = np.random.default_rng(seed)

    # ----------------------------------------------------------
    # Observed mean trajectory (for predictive band overlay)
    # ----------------------------------------------------------
    mean_obs = (data_df.groupby(time_col, as_index=False)[obs_y_col].mean()
                .sort_values(time_col))
    t_obs = mean_obs[time_col].to_numpy(float)
    y_obs = mean_obs[obs_y_col].to_numpy(float)

    # Standard time grid for predictive bands
    t_grid = np.arange(0.0, 150.0 + 0.05 / 2, 0.05)

    # ----------------------------------------------------------
    # A) SALib Sobol sensitivity on metrics
    # ----------------------------------------------------------
    param_names = list(bounds.keys())
    problem = {
        "num_vars": len(param_names),
        "names": param_names,
        "bounds": [list(bounds[p]) for p in param_names],
    }

    # Saltelli sampling
    X = saltelli.sample(problem, n_sobol_base, calc_second_order=False)
    # Evaluate model
    metrics_list = []
    for row in X:
        pars = dict(zip(param_names, row))
        # Simulate with these params
        with quiet():
            sim = sim_func(pars)
        if sim is None or sim.empty:
            metrics_list.append({k: np.nan for k in ["dynamic_range", "t50", "t10_90", "overshoot", "auc", "y_end"]})
            continue

        # Interpolate to a consistent grid for stable metrics
        tt = sim["time"].to_numpy(float)
        yy = sim["Y"].to_numpy(float)
        # sanitize
        m = np.isfinite(tt) & np.isfinite(yy)
        tt, yy = tt[m], yy[m]
        if tt.size < 5:
            metrics_list.append({k: np.nan for k in ["dynamic_range", "t50", "t10_90", "overshoot", "auc", "y_end"]})
            continue

        # Compute metrics on (tt,yy)
        metrics_list.append(metric_summary(yy, tt))

    metrics_df = pd.DataFrame(metrics_list)
    metrics_df.insert(0, "plasmid", plasmid)
    metrics_df.insert(0, "model", model_name)
    metrics_df.to_csv(OUT_PARAM / f"{model_name}_{plasmid}_sobol_samples_metrics.csv", index=False)

    # Sobol indices per metric
    sobol_rows = []
    for metric in ["dynamic_range", "t50", "t10_90", "overshoot", "auc", "y_end"]:
        Y = metrics_df[metric].to_numpy(float)
        ok = np.isfinite(Y)
        if ok.sum() < max(50, 0.3 * len(Y)):
            continue

        Si = sobol.analyze(problem, Y[ok], calc_second_order=False, print_to_console=False)
        for i, p in enumerate(param_names):
            sobol_rows.append({
                "model": model_name,
                "plasmid": plasmid,
                "metric": metric,
                "param": p,
                "S1": float(Si["S1"][i]),
                "S1_conf": float(Si["S1_conf"][i]),
                "ST": float(Si["ST"][i]),
                "ST_conf": float(Si["ST_conf"][i]),
            })

    sobol_df = pd.DataFrame(sobol_rows)
    sobol_df.to_csv(OUT_PARAM / f"{model_name}_{plasmid}_sobol_indices.csv", index=False)

    # Plot Sobol ST (total-order) per metric
    if not sobol_df.empty:
        for metric in sobol_df["metric"].unique():
            sub = sobol_df.query("metric == @metric").sort_values("ST", ascending=True)
            plt.figure(figsize=(7, 4))
            plt.barh(sub["param"], sub["ST"])
            plt.xlabel("Sobol ST (total-order)")
            plt.title(f"{model_name} {plasmid} – Sensitivity (metric={metric})")
            saveplot(OUT_PLOTS / f"{model_name}_{plasmid}_sobol_ST_{metric}.pdf")

    # ----------------------------------------------------------
    # B) Bootstrap parameter draws + predictive bands
    # ----------------------------------------------------------
    boot_rows = []
    boot_curves = []

    # Choose bootstrap scheme:
    # - half-times: normal(mean, se) if available
    # - K,n: lognormal around fitted with cv=0.20
    # - alpha: normal if you have per-plasmid spread; else lognormal cv=0.10
    # - delay: normal(mean, 0.5h) clipped [0,25]
    for b in range(int(n_boot)):
        pars_draw = {}
        for p in param_names:
            c = fitted[p]
            se = fitted.get(p + "_se", None)

            if p in ("t_up", "t_down"):
                pars_draw[p] = bootstrap_params(rng, c, se=se, kind="normal", clip_lo=1e-6)
            elif p in ("K", "n"):
                pars_draw[p] = bootstrap_params(rng, c, kind="lognormal", cv=0.20, clip_lo=1e-9)
            elif p == "alpha":
                # often global; keep conservative
                pars_draw[p] = bootstrap_params(rng, c, kind="lognormal", cv=0.10, clip_lo=1e-9)
            elif p.startswith("delay"):
                pars_draw[p] = float(np.clip(rng.normal(loc=float(c), scale=0.5), 0.0, 25.0))
            else:
                pars_draw[p] = bootstrap_params(rng, c, kind="lognormal", cv=0.20, clip_lo=1e-9)

        with quiet():
            sim = sim_func(pars_draw)

        if sim is None or sim.empty:
            continue

        tt = sim["time"].to_numpy(float)
        yy = sim["Y"].to_numpy(float)
        m = np.isfinite(tt) & np.isfinite(yy)
        tt, yy = tt[m], yy[m]
        if tt.size < 5:
            continue

        # Interpolate onto t_grid for band aggregation
        y_grid = np.interp(t_grid, tt, yy, left=np.nan, right=np.nan)
        boot_curves.append(y_grid)

        met = metric_summary(yy, tt)
        row = {"model": model_name, "plasmid": plasmid, "boot": b, **pars_draw, **met}
        boot_rows.append(row)

    boot_df = pd.DataFrame(boot_rows)
    boot_df.to_csv(OUT_PARAM / f"{model_name}_{plasmid}_bootstrap_params_and_metrics.csv", index=False)

    # Parameter histograms
    if not boot_df.empty:
        for p in param_names:
            plt.figure(figsize=(6, 4))
            plt.hist(boot_df[p].dropna().to_numpy(float), bins=40)
            plt.xlabel(p)
            plt.ylabel("count")
            plt.title(f"{model_name} {plasmid} – bootstrap parameter distribution")
            saveplot(OUT_PLOTS / f"{model_name}_{plasmid}_boot_param_{p}.pdf")

        # Metric histograms
        for metric in ["dynamic_range", "t50", "t10_90", "overshoot", "auc", "y_end"]:
            if metric not in boot_df.columns:
                continue
            plt.figure(figsize=(6, 4))
            plt.hist(boot_df[metric].dropna().to_numpy(float), bins=40)
            plt.xlabel(metric)
            plt.ylabel("count")
            plt.title(f"{model_name} {plasmid} – bootstrap metric distribution")
            saveplot(OUT_PLOTS / f"{model_name}_{plasmid}_boot_metric_{metric}.pdf")

    # Predictive bands
    if len(boot_curves) > 50:
        B = np.vstack(boot_curves)  # shape: n_boot x len(t_grid)
        qlo = np.nanpercentile(B, 2.5, axis=0)
        q50 = np.nanpercentile(B, 50.0, axis=0)
        qhi = np.nanpercentile(B, 97.5, axis=0)

        band = pd.DataFrame({
            "model": model_name,
            "plasmid": plasmid,
            "time": t_grid,
            "Y_q2p5": qlo,
            "Y_q50": q50,
            "Y_q97p5": qhi,
        })
        band.to_csv(OUT_PARAM / f"{model_name}_{plasmid}_predictive_band.csv", index=False)

        # Plot band + observed mean points
        plt.figure(figsize=(7, 4))
        plt.fill_between(t_grid, qlo, qhi, alpha=0.3, label="95% predictive band")
        plt.plot(t_grid, q50, linewidth=2, label="median prediction")
        plt.scatter(t_obs, y_obs, s=18, alpha=0.9, label="observed mean")
        plt.xlabel("Time (h)")
        plt.ylabel("mCherry (fc.cherry)")
        plt.title(f"{model_name} {plasmid} – predictive band vs observed mean")
        plt.legend()
        saveplot(OUT_PLOTS / f"{model_name}_{plasmid}_predictive_band_overlay.pdf")

    return sobol_df, boot_df


# ---------------------------------------------------------------------
# Model-specific simulators
# ---------------------------------------------------------------------
def make_rev_simulator(fitted_row: dict):
    """
    Returns a function sim(pars_draw)->df(time,R,Y) for REV.
    Uses step2.simulate_ode with fitted R0/Y0 from data.

    Args:
        fitted_row (dict): Dictionary containing fitted parameters including R0 and Y0.

    Returns:
        callable: Function that takes parameter draws and returns simulation DataFrame.
    """

    def _sim(pars_draw: dict):
        """
        Simulate REV model with given parameter draws.

        Args:
            pars_draw (dict): Dictionary of parameter draws for simulation.

        Returns:
            pd.DataFrame: Simulation results with columns time, R, Y.
        """
        R0 = float(fitted_row["R0"])
        Y0 = float(fitted_row["Y0"])

        pars = {
            "t_down": float(pars_draw["t_down"]),
            "K": float(pars_draw["K"]),
            "n": float(pars_draw["n"]),
            "alpha": float(pars_draw["alpha"]),
        }
        delay = float(pars_draw["delay_rev"])
        return step2.simulate_ode(R0=R0, Y0=Y0, pars=pd.Series(pars), tmax=150.0, step=0.05, delay=delay)

    return _sim


def make_kd_simulator(fitted_row: dict):
    """
    Returns a function sim(pars_draw)->df(time,R,Y) for KD.

    Args:
        fitted_row (dict): Dictionary containing fitted parameters.

    Returns:
        callable: Function that takes parameter draws and returns simulation DataFrame.
    """

    def _sim(pars_draw: dict):
        """
        Simulate KD model with given parameter draws.

        Args:
            pars_draw (dict): Dictionary of parameter draws for simulation.

        Returns:
            pd.DataFrame: Simulation results with columns time, R, Y.
        """
        pars = {
            "t_up": float(pars_draw["t_up"]),
            "K": float(pars_draw["K"]),
            "n": float(pars_draw["n"]),
            "alpha": float(pars_draw["alpha"]),
        }
        delay = float(pars_draw["delay_kd"])
        return step3.simulate_ode(R0=0.0, Y0=0.0, pars=pd.Series(pars), tmax=150.0, step=0.05, delay=delay)

    return _sim


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():

    # Load fitted parameter tables
    up, down, hill, alpha, alpha_mode, alpha_val, delays_rev, delays_kd = load_fitted_tables()

    # Build datasets
    rev = build_rev_dataset()
    kd = build_kd_dataset()

    if rev.empty:
        raise RuntimeError("REV dataset is empty. Check nfc/timecourse paths and step2 loader.")

    # Normalize plasmid labels
    rev["plasmid"] = _normalize_plasmid_labels(rev["plasmid"])
    kd["plasmid"] = _normalize_plasmid_labels(kd["plasmid"])

    # Ensure numeric
    for df in (rev, kd):
        df["time"] = df["time"].astype(float)

    # Build per-plasmid fitted center values and SEs
    # Half-times SE are from Step 1a/1b outputs.
    up = up.rename(columns={"halftime": "t_up"})
    down = down.rename(columns={"halftime": "t_down"})

    # unify plasmid key
    up["plasmid"] = up["plasmid"].astype(str)
    down["plasmid"] = down["plasmid"].astype(str)
    hill["plasmid"] = hill["plasmid"].astype(str)
    delays_rev["plasmid"] = delays_rev["plasmid"].astype(str)
    delays_kd["plasmid"] = delays_kd["plasmid"].astype(str)

    # Merge a fitted table per model
    fitted_core = (hill[["plasmid", "K", "n"]]
                   .merge(up[["plasmid", "t_up"] + (["se"] if "se" in up.columns else [])], on="plasmid", how="left",
                          suffixes=("", "_tup"))
                   .merge(down[["plasmid", "t_down"] + (["se"] if "se" in down.columns else [])], on="plasmid",
                          how="left", suffixes=("", "_tdown"))
                   .merge(delays_rev[["plasmid", "d_rev"]].rename(columns={"d_rev": "delay_rev"}), on="plasmid",
                          how="left")
                   .merge(delays_kd[["plasmid", "d_rev"]].rename(columns={"d_rev": "delay_kd"}), on="plasmid",
                          how="left")
                   )

    # Attach alpha
    if alpha_mode == "global":
        fitted_core["alpha"] = float(alpha_val)
    else:
        fitted_core = fitted_core.merge(alpha[["plasmid", "alpha"]], on="plasmid", how="left")

    # Rename SE columns cleanly
    # up.se and down.se are both named "se" from their csv; after merges they may collide.
    # We normalize them to t_up_se and t_down_se if possible.
    if "se" in fitted_core.columns:
        # ambiguous; keep as t_up_se if t_up exists and t_down se got renamed
        fitted_core = fitted_core.rename(columns={"se": "t_up_se"})
    if "se_tdown" in fitted_core.columns:
        fitted_core = fitted_core.rename(columns={"se_tdown": "t_down_se"})
    if "se_tup" in fitted_core.columns:
        fitted_core = fitted_core.rename(columns={"se_tup": "t_up_se"})

    # Save centers used
    fitted_core.to_csv(OUT_PARAM / "fitted_parameter_centers.csv", index=False)

    all_sobol = []
    all_boot = []

    # Run per plasmid
    for pl in sorted(set(fitted_core["plasmid"].dropna().astype(str))):
        row = fitted_core.query("plasmid == @pl")
        if row.empty:
            continue
        row = row.iloc[0].to_dict()

        # Pull observed ICs for REV from time==0 mean (like GOF)
        rev_pl = rev.query("plasmid == @pl").copy()
        kd_pl = kd.query("plasmid == @pl").copy()

        if not rev_pl.empty:
            mean_fc = rev_pl.groupby("time", as_index=False)[["fc.cherry", "norm.bfp"]].mean()
            t0 = mean_fc.query("time == 0")
            if not t0.empty:
                R0 = float(t0["norm.bfp"].iloc[0])
                Y0 = float(t0["fc.cherry"].iloc[0])
            else:
                R0, Y0 = 0.0, 1.0
        else:
            R0, Y0 = 0.0, 1.0

        # -----------------------------------------------------------------
        # REV
        # -----------------------------------------------------------------
        if not rev_pl.empty and np.isfinite(_safe_float(row.get("t_down"))) and np.isfinite(
                _safe_float(row.get("delay_rev"))):
            fitted = {
                "t_down": float(row["t_down"]),
                "t_down_se": float(row.get("t_down_se", np.nan)),
                "K": float(row["K"]),
                "n": float(row["n"]),
                "alpha": float(row["alpha"]),
                "delay_rev": float(row["delay_rev"]),
            }
            fitted_row = {"R0": R0, "Y0": Y0}

            bounds = {
                "t_down": make_bounds(fitted["t_down"], se=fitted.get("t_down_se"), rel=0.35, lo=1e-4),
                "K": make_bounds(fitted["K"], se=None, rel=0.35, lo=1e-9),
                "n": make_bounds(fitted["n"], se=None, rel=0.35, lo=1e-6),
                "alpha": make_bounds(fitted["alpha"], se=None, rel=0.25, lo=1e-9),
                "delay_rev": (0.0, 25.0),
            }
            sim_func = make_rev_simulator(fitted_row)

            sob, boot = run_one_model(
                model_name="REV",
                plasmid=pl,
                data_df=rev_pl,
                fitted=fitted,
                bounds=bounds,
                sim_func=sim_func,
                n_sobol_base=ARGS.sobol_base,
                n_boot=ARGS.n_boot,
                seed=ARGS.seed
            )

            if sob is not None and not sob.empty:
                all_sobol.append(sob)
            if boot is not None and not boot.empty:
                all_boot.append(boot)

        # -----------------------------------------------------------------
        # KD
        # -----------------------------------------------------------------
        if not kd_pl.empty and np.isfinite(_safe_float(row.get("t_up"))) and np.isfinite(
                _safe_float(row.get("delay_kd"))):
            fitted = {
                "t_up": float(row["t_up"]),
                "t_up_se": float(row.get("t_up_se", np.nan)),
                "K": float(row["K"]),
                "n": float(row["n"]),
                "alpha": float(row["alpha"]),
                "delay_kd": float(row["delay_kd"]),
            }
            bounds = {
                "t_up": make_bounds(fitted["t_up"], se=fitted.get("t_up_se"), rel=0.35, lo=1e-4),
                "K": make_bounds(fitted["K"], se=None, rel=0.35, lo=1e-9),
                "n": make_bounds(fitted["n"], se=None, rel=0.35, lo=1e-6),
                "alpha": make_bounds(fitted["alpha"], se=None, rel=0.25, lo=1e-9),
                "delay_kd": (0.0, 25.0),
            }
            sim_func = make_kd_simulator({})

            sob, boot = run_one_model(
                model_name="KD",
                plasmid=pl,
                data_df=kd_pl,
                fitted=fitted,
                bounds=bounds,
                sim_func=sim_func,
                n_sobol_base=512,
                n_boot=2000,
                seed=17
            )
            if sob is not None and not sob.empty:
                all_sobol.append(sob)
            if boot is not None and not boot.empty:
                all_boot.append(boot)

    # Aggregate outputs
    if all_sobol:
        sob_all = pd.concat(all_sobol, ignore_index=True)
        sob_all.to_csv(OUT_PARAM / "sobol_indices_all.csv", index=False)

        # Global plot: average ST per param across plasmids, per model+metric
        for (model, metric), sub in sob_all.groupby(["model", "metric"]):
            agg = (sub.groupby("param", as_index=False)["ST"].mean()
                   .sort_values("ST", ascending=True))
            plt.figure(figsize=(7, 4))
            plt.barh(agg["param"], agg["ST"])
            plt.xlabel("Mean Sobol ST across plasmids")
            plt.title(f"{model} – global sensitivity summary (metric={metric})")
            saveplot(OUT_PLOTS / f"{model}_sobol_ST_mean_{metric}.pdf")

    if all_boot:
        boot_all = pd.concat(all_boot, ignore_index=True)
        boot_all.to_csv(OUT_PARAM / "bootstrap_params_and_metrics_all.csv", index=False)

    print("[DONE] Saved outputs to:")
    print("  -", OUT_PARAM)
    print("  -", OUT_PLOTS)


if __name__ == "__main__":
    main()

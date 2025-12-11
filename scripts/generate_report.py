#!/usr/bin/env python3
"""
Generate a comprehensive reStructuredText report for the CasTuner analysis.

Behavior:
- On first run:
    - Parse console.txt to extract stable metadata (NFC background, gating info,
      KD/Rev file counts, etc.).
    - Save parsed metadata to parameters/run_metadata.json.
    - Copy console.txt to report/console_snapshot.txt for archival.
- On subsequent runs:
    - Use only parameters/run_metadata.json.
    - console.txt is NOT required and NOT parsed again.

- Always:
    - Scan plots/ and parameters/ for outputs.
    - Summarize key CSV tables.
    - Emit report/report_auto.rst.
"""

from __future__ import annotations
import json
import re
from pathlib import Path
from textwrap import dedent
from typing import Dict, List, Tuple, Optional

import pandas as pd


ROOT = Path(".").resolve()
PLOTS_DIR = ROOT / "plots"
PARAM_DIR = ROOT / "parameters"
REPORT_DIR = ROOT / "report"
REPORT_PATH = REPORT_DIR / "report_auto.rst"

CONSOLE_LOG = ROOT / "console.txt"
RUN_METADATA_PATH = PARAM_DIR / "run_metadata.json"
CONSOLE_SNAPSHOT = REPORT_DIR / "console_snapshot.txt"


# ---------------------------------------------------------------------
# One-time metadata extraction from console.txt
# ---------------------------------------------------------------------
def parse_console_log_raw(text: str) -> Dict[str, str]:
    """
    Extract *stable* metadata from raw console log text.
    We keep this conservative and only record things that won't change
    across re-runs unless the underlying data changes.
    """
    summaries: Dict[str, str] = {}

    # NFC background line
    m = re.search(
        r"\[NFC\]\s*mBFP_neg=([0-9.eE+-]+),\s*mmCherry_neg=([0-9.eE+-]+)",
        text,
    )
    if m:
        bfp, mch = m.groups()
        summaries["nfc"] = (
            f"Negative control fluorescence background estimated from NFC "
            f"samples: median BFP ≈ {bfp}, median mCherry ≈ {mch}."
        )

    # Count KD and Rev time-course files
    kd_count = len(re.findall(r"timecourse_.*_KD_.*\.fcs", text))
    rev_count = len(re.findall(r"timecourse_.*_Rev_.*\.fcs", text))
    if kd_count or rev_count:
        pieces = []
        if kd_count:
            pieces.append(f"{kd_count} KD time-course files")
        if rev_count:
            pieces.append(f"{rev_count} Rev time-course files")
        summaries["timecourses"] = (
            "The pipeline processed "
            + " and ".join(pieces)
            + " after gating and singlet selection."
        )

    # Gating info
    if "[gate] boundary" in text or "[gate] singlet" in text:
        summaries["gating"] = (
            "Each FCS file was gated on FSC/SSC to remove debris and on "
            "FSC-H vs FSC-A / SSC-H vs SSC-A to select singlets. "
            "All medians are computed on this singlet population."
        )

    return summaries


def ensure_run_metadata() -> Dict[str, str]:
    """
    Ensure parameters/run_metadata.json exists and return its content.

    - If the JSON exists: just load and return it.
    - If not and console.txt exists:
        * Parse console.txt once, write JSON.
        * Copy console.txt to report/console_snapshot.txt.
    - If neither exists: return empty metadata dict.
    """
    PARAM_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    if RUN_METADATA_PATH.exists():
        try:
            return json.loads(RUN_METADATA_PATH.read_text())
        except Exception:
            # If corrupted, fall back to empty dict
            return {}

    if not CONSOLE_LOG.exists():
        return {}

    raw = CONSOLE_LOG.read_text(errors="ignore")
    meta = parse_console_log_raw(raw)

    # Persist metadata
    RUN_METADATA_PATH.write_text(json.dumps(meta, indent=2))

    # Snapshot console.txt once for archival (but not used by report logic)
    try:
        CONSOLE_SNAPSHOT.write_text(raw)
    except Exception:
        pass

    return meta


# ---------------------------------------------------------------------
# Helpers for CSV summaries
# ---------------------------------------------------------------------
def safe_read_csv(path: Path) -> Optional[pd.DataFrame]:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path)
    except Exception:
        return None


def summarize_half_times(path_up: Path, path_down: Path) -> str:
    df_up = safe_read_csv(path_up)
    df_down = safe_read_csv(path_down)
    if df_up is None or df_down is None:
        return (
            "Half-times for up- and downregulation were estimated per plasmid "
            "from KD and Rev time-courses."
        )

    def basic_stats(df: pd.DataFrame, label: str) -> str:
        ht_cols = [c for c in df.columns if "half" in c.lower() or "t50" in c.lower()]
        if not ht_cols:
            return f"{label}: {len(df)} rows were fitted."
        col = ht_cols[0]
        return (
            f"{label}: {len(df)} plasmid-level fits with median {col} "
            f"≈ {df[col].median():.2f} h (min {df[col].min():.2f}, "
            f"max {df[col].max():.2f})."
        )

    return (
        basic_stats(df_up, "Upregulation (KD recovery)") + " " +
        basic_stats(df_down, "Downregulation (Rev decay)")
    )


def summarize_hill_parameters(path: Path) -> str:
    df = safe_read_csv(path)
    if df is None:
        return (
            "Hill functions were fitted to steady-state dose-response data for "
            "all constructs."
        )
    n_rows = len(df)
    cols = [c for c in df.columns if c.lower() in {"k", "n", "tup", "t_up", "alpha"}]
    col_str = ", ".join(cols) if cols else "multiple kinetic and shape parameters"
    return (
        f"Hill curves were fitted for {n_rows} construct/condition combinations, "
        f"estimating {col_str} from steady-state data."
    )


def summarize_delays(path_derep: Path, path_rep: Path) -> str:
    df_de = safe_read_csv(path_derep)
    df_re = safe_read_csv(path_rep)
    if df_de is None or df_re is None:
        return (
            "Time delays for derepression and repression were scanned on a grid "
            "to minimize ODE model errors."
        )

    def _summary(df: pd.DataFrame, label: str) -> str:
        num = df.select_dtypes("number")
        if num.empty:
            return f"{label} delays were estimated for {len(df)} plasmids."
        col = num.columns[0]
        return (
            f"{label} delays: median {col} ≈ {num[col].median():.2f} h "
            f"(range {num[col].min():.2f}–{num[col].max():.2f} h; n={len(df)})."
        )

    return _summary(df_de, "Derepression") + " " + _summary(df_re, "Repression")


def summarize_noise(path_ts: Path, path_hier: Path) -> str:
    df_ts = safe_read_csv(path_ts)
    df_h = safe_read_csv(path_hier)
    if df_ts is None or df_h is None:
        return (
            "Single-cell noise (CV²) was quantified from flow cytometry time-series, "
            "and a hierarchical model was used to decompose intrinsic vs extrinsic "
            "noise."
        )
    return (
        f"Noise analysis included {len(df_ts)} time-series entries and "
        f"{len(df_h)} hierarchical fits. CV² was evaluated as a function of "
        "mean BFP/mCherry across constructs."
    )


def summarize_candidates(path: Path) -> Tuple[str, str]:
    """
    Return (short_summary, table_rst) for candidate_selection_top10.csv
    """
    df = safe_read_csv(path)
    if df is None or df.empty:
        return (
            "The final design stage ranked tuner designs and selected the top 10 "
            "candidates.",
            "",
        )

    # Pick likely ID and score columns
    id_col = None
    score_col = None
    for c in df.columns:
        cl = c.lower()
        if id_col is None and any(k in cl for k in ["id", "plasmid", "construct", "name"]):
            id_col = c
        if score_col is None and any(k in cl for k in ["score", "rank", "composite", "objective"]):
            score_col = c

    if id_col and score_col and id_col != score_col:
        cols_to_show = [id_col, score_col]
    else:
        cols_to_show = list(df.columns[:3])

    short_summary = (
        f"The design space scan yielded {len(df)} top-ranking tuner designs; "
        "the table below lists the top 10 by composite score (or equivalent ranking "
        "metric)."
    )

    sub = df[cols_to_show].head(10)

    col_widths = [
        max(len(str(x)) for x in [col] + sub[col].astype(str).tolist())
        for col in cols_to_show
    ]

    def sep(char: str = "=") -> str:
        return " ".join(char * w for w in col_widths)

    lines = []
    lines.append(sep("="))
    header = " ".join(col.ljust(w) for col, w in zip(cols_to_show, col_widths))
    lines.append(header)
    lines.append(sep("="))
    for _, row in sub.iterrows():
        lines.append(
            " ".join(str(v).ljust(w) for v, w in zip(row.tolist(), col_widths))
        )
    lines.append(sep("="))

    table_rst = "\n".join(lines)
    return short_summary, table_rst


# ---------------------------------------------------------------------
# Plot categorisation
# ---------------------------------------------------------------------
def list_plots() -> Dict[str, List[Path]]:
    """
    Group plots into logical categories based on filename patterns.
    Returned paths are relative to project root.
    """
    categories: Dict[str, List[Path]] = {
        "hill": [],
        "kd_fits": [],
        "rev_fits": [],
        "kd_ode": [],
        "rev_delay_hill": [],
        "mae": [],
        "noise": [],
        "goodness_of_fit": [],
        "design_space": [],
        "other": [],
    }

    if not PLOTS_DIR.exists():
        return categories

    for p in PLOTS_DIR.rglob("*"):
        if not p.is_file():
            continue
        rel = p.relative_to(ROOT)
        name = p.name
        nlow = name.lower()

        if "goodness_of_fit" in rel.parts:
            categories["goodness_of_fit"].append(rel)
            continue
        if "noise_kinetics" in rel.parts:
            categories["noise"].append(rel)
            continue

        if nlow.startswith("hill_") or "hill" in nlow:
            categories["hill"].append(rel)
        elif nlow.startswith("kd_") and "fitting" in nlow:
            categories["kd_fits"].append(rel)
        elif nlow.startswith("rev_") and "fitting" in nlow:
            categories["rev_fits"].append(rel)
        elif nlow.startswith("kd_ode"):
            categories["kd_ode"].append(rel)
        elif nlow.startswith("rev_mcherry") or nlow.startswith("rev_tagbfp"):
            categories["rev_delay_hill"].append(rel)
        elif "mae_" in nlow:
            categories["mae"].append(rel)
        elif name == "design_space_map.pdf":
            categories["design_space"].append(rel)
        else:
            categories["other"].append(rel)

    return categories


def plots_to_rst(label: str, files: List[Path]) -> str:
    if not files:
        return ""
    lines = [label, "-" * len(label), ""]
    lines.append("The following plots are available for this section:")
    lines.append("")
    for f in sorted(files):
        lines.append(f"* ``{f.as_posix()}``")
    lines.append("")
    # Embed a few key examples
    for f in sorted(files)[:4]:
        lines.append(f".. image:: {f.as_posix()}")
        lines.append("   :width: 500px")
        lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------
def build_report() -> str:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    # One-time metadata extraction from console.txt (or load cached JSON)
    run_meta = ensure_run_metadata()

    plot_cats = list_plots()

    half_time_txt = summarize_half_times(
        PARAM_DIR / "half_times_upregulation.csv",
        PARAM_DIR / "half_times_downregulation.csv",
    )
    hill_txt = summarize_hill_parameters(PARAM_DIR / "Hill_parameters.csv")
    delay_txt = summarize_delays(
        PARAM_DIR / "delays_derepression.csv",
        PARAM_DIR / "delays_repression.csv",
    )
    noise_txt = summarize_noise(
        PARAM_DIR / "single_cell_noise_timeseries.csv",
        PARAM_DIR / "single_cell_noise_hierarchical.csv",
    )
    cand_summary, cand_table = summarize_candidates(
        PARAM_DIR / "candidate_selection_top10.csv"
    )

    intro = dedent(
        """
        CasTuner Analysis Pipeline Report
        =================================

        **Project:** CasTuner – Model-Driven Design of Analog Gene Tuners  
        **Author:** Abhinav Mishra

        Abstract
        --------
        This report was generated automatically from the outputs of the CasTuner
        analysis pipeline. The workflow starts from raw flow cytometry (FCS) files,
        performs robust gating and background subtraction, fits kinetic and
        dose-response models, simulates repression dynamics via ODEs, quantifies
        single-cell noise, and explores a model-driven design space to prioritize
        analog gene tuner constructs.

        """
    )

    # Data preprocessing / gating section using run_metadata
    preproc_lines = ["Data Preprocessing and Gating", "-----------------------------", ""]
    if "nfc" in run_meta:
        preproc_lines.append(f"* {run_meta['nfc']}")
    if "gating" in run_meta:
        preproc_lines.append(f"* {run_meta['gating']}")
    if "timecourses" in run_meta:
        preproc_lines.append(f"* {run_meta['timecourses']}")
    preproc_lines.append(
        "\nAll subsequent analyses use background-subtracted median fluorescence "
        "intensities for tagBFP and mCherry on the singlet population.\n"
    )
    preproc = "\n".join(preproc_lines)

    # Kinetic fitting
    kinetic = dedent(
        f"""
        1. Kinetic Fitting
        ------------------

        The first stage of the pipeline extracts basic kinetic parameters for each
        construct from KD (repressor induction) and Rev (repressor removal)
        time-courses.

        - KD (upregulation half-time): fits the recovery of the BFP reporter after dTAG
          withdrawal to estimate protein turnover / degradation dynamics.
        - Rev (downregulation half-time): fits the decay of BFP after dTAG addition.
        - Dose-response: fits Hill curves to day-4 steady-state measurements to extract
          repression strength (K) and cooperativity (n).

        Summary of fitted parameters:

        * {half_time_txt}
        * {hill_txt}

        The following plot groups illustrate per-construct kinetic fits and
        dose-response behavior across CasRx, dCas9 and chromatin-modifying fusions.
        """
    )

    kinetic_plots = "\n\n".join(
        [
            plots_to_rst("Dose-response (Hill fits)", plot_cats["hill"]),
            plots_to_rst("KD time-course fits", plot_cats["kd_fits"]),
            plots_to_rst("Rev time-course fits", plot_cats["rev_fits"]),
        ]
    )

    # ODE simulations and delays
    ode = dedent(
        f"""
        2. ODE Simulations and Delay Estimation
        ---------------------------------------

        Using the fitted kinetic parameters, the pipeline simulates the full CasTuner
        system via ordinary differential equations (ODEs). Separate state variables
        track tagBFP (degron-reporter) and mCherry (tuned gene) dynamics under KD and
        Rev conditions.

        For each construct, we scan delay parameters (e.g., chromatin remodeling lags)
        and select the value that minimizes the mean absolute error (MAE) between ODE
        trajectories and experimental time-series.

        Summary of delay estimates:

        * {delay_txt}

        The plots below show simulated vs observed trajectories and MAE-based
        diagnostics.
        """
    )

    ode_plots = "\n\n".join(
        [
            plots_to_rst("KD ODE trajectories", plot_cats["kd_ode"]),
            plots_to_rst("Rev delay & Hill overlays", plot_cats["rev_delay_hill"]),
            plots_to_rst("MAE diagnostics", plot_cats["mae"]),
        ]
    )

    # Noise section
    noise_section = dedent(
        f"""
        3. Single-Cell Noise Analysis
        -----------------------------

        Single-cell flow cytometry measurements are used to quantify expression
        noise (CV²) for each construct across conditions, and to decompose intrinsic
        and extrinsic contributions using a hierarchical model.

        {noise_txt}

        The noise plots illustrate how variability scales with mean expression and
        how noise depends on kinetic parameters such as half-times, delay and Hill
        coefficients (K, n).
        """
    )

    noise_plots = plots_to_rst("Noise and kinetic correlations", plot_cats["noise"])

    # Design space mapping & candidates
    design_section = dedent(
        f"""
        4. Design Space Mapping and Candidate Selection
        ------------------------------------------------

        Finally, the pipeline performs a model-driven scan of the tuner design space.
        By sampling combinations of key parameters (K, n, t_up, α) and propagating
        them through the ODE model, we obtain predicted dynamic range and response
        times (t_50) for a large ensemble of virtual designs.

        Experimentally characterized constructs are overlaid on this design space to
        benchmark achievable performance and to identify promising regions for new
        designs.

        {cand_summary}

        The plot below shows the global design space landscape, highlighting the
        trade-off between dynamic range and response speed. The accompanying table
        lists the top 10 prioritized tuner designs.
        """
    )

    design_plots = plots_to_rst("Design space summary", plot_cats["design_space"])

    if cand_table:
        design_table = "\nTop 10 candidate designs\n~~~~~~~~~~~~~~~~~~~~~~~~~\n\n" + cand_table + "\n"
    else:
        design_table = ""

    # Goodness-of-fit and other plots as an appendix section
    gof_section = plots_to_rst(
        "Goodness-of-fit diagnostics",
        plot_cats["goodness_of_fit"],
    )

    other_section = plots_to_rst("Additional plots", plot_cats["other"])

    parts = [
        intro,
        preproc,
        kinetic,
        kinetic_plots,
        ode,
        ode_plots,
        noise_section,
        noise_plots,
        design_section,
        design_plots,
        design_table,
        gof_section,
        other_section,
    ]

    return "\n\n".join(p for p in parts if p.strip())


def main() -> None:
    rst = build_report()
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(rst)
    print(f"[report] wrote {REPORT_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

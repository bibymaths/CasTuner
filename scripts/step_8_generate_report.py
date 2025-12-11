#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step 8 – Aggregated report & artifact sweep for CasTuner Python port.

Responsibilities
----------------
1. Collect all CSV parameter tables in `parameters/`.
2. Collect all plots under `plots/` (recursively).
   - Mirror them into `report/plots/` keeping subfolder structure.
   - For each PDF, try to create a PNG mirror using ImageMagick `convert`
     (if available).
3. Create a single PDF report at `report/CasTuner_summary_report.pdf` that:
   - Describes the analysis pipeline and main "knobs" (gating, models).
   - Summarizes key kinetic, Hill, ODE and noise parameters.
   - Summarizes design-space scan and selected top candidates.
   - Embeds a few key plots (if their PNG mirrors exist).

Assumptions
-----------
- You run this from the project root (same level as `scripts/`, `parameters/`,
  `plots/`, `fcs_files/`, etc.).
- The earlier pipeline steps have been executed, i.e. the CSVs and plots exist.
- `reportlab` is installed (`uv add reportlab` if needed).
- ImageMagick `convert` is available for PDF→PNG conversion; otherwise PNG
  mirroring is gracefully skipped with a warning.
"""

import sys
import shutil
import subprocess
import textwrap
from pathlib import Path
from typing import Optional,List

import numpy as np
import pandas as pd

from reportlab.lib.pagesizes import A4
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    PageBreak,
    Image,
    Table,
    TableStyle,
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors


# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------
# Project root is one level above scripts/
PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Try to read paths from config.yaml; fall back to "parameters" and "plots"
PARAM_PATH = PROJECT_ROOT / "parameters"
PLOTS_PATH = PROJECT_ROOT / "plots"

try:
    import yaml  # PyYAML
    cfg_path = PROJECT_ROOT / "config.yaml"
    if cfg_path.exists():
        cfg = yaml.safe_load(cfg_path.read_text())
        if "paths" in cfg:
            paths_cfg = cfg["paths"]
            if "parameters" in paths_cfg:
                PARAM_PATH = PROJECT_ROOT / paths_cfg["parameters"]
            if "plots" in paths_cfg:
                PLOTS_PATH = PROJECT_ROOT / paths_cfg["plots"]
except Exception as e:
    print(f"[warn] Could not parse config.yaml for paths: {e}")
    print("[warn] Falling back to ./parameters and ./plots relative to project root.")

print(f"[info] Using PARAM_PATH = {PARAM_PATH}")
print(f"[info] Using PLOTS_PATH = {PLOTS_PATH}")

REPORT_ROOT = PROJECT_ROOT / "report"
REPORT_PLOTS = REPORT_ROOT / "plots"
REPORT_TABLES = REPORT_ROOT / "tables"
REPORT_FILE = REPORT_ROOT / "CasTuner_summary_report.pdf"


# -------------------------------------------------------------------------
# Small utilities
# -------------------------------------------------------------------------
def ensure_dirs() -> None:
    REPORT_ROOT.mkdir(parents=True, exist_ok=True)
    REPORT_PLOTS.mkdir(parents=True, exist_ok=True)
    REPORT_TABLES.mkdir(parents=True, exist_ok=True)


def run_convert(pdf: Path, png: Path) -> bool:
    """
    Try to convert a single PDF into PNG using ImageMagick's `convert`.
    Returns True on success, False on failure or if convert is missing.
    """
    try:
        subprocess.run(
            ["convert", "-density", "300", str(pdf), "-quality", "95", str(png)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return True
    except FileNotFoundError:
        print("[warn] ImageMagick 'convert' not found – skipping PNG creation.")
        return False
    except subprocess.CalledProcessError:
        print(f"[warn] Failed to convert {pdf} → {png}")
        return False


def mirror_plots() -> None:
    """
    Mirror all plots from `plots/` into `report/plots/` and create PNG mirrors
    for each PDF (if possible).
    """
    print("[info] Mirroring plots into report/plots ...")
    for pdf in PLOTS_PATH.rglob("*.pdf"):
        rel = pdf.relative_to(PLOTS_PATH)
        target_pdf = REPORT_PLOTS / rel
        target_pdf.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(pdf, target_pdf)

        png_target = target_pdf.with_suffix(".png")
        if not png_target.exists():
            ok = run_convert(target_pdf, png_target)
            if ok:
                print(f"[info] Created PNG mirror: {png_target}")


def mirror_tables() -> None:
    """
    Mirror all CSVs from `parameters/` into `report/tables/`.
    """
    print("[info] Mirroring tables into report/tables ...")
    for csv in PARAM_PATH.glob("*.csv"):
        target_csv = REPORT_TABLES / csv.name
        shutil.copy2(csv, target_csv)


def load_csv(name: str) -> Optional[pd.DataFrame]:
    """
    Safe CSV loader from `parameters/`. Returns None if missing.
    """
    path = PARAM_PATH / name
    if not path.exists():
        print(f"[warn] Missing {name} – skipping its section.")
        return None
    df = pd.read_csv(path)
    if df.empty:
        print(f"[warn] {name} is empty.")
    return df


# -------------------------------------------------------------------------
# Biological summaries – lightweight but explicit
# -------------------------------------------------------------------------
def summarize_half_times(df: pd.DataFrame, label: str) -> str:
    """
    Generic summary for up/down half-times table.
    Expects columns: plasmid, halftime (or t_up / t_down).
    """
    cols = df.columns
    if "halftime" in cols:
        col = "halftime"
    elif "t_up" in cols:
        col = "t_up"
    elif "t_down" in cols:
        col = "t_down"
    else:
        return f"No recognizable half-time column found in {label}."

    df_use = df.dropna(subset=[col])
    if df_use.empty:
        return f"No finite half-times available for {label}."

    n = len(df_use)
    v = df_use[col].to_numpy(float)
    v_min, v_max = float(np.min(v)), float(np.max(v))
    v_med = float(np.median(v))

    fastest_row = df_use.iloc[int(np.argmin(v))]
    slowest_row = df_use.iloc[int(np.argmax(v))]
    fastest_pl = fastest_row.get("plasmid", "NA")
    slowest_pl = slowest_row.get("plasmid", "NA")

    return textwrap.dedent(
        f"""
        {label} ({n} constructs) show half-times ranging from {v_min:.2f} h
        to {v_max:.2f} h (median {v_med:.2f} h).
        The fastest response is observed for {fastest_pl} (~{v_min:.2f} h),
        while the slowest tuner is {slowest_pl} (~{v_max:.2f} h).
        Together, this defines the dynamic range of how quickly degron–Cas
        systems approach their new steady state after a perturbation.
        """
    ).strip()


def summarize_hill(df: pd.DataFrame) -> str:
    """
    Summary for Hill_parameters.csv table.
    Expects: plasmid, K or k, n.
    """
    df = df.copy()
    if "k" not in df.columns and "K" in df.columns:
        df = df.rename(columns={"K": "k"})
    if "k" not in df.columns or "n" not in df.columns:
        return "Hill_parameters.csv must contain columns 'K/k' and 'n'."

    n_rows = len(df)
    k_vals = df["k"].to_numpy(float)
    n_vals = df["n"].to_numpy(float)

    k_min, k_max = float(np.min(k_vals)), float(np.max(k_vals))
    n_min, n_max = float(np.min(n_vals)), float(np.max(n_vals))
    n_med = float(np.median(n_vals))

    steep_row = df.iloc[int(np.argmax(n_vals))]
    graded_row = df.iloc[int(np.argmin(n_vals))]

    steep_pl = steep_row.get("plasmid", "NA")
    steep_n = float(steep_row["n"])

    graded_pl = graded_row.get("plasmid", "NA")
    graded_n = float(graded_row["n"])

    return textwrap.dedent(
        f"""
        Steady-state dose–response fits (Hill curves) were obtained for {n_rows}
        constructs. Effective midpoints K span approximately {k_min:.3g}–{k_max:.3g},
        indicating how much dTAG-13 is required to tune the system into its
        responsive regime.

        Hill exponents n range from {n_min:.2f} to {n_max:.2f} (median {n_med:.2f}).
        The steepest, most switch-like response is seen in {steep_pl} (n ≈ {steep_n:.2f}),
        whereas {graded_pl} (n ≈ {graded_n:.2f}) behaves more like a graded analog tuner.
        Together, these parameters quantify how sharply degron–Cas constructs
        transition from "off" to "on" as the degrader dose increases.
        """
    ).strip()


def summarize_alpha(df: pd.DataFrame) -> str:
    if "alpha" not in df.columns:
        # try some common alt column names
        for alt in ("alphamcherry", "alpha_mcherry", "a"):
            if alt in df.columns:
                df = df.rename(columns={alt: "alpha"})
                break
    if "alpha" not in df.columns:
        return "Could not find an 'alpha' column in alphamcherry.csv."

    vals = df["alpha"].dropna().to_numpy(float)
    if vals.size == 0:
        return "No finite alpha values found."

    n = len(vals)
    a_min, a_max = float(np.min(vals)), float(np.max(vals))
    a_med = float(np.median(vals))

    t_half = np.log(2.0) / vals
    t_min, t_max = float(np.min(t_half)), float(np.max(t_half))

    return textwrap.dedent(
        f"""
        The mCherry degradation / dilution rate alpha was fitted from Rev
        time courses for {n} instances. Alpha values span {a_min:.3g}–{a_max:.3g} h⁻¹
        (median {a_med:.3g} h⁻¹), corresponding to effective reporter half-lives
        of roughly {t_min:.2f}–{t_max:.2f} hours. This sets how quickly the
        reporter output can track changes in repressor activity.
        """
    ).strip()


def summarize_delays(df_de: Optional[pd.DataFrame],
                     df_kd: Optional[pd.DataFrame]) -> str:
    parts: List[str] = []
    if df_de is not None and not df_de.empty:
        col = "d_rev" if "d_rev" in df_de.columns else df_de.columns[-1]
        vals = df_de[col].dropna().to_numpy(float)
        if vals.size:
            parts.append(
                f"Derepression delays (REV) range from {vals.min():.2f}–{vals.max():.2f} h."
            )
    if df_kd is not None and not df_kd.empty:
        col = "d_rev" if "d_rev" in df_kd.columns else df_kd.columns[-1]
        vals = df_kd[col].dropna().to_numpy(float)
        if vals.size:
            parts.append(
                f"Repression delays (KD) range from {vals.min():.2f}–{vals.max():.2f} h."
            )
    if not parts:
        return "No delay information was available for derepression or repression."
    return (
        "The ODE models required an explicit onset delay to align simulations with\n"
        + "the experimental time courses. "
        + " ".join(parts)
        + " These delays capture upstream processes such as degrader uptake, degron\n"
        + "processing and the time required until the effective nuclear repressor\n"
        + "concentration begins to change."
    )


def summarize_noise(noise_ts: Optional[pd.DataFrame],
                    noise_hier: Optional[pd.DataFrame]) -> str:
    if noise_ts is None or noise_ts.empty:
        return "Single-cell noise tables were not available."

    df = noise_ts.copy()
    # Expect: mean_BFP, cv2_BFP, mean_mCherry, cv2_mCherry, plasmid, exp, time
    text = []

    for ch, mean_col, cv_col in [
        ("BFP", "mean_BFP", "cv2_BFP"),
        ("mCherry", "mean_mCherry", "cv2_mCherry"),
    ]:
        if mean_col not in df.columns or cv_col not in df.columns:
            continue
        hi = df.groupby("plasmid")[cv_col].median().sort_values(ascending=False)
        lo = df.groupby("plasmid")[cv_col].median().sort_values(ascending=True)
        if hi.empty:
            continue

        noisy_pl = hi.index[0]
        quiet_pl = lo.index[0]
        noisy_cv = float(hi.iloc[0])
        quiet_cv = float(lo.iloc[0])

        text.append(
            f"For {ch}, median CV² across time and replicates spans roughly "
            f"{lo.iloc[0]:.3g}–{hi.iloc[-1]:.3g}. The noisiest construct "
            f"is {noisy_pl} (CV² ≈ {noisy_cv:.3g}), whereas {quiet_pl} shows "
            f"the quietest expression (CV² ≈ {quiet_cv:.3g})."
        )

    if not text:
        return "Could not extract channel-specific noise summaries."

    prefix = (
        "Single-cell flow cytometry events were propagated through the same gates\n"
        "as the kinetic analysis, and CV² was computed for both repressor (BFP)\n"
        "and reporter (mCherry). This quantifies how noisy each construct is at\n"
        "the single-cell level, beyond mean-level behavior.\n"
    )
    return prefix + "\n".join(text)


def summarize_design_space(design_df: Optional[pd.DataFrame],
                           top10_df: Optional[pd.DataFrame]) -> str:
    if design_df is None or design_df.empty:
        return "Design-space scan results were not found."

    txt = []

    for col in ["dynamic_range", "t50", "overshoot"]:
        if col not in design_df.columns:
            continue
        vals = design_df[col].dropna().to_numpy(float)
        if not vals.size:
            continue
        txt.append(
            f"{col} spans {vals.min():.3g}–{vals.max():.3g} "
            f"(median {np.median(vals):.3g})."
        )

    if top10_df is not None and not top10_df.empty:
        # Try to detect 'score' and 'design_id' or similar
        score_col = "Score" if "Score" in top10_df.columns else None
        if score_col is None:
            for c in top10_df.columns:
                if c.lower().startswith("score"):
                    score_col = c
                    break
        id_col = "design_id" if "design_id" in top10_df.columns else None
        if id_col is None and "plasmid" in top10_df.columns:
            id_col = "plasmid"

        if score_col is not None and id_col is not None:
            best = top10_df.sort_values(score_col, ascending=False).iloc[0]
            best_label = best[id_col]
            best_score = float(best[score_col])
            txt.append(
                f"The top-ranked candidate according to the dynamic-range/"
                f"response-time score is {best_label} (Score ≈ {best_score:.3g})."
            )

    intro = (
        "Using the measured parameter ranges as priors, a synthetic design space\n"
        "of repression ODE models was explored. For each sampled parameter set,\n"
        "the reporter trajectory Y(t) was summarized into dynamic range, response\n"
        "times (t50) and overshoot. This allows ranking hypothetical designs that\n"
        "would be difficult to realize or test exhaustively in the wet lab.\n"
    )

    if not txt:
        return intro + "Summary statistics could not be extracted from the table."
    return intro + " ".join(txt)


# -------------------------------------------------------------------------
# Experimental knobs: we try to introspect from your step scripts
# -------------------------------------------------------------------------
def build_knobs_paragraph() -> str:
    """
    Describe gating + model knobs by importing your step scripts.
    If imports fail, fall back to a generic description.
    """
    # try to import from scripts/
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
    try:
        import step_1a_fit_upregulation as s1a  # type: ignore
    except Exception:
        return textwrap.dedent(
            """
            All analyses share a consistent gating strategy and modeling
            backbone. Forward- and side-scatter channels define a coarse
            rectangle gate to remove debris, followed by an FSC-H/FSC-A
            singlet gate to exclude doublets. TagBFP (BV421-A) serves as the
            repressor proxy, and mCherry (PE-A) serves as the reporter.

            Kinetic parameters are obtained from exponential rise/decay fits
            to normalized BFP trajectories, while steady-state dose–response
            curves are captured by Hill functions. These parameters then feed
            into a two-variable ODE system (repressor R, reporter Y) used to
            estimate delays and to explore a broader design space.
            """
        ).strip()

    # If import worked, use the actual constants
    try:
        fsc_min = float(s1a.BOUND_MIN[s1a.CH_FSC_A])
        fsc_max = float(s1a.BOUND_MAX[s1a.CH_FSC_A])
        ssc_min = float(s1a.BOUND_MIN[s1a.CH_SSC_A])
        ssc_max = float(s1a.BOUND_MAX[s1a.CH_SSC_A])
        ratio_lo = float(s1a.SINGLET_RATIO_LOW)
        ratio_hi = float(s1a.SINGLET_RATIO_HIGH)
    except Exception:
        fsc_min = 0.4e5
        fsc_max = 2.0e5
        ssc_min = 0.2e5
        ssc_max = 1.3e5
        ratio_lo = 0.85
        ratio_hi = 1.15

    return textwrap.dedent(
        f"""
        Across all steps, the same analysis "knobs" are used to keep the
        quantitative picture consistent:

        • Boundary gate:
            FSC-A in [{fsc_min:.1f}, {fsc_max:.1f}],
            SSC-A in [{ssc_min:.1f}, {ssc_max:.1f}].
          This removes debris and very large aggregates.

        • Singlet gate:
            FSC-H/FSC-A constrained to [{ratio_lo:.2f}, {ratio_hi:.2f}],
            enriching for single cells.

        • Channels:
            tagBFP = BV421-A (repressor proxy),
            mCherry = PE-A (reporter output).

        • Kinetic models:
            – Up-regulation and down-regulation half-times are obtained from
              single-exponential fits to min–max normalized BFP trajectories.
            – At steady state (day 4), dTAG-13 dose–response is summarized
              by a Hill function (K, n) relating normalized BFP to reporter
              fold-change.

        • ODE backbone:
            A two-variable ODE system describes repressor R(t) and reporter
            Y(t), with synthesis/decay terms parameterized by the fitted
            half-times, Hill parameters and reporter decay constant alpha.
            Explicit onset delays align the simulations with the observed
            time courses.

        These knobs define a compact but biologically interpretable model
        of degron–Cas tuners as analog input–output devices.
        """
    ).strip()


# -------------------------------------------------------------------------
# Report building
# -------------------------------------------------------------------------
def add_heading(story, text, styles, lvl=1):
    size = 18 if lvl == 1 else 14
    style = ParagraphStyle(
        name=f"Heading{lvl}",
        parent=styles["Heading1"],
        fontSize=size,
        leading=size + 2,
        spaceAfter=6,
    )
    story.append(Paragraph(text, style))
    story.append(Spacer(1, 0.3 * cm))


def add_paragraph(story, text, styles):
    if isinstance(text, (list, tuple)):
        text = " ".join(str(t) for t in text if t is not None)
    else:
        text = str(text)

    for para in text.split("\n\n"):
        para = para.strip()
        if not para:
            continue
        story.append(Paragraph(para, styles["BodyText"]))
        story.append(Spacer(1, 0.2 * cm))


def add_plot_if_available(
    story, relative_plot: str, caption: str, styles, width_cm: float = 14.0
):
    """
    Embed a PNG version of a plot into the report if present in report/plots.
    `relative_plot` is relative to `plots/`, e.g. "noise_kinetics/noise_vs_mean_BFP.pdf".
    """
    png_path = REPORT_PLOTS / relative_plot
    png_path = png_path.with_suffix(".png")
    if not png_path.exists():
        return
    img = Image(str(png_path), width=width_cm * cm, height=(width_cm * 0.65) * cm)
    story.append(img)
    story.append(Spacer(1, 0.15 * cm))
    story.append(Paragraph(caption, styles["Italic"]))
    story.append(Spacer(1, 0.4 * cm))


def build_report():
    ensure_dirs()
    mirror_plots()
    mirror_tables()

    styles = getSampleStyleSheet()
    # Slightly denser body text
    styles["BodyText"].leading = 12

    story = []

    # ------------------------------------------------------------------
    # Title page
    # ------------------------------------------------------------------
    add_heading(story, "CasTuner Python Port – Integrated Quantitative Report", styles, lvl=1)
    add_paragraph(
        story,
        (
            "This report aggregates all analysis artifacts from the CasTuner "
            "Python pipeline: kinetic fits, dose–response models, ODE-based "
            "delay estimation, single-cell noise quantification and model-driven "
            "design-space exploration. The goal is to provide a compact but "
            "biologically interpretable summary of how degron–Cas constructs "
            "behave as analog gene tuners."
        ),
        styles,
    )
    story.append(PageBreak())

    # ------------------------------------------------------------------
    # Section: Experimental & modeling knobs
    # ------------------------------------------------------------------
    add_heading(story, "1. Experimental and Modeling Knobs", styles, lvl=1)
    knobs_text = build_knobs_paragraph()
    add_paragraph(story, knobs_text, styles)

    # ------------------------------------------------------------------
    # Section: Kinetic half-times
    # ------------------------------------------------------------------
    up = load_csv("half_times_upregulation.csv")
    down = load_csv("half_times_downregulation.csv")

    add_heading(story, "2. Kinetic Half-times (Up- and Down-regulation)", styles, lvl=1)

    if up is not None:
        add_heading(story, "2.1 Up-regulation after dTAG-13 withdrawal", styles, lvl=2)
        add_paragraph(story, summarize_half_times(up, "Up-regulation half-times"), styles)

    if down is not None:
        add_heading(story, "2.2 Down-regulation after dTAG-13 addition", styles, lvl=2)
        add_paragraph(story, summarize_half_times(down, "Down-regulation half-times"), styles)

    # ------------------------------------------------------------------
    # Section: Hill parameters
    # ------------------------------------------------------------------
    hill = load_csv("Hill_parameters.csv")
    if hill is not None:
        add_heading(story, "3. Steady-state Dose–Response (Hill Parameters)", styles, lvl=1)
        add_paragraph(story, summarize_hill(hill), styles)
        add_plot_if_available(
            story,
            "noise_kinetics/hill_K_vs_n.pdf",
            "Hill K vs n per plasmid, highlighting analog vs switch-like tuners.",
            styles,
        )

    # ------------------------------------------------------------------
    # Section: ODE parameters – alpha & delays
    # ------------------------------------------------------------------
    alpha_df = load_csv("alphamcherry.csv")
    delays_de = load_csv("delays_derepression.csv")
    delays_kd = load_csv("delays_repression.csv")

    add_heading(story, "4. ODE Parameters and Effective Delays", styles, lvl=1)

    if alpha_df is not None:
        add_heading(story, "4.1 Reporter Degradation (alpha)", styles, lvl=2)
        add_paragraph(story, summarize_alpha(alpha_df), styles)

    add_heading(story, "4.2 Onset Delays in Derepression and Repression", styles, lvl=2)
    add_paragraph(story, summarize_delays(delays_de, delays_kd), styles)

    # ------------------------------------------------------------------
    # Section: Single-cell noise
    # ------------------------------------------------------------------
    noise_ts = load_csv("single_cell_noise_timeseries.csv")
    noise_hier = load_csv("single_cell_noise_hierarchical.csv")

    add_heading(story, "5. Single-cell Noise Structure", styles, lvl=1)
    add_paragraph(story, summarize_noise(noise_ts, noise_hier), styles)

    add_plot_if_available(
        story,
        "noise_kinetics/noise_vs_mean_BFP.pdf",
        "Noise vs mean for BFP across time points and constructs.",
        styles,
    )
    add_plot_if_available(
        story,
        "noise_kinetics/noise_vs_mean_mCherry.pdf",
        "Noise vs mean for mCherry, highlighting output-level variability.",
        styles,
    )

    # ------------------------------------------------------------------
    # Section: Design-space scan & candidate selection
    # ------------------------------------------------------------------
    design = load_csv("design_space_scan_repression.csv")
    top10 = load_csv("candidate_selection_top10.csv")

    add_heading(story, "6. Model-driven Design Space and Candidate Tuners", styles, lvl=1)
    add_paragraph(story, summarize_design_space(design, top10), styles)

    add_plot_if_available(
        story,
        "design_space_map.pdf",
        "Design-space map: dynamic range vs response time with top-ranked designs.",
        styles,
        width_cm=15.0,
    )

    # Optionally show the top-10 table (trimmed)
    if top10 is not None and not top10.empty:
        add_paragraph(
            story,
            "Table 1 – Top-ranked candidate designs (truncated to first 6 columns):",
            styles,
        )
        cols = list(top10.columns[:6])
        data = [cols] + [list(map(str, row)) for row in top10[cols].values.tolist()]
        table = Table(data, hAlign="LEFT")
        table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                    ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
                    ("FONTSIZE", (0, 0), (-1, -1), 8),
                ]
            )
        )
        story.append(table)
        story.append(Spacer(1, 0.4 * cm))

    # ------------------------------------------------------------------
    # Section: Goodness-of-fit diagnostics (optional figure)
    # ------------------------------------------------------------------
    add_heading(story, "7. Goodness-of-fit Diagnostics", styles, lvl=1)
    add_paragraph(
        story,
        (
            "The fitted parameters are reused to reconstruct dose–response and "
            "time-course trajectories. Observed vs predicted comparisons (R², MAE) "
            "provide a sanity check that the mechanistic models capture the main "
            "features of the data without overfitting.",
        ),
        styles,
    )
    add_plot_if_available(
        story,
        "goodness_of_fit/gof_hill_fc_obs_vs_pred.pdf",
        "Observed vs predicted fold-change for Hill fits at steady state.",
        styles,
    )

    # ------------------------------------------------------------------
    # Final remarks
    # ------------------------------------------------------------------
    story.append(PageBreak())
    add_heading(story, "8. Outlook", styles, lvl=1)
    add_paragraph(
        story,
        (
            "Taken together, these analyses position the CasTuner constructs as "
            "controllable analog gene-tuning modules. The combination of kinetic "
            "half-times, steady-state dose–response, explicit delays and "
            "single-cell noise fingerprints provides a quantitative design space "
            "for proposing new degron architectures, guide configurations or "
            "promoter contexts.\n\n"
            "This report can be regenerated automatically whenever new FCS data "
            "or additional constructs are added, ensuring that modeling, noise "
            "analysis and design-space exploration stay tightly coupled to the "
            "experimental pipeline."
        ),
        styles,
    )

    # ------------------------------------------------------------------
    # Build PDF
    # ------------------------------------------------------------------
    print(f"[info] Writing report to {REPORT_FILE}")
    doc = SimpleDocTemplate(str(REPORT_FILE), pagesize=A4)
    doc.build(story)


if __name__ == "__main__":
    build_report()
    print("[info] Done.")
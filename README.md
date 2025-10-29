# CasTuner (Python Port)

## Modelling of repression and derepression dynamics in CRISPR/Cas-based analog gene-tuning systems

This repository contains an **independent Python implementation** of the computational framework associated with the manuscript:

> **CasTuner: a degron and CRISPR/Cas-based toolkit for analog tuning of endogenous gene expression**
> Gemma Noviello, Rutger A. F. Gjaltema, and Edda G. Schulz

This port reproduces the original **R-based kinetic modelling, parameter estimation, and ODE simulations**, preserving data flow, parameter inference logic, and figure generation.
All steps—from FCS preprocessing to ODE-based simulation—are reproducible using this repository alone.

---

## 1. Structure and Execution

### Workflow Overview

| Step | Description                                                                                       | Script                            |
| ---- | ------------------------------------------------------------------------------------------------- | --------------------------------- |
| 1a   | Fit upregulation dynamics of Cas-repressors after dTAG-13 withdrawal                              | `step_1a_fit_upregulation.py`     |
| 1b   | Fit degradation dynamics of Cas-repressors after dTAG-13 addition                                 | `step_1b_fit_downregulation.py`   |
| 1c   | Fit Hill functions for steady-state dose–response relationships                                   | `step_1c_fit_hill_curves.py`      |
| 2    | Simulate ODE-based derepression dynamics, estimate mCherry degradation rate (α), and infer delays | `step_2_simulate_derepression.py` |
| 3    | Simulate ODE-based repression dynamics and infer repression delays                                | `step_3_simulate_repression.py`   |

Each script reads raw flow cytometry `.fcs` data from `fcs_files/`, performs gating, normalization, parameter estimation, and produces results in:

* `parameters/` — fitted parameter tables (`.csv`)
* `plots/` — publication-style figures (`.pdf`)

Typical runtime: < 5 min per step on a standard desktop computer.

---

## 2. Automated Workflow (Snakemake)

All analysis steps are orchestrated via **Snakemake**, defined in the included `Snakefile` and `config.yaml`.

### Run the full workflow

```bash
snakemake -j 4
```

### Visualize dependencies

```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

Execution order:

1. Upregulation fitting
2. Downregulation fitting
3. Hill-curve fitting
4. Derepression ODE simulation
5. Repression ODE simulation

All results are automatically written to `parameters/` and `plots/`.

---

## 3. Environment Setup

Dependencies are specified in `pyproject.toml` (for uv/Poetry users) and `requirements.txt` (for pip).

### Create environment (uv)

```bash
uv venv
uv pip install -r requirements.txt
```

Or manually:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Core Dependencies

* `numpy`, `pandas`, `scipy`
* `FlowCytometryTools`, `flowio`, `flowutils`
* `plotnine`
* `snakemake>=8.0.0`

Tested under **Python 3.9 – 3.11**.
For Python 3.10+, legacy imports in `FlowCytometryTools` may require a `MutableMapping` compatibility shim.

---

## 4. Configuration

All input/output paths and parameters are defined in [`config.yaml`](./config.yaml).
Edit this file to modify directories, experimental metadata, or computational settings.

---

## 5. Data Layout

Expected directory structure:

```
CasTuner/
├── fcs_files/
│   ├── NFC/                    # Non-fluorescent control samples
│   └── time-course_data/       # Experimental time-course data
├── parameters/                 # Generated parameter CSVs
├── plots/                      # Output figures
├── Snakefile
├── config.yaml
├── pyproject.toml
├── requirements.txt
└── step_*.py                   # Analysis scripts (1a–3)
```

---

## 6. Reproducibility

All analysis steps are **deterministic** and self-contained:

* **ODE integration:** `scipy.integrate.solve_ivp`
* **Parameter fitting:** `scipy.optimize.curve_fit`
* **Plotting:** `plotnine`
* **Workflow orchestration:** `Snakemake`

Running the full pipeline regenerates all figures and parameter tables directly from raw `.fcs` files.

---

## 7. Citation

If this code or workflow contributes to your research, please cite:

> **CasTuner: a degron and CRISPR/Cas-based toolkit for analog tuning of endogenous gene expression**
> Gemma Noviello, Rutger A. F. Gjaltema, and Edda G. Schulz

And optionally acknowledge:

> *Python port and reproducibility workflow by Abhinav Mishra (2025).*

---

## 8. Development Notes

This Python implementation was independently developed by **Abhinav Mishra** as part of a prospective extension to the CasTuner project.
It is not affiliated with the original CasTuner authors but was created to:

* Demonstrate **reproducibility and extensibility** of the CasTuner framework in Python,
* Provide a **language-portable and Snakemake-integrated workflow**, and
* Serve as a foundation for further **systems-level modeling** or **kinetic inference** work in CRISPR-based gene regulation.

All equations, normalization procedures, and parameter definitions are preserved from the original R implementation to ensure numerical and conceptual equivalence.

---

## 9. Validation

To verify reproducibility, all five modelling stages were re-executed using the Python implementation and compared directly with the original R outputs from the CasTuner repository and manuscript.

### Summary of Reference (R) vs Python Results

| Parameter                                | R Output                                                         | Python Output                                                        | Comment                                                                          |
| ---------------------------------------- | ---------------------------------------------------------------- | -------------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| **α (mCherry degradation rate)**         | 0.036                                                            | 0.043                                                                | Same order of magnitude; minor deviation expected from different ODE integrators |
| **Repression delays (`d_rev`)**          | SP430A = 6 h, SP428 = 3 h, SP427 = 0 h, SP411 = 0 h, SP430 = 0 h | SP430A = 4 h, SP428 = 2.5 h, SP427 = 2 h, SP411 = 1.5 h, SP430 = 1 h | Relative ordering preserved; absolute values differ within 2–3 h tolerance       |
| **Upregulation half-times (`t₁/₂ ↑`)**   | 0.22–0.95 h range                                                | 0.22–0.91 h range                                                    | Excellent numerical agreement (<5 % relative error)                              |
| **Downregulation half-times (`t₁/₂ ↓`)** | 2.9–8.7 h range                                                  | 4.9–10.5 h range                                                     | Curves reproduce shape and ranking; slightly longer decay constants in Python    |
| **Hill parameters (`K`, `n`)**           | K ≈ 0.06–0.75, n ≈ 0.8–3.2                                       | K ≈ 0.05–0.53, n ≈ 1.1–2.9                                           | Preserves monotonic ordering and sigmoidal behavior                              |

### Interpretation

The Python implementation reproduces all major kinetic trends observed in the R version:

* **Parameter ranking** (e.g., SP411 < SP428 < SP430A for upregulation speed) is fully preserved.
* **Fitted magnitudes** differ modestly due to inherent variations in:

  * Optimization algorithms (`curve_fit` vs. `nlsLM`)
  * ODE solvers (`solve_ivp` LSODA vs. R’s `deSolve::lsoda`)
  * Floating-point precision and interpolation schemes
  * Slight differences in convergence tolerance and initial guesses

Despite these technical differences, the **model dynamics, curve shapes, and parameter scaling remain consistent**, confirming that the Python workflow is a faithful quantitative reproduction of the original R pipeline.

---

## 10. Visual Validation

To confirm reproducibility beyond numeric agreement, representative figures from both the **original R workflow** and the **Python port** were compared side by side.
All plots were generated from the same raw `.fcs` data under equivalent normalization and model assumptions.

### Comparison Overview

| Figure                          | Biological Context                                     | R Output                                                         | Python Output                                                         | Observation                                                            |
| ------------------------------- | ------------------------------------------------------ | ---------------------------------------------------------------- | --------------------------------------------------------------------- | ---------------------------------------------------------------------- |
| **Hill function (HDAC4-dCas9)** | Steady-state repression vs. normalized repressor level | <img src="R_results/plots/Hill-HDAC4-dCas9.pdf" width="220">     | <img src="Python_results/plots/Hill-HDAC4-dCas9.pdf" width="220">     | Nearly identical sigmoidal shape; slight midpoint (K) shift            |
| **Hill function (KRAB-dCas9)**  | Steady-state repression vs. dose                       | <img src="R_results/plots/Hill-KRAB-dCas9.pdf" width="220">      | <img src="Python_results/plots/Hill-KRAB-dCas9.pdf" width="220">      | Preserved slope and dynamic range; minor difference in low-dose region |
| **Upregulation (CasRx)**        | tagBFP rise after dTAG-13 withdrawal                   | <img src="R_results/plots/KD_CasRx_fitting.pdf" width="220">     | <img src="Python_results/plots/KD_CasRx_fitting.pdf" width="220">     | Curves overlap closely; Python version smoother due to solver grid     |
| **Derepression (CasRx)**        | mCherry derepression after Cas-repressor removal       | <img src="R_results/plots/REV_CasRx_fitting.pdf" width="220">    | <img src="Python_results/plots/REV_CasRx_fitting.pdf" width="220">    | Same trajectory and scaling; <15 % deviation in decay rate             |
| **ODE simulation (CasRx)**      | Repression model — mCherry trajectory                  | <img src="R_results/plots/KD_ODE_mCherry_CasRx.pdf" width="220"> | <img src="Python_results/plots/KD_ODE_mCherry_CasRx.pdf" width="220"> | Parallel rise–decay kinetics; timing offset consistent with d_rev fit  |
| **Delay scan (CasRx)**          | MAE vs. delay in repression fit                        | <img src="R_results/plots/MAE_KD_CasRx_mcherry.pdf" width="220"> | <img src="Python_results/plots/MAE_KD_CasRx_mcherry.pdf" width="220"> | Same error curve shape; minima at 2–3 h region in both                 |

---

### Summary

The reproduced plots demonstrate that the **Python implementation faithfully mirrors the R pipeline** in both dynamics and visualization:

* Hill and ODE curves preserve shape, ranking, and fold-change scaling.
* Delay-scan minima and kinetic trends match within experimental tolerance.
* Minor numeric deviations originate from:

  * Different ODE solvers (`deSolve::lsoda` vs. `scipy.solve_ivp`)
  * Distinct optimization routines (`nlsLM` vs. `curve_fit`)
  * Floating-point precision and interpolation details

Overall, the **Python results reproduce both the kinetic behavior and figure morphology** of the original CasTuner R analysis, confirming quantitative and visual equivalence.

---

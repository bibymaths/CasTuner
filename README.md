# CasTuner

## Modelling of repression and derepression dynamics in CRISPR/Cas-based analog gene tuning systems

This repository contains the **Python implementation** of all computational steps associated with the manuscript:

> **CasTuner: a degron and CRISPR/Cas-based toolkit for analog tuning of endogenous gene expression**
> Gemma Noviello, Rutger A. F. Gjaltema, and Edda G. Schulz

The Python version reproduces the kinetic modelling and ODE simulations originally implemented in R, maintaining identical data flow, parameter inference, and figure generation.
All components—data, scripts, configuration, and workflow—are fully accessible within this repository.

---

## 1. Structure and Execution

### Workflow overview

| Step | Description                                                                                            | Script                            |
| ---- | ------------------------------------------------------------------------------------------------------ | --------------------------------- |
| 1a   | Estimate upregulation dynamics of Cas-repressors after dTAG-13 withdrawal                              | `step_1a_fit_upregulation.py`     |
| 1b   | Estimate degradation dynamics of Cas-repressors after dTAG-13 addition                                 | `step_1b_fit_downregulation.py`   |
| 1c   | Fit steady-state Hill parameters for dose–response relationships                                       | `step_1c_fit_hill_curves.py`      |
| 2    | Simulate derepression ODE models, estimate mCherry degradation rate (α), and infer derepression delays | `step_2_simulate_derepression.py` |
| 3    | Simulate repression ODE models and infer repression delays                                             | `step_3_simulate_repression.py`   |

Each script reads raw flow-cytometry data from `fcs_files/`, performs parameter estimation, and produces results in:

* `parameters/` – numerical parameter estimates (`.csv`)
* `plots/` – generated figures (`.pdf`)

Typical runtime: < 5 minutes per step on a standard desktop system.

---

## 2. Automated Orchestration

All steps are orchestrated via **Snakemake** using the included `Snakefile` and `config.yaml`.

### Run the complete workflow

```bash
snakemake -j 4
```

### Inspect dependency graph

```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

Snakemake ensures the sequential execution:

1. Upregulation fitting
2. Downregulation fitting
3. Hill-curve fitting
4. Derepression simulation
5. Repression simulation

Results are written automatically to `parameters/` and `plots/`.

---

## 3. Environment and Dependencies

The project is defined through `pyproject.toml` (for uv) and `requirements.txt` (for pip users).

### Create the environment (uv)

```bash
uv venv
uv pip install -r requirements.txt
```

or equivalently:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

**Primary dependencies:**

* `pandas`, `numpy`, `scipy`
* `FlowCytometryTools`, `flowio`, `flowutils`
* `plotnine`
* `snakemake>=8.0.0`

---

## 4. Configuration

All user-defined paths and parameters are centralized in [`config.yaml`](./config.yaml).
Modify it to adjust input folders, output paths, or computational options.

---

## 5. Data Access

Expected directory layout:

```
CasTuner/
├── fcs_files/
│   ├── NFC/                    # Non-fluorescent control
│   └── time-course_data/       # Experimental flow cytometry data
├── parameters/                 # Fitted parameter CSVs
├── plots/                      # Output figures
├── Snakefile
├── config.yaml
├── pyproject.toml
├── requirements.txt
└── step_*.py                   # Analysis scripts
```

---

## 6. Reproducibility

All numerical steps are deterministic and self-contained:

* ODE integration: `scipy.integrate.odeint` / `solve_ivp`
* Parameter fitting: `scipy.optimize.curve_fit`
* Plotting: `plotnine`
* Workflow control: `Snakemake`

Results can be regenerated in full from raw `.fcs` data using only the provided scripts and configuration.

---

## 7. Citation

If this code or workflow contributes to your research, please cite:

> **CasTuner: a degron and CRISPR/Cas-based toolkit for analog tuning of endogenous gene expression**
> Gemma Noviello, Rutger A. F. Gjaltema, and Edda G. Schulz

---

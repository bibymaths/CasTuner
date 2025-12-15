# Data & Experiments

## Overview

The pipeline is driven by **flow cytometry time courses and dose–response experiments** for different CasTuner
constructs and cell lines.

The key data sources are:

- `.fcs` files with:
    - **non-fluorescent controls (NFC)**,
    - **time-course** experiments (repression/derepression),
    - **dose–response** experiments (dTAG titration).

All input paths are configured in `config.yaml`.

## Experimental conditions

### 1. Degron-controlled repressors

Cells express:

- FKBP12^F36V–hHDAC4–dCas9 (CasTuner transcriptional),
- FKBP12^F36V–CasRx (CasTuner post-transcriptional),
- or comparison systems (KRAB-dCas9, KRAB-Split-dCas9).

Each construct includes a **tBFP** reporter to quantify repressor abundance.

### 2. Endogenous reporters

Endpoints are endogenous genes tagged with fluorophores, e.g.:

- **Esrrb-P2A-mCherry** in mouse ESCs,
- **STAG2–EGFP** in HeLa cells,
- **Nanog-P2A-mCherry** in mouse ESCs.

Readouts:

- mCherry intensity → NANOG/ESRRB levels,
- EGFP intensity → STAG2 levels.

### 3. Time-course experiments

Two core dynamic experiments:

- **Repression (upregulation of repressor)**  
  Cells are shifted from high dTAG (repressor degraded) to lower dTAG.  
  Measured over ~6 days:
    - tBFP (repressor) trajectories,
    - mCherry/EGFP (target) trajectories.

- **Derepression (degradation of repressor)**  
  After sustained repression (low dTAG), dTAG is re-added.  
  Measured over ~6 days:
    - rapid loss of tBFP,
    - slower recovery of mCherry/EGFP.

These experiments constrain:

- **half-times** of degron-repressed proteins,
- **delays** between repressor change and target response.

### 4. Dose–response experiments

Cells are kept long enough at different dTAG doses to reach **steady state**:

- For each dose:
    - quantify repressor levels (tBFP),
    - quantify target output (mCherry or EGFP),
    - build **dose–response curves** (normalized BFP vs fold-change in target).

The pipeline fits Hill functions to these curves.

### 5. Single-cell noise and hierarchical structure

The same `.fcs` data also gives:

- full **single-cell distributions** of BFP and mCherry/EGFP,
- from which we derive per-condition:
    - mean,
    - variance,
    - **CV²**.

A simple hierarchical model summarises:

- construct-level noise,
- across time points and replicates.

### 6. Design-space scanning

Finally, measured parameters from the above experiments seed a **simulated design space**:

- we define realistic ranges for K, n, half-times, delays and degradation rate α,
- simulate the ODE model across this space,
- compute summary metrics (dynamic range, t₅₀, t₁₀–₉₀, overshoot),
- and rank designs.

All intermediate tables are written to `parameters/` and mirrored into the report.
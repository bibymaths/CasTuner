# CasTuner - Python Implementation
 
<p align="center">
  <img src="assets/logo.svg" alt="CasTuner Logo" width="400">
</p>

[![Deploy MkDocs to GitHub Pages](https://github.com/bibymaths/CasTuner/actions/workflows/deploy.yml/badge.svg)](https://github.com/bibymaths/CasTuner/actions/workflows/deploy.yml) 
[![Create Release](https://github.com/bibymaths/CasTuner/actions/workflows/release.yml/badge.svg)](https://github.com/bibymaths/CasTuner/actions/workflows/release.yml) 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/python-3.9%20|%203.10%20|%203.11-blue.svg)](https://www.python.org/downloads/)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17905630.svg)](https://doi.org/10.5281/zenodo.17905630)

## Quantitative modelling of analog gene tuning

CasTuner combines a degron-controlled CRISPR/Cas repressor with single-cell flow cytometry to **tune endogenous gene expression in an analog way** rather than as a simple ON/OFF switch.

This repository contains an **independent Python implementation** of the CasTuner analysis framework. It reproduces the original R-based pipeline for:

- kinetic fitting of degron–repressor dynamics,
- dose–response (Hill) curve fitting,
- ODE-based modelling of repression and derepression,
- single-cell noise analysis,
- model-driven design-space exploration, and
- automated reporting.

The focus of this documentation is not only *how* to run the code, but *why* each step exists and **what biological question it answers**.

## What problem does this solve?

Many processes in development and cell fate decisions are **dose-dependent**: small changes in the amount of a transcription factor can flip decisions such as pluripotency vs differentiation or dosage compensation. CasTuner was designed to measure and control such processes by:

- tuning Cas-derived repressors with a ligand-controlled degron,
- reading out endogenous reporters (e.g. Esrrb-mCherry, Nanog-mCherry, STAG2-EGFP),
- quantifying dynamics, dose–response and noise at single-cell resolution.

This Python port turns that biological framework into a **reproducible, Snakemake-driven workflow** that can be extended to new genes, constructs or model variants.

## How the documentation is organised

- **Biological Problem** – the conceptual and experimental context.
- **Data & Experiments** – what was measured and how it is organised.
- **Analysis Pipeline** – step-by-step explanation of each script and model:
  - Step 1 – Kinetic Fits (up/down half-times)
  - Step 2 – Dose–Response (Hill curves)
  - Step 3 – ODE Repression/Derepression
  - Step 4 – Single-cell Noise
  - Step 5 – Design Space Scan
  - Step 6 – Goodness of Fit
  - Step 7 – Design Selection
- **Reproducibility & Code** – environments, Snakemake, validation vs the original R workflow.
- **References** – primary CasTuner paper and related work.

For a quick practical overview (installation, Snakemake usage, validation tables), see the project-level `README.md`.

## Launch Binder

Note: It may take a few minutes to launch the environment on Binder. 

1. Notebook: Quickstart with Toy Dataset

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/bibymaths/CasTuner/main?urlpath=%2Fdoc%2Ftree%2Fnotebooks%2Fexample_notebook_00_quickstart_subset.ipynb)

2. Notebook: ODE Simulation with Toy Dataset

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/bibymaths/CasTuner/main?urlpath=%2Fdoc%2Ftree%2Fnotebooks%2Fexample_notebook_01_ode_steps_and_gof_subset.ipynb)

3. Notebook: Real FCS data analysis

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/bibymaths/CasTuner/main?urlpath=%2Fdoc%2Ftree%2Fnotebooks%2Fexample_notebook_02_real_fcs_subset.ipynb)
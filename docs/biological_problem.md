# Biological Problem

## Dose-dependent control of gene expression

Many cellular processes are **not** simple ON/OFF switches. They depend on *how much* of a given protein or RNA is
present:

- dosage compensation of X-linked genes,
- haploinsufficiency,
- fine-tuned thresholds in pluripotency,
- graded responses to morphogen gradients.

In these settings, changing a transcription factor by 20–30% can be enough to:

- keep cells pluripotent,
- push them into a specific lineage,
- or trigger compensation mechanisms on entire chromosomes.

To understand such systems, you need tools that:

1. **Tune endogenous gene expression** within a physiological range,
2. Do so **homogeneously** across cells (analog, not digital),
3. Are **reversible** and **quantitatively characterisable**.

## CasTuner: analog tuning instead of digital switching

CasTuner is built around degron-controlled Cas-based repressors:

- A FKBP12^F36V degron domain that makes the repressor level tunable by **dTAG-13**.
- Repressors such as **hHDAC4–dCas9** (transcriptional) or **CasRx** (post-transcriptional).
- Endogenous reporters tagged with fluorescent proteins (e.g. Esrrb-P2A-mCherry, STAG2–EGFP, Nanog–P2A-mCherry).

By titrating dTAG-13, you titrate **repressor abundance**, which in turn titrates **target gene expression**.
Importantly:

- KRAB-based systems tend to behave **digitally** (bimodal: ON vs OFF).
- hHDAC4–dCas9 and CasRx behave **analogically**: intermediate repressor levels lead to stable intermediate expression
  levels at single-cell resolution.

This makes CasTuner an ideal system to:

- map **dose–response curves** between transcription factors (e.g. NANOG, OCT4) and their targets,
- measure **kinetic responses** to degrader addition/withdrawal,
- quantify **single-cell noise** in analog repression,
- and **design new constructs** with desired dynamic range, speed and noise.

## Role of this Python port

The original CasTuner analysis was implemented in R. The Python port:

- re-implements the full modelling and analysis logic,
- uses **Snakemake** for workflow management,
- exposes all steps as **transparent, inspectable scripts**,
- and produces a final **summary report** linking:

    - experimental time courses,
    - fitted parameters (half-times, Hill K and n, delays, degradation rates),
    - ODE simulations,
    - noise statistics,
    - and a model-driven **design space map**.

The rest of this documentation walks through **how each step of the pipeline answers a concrete biological question**.
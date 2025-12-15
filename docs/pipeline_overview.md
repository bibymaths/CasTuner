# Analysis Pipeline Overview

The CasTuner Python port is structured as a **seven-step pipeline**, orchestrated via Snakemake:

1. **Kinetic fitting** of up- and down-regulation half-times (Step 1a–c),
2. **Hill curve** fitting for steady-state dose–response (Step 1c),
3. **ODE-based derepression and repression simulations** (Steps 2–3),
4. **Single-cell noise analysis** (Step 4a),
5. **Model-driven design-space scan** (Step 4b),
6. **Goodness-of-fit diagnostics** (Step 6),
7. **Design selection & mapping** (Step 7).

A final step collects all tables and plots into a single PDF report.

## High-level data flow

1. **Read and gate flow cytometry `.fcs` files**  
   All steps use consistent:
    - FSC/SSC gating,
    - singlet gating,
    - NFC-based background subtraction.

2. **Compute summary statistics per condition**
    - medians for BFP/mCherry/EGFP,
    - CV² for noise analysis,
    - dose-normalised readouts.

3. **Fit phenomenological models**
    - Exponential functions → **half-times** of BFP rise/decay,
    - Hill functions → **K** and **n** for steady-state dose–response.

4. **Fit mechanistic ODE model**
    - one model for **derepression**,
    - one model for **repression**,
    - shared parameters: Hill K, n and mCherry degradation rate α,
    - inferred delays for REV and KD.

5. **Propagate uncertainty & explore design space**
    - sample designs from measured parameter ranges,
    - simulate trajectories,
    - score candidates by dynamic range + speed.

6. **Summarise**
    - integrated plots:
        - kinetics vs noise,
        - K vs n,
        - delays vs half-times,
    - goodness-of-fit plots (obs vs pred),
    - final **design space map** and **top-10 candidate table**.

7. **Generate automated report**
    - Compile all results into a single PDF.
    - Sections for each analysis step with methods, results, interpretations.
    - Tables and plots embedded inline.

8. **Sensitivity & uncertainty analysis**
    - Quantify parameter sensitivities and output uncertainties.
    - Generate additional plots and tables.

Each of the following pages explains one step in more depth: what biological question it answers, what inputs it uses,
which models it fits, and how to interpret the outputs.
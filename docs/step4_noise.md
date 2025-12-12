# Step 4 – Single-Cell Noise & Hierarchical Analysis

**Script:**

- `step_4a_single_cell_hierarchical_noise.py`

## Biological question

> How noisy is analog tuning at the single-cell level, and how is this noise structured across time and constructs?

Analog tuning is only useful if intermediate levels are **stable and homogeneous** across cells. This step quantifies noise in:

- repressor levels (BFP),
- target levels (mCherry/EGFP).

## Inputs

- Time-course `.fcs` files,
- NFC controls,
- Paths from `config.yaml`,
- Previously defined gating strategy (reused).

## Method

1. **Per-cell data extraction**
   - apply FSC/SSC and singlet gates,
   - subtract NFC background,
   - keep BFP and mCherry/EGFP for each cell.

2. **Aggregate per condition**
   For each (construct, experiment, replicate, time):

   - compute:
     - mean,
     - variance,
     - coefficient of variation squared **CV²**:
       \[
       \text{CV}^2 = \frac{\sigma^2}{\mu^2}
       \]

3. **Hierarchical summary**
   - treat per-replicate CV² estimates as observations,
   - fit a simple hierarchical (normal–normal) model:
     - each construct has a latent “true” CV²,
     - replicates scatter around that with measurement noise.
   - summarise:
     - construct-level mean CV²,
     - uncertainty (e.g. SE, CI).

4. **Export**
   - `single_cell_noise_timeseries.csv` – raw time-series noise metrics,
   - `single_cell_noise_hierarchical.csv` – hierarchical summaries.

## Outputs

- Time-resolved noise trajectories,
- Per-construct noise summaries,
- Plots of noise vs mean BFP/mCherry.

## How to interpret

- Lower CV² at intermediate expression levels → **clean analog tuning**.
- Elevated CV²:
  - may indicate intrinsic stochasticity,
  - or subpopulations and hidden bistability.

In combination with Hill parameters and ODE fits, noise analysis helps answer:

- whether analog tuning is achieved without sacrificing stability,
- whether certain parameter regimes (fast vs slow, steep vs shallow) correlate with higher noise.

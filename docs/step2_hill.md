# Step 2 – Steady-State Dose–Response (Hill Curves)

**Script:**

- `step_1c_fit_hill_curves.py`

## Biological question

> Given a fixed construct, how does steady-state target expression depend on degron-controlled repressor abundance?

This defines the **input–output transfer function** of CasTuner:

- input: repressor abundance (from tBFP),
- output: target fold-change (mCherry or EGFP).

We model this with a Hill function to capture non-linearity and potential cooperativity.

## Inputs

- Dose–response `.fcs` data:
    - cells grown at different dTAG-13 concentrations until steady state (e.g. 4 days),
- NFC controls for background subtraction,
- Metadata describing which files correspond to which:
    - construct,
    - guide set (targeting vs NTC),
    - dTAG dose.

## Method

1. **Gating and background subtraction**
    - same gating strategy as in Step 1,
    - compute NFC-subtracted medians for BFP and mCherry/EGFP.

2. **Normalisation**
    - **Repressor axis (input)**:
        - per construct, apply min–max normalisation to BFP:
          $$
          x = \frac{BFP - BFP_\text{min}}{BFP_\text{max}-BFP_\text{min}}
          $$
            - this gives a dose-scaled “repressor activity” between 0 and 1.
    - **Target axis (output)**:
        - compute fold-change in mCherry/EGFP relative to NTC controls at high dTAG.

3. **Hill function fitting**

We fit:

$$
\text{FC}(x) = y_\text{min} +
\frac{y_\text{max} - y_\text{min}}{1 + \left(\frac{x}{K}\right)^n}
$$

where:

- `x` – normalized BFP (0–1),
- `K` – effective *midpoint* (dose at half-maximal effect),
- `n` – Hill coefficient (steepness, cooperativity),
- `y_min`, `y_max` – asymptotic fold-change bounds.

Parameters are estimated by non-linear least squares (`scipy.optimize.curve_fit`).

4. **Diagnostics**

- scatter plots of observed vs fitted fold-change,
- per-construct curves overlaid on data points,
- R² or MAE summaries.

## Outputs

- `parameters/Hill_parameters.csv` with columns:
    - construct, assay type,
    - `K`, `n`, `y_min`, `y_max`,
    - fit diagnostics.
- Plots in `plots/hill_curves/`.

## How to interpret

- **K (midpoint)**:
    - small K → construct responds strongly already at low repressor levels,
    - large K → construct needs higher repressor levels to repress strongly.

- **n (Hill coefficient)**:
    - `n ≈ 1` → graded, almost linear response,
    - `n >> 1` → switch-like behavior,
    - `n < 1` → shallower-than-linear response.

These Hill parameters are *reused directly* in the ODE model to connect repressor trajectories to mCherry dynamics.

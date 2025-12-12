# Step 6 – Goodness of Fit Diagnostics

**Script:**

- `step_6_goodness_of_fit.py`

## Biological question

> Do our fitted models (Hill, ODE) actually describe the data well enough to be trusted for design?

We want to check:

- how observed steady-state fold-changes compare to Hill predictions,
- how observed time courses compare to ODE simulations.

## Inputs

- `parameters/Hill_parameters.csv`
- `parameters/half_times_upregulation.csv`
- `parameters/half_times_downregulation.csv`
- `parameters/alphamcherry.csv`
- `parameters/delays_derepression.csv`
- `parameters/delays_repression.csv`
- Derived tables used earlier for dose–response and time-course fits.

## Method

1. **Rebuild predicted curves**

   - For each construct:
     - use Hill parameters to compute predicted fold-change across the observed BFP range,
     - use ODE parameters to re-simulate mCherry time courses.

2. **Overlay with observed data**

   - scatter plots of:
     - observed vs predicted fold-change (Hill),
     - observed vs predicted mCherry at each time point (ODE).

3. **Compute fit statistics**

   For each construct:

   - R² (coefficient of determination),
   - MAE or RMSE,
   - possibly residual diagnostics.

4. **Summarise**

   - per-construct summary table of goodness-of-fit,
   - global histogram/boxplots of fit quality.

## Outputs

- `plots/goodness_of_fit/gof_hill_fc_obs_vs_pred.pdf`
- Additional GOF plots for ODE vs data (depending on script options).

## How to interpret

- High R² and low MAE → model is capturing the essential behaviour.
- Systematic deviations:
  - might reveal missing processes (e.g. extra delays, feedback),
  - or parameter regimes where the simple Hill + first-order decay is insufficient.

This step is a sanity check before we take the design-space scan and design ranking too seriously.
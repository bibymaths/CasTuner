## Step 8 — Sensitivity, Uncertainty, and Distributions

Script: `scripts/step_8_sensitivity_uncertainty_and_distributions.py`

This step performs **global sensitivity analysis (Sobol)** and **uncertainty quantification (bootstrap)** for each
plasmid under both ODE modes:

* **REV** (derepression)
* **KD** (repression)

It quantifies which fitted parameters drive each performance metric, and how parameter uncertainty propagates into
predicted trajectories.

---

### What this step does

For each **(model, plasmid)** pair:

1. **Loads fitted parameter centers** from previous steps (half-times, Hill params, alpha, delays).
2. **Builds observed timecourse datasets**

    * REV dataset via functions imported from `step_2_simulate_derepression.py`
    * KD dataset via functions imported from `step_3_simulate_repression.py`
3. Runs **Sobol sensitivity (total-order ST)** over key metrics:

    * `dynamic_range`, `t50`, `t10_90`, `overshoot`, `auc`, `y_end`
4. Runs **bootstrap sampling** of parameters to generate:

    * parameter/metric distributions
    * **95% predictive bands** over time (with observed mean overlay)
5. Aggregates outputs into global summary CSVs and global mean-ST plots.

---

### Inputs

Expected files in `--params-dir` (default: `parameters/`):

* `half_times_upregulation.csv` *(plasmid, halftime, se)* → used as `t_up`
* `half_times_downregulation.csv` *(plasmid, halftime, se)* → used as `t_down`
* `Hill_parameters.csv` *(plasmid, K or k, n)*
* `alphamcherry.csv`

    * either one-row **global** `alpha`, or per-plasmid `(plasmid, alpha)`
* `delays_derepression.csv` *(plasmid, d_rev)*
* `delays_repression.csv` *(plasmid, d_rev)*

And dataset roots provided via CLI:

* `--nfc-dir` (default: `fcs_files/NFC`)
* `--timecourse-dir` (default: `fcs_files/time-course_data`)

---

### Outputs

Written under:

* Parameters: `parameters/uq_sensitivity/`
* Plots: `plots/uq_sensitivity/`

Per plasmid + model (examples shown):

**CSV outputs**

* `{MODEL}_{PLASMID}_sobol_samples_metrics.csv`
* `{MODEL}_{PLASMID}_sobol_indices.csv`
* `{MODEL}_{PLASMID}_bootstrap_params_and_metrics.csv`
* `{MODEL}_{PLASMID}_predictive_band.csv` *(if enough valid boot curves)*

**Plots**

* `{MODEL}_{PLASMID}_sobol_ST_{metric}.pdf` (bar chart)
* `{MODEL}_{PLASMID}_boot_param_{param}.pdf` (hist)
* `{MODEL}_{PLASMID}_boot_metric_{metric}.pdf` (hist)
* `{MODEL}_{PLASMID}_predictive_band_overlay.pdf` (band + observed mean)

Global aggregated outputs:

* `fitted_parameter_centers.csv`
* `sobol_indices_all.csv`
* `bootstrap_params_and_metrics_all.csv`
* `{MODEL}_sobol_ST_mean_{metric}.pdf` (mean ST across plasmids)

---

### Method details

**Sobol sampling**

* Uses SALib Saltelli sampling with `calc_second_order=False`.
* Parameter bounds are built from fitted centers:

    * half-times: `center ± 2*se` if `se` exists, else `± rel`
    * delays: fixed bounds `[0, 25]` hours
    * other params: relative bounds around fitted center (defaults in code).

**Bootstrap uncertainty**

* `t_up`, `t_down`: Normal(mean, se) (clipped positive) when SE available
* `K`, `n`: lognormal perturbation (default CV=0.20)
* `alpha`: lognormal perturbation (default CV=0.10)
* delays: Normal(mean, 0.5h) clipped to `[0, 25]`

**Metrics**
Computed from simulated `Y(t)` using:

* dynamic range, response times (t10/t50/t90), overshoot vs tail mean, AUC, end value.

**Parallelism**

* Parallelizes **across plasmids** using `joblib.Parallel(n_jobs=-1, backend="multiprocessing")`.

---


### How to interpret the outputs

* **Sobol ST plots** tell you which parameters dominate a metric:

    * e.g., if `n` and `K` have high ST for `dynamic_range`, dose-response shape controls amplitude.
* **Predictive band overlays** show whether the model’s uncertainty envelope is compatible with the observed mean
  trajectory.
* **Global mean ST plots** are the most report-friendly “one-glance” summary across constructs.

---

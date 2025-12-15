# Step 5 – Model-Driven Design Space Scan

**Script:**

- `step_4b_model_driven_design_scan.py`

## Biological question

> Given the parameter ranges we observe experimentally, which **hypothetical designs** would give us the best
> behaviour (dynamic range, speed, minimal overshoot)?

This step turns measurements into a **forward model** for exploring what CasTuner designs are even possible.

## Inputs

- `parameters/half_times_upregulation.csv`
- `parameters/half_times_downregulation.csv`
- `parameters/Hill_parameters.csv`
- `parameters/alphamcherry.csv`
- `parameters/delays_derepression.csv`
- `parameters/delays_repression.csv`

From these, the script derives realistic ranges for:

- Hill midpoint `K`,
- Hill coefficient `n`,
- half-times,
- α,
- delays.

## Method

1. **Define parameter ranges**

   For each parameter (e.g. K, n, t₁/₂, α, Δt), extract:

    - minimum,
    - maximum,
    - possibly trimmed quantiles to avoid outliers.

2. **Sample synthetic designs**

   Options include:

    - grid sampling,
    - Latin hypercube,
    - or a random sample within hyper-rectangular bounds.

   Each synthetic design = one parameter combination:

   $$
   \theta = (K, n, \alpha, t_{1/2,\uparrow}, t_{1/2,\downarrow}, \Delta t_\text{rev}, \Delta t_\text{kd})
   $$

3. **Simulate trajectories**

   For each \(\theta\):

    - simulate ODE during a repression protocol (or derepression, depending on design),
    - record `Y(t)` over a fixed window.

4. **Compute performance metrics**

   From each trajectory, extract:

    - dynamic range (max/min steady-state),
    - time to 50% change (**t₅₀**),
    - rise/decay times (e.g. t₁₀–₉₀),
    - overshoot or undershoot,
    - residual steady-state error.

   These are combined into one or multiple **scores**.

5. **Score & rank**

    - assign composite scores (e.g. “high dynamic range + low t₅₀”),
    - keep all designs but mark top performers.

## Outputs

- `parameters/design_space_scan_repression.csv` containing:
    - sampled parameter sets,
    - derived performance metrics,
    - ranking scores.

These results feed directly into Step 7 (design selection & map).

## How to interpret

This is where you can ask **“what if”** questions:

- *What parameter regimes give fast yet gentle repression?*
- *Is there a trade-off between dynamic range and noise or speed?*
- *Are our current constructs anywhere near the optimal region?*

Because the scan is based on **experimentally grounded ranges**, it remains biologically realistic while still exploring
combinations that may not yet exist in the current plasmid library.

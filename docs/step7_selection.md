# Step 7 – Design Selection & Mapping

**Script:**

- `step_7_design_selection_and_map.py`

## Biological question

> Among all real and simulated designs, which ones are the **best candidates** for experimental follow-up?

Using results from the design-space scan, we construct:

- a **design space map**,
- a **ranked shortlist** of top designs.

## Inputs

- `parameters/design_space_scan_repression.csv`
- `parameters/Hill_parameters.csv`
- `parameters/half_times_upregulation.csv`
- `parameters/alphamcherry.csv`
- `plots/goodness_of_fit/gof_hill_fc_obs_vs_pred.pdf` (to ensure we only consider well-behaved models).

## Method

1. **Filter designs**

   - optionally remove poor-fitting constructs or parameter sets with bad GOF,
   - focus on designs that:
     - operate in a desired dynamic-range window,
     - avoid extreme delays or unrealistic half-times.

2. **Define scoring criteria**

   Example components:

   - dynamic range score (higher is better),
   - speed score (lower t₅₀ is better),
   - smoothness/overshoot penalties,
   - consistency with experimentally observed behaviour.

   Composite score `S` can be a weighted sum or multi-objective ranking.

3. **Rank designs**

   - assign scores to real constructs and simulated parameter sets,
   - rank by `S` descending.

4. **Map design space**

   - visualise designs in 2D or 3D projections, e.g.:
     - K vs n coloured by performance,
     - t₅₀ vs dynamic range,
     - delay vs half-time.

## Outputs

- `plots/design_space_map.pdf` – a global map of design space.
- `parameters/candidate_selection_top10.csv` – table of top-ranked designs.

Each row of `candidate_selection_top10.csv` typically includes:

- construct ID or synthetic design label,
- key parameter values (K, n, half-times, delays, α),
- performance metrics and final score.

## How to interpret

This is the **actionable output** for wet-lab work:

- which constructs to build/test next,
- which parameter regimes to try to realise via new guide designs, linkers or degron variants.

It also gives a compact, visual summary of:

- where current constructs sit in the space of all plausible CasTuner designs,
- whether there is unused “room” for faster, cleaner, or more dynamic tuning.
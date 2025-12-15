# Step 3 – ODE Modelling of Derepression and Repression

**Scripts:**

- `step_2_simulate_derepression.py`
- `step_3_simulate_repression.py`

## Biological question

> Can a simple ODE model explain the observed dynamics of derepression and repression, and what delays and decay rates
> does it imply?

We want a **mechanistic** description of:

- how target gene expression responds over time to changes in repressor,
- using parameters already measured (half-times, Hill K and n).

## Inputs

- `parameters/half_times_upregulation.csv`
- `parameters/half_times_downregulation.csv`
- `parameters/Hill_parameters.csv`
- Time-course `.fcs` data for:
    - **derepression** (repressor going down),
    - **repression** (repressor going up),
- Configuration of which construct/experiment to treat as “reference” for α.

## Model

We model:

- `R(t)` – repressor level (e.g. BFP proxy),
- `Y(t)` – mCherry output (target gene).

### 1. Repressor dynamics

From Step 1, we have effective half-times `t₁/₂↑` and `t₁/₂↓`. We approximate:

- For upregulation:
  \[
  \frac{dR}{dt} = k_\text{up} (R_\text{max} - R)
  \]
- For downregulation:
  \[
  \frac{dR}{dt} = -k_\text{down} R
  \]

with `k = ln(2) / t₁/₂`.

In practice, `R(t)` is often reconstructed directly from data or parameterised implicitly.

### 2. Target dynamics

We assume:

- mCherry production is a Hill function of repressor:
  \[
  f(R) = \frac{1}{1 + (R/K)^n}
  \]
- mCherry decays with rate α.

Basic ODE:

\[
\frac{dY}{dt} = \beta f(R(t - \Delta t)) - \alpha Y
\]

where:

- `Δt` – delay between repressor change and effect on target,
- `β` – effective production rate.

α is shared across constructs and estimated from a reference time course.

## Derepression (Step 2)

- Use derepression time-course:
    - repressor drops after dTAG addition,
    - target recovers.

Procedure:

1. Fix Hill K and n from the dose–response fit.
2. Choose a reference construct/time course to fit **α (mCherry degradation rate)**.
3. For each construct:
    - integrate ODE for candidate delays Δt over a time grid,
    - compare simulated Y(t) to measured fold-change curve,
    - pick Δt that minimises MAE.
4. Save per-construct delays for derepression.

Outputs:

- `parameters/alphamcherry.csv`
- `parameters/delays_derepression.csv`

## Repression (Step 3)

Similar logic, but now:

- repressor rises after dTAG withdrawal,
- target decays.

Procedure:

1. Fix α and Hill parameters.
2. For each construct, scan Δt to best match mCherry time course under repression.
3. Save per-construct delays.

Outputs:

- `parameters/delays_repression.csv`
- ODE fit plots (obs vs simulated) in `plots/ode/`.

## How to interpret

- α (mCherry degradation rate) sets the **baseline turnover** of the reporter.
- Larger `Δt` values indicate:
    - additional layers between repressor and reporter (e.g. chromatin or network delays),
    - or slower establishment of the effective repression state.

Comparing Δt for derepression vs repression helps distinguish:

- immediate mechanisms (e.g. CasRx RNA degradation),
- slower chromatin rewiring (e.g. hHDAC4 vs KRAB),
- and potential short-term epigenetic memory.
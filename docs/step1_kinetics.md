# Step 1 – Kinetic Fits (Upregulation & Downregulation)

**Script(s):**

- `step_1a_fit_upregulation.py`
- `step_1b_fit_downregulation.py`
- `step_1c_fit_hill_curves.py` (dose–response part is covered in the next page)

## Biological question

> How fast do degron-controlled Cas repressors turn **ON** and **OFF** in response to dTAG changes?

This is captured by:

- **upregulation half-time** `t₁/₂↑`: how quickly tBFP (repressor) rises after removing dTAG,
- **downregulation half-time** `t₁/₂↓`: how quickly tBFP decays after adding dTAG.

These half-times set the **timescale** for gene repression/derepression.

## Inputs

- Gated FCS data for:
    - **upregulation** time courses (dTAG withdrawal),
    - **downregulation** time courses (dTAG addition),
- NFC samples for background subtraction (from `config.yaml`),
- Experimental metadata (mapping file names to constructs/time points).

## Method

1. **Gating and preprocessing**
    - boundary gate (FSC/SSC),
    - singlet gate (FSC-H/FSC-A),
    - remove NFC background,
    - compute **median tBFP** per (construct, experiment, time).

2. **Normalisation**
    - for each construct, rescale BFP into \[0, 1]:
        - 0 → lowest observed level (fully degraded),
        - 1 → plateau level (steady state).

3. **Exponential kinetic fits**
    - For upregulation (`t₁/₂↑`):
        - fit an exponential rise function:
          \[
          BFP(t) = 1 - \exp(-k_{\text{up}} t)
          \]
    - For downregulation (`t₁/₂↓`):
      \[
      BFP(t) = \exp(-k_{\text{down}} t)
      \]
    - Convert rate constants to half-times:
      \[
      t_{1/2} = \frac{\ln 2}{k}
      \]

4. **Quality checks**
    - keep only constructs with enough time points and good fit quality,
    - plot fits overlaid with data for visual inspection.

## Outputs

- `parameters/half_times_upregulation.csv`
- `parameters/half_times_downregulation.csv`
- Diagnostic plots for each construct in `plots/kinetics/` (mirrored into the final report).

Each row includes:

- construct ID,
- experiment/replicate,
- estimated `t₁/₂↑` or `t₁/₂↓`,
- fit error metrics.

## How to interpret

- Small `t₁/₂↑` → repressor accumulates quickly after dTAG removal.
- Small `t₁/₂↓` → repressor is removed quickly after dTAG re-addition.
- Comparing constructs tells you:
    - which design gives **fast control**,
    - whether different systems (hHDAC4–dCas9, CasRx, KRAB) differ primarily at the degron level or in downstream
      chromatin/RNA processing.

These half-times are later used as fixed or initial parameters in the ODE model.
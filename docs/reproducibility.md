# Reproducibility & Code

## Philosophy

This Python port was written to be:

- **transparent** – every step is in a plain Python script,
- **reproducible** – Snakemake encodes the full dependency graph,
- **extensible** – you can plug in new models, new datasets, or new constructs without rewriting everything.

## Environments

You can reproduce the full analysis using:

- `uv` or `pip` with `requirements.txt`,
- Python 3.9–3.11.

All versions are pinned in `pyproject.toml` / `requirements.txt`.

## Workflow (Snakemake)

The core analysis is orchestrated by the `Snakefile`:

- `snakemake -j 4`  
  runs the full pipeline and generates:

    - intermediate parameter tables in `parameters/`,
    - diagnostic plots in `plots/`,
    - final integrated report in `report/`.

The rules:

- encode dependencies between steps (1–7),
- enforce that later steps only run once earlier parameters are available,
- clean up intermediate directories automatically if configured that way.

You can visualise the DAG with:

```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

## Determinism

The pipeline is designed to be **deterministic** given the same inputs:

* ODE integration uses `scipy.integrate.solve_ivp` with fixed settings,
* optimisation uses `scipy.optimize.curve_fit` with explicit initial guesses,
* plotting is scripted and does not involve random jitter.

Any remaining variability (e.g. due to local BLAS/OpenMP threading) should be minor and not affect qualitative
conclusions.

## Validation vs original R workflow

To confirm that the Python port faithfully reproduces the original CasTuner analysis, the following was compared:

* **mCherry degradation rate α**
* **Repression/derepression delays**
* **Up- and down-regulation half-times**
* **Hill parameters K and n**, including their ordering across constructs
* Representative **Hill curves, time courses and ODE trajectories**

Results:

* parameter *rankings* and curve shapes match very closely,
* absolute values differ modestly, consistent with:

    * different ODE solvers (`deSolve::lsoda` vs `solve_ivp`),
    * different non-linear optimisers and tolerance settings.

In short, the Python implementation captures the same **biological trends and quantitative behaviour** as the original R
code, while being easier to integrate into modern Python-based analysis pipelines.

## Extending the pipeline

You can:

* add new steps as additional scripts + Snakemake rules,
* plug in alternative models (e.g. multi-state chromatin models) by reusing:

    * the data loading,
    * gating,
    * and plotting utilities already present.

The documentation for each step should give you enough context to do this without breaking the global narrative.



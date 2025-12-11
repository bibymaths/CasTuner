from snakemake.io import touch

configfile: "config.yaml"

# --------------------------------------------------------------------------
# Paths & helpers
# --------------------------------------------------------------------------
OUT_PLOTS  = config["paths"]["plots"]
OUT_PARAMS = config["paths"]["parameters"]

SCRIPTS = "scripts"
UV = "uv run python"


def script(name: str) -> str:
    """Return full path to a script inside the scripts/ directory."""
    return f"{SCRIPTS}/{name}"


# --------------------------------------------------------------------------
# All Targets
# --------------------------------------------------------------------------
rule all:
    input:
        f"{OUT_PARAMS}/half_times_upregulation.csv",
        f"{OUT_PARAMS}/half_times_downregulation.csv",
        f"{OUT_PARAMS}/Hill_parameters.csv",
        f"{OUT_PARAMS}/alphamcherry.csv",
        f"{OUT_PARAMS}/delays_derepression.csv",
        f"{OUT_PARAMS}/delays_repression.csv",
        f"{OUT_PARAMS}/single_cell_noise_timeseries.csv",
        f"{OUT_PARAMS}/single_cell_noise_hierarchical.csv",
        f"{OUT_PARAMS}/design_space_scan_repression.csv",
        f"{OUT_PLOTS}/noise_kinetics/noise_vs_mean_BFP.pdf",
        f"{OUT_PLOTS}/goodness_of_fit/gof_hill_fc_obs_vs_pred.pdf",
        f"{OUT_PLOTS}/design_space_map.pdf",
        f"{OUT_PARAMS}/candidate_selection_top10.csv"

# --------------------------------------------------------------------------
# Setup
# --------------------------------------------------------------------------
rule _mkdirs:
    output:
        touch(f"{OUT_PARAMS}/.ok"),
        touch(f"{OUT_PLOTS}/.ok")
    shell:
        f"mkdir -p {OUT_PARAMS} {OUT_PLOTS}"

# --------------------------------------------------------------------------
# Step 1: Kinetic Fitting (forced order inside Step 1)
# --------------------------------------------------------------------------
rule fit_upregulation:
    input:
        f"{OUT_PARAMS}/.ok"
    output:
        f"{OUT_PARAMS}/half_times_upregulation.csv"
    shell:
        "{uv} {script} && test -s {{output}}".format(
            uv=UV,
            script=script("step_1a_fit_upregulation.py"),
        )


rule fit_downregulation:
    input:
        f"{OUT_PARAMS}/.ok",
        # force this after fit_upregulation
        f"{OUT_PARAMS}/half_times_upregulation.csv"
    output:
        f"{OUT_PARAMS}/half_times_downregulation.csv"
    shell:
        "{uv} {script} && test -s {{output}}".format(
            uv=UV,
            script=script("step_1b_fit_downregulation.py"),
        )


rule fit_hill_curves:
    input:
        f"{OUT_PARAMS}/.ok",
        # force this after fit_downregulation
        f"{OUT_PARAMS}/half_times_downregulation.csv"
    output:
        f"{OUT_PARAMS}/Hill_parameters.csv"
    shell:
        "{uv} {script} && test -s {{output}}".format(
            uv=UV,
            script=script("step_1c_fit_hill_curves.py"),
        )

# --------------------------------------------------------------------------
# Step 2 & 3: ODE Simulation (Estimating Delays & Alpha)
# --------------------------------------------------------------------------
rule simulate_derepression:
    input:
        f"{OUT_PARAMS}/half_times_upregulation.csv",
        f"{OUT_PARAMS}/half_times_downregulation.csv",
        f"{OUT_PARAMS}/Hill_parameters.csv",
        f"{OUT_PARAMS}/.ok",
        f"{OUT_PLOTS}/.ok",
        lambda wc: config['paths']['nfc_dir'],
        lambda wc: config['paths']['timecourse_dir']
    output:
        f"{OUT_PARAMS}/alphamcherry.csv",
        f"{OUT_PARAMS}/delays_derepression.csv"
    shell:
        (
            "PLOTNINE_VERBOSE=0 {uv} {script} && "
            "test -s {{output[0]}} && test -s {{output[1]}}"
        ).format(
            uv=UV,
            script=script("step_2_simulate_derepression.py"),
        )


rule simulate_repression:
    input:
        f"{OUT_PARAMS}/half_times_upregulation.csv",
        f"{OUT_PARAMS}/half_times_downregulation.csv",
        f"{OUT_PARAMS}/Hill_parameters.csv",
        f"{OUT_PARAMS}/alphamcherry.csv",
        f"{OUT_PARAMS}/.ok",
        f"{OUT_PLOTS}/.ok",
        lambda wc: config['paths']['nfc_dir'],
        lambda wc: config['paths']['timecourse_dir']
    output:
        f"{OUT_PARAMS}/delays_repression.csv"
    shell:
        (
            "PLOTNINE_VERBOSE=0 {uv} {script} && "
            "test -s {{output}}"
        ).format(
            uv=UV,
            script=script("step_3_simulate_repression.py"),
        )

# --------------------------------------------------------------------------
# Step 4a: Single-Cell Noise Analysis
# --------------------------------------------------------------------------
rule analyze_single_cell_noise:
    input:
        f"{OUT_PARAMS}/.ok",
        # force this after simulate_repression (Step 3)
        f"{OUT_PARAMS}/delays_repression.csv",
        lambda wc: config['paths']['nfc_dir'],
        lambda wc: config['paths']['timecourse_dir']
    output:
        f"{OUT_PARAMS}/single_cell_noise_timeseries.csv",
        f"{OUT_PARAMS}/single_cell_noise_hierarchical.csv"
    shell:
        (
            "{uv} {script} && "
            "test -s {{output[0]}} && test -s {{output[1]}}"
        ).format(
            uv=UV,
            script=script("step_4a_single_cell_hierarchical_noise.py"),
        )

# --------------------------------------------------------------------------
# Step 4b: Model-Driven Design Scan
# --------------------------------------------------------------------------
rule scan_design_space:
    input:
        f"{OUT_PARAMS}/half_times_upregulation.csv",
        f"{OUT_PARAMS}/half_times_downregulation.csv",
        f"{OUT_PARAMS}/Hill_parameters.csv",
        f"{OUT_PARAMS}/alphamcherry.csv",
        f"{OUT_PARAMS}/delays_derepression.csv",
        f"{OUT_PARAMS}/delays_repression.csv",
        # force this after analyze_single_cell_noise (Step 4a)
        f"{OUT_PARAMS}/single_cell_noise_timeseries.csv"
    output:
        f"{OUT_PARAMS}/design_space_scan_repression.csv"
    shell:
        "{uv} {script} && test -s {{output}}".format(
            uv=UV,
            script=script("step_4b_model_driven_design_scan.py"),
        )

# --------------------------------------------------------------------------
# Step 5: Visualization of Noise & Kinetics
# --------------------------------------------------------------------------
rule visualize_noise_kinetics:
    input:
        f"{OUT_PARAMS}/single_cell_noise_timeseries.csv",
        f"{OUT_PARAMS}/single_cell_noise_hierarchical.csv",
        f"{OUT_PARAMS}/design_space_scan_repression.csv",
        f"{OUT_PARAMS}/Hill_parameters.csv",
        f"{OUT_PARAMS}/half_times_upregulation.csv",
        f"{OUT_PARAMS}/half_times_downregulation.csv",
        f"{OUT_PARAMS}/delays_derepression.csv",
        f"{OUT_PARAMS}/delays_repression.csv"
    output:
        f"{OUT_PLOTS}/noise_kinetics/noise_vs_mean_BFP.pdf"
    shell:
        (
            "PLOTNINE_VERBOSE=0 {uv} {script} && "
            "test -s {{output}}"
        ).format(
            uv=UV,
            script=script("step_5_visualize_noise_and_kinetics.py"),
        )

# --------------------------------------------------------------------------
# Step 6: Goodness of Fit Diagnostics
# --------------------------------------------------------------------------
rule goodness_of_fit:
    input:
        f"{OUT_PARAMS}/Hill_parameters.csv",
        f"{OUT_PARAMS}/half_times_upregulation.csv",
        f"{OUT_PARAMS}/half_times_downregulation.csv",
        f"{OUT_PARAMS}/alphamcherry.csv",
        f"{OUT_PARAMS}/delays_derepression.csv",
        f"{OUT_PARAMS}/delays_repression.csv",
        # force this after visualize_noise_kinetics (Step 5)
        f"{OUT_PLOTS}/noise_kinetics/noise_vs_mean_BFP.pdf"
    output:
        f"{OUT_PLOTS}/goodness_of_fit/gof_hill_fc_obs_vs_pred.pdf"
    shell:
        (
            "PLOTNINE_VERBOSE=0 {uv} {script} && "
            "test -s {{output}}"
        ).format(
            uv=UV,
            script=script("step_6_goodness_of_fit.py"),
        )

# --------------------------------------------------------------------------
# Step 7: Final Design Selection & Mapping
# --------------------------------------------------------------------------
rule select_designs:
    input:
        f"{OUT_PARAMS}/design_space_scan_repression.csv",
        f"{OUT_PARAMS}/Hill_parameters.csv",
        f"{OUT_PARAMS}/half_times_upregulation.csv",
        f"{OUT_PARAMS}/alphamcherry.csv",
        # force this after goodness_of_fit (Step 6)
        f"{OUT_PLOTS}/goodness_of_fit/gof_hill_fc_obs_vs_pred.pdf"
    output:
        f"{OUT_PLOTS}/design_space_map.pdf",
        f"{OUT_PARAMS}/candidate_selection_top10.csv"
    shell:
        (
            "{uv} {script} && "
            "test -s {{output[0]}}"
        ).format(
            uv=UV,
            script=script("step_7_design_selection_and_map.py"),
        )

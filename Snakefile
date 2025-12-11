from snakemake.io import touch
configfile: "config.yaml"

OUT_PLOTS  = config["paths"]["plots"]
OUT_PARAMS = config["paths"]["parameters"]
SCRIPTS = "scripts"

# -----------------------------------------------------------------------------
# All Targets
# -----------------------------------------------------------------------------
rule all:
    input:
        # Existing outputs
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

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
rule _mkdirs:
    output:
        touch(f"{OUT_PARAMS}/.ok"),
        touch(f"{OUT_PLOTS}/.ok")
    shell:
        "mkdir -p {OUT_PARAMS} {OUT_PLOTS}"

# -----------------------------------------------------------------------------
# Step 1: Kinetic Fitting
# -----------------------------------------------------------------------------
rule fit_upregulation:
    input: f"{OUT_PARAMS}/.ok"
    output: f"{OUT_PARAMS}/half_times_upregulation.csv"
    shell: f"uv run python {SCRIPTS}/step_1a_fit_upregulation.py"
            "test -s {output}"

rule fit_downregulation:
    input: f"{OUT_PARAMS}/.ok"
    output: f"{OUT_PARAMS}/half_times_downregulation.csv"
    shell: f"uv run python {SCRIPTS}/step_1b_fit_downregulation.py "
            "test -s {output}"

rule fit_hill_curves:
    input: f"{OUT_PARAMS}/.ok"
    output: f"{OUT_PARAMS}/Hill_parameters.csv"
    shell: f"uv run python {SCRIPTS}/step_1c_fit_hill_curves.py "
            "test -s {output}"

# -----------------------------------------------------------------------------
# Step 2 & 3: ODE Simulation (Estimating Delays & Alpha)
# -----------------------------------------------------------------------------
rule simulate_derepression:
    input:
        f"{OUT_PARAMS}/half_times_upregulation.csv",
        f"{OUT_PARAMS}/half_times_downregulation.csv",
        f"{OUT_PARAMS}/Hill_parameters.csv",
        f"{OUT_PARAMS}/.ok",
        f"{OUT_PLOTS}/.ok",
        # Config inputs trigger rebuild if paths change
        lambda wc: config['paths']['nfc_dir'],
        lambda wc: config['paths']['timecourse_dir']
    output:
        f"{OUT_PARAMS}/alphamcherry.csv",
        f"{OUT_PARAMS}/delays_derepression.csv"
    shell:
        "PLOTNINE_VERBOSE=0 uv run python step_2_simulate_derepression.py && "
        "test -s {output[0]} && test -s {output[1]}"

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
    output: f"{OUT_PARAMS}/delays_repression.csv"
    shell:
        "PLOTNINE_VERBOSE=0 uv run python step_3_simulate_repression.py "
        "test -s {output}"

# -----------------------------------------------------------------------------
# Step 4a: Single-Cell Noise Analysis
# -----------------------------------------------------------------------------
rule analyze_single_cell_noise:
    input:
        f"{OUT_PARAMS}/.ok",
        lambda wc: config['paths']['nfc_dir'],
        lambda wc: config['paths']['timecourse_dir']
    output:
        f"{OUT_PARAMS}/single_cell_noise_timeseries.csv",
        f"{OUT_PARAMS}/single_cell_noise_hierarchical.csv"
    shell:
        f"uv run python {SCRIPTS}/step_4a_single_cell_hierarchical_noise.py && "
        "test -s {output[0]} && test -s {output[1]}"

# -----------------------------------------------------------------------------
# Step 4b: Model-Driven Design Scan
# -----------------------------------------------------------------------------
rule scan_design_space:
    input:
        f"{OUT_PARAMS}/half_times_upregulation.csv",
        f"{OUT_PARAMS}/half_times_downregulation.csv",
        f"{OUT_PARAMS}/Hill_parameters.csv",
        f"{OUT_PARAMS}/alphamcherry.csv",
        f"{OUT_PARAMS}/delays_derepression.csv",
        f"{OUT_PARAMS}/delays_repression.csv"
    output:
        f"{OUT_PARAMS}/design_space_scan_repression.csv"
    shell:
        f"uv run python {SCRIPTS}/step_4b_model_driven_design_scan.py "
        "test -s {output}"

# -----------------------------------------------------------------------------
# Step 5: Visualization of Noise & Kinetics
# -----------------------------------------------------------------------------
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
        # Use one file as a representative timestamp for the directory
        f"{OUT_PLOTS}/noise_kinetics/noise_vs_mean_BFP.pdf"
    shell:
        "PLOTNINE_VERBOSE=0 uv run python step_5_visualize_noise_and_kinetics.py "
        "test -s {output}"

# -----------------------------------------------------------------------------
# Step 6: Goodness of Fit Diagnostics
# -----------------------------------------------------------------------------
rule goodness_of_fit:
    input:
        f"{OUT_PARAMS}/Hill_parameters.csv",
        f"{OUT_PARAMS}/half_times_upregulation.csv",
        f"{OUT_PARAMS}/half_times_downregulation.csv",
        f"{OUT_PARAMS}/alphamcherry.csv",
        f"{OUT_PARAMS}/delays_derepression.csv",
        f"{OUT_PARAMS}/delays_repression.csv"
    output:
        # Use one file as a representative timestamp for the directory
        f"{OUT_PLOTS}/goodness_of_fit/gof_hill_fc_obs_vs_pred.pdf"
    shell:
        "PLOTNINE_VERBOSE=0 uv run python step_6_goodness_of_fit.py "
        "test -s {output}"

# -----------------------------------------------------------------------------
# Step 7: Final Design Selection & Mapping
# -----------------------------------------------------------------------------
rule select_designs:
    input:
        f"{OUT_PARAMS}/design_space_scan_repression.csv",
        f"{OUT_PARAMS}/Hill_parameters.csv",
        f"{OUT_PARAMS}/half_times_upregulation.csv",
        f"{OUT_PARAMS}/alphamcherry.csv"
    output:
        f"{OUT_PLOTS}/design_space_map.pdf",
        f"{OUT_PARAMS}/candidate_selection_top10.csv"
    shell:
        f"uv run python {SCRIPTS}/step_7_design_selection_and_map.py" 
        "test -s {output[0]} && test -s {output[1]}"
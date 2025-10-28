from snakemake.io import touch
configfile: "config.yaml"

OUT_PLOTS  = config["paths"]["plots"]
OUT_PARAMS = config["paths"]["parameters"]

rule all:
    input:
        f"{OUT_PARAMS}/half_times_upregulation.csv",
        f"{OUT_PARAMS}/half_times_downregulation.csv",
        f"{OUT_PARAMS}/Hill_parameters.csv",
        f"{OUT_PARAMS}/alphamcherry.csv",
        f"{OUT_PARAMS}/delays_derepression.csv",
        f"{OUT_PARAMS}/delays_repression.csv",

rule _mkdirs:
    output:
        touch(f"{OUT_PARAMS}/.ok"),
        touch(f"{OUT_PLOTS}/.ok")
    shell:
        "mkdir -p {OUT_PARAMS} {OUT_PLOTS}"

rule fit_upregulation:
    input: f"{OUT_PARAMS}/.ok"
    output: f"{OUT_PARAMS}/half_times_upregulation.csv"
    shell: "uv run python step_1a_fit_upregulation.py && test -s {output}"

rule fit_downregulation:
    input: f"{OUT_PARAMS}/.ok"
    output: f"{OUT_PARAMS}/half_times_downregulation.csv"
    shell: "uv run python step_1b_fit_downregulation.py && test -s {output}"

rule fit_hill_curves:
    input: f"{OUT_PARAMS}/.ok"
    output: f"{OUT_PARAMS}/Hill_parameters.csv"
    shell: "uv run python step_1c_fit_hill_curves.py && test -s {output}"

rule simulate_derepression:
    input:
        f"{OUT_PARAMS}/half_times_upregulation.csv",
        f"{OUT_PARAMS}/half_times_downregulation.csv",
        f"{OUT_PARAMS}/Hill_parameters.csv",
        f"{OUT_PARAMS}/.ok",
        f"{OUT_PLOTS}/.ok",
        lambda wc: config['paths']['nfc_dir'],
        lambda wc: config['paths']['timecourse_dir'],
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
        lambda wc: config['paths']['timecourse_dir'],
    output: f"{OUT_PARAMS}/delays_repression.csv"
    shell:
        "PLOTNINE_VERBOSE=0 uv run python step_3_simulate_repression.py && test -s {output}"

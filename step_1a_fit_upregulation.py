import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from plotnine import (
    ggplot,
    aes,
    geom_point,
    geom_line,
    labs,
    scale_y_continuous,
)
from utils import (
    load_nfc_background,
    load_fcs_data,
    theme_castuner,
    color_bfp
)

# Optional: define or replace with your actual mapping
plasmid_map = {}  # e.g., {"pX": "Plasmid_X"}

# Exponential rise-to-max from 0 to 1 with half-time t1/2
def exp_fit_up(t, t1_2):
    yf, y0 = 1.0, 0.0
    return yf + (y0 - yf) * np.exp(-t * (np.log(2) / t1_2))

def main():
    OUT_PATH = "plots"
    PARAM_PATH = "parameters"
    os.makedirs(OUT_PATH, exist_ok=True)
    os.makedirs(PARAM_PATH, exist_ok=True)

    # Load background controls
    mBFP_neg, mmCherry_neg = load_nfc_background()

    # Load time-course data (expects directory to exist)
    medianexp = load_fcs_data("fcs_files/time-course_data", mBFP_neg, mmCherry_neg)

    # Parse filenames: format: _ _ plasmid _ exp _ rep _ time _ _
    parts = medianexp["filename"].str.split("_", expand=True)
    medianexp["plasmid"] = parts[2]
    medianexp["exp"] = parts[3]
    medianexp["time"] = pd.to_numeric(parts[5], errors="coerce")

    # Keep only KD (upregulation) experiment rows with valid time
    KD = medianexp[(medianexp["exp"] == "KD") & (~medianexp["time"].isna())].copy()

    # Min-max scale tagBFP (BV421-A) per plasmid using early and late time windows
    mean_final = (
        KD[KD["time"] > 10]
        .groupby("plasmid")["BV421-A"]
        .mean()
        .reset_index(name="mean_final")
    )
    mean_init = (
        KD[KD["time"] == 0]
        .groupby("plasmid")["BV421-A"]
        .mean()
        .reset_index(name="mean_init")
    )
    KD = KD.merge(mean_final, on="plasmid", how="left").merge(mean_init, on="plasmid", how="left")

    # Avoid div-by-zero; if mean_final == mean_init, set denom to NaN (row dropped by fitting)
    denom = (KD["mean_final"] - KD["mean_init"]).replace(0, np.nan)
    KD["norm.bfp"] = (KD["BV421-A"] - KD["mean_init"]) / denom

    plasmids = KD["plasmid"].dropna().unique()
    half_times = []

    for plasmid in plasmids:
        data = KD[(KD["plasmid"] == plasmid) & (~KD["norm.bfp"].isna())]
        if data.empty:
            print(f"Skipping {plasmid}: no valid normalized points.")
            continue

        # Fit exponential rise (initial guess t1/2 ~ 0.8 hours as in your script)
        try:
            popt, pcov = curve_fit(
                exp_fit_up,
                data["time"].to_numpy(dtype=float),
                data["norm.bfp"].to_numpy(dtype=float),
                p0=[0.8],
                bounds=(0, np.inf),
                maxfev=10000,
            )
            t1_2 = float(popt[0])
            se = float(np.sqrt(np.diag(pcov))[0]) if pcov.size else np.nan

            # Optional name mapping for output consistency
            param_plasmid_name = plasmid_map.get(plasmid, plasmid)
            half_times.append({"plasmid": param_plasmid_name, "halftime": t1_2, "se": se})

            # Build a smooth curve (replacement for geom_function/stat_function)
            t_min = float(data["time"].min())
            t_max = float(data["time"].max())
            curve_df = pd.DataFrame({"time": np.linspace(t_min, t_max, 200)})
            curve_df["norm_bfp"] = exp_fit_up(curve_df["time"], t1_2)

            # Plot points + fitted curve
            p = (
                ggplot(data, aes("time", "norm.bfp"))
                + geom_point(size=0.4, alpha=0.7, color=color_bfp)
                + geom_line(
                    data=curve_df,
                    mapping=aes(x="time", y="norm_bfp"),
                    color="black",
                )
                + labs(x="Time (hours)", y="tagBFP (% of final)")
                + scale_y_continuous(breaks=[0, 0.25, 0.5, 0.75, 1], limits=[-0.15, 1.2])
                + theme_castuner()
            )
            p.save(
                os.path.join(OUT_PATH, f"KD_{plasmid}_fitting.pdf"),
                width=1.5 * 1.618,
                height=1.5,
                units="in",
            )

        except RuntimeError:
            print(f"Could not fit {plasmid}")
        except Exception as e:
            print(f"Error fitting {plasmid}: {e}")

    # Save estimated half-times
    half_time_df = pd.DataFrame(half_times)
    half_time_df.to_csv(os.path.join(PARAM_PATH, "half_times_upregulation.csv"), index=False)
    print("Step 1a: Upregulation fitting complete.")
    if not half_time_df.empty:
        print(half_time_df)

if __name__ == "__main__":
    main()

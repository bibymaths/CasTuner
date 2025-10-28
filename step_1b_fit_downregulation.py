import pandas as pd
import numpy as np
import os
from scipy.optimize import curve_fit
from plotnine import (
    ggplot,
    aes,
    geom_point,
    labs,
    scale_y_continuous,
    geom_function,
)
from utils import (
    load_nfc_background,
    load_fcs_data,
    theme_castuner,
    color_bfp,
    plasmid_map
)


# Define the fitting function for downregulation
def exp_fit_down(t, t1_2):
    """Exponential decay from max (1) to min (0)."""
    yf = 0
    y0 = 1
    return yf + (y0 - yf) * np.exp(-t * (np.log(2) / t1_2))


def main():
    OUT_PATH = "plots"
    PARAM_PATH = "parameters"
    os.makedirs(OUT_PATH, exist_ok=True)
    os.makedirs(PARAM_PATH, exist_ok=True)

    # Load background
    mBFP_neg, mmCherry_neg = load_nfc_background()

    # Load time-course data
    medianexp = load_fcs_data("fcs_files/time-course_data", mBFP_neg, mmCherry_neg)

    # Parse filenames
    parts = medianexp['filename'].str.split('_', expand=True)
    medianexp['plasmid'] = parts[2]
    medianexp['exp'] = parts[3]
    medianexp['time'] = pd.to_numeric(parts[5])

    # Filter for REV (derepression/downregulation) experiment
    REV = medianexp[medianexp['exp'] == 'Rev'].copy()

    # Min-max scaling of bfp
    # Note: R script uses mean.final as min and mean.init as max
    mean_final = REV[REV['time'] > 10].groupby('plasmid')['BV421-A'].mean().reset_index(name='mean_final')
    mean_init = REV[REV['time'] == 0].groupby('plasmid')['BV421-A'].mean().reset_index(name='mean_init')

    REV = pd.merge(REV, mean_final, on='plasmid', how='left')
    REV = pd.merge(REV, mean_init, on='plasmid', how='left')

    REV['norm.bfp'] = (REV['BV421-A'] - REV['mean_final']) / (REV['mean_init'] - REV['mean_final'])

    plasmids = REV['plasmid'].unique()
    half_times = []

    for plasmid in plasmids:
        data = REV[REV['plasmid'] == plasmid]

        # Fit data
        try:
            popt, pcov = curve_fit(exp_fit_down, data['time'], data['norm.bfp'], p0=[0.1])
            t1_2 = popt[0]
            se = np.sqrt(np.diag(pcov))[0]

            param_plasmid_name = plasmid_map.get(plasmid, plasmid)
            half_times.append({'plasmid': param_plasmid_name, 'halftime': t1_2, 'se': se})

            # Plot
            p = (
                    ggplot(data, aes('time', 'norm.bfp'))
                    + geom_point(size=0.4, alpha=0.7, color=color_bfp)
                    + geom_function(
                fun=lambda t: exp_fit_down(t, t1_2),
                color="black"
            )
                    + labs(x="Time (hours)", y="tagBFP (% of final)")
                    + scale_y_continuous(breaks=[0, 0.25, 0.5, 0.75, 1], limits=[-0.15, 1.4])
                    + theme_castuner()
            )
            p.save(
                os.path.join(OUT_PATH, f"REV_{plasmid}_fitting.pdf"),
                width=1.5 * 1.618,
                height=1.5,
                units="in"
            )
        except RuntimeError:
            print(f"Could not fit {plasmid}")

    # Save parameters
    half_time_df = pd.DataFrame(half_times)
    half_time_df.to_csv(
        os.path.join(PARAM_PATH, "half_times_downregulation.csv"),
        index=False
    )
    print("Step 1b: Downregulation fitting complete.")
    print(half_time_df)


if __name__ == "__main__":
    main()
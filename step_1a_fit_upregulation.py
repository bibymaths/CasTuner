import pandas as pd
import numpy as np
import os
from scipy.optimize import curve_fit
from plotnine import (
    ggplot,
    aes,
    geom_point,
    geom_line,
    labs,
    scale_y_continuous,
    geom_function,
)
from utils import (
    load_nfc_background,
    load_fcs_data,
    theme_castuner,
    color_bfp
)


# Define the fitting function for upregulation
def exp_fit_up(t, t1_2):
    """Exponential rise to max (1) from min (0)."""
    yf = 1
    y0 = 0
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
    # Filename format: _ _ plasmid _ exp _ rep _ time _ _
    parts = medianexp['filename'].str.split('_', expand=True)
    medianexp['plasmid'] = parts[2]
    medianexp['exp'] = parts[3]
    medianexp['time'] = pd.to_numeric(parts[5])

    # Filter for KD (knockdown/upregulation) experiment
    KD = medianexp[medianexp['exp'] == 'KD'].copy()

    # Min-max scaling of bfp
    mean_final = KD[KD['time'] > 10].groupby('plasmid')['BV421-A'].mean().reset_index(name='mean_final')
    mean_init = KD[KD['time'] == 0].groupby('plasmid')['BV421-A'].mean().reset_index(name='mean_init')

    KD = pd.merge(KD, mean_final, on='plasmid', how='left')
    KD = pd.merge(KD, mean_init, on='plasmid', how='left')

    KD['norm.bfp'] = (KD['BV421-A'] - KD['mean_init']) / (KD['mean_final'] - KD['mean_init'])

    plasmids = KD['plasmid'].unique()
    half_times = []

    for plasmid in plasmids:
        data = KD[KD['plasmid'] == plasmid]

        # Fit data
        try:
            popt, pcov = curve_fit(exp_fit_up, data['time'], data['norm.bfp'], p0=[0.8])
            t1_2 = popt[0]
            se = np.sqrt(np.diag(pcov))[0]

            # Map plasmid name for consistency with R script outputs
            param_plasmid_name = plasmid_map.get(plasmid, plasmid)
            half_times.append({'plasmid': param_plasmid_name, 'halftime': t1_2, 'se': se})

            # Plot
            p = (
                    ggplot(data, aes('time', 'norm.bfp'))
                    + geom_point(size=0.4, alpha=0.7, color=color_bfp)
                    + geom_function(
                fun=lambda t: exp_fit_up(t, t1_2),
                color="black"
            )
                    + labs(x="Time (hours)", y="tagBFP (% of final)")
                    + scale_y_continuous(breaks=[0, 0.25, 0.5, 0.75, 1], limits=[-0.15, 1.2])
                    + theme_castuner()
            )
            p.save(
                os.path.join(OUT_PATH, f"KD_{plasmid}_fitting.pdf"),
                width=1.5 * 1.618,
                height=1.5,
                units="in"
            )
        except RuntimeError:
            print(f"Could not fit {plasmid}")

    # Save parameters
    half_time_df = pd.DataFrame(half_times)
    half_time_df.to_csv(
        os.path.join(PARAM_PATH, "half_times_upregulation.csv"),
        index=False
    )
    print("Step 1a: Upregulation fitting complete.")
    print(half_time_df)


if __name__ == "__main__":
    main()
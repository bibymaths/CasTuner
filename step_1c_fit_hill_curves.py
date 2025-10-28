import pandas as pd
import numpy as np
import os
from scipy.optimize import curve_fit
from plotnine import (
    ggplot,
    aes,
    geom_point,
    labs,
    coord_cartesian,
    geom_function,
)
from utils import (
    load_nfc_background,
    load_fcs_data,
    theme_castuner,
    plasmid_map
)


# Define the Hill function
def hill_func(R, K, n):
    """Inhibitory Hill function."""
    return (K ** n) / (K ** n + R ** n)


# dTAG concentration mapping
dtag_map = {
    "1": 0, "2": 0.5, "3": 1, "4": 2, "5": 3, "6": 5,
    "7": 8, "8": 10, "9": 25, "10": 50, "11": 100, "12": 500
}


def load_dose_response_data(path, mBFP_neg, mmCherry_neg):
    """Helper to load one replicate of dose-response data."""
    medianexp = load_fcs_data(path, mBFP_neg, mmCherry_neg)

    # Filename format: plasmid _ guide _ dTAG _ .fcs
    parts = medianexp['filename'].str.split('_', expand=True)
    medianexp['plasmid'] = parts[0]
    medianexp['guide'] = parts[1]

    # Extract dTAG number (e.g., "dTAG1", "dTAG10")
    dtag_part = parts[2].str.extract(r'dTAG(\d+)')[0]
    medianexp['dTAG'] = dtag_part.map(dtag_map).astype(float)
    return medianexp


def main():
    OUT_PATH = "plots"
    PARAM_PATH = "parameters"
    os.makedirs(OUT_PATH, exist_ok=True)
    os.makedirs(PARAM_PATH, exist_ok=True)

    mBFP_neg, mmCherry_neg = load_nfc_background()

    # Load 3 replicates
    rep1 = load_dose_response_data("fcs_files/dose_response_data/R1", mBFP_neg, mmCherry_neg)
    rep2 = load_dose_response_data("fcs_files/dose_response_data/R2", mBFP_neg, mmCherry_neg)
    rep3 = load_dose_response_data("fcs_files/dose_response_data/R3", mBFP_neg, mmCherry_neg)

    rep1['rep'] = 1
    rep2['rep'] = 2
    rep3['rep'] = 3

    d4 = pd.concat([rep1, rep2, rep3])

    # Fix plasmid names (R script used 430, 428, etc.)
    d4['plasmid'] = d4['plasmid'].map({
        "SP430": "430", "SP428": "428", "SP430ABA": "430ABA",
        "SP427": "427", "SP411": "411"
    }).fillna(d4['plasmid'])

    # Calculate mCherry fold-change
    meanNTC = d4[d4['guide'] == 'N'].groupby(['plasmid', 'dTAG'])['PE-A'].mean().reset_index(name='meanNTC')
    d4 = pd.merge(d4, meanNTC, on=['plasmid', 'dTAG'], how='left')
    d4['fc'] = d4['PE-A'] / d4['meanNTC']

    # Min-max scaling of bfp (0 to 1)
    max_bfp = d4[d4['dTAG'] == 0].groupby(['plasmid', 'guide'])['BV421-A'].mean().reset_index(name='max.bfp')
    d4 = pd.merge(d4, max_bfp, on=['plasmid', 'guide'], how='left')

    min_bfp = d4.groupby(['plasmid', 'guide'])['BV421-A'].min().reset_index(name='min.bfp')
    d4 = pd.merge(d4, min_bfp, on=['plasmid', 'guide'], how='left')

    d4['norm.bfp'] = (d4['BV421-A'] - d4['min.bfp']) / (d4['max.bfp'] - d4['min.bfp'])

    # Filter for targeting guides
    d4g = d4[d4['guide'] == 'G'].copy()

    plasmids_data = {
        "SP430": d4g[d4g['plasmid'] == "430"],
        "SP411": d4g[d4g['plasmid'] == "411"],
        "SP427": d4g[d4g['plasmid'] == "427"],
        "SP428": d4g[d4g['plasmid'] == "428"],
        "SP430A": d4g[d4g['plasmid'] == "430ABA"],  # Note: 430ABA maps to SP430A
    }

    parameters = []

    for name, data in plasmids_data.items():
        if data.empty:
            print(f"No data for {name}")
            continue

        R = data['norm.bfp'].dropna()
        y = data['fc'].loc[R.index]

        try:
            popt, _ = curve_fit(hill_func, R, y, p0=[0.1, 1])
            K, n = popt

            parameters.append({'plasmid': name, 'K': K, 'n': n})

            # Plot
            p = (
                    ggplot(data, aes('norm.bfp', 'fc'))
                    + geom_point(size=0.8, alpha=0.4, color="grey")
                    + geom_function(
                fun=lambda r: hill_func(r, K, n),
                color="black"
            )
                    + labs(x="Normalized repressor level", y="Normalized reporter level")
                    + coord_cartesian(ylim=[0, 1.1])
                    + theme_castuner()
            )
            p.save(
                os.path.join(OUT_PATH, f"Hill_{name}.pdf"),
                width=1.5 * 1.618,
                height=1.5,
                units="in"
            )
        except RuntimeError:
            print(f"Could not fit Hill curve for {name}")

    # Save parameters
    param_df = pd.DataFrame(parameters)
    param_df.to_csv(
        os.path.join(PARAM_PATH, "Hill_parameters.csv"),
        index=False
    )
    print("Step 1c: Hill curve fitting complete.")
    print(param_df)


if __name__ == "__main__":
    main()
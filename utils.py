import pandas as pd
import numpy as np
import os
from FlowCytometryTools import FCMeasurement, ThresholdGate
from plotnine import (
    theme_classic,
    theme,
    element_text,
    element_rect,
    element_blank,
    unit
)
from palettable.npg.colors import NpgColors

# Define color palettes
mypal = NpgColors(2).hex_colors
color_bfp = "#4DBBD5FF"
color_mcherry = "#E64B35FF"

# Define plasmid mapping for parameters
# R scripts use 'SP430A' for 'SP430ABA' data
plasmid_map = {
    "SP430": "SP430",
    "SP428": "SP428",
    "SP430ABA": "SP430A",
    "SP427": "SP427",
    "SP411": "SP411",
}


def theme_castuner():
    """A plotnine theme replicating the R project's ggplot2 theme."""
    return theme_classic() + theme(
        legend_text=element_text(size=6),
        panel_border=element_rect(color="black", fill=np.nan, size=0.5),
        axis_line=element_blank(),
        axis_text=element_text(size=6, color="black"),
        axis_text_x=element_text(va="center", color="black"),
        axis_text_y=element_text(va="center", color="black"),
        axis_title=element_text(size=6),
        strip_text=element_text(size=6, color="black"),
        strip_background=element_blank(),
        legend_title=element_blank(),
        figure_size=(1.5 * 1.618, 1.5)  # Use save(..., width=, height=) instead
    )


def load_nfc_background(path="fcs_files/NFC"):
    """Loads non-fluorescent control files and calculates mean background."""
    fcs_files = [os.path.join(path, f) for f in os.listdir(path) if f.endswith('.fcs')]

    # Use the 3rd file for gating, as in the R script
    sample = FCMeasurement(ID="NFC_Gate", datafile=fcs_files[2])

    # Gating
    gate_bound = ThresholdGate(
        {'FSC-A': (0.4e5, 2e5), 'SSC-A': (0.20e5, 1.3e5)},
        region='in',
        name='Boundary'
    )
    gate_singlet = ThresholdGate(
        {'FSC-A': (0, np.inf), 'FSC-H': (0, np.inf)},
        region='in',
        name='Singlet'
    )  # A simple singlet gate placeholder, R uses .singletGate

    gated_sample = sample.gate(gate_bound).gate(gate_singlet)

    # Apply gates to all samples
    medians = []
    for f in fcs_files:
        s = FCMeasurement(ID=os.path.basename(f), datafile=f)
        s_gated = s.gate(gate_bound).gate(gate_singlet)
        medians.append(s_gated.data[['BV421-A', 'PE-A']].median())

    median_df = pd.DataFrame(medians)

    # Get mean of first 3 samples
    mBFP_neg = median_df.iloc[0:3]['BV421-A'].mean()
    mmCherry_neg = median_df.iloc[0:3]['PE-A'].mean()

    return mBFP_neg, mmCherry_neg


def load_fcs_data(path, mBFP_neg, mmCherry_neg):
    """Loads, gates, and processes FCS data from a directory."""
    fcs_files = [os.path.join(path, f) for f in os.listdir(path) if f.endswith('.fcs')]

    if not fcs_files:
        print(f"No FCS files found in {path}")
        return pd.DataFrame()

    # Use the 3rd file for gating
    sample = FCMeasurement(ID="GateSample", datafile=fcs_files[2])
    gate_bound = ThresholdGate(
        {'FSC-A': (0.4e5, 2e5), 'SSC-A': (0.20e5, 1.3e5)},
        region='in',
        name='Boundary'
    )
    gate_singlet = ThresholdGate(
        {'FSC-A': (0, np.inf), 'FSC-H': (0, np.inf)},
        region='in',
        name='Singlet'
    )
    gated_sample = sample.gate(gate_bound).gate(gate_singlet)

    # Process all files
    all_data = []
    for f in fcs_files:
        s = FCMeasurement(ID=os.path.basename(f), datafile=f)
        s_gated = s.gate(gate_bound).gate(gate_singlet)

        # Background subtraction
        data_bgsub = s_gated.data.copy()
        data_bgsub['BV421-A'] -= mBFP_neg
        data_bgsub['PE-A'] -= mmCherry_neg

        median_exp = data_bgsub[['BV421-A', 'PE-A']].median()
        median_exp['filename'] = os.path.basename(f)
        all_data.append(median_exp)

    return pd.DataFrame(all_data)
#!/usr/bin/env python3

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

###############################################################################
# 1. SETUP
###############################################################################

WORKING_DIR = "../output/atacseqintersects/xrseq"
OUTPUT_DIR = os.path.join(WORKING_DIR, "atacseqquartiles")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.chdir(WORKING_DIR)

###############################################################################
# 2. GLOBAL SETTINGS
###############################################################################

NUM_BINS = 200
PEAK_BIN = 100
TIME_POINTS = ["5min", "1h", "8h", "24h", "48h"]
QUARTILES = ["0_25", "25_50", "50_75", "75_100"]
STRAINS = ["wt", "csb", "xpc"]
DAMAGES = ["64", "CPD"]
SMOOTH_WINDOW = 11  # smoothing window size

QUARTILE_COLORS = {
    "0_25": "#c9e4c5",
    "25_50": "#a4cf9c",
    "50_75": "#6ca86e",
    "75_100": "#3f7d4e"
}

plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'axes.facecolor': 'white',
    'axes.edgecolor': 'black',
    'axes.linewidth': 1.5
})

###############################################################################
# 3. HELPER FUNCTIONS
###############################################################################

def load_readcount(strain, timepoint, damage):
    filename = f"{strain}_{timepoint}_{damage}_readCount.txt"
    with open(filename, 'r') as f:
        return float(f.readline().strip())

def load_quartile_bed(strain, timepoint, damage, quartile):
    filename = f"{strain}_{timepoint}_{damage}_{quartile}.bed"
    return pd.read_csv(filename, sep='\t', header=None)

def compute_stats(strain, timepoint, damage, quartile):
    try:
        total_reads = load_readcount(strain, timepoint, damage)
        df = load_quartile_bed(strain, timepoint, damage, quartile)
    except FileNotFoundError:
        raise

    df['RPM'] = df[5] / (total_reads / 1e6)
    rpm_by_bin = df.groupby(4)['RPM'].agg(['mean', 'count', 'std']).reindex(range(1, NUM_BINS + 1), fill_value=0)
    rpm_by_bin['sem'] = rpm_by_bin['std'] / rpm_by_bin['count'].replace(0, np.nan).pow(0.5)
    rpm_by_bin['ci95'] = 1.96 * rpm_by_bin['sem']
    rpm_by_bin = rpm_by_bin.fillna(0)

    rpm_by_bin['mean'] = rpm_by_bin['mean'].rolling(window=SMOOTH_WINDOW, center=True, min_periods=1).mean()
    rpm_by_bin['ci95'] = rpm_by_bin['ci95'].rolling(window=SMOOTH_WINDOW, center=True, min_periods=1).mean()
    return rpm_by_bin

###############################################################################
# 4. PLOTTING FUNCTION
###############################################################################

def plot_enrichment(strain, damage):
    fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(20, 5), sharey=True)
    fig.suptitle(f"{strain.upper()} {damage} XR-Seq Quartile Enrichment (±95% CI)", fontsize=15)

    for ax, timepoint in zip(axes, TIME_POINTS):
        for quartile in QUARTILES:
            try:
                stats = compute_stats(strain, timepoint, damage, quartile)
                x = np.arange(1, NUM_BINS + 1)
                mean = stats['mean'].values
                ci = stats['ci95'].values

                ax.plot(x, mean, label=quartile.replace("_", "-") + "%", color=QUARTILE_COLORS[quartile], linewidth=2)
                ax.fill_between(x, mean - ci, mean + ci, color=QUARTILE_COLORS[quartile], alpha=0.2)
            except FileNotFoundError:
                print(f"Missing: {strain}_{timepoint}_{damage}_merged_{quartile}.bed")
                continue

        ax.set_title(timepoint)
        ax.set_xlim(1, NUM_BINS)
        ax.set_ylim(0.15, 2.25)
        ax.set_xticks([1, PEAK_BIN, NUM_BINS])
        ax.set_xticklabels(["+1kb", "Peak", "-1kb"])
        ax.set_xlabel("Peak-centered bins", fontsize=12)
        ax.label_outer()

    axes[0].set_ylabel("RPM (Reads Per Million)", fontsize=12)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=4, frameon=False, fontsize=11)

    plt.tight_layout()
    plt.subplots_adjust(top=0.85, bottom=0.15)
    outname = os.path.join(OUTPUT_DIR, f"{strain}_{damage}_atacseq_quartiles.pdf")
    plt.savefig(outname, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outname}")

###############################################################################
# 5. MAIN
###############################################################################

if __name__ == "__main__":
    for strain in STRAINS:
        for damage in DAMAGES:
            plot_enrichment(strain, damage)

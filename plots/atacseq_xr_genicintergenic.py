#!/usr/bin/env python3

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

###############################################################################
# 1. SETUP
###############################################################################

WORKING_DIR = "../output/atacseqintersects/xrseq"
OUTPUT_DIR = os.path.join(WORKING_DIR, "atacseqgenicintergenic)
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.chdir(WORKING_DIR)

###############################################################################
# 2. GLOBAL SETTINGS
###############################################################################

NUM_BINS = 200
PEAK_BIN = 100
TIME_POINTS = ["5min", "1h", "8h", "24h", "48h"]
STRAINS = ["wt", "csb", "xpc"]
DAMAGES = ["64", "CPD"]
CATEGORIES = ["genic", "intergenic"]
CATEGORY_COLORS = {"genic": "black", "intergenic": "gray"}
SMOOTH_WINDOW = 11

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
    fname = f"{strain}_{timepoint}_{damage}_merged_simulation_readCount.txt"
    if not os.path.exists(fname):
        print(f"Missing read count: {fname}")
        return None
    with open(fname) as f:
        return float(f.readline().strip())

def load_bed_file(strain, timepoint, damage, category):
    suffix = "promotergenic" if category == "genic" else "intergenic"
    fname = f"{strain}_{timepoint}_{damage}_merged_simulation_{suffix}.bed"
    if not os.path.exists(fname):
        print(f"Missing bed: {fname}")
        return pd.DataFrame()
    return pd.read_csv(fname, sep="\t", header=None)

def compute_smoothed_rpm_and_ci(strain, timepoint, damage, category):
    total_reads = load_readcount(strain, timepoint, damage)
    df = load_bed_file(strain, timepoint, damage, category)

    if df.empty or total_reads is None:
        return np.zeros(NUM_BINS), np.zeros(NUM_BINS)

    df["RPM"] = df[5] / (total_reads / 1e6)

    # Compute RPM values for each bin
    grouped = df.groupby(4)["RPM"]
    mean = grouped.mean().reindex(range(1, NUM_BINS + 1), fill_value=0)
    std = grouped.std().reindex(range(1, NUM_BINS + 1), fill_value=0)
    count = grouped.count().reindex(range(1, NUM_BINS + 1), fill_value=1)
    sem = std / np.sqrt(count)
    ci95 = 1.96 * sem

    # Smooth both mean and CI
    mean_smooth = mean.rolling(window=SMOOTH_WINDOW, center=True, min_periods=1).mean()
    ci_smooth = ci95.rolling(window=SMOOTH_WINDOW, center=True, min_periods=1).mean()

    return mean_smooth.values, ci_smooth.values

###############################################################################
# 4. PLOTTING FUNCTION
###############################################################################

def plot_genic_vs_intergenic(strain, damage, y_limit):
    fig, axes = plt.subplots(1, 5, figsize=(20, 5), sharey=(strain in ["wt", "csb"]))
    fig.suptitle(f"{strain.upper()} {damage} ATAC Genic vs Intergenic (Smoothed ±95% CI)", fontsize=16)

    for ax, tp in zip(axes, TIME_POINTS):
        for category in CATEGORIES:
            mean, ci = compute_smoothed_rpm_and_ci(strain, tp, damage, category)
            x = np.arange(1, NUM_BINS + 1)
            ax.plot(x, mean, label=category, color=CATEGORY_COLORS[category], linewidth=2)
            ax.fill_between(x, mean - ci, mean + ci, color=CATEGORY_COLORS[category], alpha=0.2)

        ax.set_title(tp)
        ax.set_xlim(1, NUM_BINS)
        ax.set_xticks([1, PEAK_BIN, NUM_BINS])
        ax.set_xticklabels(['-1kb', 'Peak', '+1kb'])
        if y_limit and strain in y_limit:
            ax.set_ylim(y_limit[strain])
        ax.label_outer()

    axes[0].set_ylabel("RPM (Smoothed ± CI)", fontsize=12)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    outpath = os.path.join(OUTPUT_DIR, f"{strain}_{damage}_atacseqgenicintergenic.pdf")
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outpath}")

###############################################################################
# 5. MAIN EXECUTION
###############################################################################

FIXED_Y_LIMITS = {"wt": (0.2, 1.5), "csb": (0.2, 1.5), "xpc": (0.2, 1.5)}

if __name__ == "__main__":
    for strain in STRAINS:
        for damage in DAMAGES:
            plot_genic_vs_intergenic(strain, damage, FIXED_Y_LIMITS)

#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

###############################################################################
# 1. SETUP
###############################################################################

WORKING_DIR = "../output/atacseqintersects/xrseq"
OUTPUT_DIR  = os.path.join(WORKING_DIR, "atacseqplots")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.chdir(WORKING_DIR)

###############################################################################
# 2. SETTINGS
###############################################################################

NUM_BINS   = 200
PEAK_BIN   = 100
TIME_POINTS = ["5min", "1h", "8h", "24h", "48h"]
STRAINS = ["wt", "csb", "xpc"]
DAMAGES = ["CPD", "64"]
SMOOTH_WINDOW = 11

# Custom color mapping
color_map = {
    "5min": "#A3B57A",   # olive green
    "1h":   "#D4B46D",   # gold
    "8h":   "#BF6B57",   # red-brown
    "24h":  "#8A5D71",   # plum
    "48h":  "#4F4D63"    # slate gray
}

###############################################################################
# 3. FUNCTIONS
###############################################################################

def load_readcount(strain, timepoint, damage):
    fname = f"{strain}_{timepoint}_{damage}_readCount.txt"
    with open(fname, 'r') as f:
        return float(f.readline().strip())

def load_bed_file(strain, timepoint, damage):
    fname = f"{strain}_{timepoint}_{damage}_allpeaks.bed"
    df = pd.read_csv(fname, sep='\t', header=None)
    return df

def compute_smoothed_stats(strain, timepoint, damage):
    df = load_bed_file(strain, timepoint, damage)
    total_reads = load_readcount(strain, timepoint, damage)
    df["RPM"] = df[5] / (total_reads / 1e6)

    grouped = df.groupby(df[4])["RPM"]
    stats = grouped.agg(["mean", "count", "std"]).reindex(range(1, NUM_BINS + 1), fill_value=0)
    stats["sem"] = stats["std"] / stats["count"].replace(0, np.nan).pow(0.5)
    stats["ci95"] = 1.96 * stats["sem"]
    stats = stats.fillna(0)

    mean_smooth = stats["mean"].rolling(window=SMOOTH_WINDOW, center=True, min_periods=1).mean()
    ci_smooth = stats["ci95"].rolling(window=SMOOTH_WINDOW, center=True, min_periods=1).mean()

    return mean_smooth.values, ci_smooth.values

###############################################################################
# 4. PLOTTING
###############################################################################

def plot_condition(strain, damage):
    plt.figure(figsize=(6, 5))
    ax = plt.gca()

    for tp in TIME_POINTS:
        try:
            mean, ci = compute_smoothed_stats(strain, tp, damage)
            x = np.arange(1, NUM_BINS + 1)

            ax.plot(x, mean, label=tp, color=color_map[tp], linewidth=2)
            ax.fill_between(x, mean - ci, mean + ci, color=color_map[tp], alpha=0.2)
        except FileNotFoundError:
            print(f"Missing file for {strain}_{tp}_{damage}")
            continue

    ax.set_xlim(1, NUM_BINS)
    ax.set_xticks([1, PEAK_BIN, NUM_BINS])
    ax.set_xticklabels(["-1kb", "peak", "+1kb"])
    ax.set_ylim(0.2, 1.3)
    ax.set_xlabel("Peak-centered bins", fontsize=12)
    ax.set_ylabel("RPM (Reads Per Million) ± CI", fontsize=12)

    ax.set_title(f"{strain.upper()} {damage} XR-Seq over ATAC-Seq peaks", fontsize=14, pad=10)
    ax.legend(loc='upper center',
              bbox_to_anchor=(0.5, -0.15),
              ncol=5,
              fontsize=10,
              frameon=False)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)

    outname = f"{strain}_{damage}_atacseqallpeaks.pdf"
    outpath = os.path.join(OUTPUT_DIR, outname)
    plt.savefig(outpath, bbox_inches='tight')
    plt.close()
    print(f"Saved: {outpath}")

###############################################################################
# 5. RUN
###############################################################################

if __name__ == "__main__":
    for strain in STRAINS:
        for damage in DAMAGES:
            plot_condition(strain, damage)

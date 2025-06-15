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

WORKING_DIR = "../output/atacseqintersects/damageseq"
OUTPUT_DIR = os.path.join(WORKING_DIR, "atacseqallpeaks")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.chdir(WORKING_DIR)

###############################################################################
# 2. GLOBAL SETTINGS
###############################################################################

NUM_BINS = 200
PEAK_BIN = 100
TIME_POINTS = ["nak", "0h", "8h", "24h", "48h"]
FINAL_PLOT_POINTS = ["0h", "8h", "24h", "48h"]
REFERENCE_TIMEPOINT = "noUV_dip"
EPIGENOMES = ["allpeaks"]
STRAINS = ["wt"]
DAMAGES = ["64", "CPD"]
colors = {
    "0h": "#85BFC4",
    "8h": "#C2665A",
    "24h": "#94617D",
    "48h": "#4B4B66"
}

SMOOTH_WINDOW = 11

###############################################################################
# 3. FUNCTIONS
###############################################################################

def load_readcount(strain, timepoint, damage):
    fname = f"DS_{strain}_damage{damage}_{timepoint}_readCount.txt"
    with open(fname) as f:
        return float(f.readline().strip())

def load_epigenome_bed(strain, timepoint, damage, mark):
    fname = f"DS_{strain}_damage{damage}_{timepoint}_{mark}.bed"
    df = pd.read_csv(fname, sep="\t", header=None)
    return df

def get_rpm_df(strain, timepoint, damage, mark):
    total_reads = load_readcount(strain, timepoint, damage)
    df = load_epigenome_bed(strain, timepoint, damage, mark)
    df["RPM"] = df[5] / (total_reads / 1e6)
    return df[[4, "RPM"]].rename(columns={4: "bin"})

def compute_two_step_norm_and_stats(strain, damage, mark):
    # Step 1: load input (noUV_dip)
    input_df = get_rpm_df(strain, REFERENCE_TIMEPOINT, damage, mark)
    input_means = input_df.groupby("bin")["RPM"].mean().reindex(range(1, NUM_BINS + 1), fill_value=0) + 1e-6

    # Step 2: load nak and normalize to input
    nak_df = get_rpm_df(strain, "nak", damage, mark)
    nak_df["RPM_norm1"] = nak_df["RPM"] / nak_df["bin"].map(input_means)
    mean_nak_norm1 = nak_df.groupby("bin")["RPM_norm1"].mean().reindex(range(1, NUM_BINS + 1), fill_value=1)

    stats_dict = {}

    for tp in FINAL_PLOT_POINTS:
        df = get_rpm_df(strain, tp, damage, mark)
        df["norm1"] = df["RPM"] / df["bin"].map(input_means)
        df["final"] = df["norm1"] / df["bin"].map(mean_nak_norm1)

        grouped = df.groupby("bin")["final"]
        stats = grouped.agg(["mean", "count", "std"]).reindex(range(1, NUM_BINS + 1), fill_value=0)
        stats["sem"] = stats["std"] / stats["count"].replace(0, np.nan).pow(0.5)
        stats["ci95"] = 1.96 * stats["sem"]
        stats_dict[tp] = stats.fillna(0)

    return stats_dict

###############################################################################
# 4. PLOTTING
###############################################################################

def plot_final_normalized(strain, damage, mark):
    print(f"Plotting: {strain}, {damage}, {mark}")
    stats_dict = compute_two_step_norm_and_stats(strain, damage, mark)

    plt.figure(figsize=(6, 5))
    ax = plt.gca()

    for i, tp in enumerate(FINAL_PLOT_POINTS):
        stats = stats_dict[tp]
        x = np.arange(1, NUM_BINS + 1)

        mean_smooth = stats["mean"].rolling(window=SMOOTH_WINDOW, center=True, min_periods=1).mean()
        ci_smooth = stats["ci95"].rolling(window=SMOOTH_WINDOW, center=True, min_periods=1).mean()

        ax.plot(x, mean_smooth, label=tp, color=colors[tp], linewidth=2)
        ax.fill_between(x, mean_smooth - ci_smooth, mean_smooth + ci_smooth, color=colors[tp], alpha=0.2)

    ax.set_xlim(1, NUM_BINS)
    ax.set_xticks([1, PEAK_BIN, NUM_BINS])
    ax.set_xticklabels(["+1kb", "peak", "-1kb"])
    ax.set_ylim(0.82, 1.05)
    ax.set_xlabel("Peak-centered bins", fontsize=12)
    ax.set_ylabel("Final normalized signal (TP / noUV_dip → then / mean(nak))", fontsize=11)

    damage_label = "6-4PP" if damage == "64" else "CPD"
    ax.set_title(f"{strain.upper()} {damage_label} over {mark}", fontsize=13)

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=10, frameon=False)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)

    outname = f"{strain}_{damage}_{mark}_atacseqallpeaks.pdf"
    plt.savefig(os.path.join(OUTPUT_DIR, outname), bbox_inches="tight")
    plt.close()
    print(f"Saved: {outname}")

###############################################################################
# 5. RUN
###############################################################################

if __name__ == "__main__":
    for strain in STRAINS:
        for damage in DAMAGES:
            for mark in EPIGENOMES:
                plot_final_normalized(strain, damage, mark)

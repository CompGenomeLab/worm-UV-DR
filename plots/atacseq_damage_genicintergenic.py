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
OUTPUT_DIR = os.path.join(WORKING_DIR, "atacseqgenicintergenic")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.chdir(WORKING_DIR)

###############################################################################
# 2. GLOBAL SETTINGS
###############################################################################

NUM_BINS = 200
PEAK_BIN = 100
STRAIN = "wt"
DAMAGES = ["64", "CPD"]
TIME_POINTS = ["0h", "8h", "24h", "48h"]
REFERENCE_TIMEPOINT = "noUV_dip"
NAK_TIMEPOINT = "nak"
CATEGORIES = ["promotergenic", "intergenic"]
CATEGORY_LABELS = {"promotergenic": "genic", "intergenic": "intergenic"}
CATEGORY_COLORS = {"promotergenic": "black", "intergenic": "gray"}
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
# 3. FUNCTIONS
###############################################################################

def load_readcount(strain, timepoint, damage):
    fname = f"DS_{strain}_damage{damage}_{timepoint}_readCount.txt"
    with open(fname) as f:
        return float(f.readline().strip())

def load_category_bed(strain, timepoint, damage, category):
    fname = f"DS_{strain}_damage{damage}_{timepoint}_{category}.bed"
    return pd.read_csv(fname, sep="\t", header=None)

def get_rpm_df(strain, timepoint, damage, category):
    total_reads = load_readcount(strain, timepoint, damage)
    df = load_category_bed(strain, timepoint, damage, category)
    df["RPM"] = df[5] / (total_reads / 1e6)
    return df[[4, "RPM"]].rename(columns={4: "bin"})

def compute_norm_stats(strain, damage, category):
    input_df = get_rpm_df(strain, REFERENCE_TIMEPOINT, damage, category)
    input_means = input_df.groupby("bin")["RPM"].mean().reindex(range(1, NUM_BINS + 1), fill_value=0) + 1e-6

    nak_df = get_rpm_df(strain, NAK_TIMEPOINT, damage, category)
    nak_df["RPM_norm1"] = nak_df["RPM"] / nak_df["bin"].map(input_means)
    mean_nak_norm1 = nak_df.groupby("bin")["RPM_norm1"].mean().reindex(range(1, NUM_BINS + 1), fill_value=1)

    stats_dict = {}
    for tp in TIME_POINTS:
        df = get_rpm_df(strain, tp, damage, category)
        df["norm1"] = df["RPM"] / df["bin"].map(input_means)
        df["final"] = df["norm1"] / df["bin"].map(mean_nak_norm1)

        grouped = df.groupby("bin")["final"]
        stats = grouped.agg(["mean", "count", "std"]).reindex(range(1, NUM_BINS + 1), fill_value=0)
        stats["sem"] = stats["std"] / stats["count"].replace(0, np.nan).pow(0.5)
        stats["ci95"] = 1.96 * stats["sem"]
        stats_dict[tp] = stats.fillna(0)

    return stats_dict

###############################################################################
# 4. PLOTTING FUNCTION
###############################################################################

def plot_genic_intergenic(damage):
    fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(20, 5), sharey=True)
    fig.suptitle(f"{STRAIN.upper()} {damage} Damage-Seq Genic vs. Intergenic", fontsize=14, y=1.05)

    stats_by_category = {cat: compute_norm_stats(STRAIN, damage, cat) for cat in CATEGORIES}
    x = np.arange(1, NUM_BINS + 1)

    for ax, tp in zip(axes, TIME_POINTS):
        for cat in CATEGORIES:
            stats = stats_by_category[cat][tp]
            mean_smooth = stats["mean"].rolling(window=SMOOTH_WINDOW, center=True, min_periods=1).mean()
            ci_smooth = stats["ci95"].rolling(window=SMOOTH_WINDOW, center=True, min_periods=1).mean()

            ax.plot(x, mean_smooth, label=CATEGORY_LABELS[cat],
                    color=CATEGORY_COLORS[cat], linewidth=2)
            ax.fill_between(x, mean_smooth - ci_smooth, mean_smooth + ci_smooth,
                            color=CATEGORY_COLORS[cat], alpha=0.2)

        ax.set_title(tp, fontsize=13)
        ax.set_xlim(1, NUM_BINS)
        ax.set_ylim(0.80, 1.1)
        ax.set_xticks([1, PEAK_BIN, NUM_BINS])
        ax.set_xticklabels(["+1kb", "Peak", "-1kb"])
        ax.set_xlabel("Peak-centered bins", fontsize=12)
        ax.label_outer()

    axes[0].set_ylabel("Normalized signal ± 95% CI", fontsize=12)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2, fontsize=11, frameon=False)

    plt.tight_layout()
    plt.subplots_adjust(top=0.88, bottom=0.2)
    outname = os.path.join(OUTPUT_DIR, f"{STRAIN}_{damage}_atacseqgenicintergenic.pdf")
    plt.savefig(outname, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outname}")

###############################################################################
# 5. EXECUTE
###############################################################################

if __name__ == "__main__":
    for damage in DAMAGES:
        plot_genic_intergenic(damage)

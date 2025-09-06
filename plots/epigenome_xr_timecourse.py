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

WORKING_DIR = "../output/epigenome/intersect/xrseq/5kb/stranded"
OUTPUT_DIR = os.path.join(WORKING_DIR, "./strandedepigenomeplots")
SIM_DIR = "../output/epigenome/intersect/xrseq/5kb/stranded/simulation"

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.chdir(WORKING_DIR)

###############################################################################
# 2. GLOBAL SETTINGS
###############################################################################

NUM_BINS = 500
PEAK_BIN = 250
TIME_POINTS = ["5min", "1h", "8h", "24h", "48h"]
EPIGENOMES = ["H3K4me3", "H3K4me1", "H3K27me3", "H3K36me3"]
STRAINS = ["wt", "csb", "xpc"]
DAMAGES = ["64", "CPD"]

color_map = {
    "5min": "#A3B57A",
    "1h":   "#D4B46D",
    "8h":   "#BF6B57",
    "24h":  "#8A5D71",
    "48h":  "#4F4D63"
}

###############################################################################
# 3. HELPER FUNCTIONS (STRANDED)
###############################################################################

def load_stranded_readcount(strain, timepoint, damage, sim=False):
    prefix = f"{strain}_{timepoint}_{damage}_merged_"
    suffix = "simulation_" if sim else ""

    f_plus = os.path.join(SIM_DIR if sim else WORKING_DIR, f"{prefix}{suffix}plus_readCount.txt")
    f_minus = os.path.join(SIM_DIR if sim else WORKING_DIR, f"{prefix}{suffix}minus_readCount.txt")

    with open(f_plus) as f1, open(f_minus) as f2:
        return float(f1.readline().strip()) + float(f2.readline().strip())

def load_stranded_bed_df(strain, timepoint, damage, mark, sim=False):
    prefix = f"{strain}_{timepoint}_{damage}_merged_"
    suffix = "simulation_" if sim else ""

    folder = SIM_DIR if sim else WORKING_DIR
    f_plus = os.path.join(folder, f"{prefix}{suffix}plus_{mark}_5kb.bed")
    f_minus = os.path.join(folder, f"{prefix}{suffix}minus_{mark}_5kb.bed")

    df_plus = pd.read_csv(f_plus, sep="\t", header=None)
    df_minus = pd.read_csv(f_minus, sep="\t", header=None)

    # Flip bins for plus strand
    df_plus[4] = NUM_BINS + 1 - df_plus[4]

    return pd.concat([df_plus, df_minus], ignore_index=True)

def calculate_rpm_stats(df, total_reads):
    df["RPM"] = df[5] / (total_reads / 1e6)
    grouped = df.groupby(df[4])
    stats = grouped["RPM"].agg(["mean", "count", "std"])
    stats = stats.reindex(range(1, NUM_BINS + 1), fill_value=np.nan)
    return stats

def get_normalized_stats(strain, timepoint, damage, mark):
    # --- Experimental ---
    df_exp = load_stranded_bed_df(strain, timepoint, damage, mark)
    total_exp_reads = load_stranded_readcount(strain, timepoint, damage)
    exp_stats = calculate_rpm_stats(df_exp, total_exp_reads)

    # --- Simulation ---
    df_sim = load_stranded_bed_df(strain, timepoint, damage, mark, sim=True)
    total_sim_reads = load_stranded_readcount(strain, timepoint, damage, sim=True)
    sim_stats = calculate_rpm_stats(df_sim, total_sim_reads)

    # --- Normalize ---
    norm_mean = exp_stats["mean"] / sim_stats["mean"]
    sem = exp_stats["std"] / np.sqrt(exp_stats["count"])
    ci95 = 1.96 * sem / sim_stats["mean"]

    norm_df = pd.DataFrame({
        "norm_mean": norm_mean,
        "ci95": ci95
    }).fillna(0)

    return norm_df

###############################################################################
# 4. PLOTTING FUNCTION
###############################################################################

def plot_condition(strain, damage, mark):
    plt.figure(figsize=(6, 5))
    ax = plt.gca()

    # Y-axis limits
    if mark == "H3K4me1":
        if strain in ["wt", "csb"]:
            ax.set_ylim(-0.05, 0.85)
        elif strain == "xpc":
            ax.set_ylim(-0.05, 1.25)
    elif mark == "H3K4me3":
        if strain in ["wt", "csb"]:
            ax.set_ylim(-0.05, 1.3)
        elif strain == "xpc":
            ax.set_ylim(-0.05, 1.6)
    elif mark == "H3K27me3":
        ax.set_ylim(-1.4, 0)
    elif mark == "H3K36me3":
        if strain in ["wt", "csb"]:
            ax.set_ylim(-0.5, 0.3)
        elif strain == "xpc":
            ax.set_ylim(0, 1.3)

    for tp in TIME_POINTS:
        try:
            stats = get_normalized_stats(strain, tp, damage, mark)

            mean_smooth = stats["norm_mean"].rolling(window=11, center=True, min_periods=1).mean()
            ci_smooth = stats["ci95"].rolling(window=11, center=True, min_periods=1).mean()

            mean_log2 = np.log2(mean_smooth.replace(0, np.nan))
            ci_log2 = ci_smooth / (mean_smooth * np.log(2))
            ci_log2 = ci_log2.replace([np.inf, -np.inf], np.nan).fillna(0)

            x = np.arange(1, NUM_BINS + 1)
            ax.plot(x, mean_log2, label=tp, color=color_map[tp], linewidth=2)
            ax.fill_between(x, mean_log2 - ci_log2, mean_log2 + ci_log2,
                            color=color_map[tp], alpha=0.2)

        except FileNotFoundError as e:
            print(f"Missing file for {strain} {tp} {damage} {mark}: {e}")
            continue
        except ZeroDivisionError as e:
            print(f"Division error for {strain} {tp} {damage} {mark}: {e}")
            continue

    ax.set_xlim(1, NUM_BINS)
    ax.axvline(x=250.5, color='gray', linestyle=':', linewidth=1)  # Faint dotted line at real peak
    ax.set_xticks([1, 250.5, NUM_BINS])
    ax.set_xticklabels(["-2.5kb", "peak center", "+2.5kb"])
    ax.set_xlabel("Peak-centered bins", fontsize=12)
    ax.set_ylabel("log2(Normalized RPM)", fontsize=12)

    damage_label = "6-4PP" if damage == "64" else "CPD"
    ax.set_title(f"{strain.upper()} {damage_label} over {mark} (stranded simnorm)", fontsize=14, pad=10)

    ax.legend(loc='upper center',
              bbox_to_anchor=(0.5, -0.15),
              ncol=5,
              fontsize=10,
              frameon=False)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)

    outname = f"{strain}_{damage}_{mark}_simnorm_stranded_log2.pdf"
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
            for mark in EPIGENOMES:
                plot_condition(strain, damage, mark)

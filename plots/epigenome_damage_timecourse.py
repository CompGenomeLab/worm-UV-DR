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

WORKING_DIR = "../output/epigenome/intersect/damageseq/5kb/stranded"
OUTPUT_DIR = os.path.join(WORKING_DIR, "/strandedepigenomeplots")
SIM_DIR = "../output/epigenome/intersect/damageseq/5kb/stranded/simulation"

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.chdir(WORKING_DIR)

###############################################################################
# 2. GLOBAL SETTINGS
###############################################################################

NUM_BINS = 500
PEAK_BIN = 250
TIME_POINTS = ["noUV_dip", "nak", "0h", "8h", "24h", "48h"]
MARKS = ["H3K4me3", "H3K4me1", "H3K27me3", "H3K36me3"]
DAMAGES = ["64", "CPD"]
STRAIN = "wt"

color_map = {
    "0h": "#85BFC4",
    "8h": "#C2665A",
    "24h": "#94617D",
    "48h": "#4B4B66"
}

###############################################################################
# 3. FUNCTIONS
###############################################################################

def load_stranded_readcount(tp, damage, sim=False):
    prefix = "simulation_" if sim else ""
    base = f"DS_{STRAIN}_damage{damage}_{tp}_merged_"
    folder = SIM_DIR if sim else WORKING_DIR
    f_plus = os.path.join(folder, f"{base}{prefix}plus_readCount.txt")
    f_minus = os.path.join(folder, f"{base}{prefix}minus_readCount.txt")
    with open(f_plus) as f1, open(f_minus) as f2:
        return float(f1.readline().strip()) + float(f2.readline().strip())

def load_stranded_bed(tp, damage, mark, sim=False):
    prefix = "simulation_" if sim else ""
    base = f"DS_{STRAIN}_damage{damage}_{tp}_merged_"
    folder = SIM_DIR if sim else WORKING_DIR
    f_plus = os.path.join(folder, f"{base}{prefix}plus_{mark}_5kb.bed")
    f_minus = os.path.join(folder, f"{base}{prefix}minus_{mark}_5kb.bed")
    df_plus = pd.read_csv(f_plus, sep="\t", header=None)
    df_minus = pd.read_csv(f_minus, sep="\t", header=None)

    # Flip bin numbers for minus strand so that both strands are 5' to 3'
    df_minus[4] = NUM_BINS + 1 - df_minus[4]

    return pd.concat([df_plus, df_minus], ignore_index=True)

def compute_mean_rpm_by_bin(df, total_reads):
    df["RPM"] = df[5] / (total_reads / 1e6)
    return df.groupby(df[4])["RPM"].mean().reindex(range(1, NUM_BINS + 1), fill_value=np.nan)

def get_normalized_profiles(damage, mark, sim=False):
    profiles = {}
    for tp in TIME_POINTS:
        try:
            total_reads = load_stranded_readcount(tp, damage, sim=sim)
            df = load_stranded_bed(tp, damage, mark, sim=sim)
            rpm = compute_mean_rpm_by_bin(df, total_reads)
            profiles[tp] = rpm
        except Exception as e:
            print(f"Missing {tp} {'simulation' if sim else 'original'}: {e}")
            profiles[tp] = pd.Series(np.nan, index=range(1, NUM_BINS + 1))
    return profiles

def get_shared_ylim(mark):
    all_values = []

    for damage in DAMAGES:
        orig_profiles = get_normalized_profiles(damage, mark, sim=False)
        sim_profiles = get_normalized_profiles(damage, mark, sim=True)

        norm_orig = {tp: orig_profiles[tp] / orig_profiles["noUV_dip"] for tp in TIME_POINTS}
        norm_sim = {tp: sim_profiles[tp] / sim_profiles["noUV_dip"] for tp in TIME_POINTS}

        nak_mean_orig = norm_orig["nak"].mean()
        nak_mean_sim = norm_sim["nak"].mean()

        for tp in ["0h", "8h", "24h", "48h"]:
            try:
                signal = norm_orig[tp] / nak_mean_orig
                sim = norm_sim[tp] / nak_mean_sim
                ratio = signal / sim
                smooth = ratio.rolling(window=11, center=True, min_periods=1).mean()
                log2_smooth = np.log2(smooth.replace(0, np.nan))
                all_values.append(log2_smooth)
            except:
                continue

    all_concat = pd.concat(all_values, axis=0)
    return all_concat.min(), all_concat.max()

###############################################################################
# 4. PLOTTING FUNCTION
###############################################################################

def plot_damage_condition(damage, mark, ylim):
    orig_profiles = get_normalized_profiles(damage, mark, sim=False)
    sim_profiles = get_normalized_profiles(damage, mark, sim=True)

    norm_orig = {tp: orig_profiles[tp] / orig_profiles["noUV_dip"] for tp in TIME_POINTS}
    norm_sim = {tp: sim_profiles[tp] / sim_profiles["noUV_dip"] for tp in TIME_POINTS}

    nak_mean_orig = norm_orig["nak"].mean()
    nak_mean_sim = norm_sim["nak"].mean()

    plt.figure(figsize=(6, 5))
    ax = plt.gca()
    ax.set_ylim(*ylim)

    for tp in ["0h", "8h", "24h", "48h"]:
        try:
            signal = norm_orig[tp] / nak_mean_orig
            sim = norm_sim[tp] / nak_mean_sim
            ratio = signal / sim

            signal_sd = (orig_profiles[tp] / orig_profiles["noUV_dip"]).std()
            sem = signal_sd / np.sqrt(len(signal.dropna()))
            sim_mean = sim.mean()
            ci_linear = 1.96 * sem / sim_mean

            smooth = ratio.rolling(window=11, center=True, min_periods=1).mean()
            log2_smooth = np.log2(smooth.replace(0, np.nan))
            ci_log2 = ci_linear / (smooth * np.log(2))
            ci_log2 = ci_log2.replace([np.inf, -np.inf], np.nan).fillna(0)

            x = np.arange(1, NUM_BINS + 1)
            ax.plot(x, log2_smooth, label=tp, color=color_map[tp], linewidth=2)
            ax.fill_between(x, log2_smooth - ci_log2, log2_smooth + ci_log2,
                            color=color_map[tp], alpha=0.2)
        except Exception as e:
            print(f"Plot error at {tp}: {e}")
            continue

    ax.set_xlim(1, NUM_BINS)
    ax.axvline(x=250.5, color='gray', linestyle=':', linewidth=1)  # Faint dotted line at real peak
    ax.set_xticks([1, 250.5, NUM_BINS])
    ax.set_xticklabels(["-2.5kb", "peak center", "+2.5kb"])
    ax.set_xlabel("Peak-centered bins", fontsize=12)
    ax.set_ylabel("log2(2-step Normalized RPM)", fontsize=12)
    ax.set_title(f"DS_wt damage{damage} over {mark} (stranded)", fontsize=14, pad=10)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=6, fontsize=10, frameon=False)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)
    outname = f"DS_wt_damage{damage}_{mark}_simnorm_stranded_log2.pdf"
    plt.savefig(os.path.join(OUTPUT_DIR, outname), bbox_inches='tight')
    plt.close()
    print(f"Saved: {outname}")

###############################################################################
# 5. RUN
###############################################################################

if __name__ == "__main__":
    ylimits_per_mark = {mark: get_shared_ylim(mark) for mark in MARKS}

    for damage in DAMAGES:
        for mark in MARKS:
            plot_damage_condition(damage, mark, ylimits_per_mark[mark])

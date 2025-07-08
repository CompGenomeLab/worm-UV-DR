import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import periodicity_helper_func as hf
import snr

output_folder="" # Directory of the output (should be same as the output directory in run_rep_dam.py)

def smooth_rolling(df, column, window_size=3):
    df = df[["quantile_label", "window", column]].copy()
    df[column] = df.groupby('quantile_label')[column].transform(lambda x: x.rolling(window=window_size, center=True, min_periods=1).mean())
    return df

def create_tt_freq(region_size=100):

    expanded_regions = hf.prepare_atac_seq(quant_num=4, file="resources/atac/c_elegans_ce11/c_elegans_ce11_atac.nucpos_org.bed", region_size=region_size, interval_length=2)
    genome_fasta_path = 'resources/ce11.fa'

    # Calculate TT frequencies
    print("Calculating TT frequencies...")
    tt_results = hf.calculate_dinucleotide_frequency_celegans(expanded_regions, genome_fasta_path, dinucleotide="TT")
    aa_results = hf.calculate_dinucleotide_frequency_celegans(expanded_regions, genome_fasta_path, dinucleotide="AA")
    aa_results["window"] = aa_results["window"] * -1

    tt_results = pd.merge(tt_results, aa_results, on=["region_id","chrom", "start", "end", "quantile", "window", "sequence_length", "sequence"], how="outer")

    # Save results
    # tt_results.to_csv('celegans_tt_frequencies.csv', index=False)

    # Print summary statistics
    print("\nTT Frequency Summary:")
    print(tt_results['TT_frequency'].describe())

    print(f"\nTotal regions processed: {len(tt_results)}")
    print(f"Average TT frequency: {tt_results['TT_frequency'].mean():.4f}")
    print(f"Total TT dinucleotides found: {tt_results['TT_count'].sum()}")

    # Print summary statistics
    print("\nAA Frequency Summary:")
    print(tt_results['AA_frequency'].describe())

    print(f"\nTotal regions processed: {len(tt_results)}")
    print(f"Average AA frequency: {tt_results['AA_frequency'].mean():.4f}")
    print(f"Total AA dinucleotides found: {tt_results['AA_count'].sum()}")

    tt_results_agg = tt_results[["quantile", "window", "TT_count", "AA_count"]].groupby(["quantile", "window"]).agg("mean").reset_index()
    tt_results_agg["window"] = tt_results_agg["window"] * 2
    tt_results_agg["total_count_TT"] = tt_results_agg["TT_count"] + tt_results_agg["AA_count"]

    return tt_results_agg

if os.path.exists(f"{output_folder}/dyad_mapped_all.csv"):
    mapped_all = pd.read_csv(f"{output_folder}/dyad_mapped_all.csv")
else:
    file_paths = hf.get_files(output_folder, f"*_nucleoatac_dyad/*_mapped.csv")
    file_paths = [x for x in file_paths if ("nak" not in x and "noUV" not in x)]

    # Read all CSV files into DataFrames
    dfs = [pd.read_csv(file, index_col=0) for file in file_paths]

    # Merge them sequentially on the specified columns
    mapped_all = dfs[0]  # Start with the first DataFrame
    for df in dfs[1:]:    # Iterate through the remaining DataFrames
        mapped_all = pd.merge(mapped_all, df, on=["quantile_label", "window"], how='outer')

    file_paths = hf.get_files(output_folder, f"*_nucleoatac_dyad_sim_norm/*_mapped.csv")
    file_paths = [x for x in file_paths if ("nak" not in x and "noUV" not in x)]

    # Read all CSV files into DataFrames
    dfs = [pd.read_csv(file, index_col=0) for file in file_paths]

    # Merge them sequentially on the specified columns
    mapped_all_sim_norm = dfs[0]  # Start with the first DataFrame
    for df in dfs[1:]:    # Iterate through the remaining DataFrames
        mapped_all_sim_norm = pd.merge(mapped_all_sim_norm, df, on=["quantile_label", "quantile", "window"], how='outer')

    file_paths = hf.get_files(output_folder, f"*_nucleoatac_dyad_double_norm/*_mapped.csv")
    # file_paths = [x for x in file_paths if ("nak" not in x and "noUV" not in x)]

    # Read all CSV files into DataFrames
    dfs = [pd.read_csv(file, index_col=0) for file in file_paths]

    # Merge them sequentially on the specified columns
    mapped_all_double_norm = dfs[0]  # Start with the first DataFrame
    for df in dfs[1:]:    # Iterate through the remaining DataFrames
        mapped_all_double_norm = pd.merge(mapped_all_double_norm, df, on=["quantile_label", "quantile", "window"], how='outer')

    mapped_all = pd.merge(mapped_all, mapped_all_sim_norm, on=["quantile_label", "window"], how='outer')
    mapped_all = pd.merge(mapped_all, mapped_all_double_norm, on=["quantile_label", "quantile", "window"], how='outer')

    tt_results_agg = create_tt_freq(region_size=80)
    tt_results_agg["window"] = tt_results_agg["window"] * -1

    mapped_all = pd.merge(mapped_all, tt_results_agg, on=["quantile", "window"], how='outer')

    mapped_all.to_csv(f"{output_folder}/dyad_mapped_all.csv", index=False)

mapped_high = mapped_all[mapped_all["quantile_label"] == "High"]

color_map_repair = {
    "5 min": "#A3B57A",
    "1 hour":   "#D4B46D",
    "8 hours":   "#BF6B57",
    "24 hours":  "#8A5D71",
    "48 hours":  "#4F4D63"
}

color_map_damage = {
    "0 hour": "#85BFC4",
    "8 hours": "#C2665A",
    "24 hours": "#94617D",
    "48 hours": "#4B4B66"
}

os.makedirs(f"{output_folder}/c_elegans_dyad", exist_ok=True)

fig = plt.figure(figsize=(8, 4))
ax1 = fig.add_subplot(111)

sns.lineplot(
    data=smooth_rolling(mapped_high, "wt_0h_CPD_ds_merged_sim"),
    x='window',
    y='wt_0h_CPD_ds_merged_sim',
    color="black",
    linestyle="--",
    label='Simulated CPD Damages (0 hour)',
    ax=ax1
)
sns.lineplot(
    data=smooth_rolling(mapped_high, "wt_0h_CPD_ds_merged"),
    x='window',
    y='wt_0h_CPD_ds_merged',
    color=color_map_damage["0 hour"],
    label='CPD Damages (0 hour)',
    ax=ax1
)

ax1.margins(x=0)
ax2 = ax1.twinx()
xticks = [min(mapped_high["window"]), 0, max(mapped_high["window"])]
xticklabels = [f'-80', 'Dyad Center', f'80']
snr.arrange_bar_dyad(ax2, x_min=-80, x_max=80, bar_width=5)
hf.format_plot(df_plot=mapped_high, ax=ax1, xticks=xticks, xticklabels=xticklabels)

lines_1, labels_1 = ax1.get_legend_handles_labels()
labels_1 = labels_1[::-1]
lines_1 = lines_1[::-1]
ax1.legend(lines_1, labels_1, title='Time After UV', loc='upper left').get_title().set_fontweight('bold')

ax1.set_ylabel("RPKM")
ax1.set_ylim(5.5,10)
ax1.set_xlabel('')
ax1.spines['bottom'].set_zorder(0)
ax2.spines['bottom'].set_zorder(0)

plt.savefig(f'{output_folder}/c_elegans_dyad/cpd_sim_vs_real.svg'.lower().replace(" ", "_"), dpi=300)
plt.close()


fig = plt.figure(figsize=(8, 4))
ax1 = fig.add_subplot(111)

sns.lineplot(
    data=smooth_rolling(mapped_high, "wt_0h_64_ds_merged_sim"),
    x='window',
    y='wt_0h_64_ds_merged_sim',
    color="black",
    linestyle="--",
    label='Simulated 64 Damages (0 hour)',
    ax=ax1
)
sns.lineplot(
    data=smooth_rolling(mapped_high, "wt_0h_64_ds_merged"),
    x='window',
    y='wt_0h_64_ds_merged',
    color=color_map_damage["0 hour"],
    label='64 Damages (0 hour)',
    ax=ax1
)

ax1.margins(x=0)
ax2 = ax1.twinx()
xticks = [min(mapped_high["window"]), 0, max(mapped_high["window"])]
xticklabels = [f'-80', 'Dyad Center', f'80']
snr.arrange_bar_dyad(ax2, x_min=-80, x_max=80, bar_width=5)
hf.format_plot(df_plot=mapped_high, ax=ax1, xticks=xticks, xticklabels=xticklabels)

lines_1, labels_1 = ax1.get_legend_handles_labels()
labels_1 = labels_1[::-1]
lines_1 = lines_1[::-1]
ax1.legend(lines_1, labels_1, title='Time After UV', loc='upper left').get_title().set_fontweight('bold')

ax1.set_ylabel("RPKM")
ax1.set_xlabel('')
ax1.set_ylim(5.5,10)
ax1.spines['bottom'].set_zorder(0)
ax2.spines['bottom'].set_zorder(0)

plt.savefig(f'{output_folder}/c_elegans_dyad/64_sim_vs_real.svg'.lower().replace(" ", "_"), dpi=300)
plt.close()


fig = plt.figure(figsize=(8, 4))
ax1 = fig.add_subplot(111)

sns.lineplot(
    data=smooth_rolling(mapped_high, "wt_0h_CPD_ds_merged_sim_norm"),
    x='window',
    y='wt_0h_CPD_ds_merged_sim_norm',
    color=color_map_damage["0 hour"],
    label='CPD Damages (0 hour)',
    ax=ax1
)

ax2 = ax1.twinx()
sns.lineplot(
    data=smooth_rolling(mapped_high, "total_count_TT"),
    x='window',
    y='total_count_TT',
    color="black",
    label='TT Count',
    linestyle="--",
    ax=ax2
)
ax2.set_ylabel("TT Count", color='black', fontsize=12)
ax2.tick_params(axis='y', labelcolor='black')

ax1.margins(x=0)
ax2.margins(x=0)
xticks = [min(mapped_high["window"]), 0, max(mapped_high["window"])]
xticklabels = [f'-80', 'Dyad Center', f'80']
snr.arrange_bar_dyad(ax2, x_min=-80, x_max=80, bar_width=5)
hf.format_plot(df_plot=mapped_high, ax=ax1, xticks=xticks, xticklabels=xticklabels)
hf.format_plot(df_plot=mapped_high, ax=ax2, xticks=xticks, xticklabels=xticklabels)

lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
ax2.legend(lines_1 + lines_2, labels_1 + labels_2, title='Time After UV', loc='upper left').get_title().set_fontweight('bold')
ax1.get_legend().remove()

ax1.set_ylabel("Simultion Normalized\nCounts")
ax2.set_ylabel("TT Count\n", rotation=270, va="center")
ax1.set_ylim(0.6,1.2)
ax1.set_xlabel('')
ax2.set_xlabel('')
ax1.spines['bottom'].set_zorder(0)
ax2.spines['bottom'].set_zorder(0)
plt.savefig(f'{output_folder}/c_elegans_dyad/cpd_vs_tt.svg'.lower().replace(" ", "_"), dpi=300)
plt.close()


fig = plt.figure(figsize=(8, 4))
ax1 = fig.add_subplot(111)

sns.lineplot(
    data=smooth_rolling(mapped_high, "wt_0h_64_ds_merged_sim_norm"),
    x='window',
    y='wt_0h_64_ds_merged_sim_norm',
    color=color_map_damage["0 hour"],
    label='(6-4)PP Damages (0 hour)',
    ax=ax1
)

ax2 = ax1.twinx()
sns.lineplot(
    data=smooth_rolling(mapped_high, "total_count_TT"),
    x='window',
    y='total_count_TT',
    color="black",
    label='TT Count',
    linestyle="--",
    ax=ax2
)
ax2.set_ylabel("TT Count", color='black', fontsize=12)
ax2.tick_params(axis='y', labelcolor='black')

ax1.margins(x=0)
ax2.margins(x=0)
xticks = [min(mapped_high["window"]), 0, max(mapped_high["window"])]
xticklabels = [f'-80', 'Dyad Center', f'80']
snr.arrange_bar_dyad(ax2, x_min=-80, x_max=80, bar_width=5)
hf.format_plot(df_plot=mapped_high, ax=ax1, xticks=xticks, xticklabels=xticklabels)
hf.format_plot(df_plot=mapped_high, ax=ax2, xticks=xticks, xticklabels=xticklabels)

lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
ax2.legend(lines_1 + lines_2, labels_1 + labels_2, title='Time After UV', loc='upper left').get_title().set_fontweight('bold')
ax1.get_legend().remove()

ax1.set_ylabel("Simultion Normalized\nCounts")
ax2.set_ylabel("TT Count\n", rotation=270, va="center")
ax1.set_ylim(0.6,1.2)
ax1.set_xlabel('')
ax2.set_xlabel('')
ax1.spines['bottom'].set_zorder(0)
ax2.spines['bottom'].set_zorder(0)
plt.savefig(f'{output_folder}/c_elegans_dyad/64_vs_tt.svg'.lower().replace(" ", "_"), dpi=300)
plt.close()



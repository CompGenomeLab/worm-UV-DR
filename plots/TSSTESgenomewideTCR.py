import os
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend suitable for scripts
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

# Functions (calculate_rpkm, prepare_plot_data, create_graph) remain unchanged

def calculate_rpkm(read_counts, total_reads_million):
    gene_length_kb_per_bin = 0.015  # 15 nucleotides = 0.015 kb
    return (read_counts / gene_length_kb_per_bin) / total_reads_million

def prepare_plot_data(df, total_reads_million):
    df['RPKM'] = calculate_rpkm(df[7], total_reads_million)  # Pass the total reads in millions
    grouped_data = df.groupby(6)['RPKM'].agg(['mean', 'std', 'count'])
    grouped_data['sem'] = grouped_data['std'] / np.sqrt(grouped_data['count'])
    grouped_data['ci95'] = 1.96 * grouped_data['sem']
    return grouped_data.reset_index()

def create_graph(sample_name, folder_path):
    files = {
        "tss_TS": os.path.join(folder_path, f"{sample_name}_div_tssTS.bed"),
        "tes_TS": os.path.join(folder_path, f"{sample_name}_div_tesTS.bed"),
        "tss_NTS": os.path.join(folder_path, f"{sample_name}_div_tssNTS.bed"),
        "tes_NTS": os.path.join(folder_path, f"{sample_name}_div_tesNTS.bed")
    }
    print(files)
    plot_data = {}
    total_reads_file = os.path.join(folder_path, f"{sample_name}_readCount.txt")
    with open(total_reads_file, 'r') as file:
        total_reads_million = float(file.readline().strip()) / 1e6

    for key, file_name in files.items():
        df = pd.read_csv(file_name, sep="\t", header=None)
        plot_data[key] = prepare_plot_data(df, total_reads_million)

    fig, axes = plt.subplots(1, 2, figsize=(8, 4), facecolor='white')
    colors = {'TS': 'darkgreen', 'NTS': 'limegreen'}

    for ax, key in zip(axes, ['tss', 'tes']):
        for strand in ['TS', 'NTS']:
            data = plot_data[f'{key}_{strand}']
            ax.plot(data[6], data['mean'], color=colors[strand])
            ax.fill_between(data[6], data['mean'] - data['ci95'], data['mean'] + data['ci95'], alpha=0.2, color=colors[strand])

        ax.set_ylim(4, 20)
        ax.set_facecolor('white')
        ax.grid(False)
        if key == 'tss':
            ax.set_xticks([1, 33, 100])
            ax.set_xticklabels(['-500bp', 'TSS', ''])
            ax.axvline(x=33, color='black', linestyle='dotted')  # Dotted line at TSS
            ax.text(99, -0.085, '+1kb', transform=ax.get_xaxis_transform(), ha='center', fontsize=12)    # Adding '+1kb' as a label at bin 28
            ax.set_ylabel('RPKM', fontsize = 16)  # Set Y-axis label for TSS
            ax.spines['right'].set_visible(False)
        else:  # 'tes'
            ax.set_xticks([1, 67, 100])
            ax.set_xticklabels(['', 'TES', '+500bp'])
            ax.axvline(x=67, color='black', linestyle='dotted')  # Dotted line at TES
            ax.text(2, -0.085, '-1kb', transform=ax.get_xaxis_transform(), ha='center', fontsize=12)    # Adding '-1kb' as a label at bin 2
            ax.yaxis.tick_right()
            ax.spines['left'].set_visible(False)  # Hide the left spine
        ax.tick_params(axis='y', labelsize=11)
        ax.tick_params(axis='x', labelsize=12)
    plt.suptitle(f"{sample_name}", fontsize=14)
    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    plt.subplots_adjust(wspace=0.016)

    ts_line = Line2D([0], [0], color=colors['TS'], lw=4, label='TS')  # Adding label for TS
    nts_line = Line2D([0], [0], color=colors['NTS'], lw=4, label='NTS')  # Adding label for NTS
    ci_patch = Patch(color='green', alpha=0.2, label='95% CI')  # CI label already present

    legend = plt.legend(handles=[ts_line, nts_line, ci_patch], ncol=3, fontsize=10, frameon=False,
                        bbox_to_anchor=(0.5, -0.15), loc='upper center')
    plt.setp(legend.get_texts(), color='black')  # Set legend text color to black
    pdf_filename = os.path.join(folder_path, f"{sample_name}_tcr.pdf")
    plt.savefig(pdf_filename)
    plt.savefig(png_filename)


if __name__ == "__main__":
    folder_path = "../output/genomewidetcr"  # Specify the folder path
    read_count_files = [f for f in os.listdir(folder_path) if f.endswith('_readCount.txt')]
    sample_names = [os.path.splitext(f)[0].replace('_readCount', '') for f in read_count_files]

    for sample_name in sample_names:
        print(f"Processing sample: {sample_name}")
        create_graph(sample_name, folder_path)

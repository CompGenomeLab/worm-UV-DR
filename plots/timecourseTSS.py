import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend suitable for scripts
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#############################################################################
# Directory where input BED/readCount files live, and where PDFs are saved
WORKING_DIR = "../output/TSSquartile/xr"
#for damage se	q, change the working directory: ../output/TSSquartile/damage 
os.chdir(WORKING_DIR)
#############################################################################


def calculate_rpkm(read_counts, total_reads_million):
    """
    Calculate RPKM for a given column of read counts.
    Adjust gene_length_kb_per_bin if each bin is 50 nt (0.05 kb), etc.
    """
    gene_length_kb_per_bin = 0.01
    return (read_counts / gene_length_kb_per_bin) / total_reads_million


def prepare_plot_data(df, total_reads_million):
    """
    Given a DataFrame with columns (6 -> bin, 7 -> read counts),
    calculate RPKM, group by bin, return mean, std, sem, ci95.
    """
    df['RPKM'] = calculate_rpkm(df[7], total_reads_million)
    grouped_data = df.groupby(6)['RPKM'].agg(['mean', 'std', 'count'])
    grouped_data['sem'] = grouped_data['std'] / np.sqrt(grouped_data['count'])
    grouped_data['ci95'] = 1.96 * grouped_data['sem']
    return grouped_data.reset_index()


def load_and_merge_quartiles(sample_name, strand, quartiles):
    """
    Load and merge the given quartiles for a sample on a particular strand.
    quartiles is a list like [(0,25), (25,50)] or [(50,75), (75,100)].
    Filenames follow:
      [strain]_[timepoint]_[damage]_[replicate]_TSS_[strain]_[start]_[end]_quartile_final_150bins_[strand].bed
    """
    parts = sample_name.split('_')
    if len(parts) < 4:
        raise ValueError(f"sample_name '{sample_name}' does not have at least 4 parts.")

    strain     = parts[0]  # e.g. "xpc"
    timepoint  = parts[1]  # e.g. "8h"
    damage     = parts[2]  # e.g. "CPD"
    replicate  = parts[3]  # e.g. "rep2"

    dataframes = []
    for (start, end) in quartiles:
        file_name = f"{strain}_{timepoint}_{damage}_{replicate}_TSS_{strain}_{start}_{end}_quartile_final_150bins_{strand}.bed"
        df_tmp = pd.read_csv(file_name, sep="\t", header=None)
        dataframes.append(df_tmp)

    # Concatenate row-wise
    df_merged = pd.concat(dataframes, ignore_index=True)
    return df_merged


def create_graph(samples_with_labels, strand, quartiles, quartile_label):
    """
    Generate & save a PDF plot for the given quartiles (e.g. 0–50 or 50–100).
    
    Args:
      samples_with_labels: ["xpc_8h_CPD_rep2:8h", ...]
      strand: "TS" or "NTS"
      quartiles: [(start, end), (start2, end2)] — e.g. [(0,25),(25,50)]
      quartile_label: "0to50_quartile" or "50to100_quartile"
    """

    # Choose colormap based on which quartiles we are plotting
    if "0to50" in quartile_label:
        cmap = plt.cm.Blues   # more bluish
    else:
        cmap = plt.cm.Reds    # more reddish
    
    plt.figure(figsize=(6, 5))
    ax = plt.gca()
    n_samples = len(samples_with_labels)

    # Use a range in the chosen colormap
    colors = cmap(np.linspace(0.3, 1, n_samples))

    for index, sample_with_label in enumerate(samples_with_labels):
        sample_name, label = sample_with_label.split(':')

        # The read-count file is assumed to be [sample_name]_readCount.txt
        total_reads_file = f"{sample_name}_readCount.txt"
        with open(total_reads_file, 'r') as f:
            total_reads_million = float(f.readline().strip()) / 1e6

        # Load & merge chosen quartiles
        df_merged = load_and_merge_quartiles(sample_name, strand, quartiles)
        plot_data = prepare_plot_data(df_merged, total_reads_million)

        ax.plot(plot_data[6], plot_data['mean'], 
                label=label, 
                color=colors[index], 
                linewidth=2)

        ax.fill_between(plot_data[6],
                        plot_data['mean'] - plot_data['ci95'],
                        plot_data['mean'] + plot_data['ci95'],
                        color=colors[index],
                        alpha=0.2)

    # Mark TSS
    ax.axvline(x=51, color='gray', linestyle='dotted')

    ax.set_ylim(1, 60)
    ax.set_xlim(1, 150)
    ax.set_xticks([1, 51, 150])
    ax.set_xticklabels(['-0.5kb', 'TSS', '+1kb'])
    ax.set_ylabel('RPKM', fontsize=18)

    # Legend
    ax.legend(loc='upper center',
              bbox_to_anchor=(0.5, -0.15),
              ncol=n_samples,
              fontsize=12,
              frameon=False)

    # Title & PDF naming
    last_arg = sys.argv[-1]
    strain   = sys.argv[-2]
    plt.title(f"{strain} Time Course {last_arg} - TSS ({strand}, {quartile_label})", 
              fontsize=16, pad=20)

    plt.tight_layout()
    plt.subplots_adjust(top=0.85, bottom=0.25)

    pdf_filename = f"{strain}_{last_arg}_{quartile_label}_TSS_{strand}.pdf"
    plt.savefig(pdf_filename, bbox_inches='tight')
    plt.close()

    print(f"Graph saved as {pdf_filename}")


if __name__ == "__main__":
   

    # The sample name:label arguments
    samples_with_labels = sys.argv[1:-2]

    # 1) Plot upper half (50–100) in Reds
    create_graph(
        samples_with_labels, 
        strand='TS',
        quartiles=[(50,75), (75,100)],
        quartile_label="50to100_quartile"
    )
    create_graph(
        samples_with_labels, 
        strand='NTS',
        quartiles=[(50,75), (75,100)],
        quartile_label="50to100_quartile"
    )

    # 2) Plot lower half (0–50) in Blues
    create_graph(
        samples_with_labels, 
        strand='TS',
        quartiles=[(0,25), (25,50)],
        quartile_label="0to50_quartile"
    )
    create_graph(
        samples_with_labels, 
        strand='NTS',
        quartiles=[(0,25), (25,50)],
        quartile_label="0to50_quartile"
    )

import os
import itertools
import pandas as pd
import numpy as np
import bioframe as bf
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import find_peaks


import scipy.signal as signal
from statsmodels.stats.multitest import multipletests

try:
    import workflow.scripts.periodicity_helper_func as hf
except:
    import periodicity_helper_func as hf

def calculate_power_spectrum(S, limit=26, normalize=True):
    N = len(S)
    S_star = []
    
    # Define the range of periods P
    # For simplicity, we use a range from 2 to N
    P_values = np.arange(2, limit)
    
    for P in P_values:
        # Calculate the sum for each period P
        sum_term = np.sum(S * np.exp(-2j * np.pi * np.arange(N) / P))
        
        # Calculate the magnitude squared and normalize by N
        S_star_P = (1 / N) * np.abs(sum_term) ** 2
        S_star.append(S_star_P)

    if normalize:
        S_star = S_star / np.mean(S_star)

    return P_values, S_star

def compute_snr(S_star, target_period=16):
    # Calculate the SNR
    max_power = S_star[target_period - 2]
    S_star_without_max = np.delete(S_star, target_period - 2)
    snr = max_power / np.mean(S_star_without_max)
    return snr

def calculate_p_value(observed_snr, permuted_snrs):
    # Convert to NumPy arrays for vectorized operations
    permuted_snrs = np.array(permuted_snrs)

    # Calculate p-value only for the observed SNR
    p_value = (np.sum(permuted_snrs > observed_snr) + 1) / (len(permuted_snrs) + 1)

    return p_value

def permutation_snr(S, n_permutations=1000, limit=26, target_period=16, xpc=False):

    # Calculate the power spectrum
    if xpc:
        normalize = True
    else:
        normalize = True

    P_values, S_star = calculate_power_spectrum(S, limit=limit, normalize=normalize)

    if target_period == "max":
        target_period = P_values[np.argmax(S_star)]

    # Calculate the SNR
    observed_snr = compute_snr(S_star, target_period)

    # Perform permutations
    permuted_snrs = []
    for _ in range(n_permutations):
        permuted_S = np.random.permutation(S)
        _, permuted_S_star = calculate_power_spectrum(permuted_S, limit)
        permuted_snr = compute_snr(permuted_S_star, target_period)
        permuted_snrs.append(permuted_snr)

    # Calculate p-values for all SNRs
    all_snr = [observed_snr] + permuted_snrs

    p_value = calculate_p_value(observed_snr, permuted_snrs)

    return observed_snr, p_value, permuted_snrs, P_values, S_star


def plot_periodogram(signal_rate, P_values, S_star, example_shuffled_data, example_P_values, example_S_star, permuted_snrs, observed_snr, observed_p_value, output, name):
    """
    Enhanced visualization that includes the power vs. period plot.
    """
    window = np.linspace(-100, 100, 201)

    plt.figure(figsize=(5,4))

    # Original Signal data
    plt.subplot(3, 2, 1)
    plt.plot(window, signal_rate, label=f"Original Signal\nObserved SNR: {observed_snr:.2f}")
    plt.title("Original Signal Data")
    plt.xlabel("Window")
    plt.ylabel("Signal")
    plt.legend()

    # Original periodogram: frequency vs. power
    plt.subplot(3, 2, 2)
    plt.plot(P_values, S_star, label="Original Power Spectrum S*(P)")
    plt.axvline(P_values[np.argmax(S_star)], color='red', linestyle='--', label="Max Power Period")
    plt.title("Original Periodogram")
    plt.xlabel('Period P')
    plt.ylabel("Power")
    plt.legend()


    plt.subplot(3, 2, 3)
    plt.plot(window, example_shuffled_data, label="Permuted Signal (Example)")
    plt.title("Permuted Signal Data Example")
    plt.xlabel("Window")
    plt.ylabel("Signal")
    plt.legend()

    # Example permuted periodogram
    plt.subplot(3, 2, 4)
    plt.plot(example_P_values, example_S_star, label="Permuted Power Spectrum S*(P)")
    plt.axvline(example_P_values[np.argmax(example_S_star)], color='red', linestyle='--', label="Max Power Period")
    plt.title("Permuted Periodogram")
    plt.xlabel('Period P')
    plt.ylabel("Power")
    plt.legend()

    # Permuted SNR distribution
    plt.subplot(3, 1, 3)
    plt.hist(permuted_snrs, bins=30, alpha=0.7, label="Permuted SNR Distribution")
    plt.axvline(observed_snr, color='red', linestyle='--', label=f"Observed SNR: {observed_snr:.2f}")
    plt.title("Permutation Test: SNR Distribution")
    plt.xlabel("SNR")
    plt.ylabel("Frequency")
    plt.annotate(f"p-value: {observed_p_value:.3f}",
                 xy=(observed_snr, 6), xytext=(observed_snr + 2, 8),
                 arrowprops=dict(facecolor='black', arrowstyle="->"),
                 fontsize=10)
    plt.legend()

    plt.tight_layout()
    plt.savefig(f"{output}/{name}_snr.png", dpi=300)

def scale_peaks(signal):

    # Find peaks
    peaks, _ = find_peaks(signal)

    # Scale each peak to a fixed height (e.g., 1)
    signal_normalized = signal.copy()
    for peak in peaks:
        if signal[peak] != 0:
            signal_normalized[peak] = max(signal)

    return signal_normalized

import numpy as np
from scipy.signal import find_peaks, convolve
from scipy.signal.windows import gaussian 

def scale_peaks_v2(signal, width=5):
    # Find peaks
    peaks, _ = find_peaks(signal, width=1)
    
    # Create a copy of the signal
    signal_widened = np.zeros_like(signal)

    # Define a Gaussian kernel for smoothing
    kernel = gaussian(width * 2, std=width / 3)  # Adjust std for sharpness
    kernel /= kernel.max()  # Normalize

    # Expand each peak using the kernel
    for peak in peaks:
        start = max(0, peak - width)
        end = min(len(signal), peak + width)
        signal_widened[start:end] += kernel[:end-start]

    # Normalize to keep the peaks at the same scale
    signal_widened *= max(signal) / max(signal_widened)

    return signal_widened


def possible_combinations(norm=False):

    cell_line = ["csb", "xpc", "wt", "NHF1"]
    experiments = ["xr", "ds"]
    conditions = ["6-4PP", "64", "CPD"]
    #replicates = ["_Rep1_", "_Rep2_", "_A_", "_B_"]
    types = ["sim", "obs"]

    #target_list_of_list = list(itertools.product(cell_line, experiments, conditions, replicates, types))
    if norm:
        target_list_of_list = list(itertools.product(cell_line, conditions, experiments))
    else:
        target_list_of_list = list(itertools.product(cell_line, conditions, experiments, types))
    return [list(combination) for combination in target_list_of_list]

#### MAIN ####
def process_snr(file, expanded_regions, output):

    print(file, "started", flush=True)

    df = bf.read_table(file, schema="bed")

    overlapped = bf.count_overlaps(expanded_regions[["chrom", "start", "end", "quantile", "window"]], df)
    name = os.path.basename(file).replace(".bed", "")
    overlapped.rename(columns={"count": name}, inplace=True)

    with open(file, 'r') as f:
        total_reads = sum(1 for _ in f)

    if "_plus" in name:
        opposite_strand = file.replace("_plus.bed", "_minus.bed")
    else:
        opposite_strand = file.replace("_minus.bed", "_plus.bed")

    with open(opposite_strand, 'r') as f_opp:
        total_reads += sum(1 for _ in f_opp)

    overlapped[f"{name}"] = hf.calculate_rpkm(
        region_df=overlapped,
        read_counts_column=f"{name}",
        total_reads=total_reads
    )


    signal_rate = overlapped[overlapped["quantile"]==max(overlapped["quantile"])][["window",name]].groupby("window").mean().reset_index()[name].values
    signal_rate_0 = signal_rate - np.mean(signal_rate)

    observed_snr, p_value, permuted_snrs, P_values, S_star = permutation_snr(signal_rate_0, n_permutations=1000)

    # Example permuted Signal
    example_shuffled_data = np.random.permutation(signal_rate_0)
    _, _, _, example_P_values, example_S_star = permutation_snr(example_shuffled_data)

    plot_periodogram(signal_rate_0, P_values, S_star, example_shuffled_data, example_P_values, example_S_star, permuted_snrs, observed_snr, p_value, output, name)

    print(file, "finished", flush=True)

def get_time_label(name):
    if ("_0_h" in name) or ("_0h_" in name):
        return '0 hour'
    elif ("_5_min" in name) or ("_5min_" in name):
        return '5 min'
    elif ("_20_min" in name) or ("_20min_" in name):
        return '20 min'
    elif ("_1_hr" in name) or ("_1h_" in name):
        return '1 hour'
    elif ("_2_hr" in name) or ("_2h_" in name):
        return '2 hours'
    elif ("_4_hr" in name) or ("_4h_" in name):
        return '4 hours'
    elif ("_8_hr" in name) or ("_8h_" in name):
        return '8 hours'
    elif ("_16_hr" in name) or ("_16h_" in name):
        return '16 hours'
    elif ("_24_hrs" in name) or ("_24h_" in name):
        return '24 hours'
    elif ("_36_hrs" in name) or ("_36h_" in name):
        return '36 hours'
    elif ("_48_hrs" in name) or ("_48h_" in name):
        return '48 hours'
    elif ("nak" in name):
        return 'Naked'
    elif ("noUV_dip" in name):
        return 'No UV dip'
    elif ("noUV" in name):
        return 'No UV'
    else:
        return 'Unknown'

def plot_periodicity(data, name, output_folder, ymax=18, xticks=[40, 80, 120, 160, 200, 240], vline=160):

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

    if "_ds_" in name:
        timepoint_color_palette = color_map_damage
    if "_xr_" in name:
        timepoint_color_palette = color_map_repair

    name_order = [key for key in timepoint_color_palette.keys() if key in data['name_label'].unique()]
    fontsize = 18

    def plot(data=data, name=name):
        fig, ax = plt.subplots(figsize=(5,4))
        sns.lineplot(x='period', y='signal', data=data, hue='name_label', hue_order= name_order, palette=timepoint_color_palette, ax=ax, linewidth=3)

        # draw a vertical line at the period of 16
        plt.axvline(x=vline, color='gray', linestyle='--')

        if not ymax=="Free":
            plt.ylim(0, ymax)
        plt.xticks(xticks)
        #ax.get_legend().remove()
        ax.tick_params(axis='both')
        plt.tight_layout()
        hf.format_plot(df_plot=data, ax=ax, xticks=xticks)


        # if "64" in name and "obs" in name:
        plt.xlabel('Period (bp)')
        plt.ylabel('Power Spectrum')
        legend = plt.legend(title='Time After UV', loc='upper left', frameon=True, framealpha=0.95, fontsize=12)
        plt.setp(legend.get_title(), fontsize=12, weight="bold")

        # else:
        #     plt.ylabel('')
        #     ax.set_yticklabels([])
        #     ax.get_legend().remove()
        #     ax.set_xticklabels([])
        #     plt.xlabel('')

        plt.savefig(f"{output_folder}/{name}.svg")

    if "xpc" in name:
        data_ds = data[data["name"].str.contains("_ds")]
        name1 = name + "_downstream"
        if len(data_ds) > 0:
            plot(data=data_ds, name=name1)
        data_us = data[data["name"].str.contains("_us")]
        name2 = name + "_upstream"
        if len(data_us) > 0:
            plot(data=data_us, name=name2)
    else:
        plot()


def plot_snr(df, name, output_folder, ymax=50):

    name_order_all = ["0 hour", "5 min", "20 min", "1 hour", "2 hours", "4 hours", "8 hours", "16 hours", "24 hours", "36 hours", "48 hours"]
    name_order_all2 = ["0 h", "5 min", "20 min", "1 h", "2 h", "4 h", "8 h", "16 h", "24 h", "36 h", "48 h"]
    name_order = [key for key in name_order_all if key in df['name_label'].unique()]
    df["name_label"] = pd.Categorical(df["name_label"], categories=name_order, ordered=True)
    df["name_label_short"] = df["name_label"].str.replace(" hour", " h").str.replace(" hs", " h")
    name_order = [key for key in name_order_all2 if key in df['name_label_short'].unique()]
    df["name_label_short"] = pd.Categorical(df["name_label_short"], categories=name_order, ordered=True)
    fontsize=18
    quantile_labels = {0: 'Low', 1: 'Mid-Low', 2: 'Mid-High', 3: 'High'}
    df['quantile_label'] = df['quantile'].map(quantile_labels)
    quantile_order = ['High', 'Mid-High', 'Mid-Low', 'Low']
    custom_palette = {
    'High': '#2b2d42',
    'Mid-High': '#8d99ae',
    'Mid-Low': '#f49cbb',
    'Low': '#ef233c'
    }

    # Add replicate information
    # df['replicate'] = df['name'].apply(
    #     lambda x: (
    #         'Rep1' if 'Rep1' in x or '_A_' in x else
    #         'Rep2' if 'Rep2' in x or '_B_' in x else
    #         'Rep3' if 'Rep3' in x else
    #         'Rep4' if 'Rep4' in x else
    #         'unknown'
    #     )
    # )

    # Line styles
    # line_styles = {'Rep1': '', 'Rep2': (4, 2), 'Rep3': (1,1), 'Rep4': (4,2,1,2)}

    def plot(df=df, name=name):
        # Plotting
        fig, ax = plt.subplots(figsize=(5,4))
        sns.lineplot(
            x='name_label_short', y='snr',
            data=df,
            hue='quantile_label',
            hue_order=quantile_order,
            palette=custom_palette,
            #style='replicate',
            #dashes=line_styles,
            ax=ax,
            linewidth=3
        )

        # Add asterisks for significance
        name_to_xpos = {label.get_text(): pos for pos, label in zip(ax.get_xticks(), ax.get_xticklabels())}
        for i, row in df.iterrows():
            x_pos = name_to_xpos[row['name_label_short']]
            y_pos = row['snr']
            # if 0.01 < row['p_value'] <= 0.05:
            #     ax.text(x_pos, y_pos, '*', fontsize=fontsize, color='black', ha='center', va='bottom')
            if row['p_value'] <= 0.01:
                ax.text(x_pos, y_pos, '*', fontsize=fontsize, color='black', ha='center', va='center')
            # elif row['p_value'] <= 0.001:
            #     ax.text(x_pos, y_pos, '***', fontsize=fontsize, color='black', ha='center', va='bottom')

        # Adjustments
        if not ymax=="Free":
            plt.ylim(0, ymax)

        # Adjust legend for replicates and quantiles
        # handles, labels = ax.get_legend_handles_labels()
        #ax.legend(handles=handles, labels=labels, title='Quantile / Replicate', fontsize=fontsize - 2, title_fontsize=fontsize)
        # ax.get_legend().remove()
        plt.tight_layout()
        hf.format_plot(df_plot=df, ax=ax, xticks=name_order)

        # if "64" in name and "obs" in name:
        ax.tick_params(axis='both')
        ax.set_ylabel("Signal to Noise Ratio")
        ax.set_xlabel("Timepoint")
        legend = plt.legend(title='Quantiles', loc='upper left', frameon=True, framealpha=0.95, fontsize=12)
        plt.setp(legend.get_title(), fontsize=12, weight="bold")
        # else:
        #     plt.ylabel('')
        #     ax.set_yticklabels([])
        #     ax.get_legend().remove()
        #     ax.set_xticklabels([])
        #     plt.xlabel('')

        plt.savefig(f"{output_folder}/{name}_snr.svg")

    if "xpc" in name:
        df_ds = df[df["name"].str.contains("_ds")]
        name1 = name + "_downstream"
        if len(df_ds) > 0:
            plot(df=df_ds, name=name1)
        df_us = df[df["name"].str.contains("_us")]
        name2 = name + "_upstream"
        if len(df_us) > 0:
            plot(df=df_us, name=name2)
    else:
        plot()


def arrange_bar_dyad(ax, x_min, x_max, bar_width=5):
    # Find the range of centers so that 0 is always a center of a yellow bar
    # Start from 0, go outwards in both directions
    centers = np.arange(0, x_max + bar_width, bar_width)
    centers = np.concatenate([centers, np.arange(-bar_width, x_min - bar_width, -bar_width)])
    centers = np.sort(centers)

    bar_y = ax.get_ylim()[0] - 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    bar_height = 0.03 * (ax.get_ylim()[1] - ax.get_ylim()[0])

    for i, center in enumerate(centers):
        # Omit bars at -5, 0, +5
        if any(np.isclose(center, omit, atol=0.1) for omit in [-5, 0, 5]):
            continue

        # Determine color: yellow if (center/10) is integer, white otherwise
        #color = '#bdb76b' if (center % 10 == 0) else '#2d4739'
        color = 'black' if (center % 10 == 0) else 'gray'

        # Calculate left and right edges
        left = center - bar_width / 2
        right = center + bar_width / 2

        # Edge bars: half-width if they would extend beyond the axis
        if left < x_min:
            left = x_min -0.15
        if right > x_max:
            right = x_max +0.5

        width = right - left
        if width <= 0:
            continue

        ax.add_patch(
            plt.Rectangle(
                (left, bar_y),
                width,
                bar_height,
                color=color,
                transform=ax.transData,
                clip_on=False,
                linewidth=0
            )
        )

        # Dashed lines at yellow bars, except at 0
        if color == 'black' and not np.isclose(center, 0, atol=0.1):
            ax.axvline(center, color='gray', linestyle='--', linewidth=0.7, alpha=0.5, zorder=0)
        # if color == '#bdb76b' and not np.isclose(center, 0, atol=0.1):
        #     ax.axvline(center, color='gray', linestyle='--', linewidth=0.7, alpha=0.5, zorder=0)

    ax.set_ylim(ax.get_ylim()[0] - bar_height, ax.get_ylim()[1])


def arrange_bar_nucleosome(ax, x_min, x_max, width1, width2, bar_color1, bar_color2, central_width):
    """
    Alternate bars of width1 and width2, but for the central bar (center=0), use central_width and omit it.
    Dashed lines are drawn at the center of each bar_color1 (width1) bar, except at 0 and the omitted central bar.
    """
    bar_y = ax.get_ylim()[0] - 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    bar_height = 0.03 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    trans = ax.transData

    # Central bar (center=0, width=central_width) is omitted
    left = -central_width / 2
    right = central_width / 2

    # Omit the central bar (center=0)
    # Bars to the right
    x = right
    is_first = True  # Start with bar1 to the right of center
    while x < x_max:
        width = width1 if is_first else width2
        color = bar_color1 if is_first else bar_color2
        center_bar = x + width / 2
        if np.isclose(center_bar, 0, atol=0.1):
            width = central_width  # Use custom width for central bar
            x += width
            is_first = not is_first
            continue
        ax.add_patch(
            plt.Rectangle(
                (x, bar_y),
                min(width, x_max - x),
                bar_height,
                color=color,
                transform=trans,
                clip_on=False,
                linewidth=0
            )
        )
        # Dashed lines at center of bar1 (width1), except at omitted central bar
        if is_first and not np.isclose(center_bar, 0, atol=0.1):
            ax.axvline(center_bar, color='gray', linestyle='--', linewidth=1.5, alpha=0.5, zorder=0)
        x += width
        is_first = not is_first

    # Bars to the left
    x = left
    is_first = True
    while x > x_min:
        width = width1 if is_first else width2
        color = bar_color1 if is_first else bar_color2
        center_bar = x - width / 2
        if np.isclose(center_bar, 0, atol=0.1):
            width = central_width
            x -= width
            is_first = not is_first
            continue
        ax.add_patch(
            plt.Rectangle(
                (max(x - width, x_min), bar_y),
                min(width, x - x_min),
                bar_height,
                color=color,
                transform=trans,
                clip_on=False,
                linewidth=0
            )
        )
        # Dashed lines at center of bar1 (width1), except at omitted central bar
        if is_first and not np.isclose(center_bar, 0, atol=0.1):
            ax.axvline(center_bar, color='gray', linestyle='--', linewidth=1.5, alpha=0.5, zorder=0)
        x -= width
        is_first = not is_first

    ax.set_ylim(ax.get_ylim()[0] - bar_height, ax.get_ylim()[1])
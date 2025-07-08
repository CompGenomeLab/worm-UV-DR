import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import periodicity_helper_func as hf
import snr


output_folder="" # Directory of the output (should be same as the output directory in run_rep_dam.py)
region = ["atac", "nucleoatac_dyad"]
region_names = ["ATAC-seq", "Nucleosome Dyad"]

def organize_data(output_folder, regex):

    file_paths = hf.get_files(output_folder, regex)
    file_paths = [x for x in file_paths if ("nak" not in x and "noUV" not in x)]

    # Read all CSV files into DataFrames
    dfs = [pd.read_csv(file, index_col=0) for file in file_paths]

    # Merge them sequentially on the specified columns
    mapped_all = dfs[0]  # Start with the first DataFrame
    for df in dfs[1:]:    # Iterate through the remaining DataFrames
        mapped_all = pd.merge(mapped_all, df[["quantile_label", "window", f"{df.columns[-1]}"]], on=["quantile_label", "window"], how='outer')

    mapped_high = mapped_all[mapped_all["quantile_label"]=="High"].drop(columns=["quantile_label"])
    df_melted = mapped_high.melt(
        id_vars=['window'],
        var_name='variable',
        value_name='value'
    )

    df_melted["name_labels"] = df_melted["variable"].apply(snr.get_time_label)

    return df_melted


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

def plot_timeseries(df_melted, cell_line, method, product, region, name, norm=""):

    name_order = ['0 hour', '5 min', '1 hour', '8 hours', '24 hours', '48 hours']

    if norm == "sim_norm" or norm == "double_norm":
        df_plot = df_melted[~df_melted["variable"].str.contains("nak") &
                        ~df_melted["variable"].str.contains("noUV") &
                        df_melted["variable"].str.contains(method) &
                        df_melted["variable"].str.contains(cell_line) &
                        df_melted["variable"].str.contains(product)]
    else:
        df_plot = df_melted[~df_melted["variable"].str.contains("nak") &
                        ~df_melted["variable"].str.contains("noUV") &
                        df_melted["variable"].str.contains(method) &
                        df_melted["variable"].str.contains(cell_line) &
                        df_melted["variable"].str.contains(product) &
                        ~df_melted["variable"].str.contains("sim")]

    # Apply rolling average smoothing
    if region == "nucleoatac_dyad":
        window_size = 3  # Adjust this value to control smoothing amount
    else:
        window_size = 5  # Adjust this value to control smoothing amount
    df_plot = df_plot.sort_values(['name_labels', 'window'])
    df_plot['value'] = df_plot.groupby('name_labels')['value'].transform(lambda x: x.rolling(window=window_size, center=True, min_periods=1).mean())


    name_order_curr = [name for name in name_order if name in df_plot["name_labels"].unique()]

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

    if method == "ds":
        timepoint_color_palette = color_map_damage
    elif method == "xr":
        timepoint_color_palette = color_map_repair


    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)
    sns.lineplot(
        data=df_plot,
        x='window',
        y='value',
        hue='name_labels',
        palette=timepoint_color_palette,
        hue_order=name_order_curr,
        ax=ax
    )

    if method == "ds" and norm == "sim_norm" and region == "atac":
        plt.ylim(0.6,1.2)
    elif method == "xr" and norm == "double_norm" and region == "atac":
        plt.ylim(0.8,5)
    elif method == "xr" and norm == "":
        plt.ylim(3,45)
    elif method == "ds" and norm == "" and region == "atac":
        plt.ylim(3.5, 6.5)
    elif method == "ds" and norm == "sim_norm" and region == "nucleoatac_dyad":
        plt.ylim(0.6,1.25)

    if region == "atac":
        ax.axvspan(380, 545, color="gray", alpha=0.3)

    ax.margins(x=0)
    if region=="nucleoatac_dyad":
        xticks = [min(df_plot["window"]), 0, max(df_plot["window"])]
        xticklabels = [f'-80', '0', f'80']
        plt.xlabel(f'Distance from\n{name.capitalize()} Center (bp)')
        arrange_bar_dyad(ax, x_min=-80, x_max=80, bar_width=5)
    else:
        xticklabels = ['-1', '0', '1']
        xticks = [min(df_plot["window"]), 0, max(df_plot["window"])]
        plt.xlabel(f'Distance from\n{name.capitalize()} Center (kb)')
        if region == "nucleoatac":
            central_width = 150
        if region == "atac":
            central_width = 120
        arrange_bar_nucleosome(
            ax, x_min=-1000, x_max=1000,
            width1=10, width2=150,
            bar_color1='black', bar_color2='gray',
            central_width=central_width
        )

    hf.format_plot(df_plot=df_plot, ax=ax, xticks=xticks, xticklabels=xticklabels)

    # legend = plt.legend(title='Time After UV', loc='upper left', frameon=True, framealpha=0.95, fontsize=12)
    # plt.setp(legend.get_title(), fontsize=12, weight="bold")
    plt.legend().remove()


    ext = "_" + norm
    if norm == "sim_norm":
        plt.ylabel('Simulation Normalized\nRPKM')
    elif norm == "double_norm":
        plt.ylabel('Simulation Normalized\nRelative Repair')
    else:
        plt.ylabel('RPKM')
        ext = ""
    ax.spines['bottom'].set_zorder(0)
    plt.savefig(f'results/c_elegans_ce11_final/c_elegans_{cell_line}_{region}{ext}/{cell_line}_{method}_{product}_{region}.svg'.lower().replace(" ", "_"), dpi=300)

    plt.close()


for i in range(len(region)):
    df_melted = organize_data(output_folder, f"*_{region[i]}/*_mapped.csv")
    for cell_line in ["csb", "xpc", "wt"]:
            for product in ["64", "CPD"]:
                method = "xr"
                plot_timeseries(df_melted, cell_line, method, product, region[i], region_names[i])
                if cell_line == "wt":
                    method = "ds"
                    plot_timeseries(df_melted, cell_line, method, product, region[i], region_names[i])

for i in range(len(region)):
    df_melted = organize_data(output_folder, f"*_{region[i]}_sim_norm/*_mapped.csv")
    for cell_line in ["csb", "xpc", "wt"]:
            for product in ["64", "CPD"]:
                method = "xr"
                plot_timeseries(df_melted, cell_line, method, product, region[i], region_names[i], norm="sim_norm")
                if cell_line == "wt":
                    method = "ds"
                    plot_timeseries(df_melted, cell_line, method, product, region[i], region_names[i], norm="sim_norm")

for i in range(len(region)):
    df_melted = organize_data(output_folder, f"*_{region[i]}_double_norm/*_mapped.csv")
    for cell_line in ["csb", "xpc", "wt"]:
            for product in ["64", "CPD"]:
                method = "xr"
                plot_timeseries(df_melted, cell_line, method, product, region[i], region_names[i], norm="double_norm")
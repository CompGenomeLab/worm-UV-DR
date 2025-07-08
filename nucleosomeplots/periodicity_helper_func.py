
import os
import glob
import pandas as pd
import numpy as np
import bioframe as bf
import matplotlib.pyplot as plt
import seaborn as sns
import re

from Bio import SeqIO

import snr

from datetime import datetime
from zoneinfo import ZoneInfo
from functools import wraps
import inspect

def timing_decorator(timezone='Europe/Istanbul'):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Get the start time
            start_time = datetime.now(ZoneInfo(timezone))

            # Inspect the function's signature to find the 'file' argument
            sig = inspect.signature(func)
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()

            # Extract the value of the 'file' argument if it exists
            file_value = bound_args.arguments.get('file', None)

            # Log the start message with the function name and file value
            if file_value:
                print(start_time.strftime('%Y-%m-%d %H:%M:%S'), f"{func.__name__} ({file_value}) started.", flush=True)
            else:
                print(start_time.strftime('%Y-%m-%d %H:%M:%S'), f"{func.__name__} started.", flush=True)

            # Execute the function
            result = func(*args, **kwargs)

            # Get the end time
            end_time = datetime.now(ZoneInfo(timezone))
            duration = end_time - start_time

            # Log the end message with the function name and file value
            if file_value:
                print(end_time.strftime('%Y-%m-%d %H:%M:%S'), f"{func.__name__} ({file_value}) finished. Duration:", duration, flush=True)
            else:
                print(end_time.strftime('%Y-%m-%d %H:%M:%S'), f"{func.__name__} finished. Duration:", duration, flush=True)

            return result
        return wrapper
    return decorator

def format_plot(df_plot, ax=None, fontsize=18, xticklabels=None, xticks="default", title=None):
    """Apply consistent formatting to matplotlib plots.
    
    Args:
        ax: matplotlib axes object (defaults to current axes)
        fontsize: base font size (default 12)
        title: plot title (optional)
    """
    if ax is None:
        ax = plt.gca()
    
    # Apply formatting
    ax.grid(False)
    ax.tick_params(axis='both', labelsize=fontsize)

    if xticks == "default":
        xticks = [min(df_plot["window"]), 0, max(df_plot["window"])]

    ax.set_xticks(xticks)
    if xticklabels:
        ax.set_xticklabels(xticklabels)
    else:
        ax.set_xticklabels(xticks)

    # Make all spines visible and solid
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1)
    
    if title:
        ax.set_title(title, fontsize=fontsize+2, pad=10)
    
    ax.xaxis.label.set_size(fontsize)
    ax.yaxis.label.set_size(fontsize)
    
    # Adjust layout
    plt.tight_layout()

def get_files(directory, pattern):
    return glob.glob(os.path.join(directory, pattern))


def prepare_atac_seq(file, quant_num=4, region_size=1000, interval_length=10):
    atacs = pd.read_table(file, header=None)
    atacs.columns=["chrom","start","end","name","signal"]
    atacs['quantile'] = pd.qcut(atacs["signal"], q=quant_num, labels=False)

    atac_regions = atacs[["chrom", "start", "end", "signal", "quantile"]].copy()
    atac_regions["name"] = "atac_" + atac_regions.index.astype(str)
    atac_regions.columns = ["chrom", "start", "end", "signal", "quantile", "name"]
    atac_regions["start_ext"] = (atac_regions["start"] + atac_regions["end"]) //2 - region_size
    atac_regions["end_ext"] = (atac_regions["start"] + atac_regions["end"]) //2 + region_size + interval_length

    expanded_regions = pd.DataFrame({
        'chrom': np.repeat(atac_regions['chrom'], (atac_regions['end_ext'] - atac_regions['start_ext']) // interval_length),
        'start': np.concatenate([np.arange(s, e, interval_length) for s, e in zip(atac_regions['start_ext'], atac_regions['end_ext'])]),
        'end': np.concatenate([np.arange(s, e, interval_length) + interval_length for s, e in zip(atac_regions['start_ext'], atac_regions['end_ext'])]),
        'name': np.repeat(atac_regions['name'], (atac_regions['end_ext'] - atac_regions['start_ext']) // interval_length),
        'midpoint': np.repeat(atac_regions['start_ext'] + (atac_regions['end_ext'] - atac_regions['start_ext']) // 2, (atac_regions['end_ext'] - atac_regions['start_ext']) // interval_length),
        'quantile': np.repeat(atac_regions['quantile'], (atac_regions['end_ext'] - atac_regions['start_ext']) // interval_length),
    }).reset_index(drop=True)

    expanded_regions["window"] = ((expanded_regions["start"] - expanded_regions["midpoint"]) // interval_length).astype(int) +1

    return expanded_regions


def calculate_dinucleotide_frequency_celegans(regions_df, genome_fasta_path, dinucleotide="TT"):
    """
    Calculate dinucleotide frequency for regions in C. elegans genome.
    
    Args:
        regions_df: DataFrame with chrom, start, end columns
        genome_fasta_path: Path to C. elegans genome FASTA file
    
    Returns:
        DataFrame with dinucleotide frequencies for each region
    """
    # Load C. elegans genome
    print("Loading C. elegans genome...")
    genome = {}
    for record in SeqIO.parse(genome_fasta_path, "fasta"):
        genome[record.id] = str(record.seq)
    
    print(f"Loaded {len(genome)} chromosomes")
    
    results = []
    
    for idx, row in regions_df.iterrows():
        chrom = row['chrom']
        start = row['start']
        end = row['end']
        
        # Get sequence for this region
        if chrom in genome:
            seq = genome[chrom][start:end+1].upper()
            
            # Calculate TT frequency
            dinucleotide_count = len(re.findall(f'(?={dinucleotide})', seq))
            # if dinucleotide_count > 0:
            #     print(seq, dinucleotide_count)
            dinucleotide_freq = dinucleotide_count / (len(seq) - 1) if len(seq) > 1 else 0
            
            results.append({
                'region_id': idx,
                'chrom': chrom,
                'start': start,
                'end': end,
                'window': row['window'],
                'quantile': row['quantile'],
                'sequence_length': len(seq),
                f'{dinucleotide}_count': dinucleotide_count,
                f'{dinucleotide}_frequency': dinucleotide_freq,
                'sequence': seq  # Optional: remove if you don't need the sequence
            })
        else:
            print(f"Warning: Chromosome {chrom} not found in genome")
            results.append({
                'region_id': idx,
                'chrom': chrom,
                'start': start,
                'end': end,
                'window': row['window'],
                'quantile': row['quantile'],
                'sequence_length': 0,
                f'{dinucleotide}_count': 0,
                f'{dinucleotide}_frequency': 0,
                'sequence': ''
            })
    
    return pd.DataFrame(results)


def calculate_rpkm(region_df, read_counts_column, total_reads):
    """
    Calculate RPKM values for the given region DataFrame.

    Parameters:
        region_df (DataFrame): A DataFrame with 'start' and 'end' columns defining region lengths.
        read_counts_column (str): The name of the column containing read counts for each region.
        total_reads (int): Total number of reads in the dataset (line count of the file).

    Returns:
        Series: A pandas Series with RPKM values.
    """
    # Region lengths in kilobases
    region_lengths_kb = (region_df["end"] - region_df["start"]) / 1000  # Length in kilobases

    # Total reads in millions
    total_reads_million = total_reads / 1_000_000  # Convert to millions

    # RPKM calculation
    rpkm_values = region_df[read_counts_column] / (region_lengths_kb * total_reads_million)
    return rpkm_values


def process_xr_ds(curr_file, expanded_regions, name):

        df = bf.read_table(curr_file, schema="bed")
        overlapped = bf.count_overlaps(expanded_regions[["chrom", "start", "end", "name", "quantile", "quantile_label", "window"]], df)

        overlapped.rename(columns={"count": name}, inplace=True)

        with open(curr_file, 'r') as f:
            total_reads = sum(1 for _ in f)

        if "_plus.bed" in curr_file:
            opposite_strand = curr_file.replace("_plus.bed", "_minus.bed")
        else:
            opposite_strand = curr_file.replace("_minus.bed", "_plus.bed")

        with open(opposite_strand, 'r') as f_opp:
            total_reads += sum(1 for _ in f_opp)

        overlapped[name] = calculate_rpkm(
            region_df=overlapped,
            read_counts_column=name,
            total_reads=total_reads
        )

        return overlapped


@timing_decorator()
def process_file(file, expanded_regions, output, interval=10, target_period=16, limit=26):

    quantile_labels = {0: 'Low', 1: 'Mid-Low', 2: 'Mid-High', 3: 'High'}
    expanded_regions['quantile_label'] =expanded_regions['quantile'].map(quantile_labels)
    quantile_order = ['High', 'Mid-High', 'Mid-Low', 'Low']
    custom_palette = {
    'High': '#2b2d42',
    'Mid-High': '#8d99ae',
    'Mid-Low': '#f49cbb',
    'Low': '#ef233c'
    }
    fontsize = 18
    name = os.path.basename(file).replace(".bed", "")
    signal_info = pd.DataFrame(columns=["name", "quantile", "snr", "p_value", "period", "signal"])

    if "csb-1_L1" in name:
        ylow, yhigh = 0, 110
    elif "Damage" in name:
        ylow, yhigh = 6, 14
    elif "WT_L1" in name or "xpc-1_L1" in name:
        ylow, yhigh = 10, 100
    elif "NHF1" in name or "CS-B" in name:
        ylow, yhigh = 0, 6

    for ext in ["", "_sim"]:
        curr_name = name + ext
        curr_file = file + ext + "_plus.bed"
        overlapped_plus = process_xr_ds(curr_file, expanded_regions, curr_name)
        overlapped_plus.loc[:, "window"] *= interval * -1

        curr_file = file + ext + "_minus.bed"
        overlapped_minus = process_xr_ds(curr_file, expanded_regions, curr_name)
        overlapped_minus.loc[:, "window"] *= interval

        overlapped = pd.concat([overlapped_plus, overlapped_minus], ignore_index=True)

        mapped_agg = overlapped[["quantile_label", "window", curr_name]].groupby(["quantile_label", "window"]).agg("mean").reset_index()
        mapped_agg.to_csv(f"{output}/{name}{ext}_mapped.csv")

        for quant in quantile_labels:
            if "xpc" in name:
                if interval == 10:
                    ds_start = 25
                elif interval == 2:
                    ds_start = 10
                signal_rate_ds = overlapped[(overlapped["quantile"]==quant) & (overlapped["window"]>ds_start*interval)][["window",curr_name]].groupby("window").mean().reset_index()[curr_name].values
                signal_rate_ds_0 = signal_rate_ds - np.mean(signal_rate_ds)
                observed_snr_ds, p_value_ds, permuted_snrs, signal_idxs_ds, signals_ds = snr.permutation_snr(signal_rate_ds_0, n_permutations=1000, limit=limit, target_period=target_period, xpc=True)
                maximum_period_ds = signal_idxs_ds[np.argmax(signals_ds)] * interval
                new_rows_ds = pd.DataFrame({"name": f"{curr_name}_ds", "quantile": quant, "snr": observed_snr_ds, "p_value": p_value_ds, "period": signal_idxs_ds, "signal": signals_ds})

                signal_rate_us = overlapped[(overlapped["quantile"]==quant) & (overlapped["window"]<-1*ds_start*interval)][["window",curr_name]].groupby("window").mean().reset_index()[curr_name].values
                signal_rate_us_0 = signal_rate_us - np.mean(signal_rate_us)
                observed_snr_us, p_value_us, permuted_snrs, signal_idxs_us, signals_us = snr.permutation_snr(signal_rate_us_0, n_permutations=1000, limit=limit, target_period=target_period, xpc=True)
                maximum_period_us = signal_idxs_us[np.argmax(signals_us)] * interval

                signals_us = signals_us / np.mean(signals_us + signals_ds)
                signals_ds = signals_ds / np.mean(signals_us + signals_ds)

                new_rows_ds = pd.DataFrame({"name": f"{curr_name}_ds", "quantile": quant, "snr": observed_snr_ds, "p_value": p_value_ds, "period": signal_idxs_ds, "signal": signals_ds})
                new_rows_us = pd.DataFrame({"name": f"{curr_name}_us", "quantile": quant, "snr": observed_snr_us, "p_value": p_value_us, "period": signal_idxs_us, "signal": signals_us})
                new_rows = pd.concat([new_rows_ds, new_rows_us], ignore_index=True)
            else:
                signal_rate = overlapped[overlapped["quantile"]==quant][["window",curr_name]].groupby("window").mean().reset_index()[curr_name].values
                signal_rate_0 = signal_rate - np.mean(signal_rate)
                observed_snr, p_value, permuted_snrs, signal_idxs, signals = snr.permutation_snr(signal_rate_0, n_permutations=1000, limit=limit, target_period=target_period)
                maximum_period = signal_idxs[np.argmax(signals)] * interval
                new_rows = pd.DataFrame({"name": curr_name, "quantile": quant, "snr": observed_snr, "p_value": p_value, "period": signal_idxs, "signal": signals})

            if signal_info.empty:
                signal_info = new_rows
            else:
                signal_info = pd.concat([signal_info, new_rows], ignore_index=False)

            if quant == max(overlapped["quantile"]):

                fig, ax = plt.subplots(figsize=(12, 6))
                fig.patch.set_facecolor('white')
                sns.lineplot(x='window', y=curr_name, data=overlapped, hue='quantile_label', 
                    hue_order=quantile_order[::-1], palette=custom_palette, ax=ax)

                ax.tick_params(axis='both', labelsize=fontsize)


                plt.xlabel("Distance from ATAC-seq center (bp)", fontsize=fontsize)
                plt.ylabel(ext.replace(".bed","").replace("_"," ").title() + " RPKM", fontsize=fontsize)
                plt.ylim(ylow, yhigh)
                plt.title(name.replace("_"," "))
                if output == "here":
                    plt.show()
                else:
                    plt.savefig(f"{output}/{name}{ext}.svg")
                plt.close()

    signal_info.to_csv(f"{output}/{name}_snr.csv")


@timing_decorator()
def process_file_sim_norm(file, expanded_regions, output, interval=10, target_period=16, limit=26):

    quantile_labels = {0: 'Low', 1: 'Mid-Low', 2: 'Mid-High', 3: 'High'}
    expanded_regions['quantile_label'] =expanded_regions['quantile'].map(quantile_labels)

    name = os.path.basename(file).replace(".bed", "")
    signal_info = pd.DataFrame(columns=["name", "quantile", "snr", "p_value", "period", "signal"])

    def process_read(expanded_regions, interval, name, file, ext):
        curr_name = name + ext
        curr_file = file + ext + "_plus.bed"
        overlapped_plus = process_xr_ds(curr_file, expanded_regions, curr_name)
        overlapped_plus.loc[:, "window"] *= interval * -1
        curr_file = file + ext + "_minus.bed"
        overlapped_minus = process_xr_ds(curr_file, expanded_regions, curr_name)
        overlapped_minus.loc[:, "window"] *= interval
        overlapped = pd.concat([overlapped_plus, overlapped_minus], ignore_index=True)

        mapped_agg = overlapped[["quantile_label", "quantile", "window", curr_name]].groupby(["quantile_label", "quantile", "window"]).agg("mean").reset_index()

        return mapped_agg

    overlapped = process_read(expanded_regions, interval, name, file, "")
    overlapped_sim = process_read(expanded_regions, interval, name, file, "_sim")
    overlapped = pd.merge(overlapped, overlapped_sim, on=["quantile_label", "quantile", "window"], how="left")
    overlapped[f"{name}_sim_norm"] = overlapped[name] / overlapped[f"{name}_sim"]
    overlapped.drop(columns=[f"{name}_sim"], inplace=True)
    overlapped.drop(columns=[f"{name}"], inplace=True)

    overlapped.to_csv(f"{output}/{name}_mapped.csv")
    curr_name = f"{name}_sim_norm"

    for quant in quantile_labels:
        if "xpc" in name:
            if interval == 10:
                ds_start = 25
            elif interval == 2:
                ds_start = 10
            signal_rate_ds = overlapped[(overlapped["quantile"]==quant) & (overlapped["window"]>ds_start*interval)][curr_name].values
            signal_rate_ds_0 = signal_rate_ds - np.mean(signal_rate_ds)
            observed_snr_ds, p_value_ds, permuted_snrs, signal_idxs_ds, signals_ds = snr.permutation_snr(signal_rate_ds_0, n_permutations=1000, limit=limit, target_period=target_period, xpc=True)
            maximum_period_ds = signal_idxs_ds[np.argmax(signals_ds)] * interval
            new_rows_ds = pd.DataFrame({"name": f"{curr_name}_ds", "quantile": quant, "snr": observed_snr_ds, "p_value": p_value_ds, "period": signal_idxs_ds, "signal": signals_ds})

            signal_rate_us = overlapped[(overlapped["quantile"]==quant) & (overlapped["window"]<-1*ds_start*interval)][curr_name].values
            signal_rate_us_0 = signal_rate_us - np.mean(signal_rate_us)
            observed_snr_us, p_value_us, permuted_snrs, signal_idxs_us, signals_us = snr.permutation_snr(signal_rate_us_0, n_permutations=1000, limit=limit, target_period=target_period, xpc=True)
            maximum_period_us = signal_idxs_us[np.argmax(signals_us)] * interval

            signals_us = signals_us / np.mean(signals_us + signals_ds)
            signals_ds = signals_ds / np.mean(signals_us + signals_ds)

            new_rows_ds = pd.DataFrame({"name": f"{curr_name}_ds", "quantile": quant, "snr": observed_snr_ds, "p_value": p_value_ds, "period": signal_idxs_ds, "signal": signals_ds})
            new_rows_us = pd.DataFrame({"name": f"{curr_name}_us", "quantile": quant, "snr": observed_snr_us, "p_value": p_value_us, "period": signal_idxs_us, "signal": signals_us})
            new_rows = pd.concat([new_rows_ds, new_rows_us], ignore_index=True)
        else:
            signal_rate = overlapped[overlapped["quantile"]==quant][curr_name].values
            signal_rate_0 = signal_rate - np.mean(signal_rate)
            observed_snr, p_value, permuted_snrs, signal_idxs, signals = snr.permutation_snr(signal_rate_0, n_permutations=1000, limit=limit, target_period=target_period)
            maximum_period = signal_idxs[np.argmax(signals)] * interval
            new_rows = pd.DataFrame({"name": curr_name, "quantile": quant, "snr": observed_snr, "p_value": p_value, "period": signal_idxs, "signal": signals})

        if signal_info.empty:
            signal_info = new_rows
        else:
            signal_info = pd.concat([signal_info, new_rows], ignore_index=False)

    signal_info.to_csv(f"{output}/{name}_snr.csv")


@timing_decorator()
def process_file_double_norm(file, expanded_regions, output, interval=10, target_period=16, limit=26):

    quantile_labels = {0: 'Low', 1: 'Mid-Low', 2: 'Mid-High', 3: 'High'}
    expanded_regions['quantile_label'] =expanded_regions['quantile'].map(quantile_labels)

    name = os.path.basename(file).replace(".bed", "")
    signal_info = pd.DataFrame(columns=["name", "quantile", "snr", "p_value", "period", "signal"])

    def process_read(expanded_regions, interval, name, file, ext):
        curr_name = name + ext
        curr_file = file + ext + "_plus.bed"
        overlapped_plus = process_xr_ds(curr_file, expanded_regions, curr_name)
        overlapped_plus.loc[:, "window"] *= interval * -1
        curr_file = file + ext + "_minus.bed"
        overlapped_minus = process_xr_ds(curr_file, expanded_regions, curr_name)
        overlapped_minus.loc[:, "window"] *= interval
        overlapped = pd.concat([overlapped_plus, overlapped_minus], ignore_index=True)

        mapped_agg = overlapped[["quantile_label", "quantile", "window", curr_name]].groupby(["quantile_label", "quantile", "window"]).agg("mean").reset_index()

        return mapped_agg

    overlapped_xr = process_read(expanded_regions, interval, name, file, "")
    overlapped_xr_sim = process_read(expanded_regions, interval, name, file, "_sim")
    overlapped_xr = pd.merge(overlapped_xr, overlapped_xr_sim, on=["quantile_label", "quantile", "window"], how="left")
    overlapped_xr[f"{name}_sim_norm"] = overlapped_xr[name] / overlapped_xr[f"{name}_sim"]
    overlapped_xr.drop(columns=[f"{name}_sim"], inplace=True)
    overlapped_xr.drop(columns=[f"{name}"], inplace=True)

    if "5min" in name or "1h" in name:
        name_ds = name.replace("5min", "0h").replace("1h", "0h").replace("_xr_", "_ds_")
        file_ds = file.replace("_xr_", "_ds_").replace("5min", "0h").replace("1h", "0h")
    else:
        name_ds = name.replace("_xr_", "_ds_")
        file_ds = file.replace("_xr_", "_ds_")

    if "xpc" in name or "csb" in name:
        name_ds = name_ds.replace("xpc", "wt").replace("csb", "wt").replace("48h", "0h").replace("8h", "0h").replace("24h", "0h")
        file_ds = file_ds.replace("xpc", "wt").replace("csb", "wt").replace("48h", "0h").replace("8h", "0h").replace("24h", "0h")

    overlapped_ds = process_read(expanded_regions, interval, name_ds, file_ds, "")
    overlapped_ds_sim = process_read(expanded_regions, interval, name_ds, file_ds, "_sim")
    overlapped_ds = pd.merge(overlapped_ds, overlapped_ds_sim, on=["quantile_label", "quantile", "window"], how="left")
    overlapped_ds[f"{name_ds}_sim_norm"] = overlapped_ds[name_ds] / overlapped_ds[f"{name_ds}_sim"]
    overlapped_ds.drop(columns=[f"{name_ds}_sim"], inplace=True)
    overlapped_ds.drop(columns=[f"{name_ds}"], inplace=True)

    overlapped = pd.merge(overlapped_xr, overlapped_ds, on=["quantile_label", "quantile", "window"], how="left")
    overlapped[f"{name}_double_norm"] = overlapped[f"{name}_sim_norm"] / overlapped[f"{name_ds}_sim_norm"]
    overlapped.drop(columns=[f"{name}_sim_norm"], inplace=True)
    overlapped.drop(columns=[f"{name_ds}_sim_norm"], inplace=True)

    overlapped.to_csv(f"{output}/{name}_mapped.csv")
    curr_name = f"{name}_double_norm"

    for quant in quantile_labels:
        if "xpc" in name:
            if interval == 10:
                ds_start = 25
            elif interval == 2:
                ds_start = 10
            signal_rate_ds = overlapped[(overlapped["quantile"]==quant) & (overlapped["window"]>ds_start*interval)][curr_name].values
            signal_rate_ds_0 = signal_rate_ds - np.mean(signal_rate_ds)
            observed_snr_ds, p_value_ds, permuted_snrs, signal_idxs_ds, signals_ds = snr.permutation_snr(signal_rate_ds_0, n_permutations=1000, limit=limit, target_period=target_period, xpc=True)
            maximum_period_ds = signal_idxs_ds[np.argmax(signals_ds)] * interval
            new_rows_ds = pd.DataFrame({"name": f"{curr_name}_ds", "quantile": quant, "snr": observed_snr_ds, "p_value": p_value_ds, "period": signal_idxs_ds, "signal": signals_ds})

            signal_rate_us = overlapped[(overlapped["quantile"]==quant) & (overlapped["window"]<-1*ds_start*interval)][curr_name].values
            signal_rate_us_0 = signal_rate_us - np.mean(signal_rate_us)
            observed_snr_us, p_value_us, permuted_snrs, signal_idxs_us, signals_us = snr.permutation_snr(signal_rate_us_0, n_permutations=1000, limit=limit, target_period=target_period, xpc=True)
            maximum_period_us = signal_idxs_us[np.argmax(signals_us)] * interval

            signals_us = signals_us / np.mean(signals_us + signals_ds)
            signals_ds = signals_ds / np.mean(signals_us + signals_ds)

            new_rows_ds = pd.DataFrame({"name": f"{curr_name}_ds", "quantile": quant, "snr": observed_snr_ds, "p_value": p_value_ds, "period": signal_idxs_ds, "signal": signals_ds})
            new_rows_us = pd.DataFrame({"name": f"{curr_name}_us", "quantile": quant, "snr": observed_snr_us, "p_value": p_value_us, "period": signal_idxs_us, "signal": signals_us})
            new_rows = pd.concat([new_rows_ds, new_rows_us], ignore_index=True)
        else:
            signal_rate = overlapped[overlapped["quantile"]==quant][curr_name].values
            signal_rate_0 = signal_rate - np.mean(signal_rate)
            observed_snr, p_value, permuted_snrs, signal_idxs, signals = snr.permutation_snr(signal_rate_0, n_permutations=1000, limit=limit, target_period=target_period)
            maximum_period = signal_idxs[np.argmax(signals)] * interval
            new_rows = pd.DataFrame({"name": curr_name, "quantile": quant, "snr": observed_snr, "p_value": p_value, "period": signal_idxs, "signal": signals})

        if signal_info.empty:
            signal_info = new_rows
        else:
            signal_info = pd.concat([signal_info, new_rows], ignore_index=False)

    signal_info.to_csv(f"{output}/{name}_snr.csv")


def process_power_spectrum_snr_plots(file_paths, output_folder, interval, target_period, norm=False):

    def combine_snr(file_paths):
        snr_all = pd.concat([pd.read_csv(file, index_col=0) for file in file_paths]).reset_index(drop=True)
        snr_all["period"] *=interval
        name_labels = {name: snr.get_time_label(name) for name in snr_all['name'].unique()}
        snr_all['name_label'] =snr_all['name'].map(name_labels)
        snr_all["name"] = snr_all["name"].apply(lambda x: x if x.endswith(("sim", "sim_ds", "sim_us")) else x + "_obs")

        return snr_all


    def periodicity_plot(snr_target, target_list, norm):

        if norm == False and not "dyad" in output_folder:
            ymax = 15
        else:
            ymax = 12

        snr.plot_periodicity(snr_target[snr_target["quantile"]==min(snr_target["quantile"])],
            output_folder=output_folder, ymax=ymax, name="_".join(target_list + ["min_quantile"]),
            vline=target_period*interval, xticks=x_list)
        snr.plot_periodicity(snr_target[snr_target["quantile"]==max(snr_target["quantile"])],
            output_folder=output_folder, ymax=ymax, name="_".join(target_list + ["max_quantile"]),
            vline=target_period*interval, xticks=x_list)


    def plot_all(snr_all, target_list, norm):
        snr_target = snr_all.copy()

        if norm == True and "sim_norm" in output_folder and "dyad" in output_folder and "xr" in target_list:
            ymax = 40
        elif "dyad" in output_folder and "xr" in target_list:
            ymax = 30
        elif "dyad" in output_folder:
            ymax = 130
        elif norm == False and not "dyad" in output_folder:
            ymax = 35
        else:
            ymax = 25

        for target in target_list:
            snr_target = snr_target[snr_target["name"].str.contains(target, na=False)]
        if len(snr_target) != 0:
            periodicity_plot(snr_target, target_list, norm)

        snr_target = snr_target[snr_target["period"] == target_period*interval]

        if len(snr_target) != 0:
            if norm == False and "sim" in target_list:
                if "ds" == target_list[-1]:
                    sim_name = "sim_ds"
                elif "us" == target_list[-1]:
                    sim_name = "sim_us"
                else:
                    sim_name = "sim"

                snr_target_sim = snr_target[snr_target["name"].str.endswith(sim_name)]
                snr_target_sim.to_csv(f"{output_folder}/{'_'.join(target_list)}_snr_bf_plot.csv", index=False)
                snr.plot_snr(snr_target_sim, output_folder=output_folder, ymax=ymax, name="_".join(target_list))

            if "obs" in target_list or not ("sim" in target_list and "obs" in target_list):
                snr_target_obs = snr_target[snr_target["name"].str.endswith("obs")]
                snr_target_obs.to_csv(f"{output_folder}/{'_'.join(target_list)}_snr_bf_plot.csv", index=False)
                snr.plot_snr(snr_target_obs, output_folder=output_folder, ymax=ymax, name="_".join(target_list))

    ##################################### Main Function ############################################

    x_list = [i*interval for i in [5,16,24]]
    if interval ==2:
        x_list=[4,10,15]

    snr_all = combine_snr(file_paths)

    target_list_of_list = snr.possible_combinations(norm)
    filt_target_list_of_list = [
        target_list for target_list in target_list_of_list
        if any(all(target in sample for target in target_list[:len(target_list)-1])
               for sample in file_paths)
    ]
    for target_list in filt_target_list_of_list:

        if (("xpc" in target_list and "ds" in target_list) or
            ("csb" in target_list and "ds" in target_list)):
            continue
        else:
            if "xpc" in target_list:
                target_list_us = target_list.copy() + ["us"]
                plot_all(snr_all, target_list_us, norm)
                target_list_ds = target_list.copy() + ["ds"]
                plot_all(snr_all, target_list_ds, norm)
            else:
                plot_all(snr_all, target_list, norm)

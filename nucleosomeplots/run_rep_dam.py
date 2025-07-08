import os
import glob

import periodicity_helper_func as hf

from concurrent.futures import ProcessPoolExecutor


def get_files(directory, pattern):
    return glob.glob(os.path.join(directory, pattern))

def get_samples(xr_ds_dir, filt=None):
    files = get_files(xr_ds_dir, "*.bed")
    samples=[]
    for file in files:
        file = file.replace("_sim_plus.bed", "").replace("_sim_minus.bed", "").replace("_plus.bed", "").replace("_minus.bed", "")
        if file not in samples:
            samples.append(file)

    filt_samples = [sample for sample in samples if filt in sample]
    return filt_samples

info = {
    "samples": ["wt", "csb", "xpc"],
    "region": ["atac", "nucleoatac_dyad"],
    "region_size": [1000, 80],
    "interval": [10, 2],
    "target_period": [16, 5],
    "period_limit": [26, 10]
}


sample_dir = "" # Directory of the xr and ds samples (e.g. workdir/merged_bed_exact_dam_site)
output_dir = "" # Directory of the output

for sample in info["samples"]:
    samples = get_samples(sample_dir, filt=sample)
    for i,region in enumerate(info["region"]):
        region_size = info["region_size"][i]
        interval = info["interval"][i]
        target_period = info["target_period"][i]
        limit = info["period_limit"][i]
        output_folder = f"{output_dir}/{sample}_{region}"
        os.makedirs(output_folder, exist_ok=True)

        quant_num = 4

        output_folder_sim_norm = output_folder + "_sim_norm"
        os.makedirs(output_folder_sim_norm, exist_ok=True)

        output_folder_double_norm = output_folder + "_double_norm"
        os.makedirs(output_folder_double_norm, exist_ok=True)

        samples_double_norm = [x for x in samples if "xr" in x]

        expanded_regions = hf.prepare_atac_seq(quant_num=quant_num, file=region, region_size=region_size, interval_length=interval)

        with ProcessPoolExecutor(max_workers=18) as executor:
            executor.map(hf.process_file, samples, [expanded_regions]*len(samples), [output_folder]*len(samples), [interval]*len(samples), [target_period]*len(samples), [limit]*len(samples))

        with ProcessPoolExecutor(max_workers=18) as executor:
            executor.map(hf.process_file_sim_norm, samples, [expanded_regions]*len(samples), [output_folder_sim_norm]*len(samples), [interval]*len(samples), [target_period]*len(samples), [limit]*len(samples))

        with ProcessPoolExecutor(max_workers=18) as executor:
            executor.map(hf.process_file_double_norm, samples_double_norm, [expanded_regions]*len(samples), [output_folder_double_norm]*len(samples), [interval]*len(samples), [target_period]*len(samples), [limit]*len(samples))

        file_paths = hf.get_files(output_folder, "*_snr.csv")
        file_paths = [x for x in file_paths if ("nak" not in x and "noUV" not in x)]
        hf.process_power_spectrum_snr_plots(file_paths, output_folder, interval, target_period)

        file_paths_sim_norm = hf.get_files(output_folder_sim_norm, "*_snr.csv")
        file_paths_sim_norm = [x for x in file_paths_sim_norm if ("nak" not in x and "noUV" not in x)]
        hf.process_power_spectrum_snr_plots(file_paths_sim_norm, output_folder_sim_norm, interval, target_period, norm=True)

        file_paths_double_norm = hf.get_files(output_folder_double_norm, "*_snr.csv")
        file_paths_double_norm = [x for x in file_paths_double_norm if ("nak" not in x and "noUV" not in x)]
        file_paths_double_norm = [x for x in file_paths_double_norm if "xr" in x]
        hf.process_power_spectrum_snr_plots(file_paths_double_norm, output_folder_double_norm, interval, target_period, norm=True)

        print("All samples are run successfully.")

        with open(f'{output_folder}/successful.txt', 'w') as f:
            f.write("All samples are run successfully.")

        with open(f'{output_folder_sim_norm}/successful.txt', 'w') as f:
            f.write("All samples are run successfully.")

        with open(f'{output_folder_double_norm}/successful.txt', 'w') as f:
            f.write("All samples are run successfully.")

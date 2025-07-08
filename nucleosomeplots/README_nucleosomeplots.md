## nucleosome_plots

This directory contains scripts and helper functions for analyzing and visualizing nucleosome positioning, periodicity, and DNA damage/repair signals in the context of ATAC-seq and NucleoATAC data. The pipeline supports quantification, normalization, periodicity analysis, and advanced plotting of nucleosome and DNA damage/repair features.

### 📜 Script Overview

- **`run_rep_dam.py`**  
  Main analysis pipeline for nucleosome. Runs region expansion, quantile binning, 
  and parallelized processing of DNA damage/repair signals across samples and regions.
  Generates normalized and double-normalized outputs, and computes power spectrum and SNR
  (signal-to-noise ratio) plots for periodicity analysis.

- **`periodicity_helper_func.py`**  
  Contains core helper functions for region expansion, dinucleotide frequency calculation, 
  quantile binning, RPKM calculation.

- **`snr.py`**  
  Implements signal processing and statistical analysis for periodicity, including power spectrum
  calculation, SNR computation, permutation testing, and advanced plotting
  (periodograms, SNR distributions, etc.).

- **`normalized.py`**  
  Aggregates, normalizes, and visualizes time-series data for nucleosome and ATAC-seq signals. Includes functions for rolling average smoothing, bar arrangement for dyad/nucleosome plots, and time-course plotting.

- **`tt_counts_dyad.py`**  
  Calculates and visualizes dinucleotide (TT/AA) frequencies around nucleosome dyads, merges and processes mapped dyad data, and generates comparative plots for simulated and real DNA damage signals.

### ⚙️ Requirements

- Python 3.x
- `pandas`, `numpy`
- `matplotlib`, `seaborn`
- `bioframe`
- `scipy`
- `statsmodels`
- `biopython`

### 📦 Outputs

- Aggregated and normalized CSV files for dyad and ATAC-seq signals
- Power spectrum and SNR plots for periodicity analysis
- Time-course and comparative plots for DNA damage/repair and nucleosome features
- Dinucleotide frequency statistics and plots

### 🧩 Usage

1. **Prepare input data**:  
   - Process ATAC-seq/NucleoATAC data and DNA damage/repair BED files as required.
2. **Configure paths**:  
   - Set input/output directories in the scripts as needed.
3. **Run the pipeline**:  
   - Use `run_rep_dam.py` to process all samples and generate outputs.
   - Use plotting scripts (`normalized.py`, `tt_counts_dyad.py`) for visualization and further analysis.

---

## Author

Cem Azgari, 2025 – Sabancı University | Adebali Lab

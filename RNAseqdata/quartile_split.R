# Load necessary library
library(data.table)

# Define the input/output directory
data_dir <- commandArgs(trailingOnly = TRUE)[1]

# Find all BED files to process in the data directory
bed_files <- list.files(path = data_dir, pattern = "_500bpaway.bed", full.names = TRUE)

# Process each file
for (bed_file in bed_files) {
  # Read the BED file
  bed_data <- fread(bed_file, header = FALSE, sep = "\t")

  # Assuming TPM values are in the 5th column
  tpm_values <- bed_data[[5]]

  # Calculate the quartile thresholds
  quartiles <- quantile(tpm_values, probs = 0:4/4, na.rm = TRUE)

  # Split data into quartiles
  quartile_names <- c("0_25", "25_50", "50_75", "75_100")
  for (i in 1:4) {
    quartile_data <- bed_data[tpm_values > quartiles[i] & tpm_values <= quartiles[i + 1]]

    # Create the output filename
    base_name <- gsub("_500bpaway.bed", "", basename(bed_file))
    file_name <- file.path(data_dir, paste0(base_name, "_", quartile_names[i], "_quartile.bed"))

    # Write each quartile to a separate file
    fwrite(quartile_data, file = file_name, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

cat("Quartile analysis completed successfully.\n")

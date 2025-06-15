# Load required library
library(data.table)

# Define the directory paths
input_dir <- "/RNAseqdirectory"
output_dir <- "/RNAseqdirectory"

# List of genotypes to process
genotypes <- c("wt", "csb", "xpc")

# Function to calculate TPM
calculate_tpm <- function(counts, gene_length) {
  tpm <- (counts / gene_length) * 1e3
  tpm <- tpm / sum(tpm) * 1e6
  return(tpm)
}

# Set threshold for filtering
threshold <- 0.1  # Adjust as needed

# Load gene data from the BED file
gene_data <- fread("/work/users/c/a/cansuk/2025JanuaryCeleganspaperworkflow/wormbaseWS295/c_elegansPRJNA13758_WS295_annotations.bed",
                   header = FALSE, sep = "\t",
                   col.names = c("chr", "start", "end", "strand", "gene_ids", "gene_name", "gene_type"))

# Debugging: Print the first few rows of the gene data
cat("Gene data preview:\n")
print(head(gene_data))

# Process each genotype
for (genotype in genotypes) {
  # Read replicate count data
  replicate1 <- fread(file.path(input_dir, paste0(genotype, "L1R1_count.txt")))
  replicate2 <- fread(file.path(input_dir, paste0(genotype, "L1R2_count.txt")))
  
  # Debugging: Print file details
  cat("\nProcessing genotype:", genotype, "\n")
  cat("Replicate1 dimensions:", dim(replicate1), "\n")
  cat("Replicate2 dimensions:", dim(replicate2), "\n")
  cat("Columns in replicate1:\n", colnames(replicate1), "\n")
  cat("Columns in replicate2:\n", colnames(replicate2), "\n")
  
  # Dynamically match the correct column for counts
  counts1_column <- grep(paste0(genotype, "L1R1_sorted.bam"), colnames(replicate1), value = TRUE)
  counts2_column <- grep(paste0(genotype, "L1R2_sorted.bam"), colnames(replicate2), value = TRUE)
  
  # Debugging: Print identified column names
  cat("Counts1 column found:", counts1_column, "\n")
  cat("Counts2 column found:", counts2_column, "\n")
  
  # Extract counts using the identified column names
  counts1 <- replicate1[[counts1_column]]
  counts2 <- replicate2[[counts2_column]]
  
  # Extract required information
  gene_ids <- replicate1$Geneid
  gene_length <- replicate1$Length
  
  # Create count matrix
  count_matrix <- data.frame(gene_ids, replicate1 = counts1, replicate2 = counts2)
  rownames(count_matrix) <- count_matrix$gene_ids
  count_matrix <- count_matrix[, -1]
  
  # Calculate TPM for both replicates
  tpm1 <- calculate_tpm(count_matrix$replicate1, gene_length)
  tpm2 <- calculate_tpm(count_matrix$replicate2, gene_length)
  
  # Combine TPMs and calculate average
  combined_tpm <- data.frame(gene_ids, tpm1, tpm2)
  combined_tpm$average_tpm <- rowMeans(combined_tpm[, c("tpm1", "tpm2")])
  
  # Filter genes based on TPM threshold
  filtered_tpm <- combined_tpm[combined_tpm$average_tpm > threshold, ]
  
  # Merge filtered TPM data with gene data
  merged_data <- merge(filtered_tpm, gene_data, by = "gene_ids")
  
  # Convert to data.table
  merged_data <- as.data.table(merged_data)
  
  # Correctly swap 'start' and 'end' for entries where the strand is '-'
  merged_data[strand == "-", c("start", "end") := .(end, start)]
  
  # Write filtered dataset to a file
  output_file <- file.path(output_dir, paste0(genotype, "_filtered_gene_data.bed"))
  write.table(merged_data[, c("chr", "start", "end", "strand", "gene_ids", "gene_name", "average_tpm")],
              file = output_file, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = FALSE, 
              quote = FALSE)
}

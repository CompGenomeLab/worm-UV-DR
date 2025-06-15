# Load necessary libraries
library(ggplot2)
library(dplyr)

# Set paths for data and output
data_path <- "../output/genebodyintersect/xr"
#for damageseq: ../output/genebodyintersect/damage
output_path <- "../output/genebodyintersect/xr/tsntsviolinplots"

# Get a list of all files ending with *_readCount.txt in the data directory
sample_files <- list.files(data_path, pattern = "_readCount.txt$", full.names = TRUE)

# Iterate over each sample file
for (sample_file in sample_files) {
  # Extract the sample name from the file name
  sample_name <- sub("_readCount.txt$", "", basename(sample_file))
  
  # Read the total reads from the read count file
  total_reads <- as.numeric(readLines(sample_file))
  
  # Function to calculate RPKM
  calculate_rpkm <- function(bed_file) {
    df <- read.table(bed_file, header = FALSE, sep = "\t", 
                     col.names = c("chr", "start", "end", "name", "score", "strand", "reads"))
    df$length <- df$end - df$start
    df$rpkm <- (df$reads / (df$length / 1000)) / (total_reads / 1e6)
    return(df[, c("name", "rpkm")]) # Return a data frame with gene names and RPKM values
  }
  
  # Read and calculate RPKM for each BED file
  bed_files <- list.files(data_path, pattern = paste0(sample_name, "_.*_TS.bed$"), full.names = TRUE)
  nts_files <- list.files(data_path, pattern = paste0(sample_name, "_.*_NTS.bed$"), full.names = TRUE)
  
  # Ensure corresponding TS and NTS files match
  bed_pairs <- data.frame(
    TS = bed_files[order(basename(bed_files))],
    NTS = nts_files[order(basename(nts_files))]
  )
  
  rpkm_list <- lapply(1:nrow(bed_pairs), function(i) {
    ts_rpkm <- calculate_rpkm(bed_pairs$TS[i])
    nts_rpkm <- calculate_rpkm(bed_pairs$NTS[i])
    quartile <- sub(".*_(\\d+_\\d+)_.*", "\\1", basename(bed_pairs$TS[i]))
    data.frame(
      Quartile = quartile,
      Strand = c(rep("Transcribed", nrow(ts_rpkm)), rep("Non-Transcribed", nrow(nts_rpkm))),
      RPKM = c(ts_rpkm$rpkm, nts_rpkm$rpkm)
    )
  })
  
  # Combine all RPKM data
  data_combined <- do.call(rbind, rpkm_list)
  
  # Factorize Quartile and Strand columns
  data_combined$Quartile <- factor(data_combined$Quartile, 
                                   levels = c("0_25", "25_50", "50_75", "75_100"),
                                   labels = c("0-25%", "25-50%", "50-75%", "75-100%"))
  data_combined$Strand <- factor(data_combined$Strand, levels = c("Transcribed", "Non-Transcribed"))
  
  # Perform paired t-test for each quartile
  quartiles <- unique(data_combined$Quartile)
  significance_results <- data.frame(Quartile = quartiles, PValue = NA)

  for (quartile in quartiles) {
    ts_values <- data_combined %>%
      filter(Quartile == quartile, Strand == "Transcribed") %>%
      pull(RPKM)
    nts_values <- data_combined %>%
      filter(Quartile == quartile, Strand == "Non-Transcribed") %>%
      pull(RPKM)
    
if (length(ts_values) > 0 && length(nts_values) > 0 && length(ts_values) == length(nts_values)) {
    test_result <- wilcox.test(ts_values, nts_values, paired = TRUE, exact = FALSE)
    significance_results[significance_results$Quartile == quartile, "PValue"] <- test_result$p.value
  } else {
    significance_results[significance_results$Quartile == quartile, "PValue"] <- NA
  }
  }

  # Add significance stars
  significance_results$Stars <- ifelse(significance_results$PValue < 0.001, "***",
                              ifelse(significance_results$PValue < 0.01, "**",
                              ifelse(significance_results$PValue < 0.05, "*", "")))
if (startsWith(sample_name, "DS_wt")) {
  strand_colors <- c("Transcribed" = "#7092a5", "Non-Transcribed" = "#6badc9")
} else {
  strand_colors <- c("Transcribed" = "#8b74a6", "Non-Transcribed" = "#c28b9a")
}
  # Generate violin plot with significance stars
  p <- ggplot(data_combined, aes(x = Quartile, y = RPKM, fill = Strand)) +
    geom_violin(trim = TRUE, color = NA, alpha = 0.8, position = position_dodge(0.75)) +
    scale_fill_manual(values = strand_colors) +
    scale_y_continuous(limits = c(0, 15)) +
    labs(title = paste0("RPKM for Transcribed and Non-Transcribed Strands Across Expression Quartiles\n", sample_name),
         x = "Expression Quartile", y = "RPKM") +
    geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.75), alpha = 0.6) +
    theme_classic() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.3), # Faint horizontal lines
      panel.grid.major.x = element_blank(),  # No vertical lines
      panel.border = element_rect(color = "black", fill = NA, size = 1), # Add frame
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 16), # Larger axis titles
      axis.text = element_text(size = 14)  # Larger axis text
    ) +
    geom_text(data = significance_results, aes(x = Quartile, y = 13, label = Stars), 
              inherit.aes = FALSE, size = 5)

  # Save the violin plot with smaller dimensions
  ggsave(file.path(output_path, paste0(sample_name, "_tsntsviolinwilcoxon.pdf")), 
         plot = p, width = 6, height = 4, device = "pdf")  # Smaller dimensions
}

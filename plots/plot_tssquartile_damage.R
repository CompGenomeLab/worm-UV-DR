library(ggplot2)
library(dplyr)

setwd('../output/TSSquartile/damage')

# Function to calculate RPKM
calculate_rpkm <- function(read_counts, total_reads_million) {
  gene_length_kb_per_bin <- 0.01  # 50 nucleotides = 0.05 kb
  (read_counts / gene_length_kb_per_bin) / total_reads_million
}

# Function to prepare plot data
prepare_plot_data <- function(df, total_reads_million) {
  df <- df %>%
    mutate(RPKM = calculate_rpkm(V8, total_reads_million)) %>%
    group_by(V7) %>%
    summarise(mean = mean(RPKM),
              std = sd(RPKM),
              count = n(),
              sem = std / sqrt(count),
              ci95 = 1.96 * sem)
  return(df)
}

# Function to create graph
create_graph <- function(sample_name, strand) {
  # Set Y-axis limits to 0-10 for all samples
  y_limits <- c(0, 25)
  
  # Define quartile labels
  quartile_labels <- c("_0_25", "_25_50", "_50_75", "_75_100")
  colors <- c("black", "blue", "green", "red")
  
  plot_data <- list()
  
  for (quartile in quartile_labels) {
    file_name <- paste0(sample_name, "_TSS_wt", quartile, "_quartile_final_150bins_", strand, ".bed")
    total_reads_file <- paste0(sample_name, "_readCount.txt")
    if (!file.exists(total_reads_file) || !file.exists(file_name)) {
      message("Skipping missing file for sample: ", sample_name, " - ", quartile)
      next
    }
    total_reads_million <- as.numeric(readLines(total_reads_file, n = 1)) / 1e6
    df <- read.table(file_name, header = FALSE, sep = "\t")
    plot_data[[quartile]] <- prepare_plot_data(df, total_reads_million)
  }
  
  # Combine data for plotting
  combined_data <- bind_rows(plot_data, .id = "quartile")
  combined_data$quartile <- factor(combined_data$quartile, levels = quartile_labels)
  
  # Check for empty data
  if (nrow(combined_data) == 0) {
    warning("No data available for sample: ", sample_name)
    return()
  }
  
  # Create plot
  p <- ggplot(combined_data, aes(x = as.numeric(V7), y = mean, color = quartile, fill = quartile)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = mean - ci95, ymax = mean + ci95), alpha = 0.2, color = NA) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_x_continuous(breaks = c(1, 51, 150), labels = c("-0.5kb", "TSS", "+1kb")) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    labs(x = NULL, y = "RPKM", title = paste(sample_name, "- TSS (", strand, ")", sep = "")) +
    theme_minimal(base_size = 15) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Save plot
  pdf_filename <- paste0(sample_name, "_TSS_", strand, ".pdf")
  ggsave(pdf_filename, plot = p, width = 6, height = 5)
  message("Graph saved as ", pdf_filename)
}

# Main: Process all samples
read_count_files <- list.files(pattern = "_readCount.txt$")
sample_names <- gsub("_readCount.txt$", "", read_count_files)

for (sample_name in sample_names) {
  message("Processing sample: ", sample_name)
  create_graph(sample_name, "NTS")  # Process NTS strand
  create_graph(sample_name, "TS")   # Process TS strand
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos = "https://cloud.r-project.org/")
}
library(pheatmap)
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer", repos = "https://cloud.r-project.org/")
}
library(RColorBrewer)

# Define a custom pastel palette
pastel_palette <-colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)


# Set working directory
setwd("../output/genelistintersect/2kbwindow")

# Output directory
plots_dir <- file.path(getwd(), "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir)

# Read count function
read_total_reads <- function(ts_file) {

  base <- sub("_2kbwindow_rpkm\\.bed$", "", basename(ts_file))
  readcount_file <- paste0(base, "_readCount.txt")
  fullpath <- file.path(getwd(), readcount_file)

  if (!file.exists(fullpath)) {
    warning("Missing read count file: ", fullpath)
    return(0)
  }
  as.numeric(readLines(fullpath))
}

# Metadata parsing from filename
extract_info <- function(fname) {
  parts <- strsplit(basename(fname), "_")[[1]]
  list(
    strain = parts[1],
    timepoint = parts[2],
    damage = parts[3],
    file = fname
  )
}
# Desired order
strain_order <- c("wt", "csb", "xpc")
damage_order <- c("64", "CPD")
timepoint_order <- c("5min", "1h", "8h", "24h", "48h")

# File list
rpkm_files <- list.files(pattern = "_2kbwindow_rpkm.bed$", full.names = TRUE)

# Build metadata
metadata <- do.call(rbind, lapply(rpkm_files, function(f) {
  info <- extract_info(f)
  data.frame(
    strain = info$strain,
    timepoint = info$timepoint,
    damage = info$damage,
    file = f,
    stringsAsFactors = FALSE
  )
}))


metadata$strain <- factor(metadata$strain, levels = strain_order)
metadata$damage <- factor(metadata$damage, levels = damage_order)
metadata$timepoint <- factor(metadata$timepoint, levels = timepoint_order)
metadata <- metadata[order(metadata$strain, metadata$damage, metadata$timepoint), ]
print(metadata$file_ts)
print(metadata$file_nts)


# Helper functions
read_bed_col <- function(file, colnum) {
  df <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  df[[colnum]]
}

normalize_to_rpm <- function(raw_counts, total_reads) {
  if (total_reads == 0) return(rep(0, length(raw_counts)))
  (raw_counts / total_reads) * 1e6
}

# Build RPKM matrix (from column 8)
combined_matrix_rpkm <- do.call(cbind, lapply(metadata$file, function(f) {
  read_bed_col(f, 5)
}))

# Optional: RPM matrix (from column 4, normalized by read count file)
read_total_reads <- function(file_path) {
  base <- sub("_2kbwindow_rpkm\\.bed$", "", basename(file_path))
  readcount_file <- paste0(base, "_readCount.txt")
  fullpath <- file.path(getwd(), readcount_file)
  if (!file.exists(fullpath)) {
    warning("Missing read count file: ", fullpath)
    return(0)
  }
  as.numeric(readLines(fullpath))
}

combined_matrix_rpm <- do.call(cbind, lapply(metadata$file, function(f) {
  raw_counts <- read_bed_col(f, 4)
  total_reads <- read_total_reads(f)
  normalize_to_rpm(raw_counts, total_reads)
}))

# Add row/column names
rownames(combined_matrix_rpkm) <- paste0("Region_", seq_len(nrow(combined_matrix_rpkm)))
rownames(combined_matrix_rpm) <- rownames(combined_matrix_rpkm)

colnames_format <- paste(metadata$strain, metadata$timepoint, metadata$damage, sep = "_")
colnames(combined_matrix_rpkm) <- colnames_format
colnames(combined_matrix_rpm) <- colnames_format

# Save matrices if needed
save(combined_matrix_rpkm, combined_matrix_rpm, file = file.path(plots_dir, "merged_combined_matrices.rda"))

# Clustering and heatmap plotting
do_hclust <- function(data_matrix, filename) {
  cor_matrix <- cor(data_matrix, method = "spearman", use = "complete.obs")
  dist_matrix <- as.dist(1 - cor_matrix)
  hclust_res <- hclust(dist_matrix, method = "ward.D2")
  pdf(file = file.path(plots_dir, filename), width = 13, height = 12)
  pheatmap(
    cor_matrix,
    color = pastel_palette,
    cluster_rows = hclust_res,
    cluster_cols = hclust_res,
    main = filename,
    fontsize_row = 14,
    fontsize_col = 14
  )
  dev.off()
}

# Generate heatmaps
do_hclust(combined_matrix_rpkm,         "hclust_2kbwindow_rpkm_.pdf")
do_hclust(log2(combined_matrix_rpkm+1), "hclust_2kbwindow_gene_rpkm_.pdf")
do_hclust(combined_matrix_rpm,          "hclust_2kbwindow_rpm_.pdf")
do_hclust(log2(combined_matrix_rpm+1),  "hclust_log2_2kbwindow_rpm_.pdf")

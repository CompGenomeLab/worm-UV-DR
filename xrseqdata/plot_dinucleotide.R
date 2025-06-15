# Set CRAN mirror explicitly
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install package if missing
if (!requireNamespace("wesanderson", quietly = TRUE)) {
  install.packages("wesanderson")
}
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(wesanderson)

# List all files with the specified pattern
# xlist <- list.files(pattern = "_get2ker.csv")
xlist <- c("wt_5min_64_merged_20senssim_24_bed_get2ker.csv", "wt_5min_64_merged_20senssim_24_get2ker.csv")
print(xlist)

# Read and process each file
for(filename in xlist) {
  data <- read.delim(filename)
  col_names <- c("base", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8", "8-9", "9-10",
                 "10-11", "11-12", "12-13", "13-14", "14-15", "15-16", "16-17", "17-18",
                 "18-19", "19-20", "20-21", "21-22", "22-23", "23-24")
  names(data) <- col_names
  custom_order <- col_names

  data <- data %>% filter(base %in% c('CC', 'CT', 'TC', 'TT'))
  data <- data %>% gather("position", "percent", 2:24)
  data$position <- factor(data$position, levels = custom_order)
  data <- data %>% mutate(percent = percent * 100, type = "genome")

  p <- ggplot(data, aes(position, percent, fill = base)) +
    geom_bar(stat = 'identity') +
    labs(x = "Dinucleotide Position", y = "% of Reads") +
    theme_tufte() +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          axis.line = element_line(colour = "black")) +
    scale_x_discrete(breaks = c("1-2", "3-4", "5-6", "7-8", "9-10","11-12",'13-14','15-16','17-18','19-20','21-22','23-24'),
                     labels = c("1-2", "", "5-6", "", "9-10","",'13-14','','17-18','','21-22','')) +
    scale_y_continuous(limits = c(0, 70), breaks = c(0, 20, 40, 60), labels = c("0", "20", "40", "60")) +
    coord_equal(ratio = 0.1) +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = wes_palette("AsteroidCity1", n = 4)) +
    labs(title = tools::file_path_sans_ext(basename(filename)))

  base_output_filename <- gsub("_get2ker.csv", "", filename)
  ggsave(paste0(base_output_filename, "_plot.png"), p)
  ggsave(paste0(base_output_filename, "_plot.pdf"), p, device = "pdf")
}

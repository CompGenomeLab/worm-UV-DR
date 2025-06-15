#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sample = args[1]

setwd('results')
library(ggplot2)
theme_set(theme_minimal())

# Read Length Distribution
rl_file = paste(sample, "_read_length_distribution.txt", sep = "")
df = as.data.frame(read.delim(rl_file))
colnames(df) = c("Length", "Count")
df$Percent = prop.table(df$Count) * 100
plot = ggplot(df, aes(x = Length, y = Percent)) +
  geom_bar(stat = "identity", fill = "deepskyblue3") +
  labs(y = "Frequency (%)", x = "Length (nt)", title = sample) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  coord_cartesian(clip = "off") +
  theme(plot.title = element_text(hjust = 0.5)) + theme_light()
ggsave(paste(sample, '_rld_graph.png', sep = ""), width = 5.54, height = 4.39)

# Monomer Analysis
mon_file = paste(sample, "_monomer_R_df.txt", sep = "")
df <- read.table(mon_file, header = FALSE)
colnames(df) <- c("Length", "Position", "Base", "Percentage")
df$Base <- factor(df$Base, levels = c("G", "A", "C", "T"), ordered = TRUE)
plot = ggplot(df, aes(x=Position, y=Percentage)) +
  geom_bar(stat = "identity", aes(fill = Base)) +
  labs(title = sample) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("G" = "purple4", "C" = "dodgerblue4", "A" = "green4", "T" = "orange")) +
  facet_grid(Length ~ .) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Nucleotide frequency (% of Total)") +
  theme(axis.text = element_text(size = 3),
        axis.title = element_text(size = 10),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        strip.text = element_text(size=5)) +
  theme_light()
ggsave(paste(sample, '_monomer_graph.png', sep = ""), width = 5.54, height = 4.39)

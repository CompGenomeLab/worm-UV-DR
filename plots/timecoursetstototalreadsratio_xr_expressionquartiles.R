library(ggplot2)
library(dplyr)
library(ggpubr)
library(purrr)

setwd('../output/genebodyintersect/xr')

strains <- c("csb", "wt", "xpc")
damages <- c("64", "CPD")
time_points <- c("5min", "1h", "8h", "24h", "48h")
quartiles <- c("0_25", "25_50", "50_75", "75_100")
quartile_labels <- c("0-25%", "25-50%", "50-75%", "75-100%")

process_bed_file <- function(bed_file) {
  df <- tryCatch({
    read.table(bed_file, header = FALSE, sep = "\t",
               col.names = c("chr", "start", "end", "name", "score", "strand", "reads"))
  }, error = function(e) {
    warning(paste("Error reading file:", bed_file, "-", e$message))
    return(NULL)
  })
  if (!is.null(df)) {
    df$length <- df$end - df$start
    df$reads[df$reads == 0] <- 0
  }
  return(df)
}

cat("Starting script...\n")

for (strain in strains) {
  for (damage in damages) {
    data_combined <- data.frame()
    
    for (time_point in time_points) {
      for (quartile in quartiles) {
        ts_file <- paste0(strain, "_", time_point, "_", damage, "_merged_", strain, "_", quartile, "_quartile_TS.bed")
        nts_file <- paste0(strain, "_", time_point, "_", damage, "_merged_", strain, "_", quartile, "_quartile_NTS.bed")
        
        if (!file.exists(ts_file) || !file.exists(nts_file)) {
          warning(paste("File not found:", ts_file, "or", nts_file))
          next
        }
        
        ts_df <- process_bed_file(ts_file)
        nts_df <- process_bed_file(nts_file)
        
        if (!is.null(ts_df) && !is.null(nts_df)) {
          total_reads <- ts_df$reads + nts_df$reads
          ratio <- ifelse(total_reads == 0, 0, ts_df$reads / total_reads)
          
          quartile_data <- data.frame(Quartile = quartile,
                                      TS_Ratio = ratio,
                                      Timepoint = time_point,
                                      Strain = strain,
                                      Damage = damage)
          data_combined <- rbind(data_combined, quartile_data)
        }
      }
    }
    
    if (nrow(data_combined) == 0) {
      warning(paste("No data for:", strain, damage))
      next
    }
    
    data_combined <- data_combined %>%
      filter(is.finite(TS_Ratio)) %>%
      mutate(
        Timepoint = factor(Timepoint, levels = time_points),
        Quartile = factor(Quartile, levels = quartiles, labels = quartile_labels)
      )
    
    quartile_indices <- setNames(1:4, quartile_labels)
    time_levels <- levels(data_combined$Timepoint)
    
    # median line positions
    median_lines <- data_combined %>%
      group_by(Timepoint, Quartile) %>%
      summarise(median_ratio = median(TS_Ratio, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        QuartileIndex = as.numeric(factor(Quartile, levels = quartile_labels)),
        TimepointIndex = as.numeric(factor(Timepoint, levels = time_points)),
        xpos = TimepointIndex + (QuartileIndex - 2.5) * 0.2
      )
    
    # Define y-axis
    if (strain %in% c("wt", "csb")) {
      y_limits <- c(0.25, 0.85)
      y_breaks <- c(0.25, 0.5, 0.75)
    } else if (strain == "xpc") {
      y_limits <- c(0.25, 1.17)
      y_breaks <- c(0.25, 0.5, 0.75, 1.0)
    } else {
      y_limits <- c(0.25, 1.25)
      y_breaks <- c(0.25, 0.5, 0.75, 1.0)
    }
    
    # Generate base plot
    p <- ggplot(data_combined, aes(x = Timepoint, y = TS_Ratio, fill = Quartile)) +
      geom_violin(trim = TRUE, color = NA, alpha = 0.9, position = position_dodge(0.75)) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5, position = position_dodge(0.75)) +
      geom_line(data = median_lines, aes(x = xpos, y = median_ratio, group = Timepoint),
                color = "black", size = 0.35, linetype = "dotted", inherit.aes = FALSE) +
      scale_fill_manual(values = rep("#B383BA", 4)) +
      scale_y_continuous(limits = y_limits, breaks = y_breaks) +
      labs(title = paste0("TS / (TS + NTS) - ", strain, " ", damage),
           x = "Timepoint", y = "TS / (TS + NTS)") +
      theme_classic() +
      theme(
        panel.grid.major.y = element_line(color = "grey80", size = 0.2),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    # Define comparisons
    adjacent_pairs <- list(
      c("0-25%", "25-50%"),
      c("25-50%", "50-75%"),
      c("50-75%", "75-100%")
    )
    adjacent_pair_tags <- purrr::map_chr(adjacent_pairs, ~ paste0(.x[1], "_", .x[2]))
    all_pairs <- combn(quartile_labels, 2, simplify = FALSE)
    
    # Stats
    comparison_results <- data_combined %>%
      group_by(Timepoint) %>%
      group_split() %>%
      purrr::map_df(function(subset) {
        tp <- unique(subset$Timepoint)
        purrr::map_dfr(all_pairs, function(pair) {
          sub <- subset %>% filter(Quartile %in% pair)
          sub$Quartile <- droplevels(sub$Quartile)
          if (length(unique(sub$Quartile)) == 2) {
            test <- t.test(TS_Ratio ~ Quartile, data = sub, var.equal = FALSE)
            data.frame(
              Timepoint = tp,
              group1 = pair[1],
              group2 = pair[2],
              p = test$p.value
            )
          } else {
            NULL
          }
        })
      }) %>%
      mutate(
        p.adj = p.adjust(p, method = "bonferroni"),
        p.signif = cut(p.adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("***", "**", "*", "ns")),
        TimepointIndex = as.numeric(factor(Timepoint, levels = time_points)),
        group1_index = quartile_indices[group1],
        group2_index = quartile_indices[group2],
        pair_tag = paste0(group1, "_", group2),
        is_adjacent = pair_tag %in% adjacent_pair_tags,
        bracket_length = abs(group1_index - group2_index),
        xmin = TimepointIndex + (group1_index - 2.5) * ifelse(pair_tag %in% adjacent_pair_tags, 0.18, 0.2),
        xmax = TimepointIndex + (group2_index - 2.5) * ifelse(pair_tag %in% adjacent_pair_tags, 0.18, 0.2)
      )
    
    # y-position offset for xpc
    y_offset <- if (strain == "xpc") 0.25 else 0
    
    comparison_results <- comparison_results %>%
      group_by(Timepoint) %>%
      arrange(is_adjacent, bracket_length, group1_index, group2_index) %>%
      mutate(
        y.position = ifelse(
          is_adjacent,
          0.75 + y_offset,
          0.75 + y_offset + (0.045 * (y_limits[2] - y_limits[1])) * row_number()
        )
      ) %>%
      ungroup()
    
    # Add stats to plot
    final_plot <- p + stat_pvalue_manual(
      data = comparison_results,
      label = "p.signif",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.002,
      size = 2,
      inherit.aes = FALSE
    )
    
    # Save
    output_file <- paste0(strain, "_", damage, "_TS_ratios_violin.pdf")
    ggsave(output_file, plot = final_plot, width = 8, height = 3)
    cat("✅ Saved plot:", output_file, "\n")
  }
}

cat("\n🎉 All plots complete.\n")

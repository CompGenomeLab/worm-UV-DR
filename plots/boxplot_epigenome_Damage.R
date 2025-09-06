#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  need <- c("data.table","ggplot2","optparse")
  for (p in need) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
  library(data.table); library(ggplot2); library(optparse)
})

opt <- OptionParser(
  description = "Damage-seq boxplots stratified by histone-mark tertiles from .tab files (mean0 column)."
) |>
  add_option(c("-i","--input-dir"),  type="character", default=".",
             help="Directory with Damage-seq .tab files (default: %default)") |>
  add_option(c("-s","--sim-dir"),    type="character", default="simulation",
             help="Directory with simulation .tab files (default: %default)") |>
  add_option(c("-o","--out-dir"),    type="character", default=".",
             help="Output directory (default: %default)") |>
  add_option(c("-t","--hist-tab"),   type="character", help="Histone .tab file (mean0 in col 5) [required]") |>
  add_option(c("-l","--hist-label"), type="character", default="H3K27me3",
             help="Histone label (default: %default)") |>
  add_option(c("--trim-top"),        type="double",  default=0.5,
             help="Trim top %% of histone signal before tertiles (default: %default)") |>
  add_option(c("--bins-label"),      type="character", default="2 kb",
             help="Bin label for axis (default: %default)") |>
  add_option(c("--norm-mode"),       type="character", default="ratio",
             help="Normalization: ratio(real/sim) or difference(real-sim) (default: %default)") |>
  add_option(c("--log2-plot"),       action="store_true", default=TRUE,
             help="Plot on log2 scale") |>
  parse_args()

stopifnot(!is.null(opt$`hist-tab`))

# Timepoints in order (0h, 8h, 24h, 48h), feel free to edit in your repo
TP_LEVELS <- c("0 hour","8 hours","24 hours","48 hours")
TP_COLORS <- c("0 hour"="#85BFC4","8 hours"="#C2665A","24 hours"="#94617D","48 hours"="#4B4B66")
TERT_LABELS <- c("Low (0–33%)","Medium (33–67%)","High (67–100%)")
BOX_WIDTH <- 0.7
POS_DODGE <- position_dodge2(width=1, padding=0, preserve="single")
EPS <- 1e-12
LOG_FLOOR <- 1e-12
PDF_W <- 7.9; PDF_H <- 5.3

read_tab_mean0 <- function(path){
  dt <- fread(path, header=FALSE, sep="\t",
              col.names=c("name","size","basesCovered","sum","mean0","mean"))
  dt$mean0
}
make_tertiles <- function(x, trim_top_pct){
  cut_trim <- as.numeric(quantile(x, probs=1-trim_top_pct/100, na.rm=TRUE))
  keep <- x <= cut_trim
  tb <- as.numeric(quantile(x[keep], probs=c(0,1/3,2/3,1), na.rm=TRUE))
  if (any(diff(tb) <= 0)) for (i in 2:4) if (tb[i] <= tb[i-1]) tb[i] <- tb[i-1] + .Machine$double.eps
  list(keep=keep, groups=cut(x, breaks=tb, include.lowest=TRUE, labels=TERT_LABELS))
}
normalize_pair <- function(real_vec, sim_vec, mode){
  n <- min(length(real_vec), length(sim_vec))
  r <- real_vec[seq_len(n)]; s <- sim_vec[seq_len(n)]
  if (mode=="ratio") list(values = r/(s + EPS)) else list(values = r - s)
}
transform_for_plot <- function(v) if (opt$`log2-plot`) log2(pmax(v, LOG_FLOOR)) else v

# ------------- FILE GLOBS (expects your standard DS naming) -------------
# Example names:
#   DS_wt_damage64_0h_merged_2kb.tab
#   DS_wt_damageCPD_24h_merged_2kb.tab
ds_files <- list.files(opt$`input-dir`,
  pattern="^DS_.*_merged_2kb\\.tab$", full.names=FALSE)
if (!length(ds_files)) stop("No Damage-seq DS_*_merged_2kb.tab files in --input-dir")

# Parse keys: damage type & time
parse_meta <- function(f){
  # DS_wt_damage64_0h_merged_2kb.tab  OR  DS_wt_damageCPD_48h_merged_2kb.tab
  p <- strsplit(f, "_", fixed=TRUE)[[1]]
  dmg <- if (grepl("damage64", f)) "Damage64" else "DamageCPD"
  tp  <- p[4] # like "0h" "8h" ...
  tp_label <- switch(tp, "0h"="0 hour","8h"="8 hours","24h"="24 hours","48h"="48 hours", tp)
  list(damage=dmg, time=tp_label)
}

# Build long dt for each damage class
build_long <- function(dmg_regex){
  sel <- ds_files[grepl(dmg_regex, ds_files)]
  if (!length(sel)) return(data.table())
  rbindlist(lapply(sel, function(f){
    m <- parse_meta(f)
    real_vals <- read_tab_mean0(file.path(opt$`input-dir`, f))
    sim_vals  <- read_tab_mean0(file.path(opt$`sim-dir`, sub("_merged_2kb\\.tab$","_merged_simulation_2kb.tab", f)))
    data.table(damage=m$damage,
               timepoint=factor(m$time, levels=TP_LEVELS),
               value=normalize_pair(real_vals, sim_vals, opt$`norm-mode`)$values)
  }))
}

h_vals <- read_tab_mean0(opt$`hist-tab`)
tt <- make_tertiles(h_vals, opt$`trim-top`)
keep_hist <- tt$keep; tert_all <- tt$groups

df <- rbind(
  build_long("damage64"), 
  build_long("damageCPD"),
  use.names=TRUE, fill=TRUE
)

# align to tertiles
n <- min(length(h_vals), nrow(df))
df <- df[seq_len(n)]
df$tertile <- tert_all[seq_len(n)][keep_hist[seq_len(n)]]
df <- df[!is.na(df$tertile)]
df$tertile   <- factor(df$tertile,   levels=TERT_LABELS)
df$timepoint <- factor(df$timepoint, levels=TP_LEVELS)

# Shared y from transformed 12.5–87.5% whiskers
bx <- as.data.table(df)[, {
  v <- transform_for_plot(value)
  data.table(ymin = quantile(v, 0.125, na.rm=TRUE), ymax = quantile(v, 0.875, na.rm=TRUE))
}, by=.(damage, tertile, timepoint)]
lo <- min(bx$ymin, na.rm=TRUE); hi <- max(bx$ymax, na.rm=TRUE); rg <- hi - lo
YLIM <- c(lo - 0.02*rg, hi + 0.10*rg)

build_boxstats <- function(d){
  dt <- as.data.table(d); dt[, v := transform_for_plot(value)]
  dt[, .(ymin=quantile(v,0.125,na.rm=TRUE),
         lower=quantile(v,0.25,na.rm=TRUE),
         middle=median(v,na.rm=TRUE),
         upper=quantile(v,0.75,na.rm=TRUE),
         ymax=quantile(v,0.875,na.rm=TRUE)),
     by=.(damage, tertile, timepoint)]
}

plot_one <- function(dmg){
  sub <- df[damage==dmg]
  bx  <- build_boxstats(sub)
  p <- ggplot(bx, aes(x=tertile, fill=timepoint)) +
    geom_boxplot(aes(ymin=ymin, lower=lower, middle=middle, upper=upper, ymax=ymax,
                     group=interaction(tertile,timepoint,drop=TRUE)),
                 stat="identity", position=POS_DODGE, width=BOX_WIDTH, outlier.shape=NA) +
    scale_fill_manual(values=TP_COLORS, drop=FALSE) +
    coord_cartesian(ylim=YLIM, expand=FALSE) +
    labs(
      title=sprintf("%s vs %s (normalized by simulation, %s)", opt$`hist-label`, dmg,
                    opt$`norm-mode`),
      x=sprintf("%s signal tertiles (%s; %.1f%% top-trim)", opt$`hist-label`, opt$`bins-label`, opt$`trim-top`),
      y=if (opt$`log2-plot`)
           sprintf("log2(Damage-seq %s, mean0 %s)", if (opt$`norm-mode`=="ratio") "(real/sim)" else "(real - sim)", opt$`bins-label`)
         else if (opt$`norm-mode`=="ratio") sprintf("Damage-seq (real / simulation), mean0 (%s)", opt$`bins-label`)
         else                                sprintf("Damage-seq (real - simulation), mean0 (%s)", opt$`bins-label`),
      fill="Time"
    ) +
    theme_bw(base_size=11) +
    theme(plot.margin=margin(t=18,r=16,b=12,l=12))
  out <- file.path(opt$`out-dir`, sprintf("%s_tertiles_box_%s_mean0_%s%s.pdf",
                                           opt$`hist-label`, dmg, opt$`norm-mode`,
                                           if (opt$`log2-plot`) "_log2" else ""))
  ggsave(out, p, width=PDF_W, height=PDF_H, device=cairo_pdf); message("Saved: ", out)
}

dir.create(opt$`out-dir`, showWarnings=FALSE, recursive=TRUE)
if (nrow(df[damage=="Damage64"])) plot_one("Damage64")
if (nrow(df[damage=="DamageCPD"])) plot_one("DamageCPD")

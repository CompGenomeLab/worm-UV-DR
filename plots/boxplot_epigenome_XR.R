#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  need <- c("data.table","ggplot2","optparse")
  for (p in need) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
  library(data.table); library(ggplot2); library(optparse)
})

# --------------------------- CLI OPTIONS ---------------------------
opt <- OptionParser(
  description = "XR-seq boxplots stratified by histone-mark tertiles from .tab files (mean0 column)."
) |>
  add_option(c("-i","--input-dir"),  type="character", default=".",
             help="Directory with XR .tab files (default: %default)") |>
  add_option(c("-s","--sim-dir"),    type="character", default="simulation",
             help="Directory with simulation .tab files (default: %default)") |>
  add_option(c("-o","--out-dir"),    type="character", default=".",
             help="Output directory (default: %default)") |>
  add_option(c("-t","--hist-tab"),   type="character", help="Histone .tab file (mean0 in col 5) [required]") |>
  add_option(c("-l","--hist-label"), type="character", default="H3K4me1",
             help="Histone label for axis/title (default: %default)") |>
  add_option(c("--trim-top"),        type="double",  default=0.5,
             help="Trim top %% of histone signal before tertiles (default: %default)") |>
  add_option(c("--bins-label"),      type="character", default="2 kb",
             help="Bin label for axis (default: %default)") |>
  add_option(c("--norm-mode"),       type="character", default="ratio",
             help="Normalization: ratio(real/sim) or difference(real-sim) (default: %default)") |>
  add_option(c("--drop-sim-zero"),   action="store_true", default=FALSE,
             help="Drop bins with sim <= 1e-12 when ratio") |>
  add_option(c("--log2-plot"),       action="store_true", default=TRUE,
             help="Plot on log2 scale (stats always on raw normalized)") |>
  add_option(c("--grouped-y"),       type="character", default="wtcsb_xpc",
             help="y-scaling: 'global' or 'wtcsb_xpc' (share y for wt+csb; xpc separate)") |>
  add_option(c("--whisk-low"),       type="double", default=0.125,
             help="Lower whisker quantile on transformed values (default: %default)") |>
  add_option(c("--whisk-high"),      type="double", default=0.875,
             help="Upper whisker quantile (default: %default)") |>
  parse_args()

stopifnot(!is.null(opt$`hist-tab`))

# --------------------------- USER CONSTANTS ---------------------------
TIME_PALETTE <- c("5 min"="#A3B57A","1 hour"="#D4B46D","8 hours"="#BF6B57",
                  "24 hours"="#8A5D71","48 hours"="#4F4D63")
TIME_LEVELS  <- names(TIME_PALETTE)
TERT_LABELS  <- c("Low (0–33%)","Medium (33–67%)","High (67–100%)")
BOX_DODGE    <- 0.85
PDF_W <- 7.9; PDF_H <- 5.3
EPS <- 1e-12
LOG_FLOOR <- 1e-12

# --------------------------- HELPERS ---------------------------
time_label <- function(t) switch(t,
  "5min"="5 min","1h"="1 hour","8h"="8 hours","24h"="24 hours","48h"="48 hours", t)

read_tab_mean0 <- function(path){
  dt <- fread(path, header=FALSE, sep="\t",
              col.names=c("name","size","basesCovered","sum","mean0","mean"))
  dt$mean0
}

make_tertiles <- function(x, trim_top_pct){
  cut_trim <- as.numeric(quantile(x, probs = 1 - trim_top_pct/100, na.rm=TRUE))
  keep <- x <= cut_trim
  tb <- as.numeric(quantile(x[keep], probs=c(0,1/3,2/3,1), na.rm=TRUE))
  if (any(diff(tb) <= 0)) for (i in 2:4) if (tb[i] <= tb[i-1]) tb[i] <- tb[i-1] + .Machine$double.eps
  list(keep=keep, groups=cut(x, breaks=tb, include.lowest=TRUE, labels=TERT_LABELS))
}

normalize_pair <- function(real_vec, sim_vec, mode, drop_sim_zero){
  n <- min(length(real_vec), length(sim_vec))
  r <- real_vec[seq_len(n)]; s <- sim_vec[seq_len(n)]
  if (mode == "ratio"){
    if (drop_sim_zero){ keep <- s > EPS; list(values = r[keep]/s[keep]) }
    else               { list(values = r/(s+EPS)) }
  } else list(values = r - s)
}

transform_for_plot <- function(v){
  if (!opt$`log2-plot`) return(v)
  log2(pmax(v, LOG_FLOOR))
}

p_to_bucket <- function(p) if (is.na(p)) "ns" else if (p < 0.01) "p<0.01" else if (p < 0.05) "p<0.05" else "ns"
p_to_stars  <- function(p) if (is.na(p)) ""   else if (p < 0.01) "**"      else if (p < 0.05) "*"       else ""

# --------------------------- LOAD HISTONE / TERTILES ---------------------------
h_vals <- read_tab_mean0(opt$`hist-tab`)
tt <- make_tertiles(h_vals, opt$`trim-top`)
keep_hist <- tt$keep; tert_all <- tt$groups

# --------------------------- XR FILES (by your standard naming) ---------------------------
xr_files <- list.files(opt$`input-dir`,
  pattern="^(wt|csb|xpc)_(5min|1h|8h|24h|48h)_(64|CPD)_merged_2kb\\.tab$", full.names=FALSE)
if (!length(xr_files)) stop("No XR *_merged_2kb.tab files found in --input-dir")

meta_of <- function(f){ p <- strsplit(f,"_",fixed=TRUE)[[1]]; list(strain=p[1], time=p[2], damage=p[3]) }

xr_dt <- rbindlist(lapply(xr_files, function(f){
  m <- meta_of(f)
  data.table(file=f,
             strain=m$strain, damage=m$damage,
             time = factor(time_label(m$time), levels=TIME_LEVELS),
             value = read_tab_mean0(file.path(opt$`input-dir`, f)))
}))

# map to simulation
sim_map <- setNames(
  file.path(opt$`sim-dir`, sub("_merged_2kb\\.tab$","_merged_simulation_2kb.tab", xr_files)),
  xr_files
)
if (!all(file.exists(sim_map))) stop("Missing simulation files in --sim-dir")

# --------------------------- BUILD LONG TABLE ---------------------------
long_all <- xr_dt[, {
  sim_vals <- read_tab_mean0(sim_map[file])
  n <- min(length(value), length(sim_vals), length(h_vals))
  norm <- normalize_pair(value[seq_len(n)], sim_vals[seq_len(n)], opt$`norm-mode`, opt$`drop-sim-zero`)
  keep <- keep_hist[seq_len(n)]
  data.table(tertile = tert_all[seq_len(n)][keep],
             timepoint = time[1],
             value = norm$values[keep],
             strain = strain[1], damage = damage[1])
}, by=.(file)]

long_all[, `:=`(
  tertile = factor(tertile, levels=TERT_LABELS),
  timepoint = factor(timepoint, levels=TIME_LEVELS)
)]

# --------------------------- Y-LIMITS (grouped or global) ---------------------------
tmp <- copy(long_all); tmp[, value_plot := transform_for_plot(value)]
bx <- tmp[, .(ymin = quantile(value_plot, opt$`whisk-low`, na.rm=TRUE),
              ymax = quantile(value_plot, opt$`whisk-high`, na.rm=TRUE)),
          by=.(strain, damage, tertile, timepoint)]

pad_top <- 0.22; pad_bot <- 0.06
y_for <- function(subset_bx){
  lo <- min(subset_bx$ymin, na.rm=TRUE); hi <- max(subset_bx$ymax, na.rm=TRUE); rg <- hi - lo
  c(lo - pad_bot*rg, hi + pad_top*rg)
}

YLIM_WTCSB <- y_for(bx[strain %in% c("wt","csb")])
YLIM_XPC   <- y_for(bx[strain == "xpc"])
YLIM_GLOBAL <- y_for(bx)

pick_ylim <- function(strain){
  if (opt$`grouped-y` == "wtcsb_xpc") if (strain %in% c("wt","csb")) YLIM_WTCSB else YLIM_XPC
  else YLIM_GLOBAL
}

# --------------------------- STATS ---------------------------
WITHIN <- long_all[, {
  out <- list(); base <- "5 min"
  later <- setdiff(levels(timepoint), base)
  for (tp in later){
    x <- value[timepoint == base]; y <- value[timepoint == tp]
    n <- min(length(x), length(y))
    if (n > 0){
      p <- tryCatch(wilcox.test(x[seq_len(n)], y[seq_len(n)],
                                paired=TRUE, exact=FALSE)$p.value, error=function(e) NA_real_)
      out[[tp]] <- data.table(timepoint=tp, compare=paste0(base," vs ",tp), p=p)
    }
  }
  rbindlist(out)
}, by=.(strain, damage, tertile)]
WITHIN[, `:=`(p_adj = p.adjust(p, "BH"), stars=vapply(p.adjust(p,"BH"), p_to_stars, ""), sig=vapply(p.adjust(p,"BH"), p_to_bucket, ""))]

# --------------------------- PLOTTER ---------------------------
build_boxstats <- function(df){
  dt <- as.data.table(df); dt[, value_plot := transform_for_plot(value)]
  dt[, .(ymin=quantile(value_plot, opt$`whisk-low`, na.rm=TRUE),
         lower=quantile(value_plot, 0.25, na.rm=TRUE),
         middle=median(value_plot, na.rm=TRUE),
         upper=quantile(value_plot, 0.75, na.rm=TRUE),
         ymax=quantile(value_plot, opt$`whisk-high`, na.rm=TRUE)),
     by=.(tertile, timepoint)]
}

make_one <- function(df, title, outfile, ylim, ann_df){
  bx <- build_boxstats(df)
  p  <- ggplot(bx, aes(x=tertile, fill=timepoint)) +
    geom_boxplot(aes(ymin=ymin, lower=lower, middle=middle, upper=upper, ymax=ymax,
                     group=interaction(tertile,timepoint, drop=TRUE)),
                 stat="identity",
                 position=position_dodge2(width=BOX_DODGE, padding=0),
                 width=0.8, outlier.shape=NA) +
    scale_fill_manual(values=TIME_PALETTE, drop=FALSE) +
    scale_y_continuous(limits=ylim, expand=c(0,0)) +
    labs(
      title=title,
      x=sprintf("%s signal tertiles (%s; %.1f%% top-trim)", opt$`hist-label`, opt$`bins-label`, opt$`trim-top`),
      y=if (opt$`log2-plot`)
           sprintf("log2(XR-seq %s, mean0 %s)", if (opt$`norm-mode`=="ratio") "(real/sim)" else "(real - sim)", opt$`bins-label`)
         else if (opt$`norm-mode`=="ratio") sprintf("XR-seq (real / simulation), mean0 (%s)", opt$`bins-label`)
         else                                sprintf("XR-seq (real - simulation), mean0 (%s)", opt$`bins-label`),
      fill="Time"
    ) +
    theme_bw(base_size=11) +
    theme(plot.margin=margin(t=18,r=16,b=12,l=12))

  if (nrow(ann_df)){
    key <- bx[, .(tertile, timepoint, ymax)]
    ann <- merge(ann_df[stars != ""], key, by=c("tertile","timepoint"), all.x=TRUE)
    rg <- diff(ylim)
    ann[, y := pmin(ymax + 0.02*rg, ylim[2] - 0.01*rg)]
    p <- p + geom_text(data=ann, aes(x=tertile, y=y, label=stars, group=timepoint),
                       position=position_dodge2(width=BOX_DODGE, padding=0),
                       vjust=0, size=4.2)
  }

  ggsave(outfile, p, width=PDF_W, height=PDF_H, device=cairo_pdf); message("Saved: ", outfile)
}

# split by (strain, damage)
dir.create(opt$`out-dir`, showWarnings=FALSE, recursive=TRUE)
spl <- split(long_all, list(long_all$strain, long_all$damage), drop=TRUE)
for (nm in names(spl)){
  df <- spl[[nm]]; s <- unique(df$strain); d <- unique(df$damage)
  ylim <- pick_ylim(s)
  ann  <- WITHIN[strain==s & damage==d, .(tertile, timepoint, stars)]
  out  <- file.path(opt$`out-dir`, sprintf("XR_%s_%s_by_%s_tertiles_2kb_simNorm%s.pdf",
                                           s, d, opt$`hist-label`, if (opt$`log2-plot`) "_log2" else ""))
  ttl  <- sprintf("XR-seq %s normalized (%s) — %s", toupper(d), opt$`norm-mode`, toupper(s))
  make_one(df, ttl, out, ylim, ann)
}

# write stats
fwrite(WITHIN[, .(strain,damage,tertile,timepoint,stars,sig)], 
       file.path(opt$`out-dir`, sprintf("XR_within_tertile_wilcoxon_%s_simNorm.csv", opt$`hist-label`)))

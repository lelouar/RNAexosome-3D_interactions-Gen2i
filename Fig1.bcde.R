# STAR METHODS - Figure 1: ChIP-seq Recruitment Profiles
# ============================================================
# Produces separate figures:
#   fig1b_venn.png
#   fig1c_profile_zcchc8_peaks.png   -- profiles at ZCCHC8 recruitment sites
#   fig1d_profile_zfc3h1_peaks.png   -- profiles at ZFC3H1 recruitment sites
#   fig1e_profile_mtr4_peaks.png     -- profiles at MTR4 recruitment sites
#   fig1c_heatmap_zcchc8sorted.png   -- ZCCHC8-sorted ±10kb: [ZCCHC8 | MTR4]
#   fig1d_heatmap_zfc3h1sorted.png   -- ZFC3H1-sorted ±10kb: [ZFC3H1 | MTR4]
#   fig1e_heatmap_mtr4sorted.png     -- MTR4-sorted  ±10kb: [MTR4 | ZFC3H1 | ZCCHC8]
#
# Color conventions (project-wide):
#   MTR4   = "#ff5733"   ZFC3H1 = "#1898dd"   ZCCHC8 = "#7bb84f"
#   Control sites = "grey60"
#
# Working directory: project root (00Kiernan_PROJECT-001/)
# Run: Rscript STAR_METHODS/Fig1_Recruitment_Profiles.R
# ============================================================

source("STAR_METHODS/Analysis_Functions.R")
source("STAR_METHODS/g2i_plotfig_lib.R")   # provides redim_matrix()
library(pheatmap)
library(eulerr)
library(ggpubr)
library(grid)
library(gridExtra)

# ============================================================
# 1. File Paths
# ============================================================

bed_mtr4 <- check_file("STAR_METHODS/data/mtr4_peaks.bed")
bed_z1   <- check_file("STAR_METHODS/data/z1_peaks.bed")
bed_z8   <- check_file("STAR_METHODS/data/z8_peaks.bed")

bw_mtr4  <- check_file("STAR_METHODS/data/mtr4.bw")
bw_z1    <- check_file("STAR_METHODS/data/z1.bw")
bw_z8    <- check_file("STAR_METHODS/data/z8.bw")

path_sgl <- check_file("STAR_METHODS/data/sglTSS.gr")

# ============================================================
# 2. Load and Format Peak Files (Ensembl chr style)
# ============================================================

message("Loading peak BED files...")
peaks_mtr4 <- addSeqinfo(import(bed_mtr4))
peaks_z1   <- addSeqinfo(import(bed_z1))
peaks_z8   <- addSeqinfo(import(bed_z8))

message(sprintf("Peaks: MTR4=%d  ZFC3H1=%d  ZCCHC8=%d",
                length(peaks_mtr4), length(peaks_z1), length(peaks_z8)))

# ============================================================
# 3. Figure 1b - Venn Diagram (3-way peak overlap)
# ============================================================

message("Computing Venn diagram...")

# Peak-based counts for diagram (standard for ChIP-seq overlap Venns)
n_mtr4    <- length(peaks_mtr4)
n_z1      <- length(peaks_z1)
n_z8      <- length(peaks_z8)
n_mtr4_z1 <- length(subsetByOverlaps(peaks_mtr4, peaks_z1))
n_mtr4_z8 <- length(subsetByOverlaps(peaks_mtr4, peaks_z8))
n_z1_z8   <- length(subsetByOverlaps(peaks_z1,   peaks_z8))
n_all     <- length(subsetByOverlaps(
               subsetByOverlaps(peaks_mtr4, peaks_z1), peaks_z8))

venn_counts <- c(
  "MTR4"               = n_mtr4 - n_mtr4_z1 - n_mtr4_z8 + n_all,
  "ZFC3H1"             = n_z1   - n_mtr4_z1 - n_z1_z8   + n_all,
  "ZCCHC8"             = n_z8   - n_mtr4_z8 - n_z1_z8   + n_all,
  "MTR4&ZFC3H1"        = n_mtr4_z1 - n_all,
  "MTR4&ZCCHC8"        = n_mtr4_z8 - n_all,
  "ZFC3H1&ZCCHC8"      = n_z1_z8   - n_all,
  "MTR4&ZFC3H1&ZCCHC8" = n_all
)
venn_counts <- pmax(venn_counts, 0)

cat("\nVenn counts:\n"); print(venn_counts)

# Fisher test: enrichment of triple co-binding at sglTSS gene promoters
# Universe = 14741 single-TSS gene TSSs
sglTSS   <- readRDS(path_sgl)
tss_pts  <- resize(sglTSS, 1L, "start")
strand(tss_pts) <- "*"

in_m <- countOverlaps(tss_pts, peaks_mtr4, maxgap = 500L) > 0
in_z1_u <- countOverlaps(tss_pts, peaks_z1,   maxgap = 500L) > 0
in_z8_u <- countOverlaps(tss_pts, peaks_z8,   maxgap = 500L) > 0

in_all3   <- in_m & in_z1_u & in_z8_u
in_any    <- in_m | in_z1_u | in_z8_u
ft_venn1b <- fisher.test(matrix(c(
  sum( in_all3),           sum(!in_all3 &  in_any),
  sum(!in_any),            0L   # no site is in none AND in all3
), nrow = 2), simulate.p.value = TRUE, B = 1e5)
# Use standard 2x2: triple-co-bound vs rest
ft_venn1b <- fisher.test(matrix(c(
  sum( in_all3 &  in_any),   sum(!in_all3 &  in_any),
  sum( in_all3 & !in_any),   sum(!in_all3 & !in_any)
), nrow = 2))

p_fisher <- ft_venn1b$p.value
message(sprintf("Venn Fisher p = %.2e  OR = %.2f", p_fisher, ft_venn1b$estimate))

# Format p-value label
fmt_p_exp <- function(p) {
  e <- floor(log10(p))
  sprintf("p < 1\u00d710^%d", e)
}

venn_fit <- euler(venn_counts, shape = "ellipse")
p_venn <- plot(venn_fit,
  fills      = list(fill = c("#ff5733", "#1898dd", "#7bb84f"), alpha = 0.5),
  labels     = list(labels = c("MTR4", "ZFC3H1", "ZCCHC8"), fontsize = 12),
  quantities = list(fontsize = 10, type = "counts"),
  edges      = list(col = "white", lwd = 1),
  main       = list(label = sprintf("Peak overlap\n%s", fmt_p_exp(p_fisher)),
                    fontsize = 10)
)

png(file.path(plots_dir, "fig1b_venn.png"), width = 600, height = 600, res = 120)
print(p_venn)
dev.off()
message("Saved: fig1b_venn.png")

# ============================================================
# 4. Generate Control Sites
# ============================================================
# Random genomic sites (UCSC style) not overlapping any peak set ±2kb.
# Used as the grey "control sites" in average profiles.

message("Generating control sites...")

generate_control_sites <- function(n, peaks_list, bw_path, seed = 42) {
  set.seed(seed)
  si   <- seqinfo(BigWigFile(bw_path))
  lens <- seqlengths(si)
  lens <- lens[!is.na(lens)]
  # Keep standard autosomes + sex chromosomes (UCSC or Ensembl style)
  is_ucsc <- any(grepl("^chr", names(lens)))
  if (is_ucsc) {
    lens <- lens[grepl("^chr([0-9]+|[XY])$", names(lens))]
  } else {
    lens <- lens[grepl("^([0-9]+|[XY])$", names(lens))]
  }
  message(sprintf("  BigWig chromosomes available: %d (%s style)",
                  length(lens), if (is_ucsc) "UCSC" else "Ensembl"))

  # Merge all peaks; match seqlevel style to BigWig for overlap filtering
  all_peaks_gr <- do.call(c, unname(peaks_list))
  if (is_ucsc && !any(grepl("^chr", seqlevels(all_peaks_gr)))) {
    seqlevels(all_peaks_gr) <- paste0("chr", seqlevels(all_peaks_gr))
  } else if (!is_ucsc && any(grepl("^chr", seqlevels(all_peaks_gr)))) {
    seqlevels(all_peaks_gr) <- sub("^chr", "", seqlevels(all_peaks_gr))
  }
  seqlevels(all_peaks_gr, pruning.mode = "coarse") <- names(lens)
  all_excl <- reduce(all_peaks_gr)

  n_cand <- n * 8L
  probs  <- as.numeric(lens) / sum(as.numeric(lens))
  chrs   <- sample(names(lens), n_cand, replace = TRUE, prob = probs)
  pos    <- mapply(function(ch) {
    max_start <- lens[ch] - 4002L
    if (is.na(max_start) || max_start < 1L) return(NA_integer_)
    sample.int(max_start, 1L) + 2001L
  }, chrs)
  ok <- !is.na(pos)
  sites <- GRanges(seqnames = chrs[ok],
                   ranges   = IRanges(pos[ok], pos[ok]))
  seqlevels(sites) <- names(lens)
  seqlengths(sites) <- lens

  sites <- sites[!overlapsAny(sites, all_excl, maxgap = 2000L)]
  head(sites, n)
}

N_SITES <- 2000L
ctrl_sites <- generate_control_sites(
  N_SITES,
  list(peaks_mtr4, peaks_z1, peaks_z8),
  bw_mtr4
)
message(sprintf("Control sites generated: %d", length(ctrl_sites)))

# ============================================================
# 5. Profile Plot Helper
# ============================================================
# Profile = absolute colMeans (no background subtraction).
# Shows all 3 proteins (colored) + all 3 at control sites (grey).
# Wilcoxon test for each protein at center bins vs control.
# Marks *** if p < 1e-6 for the primary protein (the "sorted" one).

EXTEND   <- 2000L
BIN_SIZE <- 50L    # 50bp bins → 80 columns for ±2kb (higher resolution profiles)
COL <- c(MTR4 = "#ff5733", ZFC3H1 = "#1898dd", ZCCHC8 = "#7bb84f")

make_chip_profile_plot <- function(spec_mats,   # named list: MTR4, ZFC3H1, ZCCHC8 at specific sites
                                   ctrl_mats,   # same proteins at control sites
                                   primary,     # protein name for Wilcoxon annotation
                                   extend, title = "") {
  n_bins      <- ncol(spec_mats[[1]])
  positions   <- seq(-extend, extend, length.out = n_bins)
  center_bins <- seq(round(n_bins * 0.45), round(n_bins * 0.55))

  # Wilcoxon: primary protein center-bin signal, specific vs control
  spec_c <- rowMeans(spec_mats[[primary]][, center_bins, drop = FALSE], na.rm = TRUE)
  ctrl_c <- rowMeans(ctrl_mats[[primary]][, center_bins, drop = FALSE], na.rm = TRUE)
  pval   <- wilcox.test(spec_c, ctrl_c, alternative = "greater")$p.value
  annot  <- if (pval < 1e-6) "***" else if (pval < 0.05) "*" else "NS"
  message(sprintf("  Wilcoxon %s: p=%.2e  %s", primary, pval, annot))

  proteins <- names(spec_mats)

  # Build long data frame
  rows <- vector("list", length(proteins) * 2)
  i <- 1L
  for (nm in proteins) {
    rows[[i]]   <- data.frame(Position = positions,
                              Signal   = colMeans(spec_mats[[nm]], na.rm = TRUE),
                              Protein  = nm, Type = "Specific", stringsAsFactors = FALSE)
    rows[[i+1]] <- data.frame(Position = positions,
                              Signal   = colMeans(ctrl_mats[[nm]], na.rm = TRUE),
                              Protein  = nm, Type = "Control",  stringsAsFactors = FALSE)
    i <- i + 2L
  }
  df <- do.call(rbind, rows)
  df$group     <- paste0(df$Protein, "_", df$Type)
  df$LegendKey <- ifelse(df$Type == "Specific", df$Protein, "Control sites")
  df$col       <- ifelse(df$Type == "Specific", COL[df$Protein], "grey60")

  col_map    <- c(COL[proteins], "Control sites" = "grey60")
  ltype_map  <- setNames(
    c(rep("solid", length(proteins)), "dashed"),
    c(proteins, "Control sites")
  )

  y_max <- max(unlist(lapply(spec_mats, function(m) colMeans(m, na.rm = TRUE))), na.rm = TRUE)

  ggplot(df, aes(x = Position, y = Signal,
                 color = LegendKey, linetype = LegendKey, group = group)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = col_map,
                       breaks = c(proteins, "Control sites")) +
    scale_linetype_manual(values = ltype_map,
                          breaks = c(proteins, "Control sites")) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
    annotate("text", x = 0, y = y_max * 1.1, label = annot,
             size = 5, fontface = "bold", hjust = 0.5) +
    theme_classic(base_size = 11) +
    labs(x  = "Distance from peak center (bp)",
         y  = "Normalized ChIP-seq reads",
         title = title, color = NULL, linetype = NULL) +
    theme(legend.position = "right",
          legend.key.width = unit(1.2, "cm"))
}

# ============================================================
# 6. Extract signal matrices for all 6 site sets (±2kb)
# ============================================================

message("Extracting signal matrices at ±2kb...")

set.seed(42)
idx_z8 <- if (length(peaks_z8) > N_SITES) sample(length(peaks_z8), N_SITES) else seq_along(peaks_z8)
idx_z1 <- if (length(peaks_z1) > N_SITES) sample(length(peaks_z1), N_SITES) else seq_along(peaks_z1)
idx_m  <- if (length(peaks_mtr4) > N_SITES) sample(length(peaks_mtr4), N_SITES) else seq_along(peaks_mtr4)

# --- at ZCCHC8 peaks (Fig 1c) ---
message("  ZCCHC8 peaks (Fig 1c)...")
mat_mtr4_at_z8 <- extract_signal_matrix(bw_mtr4, peaks_z8[idx_z8], BIN_SIZE, EXTEND)
mat_z1_at_z8   <- extract_signal_matrix(bw_z1,   peaks_z8[idx_z8], BIN_SIZE, EXTEND)
mat_z8_at_z8   <- extract_signal_matrix(bw_z8,   peaks_z8[idx_z8], BIN_SIZE, EXTEND)

# --- at ZFC3H1 peaks (Fig 1d) ---
message("  ZFC3H1 peaks (Fig 1d)...")
mat_mtr4_at_z1 <- extract_signal_matrix(bw_mtr4, peaks_z1[idx_z1], BIN_SIZE, EXTEND)
mat_z1_at_z1   <- extract_signal_matrix(bw_z1,   peaks_z1[idx_z1], BIN_SIZE, EXTEND)
mat_z8_at_z1   <- extract_signal_matrix(bw_z8,   peaks_z1[idx_z1], BIN_SIZE, EXTEND)

# --- at MTR4 peaks (Fig 1e) ---
message("  MTR4 peaks (Fig 1e)...")
mat_mtr4_at_m  <- extract_signal_matrix(bw_mtr4, peaks_mtr4[idx_m], BIN_SIZE, EXTEND)
mat_z1_at_m    <- extract_signal_matrix(bw_z1,   peaks_mtr4[idx_m], BIN_SIZE, EXTEND)
mat_z8_at_m    <- extract_signal_matrix(bw_z8,   peaks_mtr4[idx_m], BIN_SIZE, EXTEND)

# --- at control sites ---
message("  Control sites...")
mat_mtr4_ctrl <- extract_signal_matrix(bw_mtr4, ctrl_sites, BIN_SIZE, EXTEND)
mat_z1_ctrl   <- extract_signal_matrix(bw_z1,   ctrl_sites, BIN_SIZE, EXTEND)
mat_z8_ctrl   <- extract_signal_matrix(bw_z8,   ctrl_sites, BIN_SIZE, EXTEND)

ctrl_mats <- list(MTR4 = mat_mtr4_ctrl, ZFC3H1 = mat_z1_ctrl, ZCCHC8 = mat_z8_ctrl)

# ============================================================
# 7. Fig 1c - Average profiles at ZCCHC8 peaks
# ============================================================

message("Plotting Fig 1c (ZCCHC8 peaks)...")
p_1c <- make_chip_profile_plot(
  spec_mats = list(MTR4 = mat_mtr4_at_z8, ZFC3H1 = mat_z1_at_z8, ZCCHC8 = mat_z8_at_z8),
  ctrl_mats = ctrl_mats,
  primary   = "ZCCHC8",
  extend    = EXTEND,
  title     = "c   ChIP-seq at ZCCHC8 recruitment sites"
)
save_plot(p_1c, "fig1c_profile_zcchc8_peaks.png", width = 6, height = 4)

# ============================================================
# 8. Fig 1d - Average profiles at ZFC3H1 peaks
# ============================================================

message("Plotting Fig 1d (ZFC3H1 peaks)...")
p_1d <- make_chip_profile_plot(
  spec_mats = list(MTR4 = mat_mtr4_at_z1, ZFC3H1 = mat_z1_at_z1, ZCCHC8 = mat_z8_at_z1),
  ctrl_mats = ctrl_mats,
  primary   = "ZFC3H1",
  extend    = EXTEND,
  title     = "d   ChIP-seq at ZFC3H1 recruitment sites"
)
save_plot(p_1d, "fig1d_profile_zfc3h1_peaks.png", width = 6, height = 4)

# ============================================================
# 9. Fig 1e - Average profiles at MTR4 peaks
# ============================================================

message("Plotting Fig 1e (MTR4 peaks)...")
p_1e <- make_chip_profile_plot(
  spec_mats = list(MTR4 = mat_mtr4_at_m, ZFC3H1 = mat_z1_at_m, ZCCHC8 = mat_z8_at_m),
  ctrl_mats = ctrl_mats,
  primary   = "MTR4",
  extend    = EXTEND,
  title     = "e   ChIP-seq at MTR4 recruitment sites"
)
save_plot(p_1e, "fig1e_profile_mtr4_peaks.png", width = 6, height = 4)

# ============================================================
# 10. Heatmaps at ±10kb (separate composite PNG per panel)
# ============================================================

message("Extracting matrices at ±10kb...")

EXTEND_HEAT   <- 10000L
BIN_SIZE_HEAT <- 100L  # 100bp bins → 200 columns input for redim_matrix

# All peaks, no subsampling (max_sites default = 3000 in extract_signal_matrix)
mat_z8_10k_z8    <- extract_signal_matrix(bw_z8,   peaks_z8,   BIN_SIZE_HEAT, EXTEND_HEAT)
mat_mtr4_10k_z8  <- extract_signal_matrix(bw_mtr4, peaks_z8,   BIN_SIZE_HEAT, EXTEND_HEAT)

mat_z1_10k_z1    <- extract_signal_matrix(bw_z1,   peaks_z1,   BIN_SIZE_HEAT, EXTEND_HEAT)
mat_mtr4_10k_z1  <- extract_signal_matrix(bw_mtr4, peaks_z1,   BIN_SIZE_HEAT, EXTEND_HEAT)

mat_mtr4_10k_m   <- extract_signal_matrix(bw_mtr4, peaks_mtr4, BIN_SIZE_HEAT, EXTEND_HEAT)
mat_z1_10k_m     <- extract_signal_matrix(bw_z1,   peaks_mtr4, BIN_SIZE_HEAT, EXTEND_HEAT)
mat_z8_10k_m     <- extract_signal_matrix(bw_z8,   peaks_mtr4, BIN_SIZE_HEAT, EXTEND_HEAT)

# Linear white→black palette.
# Contrast is driven by log-spaced breaks in make_hm_grob, not by a skewed palette.
col_grey <- colorRampPalette(c("white", "black"))(100)

winsorize <- function(m, p = 0.99) {
  thresh <- quantile(m, p, na.rm = TRUE)
  m[m > thresh] <- thresh
  m
}

# Display dimensions coherent with PNG pixel budgets:
#   2-col PNG (800×1200 px): each heatmap ~350px wide × 1050px tall
#   3-col PNG (1200×1200 px): each heatmap ~370px wide × 1050px tall
# → target_height = 1000 rows  (1px ≈ 1 site row at 120 dpi)
#   target_width  = 100 cols   (≈ 3.5px per 200bp bin; input has 200 cols)
TARGET_H <- 1000L
TARGET_W <- 100L
N_CORE   <- 14L

rdim <- function(mat_sorted) {
  redim_matrix(winsorize(mat_sorted),
               target_height  = min(TARGET_H, nrow(mat_sorted) - 1L),
               target_width   = TARGET_W,
               summary_func   = function(x) mean(x, na.rm = TRUE),
               n_core         = N_CORE)
}

make_hm_grob <- function(mat, title = "") {
  mat  <- pmax(mat, 0)          # clip negatives from input normalization
  vmax <- max(mat, na.rm = TRUE)
  if (vmax <= 0) vmax <- 1e-6
  # Log-spaced breaks: same 100 colors but mapped logarithmically over [0, vmax].
  # Low values (background) stay white; mid-high signal transitions rapidly to dark.
  breaks <- c(-1e-9, expm1(seq(0, log1p(vmax), length.out = 100)))
  pheatmap::pheatmap(mat,
    cluster_rows = FALSE, cluster_cols = FALSE,
    show_rownames = FALSE, show_colnames = FALSE,
    color        = col_grey,
    breaks       = breaks,
    border_color = NA, main = title,
    silent       = TRUE
  )$gtable
}

# --- Fig 1c heatmap: ZCCHC8-sorted | [ZCCHC8 | MTR4] ---
message("Heatmap Fig 1c (ZCCHC8-sorted)...")
ord_z8 <- order(rowMeans(mat_z8_10k_z8, na.rm = TRUE), decreasing = TRUE)
g_z8   <- make_hm_grob(rdim(mat_z8_10k_z8[ord_z8, ]),   "ZCCHC8 at ZCCHC8 peaks")
g_m_z8 <- make_hm_grob(rdim(mat_mtr4_10k_z8[ord_z8, ]), "MTR4 at ZCCHC8 peaks")

png(file.path(plots_dir, "fig1c_heatmap_zcchc8sorted.png"),
    width = 800, height = 1200, res = 120)
grid.arrange(g_z8, g_m_z8, ncol = 2,
  top = textGrob("c  ZCCHC8 recruitment sites (\u00b110kb)",
                 gp = gpar(fontsize = 11, fontface = "bold")))
dev.off()
message("Saved: fig1c_heatmap_zcchc8sorted.png")

# --- Fig 1d heatmap: ZFC3H1-sorted | [ZFC3H1 | MTR4] ---
message("Heatmap Fig 1d (ZFC3H1-sorted)...")
ord_z1 <- order(rowMeans(mat_z1_10k_z1, na.rm = TRUE), decreasing = TRUE)
g_z1   <- make_hm_grob(rdim(mat_z1_10k_z1[ord_z1, ]),   "ZFC3H1 at ZFC3H1 peaks")
g_m_z1 <- make_hm_grob(rdim(mat_mtr4_10k_z1[ord_z1, ]), "MTR4 at ZFC3H1 peaks")

png(file.path(plots_dir, "fig1d_heatmap_zfc3h1sorted.png"),
    width = 800, height = 1200, res = 120)
grid.arrange(g_z1, g_m_z1, ncol = 2,
  top = textGrob("d  ZFC3H1 recruitment sites (\u00b110kb)",
                 gp = gpar(fontsize = 11, fontface = "bold")))
dev.off()
message("Saved: fig1d_heatmap_zfc3h1sorted.png")

# --- Fig 1e heatmap: MTR4-sorted | [MTR4 | ZFC3H1 | ZCCHC8] ---
message("Heatmap Fig 1e (MTR4-sorted)...")
ord_m  <- order(rowMeans(mat_mtr4_10k_m, na.rm = TRUE), decreasing = TRUE)
g_m    <- make_hm_grob(rdim(mat_mtr4_10k_m[ord_m, ]), "MTR4 at MTR4 peaks")
g_z1_m <- make_hm_grob(rdim(mat_z1_10k_m[ord_m, ]),  "ZFC3H1 at MTR4 peaks")
g_z8_m <- make_hm_grob(rdim(mat_z8_10k_m[ord_m, ]),  "ZCCHC8 at MTR4 peaks")

png(file.path(plots_dir, "fig1e_heatmap_mtr4sorted.png"),
    width = 1200, height = 1200, res = 120)
grid.arrange(g_m, g_z1_m, g_z8_m, ncol = 3,
  top = textGrob("e  MTR4 recruitment sites (±10kb)",
                 gp = gpar(fontsize = 11, fontface = "bold")))
dev.off()
message("Saved: fig1e_heatmap_mtr4sorted.png")

message("\n=== Figure 1 complete ===")
message("Outputs in: ", plots_dir)

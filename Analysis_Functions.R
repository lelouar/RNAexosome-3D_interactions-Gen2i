# STAR METHODS - Shared Helper Functions
# Core functions for signal extraction, genomic manipulation, and statistical analysis.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(RColorBrewer)
})

# --- Utility ---

`%ni%` <- Negate('%in%')

check_file <- function(path, label = path) {
  if (!file.exists(path)) stop("File not found: ", label, "\n  -> ", path)
  invisible(path)
}

# Create output directory
plots_dir <- "STAR_METHODS/plots"
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, filename, width = 7, height = 7, dpi = 150) {
  out <- file.path(plots_dir, filename)
  ggsave(out, plot = p, width = width, height = height, dpi = dpi)
  message("Saved: ", out)
  invisible(out)
}

save_png <- function(expr, filename, width = 500, height = 500) {
  out <- file.path(plots_dir, filename)
  png(out, width = width, height = height)
  eval(expr)
  dev.off()
  message("Saved: ", out)
  invisible(out)
}

# --- GRanges Manipulation ---

#' Add Sequence Information and Filter to HeLa Chromosomes
addSeqinfo <- function(my.gr, chrTYPE = "Ensembl") {
  hela.chr <- c(as.character(1:22), "X")
  seqlevelsStyle(my.gr) <- "UCSC"
  valid <- paste0("chr", hela.chr)
  valid <- valid[valid %in% seqlevels(my.gr)]
  my.gr <- keepSeqlevels(my.gr, valid, pruning.mode = "coarse")
  seqlevelsStyle(my.gr) <- chrTYPE
  return(my.gr)
}

# Fisher test across N quantile bins: returns log2 odds-ratio and significance per bin
fisher_ntiles <- function(DF = NULL, COLN = NULL, ROWN = NULL, NTILE = 4, NAME = NULL) {
  fish.df <- NULL
  for (quart in 1:NTILE) {
    myDF <- as.data.frame(DF %>%
      group_by(!!as.name(COLN), !!as.name(ROWN)) %>%
      summarize(tot = n(), .groups = "drop"))
    pq   <- sum(myDF$tot[myDF[[COLN]] == quart  & myDF[[ROWN]] == 1], na.rm = TRUE)
    pnq  <- sum(myDF$tot[myDF[[COLN]] != quart  & myDF[[ROWN]] == 1], na.rm = TRUE)
    npq  <- sum(myDF$tot[myDF[[COLN]] == quart  & myDF[[ROWN]] == 0], na.rm = TRUE)
    npnq <- sum(myDF$tot[myDF[[COLN]] != quart  & myDF[[ROWN]] == 0], na.rm = TRUE)
    ft   <- fisher.test(matrix(c(pq, pnq, npq, npnq), nrow = 2))
    lfc  <- round(log2(ft$estimate), 2)
    sign <- ifelse(ft$p.value < 0.001, "***", ifelse(ft$p.value < 0.01, "**", ifelse(ft$p.value < 0.05, "*", "NS")))
    fish.df <- rbind.data.frame(fish.df,
      data.frame(XNAME = as.character(quart), pv = formatC(ft$p.value, format = "e", digits = 2),
                 lfcb = lfc, sign = sign, YNAME = NAME, stringsAsFactors = FALSE))
  }
  fish.df$intSign <- sign(fish.df$lfcb)
  fish.df <- fish.df %>% mutate(lfc = ifelse(abs(lfcb) > 2, intSign * 2, intSign * abs(lfcb)))
  fish.df
}

#' Extend Genomic Ranges Relative to Strand
extendGR <- function(x, upstream = 0, downstream = 0) {
  if (any(strand(x) == "*")) strand(x)[strand(x) == "*"] <- "+"
  on_plus <- strand(x) == "+"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end   <- end(x)   + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

#' Reduce (shrink) Genomic Ranges Relative to Strand
reduceGR <- function(x, upstream = 0, downstream = 0) {
  if (any(strand(x) == "*")) strand(x)[strand(x) == "*"] <- "+"
  on_plus <- strand(x) == "+"
  new_start <- start(x) + ifelse(on_plus, upstream, downstream)
  new_end   <- end(x)   - ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

#' Overlap helper: returns subsets of gr1 and gr2 with reciprocal overlaps
myOverlaps <- function(n1 = "n1", n2 = "n2", gr1.GR, gr2.GR,
                       gap = 1, ign.str = TRUE) {
  if (ign.str) {
    strand(gr1.GR) <- "*"
    strand(gr2.GR) <- "*"
  }
  gr1xgr2 <- unique(subsetByOverlaps(gr1.GR, gr2.GR, maxgap = gap))
  gr2xgr1 <- unique(subsetByOverlaps(gr2.GR, gr1.GR, maxgap = gap))
  gr1nogr2 <- gr1.GR[gr1.GR %ni% gr1xgr2]
  gr2nogr1 <- gr2.GR[gr2.GR %ni% gr2xgr1]
  list(gr1xgr2 = gr1xgr2, gr2xgr1 = gr2xgr1,
       gr1nogr2 = gr1nogr2, gr2nogr1 = gr2nogr1)
}

# --- BigWig Signal Extraction ---

#' Extract a signal matrix from a BigWig file
#'
#' Approach mirrored from g2i_serializeAnno.R::extractBigWigProfile:
#'   - import.bw(..., as = "NumericList") for fully vectorized import
#'   - mapSeqlevels + renameSeqlevels for robust style matching
#'   - data.table column-block means for fast binning
#'
#' @param bw_path  Path to BigWig file (must exist)
#' @param regions  GRanges of sites (strand-aware)
#' @param bin_size Resolution in bp (default 50)
#' @param extend   Extension around center in bp (default 2000)
#' @param max_sites Max number of sites (subsampled if larger)
#' @return Matrix [n_sites x n_bins]
extract_signal_matrix <- function(bw_path, regions, bin_size = 50,
                                  extend = 2000, max_sites = 3000) {
  check_file(bw_path)

  if (length(regions) > max_sites) {
    set.seed(42)
    regions <- regions[sample(length(regions), max_sites)]
  }

  n_bins    <- floor(2L * as.integer(extend) / bin_size)
  win_width <- 2L * as.integer(extend)

  # 1. Fixed-width windows centered on each region
  #    resize(center, 2*extend, "center") gives exactly win_width bp
  centers    <- resize(regions, 1L, "center")
  windows    <- resize(centers, win_width, "center")
  minus_idx  <- which(as.character(strand(regions)) == "-")

  # 2. Match seqlevel style to BigWig (mapSeqlevels is robust vs seqlevelsStyle<-)
  bw_lvls  <- seqlevels(BigWigFile(bw_path))
  bw_style <- if (any(grepl("^chr", bw_lvls))) "UCSC" else "NCBI"
  new_lvls <- mapSeqlevels(seqlevels(windows), bw_style)
  new_lvls <- new_lvls[!is.na(new_lvls)]
  windows_bw <- if (length(new_lvls) > 0) renameSeqlevels(windows, new_lvls) else windows
  windows_bw <- trim(windows_bw)

  # 3. Keep only windows that were not trimmed at chromosome boundaries
  valid <- which(width(windows_bw) == win_width)
  mat   <- matrix(0, nrow = length(regions), ncol = n_bins)
  if (length(valid) == 0) return(mat)

  # 4. Vectorized import: import.bw returns one numeric vector per window
  #    (one element of the NumericList per element of `which`)
  sig <- tryCatch(
    suppressWarnings(import.bw(bw_path, which = windows_bw[valid], as = "NumericList")),
    error = function(e) { warning("BigWig import error: ", e$message); NULL }
  )
  if (is.null(sig) || length(sig) == 0) return(mat)

  # 5. Build raw matrix [n_valid x win_width] via do.call(rbind)
  raw_mat <- do.call(rbind, as.list(sig))

  # 6. Strand reversal for minus-strand regions (vectorized row reversal)
  minus_valid <- which(valid %in% minus_idx)
  if (length(minus_valid) > 0)
    raw_mat[minus_valid, ] <- raw_mat[minus_valid, rev(seq_len(win_width)), drop = FALSE]

  # 7. Bin columns with data.table (vectorized block means, same as extractBigWigProfile)
  dt     <- data.table(raw_mat)
  n_cols <- ncol(dt)
  binned <- as.matrix(
    dt[, lapply(seq_len(n_bins), function(b) {
      cols <- seq(round((b - 1) * n_cols / n_bins) + 1L, round(b * n_cols / n_bins))
      cols <- cols[cols >= 1L & cols <= n_cols]
      rowMeans(.SD[, cols, with = FALSE], na.rm = TRUE)
    })]
  )

  mat[valid, ] <- binned
  mat
}

# --- Profile Normalization ---

#' Background (flank) normalization for average profiles.
#'
#' Scales every track multiplicatively so that its peripheral signal (first and
#' last `n` bins, i.e. the flanks far from the peak center) matches the
#' reference track's peripheral signal.  Preserves relative peak-to-background
#' ratios — suitable for comparing co-recruitment profiles.
#'
#' @param profile_list Named list of numeric vectors (one per ChIP track).
#' @param ref_name     Name of the reference track (default "MTR4").
#' @param n            Number of flank bins used to estimate background.
#' @return Rescaled list with the same structure.
correct_to_baseline <- function(profile_list, ref_name = "MTR4", n = 10) {
  baselines <- sapply(profile_list, function(v) mean(v[seq_len(n)], na.rm = TRUE))
  ref       <- baselines[[ref_name]]
  lapply(profile_list, function(v) v * (ref / mean(v[seq_len(n)], na.rm = TRUE)))
}

#' Z-score normalization per genomic site (mirrored from g2i_plotfig_lib.R).
#'
#' Applies R's `scale()` row-wise to a signal matrix: each site's profile is
#' converted to mean=0 / SD=1.  After column-averaging, the profile shows mean
#' deviation from each site's local background — comparable across tracks with
#' very different absolute signals.
#'
#' @param mat Numeric matrix [n_sites x n_bins].
#' @return Scaled matrix with same dimensions.
scale_per_site <- function(mat) {
  t(apply(mat, 1, scale))
}

# --- Profile Plotting ---

#' Plot Average Signal Profiles
plot_avg_profiles <- function(profile_list, extend, colors = NULL, title = "") {
  df <- do.call(rbind, lapply(names(profile_list), function(nm) {
    data.frame(
      Position = seq(-extend, extend, length.out = length(profile_list[[nm]])),
      Signal   = profile_list[[nm]],
      Group    = nm,
      stringsAsFactors = FALSE
    )
  }))
  p <- ggplot(df, aes(x = Position, y = Signal, color = Group)) +
    geom_line(linewidth = 0.8) +
    theme_classic(base_size = 12) +
    labs(x = "Distance from center (bp)", y = "Normalized Signal", title = title) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")
  if (!is.null(colors)) p <- p + scale_color_manual(values = colors)
  return(p)
}

# --- APA Heatmap Plotting ---

APA_COLORS <- colorRampPalette(c("white", "orange", "firebrick1"))(100)

#' Plot APA matrix using pheatmap
plot_apa <- function(mat, title = "APA", col = APA_COLORS, breaks = NULL,
                     show_rownames = FALSE, show_colnames = FALSE) {
  library(pheatmap)
  if (is.null(breaks)) {
    lim <- quantile(mat, 0.99, na.rm = TRUE)
    breaks <- seq(0, max(lim, 0.01), length.out = 101)
  }
  pheatmap::pheatmap(mat,
    cluster_rows  = FALSE, cluster_cols = FALSE,
    main          = title,
    color         = col,
    breaks        = breaks,
    border_color  = NA,
    show_rownames = show_rownames,
    show_colnames = show_colnames
  )
}

# --- Fisher Exact Test on Quintiles ---

#' Fisher test across NTILE quantiles
#' @param DF data.table with columns COLN (quantile) and ROWN (binary feature)
fisher_ntiles <- function(DF, COLN, ROWN, NTILE = 5, NAME = NULL) {
  fish.df <- NULL
  for (quart in 1:NTILE) {
    pq    <- nrow(DF[get(COLN) == quart & get(ROWN) == 1, ])
    pnq   <- nrow(DF[get(COLN) != quart & get(ROWN) == 1, ])
    npq   <- nrow(DF[get(COLN) == quart & get(ROWN) == 0, ])
    npnq  <- nrow(DF[get(COLN) != quart & get(ROWN) == 0, ])
    f_mat <- matrix(c(pq, pnq, npq, npnq), nrow = 2)
    ft    <- fisher.test(f_mat)
    mypv  <- formatC(ft$p.value, format = "e", digits = 2)
    lfc   <- round(log2(ft$estimate), 2)
    sign  <- ifelse(ft$p.value < 0.05,
               ifelse(ft$p.value < 0.01,
                 ifelse(ft$p.value < 0.001, "***", "**"), "*"), "NS")
    fish.df <- rbind(fish.df,
      data.frame(XNAME = as.character(quart), pv = mypv,
                 lfcb = lfc, sign = sign, YNAME = NAME,
                 stringsAsFactors = FALSE))
  }
  fish.df$intSign <- sign(fish.df$lfcb)
  fish.df <- fish.df %>%
    mutate(lfc = ifelse(abs(lfcb) > 2, intSign * 2, intSign * abs(lfcb)))
  fish.df
}

#' Plot Fisher heatmap (tile plot with significance labels)
fisher_plot <- function(DF, main = "", range = c(-2, 2),
                        xl = "Quintile", yl = "Feature") {
  colfunc <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
  DF$lfc <- pmax(range[1], pmin(range[2], DF$lfc))
  ggplot(DF, aes(x = XNAME, y = YNAME, fill = lfc)) +
    geom_tile(col = "black") +
    theme_minimal(base_size = 12) +
    scale_x_discrete(position = "top") +
    labs(x = xl, y = yl, fill = "Log2 OR", title = main) +
    geom_text(aes(label = sign), size = 5, fontface = "bold") +
    scale_fill_gradientn(colours = colfunc(10),
                         limits = range,
                         breaks = range,
                         labels = as.character(range)) +
    theme(axis.text = element_text(face = "bold"))
}

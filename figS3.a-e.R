# STAR METHODS - Figure S3: Venn Diagrams DEG vs PROMPTs
# ============================================================
# Reproduces Figure S3 panels A-E:
#   Overlap between MTR4-regulated genes (siMTR4/siCon DESeq2)
#   and Upregulated PROMPTs from siZFC3H1, siMTR4, siZCCHC8 KD.
#
# PROMPTS are defined as ncRNAs (normR, bs=500, FDR<0.05)
# whose genomic midpoint falls within 2 kb upstream of a
# single-TSS gene (sglTSS catalogue).
#
# Working directory: project root (00Kiernan_PROJECT-001/)
# Run: Rscript FigS3_Venn_DEG_PROMPTS.R
# ============================================================

source("Analysis_Functions.R")
suppressPackageStartupMessages(library(eulerr))
suppressPackageStartupMessages(library(ggpubr))

# ============================================================
# 0. Paths
# ============================================================

path_tx   <- check_file("data/matrecap_genes.rds")
path_sgl  <- check_file("data/sglTSS.gr")
path_lZ1  <- check_file("data/prompts_Z1kd.rds")
path_lM4  <- check_file("data/prompts_MTR4kd.rds")
path_lZ8  <- check_file("data/prompts_Z8kd.rds")

# ============================================================
# 1. Load data
# ============================================================

message("Loading data...")

# Gene expression data with DESeq2 results (MTR4 and ZFC3H1 KD)
# Columns: MTR4_DEREG_ADJ (padj<0.05), MTR4_LFC (log2FC siCon/siMTR4),
#          Z1_DEREG_ADJ, Z1_LFC
# NOTE: LFC convention: LFC > 0 = higher in siCon = downregulated in KD
#                       LFC < 0 = lower in siCon = upregulated in KD
tx.gr <- readRDS(path_tx)

# Single-TSS gene catalogue (unique promoter regions, ±2 kb windows)
sglTSS <- readRDS(path_sgl)

# ncRNA/PROMPTS detected by normR (binsize=500, FDR<0.05) in each KD condition
# These are genomic intervals upregulated when the respective factor is knocked down
lncZ1 <- readRDS(path_lZ1)   # siZFC3H1 KD PROMPTs
lncM4 <- readRDS(path_lM4)   # siMTR4 KD PROMPTs
lncZ8 <- readRDS(path_lZ8)   # siZCCHC8 KD PROMPTs

# ============================================================
# 2. Define PROMPT sets
# ============================================================
# A gene is associated with a PROMPT if its TSS (from sglTSS catalogue)
# is within 2 kb of an upregulated ncRNA detected in that KD condition.
# The PROMPT "circle" in the Venn = sglTSS sites near a KD ncRNA.

message("Defining PROMPT sets...")

tss_pts <- resize(sglTSS, 1L, "start")   # TSS single points
strand(tss_pts) <- "*"
strand(lncZ1)   <- "*"
strand(lncM4)   <- "*"
strand(lncZ8)   <- "*"

# TSS sites that have a ncRNA within 2 kb (the PROMPT catalogue)
promptsZ1  <- unique(subsetByOverlaps(tss_pts, lncZ1, maxgap = 2000L))
promptsMTR4 <- unique(subsetByOverlaps(tss_pts, lncM4, maxgap = 2000L))
promptsZ8  <- unique(subsetByOverlaps(tss_pts, lncZ8, maxgap = 2000L))

message(sprintf("PROMPTs: siZFC3H1=%d  siMTR4=%d  siZCCHC8=%d",
                length(promptsZ1), length(promptsMTR4), length(promptsZ8)))

# ============================================================
# 3. Gene classification flags
# ============================================================
# Each gene in tx.gr is flagged:
#   isDEG  - differentially expressed (siMTR4, padj<0.05)
#   isUP   - upregulated in siMTR4/siCon (LFC < 0, higher in KD)
#   isDOWN - downregulated in siMTR4/siCon (LFC > 0, higher in WT)
#   hasZ1  - gene TSS is in the siZFC3H1 PROMPT set
#   hasMTR4 - gene TSS is in the siMTR4 PROMPT set
#   hasZ8  - gene TSS is in the siZCCHC8 PROMPT set

tx_tss <- resize(tx.gr, 1L, "start")
strand(tx_tss) <- "*"

isDEG  <- tx.gr$MTR4_DEREG_ADJ %in% TRUE
isUP   <- isDEG & (tx.gr$MTR4_LFC < 0)    # upregulated in siMTR4
isDOWN <- isDEG & (tx.gr$MTR4_LFC > 0)    # downregulated in siMTR4

hasZ1   <- countOverlaps(tx_tss, promptsZ1,   maxgap = 500L) > 0
hasMTR4 <- countOverlaps(tx_tss, promptsMTR4, maxgap = 500L) > 0
hasZ8   <- countOverlaps(tx_tss, promptsZ8,   maxgap = 500L) > 0

message(sprintf("DEGs: all=%d  UP=%d  DOWN=%d", sum(isDEG), sum(isUP), sum(isDOWN)))
message(sprintf("Genes near PROMPT: Z1=%d  MTR4=%d  Z8=%d",
                sum(hasZ1), sum(hasMTR4), sum(hasZ8)))

# ============================================================
# 4. Compute 3-way Venn counts and Fisher p-values
# ============================================================

# Returns named vector of 7 exclusive-region counts for eulerr::euler()
# and the Fisher p-value for set A enrichment vs {B or C}
venn3_counts <- function(A, B, C, labels = c("A", "B", "C")) {
  n <- c(
    sum( A & !B & !C),    # A only
    sum(!A &  B & !C),    # B only
    sum(!A & !B &  C),    # C only
    sum( A &  B & !C),    # A & B
    sum( A & !B &  C),    # A & C
    sum(!A &  B &  C),    # B & C
    sum( A &  B &  C)     # A & B & C
  )
  names(n) <- c(
    labels[1],
    labels[2],
    labels[3],
    paste0(labels[1], "&", labels[2]),
    paste0(labels[1], "&", labels[3]),
    paste0(labels[2], "&", labels[3]),
    paste0(labels[1], "&", labels[2], "&", labels[3])
  )

  # Fisher test: enrichment of A in {B | C} vs neither
  in_prompt  <- B | C
  ft <- fisher.test(matrix(c(
    sum( A &  in_prompt), sum(!A &  in_prompt),
    sum( A & !in_prompt), sum(!A & !in_prompt)
  ), nrow = 2))
  attr(n, "fisher_p") <- ft$p.value
  attr(n, "fisher_or") <- ft$estimate
  n
}

cnt_A <- venn3_counts(isDEG,  hasZ1, hasMTR4, c("DEGs",    "siZFC3H1", "siMTR4"))
cnt_B <- venn3_counts(isUP,   hasZ1, hasMTR4, c("UP_KD",   "siZFC3H1", "siMTR4"))
cnt_C <- venn3_counts(isDOWN, hasZ1, hasMTR4, c("DOWN_KD", "siZFC3H1", "siMTR4"))
cnt_D <- venn3_counts(isUP,   hasMTR4, hasZ8, c("UP_KD",   "siMTR4",   "siZCCHC8"))
cnt_E <- venn3_counts(isDOWN, hasMTR4, hasZ8, c("DOWN_KD", "siMTR4",   "siZCCHC8"))

# Print counts for verification
fmt_venn <- function(cnt, lbl) {
  p <- attr(cnt, "fisher_p")
  or <- attr(cnt, "fisher_or")
  cat(sprintf("\n--- %s ---  Fisher p=%.2e  OR=%.2f\n", lbl, p, or))
  print(cnt)
}
fmt_venn(cnt_A, "Panel A: All DEGs x siZFC3H1 x siMTR4")
fmt_venn(cnt_B, "Panel B: UP in siMTR4 x siZFC3H1 x siMTR4")
fmt_venn(cnt_C, "Panel C: DOWN in siMTR4 x siZFC3H1 x siMTR4")
fmt_venn(cnt_D, "Panel D: UP in siMTR4 x siMTR4 x siZCCHC8")
fmt_venn(cnt_E, "Panel E: DOWN in siMTR4 x siMTR4 x siZCCHC8")

# ============================================================
# 5. Plot Venn diagrams
# ============================================================

# Colors matching project conventions (paxt_next_database.R)
COL_DEG_ALL  <- "#AAAAAA"    # light grey  - all DEGs
COL_DEG_UP   <- "#555555"    # dark grey   - upregulated in KD
COL_DEG_DOWN <- "#CCCCCC"    # white-grey  - downregulated in KD
COL_Z1       <- "#1898dd"    # blue        - siZFC3H1 PROMPTs
COL_MTR4     <- "#ff5733"    # orange-red  - siMTR4 PROMPTs
COL_Z8       <- "#7bb84f"    # green       - siZCCHC8 PROMPTs

fmt_p <- function(p) {
  e <- floor(log10(p))
  sprintf("p < 1\u00d710^%d", e)
}

make_venn_plot <- function(cnt, fills, labels, title = "") {
  p_val <- attr(cnt, "fisher_p")
  fit   <- euler(cnt, shape = "ellipse")
  plot(fit,
       fills      = list(fill = fills, alpha = 0.55),
       labels     = list(labels = labels, fontsize = 9),
       quantities = list(fontsize = 8, type = "counts"),
       edges      = list(col = "white", lwd = 1),
       main       = list(label = sprintf("%s\n%s", title, fmt_p(p_val)),
                         fontsize = 9)
  )
}

message("\nGenerating Venn diagrams...")

p_A <- make_venn_plot(cnt_A,
  fills  = c(COL_DEG_ALL, COL_Z1, COL_MTR4),
  labels = c("All DEGs", "siZFC3H1\nPROMPTs", "siMTR4\nPROMPTs"),
  title  = "A")

p_B <- make_venn_plot(cnt_B,
  fills  = c(COL_DEG_UP, COL_Z1, COL_MTR4),
  labels = c("UP in KD", "siZFC3H1\nPROMPTs", "siMTR4\nPROMPTs"),
  title  = "B")

p_C <- make_venn_plot(cnt_C,
  fills  = c(COL_DEG_DOWN, COL_Z1, COL_MTR4),
  labels = c("DOWN in KD", "siZFC3H1\nPROMPTs", "siMTR4\nPROMPTs"),
  title  = "C")

p_D <- make_venn_plot(cnt_D,
  fills  = c(COL_DEG_UP, COL_MTR4, COL_Z8),
  labels = c("UP in KD", "siMTR4\nPROMPTs", "siZCCHC8\nPROMPTs"),
  title  = "D")

p_E <- make_venn_plot(cnt_E,
  fills  = c(COL_DEG_DOWN, COL_MTR4, COL_Z8),
  labels = c("DOWN in KD", "siMTR4\nPROMPTs", "siZCCHC8\nPROMPTs"),
  title  = "E")

# ============================================================
# 6. Save output
# ============================================================

out <- file.path(plots_dir, "figS3_venn_deg_prompts.png")
png(out, width = 1400, height = 900, res = 120)
gridExtra::grid.arrange(
  p_A, p_B, p_C,
  p_D, p_E, grid::nullGrob(),
  ncol = 3,
  top = grid::textGrob(
    "DEG vs PROMPT overlap (siMTR4/siCon)",
    gp = grid::gpar(fontsize = 13, fontface = "bold")
  )
)
dev.off()
message("Saved: ", out)

# Individual panels for compositing
for (nm in c("A","B","C","D","E")) {
  p <- get(paste0("p_", nm))
  fn <- sprintf("figS3_venn_panel%s.png", nm)
  png(file.path(plots_dir, fn), width = 480, height = 480, res = 120)
  print(p)
  dev.off()
  message("Saved: ", fn)
}

message("\n=== FigS3 complete ===")

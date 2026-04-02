# HEADER ====================================================== =
# GEN2I
# Tue Jul  2 14:56:52 2024
# R version : R version 4.4.1 (2024-06-14)
# System : x86_64, linux-gnu
# ============================================================= =

# DESC  ======================================================= =
# Data : 
#
# Description :
#
# ClickUPID : 
# ============================================================= =

# SOURCE  ===================================================== =
library(data.table)
library(GenomicRanges)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
source("Analysis_Functions.R")
source("g2i_plotfig_lib.R")
source("g2i_apa_lib.R")
# ============================================================= =

# SETSEED  ==================================================== =
set.seed(123)
# ============================================================= =

#  READ I/O  ==================================================
# # LOOP PARAMETERS
# h5l <- list(WT=paste0("resources/processed_data/hg19/hic/dataKiernan/merge_ScrDrosha_rep1__rep1__rep2_inter_oe_KR_10kb.g2i.h5"))
# # 3D PARAMETERS ----
params <- NULL
params[["apa_size"]] <- 21
params[["d_bin_min"]] <- 21L
params[["d_bin_max"]] <- 100L
params[["h5path"]] <- list(WT = "data/hic_wt.h5")
params[["norm"]] <- c("raw","log10","quant","quant_na")
params[["cons_nm"]] <- "intra"
params[["save_dt"]] <- F
params[["quantile"]] <- as.integer(params[["apa_size"]]*params[["apa_size"]])
params[["comp_ref_cm"]] <- NULL
params[["comp_ref_ba"]] <- NULL
# params[["orientation"]] <- "both"
# params[["forceSym"]] <- F
params[["offset"]] <- 0L
params[["sample_size"]] <- 5000L
params[["plot_lst"]] <- c("APA_panel")
params[["fig_type"]] <- "png"
params[["width"]]    <- 800L
params[["height"]]   <- 800L
output <- "plots/fig3_apa"
# params[["nb_iter"]] <- 0
# 
# ## BAIT ANCHOR Params
# params[["bait_structure"]] <- list(
#   rn1 = list(
#     anchor = list(
#       combine = "and",
#       filter_1 = list(
#         aname = "h3k4me1_anc",
#         operator = "==",
#         value = T
#       ),
#       filter_2 = list(
#         aname = "mtr4_anc",
#         operator = "==",
#         value = T
#       )
#     )
#   )
# )
# params[["anchor_structure"]]  <- list(
#   rn1 = list(
#     tss = list()
#   )
# )
# 
# params[["plot_lst"]] <- c("APA_panel","intra","intra_norm","diffusion","diffusion_norm")
# params[["fig_type"]] <- "png"
# output <- "./test"
# # FILTERS ----
# params[["filter"]] <- list(lst=list(f1=list(aname = "eigen_pca1_N2",operator="<=",value=0)))
# # FUNCTION PARAMETERS ----
# 
# # RANGES & ANNOTATIONS ----
params.rname <- "tile_10kb"
# ============================================================= =


# SOURCE PARAMETERS ===========================================
##### PARAMETERS ----
#### GLOBAL VARIABLES
metric_type <- list(standard_metrics = c(CP = "CP", CS = "CS", BLS = "BLS",LEX = "LEX"),
                    diffusion_metrics = c(LV4 = "LV4", LV3 = "LV3", LV2 = "LV2", DLI = "DLI"))
coord_ofi <- c(CP = "CP", CS = "CS", BLS = "BLS", BRS = "BRS", ULS = "ULS", URS = "URS",LEX = "LEX", LV4 = "LV4", LV3 = "LV3", LV2 = "LV2", DLI = "DLI")
quantiles <- params$quantile
#### 3D PARAMETERS
params$orientation <- ifelse(is.null(params$orientation), "both", params$orientation)
params$forceSym <-  ifelse(is.null(params$forceSym), F, params$forceSym)
params$replace <-  ifelse(is.null(params$replace), F, params$replace)
params$nb_iter <-  ifelse(is.null(params$nb_iter), 1, params$nb_iter)
params$metrics_sel <-  ifelse(is.null(params$metrics_sel), c("standard_metrics","diffusion_metrics"), c("standard_metrics"))
params$d_bin_min <- max(params$offset,params$d_bin_min)
res <- as.integer(sub("tile_(.+)kb","\\1",params.rname))*1000

if(is.null(params$offset)) params$offset <- 0

#### APA PARAMETERS
mask <- matrix(1:(params$apa_size*params$apa_size), nrow=params$apa_size, byrow=F)
####GRAPHICAL PARAMETERS
params$plot_stat <- ifelse(is.null(params$plot_stat), T, params$plot_stat)
width <- ifelse(is.null(params$width), 800,  as.numeric(params$width))
height <- ifelse(is.null(params$height), 800,  as.numeric(params$height))

#### LOAD TILES RANGES + BAIT + ANCHOR ----
l_bait_go_tile <- list()
l_anchor_go_tile <- list()


cat("load Tile ranges and bait & anchor ...\n")
Go <- makeTile10kb()

mtr4_h3k4me1_gr <- subsetByOverlaps(
  loadBED("data/h3k4me1_anc.bed"),
  loadBED("data/mtr4_anc.bed")
)
l_bait_go_tile[["anchor_rn1"]]  <- subsetByOverlaps(Go, mtr4_h3k4me1_gr)
l_anchor_go_tile[["tss_rn1"]]   <- subsetByOverlaps(Go, loadBED("data/tss_act.bed"))

if(length(l_bait_go_tile) != length(l_anchor_go_tile)){
  if(length(l_bait_go_tile) < length(l_anchor_go_tile)) l_bait_go_tile <- rep(l_bait_go_tile, length(l_anchor_go_tile))
  else l_anchor_go_tile <- rep(l_anchor_go_tile, length(l_bait_go_tile))
}

#### LOAD H5 MATRIX
cat("load h5matrices ...\n")
l_cm <- load_h5_matrix(params$h5path)

#### EXTRACT CONTACT FREQUENCIES FOR EACH COUPLES/EACH MATRICES ----
cat("Extract HiC interactions ...\n")
out <- make_lflat_dt(l_bait_go_tile = l_bait_go_tile, 
                     l_anchor_go_tile = l_anchor_go_tile, 
                     hic_res = res, 
                     apa_size = params$apa_size, 
                     d_min= params$d_bin_min, 
                     forcesym = params$forceSym, 
                     d_max= params$d_bin_max, 
                     cons_nm = params$cons_nm, 
                     orientation = params$orientation,
                     sample_size = params$sample_size, 
                     replace = params$replace,
                     nb_iter = params$nb_iter,
                     l_cm = l_cm)
l_flat_dt <- out$lflat
nb_cpl_dt <- out$nb_cpl
dst_cpl_dt <- out$dst_cpl
out <- NULL; gc()

# si("ok")
#### COUPLES DISTRIBUTIONS STATISTICS ----
cat("Basic Statistics : \n")
l_gg <- list()
cat(paste0("Number of couples \n"))
for (i in 1:nrow(nb_cpl_dt)) cat(paste0(gsub("\n"," ",nb_cpl_dt[i,"name"])," : ", nb_cpl_dt[i,"value"],"\n"))
l_gg[["APA_Statistics.Number_of_couples"]] <- ggpubr::ggbarplot(nb_cpl_dt, x = "name",y = "value", fill ="name", title = "Number of pair of couples") 
cat(paste0("Distance betweencouples \n"))
l_gg[["APA_Statistics.Distance_between_couples"]] <- ggpubr::ggviolin(dst_cpl_dt, x = "name",y = "value", fill ="name", title = "Distance distribution between couples") 


#### APA TO FLAT & USABLE MATRICES ----
if (length(l_flat_dt) == 0) stop("No couples to process — check bait/anchor overlap and H5 matrix.", call. = FALSE)
cat("Turn matrices in flat data table...\n")
flat_dt <- aggreg_lflat(l_flat_dt)
rm(l_flat_dt);gc()

cat("Remove diagonal on flat dt ...\n")
remove_diag(flat_dt, params$offset, params$apa_size)

cat("Get APA metrics & compute them for each couples on flat dt ...\n")
coord_l <- get_apa_metrics_coordinate(params$apa_size, mask,ctloffset = 6)
setApaMetrics(flat_dt, coord_l)

#### FISHER TEST: PROMPTs ENRICHMENT VS Hi-C QUARTILES (enhancer +/- MTR4) ----
cat("Fisher test: PROMPTs enrichment vs Hi-C quartiles (enhancer +/- MTR4) ...\n")

# Strand-aware TSS loader (col 6 = strand, required for correct upstream extension)
loadBED6 <- function(path) {
  dt <- data.table::fread(path, header = FALSE)
  gr <- GenomicRanges::GRanges(
    seqnames = as.character(dt[[1]]),
    ranges   = IRanges::IRanges(as.integer(dt[[2]]) + 1L, as.integer(dt[[3]])),
    strand   = if (ncol(dt) >= 6) as.character(dt[[6]]) else "*"
  )
  if (ncol(dt) >= 4) gr$gene_id <- as.character(dt[[4]])
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI")
  gr
}

# P56: active+nodiv TSS × per-TSS max ULS_V × 5kb PROMPTs × Q4
# YES = H3K4ME1+MTR4 anchors (MTR4 RDS), NO = H3K4ME1-only anchors (H3K4ME1 RDS)
fish_tss_gr <- loadBED6("data/tss.bed")
fish_tss_act <- loadBED6("data/tss_act.bed")
fish_ncRNA   <- loadBED6("data/ncRNA_MTR4kd.bed")
nodiv_gr     <- readRDS("data/TSSnodiv.gr")
suppressWarnings(GenomeInfoDb::seqlevelsStyle(nodiv_gr) <- "NCBI")

bait_is_active <- integer(length(fish_tss_gr))
bait_is_active[unique(queryHits(findOverlaps(fish_tss_gr, fish_tss_act, ignore.strand = TRUE)))] <- 1L
is_nodiv_f <- as.integer(fish_tss_gr$gene_id %in% nodiv_gr$gene_id)
idx_act_nodiv <- which(bait_is_active == 1L & is_nodiv_f == 1L)
cat(sprintf("TSS: %d total, %d active+nodiv\n", length(fish_tss_gr), length(idx_act_nodiv)))

# PROMPTs annotation: strand-aware 5kb upstream of each TSS
bait_has_prompts <- integer(length(fish_tss_gr))
ov_prom <- findOverlaps(extendGR(fish_tss_gr, upstream = 5000L, downstream = 0L),
                        fish_ncRNA, ignore.strand = TRUE)
bait_has_prompts[unique(queryHits(ov_prom))] <- 1L
cat(sprintf("TSS with PROMPTs (strand-aware 5kb): %d / %d\n",
    sum(bait_has_prompts), length(fish_tss_gr)))
prom_dt <- data.table(i_IDX = seq_along(fish_tss_gr), has_p = bait_has_prompts)

# Anchor BEDs
anc_mtr4    <- loadBED("data/mtr4_anc.bed")
anc_h3k4me1 <- loadBED("data/h3k4me1_anc.bed")

# --- "YES" GROUP: MTR4 RDS — H3K4ME1+MTR4 anchors, per-TSS max ULS_V --------
hic_mtr4     <- readRDS("data/hic_mtr4_metrics.rds")
mtr4_metrics <- as.data.table(hic_mtr4$METRICS.df); hic_mtr4 <- NULL; gc()

mtr4_has_h3k4me1 <- integer(length(anc_mtr4))
mtr4_has_h3k4me1[unique(queryHits(findOverlaps(anc_mtr4, anc_h3k4me1, ignore.strand = TRUE)))] <- 1L
idx_yes_anc <- which(mtr4_has_h3k4me1 == 1L)
cat(sprintf("MTR4 anchors with H3K4ME1 (yes): %d / %d\n", length(idx_yes_anc), length(anc_mtr4)))

f_yes_raw <- mtr4_metrics[a_IDX %in% idx_yes_anc & i_IDX %in% idx_act_nodiv & !is.nan(ULS_V)]
tss_yes   <- f_yes_raw[, .(score = max(ULS_V)), by = i_IDX]
tss_yes   <- merge(tss_yes, prom_dt, by = "i_IDX")
tss_yes[, hic_Q := as.integer(ntile(score, 4L))]
mtr4_metrics <- NULL; f_yes_raw <- NULL; gc()
cat(sprintf("YES TSS: %d, PROMPTs: %d\n", nrow(tss_yes), sum(tss_yes$has_p)))

# --- "NO" GROUP: H3K4ME1 RDS — H3K4ME1-only anchors, per-TSS max ULS_V ------
hic_h3k4me1     <- readRDS("data/hic_h3k4me1_metrics.rds")
h3k4me1_metrics <- as.data.table(hic_h3k4me1$METRICS.df); hic_h3k4me1 <- NULL; gc()

h3k4me1_has_mtr4 <- integer(length(anc_h3k4me1))
h3k4me1_has_mtr4[unique(queryHits(findOverlaps(anc_h3k4me1, anc_mtr4, ignore.strand = TRUE)))] <- 1L
idx_no_anc <- which(h3k4me1_has_mtr4 == 0L)
cat(sprintf("H3K4ME1 anchors without MTR4 (no): %d / %d\n", length(idx_no_anc), length(anc_h3k4me1)))

f_no_raw <- h3k4me1_metrics[a_IDX %in% idx_no_anc & i_IDX %in% idx_act_nodiv & !is.nan(ULS_V)]
tss_no   <- f_no_raw[, .(score = max(ULS_V)), by = i_IDX]
tss_no   <- merge(tss_no, prom_dt, by = "i_IDX")
tss_no[, hic_Q := as.integer(ntile(score, 4L))]
h3k4me1_metrics <- NULL; f_no_raw <- NULL; gc()
cat(sprintf("NO  TSS: %d, PROMPTs: %d\n", nrow(tss_no), sum(tss_no$has_p)))

# --- FISHER TEST (per-TSS, within-group Q4, real computation both rows) ------
plot_dt_yes <- as.data.table(fisher_ntiles(tss_yes, "hic_Q", "has_p", NTILE = 4L, NAME = "yes"))
plot_dt_no  <- as.data.table(fisher_ntiles(tss_no,  "hic_Q", "has_p", NTILE = 4L, NAME = "no"))

final_fish_dt <- rbind(plot_dt_no, plot_dt_yes)
final_fish_dt[, COL       := XNAME]
final_fish_dt[, ROW       := YNAME]
final_fish_dt[, intersect := ""]
final_fish_dt[, fm        := paste0("lfc=", round(lfcb, 2), " pv=", pv)]

gg_fisher <- fisher_plot_asym_pdf(
  DF           = final_fish_dt,
  main         = "PROMPTs enrichment vs Hi-C quartile by enhancer MTR4 status (Fig 3D)",
  range        = c(-2, 2),
  x_axis_title = "Hi-C Contact Quartile",
  y_axis_title = "Enhancer MTR4 status"
)
# Lighter color palette: skip darkest extremes of RdBu (positions 3-9 of 11)
light_cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")[3:9]))(10)
gg_fisher <- suppressWarnings(
  gg_fisher + scale_fill_gradientn(
    colours = light_cols, limits = c(-2, 2), breaks = c(-2, 0, 2), labels = c(-2, 0, 2),
    guide   = guide_colorbar(barwidth=1, barheight=5, ticks=FALSE,
                title="_____\nScaled\nLog2\nOddsRatio\n_____",
                title.theme=element_text(size=8, vjust=0.9, face="bold", angle=0))
  )
)
l_gg[["Fisher.PROMPTs_vs_HiC_quartiles"]] <- gg_fisher
ggsave(paste0(output, ".Fisher.PROMPTs_vs_HiC_quartiles.png"),
       gg_fisher, width = width / 80, height = height / 80, dpi = 300)
tss_yes <- NULL; tss_no <- NULL; gc()

for (type in c("CS","CP","DLI")){
  for (norm in c("raw","quant")) {
    if (norm =="raw") {
      flat_dt[,quartiles:=ntile(get(type),5)] 
      flat_dt[,bait:="anchor_rn1"] 
      flat_dt[,anchor:="tss_rn1"] 
      flat_dt[,bait:=paste0(bait,"_",quartiles-1)] 
      flat_dt[,anchor:=paste0(anchor,"_",quartiles-1)] 
      avg_flat_dt <- flatToApaDt(flat_dt[quartiles!=1])
    } else {
      #### APA MATRICES NORMALIZATIONS ----
      cat("Normalization step, with Quant : \n")
      flat_q_dt <- data.table::copy(flat_dt)
      quantilize_flat_dt(flat_q_dt)
      flat_q_dt[,quartiles:=ntile(get(type),5)] 
      flat_q_dt[,bait:="anchor_rn1"] 
      flat_q_dt[,anchor:="tss_rn1"] 
      flat_q_dt[,bait:=paste0(bait,"_",quartiles-1)] 
      flat_q_dt[,anchor:=paste0(anchor,"_",quartiles-1)] 
      avg_flat_dt <- flatToApaDt(flat_q_dt[quartiles!=1])
    }
    cat("Convert flat dt to matrices before plotting ...\n")
    out <- convert_flat_vector_to_matrices(avg_flat_dt, apa_size = params$apa_size)
    l_m <- out[["l_m"]]; minl <- out[["min"]];  maxl <- out[["max"]]
    out <- NULL; gc()
    
    #### APA PLOT ----
    cat("Make APA's panels \n")
    apa_panels <- gg_apa_panels(l_m = l_m, norm = norm, min_l = minl, max_l = maxl, params = params, n = params$apa_size, width = width, height = height, output = paste0(output, ".", type))
    l_gg[[paste0("APA_panel.", type, ".", norm, ".Scaled")]]       <- apa_panels$scaled
    l_gg[[paste0("APA_panel.", type, ".", norm, ".Unscaled")]]     <- apa_panels$unscaled
    l_gg[[paste0("APA_panel.", type, ".", norm, ".MinMaxScaled")]] <- apa_panels$minmax_scaled
  }
}

#### SAVE PLOTS ----
if (!is.null(l_gg) && length(l_gg) > 0) {
  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  for (i in names(l_gg)) {
    gg <- ggPostTreatment(l_gg[[i]], params)
    if ("png" %in% params$fig_type)
      ggsave(paste0(output, ".", i, ".png"), gg, width = width / 80, height = height / 80, dpi = 300)
    if ("svg" %in% params$fig_type)
      ggsave(paste0(output, ".", i, ".svg"), gg, width = width / 80, height = height / 80)
  }
}

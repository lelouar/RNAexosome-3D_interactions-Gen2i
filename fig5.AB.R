quiet <- function(...) {suppressWarnings(suppressMessages(...))}

# Libraries loading ----
cat("Loading library:... ")
quiet(library(GenomicRanges))
quiet(library(dplyr))
quiet(library(RColorBrewer))
quiet(library(data.table))
quiet(library(ggplot2))
quiet(library(ggpubr))
quiet(library(officer))
quiet(library(rvg))
quiet(library(svglite))
quiet(library(ggsci))

# SOURCE APA AND GRAPHICAL REALTED LIBRARIES ----
source("g2i_plotfig_lib.R")
source("g2i_apa_lib.R")

library(ggthemes)
theme_set(theme_tufte())
theme_update(plot.title = element_text(hjust = 0.5))
theme_update(plot.subtitle = element_text(hjust = 0.5))

#Load Required Libraries
quiet(library(gridExtra))

cat(" => ok\n")

# ---- Parameters (hardcoded from config/src/yaml/fig.yaml) ----
knitr::opts_knit$set(root.dir = here::here())
output <- "plots/fig5_apa_top5p"

params.rname   <- "tile_10kb"
params.aname   <- "chr_bin"

params <- list(
  anchor_structure = list(
    tss_act = list(
      tss = list(
        filter_1 = list(aname = "active", operator = "==", value = TRUE)
      )
    )
  ),
  bait_structure = list(
    rad21 = list(
      rad21 = list(
        filter_1 = list(aname = "rad21_mtr4_zscore_100tile", operator = ">=", value = 95L)
      )
    )
  ),
  h5path = list(
    WT     = "data/hic_wt.h5",
    MTR4KD = "data/hic_mtr4kd.h5"
  ),
  apa_size    = 21L,
  offset      = 0L,
  d_bin_min   = 21L,
  d_bin_max   = 100L,
  cons_nm     = "tad_rao_hela",
  orientation = "both",
  forceSym    = FALSE,
  replace     = TRUE,
  nb_iter     = 50L,
  sample_size = 5000L,
  quantile    = 441L,
  norm        = c("raw", "quant", "quant_na"),
  metrics_sel = c("standard_metrics", "diffusion_metrics"),
  plot_stat   = TRUE,
  save_dt     = FALSE,
  plot_lst    = c("APA_panel"),
  fig_type    = "png",
  width       = 800L,
  height      = 800L,
  sqlv        = "NCBI",
  threads     = 28L
)

data.table::setDTthreads(threads = params$threads)
width  <- as.numeric(params$width)
height <- as.numeric(params$height)


skip_go <- TRUE
  # Set option scipen for resolution
  options(scipen=100)
  
  #### PARAMETERS ----
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

  # Tile 10kb genome (chr_bin used as H5 matrix row/col index)
  Go <- makeTile10kb()

  # --- Bait: RAD21 peaks in top 5% by MTR4 enrichment (rad21_mtr4_zscore_100tile >= 95) ---
  rad21_gr <- loadBED("data/rad21.bed")
  mtr4_bw  <- "data/mtr4.bw"
  rad21_ext <- suppressWarnings(GenomicRanges::trim(GenomicRanges::resize(rad21_gr, 2001L, fix = "center")))
  bw_gr     <- rtracklayer::import.bw(mtr4_bw, which = rad21_ext)
  ol        <- GenomicRanges::findOverlaps(rad21_ext, bw_gr)
  mtr4_mean <- numeric(length(rad21_ext))
  if (length(ol) > 0) {
    scores <- tapply(bw_gr$score[S4Vectors::subjectHits(ol)], S4Vectors::queryHits(ol), mean)
    mtr4_mean[as.integer(names(scores))] <- as.numeric(scores)
  }
  mtr4_z    <- as.numeric(scale(mtr4_mean))
  mtr4_tile <- dplyr::ntile(mtr4_z, 100L)
  Go_bait   <- rad21_gr[mtr4_tile >= 95L]
  l_bait_go_tile[["rad21"]] <- IRanges::subsetByOverlaps(Go, Go_bait)

  # --- Anchor: active TSS (tss_act.hg19.bed — all rows are active) ---
  tss_act_gr <- loadBED("data/tss_act.bed")
  l_anchor_go_tile[["tss_tss_act"]] <- IRanges::subsetByOverlaps(Go, tss_act_gr)

  rm(rad21_gr, rad21_ext, bw_gr, ol, mtr4_mean, mtr4_z, mtr4_tile, Go_bait, tss_act_gr)
  
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
  
  #### COUPLES DISTRIBUTIONS STATISTICS ----
  cat("Basic Statistics : \n")
  l_gg <- list()
  cat(paste0("Number of couples \n"))
  for (i in 1:nrow(nb_cpl_dt)) cat(paste0(gsub("\n"," ",nb_cpl_dt[i,"name"])," : ", nb_cpl_dt[i,"value"],"\n"))
  l_gg[["APA_Statistics.Number_of_couples"]] <- ggpubr::ggbarplot(nb_cpl_dt, x = "name",y = "value", fill ="name", title = "Number of pair of couples") 
  cat(paste0("Distance betweencouples \n"))
  l_gg[["APA_Statistics.Distance_between_couples"]] <- ggpubr::ggviolin(dst_cpl_dt, x = "name",y = "value", fill ="name", title = "Distance distribution between couples") 
  
  #### APA TO FLAT & USABLE MATRICES ----
  if (length(l_flat_dt) == 0){
    cat("length(l_flat_dt) == 0, there no couple to process ...")
    gg <- NULL
  } else {
    
    cat("Turn matrices in flat data table...\n")
    flat_dt <- aggreg_lflat(l_flat_dt)
    rm(l_flat_dt);gc()
    
    cat("Remove diagonal on flat dt ...\n")
    remove_diag(flat_dt, params$offset, params$apa_size)
    
    cat("Get APA metrics & compute them for each couples on flat dt ...\n")
    coord_l <- get_apa_metrics_coordinate(params$apa_size, mask,ctloffset = 6)
    setApaMetrics(flat_dt, coord_l)
    
    #### APA MATRICES NORMALIZATIONS ----
    for (norm in params$norm){ # RAW HAS TO BE FIRST - LOG SECOND AND QUANT HAS TO BE LAST - see copy(dt) !!!!
      if (norm == "raw"){
        cat("Normalization step, starting with Raw : \n")
        avg_flat_dt <- flatToApaDt(flat_dt)
      } else if (norm == "quant"){
        cat("Normalization step, with Quant : \n")
        flat_q_dt <- data.table::copy(flat_dt)
        quantilize_flat_dt(flat_q_dt)
        avg_flat_dt <- flatToApaDt(flat_q_dt)
      } else if (norm == "quant_na") {
        cat("Normalization step, with Quant NA removed : \n")
        flat_q_dt <- data.table::copy(flat_dt)
        quantilize_na_flat_dt(flat_q_dt)
        avg_flat_dt <- flatToApaDt(flat_q_dt)
      } else {
        cat("Normalization step, with Quant NA removed : \n")
        avg_flat_dt <- flatToApaDt(transform_flat_dt(data.table::copy(flat_dt), norm = norm, rangeMinMax = params$rangeMinMax,metrics_names = names(coord_l)))
      }
      cat("Convert flat dt to matrices before plotting ...\n")
      out <- convert_flat_vector_to_matrices(avg_flat_dt, apa_size = params$apa_size)
      l_m <- out[["l_m"]]; minl <- out[["min"]];  maxl <- out[["max"]]
      out <- NULL; gc()
      
      saveRDS(object = l_m, paste0(output,".matrix.",norm,".apards"))
      
      #### APA PLOT ----
      cat("Make APA's panels \n")
      apa_panels <- gg_apa_panels(l_m = l_m, norm = norm, min_l = minl, max_l = maxl, params = params, n = params$apa_size, width = width , height = height, output = output)
      l_gg[[paste0("APA_panel.", norm, ".Scaled")]]       <- apa_panels$scaled
      l_gg[[paste0("APA_panel.", norm, ".Unscaled")]]     <- apa_panels$unscaled
      l_gg[[paste0("APA_panel.", norm, ".MinMaxScaled")]] <- apa_panels$minmax_scaled
      
      #### STATISTICAL PLOT ----
      # Compute Each sub matrices of interest
      cat("Prepare flat dt for statistical analyzes \n")
      # If nb_iter > 1 it means that we performed sampling nb_iter times and we want to observe statistics based on this iterated sampling.
      if (params$nb_iter > 1 & grepl("quant",norm)){# quant with sampling
        mean_flat_dt <- flat_q_dt[, lapply(.SD, function(x) {m = mean(x,na.rm=T);m[!is.finite(m)] <- 0;m}), .SDcols = names(flat_q_dt) %in% names(coord_ofi), by = .(group, bait, anchor, cm)]
        dt <- get_apa_dt(apa_dt = mean_flat_dt, coord_l = coord_l, norm = norm,  save_dt = params$save_dt, output = output, params = params)
      } else if (params$nb_iter > 1) { # raw or whatever but with sampling
        mean_flat_dt <- flat_dt[, lapply(.SD, function(x) {m = mean(x,na.rm=T);m[!is.finite(m)] <- 0;m}), .SDcols = names(flat_dt) %in% names(coord_ofi), by = .(group, bait, anchor, cm)]
        dt <- get_apa_dt(apa_dt = mean_flat_dt, coord_l = coord_l, norm = norm,  save_dt = params$save_dt, output = output, params = params)
      } else if (grepl("quant",norm)){ # No sampling, quant
        dt <- get_apa_dt(apa_dt = flat_q_dt, coord_l = coord_l, norm = norm,  save_dt = params$save_dt, output = output, params = params)
      } else { # No sampling, raw or other
        dt <- get_apa_dt(apa_dt = flat_dt, coord_l = coord_l, norm = norm,  save_dt = params$save_dt, output = output, params = params)
      }
      # In all cases, compute normalisation for all APA metrics of interest
      
      cat("APA metric statistics and condition or bait/anchor comparisons \n")
      
      if (params$plot_stat){
        if(is.null(params$compare$method)) params$compare$method = "wilcox.test"
        
        for (m_type in params$metrics_sel) {
          for (nrm in c("Normalized","Unnormalized")) {
            
            if (m_type == "standard_metrics" & nrm == "Normalized") {
              col_sel <- c("CP.norm", "CS.norm", "BLS.norm")
              comp_metrics <- combn(x = as.character(unique(dt[variable %in% col_sel, ]$variable)), m = 2, simplify = F)
              select <- unlist(lapply(comp_metrics, function(x) "CP.norm" %in% x))
              comp_metrics <- comp_metrics[select]
            } else if (m_type == "standard_metrics" & nrm == "Unnormalized") {
              col_sel <- c("CP", "CS", "BLS", "ctrl")
              comp_metrics <- combn(x = as.character(unique(dt[variable %in% col_sel, ]$variable)), m = 2, simplify = F)
              select <- unlist(lapply(comp_metrics, function(x) "ctrl" %in% x))
              comp_metrics <- comp_metrics[select]
            } else if (m_type == "diffusion_metrics" & nrm == "Unnormalized") {
              col_sel <-  c("CP", "LV2", "LV3", "LV4","DLI")
              comp_metrics <- combn(x = as.character(unique(dt[variable %in% col_sel, ]$variable)), m = 2, simplify = F)
              select <- unlist(lapply(comp_metrics, function(x) "CP" %in% x))
              comp_metrics <- comp_metrics[select]
            } else if (m_type == "diffusion_metrics" & nrm == "Normalized") {
              col_sel <-c("CP.norm", "LV2.norm", "LV3.norm", "LV4.norm","DLI.norm")
              comp_metrics <- combn(x = as.character(unique(dt[variable %in% col_sel, ]$variable)), m = 2, simplify = F)
              select <- unlist(lapply(comp_metrics, function(x) "CP.norm" %in% x))
              comp_metrics <- comp_metrics[select]
            }
            dt_noOutli <- data.table::copy(dt)
            dt_noOutli[, value := remove_outliers(value), by = .(cm, c_name, variable)]
            for (which_dt in c("dt","dt_noOutli")) {
              outli_type <- ifelse(which_dt=="dt","With_outlier","No_outliers")
              mdt <- get(which_dt)
              gg <- ggboxplot(mdt[variable %in% col_sel, ], x="variable", y="value", col="cm", fill = "c_name" ) + scale_fill_brewer(palette = 2)
              by_type <- ifelse(length(unique(mdt$cm)) > length(unique(mdt$c_name)), "cm", "c_name")
              
              l_gg[[paste0("Intra_group_comparisons.",m_type,".",nrm,".",outli_type,".",norm)]] <- gg + facet_wrap(paste0(".~",by_type)) + 
                stat_compare_means(aes(label = ..p.signif..),comparisons = comp_metrics,paired = T,method = params$compare$method) +
                ggtitle(paste0("Comparisons of APA standard metrics by : ",ifelse(by_type=="cm","HiC Conditions","Bait/Anchor couples"),"\n",nrm," signal"))
              
              l_gg[[paste0("Intra_group_comparisons.",m_type,".",nrm,"_no_pval.",outli_type,".",norm)]] <- gg + facet_wrap(paste0(".~",by_type)) +
                ggtitle(paste0("Comparisons of APA standard metrics by : ",ifelse(by_type=="cm","HiC Conditions","Bait/Anchor couples"),"\n",nrm," signal"))  
              
            }
          }
        }
        
        
        rm(dt_noOutli)
        
        
        # INTER GROUP COMPARISONS
        for (m_type in params$metrics_sel) {
          
          metric_unorm <- c(names(metric_type[[m_type]]),if(m_type=="standard_metrics") "ctrl" else NULL)
          metric_norm <- paste0(names(metric_type[[m_type]]),".norm")
          # plot normalised value
          dt.norm <- dt[variable %in% metric_norm, ]
          dt.unorm <- dt[variable %in% metric_unorm, ]
          
          for(normalize in c("Normalized", "Unnormalized")){
            if(normalize == "Normalized") dtn <- dt.norm
            else dtn <- dt.unorm
            
            if(!is.null(params$comp_ref_ba)){
              ref_comp_ba <- unique(grep(paste0("^",params$comp_ref_ba,"_vs_.+|.+_vs_",params$comp_ref_ba), dt$c_name,perl = T,value = T))
              if (!length(ref_comp_ba)) ref_comp_ba <- NULL # We can't grep anything
            } else {
              ref_comp_ba <- NULL
            }
            if(!is.null(params$comp_ref_cm)){
              ref_comp_cm <- unique(grep(params$comp_ref_cm, dt$cm,value = T))
              if (!length(ref_comp_cm)) ref_comp_cm <- NULL # We can't grep anything
            } else {
              ref_comp_cm <- NULL
            }
            # We don't use comp_ref here as when inputting multiple bait/anchor, the idea is to compare everyone
            # Not obvioulsy the case in multiple h5 matrices when we often have a WT reference.
            if(length(unique(dt$c_name))>1){
              myComp <- combn(x = as.character(unique(dt$c_name)), m = 2, simplify = F)
              for(cm.nm in dtn[, unique(cm)]){
                gg <- ggboxplot(dtn[cm == cm.nm], x="c_name", y="value", fill = "c_name") + scale_fill_brewer(palette = 2)+ xlab(NULL)
                gg <- gg + facet_wrap(cm~variable) + theme(axis.text.x=element_blank()) +
                  ggtitle(paste0("Comparisons of bait/anchor \n", normalize," ", cm.nm, " signal"))
                if (!is.null(ref_comp_ba)) gg_pv <- gg + stat_compare_means(aes(label = ..p.signif..), ref.group  = ref_comp_ba, paired = F, method = params$compare$method)
                else gg_pv <- gg + stat_compare_means(aes(label = ..p.signif..), comparisons = myComp, paired = F, method = params$compare$method)
                l_gg[[paste0("Bait_anchor_comparisons.",m_type,".",normalize,".HiC_", cm.nm,".With_outliers.",norm)]] <- gg_pv
                l_gg[[paste0("Bait_anchor_comparisons.",m_type,".",normalize,".HiC_", cm.nm,".With_outliers_no_pval.",norm)]] <- gg
              }
            }
            
            if(length(unique(dtn$cm))>1){
              myComp <- combn(x = as.character(unique(dtn$cm)), m = 2, simplify = F)
              for(ba_name in dtn[, unique(c_name)]){
                gg <- ggboxplot(dtn[ba_name == c_name], x="cm", y="value", col = "cm")+ xlab(NULL)
                gg <- gg + facet_wrap(.~variable) +
                  theme(axis.text.x=element_blank()) +
                  ggtitle(paste0("Comparisons of HiC conditions over ",ba_name," couples \n",normalize," HiC signals"))
                if (!is.null(ref_comp_cm)) gg_pv <- gg + stat_compare_means(aes(label = ..p.signif..),ref.group  = ref_comp_cm,paired = T,method = params$compare$method)
                else gg_pv <- gg + stat_compare_means(aes(label = ..p.signif..),comparisons = myComp,paired = T,method = params$compare$method)
                l_gg[[paste0("HiC_Conditions_comparisons.",m_type,".",normalize,".", ba_name,".With_outliers.",norm)]] <- gg_pv
                l_gg[[paste0("HiC_Conditions_comparisons.",m_type,".",normalize,".", ba_name,".With_outliers_no_pval.",norm)]] <- gg
              }
            }
            
            #remove outliers and replot
            dtn[, value := remove_outliers(value), by = .(cm, c_name, variable)]
            
            if(length(unique(dtn$c_name))>1){
              myComp <- combn(x = as.character(unique(dtn$c_name)), m = 2, simplify = F)
              for(cm.nm in dtn[, unique(cm)]){
                gg <- ggboxplot(dtn[cm == cm.nm], x="c_name", y="value", fill = "c_name") + scale_fill_brewer(palette = 2) + xlab(NULL)
                gg <- gg + facet_wrap(cm~variable) +
                  theme(axis.text.x=element_blank()) +
                ggtitle(paste0("Comparisons of bait/anchor \n",normalize," ",cm.nm,"HiC signals\n","Outliers removed")) + xlab(NULL)
                if (!is.null(ref_comp_ba)) gg_pv <- gg + stat_compare_means(aes(label = ..p.signif..),ref.group  = ref_comp_ba,paired = F,method = params$compare$method)
                else gg_pv <- gg + stat_compare_means(aes(label = ..p.signif..),comparisons = myComp,paired = F,method = params$compare$method)
                
                l_gg[[paste0("Bait_anchor_comparisons.",m_type,".",normalize,".HiC_", cm.nm,".No_outliers_no_pval.",norm)]] <- gg
                l_gg[[paste0("Bait_anchor_comparisons.",m_type,".",normalize,".HiC_", cm.nm,".No_outliers.",norm)]] <- gg_pv
              }
            }
            
            if(length(unique(dtn$cm))>1){
              myComp <- combn(x = as.character(unique(dtn$cm)), m = 2, simplify = F)
              for(ba_name in dtn[, unique(c_name)]){
                gg <- ggboxplot(dtn[ba_name == c_name], x="cm", y="value", col = "cm") + xlab(NULL)
                gg <- gg + facet_wrap(.~variable) +
                  theme(axis.text.x=element_blank()) +
                  ggtitle(paste0("Comparisons of HiC conditions over ",ba_name," couples \n",normalize," HiC signal\n","Outliers removed"))
                if (!is.null(ref_comp_cm)) gg_pv <- gg + stat_compare_means(aes(label = ..p.signif..),ref.group  = ref_comp_cm,paired = F,method = params$compare$method)
                else gg_pv <- gg + stat_compare_means(aes(label = ..p.signif..),comparisons = myComp,paired = F,method = params$compare$method)
                
                l_gg[[paste0("HiC_Conditions_comparisons.",m_type,".",normalize,".", ba_name,".No_outliers.",norm)]] <- gg_pv
                l_gg[[paste0("HiC_Conditions_comparisons.",m_type,".",normalize,".", ba_name,".No_outliers_no_pval.",norm)]] <- gg
              }
            }
          }
        }
      }
    }
    gg <- l_gg
  }
# SAVING PLOT ----
# if(is.ggplot(gg))
# gg <- list(ONEFIG = gg)
if ( ! skip_go ) {
  cat(str(Go_dt))
  if(is.data.table(Go_dt)){
    if("csv" %in% params$fig_type)
      fwrite(Go_dt[, lapply(.SD, as.character)], paste0(output,".data.csv"), sep = "\t")
  }
}
if(!is.null(gg)){
  l_gg <- gg
  if (length(l_gg) > 1){
    for(i in names(l_gg)){
      if (grepl("CrossTalk",i)) {
        saveRDS(l_gg[[i]], paste0(output,".",i,".ct"))
      } else {
        gg <- l_gg[[i]]
        ## add x/y lim  ----
        gg <- ggPostTreatment(gg, params)
        
        if("rds" %in% params$fig_type)
          saveRDS(gg, paste0(output,".",i,".dyn"))
        if("svg" %in% params$fig_type)
          ggsave(paste0(output,".",i,".svg"), gg, width = width/80, height = height/80)
        if("png" %in% params$fig_type)
          ggsave(paste0(output,".",i,".png"), gg, width = width/80, height = height/80)
        if("pptx" %in% params$fig_type)
          pptx <- read_pptx() %>%
            add_slide() %>%
            ph_with(value = dml(ggobj = l_gg[[i]]), location = ph_location(width = 9, height = 4.95))%>%
            print(target = file.path(paste0(output,".",i,".pptx")))
      }
    }
  } else {
    i <- names(l_gg) # Keep names for html dongles
    if (grepl("CrossTalk",i)) {
      saveRDS(l_gg[[i]], paste0(output,".",i,".ct"))
    } else {
      gg <- ggPostTreatment(l_gg[[1]], params)
      if("rds" %in% params$fig_type)
        saveRDS(gg, paste0(output,".",i,".dyn"))
      if("svg" %in% params$fig_type)
        ggsave(paste0(output,".",i,".svg"), gg, width = width/80, height = height/80)
      if("png" %in% params$fig_type)
        ggsave(paste0(output,".",i,".png"), gg, width = width/80, height = height/80, dpi = 300)
      if("pptx" %in% params$fig_type)
        pptx <- read_pptx() %>%
          add_slide() %>%
          ph_with(value = dml(ggobj = l_gg[[i]]), location = ph_location(width = 9, height = 4.95))%>%
          print(target = file.path(paste0(output,".",i,".pptx")))
    }
  }
}

source("g2i.R")

###  APA PREPROCESS ----
#' Title
#'
#' @param h5path
#'
#' @return
#' @export
#'
#' @examples
load_h5_matrix <- function(h5path){
  l_cm <- list()
  # Load H5 matrix
  if (!is.null(h5path)) {
    for(h5.nm in names(h5path)){
      # get h5 name info by h5path
      cm_info <- data.table::data.table(rhdf5::h5ls(h5path[[h5.nm]]))
      if (! nrow(cm_info[otype=="H5I_GROUP","name"])) l_cm[[h5.nm]] <- sapply(cm_info$name, function(nm) { HDF5Array::HDF5Array(h5path[[h5.nm]], nm)})
      else l_cm[[h5.nm]] <- sapply(as.vector(unlist(cm_info[otype=="H5I_GROUP","name"])), function(nm) {print(nm);HDF5Array::TENxMatrix(h5path[[h5.nm]], nm)})
    }
  } else {
    # TODO : put loop here to load cm in l_cm
    mat_nm <- paste0(cond,"_",rep,"_",norm,"_",as.integer(res))
    filt_cm <- filterCmDB(cmname=list(in_=list(value=mat_nm,strict=T)))
    l_cm[["to_dev"]] <- dumpCM(filt_cm)
  }
  l_cm
}



#' Title
#'
#' @param Go
#' @param anchor
#'
#' @return
#' @export
#'
#' @examples
load_anchor <- function(Go, anchor){
  # Load anchor genome
  l_anchor_go_tile <- list()
  for(anchor.nm in anchor){
    l_anchor_go_tile[[anchor.nm]] <- subset(Go, mcols(Go)[[anchor.nm]][[anchor.nm]] == T)
  }
  l_anchor_go_tile
}

#' Title
#'
#' @param Go
#' @param bait
#'
#' @return
#' @export
#'
#' @examples
load_bait <- function(Go, bait){
  # Load bait genome
  l_bait_go_tile <- list()
  for(bait.nm in bait){
    l_bait_go_tile[[bait.nm]] <- subset(Go, mcols(Go)[[bait.nm]][[bait.nm]] == T)
  }
  l_bait_go_tile
}



### EXTRACT SUBMATRICES ----

#' Title
#'
#' @param l_bait_go_tile
#' @param l_anchor_go_tile
#' @param hic_res
#' @param apa_size
#' @param d_min
#' @param d_max
#' @param cons_nm
#' @param l_cm
#' @param orientation 
#' @param control Integer(1). Shift in bp of the control position from Anchor and Bait ranges.
#' @param forcesym Should we force the Anchor/Bait to be symetric ? See g2ismk::plot_apa_smk. Default = F.
#' @param sample_size Size of sampling. Default NULL (No sampling perfomred, all couples are used)
#' @param replace Should we use replace argument when usiung the sample() function ? Default = F
#' @param nb_iter Number of sampling iteration performed on whole Anchor/Bait interaction DT. By default, nb_iteratio = 1 and no sampling is performed.
#'
#' @return
#' @export
#'
#' @examples
make_lflat_dt <- function(l_bait_go_tile = NULL, l_anchor_go_tile = NULL,
                          hic_res = NULL, apa_size = NULL, d_min = NULL, d_max = NULL, cons_nm = "intra",
                          orientation = NULL, control = NULL,forcesym = F, sample_size = NULL,replace = F,
                          nb_iter = 1,
                          l_cm = NULL){
  l_flat_dt <- NULL
  nb_cpl_dt <- NULL
  dst_cpl_dt <- NULL
  out <- list()
  
  for(idx in 1:length(l_bait_go_tile)){
    bait.nm <- names(l_bait_go_tile)[idx]
    anchor.nm <- names(l_anchor_go_tile)[idx]
    if (! all(length(l_bait_go_tile[[idx]]) | length(l_anchor_go_tile[[idx]]))) stop("One of bait or anchor is null. You can't perform an APA with one of them NULL, check your input data and filters.", call. = FALSE)
    
    pw_dt <- pairwiseInteractionDT(anc_go = l_anchor_go_tile[[anchor.nm]],
                                   bait_go =  l_bait_go_tile[[bait.nm]],
                                   hic_res = hic_res,
                                   apa_size =  apa_size,
                                   d_min = d_min,
                                   d_max = d_max,
                                   forceSym = forcesym,
                                   cons_nm = cons_nm,
                                   sym = orientation)
    if (!is.null(sample_size)) {
      nr_pw_dt <- nrow(pw_dt)
      smpl_s <- nb_iter * floor(ifelse(sample_size<1, nr_pw_dt*sample_size,sample_size)) # Allow for percentage or absolute sample value, but use floor to avoid non integer values for sampling (error when doing groups)
      if (nrow(pw_dt) > smpl_s & !replace){
        pw_dt <- pw_dt[sample(.N,smpl_s)]
        pw_dt$key <- 1:nrow(pw_dt)
      } else if (nrow(pw_dt) < smpl_s & !replace){
        stop("The number of bait/anchor couples is inferior to sample size and can't be sampled with replace = FALSE. Try a different sample value or allow for replacing.", call. = FALSE)
      } else {
        pw_dt <- pw_dt[sample(.N,smpl_s,replace = T)] 
        pw_dt$key <- 1:nrow(pw_dt)
      }
      pw_dt$group <- sample(rep(1:nb_iter,smpl_s/nb_iter))
    } else {
      pw_dt$group <- 1
    }
    
    if( !is.null(control) ) {
      anc_ctl <- GenomicRanges::shift(l_anchor_go_tile[[anchor.nm]], control)
      bait_ctl <- GenomicRanges::shift(l_bait_go_tile[[anchor.nm]], control)
      # Take care of oob ranges if any
      a_oob <- GenomicRanges:::get_out_of_bound_index(anc_ctl)
      b_oob <- GenomicRanges:::get_out_of_bound_index(bait_ctl)
      if( !is.null(a_oob)) anc_ctl <- anc_ctl[-a_oob]
      if( !is.null(b_oob)) bait_ctl <- bait_ctl[-b_oob]
      
      anc_ctl$chr_bin <- anc_ctl$chr_bin + floor(control/res)
      bait_ctl$chr_bin <- bait_ctl$chr_bin + floor(control/res)
      pw_ctl_dt <- pairwiseInteractionDT(anc_go = anc_ctl,
                                         bait_go =  bait_ctl,
                                         hic_res = hic_res,
                                         apa_size =  apa_size,
                                         d_min = d_min,
                                         d_max = d_max,
                                         forceSym = forcesym,
                                         cons_nm = cons_nm,
                                         sym = orientation)
      if (!is.null(sample_size)) pw_ctl_dt <- pw_ctl_dt[sample(.N,sample_size)]
      
    }
    nb_cpl_dt <- rbind.data.frame(nb_cpl_dt, cbind.data.frame(name = paste0(bait.nm,"\nVS\n",anchor.nm), value = nrow(pw_dt)))
    dst_cpl_dt <- rbind.data.frame(dst_cpl_dt, cbind.data.frame(name = paste0(bait.nm,"\nVS\n",anchor.nm), value = (pw_dt$dst_bin)*hic_res))
    for(cm.nm in names(l_cm)){
      if (nrow(pw_dt) > 1){
        l_flat_dt[[paste0(bait.nm,"\nVS\n",anchor.nm)]][[cm.nm]] <- dumpFlatInteractionsMatrices(pairwise_dt = pw_dt, cm=l_cm[[cm.nm]], intra = T, n = apa_size)
        if (!is.null(control)) {
          l_flat_dt[[paste0(bait.nm,"\nVS\n",anchor.nm,"_control")]][[cm.nm]] <- dumpFlatInteractionsMatrices(pairwise_dt = pw_ctl_dt, cm=l_cm[[cm.nm]], intra = T, n = apa_size)
        }
      }
    }
  }
  return(list(lflat = l_flat_dt, nb_cpl = nb_cpl_dt, dst_cpl = dst_cpl_dt))
}


#' Create pairwise bin/bin interactions data frame
#' to use to dump HiC contacts
#'
#' @param anc_go Tile genome of anchor
#' @param bait_go Tile genome of Bait
#' @param apa_size Size of APA. Default 41
#' @param constraint_go Constraint GO. TAD for instance or reduce(tile_genome) for all interactions
#' @param hic_res
#' @param d_min distance minimum between bait & anc in bin
#' @param d_max distance max between bait & anc in bin
#'
#' @return
#' @export
#'
#' @examples
#' pw_dt <- pairwiseInteractionDT(anc_go = anchor_go_tile,bait_go =  bait_go_tile,hic_res =res ,apa_size =  41,constraint_go = reduce(genome_t_gr))
pairwiseInteractionDT <- function(anc_go=NULL,
                                  bait_go=NULL,
                                  apa_size=41,
                                  d_min=4,
                                  d_max=100,
                                  forceSym = F,
                                  # tile_go=NULL,
                                  cons_nm=NULL,
                                  sym = F,
                                  hic_res=1000){
  if (!any(cons_nm == "inter", cons_nm == "intra")) {
    # Load constraint if any (like a TAD)
    constraint_go <- loadConstraint(cons_nm)
  } else {
    # chromosomes collection
    chrs <- unique(c(as.vector(seqnames(bait_go)),as.vector(seqnames(anc_go))))
    Go_cons <- NULL
    for (chr in chrs) {
      min_anc <- min(start(subset(anc_go, seqnames == chr)))
      max_anc <- max(end(subset(anc_go, seqnames == chr)))
      min_bait <- min(start(subset(bait_go, seqnames == chr)))
      max_bait <- max(end(subset(bait_go, seqnames == chr)))
      Go_cons <- rbind(Go_cons, data.frame(start = min(min_anc,min_bait), end = max(max_anc,max_bait) ,seqnames = chr))
    }
    
    constraint_go <- GenomicRanges::GRanges(Go_cons)
  }
  names(constraint_go) <- NULL
  #Bait feature
  bait_r_go <- IRanges::resize(bait_go, 1, fix="center") #resize the feature to do the overlap on domains without multi hits
  bait_r_go$cb_b <- ceiling((GenomicRanges::start(bait_r_go)) / hic_res )
  bait_r_go$l_idx_b <- 1:length(bait_r_go)
  mol_dt <- data.table::data.table( data.frame(IRanges::mergeByOverlaps(constraint_go, bait_r_go) ))
  mol_dt$filt <- paste0(mol_dt$constraint_go.seqnames, "_", mol_dt$constraint_go.start)
  mol_b_dt <- mol_dt[ , c( "bait_r_go.name", "chr_bin", "filt", "bait_r_go.seqnames", "bait_r_go.strand" ) ]
  colnames(mol_b_dt) <- c("bait_name", "bait_chrom_bin", "constraint_id", "bait_chrom", "bait_strand")
  data.table::setkey(mol_b_dt, constraint_id)
  #Anchoring feature
  anc_r_go <- IRanges::resize(anc_go, 1, fix="center") #resize the feature to do the overlap on domains without multi hits
  anc_r_go$cb_b <- ceiling((GenomicRanges::start(anc_r_go)) / hic_res )
  anc_r_go$l_idx_b <- 1:length(anc_r_go)
  mol_dt <- data.table::data.table( data.frame(IRanges::mergeByOverlaps(constraint_go, anc_r_go) ))
  mol_dt$filt <- paste0(mol_dt$constraint_go.seqnames, "_", mol_dt$constraint_go.start)
  mol_a_dt <- mol_dt[ , c( "anc_r_go.name", "chr_bin", "filt", "anc_r_go.seqnames", "anc_r_go.strand" ) ]
  colnames(mol_a_dt) <- c("anc_name", "anc_chrom_bin", "constraint_id", "anc_chrom", "anc_strand")
  data.table::setkey(mol_a_dt, constraint_id)
  # Merge
  inter_dt <- mol_b_dt[mol_a_dt, allow.cartesian=TRUE,nomatch=0]
  
  if(nrow(inter_dt) == 0) return(NULL)
  
  if (forceSym) { # Useful when doing leftArms vs Right Arms and one wants to have both orientation in the APA (for a given couple)
    sym_dt <- data.table::copy(inter_dt)
    bait_v <- grep("bait",colnames(sym_dt),value = T)
    anc_v <- sub("bait","anc",bait_v)
    data.table::setnames(sym_dt,bait_v,sub("bait","new",bait_v))
    data.table::setnames(sym_dt,anc_v,bait_v)
    data.table::setnames(sym_dt,sub("bait","new",bait_v),sub("bait","anc",bait_v))
    inter_dt <- rbind(sym_dt,inter_dt)
  }
  
  # Mark the one to revert
  inter_dt$rev <- ifelse(inter_dt$anc_chrom_bin < inter_dt$bait_chrom_bin, T, F)
  # Rmove reverse if bait == anchor (duplicated, we don't need it)
  if(sym == "5p3p"){
    inter_dt <- inter_dt[inter_dt$rev==F,]
  } else if (sym == "3p5p") {
    inter_dt <- inter_dt[inter_dt$rev==T,]
  } # else we take both
  # Add distance between int and anc in bin
  inter_dt$dst_bin <- abs(inter_dt$bait_chrom_bin - inter_dt$anc_chrom_bin)
  # Add a key
  inter_dt$key <- 1:nrow(inter_dt)
  # Filter distnace min and max
  inter_dt <- subset(inter_dt, dst_bin >= d_min & dst_bin <= d_max)
  # RETURN DATA.TABLE
  inter_dt
}





#' Dump all interaction matrices
#' Stored as flat data table
#' One row stand for a complete nxn matrix slice
#' centered on a given bait/anchor interaction
#'
#' @param pairwise_dt Output of pairwiseInteractionDT.
#' @param intra Is the interaction intrachrom only ? Default = TRUE
#' @param n Size of the matrices. Default = 41
#' @param cm
#'
#' @return
#' @export
#'
#' @examples
#' itr_mat <- dumpInteractionMatrix(pairwise_dt = pw_filt_dt, intra = T, n = n)

dumpFlatInteractionsMatrices <- function (pairwise_dt = NULL,cm = NULL, intra = T, n = 41) {
  if (intra) {
    stack_dt <- NULL
    # We can directly check either bait or anc seqnames as they are the same in intra
    for (chrom in unique(pairwise_dt$bait_chrom)){# For each chromosomes (might be parallelized later)
      pw_chr_dt <- subset(pairwise_dt, bait_chrom == chrom)
      rv_l <- list()
      for (reverse in unique(pw_chr_dt$rev)){ # For couples that has to be reversed or not
        print(paste0(chrom ," ",reverse))
        pw_chr_rev_dt <- subset(pw_chr_dt, rev == reverse)
        rv_l[[paste0("to rev ",reverse)]] <- nrow(pw_chr_rev_dt)
        
        # Subset by chrom
        chr_cm <- cm[[chrom]]
        # Launch fast algorithm
        bb <- pw_chr_rev_dt$bait_chrom_bin
        ba <- pw_chr_rev_dt$anc_chrom_bin
        # Make all bin by repeat
        ba_r <- rep(ba, each=n*n)
        bb_r <- rep(bb,each=n*n)
        # Make mask for extending bins
        mask <- -((n-1)/2):((n-1)/2)
        # apply it to Anch and Bait
        a_mask <- rep(rep(mask, each=n),nrow(pw_chr_rev_dt))
        b_mask <- rep(rep(mask, times=n),nrow(pw_chr_rev_dt))
        # Add to ba & bb
        ba_mask <- ba_r + a_mask
        bb_mask <- bb_r + b_mask
        
        idx <- ceiling(n*n/2)
        j <-  (idx-1)%%n + 1
        i <-  (idx-1)%/%n + 1
        
        bin_i <- i - (n+1)/2 + 20 # The +20 comes from na extend
        bin_j <- j - (n+1)/2 + 20 # The +20 comes from na extend
        itr_vec <- chr_cm[((ba_mask + bin_i)-1) * dim(chr_cm)[1] + bb_mask + bin_j]
        # print(paste0("Rev = ",reverse))
        # obj_s <- object.size(itr_vec)
        # print(paste0("Size vector = ",round(obj_s/1e6,3), " Mb"))
        if (reverse) {
          itr_vec <- rev(itr_vec)
          nkey <- rev(pw_chr_rev_dt$key)
        } else {
          nkey <- pw_chr_rev_dt$key
        }
        itr_mat <- matrix(itr_vec, nrow = nrow(pw_chr_rev_dt),ncol=n*n,byrow = T,dimnames = list(paste0("mat",1:nrow(pw_chr_rev_dt)),paste0("pix",1:(n*n))))
        itr_matrices_dt <- data.table::data.table(dtkey=nkey, itr_mat)
        stack_dt <- rbind(stack_dt, itr_matrices_dt)
      }
      print(paste0(" Chrom : ", chrom , " Nb couple : ", nrow(pw_chr_dt)))#, " - ", paste0(names(rv_l)[1], " : ",rv_l[[1]], " / ", paste0(names(rv_l)[2], " : ",rv_l[[2]]))))
      
      gc()
    }
  }
  # # Join tables back to have access to bin info
  data.table::setkey(stack_dt, dtkey)
  data.table::setkey(pairwise_dt, key)
  flat_all_chaos <- pairwise_dt[stack_dt]
  flat_all_chaos
}



### POSTPROCESS EXTRACTED SUBMATRICES ----

#' Title
#'
#' @param l_flat_dt
#'
#' @return
#' @export
#'
#' @examples
aggreg_lflat <- function(l_flat_dt){
  
  for(i in names(l_flat_dt)){
    l_flat_dt[[i]] <- rbindlist(l_flat_dt[[i]], idcol = "cm")
  }
  
  flat_dt <- rbindlist(l_flat_dt, idcol = "couple_name")
  flat_dt <- tidyr::separate(flat_dt, "couple_name", into = c("bait","anchor"),sep = "\nVS\n")
  setDT(flat_dt)
}


#' Title
#'
#' @param flat_dt
#' @param offset
#'
#' @return
#' @export
#'
#' @examples
remove_diag <- function(flat_dt, offset,apa_size){
  dst_rm_offdiag <- unique(subset(flat_dt, dst_bin < (apa_size+offset-1))$dst_bin)
  
  for(i in dst_rm_offdiag){
    # i_rm <-  i - offset
    i_rm <-  abs(i - apa_size - (offset - 1))
    msk <- mask[(apa_size-i_rm):apa_size,1:(i_rm+1)]
    idx <- msk[lower.tri(msk, T)]
    if(length(idx>0))
      flat_dt[dst_bin == i, (paste0("pix", idx)) := NA]
  }
  flat_dt
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
quant_mat <- function(x) {
  todel <- which(x==0)
  if(length(todel)){
    ok <- ntile(x[-todel], quantiles)
    # x[todel] <- NA
    x[-todel] <- ok #+ length(which(is.na(x)))
    x[which(is.na(x))] <- 0
  } else {
    x <-  ntile(x,quantiles)
  }
  unlist(x)
}

#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
quant_mat_with_na <- function(x) { 
  todel <- which(x==0)
  if(length(todel)){
    ok <- ntile(x[-todel], quantiles)
    x[todel] <- NA
    x[-todel] <- ok + length(which(is.na(x)))
    x[which(is.na(x))] <- 0
  } else {
    x <-  ntile(x,quantiles)
  }
  unlist(x)
}


#' Title
#'
#' @param working_dt
#' @param norm
#' @param rangeMinMax
#' @param metrics_names
#'
#' @return
#' @export
#'
#' @examples
make_avg_flat <- function(working_dt, norm, rangeMinMax, metrics_names){
  pix_met_col <- colnames(working_dt)[names(working_dt) %like% paste0(c("pix",metrics_names),collapse = "|")] # Transform pixel and metrics column
  pix_col <- colnames(working_dt)[names(working_dt) %like% "pix"] # Transform pixel and metrics column
  
  if (norm == "quant") {
    cat("Quantile normalization \n")
    flat_qtmp_dt <- t(apply(working_dt[,.SD, .SDcols = pix_col], 1, quant_mat))
    colnames(flat_qtmp_dt) <- paste0("pix",1:ncol(flat_qtmp_dt))
    tl <- names(working_dt)[names(working_dt) %like% "pix"]
    working_dt[,(tl):=NULL]
    working_dt <- cbind(working_dt,flat_qtmp_dt)
    setApaMetrics(working_dt, coord_l)
    
  } else if(norm == "log2") {
    cat("log2 transformation\n")
    avg_flat_dt <- working_dt[, lapply(.SD, function(x){sign(x)*log2(abs(x))}), .SDcols = pix_met_col][, lapply(.SD, function(x) {m = mean(x,na.rm=T);m[!is.finite(m)] <- 0;m}),
                                                                                                       .SDcols = names(working_dt) %like% "pix", by = .(bait, anchor, cm)]
    
  } else if (norm == "log10"){
    cat("log10 transformation\n")
    avg_flat_dt <- working_dt[, lapply(.SD, function(x){sign(x)*log10(abs(x))}), .SDcols = pix_met_col][, lapply(.SD, function(x) {m = mean(x,na.rm=T);m[!is.finite(m)] <- 0;m}),
                                                                                                        .SDcols = names(working_dt) %like% "pix", by = .(bait, anchor, cm)]
    
  } else if(norm == "log") {
    cat("log transformation\n")
    avg_flat_dt <- working_dt[, lapply(.SD, function(x){sign(x)*log1(abs(x))}), .SDcols = pix_met_col][, lapply(.SD, function(x) {m = mean(x,na.rm=T);m[!is.finite(m)] <- 0;m}),
                                                                                                       .SDcols = names(working_dt) %like% "pix", by = .(bait, anchor, cm)]
    
  } else if(norm == "alog") {
    cat("alog transformation\n")
    avg_flat_dt <- working_dt[, lapply(.SD, function(x){ifelse(x < -1 | x > 1, sign(x)*log2(abs(x)), x)}), .SDcols = pix_met_col][, lapply(.SD, function(x) {m = mean(x,na.rm=T);m[!is.finite(m)] <- 0;m}),
                                                                                                                                  .SDcols = names(working_dt) %like% "pix", by = .(bait, anchor, cm)]
    
  } else if(!is.null(rangeMinMax)){
    cat("RangeMinMax transformation\n")
    avg_flat_dt <- working_dt[, lapply(.SD, function(x){rangeMinMax(x)}), .SDcols = pix_met_col][, lapply(.SD, function(x) {m = mean(x,na.rm=T);m[!is.finite(m)] <- 0;m}),
                                                                                                 .SDcols = names(working_dt) %like% "pix", by = .(bait, anchor, cm)]
    
  } else {
    cat("No transformation\n")
    avg_flat_dt <- working_dt[, lapply(.SD, function(x) {m = mean(x,na.rm=T);m[!is.finite(m)] <- 0;m}),
                              .SDcols = names(working_dt) %like% "pix", by = .(bait, anchor, cm)]
  }
  
  list(w_dt = working_dt, avg_dt = avg_flat_dt)
}

#' Title
#'
#' @param avg_flat_dt
#' @param apa_size
#'
#' @return
#' @export
#'
#' @examples
convert_flat_vector_to_matrices <- function(avg_flat_dt, apa_size = n){
  out <- list()
  l_m <- list()
  min = Inf
  max = -Inf
  for (i in 1:nrow(avg_flat_dt)){
    bait.nm <- avg_flat_dt[i, bait]
    anc.nm <- avg_flat_dt[i, anchor]
    cm.nm <- avg_flat_dt[i, cm]
    m <- unlist(avg_flat_dt[bait == bait.nm & anchor == anc.nm & cm == cm.nm, .SD, .SDcols = !c("bait","anchor", "cm")])
    m <- matrix(m, apa_size, apa_size, byrow = F)
    m2t <- matrix2tibble(m, apa_size, apa_size)
    names(m2t) <- c(paste0(bait.nm,".Bait") , paste0(anc.nm,".Anchor"), "values")  # Bait  Anchor counts
    min <- min(min(m2t$values), min)
    max <- max(max(m2t$values), max)
    l_m[[cm.nm]][[bait.nm]][[anc.nm]] <- m2t
  }
  out[["min"]] <- min
  out[["max"]] <- max
  out[["l_m"]] <- l_m
  out
}





### APA GRAPHICS ----
#' Title
#'
#' @param l_m
#' @param norm
#' @param params
#' @param n
#' @param width
#' @param height
#' @param output
#' @param min_l
#' @param max_l
#'
#' @return
#' @export
#'
#' @examples
gg_apa_panels <- function(l_m, norm, params,n ,width ,height, output, min_l, max_l, palette="lajolla", minmax = NULL){
  sl_gg <- list()
  # if (as.numeric(gsub("\\.","",R.version$minor)) > 10) direction <-  1 else direction <- -1
  direction <- -1
  rmax <- 0
  if (!is.null(minmax)) {
    rmax <- minmax
  } else {
    for(cm.nm in names(l_m)){
      for(bait.nm in names(l_m[[cm.nm]])){
        for (anc.nm in names(l_m[[cm.nm]][[bait.nm]])){
          val <- l_m[[cm.nm]][[bait.nm]][[anc.nm]]$values
          rmax <- max(rmax,abs(max(val) - min(val)))
        }
      }
    }
  }
  
  for(cm.nm in names(l_m)){
    for(bait.nm in names(l_m[[cm.nm]])){
      for (anc.nm in names(l_m[[cm.nm]][[bait.nm]])){
        gg <- ggplot2::ggplot(l_m[[cm.nm]][[bait.nm]][[anc.nm]], ggplot2::aes_string(x=paste0(anc.nm,".Anchor"), y=paste0(bait.nm,".Bait"), fill="values")) +
          ggplot2::geom_tile() +  ggtitle(paste0("Condition : ",cm.nm,"\n ",
                                                 "HiC metrics : ",norm)) +
          ggplot2::theme(axis.text.y = element_blank(),axis.text.x = element_blank()) +
          scale_y_discrete(breaks=(n+1)/2) +
          scale_x_discrete(breaks=(n+1)/2)
        
        ggs <- gg + ggplot2::scale_fill_gradientn(colors = scico::scico(6, palette = palette,direction = direction), 
                                                  breaks = round(seq(min_l, max_l, length.out = 5),digits = 3),
                                                  limits = c(min_l, max_l))
        ggns <- gg + ggplot2::scale_fill_gradientn(colors = scico::scico(100, palette = palette,direction = direction))
        
        sl_gg[["scaled"]][[paste0("APA","_", cm.nm, "_", bait.nm,"_",anc.nm,"_",norm)]] <- ggs
        sl_gg[["unscaled"]][[paste0("APA","_", cm.nm, "_", bait.nm,"_",anc.nm,"_",norm)]] <- ggns
        
        
        # Cut extreme values in +/- rmax/2
        
        nmat <- l_m[[cm.nm]][[bait.nm]][[anc.nm]]
        mmv <- mean(c(max(nmat$values),min(nmat$values)))
        upperCut <- nmat$values >= (mmv  + (rmax/2))
        nmat$values[upperCut] <- mmv  + (rmax/2)
        lowerCut <- nmat$values < (mmv  - (rmax/2))
        nmat$values[lowerCut] <- mmv  - (rmax/2)
        
        gg <- ggplot2::ggplot(nmat, ggplot2::aes_string(x=paste0(anc.nm,".Anchor"), y=paste0(bait.nm,".Bait"), fill="values")) +
          ggplot2::geom_tile() +  ggtitle(paste0("Condition : ",cm.nm,"\n ",
                                                 "HiC metrics : ",norm)) +
          ggplot2::theme(axis.text.y = element_blank(),axis.text.x = element_blank()) +
          scale_y_discrete(breaks=(n+1)/2) +
          scale_x_discrete(breaks=(n+1)/2)
        
        ggmms <- gg + ggplot2::scale_fill_gradientn(colors = scico::scico(100, palette = palette,direction = direction), breaks = round(seq((mmv-rmax/2), (mmv+rmax/2), length.out = 4), 3),
                                                    limits = c((mmv-rmax/2), ceiling((mmv+rmax/2))))
        
        
        sl_gg[["minmax_scaled"]][[paste0("APA","_", cm.nm, "_", bait.nm,"_",anc.nm,"_",norm)]] <- ggmms
      }
    }
  }
  nb_rc <- ceiling(sqrt(max(nrow(avg_flat_dt), length(l_m))))
  ggs <- cowplot::plot_grid(plotlist = sl_gg$scaled,
                            nrow = nb_rc,
                            ncol = nb_rc)
  ggns <- cowplot::plot_grid(plotlist = sl_gg$unscaled,
                             nrow = nb_rc,
                             ncol = nb_rc)
  ggmms <- cowplot::plot_grid(plotlist = sl_gg$minmax_scaled,
                              nrow = nb_rc,
                              ncol = nb_rc)
  fig.size <- max(width, height) * 1.5
  i <- "APA_panel"
  if ("APA_panel" %in% params$plot_lst){
    if("svg" %in% params$fig_type){
      ggsave(paste0(output,".APA_panel.",norm,".Scaled.svg"), ggs,width = fig.size*nb_rc, height = fig.size*nb_rc, units = "px", dpi = 200)
      ggsave(paste0(output,".APA_panel.",norm,".Unscaled.svg"), ggns,width = fig.size*nb_rc, height = fig.size*nb_rc, units = "px", dpi = 200)
      ggsave(paste0(output,".APA_panel.",norm,".MinMaxScaled.svg"), ggmms,width = fig.size*nb_rc, height = fig.size*nb_rc, units = "px", dpi = 200)
      
    }
    if("png" %in% params$fig_type){
      ggsave(paste0(output,".APA_panel.",norm,".Scaled.png"), ggs,width = fig.size*nb_rc, height = fig.size*nb_rc, units = "px", dpi = 200)
      ggsave(paste0(output,".APA_panel.",norm,".Unscaled.png"), ggns,width = fig.size*nb_rc, height = fig.size*nb_rc, units = "px", dpi = 200)
      ggsave(paste0(output,".APA_panel.",norm,".MinMaxScaled.png"), ggmms,width = fig.size*nb_rc, height = fig.size*nb_rc, units = "px", dpi = 200)
    }
    if("pptx" %in% params$fig_type){ # See heatmap panel for Pradere where we remove the raster issue when doing pptx with multiple plots
      pptx <- read_pptx() %>%
        add_slide() %>%
        ph_with(value = dml(ggobj = ggs), location = ph_location("body", left = 1, top = 1, width = 9, height = 4.95))%>%
        print(target = file.path(paste0(output,".APA_panel.",norm,".Scaled.pptx")))
      pptx <- read_pptx() %>%
        add_slide() %>%
        ph_with(value = dml(ggobj = ggns), location = ph_location("body", left = 1, top = 1, width = 9, height = 4.95))%>%
        print(target = file.path(paste0(output,".APA_panel.",norm,".Unscaled.pptx")))
      pptx <- read_pptx() %>%
        add_slide() %>%
        ph_with(value = dml(ggobj = ggmms), location = ph_location("body", left = 1, top = 1, width = 9, height = 4.95))%>%
        print(target = file.path(paste0(output,".APA_panel.",norm,".MinMaxScaled.pptx")))
    }
  }
  invisible(list(scaled = ggs, unscaled = ggns, minmax_scaled = ggmms))
}


### UTILS ----
h52dense <- function(h5m, chr = NULL, naOffset = 20) {
  
  if (is.null(chr)) stop("Parameter chr is null. You have to set a value to extract a given submatrix.", call. = FALSE)
  
  m <- as.matrix(as.matrix(h5m[[as.character(chr)]]))
  max_bin <- dim(m)[1] - naOffset # Remove NA's extend use in APA
  m <- m[(naOffset+1):max_bin, (naOffset+1):max_bin] # Extract correct submatrix
  m
}

# Obtain the genome tiled indices for this GRanges at a given resolution
getTiledIndices <- function(gr, fix = "center", resolution = NULL) {
  checkmate::assert_class(gr, "GRanges")
  checkmate::assert_number(resolution)
  checkmate::assert_int(resolution)
  checkmate::assert_character(fix,max.len = 1)
  ceiling(GenomicRanges::start(GenomicRanges::resize(gr,1,fix))/resolution)
}


get_apa_metrics_coordinate <- function(apa_size, mask, ctloffset=6, square_size) {
  
  # Check offset consistency
  if(floor(apa_size/2) - (ctloffset + 1) < 0){
    ctloffset <- floor(apa_size/2) - 1
  }
  
  # Coordinate of central pixel
  CP <- (apa_size+1)/2
  CP <- mask[CP, CP]
  
  # Coordinate of 3x3 square centered on central pixel
  m_size <- (apa_size-1)/2 # Middle pixel value - 1
  i_CS <- m_size:(m_size+2)
  j_CS <- m_size:(m_size+2)
  CS <- as.vector(mask[i_CS, j_CS])
  
  # Coordinate of Upper left 3x3 square from 3x3 square at center
  i_ULS <- i_CS-ctloffset
  j_ULS <- j_CS-ctloffset
  ULS <- as.vector(mask[i_ULS, j_ULS])
  
  # Coordinate of Upper Right square starting from central pixel (Compartment)
  i_URS <- 1:floor(apa_size/2)
  j_URS <- (ceiling(apa_size/2)+1):apa_size
  URS <- as.vector(mask[i_URS, j_URS])
  
  # Coordinate of Bottom Right 3x3 square from 3x3 square at center (Compartment)
  i_BRS <- i_CS+ctloffset
  j_BRS <- j_CS+ctloffset
  BRS <- as.vector(mask[i_BRS, j_BRS])
  
  # Coordinate of  Bottom Left square starting from central pixel (TAD)
  i_BLS <- (ceiling(apa_size/2)+1):apa_size
  j_BLS <- 1:floor(apa_size/2)
  BLS <- as.vector(mask[i_BLS, j_BLS])
  
  # Loop extrusion, left side from central pixel
  i_LEX <- ceiling(apa_size/2)
  j_LEX <- (ceiling(apa_size/2)-8):(ceiling(apa_size/2)-1)
  LEX <- as.vector(mask[i_LEX, j_LEX])
  
  
  # Outer perimeters squares
  lv2_coord <- floor(apa_size/2) : (floor(apa_size/2) + 2)
  lv3_coord <- (floor(apa_size/2)-1) : (floor(apa_size/2) + 3)
  lv4_coord <- (floor(apa_size/2)-2) : (floor(apa_size/2) + 4)
  dli_9x9 <- (floor(apa_size/2) - 3) : (floor(apa_size/2) + 5)
  
  LV2 <- mask[lv2_coord, lv2_coord][-5]
  LV3 <- mask[lv3_coord, lv3_coord][-c(7:9, 12:14, 17:19)]
  LV4 <- mask[lv4_coord, lv4_coord][-c(9:13, 16:20, 23:27, 30:34,37:41)]
  DLI <- mask[dli_9x9,dli_9x9][-41]
  
  return(list(CP = CP, CS = CS, BLS = BLS, BRS = BRS, ULS = ULS, URS = URS, LEX = LEX,LV4 = LV4, LV3 = LV3, LV2 = LV2, DLI = DLI))
}

get_apa_dt <- function(apa_dt, coord_l, save_dt, norm, params, output) {
  
  apa_dt[,c_name:=ifelse(bait == anchor, paste0(bait,"_homot"), paste0(bait,"_vs_",anchor))]
  
  apa_dt[ , ctrl := (ULS + BRS)/2]
  
  apa_dt[,CP.norm :=(CP+1)/(ctrl+1)]
  apa_dt[,CS.norm :=(CS+1)/(ctrl+1)]
  apa_dt[,BLS.norm :=(BLS+1)/(ctrl+1)]
  apa_dt[,ULS.norm :=(ULS+1)/(ctrl+1)]
  apa_dt[,BRS.norm :=(BRS+1)/(ctrl+1)]
  apa_dt[,LEX.norm :=(LEX+1)/(CP+1)]
  apa_dt[,LV2.norm :=(LV2+1)/(ctrl+1)]
  apa_dt[,LV3.norm :=(LV3+1)/(ctrl+1)]
  apa_dt[,LV4.norm :=(LV4+1)/(ctrl+1)]
  apa_dt[,DLI.norm :=(DLI+1)/(ctrl+1)]
  
  
  flat_dt.tmp <- apa_dt[, .(cm, c_name,
                            CP, CS, BLS, ULS,  BRS, LEX,
                            CP.norm, CS.norm, BLS.norm, ULS.norm, BRS.norm,LEX.norm,
                            LV2, LV3, LV4,DLI,
                            LV2.norm, LV3.norm, LV4.norm,DLI.norm,  ctrl)]
  dt <- data.table::melt(flat_dt.tmp, id.vars = c("cm", "c_name"))
  if ( save_dt == T ) saveRDS(dt, paste0(output,"_apa_metric_",norm,"_dt.rds"))
  dt
}

#' Set APA metrics
#'
#' Set APA metrics and set it as column of flat data table.
#' It will modify directly the `apa_dt` data.table passed as reference.
#'
#' @param apa_dt A flat `data.table`.
#' @param coord_l A `list(n)`. List of metrics coordinates.
#'
#' @return A `data.table` with metrics APA.
#'
#' @export
setApaMetrics <- function(apa_dt, coord_l) {
  for (i in names(coord_l)){
    tmp <- apa_dt[, rowMeans(.SD, na.rm = T), .SDcols = paste0("pix",coord_l[[i]])]
    apa_dt[,(i):=tmp]
  }
}

############***################ #
############################### #



############################ #
#
# make color scale (tricolor with center value)
#
# return : list(breaks.val = breaks.values, breaks.col = breaks.colors)
############################ #




# NxN matrix to ggplotable heatmap data frame

#' Title
#'
#' @param m
#' @param n_r
#' @param n_c
#'
#' @return
#' @export
#'
#' @examples
matrix2tibble <- function(m, n_r=41, n_c = 41) {
  colnames(m) <- paste0(seq(1,n_c))
  rownames(m) <- paste0(seq(1,n_r))
  m %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Bait") %>%
    tidyr::pivot_longer(-c(Bait), names_to = "Anchor", values_to = "counts") %>%
    dplyr::mutate(Anchor= forcats::fct_relevel(Anchor,colnames(m))) %>%
    dplyr::mutate(Bait= forcats::fct_relevel(Bait,rev(rownames(m))))
}





# reduce matrix size, using a summarizing function (default, mean)
redim_matrix <- function(
    mat,
    target_height = 100,
    target_width = 100,
    summary_func = function(x) mean(x, na.rm = TRUE),
    output_type = 0.0, #vapply style
    n_core = 1 # parallel processing
) {
  
  if(target_height > nrow(mat) | target_width > ncol(mat)) {
    stop("Input matrix must be bigger than target width and height.")
  }
  
  seq_height <- round(seq(1, nrow(mat), length.out = target_height + 1))
  seq_width  <- round(seq(1, ncol(mat), length.out = target_width  + 1))
  
  # complicated way to write a double for loop
  do.call(rbind, parallel::mclapply(seq_len(target_height), function(i) { # i is row
    vapply(seq_len(target_width), function(j) { # j is column
      summary_func(
        mat[
          seq(seq_height[i], seq_height[i + 1]),
          seq(seq_width[j] , seq_width[j + 1] )
        ]
      )
    }, output_type)
  }, mc.cores = n_core))
}




#' PLot position of APA control area of interest
#'
#' @param n `integer(1)` Size of APA
#' @param coordinate `list(n)` List of coordinates
#'
#' @return A `pheatmap` plot
#' @export
plotApaCtlPosition <- function(n, coordinate){
  m <- matrix(rep(0,n*n),n,n)
  for (i in seq_along(coordinate)) m[coordinate[[i]]] <- i
  pheatmap::pheatmap(m,cluster_rows = F,cluster_cols = F)
}


flatToApaDt <- function(working_dt,zero2NA = F) {
  tosel <- colnames(working_dt)[colnames(working_dt)%like% "pix"]
  if (zero2NA) working_dt[,(tosel):=lapply(.SD, function(x) ifelse(x==0,NA,x)),.SDcols = tosel]
  apa_dt <- working_dt[, lapply(.SD, function(x) {m = mean(x,na.rm=T);m[!is.finite(m)] <- 0;m}),
                       .SDcols = names(working_dt) %like% "pix", by = .(bait, anchor, cm)]
  apa_dt
}


quantilize_flat_dt <- function(working_dt) {
  cat("Quantile normalization \n")
  pix_col <- colnames(working_dt)[names(working_dt) %like% "pix"] # Transform pixel and metrics column
  flat_qtmp_dt <- t(apply(working_dt[,.SD, .SDcols = pix_col], 1, quant_mat))
  colnames(flat_qtmp_dt) <- paste0("pix",1:ncol(flat_qtmp_dt))
  tl <- names(working_dt)[names(working_dt) %like% "pix"]
  working_dt[,(tl):=NULL]
  working_dt[, (pix_col) := as.data.table(flat_qtmp_dt)]
  setApaMetrics(working_dt, coord_l)
}


quantilize_na_flat_dt <- function(working_dt) {
  cat("Quantile normalization \n")
  pix_col <- colnames(working_dt)[names(working_dt) %like% "pix"] # Transform pixel and metrics column
  flat_qtmp_dt <- t(apply(working_dt[,.SD, .SDcols = pix_col], 1, quant_mat_with_na))
  colnames(flat_qtmp_dt) <- paste0("pix",1:ncol(flat_qtmp_dt))
  tl <- names(working_dt)[names(working_dt) %like% "pix"]
  working_dt[,(tl):=NULL]
  working_dt[, (pix_col) := as.data.table(flat_qtmp_dt)]
  setApaMetrics(working_dt, coord_l)
}



transform_flat_dt <- function(working_dt, norm, rangeMinMax, metrics_names){
  pix_met_col <- colnames(working_dt)[names(working_dt) %like% paste0(c("pix",metrics_names),collapse = "|")] # Transform pixel and metrics column
  if(norm == "log2") {
    cat("log2 transformation\n")
    for (j in pix_met_col) set(working_dt, j = j, value = sign(working_dt[[j]])*log2(abs(working_dt[[j]])))
    
  } else if (norm == "log10"){
    cat("log10 transformation\n")
    for (j in pix_met_col) set(working_dt, j = j, value = sign(working_dt[[j]])*log10(abs(working_dt[[j]])))
    
  } else if(norm == "log") {
    cat("log transformation\n")
    for (j in pix_met_col) set(working_dt, j = j, value = sign(working_dt[[j]])*log1(abs(working_dt[[j]])))
    
  } else if(norm == "alog") {
    cat("alog transformation\n")
    for (j in pix_met_col) set(working_dt, j = j, value = ifelse(working_dt[[j]] < -1 | working_dt[[j]] > 1, sign(working_dt[[j]])*log2(abs(working_dt[[j]])),working_dt[[j]]))
    
  } else if(!is.null(rangeMinMax)){
    cat("RangeMinMax transformation\n")
    for (j in pix_met_col) set(working_dt, j = j, value = rangeMinMax(working_dt[[j]]))
    
  } else {
    cat("No transformation\n")
  }
  working_dt
}


#TODO : GET OR SET K-th diagonal of a matrix
# Input : mat = matrix
#         k = the k-th diagonal (negative for lower diagonals)
# Output : k-th diagonal

getKDiag <- function(mat=NULL,k=0){ # k=0 the first diagonal
  for(i in k){
    print(as.vector(mat[col(mat) - row(mat) == i]))
  }
}

'getKDiag<-' <- function(mat,k,value)
{
  UseMethod('getKDiag<-',mat)
}

'getKDiag<-.matrix' <- function(mat,k,value)
{
  for(i in k){
    mat[col(mat) - row(mat) == i] <- value
  }
  return(mat)
}
#getKDiag(ok,c(1:2))
#getKDiag(ok,-1:1) <- NA

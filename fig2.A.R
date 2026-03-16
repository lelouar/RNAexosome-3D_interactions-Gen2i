# Author: Heurteau Alexandre
# Date: Fri Jul 26 14:27:32 2019
#####################################################################################--
#          DESCRIPTION  ----
#####################################################################################--

# THIS SCRIPT AIMS TO

#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-
source("lib/lrc_lib.R")
source("lib/paxt_next_database.R")

#####################################################################################-
#
#          RUN  ----
#
#####################################################################################-


MTR4.l <- list(MTR4_OLD="mtr4.hg19.gr",MTR4_NEW="mtr4N.hg19.gr")
MTR4.bwl <- list(MTR4_OLD=mtr4.hg19.bwp,MTR4_NEW=mtr4N.hg19.bwp)
Z1.l <- list(Z1="z1.hg19.gr",Z1shMTR4="z1.shMTR4.hg19.gr")
Z1.bwl <- list(Z1=z1.hg19.bwp,Z1shMTR4=z1.shMTR4.hg19.bwp)
Z8.l <- list(Z8_OLD="z8.hg19.gr",Z8_NEW="z8N.hg19.gr")
Z8.bwl <- list(Z8_OLD=z8.hg19.bwp,Z8_NEW=z8N.hg19.bwp)

RAD21.l <- list(RAD21_OLD="rad21.hg19.gr",RAD21_NEW="rad21.shC.hg19.gr")
RAD21.bwl <- list(RAD21shC=rad21.shC.hg19.bwp,RAD21shMTR4=rad21.shMTR4.hg19.bwp)

GNME <- "hg19"

gg <- NULL

params <- NULL
params$type <- "ChIP" # "NCRna" or "ChIP"
# FIGURE 1 C-D-E - WITH Z8N and MTR4N & SINGLE TSS GENES----------------------------------------------------------
if(params$type == "NCRna") {
  TSS_D_S.l <- list(TSS_SGL="tss.hg19.sglTSS.pm2kb.gr")
  
  TSS.chr  <-  names(TSS_D_S.l)[1]
  MTR4.chr  <-  names(MTR4.l)[2]
  Z8.chr  <-  names(Z8.l)[2]
  
  TSS.hg19.gr <- resize(genes.hg19.sglTSS.pm2kb.gr,1,"start")
  TES.hg19.gr <- resize(genes.hg19.sglTSS.pm2kb.gr,1,"end")
  
  prom2Kb.hg19.gr <- extendGR(TSS.hg19.gr,2000,0)
  
  promptsZ1.hg19.gr <- myOverlaps("ncRNAZ1","PROM",ncRNA_Z1kd.hg19.gr,prom2Kb.hg19.gr)$gr1xgr2
  promptsMTR4.hg19.gr <- myOverlaps("ncRNAMTR4","PROM",ncRNA_MTR4kd.hg19.gr,prom2Kb.hg19.gr)$gr1xgr2
  promptsZ8.hg19.gr <- myOverlaps("ncRNAZ8","PROM",ncRNA_Z8kd.hg19.gr,prom2Kb.hg19.gr)$gr1xgr2
  
  rad21RK.hg19.gr <- get(RAD21.l$RAD21_NEW)
  
  
  MTR4.cgr <- MTR4.l[[MTR4.chr]]
  MTR4.hg19.gr <- get(MTR4.cgr)
  MTR.nm <- toupper(gsub("(.*).hg19.gr","\\1",MTR4.cgr))
  Z8.cgr <- Z8.l[[Z8.chr]]
  Z8.nm <- toupper(gsub("(.*).hg19.gr","\\1",Z8.cgr))
  
  
  MTR4xZ1.hg19.gr <- subsetByOverlaps(z1.hg19.gr,get(MTR4.cgr),maxgap=500)
  MTR4xZ8.hg19.gr <- subsetByOverlaps(get(Z8.cgr),get(MTR4.cgr),maxgap=500)
  
  Col.l <- list(MTR4="#ff5733",Z1="#1898dd",Z8="#7bb84f")
  
  
  myWhole.l <- c("z1.hg19.gr",MTR4.cgr,Z8.cgr,"MTR4xZ1.hg19.gr","MTR4xZ8.hg19.gr","ctcf.hg19.gr","rad21RK.hg19.gr",
                 "promptsZ1.hg19.gr","promptsZ8.hg19.gr","promptsMTR4.hg19.gr",
                 "ncRNA_Z1kd.hg19.gr","ncRNA_Z8kd.hg19.gr","ncRNA_MTR4kd.hg19.gr",
                 "TSS.hg19.gr","enhCAGE.hg19.gr","h3k4me1.hg19.gr")
  
  
  xtVec <- c(T,T,T,T,
             T,T,T,
             F,F,F,
             F,F,F,
             T,F,T)
  
  xtend <- 500
  genome <- "hg19"
  celltype <- "hela"
  
  all.gr <- trim(GenomicRanges::reduce(do.call("c",sapply(1:length(myWhole.l),function(X){
    gr <- get(myWhole.l[X])
    mcols(gr) <- NULL
    if(xtVec[X]){print("REDUCE"); ret_res <- addSeqinfo(resize(gr,1,"center")+xtend,genome,celltype);strand(ret_res) <- "*"}
    else{ret_res <- addSeqinfo(gr,genome,celltype);strand(ret_res) <- "*"}
    return(ret_res)}))))
  
  apa_DF <- data.frame(IDX=1:length(all.gr),ALL=1)
  for(feat in myWhole.l){
    nm <- toupper(gsub(paste0("(.*).",genome,".gr"),"\\1",feat))
    apa_DF[,nm] <- 0
    gr <- get(feat)
    mcols(gr) <- NULL
    gr <- addSeqinfo(gr,genome,celltype)
    strand(gr) <- "*"
    IDX <- unique(findOverlaps(gr,all.gr)@to)
    apa_DF[IDX,nm] <- 1
  }
  
  genesBody.xt.hg19.gr <- addSeqinfo(genesBody.hg19.gr ,"hg19","hela")
  mcols(genesBody.xt.hg19.gr) <- NULL
  seqlevelsStyle(genesBody.xt.hg19.gr) <- "Ensembl"
  
  xtcomp.l <- c("genesBody.xt.hg19.gr")
  
  for(feat in xtcomp.l){
    nm <- toupper(gsub(paste0("(.*).xt.",genome,".gr"),"\\1",feat))
    apa_DF[,nm] <- 0
    gr <- get(feat)
    mcols(gr) <- NULL
    gr <- addSeqinfo(gr,genome,celltype)
    strand(gr) <- "*"
    IDX <- unique(findOverlaps(gr,all.gr)@to)
    apa_DF[IDX,nm] <- 1
  }
  
  
  
  
  eRNA_MTR4.hg19.gr <- myOverlaps("ncRNAMTR4","PROM",ncRNA_MTR4kd.hg19.gr,enhCAGE.hg19.gr)$gr1xgr2
  eRNA_Z1.hg19.gr <- myOverlaps("ncRNAZ1","PROM",ncRNA_Z1kd.hg19.gr,enhCAGE.hg19.gr)$gr1xgr2
  eRNA_Z8.hg19.gr <- myOverlaps("ncRNAZ8","PROM",ncRNA_Z8kd.hg19.gr,enhCAGE.hg19.gr)$gr1xgr2
  
  
  myWhole.l <- c("z1.hg19.gr",MTR4.cgr,Z8.cgr,"ctcf.hg19.gr","rad21RK.hg19.gr",
                 "promptsZ1.hg19.gr","promptsZ8.hg19.gr","promptsMTR4.hg19.gr",
                 "ncRNA_Z1kd.hg19.gr","ncRNA_Z8kd.hg19.gr","ncRNA_MTR4kd.hg19.gr",
                 "eRNA_MTR4.hg19.gr","eRNA_Z1.hg19.gr","eRNA_Z8.hg19.gr",
                 "TSS.hg19.gr","enhCAGE.hg19.gr","h3k4me1.hg19.gr")
  
  
  
  xtVec <- c(T,T,T,T,
             T,T,T,
             F,F,F,
             F,F,F,
             F,F,F,
             T,F,T)
  
  xtend <- 500
  genome <- "hg19"
  celltype <- "hela"
  
  all.gr <- trim(GenomicRanges::reduce(do.call("c",sapply(1:length(myWhole.l),function(X){
    gr <- get(myWhole.l[X])
    mcols(gr) <- NULL
    if(xtVec[X]){print("REDUCE"); ret_res <- addSeqinfo(resize(gr,1,"center")+xtend,genome,celltype);strand(ret_res) <- "*"}
    else{ret_res <- addSeqinfo(gr,genome,celltype);strand(ret_res) <- "*"}
    return(ret_res)}))))
  
  enhapa_DF <- data.frame(IDX=1:length(all.gr),ALL=1)
  for(feat in myWhole.l){
    nm <- toupper(gsub(paste0("(.*).",genome,".gr"),"\\1",feat))
    enhapa_DF[,nm] <- 0
    gr <- get(feat)
    mcols(gr) <- NULL
    gr <- addSeqinfo(gr,genome,celltype)
    strand(gr) <- "*"
    IDX <- unique(findOverlaps(gr,all.gr)@to)
    enhapa_DF[IDX,nm] <- 1
  }
  
  enhapa_dt <- data.table::data.table(enhapa_DF)
  
  
  
  
  apa_dts <- data.table::data.table(apa_DF)
  # IF tes or tss you cant be gene body
  vv <- apa_dts$GENESBODY
  
  
  feats <- c("Z1","MTR4","Z8")
  
  for(featv in feats){
    
    apa_dt <- copy(apa_dts)
    
    feat <- toupper(paste0("ncRNA_",featv,"kd"))
    
    colname <- paste0("PROMPTS", featv)
    apa_dt$INTERGENIC <- 0
    apa_dt[,GENESBODY:=ifelse(.SD[[colname]]==1 ,0,GENESBODY)]
    apa_dt[, INTERGENIC := ifelse(GENESBODY == 0 & .SD[[colname]] == 0 & .SD[[feat]]==1, 1, 0)]
    
    mydf <- NULL
    tot_feat <- nrow(apa_dt[get(feat),])
    
    
    tst.l <- toupper(c(paste0("prompts",featv),
                       "genesBody","intergenic"))
    
    
    for(tst_feat in tst.l){
      
      tot_tst <- nrow(apa_dt[get(feat) & get(tst_feat),])
      percF <- tot_tst*100/tot_feat
      mydf <- rbind(mydf,
                    rbind.data.frame(cbind.data.frame(nm=tst_feat, total = tot_tst, perc=percF, type="PROT_%")))
      
    }
    
    tot_sub <- nrow(enhapa_dt[get(paste0("ERNA_",featv))==1,])
    percF <- tot_sub*100/tot_feat
    mydf <- rbind(mydf,
                  rbind.data.frame(cbind.data.frame(nm="eRNAs", total = tot_sub, perc=percF, type="PROT_%")))
    
    
    tss_dt <- apa_dt[TSS==1,]
    tot_tss <- nrow(tss_dt)
    tot_sub <- nrow(tss_dt[get(feat)==1,])
    percF <- tot_sub*100/tot_tss
    mydf <- rbind(mydf,
                  rbind.data.frame(cbind.data.frame(nm="PERC_TSS_NCRNA", total = tot_sub, perc=percF, type="PROT_%")))
    
    
    
    p <- ggbarplot(
      mydf,
      x = "nm",
      y = "perc",
      fill = Col.l[[featv]],
      position = position_dodge(0.9),
      orientation = "horiz",
      width = 0.4,
      palette = "jco"
    ) +
      # ajout du total au bout de chaque barre
      geom_text(
        aes(label = total),
        position = position_dodge(width = 0.9),
        hjust = -0.2,     # décale légèrement à droite du bout de la barre
        size = 3.5
      )
    
    # affichage propre
    gg[[paste0("Barplot.NCRNA_Percentage.ncRNA_",featv)]] <- print(ggpar(p, xlim = c(0, 105), rotate = TRUE))
    
  }
} else {
  
  # FIGURE 1 C-D-E - BARPLOT OF CHIPSEQ PROTEINS vs TSS,TES,ENH,H3K4ME1,GENEbODY, INTERGENIC ----------------------------------------------------------
  TSS_D_S.l <- list(TSS_SGL="tss.hg19.sglTSS.pm2kb.gr")
  
  TSS.chr  <-  names(TSS_D_S.l)[1]
  MTR4.chr  <-  names(MTR4.l)[2]
  Z8.chr  <-  names(Z8.l)[2]
  
  TSS.hg19.gr <- resize(genes.hg19.sglTSS.pm2kb.gr,1,"start")
  TES.hg19.gr <- resize(genes.hg19.sglTSS.pm2kb.gr,1,"end")
  
  prom2Kb.hg19.gr <- extendGR(TSS.hg19.gr,2000,0)
  
  promptsZ1.hg19.gr <- myOverlaps("ncRNAZ1","PROM",ncRNA_Z1kd.hg19.gr,prom2Kb.hg19.gr)$gr1xgr2
  promptsMTR4.hg19.gr <- myOverlaps("ncRNAMTR4","PROM",ncRNA_MTR4kd.hg19.gr,prom2Kb.hg19.gr)$gr1xgr2
  promptsZ8.hg19.gr <- myOverlaps("ncRNAZ8","PROM",ncRNA_Z8kd.hg19.gr,prom2Kb.hg19.gr)$gr1xgr2
  
  rad21RK.hg19.gr <- get(RAD21.l$RAD21_NEW)
  
  
  MTR4.cgr <- MTR4.l[[MTR4.chr]]
  MTR4.hg19.gr <- get(MTR4.cgr)
  MTR.nm <- toupper(gsub("(.*).hg19.gr","\\1",MTR4.cgr))
  Z8.cgr <- Z8.l[[Z8.chr]]
  Z8.nm <- toupper(gsub("(.*).hg19.gr","\\1",Z8.cgr))
  
  
  MTR4xZ1.hg19.gr <- subsetByOverlaps(z1.hg19.gr,get(MTR4.cgr),maxgap=500)
  MTR4xZ8.hg19.gr <- subsetByOverlaps(get(Z8.cgr),get(MTR4.cgr),maxgap=500)
  
  Col.l <- list(MTR4N="#ff5733",Z1="#1898dd",Z8N="#7bb84f",MTR4_inter_Z1="#45b",MTR4_inter_Z8="#bd8841")
  
  
  myWhole.l <- c("z1.hg19.gr",MTR4.cgr,Z8.cgr,"MTR4xZ1.hg19.gr","MTR4xZ8.hg19.gr","ctcf.hg19.gr","rad21RK.hg19.gr",
                 "promptsZ1.hg19.gr","promptsZ8.hg19.gr","promptsMTR4.hg19.gr",
                 "ncRNA_Z1kd.hg19.gr","ncRNA_Z8kd.hg19.gr","ncRNA_MTR4kd.hg19.gr",
                 "TSS.hg19.gr","enhCAGE.hg19.gr","h3k4me1.hg19.gr")
  
  
  xtVec <- c(T,T,T,T,
             T,T,T,
             F,F,F,
             F,F,F,
             T,F,T)
  
  xtend <- 500
  genome <- "hg19"
  celltype <- "hela"
  
  all.gr <- trim(GenomicRanges::reduce(do.call("c",sapply(1:length(myWhole.l),function(X){
    gr <- get(myWhole.l[X])
    mcols(gr) <- NULL
    if(xtVec[X]){print("REDUCE"); ret_res <- addSeqinfo(resize(gr,1,"center")+xtend,genome,celltype);strand(ret_res) <- "*"}
    else{ret_res <- addSeqinfo(gr,genome,celltype);strand(ret_res) <- "*"}
    return(ret_res)}))))
  
  apa_DF <- data.frame(IDX=1:length(all.gr),ALL=1)
  for(feat in myWhole.l){
    nm <- toupper(gsub(paste0("(.*).",genome,".gr"),"\\1",feat))
    apa_DF[,nm] <- 0
    gr <- get(feat)
    mcols(gr) <- NULL
    gr <- addSeqinfo(gr,genome,celltype)
    strand(gr) <- "*"
    IDX <- unique(findOverlaps(gr,all.gr)@to)
    apa_DF[IDX,nm] <- 1
  }
  
  
  comp.l <- c("TES.hg19.gr")
  xtV <- c(T)
  for(ft in comp.l[xtV]){
    nm <- gsub("(.*).hg19.gr","\\1",ft)
    ft.gr <- get(ft)
    mcols(ft.gr) <- NULL
    assign(paste0(nm,".xt.hg19.gr"),trim(addSeqinfo(resize(ft.gr,1,"center") + xtend,"hg19","hela")))
  }
  
  genesBody.xt.hg19.gr <- addSeqinfo(genesBody.hg19.gr ,"hg19","hela")
  mcols(genesBody.xt.hg19.gr) <- NULL
  seqlevelsStyle(genesBody.xt.hg19.gr) <- "Ensembl"
  
  xtcomp.l <- c("TES.xt.hg19.gr", "genesBody.xt.hg19.gr")
  
  for(feat in xtcomp.l){
    nm <- toupper(gsub(paste0("(.*).xt.",genome,".gr"),"\\1",feat))
    apa_DF[,nm] <- 0
    gr <- get(feat)
    mcols(gr) <- NULL
    gr <- addSeqinfo(gr,genome,celltype)
    strand(gr) <- "*"
    IDX <- unique(findOverlaps(gr,all.gr)@to)
    apa_DF[IDX,nm] <- 1
  }
  
  
  
  apa_dt <- data.table::data.table(apa_DF)
  # IF tes or tss you cant be gene body
  apa_dt[,GENESBODY:=ifelse((TSS==1 & GENESBODY==1)|(TES==1 & GENESBODY==1),0,GENESBODY)]
  apa_dt[,TES:=ifelse((TES==1 & TSS==1),0,TES)]
  apa_dt$INTERGENIC <- 0
  apa_dt[,INTERGENIC:=ifelse(GENESBODY==1 | TSS==1 | TES == 1,0,1)]
  
  
  tst.l <- toupper(c("h3k4me1","TSS","TES",
                     "enhCAGE",
                     "genesBody","intergenic"))
  
  
  
  apa_dt[,MTR4_inter_Z1:= Z1 & MTR4N]
  apa_dt[,MTR4_inter_Z8:= Z8N & MTR4N]
  
  
  
  apa_dt[MTR4N==1]
  
  feat.l <- c("MTR4N","Z1","Z8N","MTR4_inter_Z1","MTR4_inter_Z8")
  
  for(feat in feat.l){
    mydf <- NULL
    tot_feat <- nrow(apa_dt[get(feat),])
    
    for(tst_feat in tst.l){
      
      tot_tst <- nrow(apa_dt[get(feat) & get(tst_feat),])
      percF <- tot_tst*100/tot_feat
      mydf <- rbind(mydf,
                    rbind.data.frame(cbind.data.frame(nm=tst_feat, total = tot_tst, perc=percF, type="PROT_%")))
      
    }
    
    
    tot_tst <- nrow(apa_dt[get(feat) & get(tst_feat),])
    
    sub_dt <- apa_dt[get(feat)==1 & INTERGENIC==1,]
    tot_sub <- nrow(sub_dt)
    tot_tst <- nrow(sub_dt[H3K4ME1==1,])
    percF <- tot_tst*100/tot_sub
    mydf <- rbind(mydf,
                  rbind.data.frame(cbind.data.frame(nm="Intergenic_H3K4me1_Enhancers", total = tot_tst, perc=percF, type="PROT_%")))
    
    
    p <- ggbarplot(
      mydf,
      x = "nm",
      y = "perc",
      fill = Col.l[[feat]],
      position = position_dodge(0.9),
      orientation = "horiz",
      width = 0.4,
      palette = "jco"
    ) +
      # ajout du total au bout de chaque barre
      geom_text(
        aes(label = total),
        position = position_dodge(width = 0.9),
        hjust = -0.2,     # décale légèrement à droite du bout de la barre
        size = 3.5
      )
    
    # affichage propre
    gg[[paste0("Barplot.ChIPseq_Percentage.",feat)]] <- print(ggpar(p, xlim = c(0, 105), rotate = TRUE))
  }
}

#### SAVE PLOTS ----
out_dir <- "STAR_METHODS/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (nm in names(gg)) {
  if (!is.null(gg[[nm]])) {
    ggsave(file.path(out_dir, paste0("fig2.A.", nm, ".png")),
           plot = gg[[nm]], width = 8, height = 6, dpi = 300)
  }
}


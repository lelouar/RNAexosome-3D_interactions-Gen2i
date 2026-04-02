# Load needed libraries ----
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(ggpubr)
# library("plot3D")
library(plyr)
library(dplyr)
library(reshape2)
library(Matrix)
library("stringr")
library("gridExtra")
require("eulerr")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
# library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(scales)
# library("seqplots")
# library(reshape)
library(RColorBrewer)
library(Vennerable)
library(trackViewer)
library(tidyverse)
library(foreach)
library(doParallel)
numCores <- detectCores()


# access to node2

#  import("/run/user/1003/gvfs/smb-share:server=node2,share=grpcuvier/Alexandre/BED_Kc167_RM/Prot/Barren_Kc167_CS_dm3.bed")

#  GLOBAL VARIABLE ---------------------------------------------------------


# GRANGES MANIPULATION ----------------------------------------------------


# To find divergent genes we have to reverse TSS minus strand


# LIFTOVER
# dm3TOdm6 <- import.chain("resources/ext_data/Database/Drosophila/Genome_info/liftover_chain/dm3ToDm6.over.chain")
# dm6TOdm3 <- import.chain("resources/ext_data/Database/Drosophila/Genome_info/liftover_chain/dm6ToDm3.over.chain")
# hg19TOhg38 <- import.chain("resources/ext_data/Database/Human/Genome_info/liftover_chain/hg19ToHg38.over.chain")
# hg38TOhg19 <- import.chain("resources/ext_data/Database/Human/Genome_info/liftover_chain/hg38ToHg19.over.chain")





liftover <- function(my.gr,chain,sqLvlstyle="Ensembl"){
  sequence <- seqlevels(my.gr)
  ## sequence is in UCSC format and we want NCBI style
  newStyle <- mapSeqlevels(sequence,"UCSC")
  newStyle <- newStyle[complete.cases(newStyle)] # removing NA cases.
  ## rename the seqlevels
  my.gr <- renameSeqlevels(my.gr,newStyle)
  my.new.gr <- unlist(liftOver(my.gr,chain))
  seqlevelsStyle(my.new.gr) = sqLvlstyle
  return(my.new.gr)
}

# # HIC LIST STORAGE
# myhiclist <- list(eagenKC167="resources/ext_data/Database/Drosophila/Kc167/dm3/hic/Eagen/GSE89112_Kc167combined_randomized.hic",
#                   wangS2="resources/ext_data/Database/Drosophila/S2R/hic/Wang_juicer/WangS2R_noFilt_Juicer.hic",
#                   helaHG19="https://hicfiles.s3.amazonaws.com/hiseq/hela/in-situ/combined.hic",
#                   helaHG19_W="resources/ext_data/Database/Human/helaS3/hg19/hic/Hela/Stocsits/GSE102884_RAW/GSM2747739_HeLa_W_2reps.hic",
#                   helaHG19_AB="resources/ext_data/Database/Human/helaS3/hg19/hic/Hela/Stocsits/GSE102884_RAW/GSM2747740_HeLa_AB_2reps.hic",
#                   helaHG19_WAB="resources/ext_data/Database/Human/helaS3/hg19/hic/Hela/Stocsits/GSE102884_RAW/GSM2747741_HeLa_WAB_2reps.hic",
#                   helaHG19_CTL="resources/ext_data/Database/Human/helaS3/hg19/hic/Hela/Stocsits/GSE102884_RAW/GSM2747738_HeLa_control_2reps.hic",
#                   ramirezS2bKD_500="resources/ext_data/Database/Drosophila/S2/dm3/hic/Ramirez/Beaf_KD/inter_500.hic",
#                   ramirezS2WT_500="resources/ext_data/Database/Drosophila/S2/dm3/hic/Ramirez/WT/inter_500.hic",
#                   ramirezS2bKD_1000="resources/ext_data/Database/Drosophila/S2/dm3/hic/Ramirez/Beaf_KD/inter_1000.hic",
#                   ramirezS2WT_1000="resources/ext_data/Database/Drosophila/S2/dm3/hic/Ramirez/WT/inter_1000.hic",
#                   tcellHG19="https://s3.amazonaws.com/hicfiles/external/goodell/tcell.hic")


`%ni%` <-  Negate('%in%')
# DM 3
kc.chr <- c("2L","2R","3L","3R","X")
s2.chr <- c("2L","2R","3L","3R","X")
sqfDM3 <- SeqinfoForBSGenome(genome="dm3")
sqfDM6 <- SeqinfoForBSGenome(genome="dm6")
dm3.gr <- GRanges(sqfDM3)
dm6.gr <- GRanges(sqfDM6)

myCHR <- c("2L","2R","3L","3R","X")


# Human
hela.chr <- c(seq(1,22,1),"X")
sqfHG19 <- sqfHG19 <- g2i:::fai_lst$NCBI$GRCh37
sqfHG38 <- sqfHG19 <- g2i:::fai_lst$NCBI$GRCh37
hg19.gr <- GRanges(sqfHG19)
hg38.gr <- GRanges(sqfHG38)




MYGENOME <- list(hg19=c(txdb=TxDb.Hsapiens.UCSC.hg19.knownGene,hela.chr=list(hela.chr),sqf=sqfHG19))# hg38=c(txdb=TxDb.Hsapiens.UCSC.hg38.knownGene,hela.chr=list(hela.chr),sqf=sqfHG38)


create <- function (...){
  tmp <- list(...)
  dir.create(paste0(...),showWarnings = F,recursive = T)
  
  return(paste0(...))
}

plotWindow <- list(PM2KB=rep(2000,2),
                   PM1KB=rep(1000,2),
                   PM500B=rep(500,2))

profcomp2_tx=function(covr,gr)
  #Compute profiles from a genome-wide coverage and a set of windows
  #covr is a genome-wide coverage (typically obtained with the coverage function)
  #gr is a GRanges object containing the windows over which to compute the profiles
  #Note that names of covr and seqnames of gr must match
{
  prof=covr[gr]
  prof[strand(gr)=='-']=lapply(prof[strand(gr)=='-'],rev)
  return(prof)
}

saveBRDS <- function(obj,s_path){
  saveRDS(obj,paste0(s_path,".rds"))
  export.bed(obj,paste0(s_path,".bed"))
}


changeSeqlevels <- function(my.gr, levels="UCSC"){
  sequence <- seqlevels(my.gr)
  ## sequence is in UCSC format and we want NCBI style
  newStyle <- mapSeqlevels(sequence,levels)
  newStyle <- newStyle[complete.cases(newStyle)] # removing NA cases.
  ## rename the seqlevels
  my.gr <- renameSeqlevels(my.gr,newStyle)
}

addSeqinfo <- function(my.gr,g_name="dm3",celltype="kc",keepSeq=T,chrTYPE="Ensembl"){
  sqf <- MYGENOME[[g_name]]$sqf
  
  chr <- as.character(unlist(MYGENOME[[g_name]][[paste0(celltype,".chr")]]))
  seqlevelsStyle(my.gr) <- "UCSC"
  seqinfo(my.gr) <- sqf[seqlevels(my.gr),]
  seqlevelsStyle(my.gr) <- "UCSC"
  seqinfo(my.gr) <- sqf[seqlevels(my.gr),]
  seqlevelsStyle(my.gr) <- "Ensembl"
  if(keepSeq) my.gr <- keepSeqlevels(my.gr,chr,pruning.mode="coarse")
  my.gr <- sortSeqlevels(my.gr)
  my.gr <- sort(my.gr,ignore.strand=T)
  seqlevelsStyle(my.gr) <- chrTYPE
  return(my.gr)
}




hg19.tiles.gr <- unlist(tileGenome(sqfHG19,tilewidth = 1000))
hg19.tiles.gr <- addSeqinfo(hg19.tiles.gr,g_name = "hg19",celltype = "hela",chrTYPE = "Ensembl")
hg19.tiles2KB.gr <- unlist(tileGenome(sqfHG19,tilewidth = 2000))
hg19.tiles2KB.gr <- addSeqinfo(hg19.tiles.gr,g_name = "hg19",celltype = "hela",chrTYPE = "Ensembl")


facToInt <- function(VEC) {
  as.numeric(as.character(VEC))
  
}

uniqueOvlp <- function(q.gr,s.gr,mxgp){
  return(q.gr[unique(queryHits(findOverlaps(q.gr,s.gr,maxgap=mxgp)))])
}

extendGR <- function(x, upstream=0, downstream=0)     
{
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}
reduceGR <- function(x, upstream=0, downstream=0)     
{
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) + ifelse(on_plus, upstream, downstream)
  new_end <- end(x) - ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}


# CREATE RANDOM GENOMIC RANGES BASED ON SEQINFO DATA
randomGR <- function(sqINFO=NULL,n=1000) {
  Reduce("c",sample(tileGenome(sqINFO,tilewidth = 200),n,F))
}


trim_vec <- function(X,perc,NA.RM=T){
  TMP.df <- data.frame(vec=X) %>% mutate(dec=ntile(vec,100))
  TMPsub.df <- subset(TMP.df,dec > perc & dec < (100-perc))
  if(NA.RM) TMPsub.df <- subset(TMPsub.df,!is.na(vec))
  RET.v <- TMPsub.df$vec
  return(RET.v)
}   


bordersGR <- function(GR) {
  s.gr <- resize(GR,1,"start")
  e.gr <- resize(GR,1,"end")
  a.gr <- c(e.gr,s.gr)
  a.gr <- sort(a.gr)
  return(list(allBORD=a.gr,leftBORD=s.gr,rightBORD=e.gr))
  
}

# CREATE BOOLEAN DATA FRAME OF A GIVEN SET OF CHIPSEQ POSITION REGARDING A TILED .gr OBJECT 
genomic_DF <- function(g_set=NULL,pattern=".hg19.gr",tiled.gr=hg19.tiles.gr,g_ver="hg19",c_type="hela") {
  l_names <- toupper(gsub(paste0("(.*)",pattern),"\\1",g_set))
  ezvec <- rep(0,length(tiled.gr))
  registerDoParallel(numCores)  # use multicore, set to the number of our cores
  as.data.frame(foreach::foreach(X=g_set,.final = function(x) setNames(x, l_names)) %dopar% {
    mygr <- get(X)
    mcols(mygr) <- NULL
    mygr <- addSeqinfo(mygr,g_ver,c_type)
    strand(mygr) <- "*"
    IDX <- findOverlaps(mygr,tiled.gr)@to
    ezvec[IDX] <- 1
    ezvec[,drop=FALSE]
  })
}

createGRcouplesAPA <- function(apa_size=41,hic_res=1000,anc.gr=NULL,int.gr=NULL,constraint.gr=NULL,SQF=NULL,MYCHR=dm3CHR){
  
  if(is.null(constraint.gr)){
    print("NO DOMAIN TO RESTRICT INTERACTIONS, USE ANOTHER FUNCTION OR PROVIDE DOMAIN GRANGES")
    return(F)
  }
  # DELETE MCOLS TO AVOID FIRST.X.BLAHBLAHBLAH
  mcols(anc.gr) <- NULL
  mcols(int.gr) <- NULL
  mcols(constraint.gr) <- NULL
  names(anc.gr) <- NULL
  names(int.gr) <- NULL
  names(constraint.gr) <- NULL
  seqlevelsStyle(anc.gr) <- "Ensembl"
  seqlevelsStyle(int.gr) <- "Ensembl"
  seqlevelsStyle(constraint.gr) <- "Ensembl"
  
  # Bin size 
  mBin <- (apa_size-1)/2
  #Interacting feature
  #Interacting feature
  feati.gr <- resize(int.gr,1,fix="center") #resize the feature to do the overlap on domains without multi hits
  feati.gr$bin_i_IDX <- ceiling((start(feati.gr))/hic_res)
  feati.gr$i_IDX <- 1:length(feati.gr)
  dfi <- data.frame(findOverlapPairs(constraint.gr,feati.gr))
  dfi$filt <- paste0(dfi$first.seqnames,"_",dfi$first.start)
  dfi <- dfi[,c("second.bin_i_IDX","second.i_IDX","filt")]
  colnames(dfi) <- c("bin_i_IDX","i_IDX","cons_info")
  
  #Anchoring feature
  feata.gr <- resize(anc.gr,1,fix="center") #resize the feature to do the overlap on domains without multi hits
  feata.gr$bin_a_IDX <- ceiling((start(feata.gr))/hic_res)
  feata.gr$a_IDX <- 1:length(feata.gr)
  dfa <- data.frame(findOverlapPairs(constraint.gr,feata.gr))
  dfa$filt <- paste0(dfa$first.seqnames,"_",dfa$first.start)
  dfa <- dfa[,c("second.bin_a_IDX","second.a_IDX","filt")]
  colnames(dfa) <- c("bin_a_IDX","a_IDX","cons_info")
  
  # Merge dataframes
  mia <- merge(dfi,dfa,by = "cons_info")
  # Remove couples that are < 15kb/binS fragment away
  if(length(which(abs(mia$bin_i_IDX-mia$bin_a_IDX)<(15000/hic_res))))
    mia <- mia[-which(abs(mia$bin_i_IDX-mia$bin_a_IDX)<(15000/hic_res)),]
  
  # Add distance between int and anc in bin 
  mia$dst_bin <- abs(mia$bin_i_IDX-mia$bin_a_IDX)
  # Add chromosome info
  mia$chr <- gsub("(.*)_.*","\\1",mia$cons_info)
  
  
  
  # remove peaks that fall less than bin+1 away from start and end on each chromosome
  MYDF <- NULL
  for(CHR in MYCHR){
    print(CHR)
    miaTMP <- subset(mia,chr==CHR)
    DIM <- as.numeric(ceiling(seqlengths(SQF[paste0("chr",CHR)])/hic_res))
    toDEL <- which(miaTMP$bin_i_IDX >= DIM-mBin | miaTMP$bin_i_IDX <= mBin | miaTMP$bin_a_IDX >= DIM-mBin | miaTMP$bin_a_IDX <= mBin  )
    if(length(toDEL)) miaTMP <- miaTMP[-toDEL,]
    MYDF <- rbind(MYDF,miaTMP)
  }
  return(MYDF)
}
createGRcouplesAPA.dst <- function(apa_size=41,hic_res=1000,anc.gr=NULL,int.gr=NULL,constraint.gr=NULL,SQF=NULL,MYCHR=dm3CHR, dmax = NULL, dmin = NULL){
  
  if(is.null(constraint.gr)){
    # chromosomes collection
    chrs <- unique(c(as.vector(seqnames(int.gr)),as.vector(seqnames(anc.gr))))
    Go_cons <- NULL
    for (chr in chrs) {
      min_anc <- min(start(subset(anc.gr, seqnames == chr)))
      max_anc <- max(end(subset(anc.gr, seqnames == chr)))
      min_bait <- min(start(subset(int.gr, seqnames == chr)))
      max_bait <- max(end(subset(int.gr, seqnames == chr)))
      Go_cons <- rbind(Go_cons, data.frame(start = min(min_anc,min_bait), end = max(max_anc,max_bait) ,seqnames = chr))
      constraint.gr <- GenomicRanges::GRanges(Go_cons)
    }
  }
  # DELETE MCOLS TO AVOID FIRST.X.BLAHBLAHBLAH
  mcols(anc.gr) <- NULL
  mcols(int.gr) <- NULL
  mcols(constraint.gr) <- NULL
  names(anc.gr) <- NULL
  names(int.gr) <- NULL
  names(constraint.gr) <- NULL
  seqlevelsStyle(anc.gr) <- "Ensembl"
  seqlevelsStyle(int.gr) <- "Ensembl"
  seqlevelsStyle(constraint.gr) <- "Ensembl"
  
  # Bin size 
  mBin <- (apa_size-1)/2
  #Interacting feature
  #Interacting feature
  feati.gr <- resize(int.gr,1,fix="center") #resize the feature to do the overlap on domains without multi hits
  feati.gr$bin_i_IDX <- ceiling((start(feati.gr))/hic_res)
  feati.gr$i_IDX <- 1:length(feati.gr)
  dfi <- data.table::data.table(data.frame(findOverlapPairs(constraint.gr,feati.gr)))
  dfi$filt <- paste0(dfi$first.seqnames,"_",dfi$first.start)
  dfi <- dfi[,c("second.bin_i_IDX","second.i_IDX","filt")]
  colnames(dfi) <- c("bin_i_IDX","i_IDX","cons_info")
  data.table::setkey(dfi, cons_info)
  
  #Anchoring feature
  feata.gr <- resize(anc.gr,1,fix="center") #resize the feature to do the overlap on domains without multi hits
  feata.gr$bin_a_IDX <- ceiling((start(feata.gr))/hic_res)
  feata.gr$a_IDX <- 1:length(feata.gr)
  dfa <- data.table::data.table(data.frame(findOverlapPairs(constraint.gr,feata.gr)))
  dfa$filt <- paste0(dfa$first.seqnames,"_",dfa$first.start)
  dfa <- dfa[,c("second.bin_a_IDX","second.a_IDX","filt")]
  colnames(dfa) <- c("bin_a_IDX","a_IDX","cons_info")
  data.table::setkey(dfa, cons_info)
  
  # Merge dataframes
  mia <-   dfi[dfa, allow.cartesian=TRUE,nomatch=0]
  # Add distance between int and anc in bin 
  mia$dst_bin <- abs(mia$bin_i_IDX-mia$bin_a_IDX)
  
  if (!is.null(dmin)) mia <- mia[dst_bin>=dmin ]
  if (!is.null(dmax)) mia <- mia[dst_bin<=dmax ]
  
  # Add chromosome info
  mia$chr <- gsub("(.*)_.*","\\1",mia$cons_info)
  
  # remove peaks that fall less than bin+1 away from start and end on each chromosome
  MYDF <- NULL
  for(CHR in MYCHR){
    print(CHR)
    miaTMP <- subset(mia,chr==CHR)
    DIM <- as.numeric(ceiling(seqlengths(SQF[CHR])/hic_res))
    toDEL <- which(miaTMP$bin_i_IDX >= DIM-mBin | miaTMP$bin_i_IDX <= mBin | miaTMP$bin_a_IDX >= DIM-mBin | miaTMP$bin_a_IDX <= mBin  )
    if(length(toDEL)) miaTMP <- miaTMP[-toDEL,]
    MYDF <- rbind(MYDF,miaTMP)
  }
  return(data.frame(MYDF))
}

# DATA IMPORT EXAMPLE

# source("https://bioconductor.org/biocLite.R")
# biocLite("trackViewer")
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(org.Hs.eg.db)
# gr <- GRanges("chr11", IRanges(122929275, 122930122), strand="-")
# trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
#                          org.Hs.eg.db,
#                          gr=gr)
# 
# library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
# mybw <- "resources/ext_data/Database/Drosophila/Kc167/dm3/protein_binding/cp190/CP190.bw"
# tss.gr <- genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
# seqlevelsStyle(tss.gr) <- "UCSC"
# 
# fox2 <- importData(mybw, format="BigWig",
#                    ranges=tss.gr)
# ok <- fox2$`resources/ext_data/Database/Drosophila/Kc167/dm3/protein_binding/cp190/CP190.bw`
# res <- ok[tss.gr]
# ok <- lapply(res,function(X){
#   colMeans(as.matrix(X))
# })
# 
# 
# fox3 <- importData(mybw, format="BigWig")
# 
# repA <- importScore(file.path(extdata, "cpsf160.repA_-.wig"),
#                     file.path(extdata, "cpsf160.repA_+.wig"),
#                     format="WIG")

dotHICplot <- function(APA.mtx){
  ok.df <- APA.mtx
  colnames(ok.df) <- 1:ncol(ok.df)
  rownames(ok.df) <- 1:nrow(ok.df)
  
  tst <- melt(ok.df)
  colnames(tst) <- c("ROW","COL","VAL")
  tst <- tst %>% mutate(QUANT=ntile(VAL,150))
  COLTST <- (colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(150))
  tst$COLOR <- COLTST[tst$QUANT]
  tst$ROW <- as.factor(tst$ROW)
  tst$COL <- as.factor(tst$COL)
  
  ggplot(tst, aes(y = rev(ROW),x =COL)) +
    # geom_tile(aes(fill = rectheat)) +         ## to get the rect filled
    geom_point(aes(colour = VAL,size=(VAL)^15),shape=20)  +    ## geom_point for circle illusion
    theme_minimal() + theme(plot.background = element_rect(fill="black"),
                            panel.grid.major = element_blank()) +
    # scale_color_discrete(COLTST)+
    scale_color_gradientn(colors=COLTST)+       ## color of the corresponding aes
    scale_size(range = c(1, 8))+             ## to tune the size of circles
    theme_tiff
}







zscore_array <- function(X) {
  M <- abs(mean(c(X[1],X[2]),na.rm=T))
  return((X[1]-X[2])/sqrt(M))  
}

zscore_df <- function(X,Y){
  M <- abs(mean(c(X,Y),na.rm=T))
  return((X-Y)/sqrt(M))  
}

hic_pipeline <- function(l_IN=l_IN,apaCPLS.df=NULL,mBin=20) {
  
  
  HIC.lmtx <- l_IN$HIC
  res <- l_IN$RES
  NM <- l_IN$NM
  
  
  # EXCTRACT HIC MATRICES BINxBIN
  myMAT <- list()
  numCores <- detectCores()
  registerDoParallel(numCores)  # use multicore, set to the number of our cores
  myMAT <- foreach::foreach(X=1:nrow(apaCPLS.df)) %dopar% {
    binANC <- apaCPLS.df[X,"bin_a_IDX"]
    binINT <- apaCPLS.df[X,"bin_i_IDX"]
    CHR <- apaCPLS.df[X,"chr"]
    if(binINT < binANC){
      TMP <- binINT
      binINT <- binANC
      binANC <- TMP
    }
    # Compute the APA matrix and return it in myMAT as a list of matrix
    feat_apa.mtx <- Matrix(HIC.lmtx[[CHR]][((binANC-mBin):(binANC+mBin)),((binINT-mBin):(binINT+mBin))],sparse=F)
    feat_apa.mtx[!is.finite(feat_apa.mtx) | feat_apa.mtx==0] <- NA
    list(a_IDX=apaCPLS.df[X,"a_IDX"],i_IDX=apaCPLS.df[X,"i_IDX"],MAT=as.matrix(feat_apa.mtx))
  }
  stopImplicitCluster()
  
  
  
  
  # CREATE DATA FRAME OF METRICS
  myMAT.l <- list()
  numCores <- detectCores()
  registerDoParallel(numCores)  # use multicore, set to the number of our cores
  myMAT.l <- foreach::foreach(X=1:length(myMAT)) %dopar% {
    a_IDX <- myMAT[[X]][["a_IDX"]]
    i_IDX <- myMAT[[X]][["i_IDX"]]
    MAT <- myMAT[[X]][["MAT"]]
    l_ <- list()
    l_[["CP_V"]] <- MAT[i_CP,j_CP]
    l_[["CS_V"]] <- MAT[i_CS,j_CS]
    l_[["ULS_V"]] <- MAT[i_ULS,j_ULS]
    l_[["URS_V"]] <- MAT[i_URS,j_URS]
    l_[["BLS_V"]] <- MAT[i_BLS,j_BLS]
    l_[["BRS_V"]] <- MAT[i_BRS,j_BRS]
    TMP.df <- data.frame(a_IDX=a_IDX,i_IDX=i_IDX)
    for(PAR in names(l_)){
      TMP.df[,PAR] <- mean(as.array(l_[[PAR]]),na.rm=T)
    }
    TMP.df
  }
  stopImplicitCluster()
  METRICS.df <- do.call("rbind",myMAT.l)
  
  
  myMAT_QUANT.l <- list()
  numCores <- detectCores()
  registerDoParallel(numCores)  # use multicore, set to the number of our cores
  
  
  myMAT_QUANT.l <- foreach::foreach(X=1:length(myMAT)) %dopar% {
    a_IDX <- myMAT[[X]][["a_IDX"]]
    i_IDX <- myMAT[[X]][["i_IDX"]]
    MAT <- myMAT[[X]][["MAT"]]
    MAT_QUANT <- matrix(ntile(MAT,500),APA_SIZE,APA_SIZE)
    l_ <- list()
    l_[["CP_V"]] <- MAT_QUANT[i_CP,j_CP]
    l_[["CS_V"]] <- MAT_QUANT[i_CS,j_CS]
    l_[["ULS_V"]] <- MAT_QUANT[i_ULS,j_ULS]
    l_[["URS_V"]] <- MAT_QUANT[i_URS,j_URS]
    l_[["BLS_V"]] <- MAT_QUANT[i_BLS,j_BLS]
    l_[["BRS_V"]] <- MAT_QUANT[i_BRS,j_BRS]
    TMP.df <- data.frame(a_IDX=a_IDX,i_IDX=i_IDX)
    for(PAR in names(l_)){
      TMP.df[,PAR] <- mean(as.array(l_[[PAR]]),na.rm=T)
    }
    TMP.df
  }
  stopImplicitCluster()
  METRICS_QUANT.df <- do.call("rbind",myMAT_QUANT.l)
  
  
  return(list(METRICS.df=METRICS.df,METRICS_QUANT.df=METRICS_QUANT.df,myMAT=myMAT))
  
  
}






find_divergent_genes <- function(tssm.gr,tssp.gr,maxgap=1000,geneID=F){
  tss.gr <- c(tssm.gr,tssp.gr)
  tssmp.gr <- tssm.gr
  strand(tssmp.gr) <- "+"
  idxTssp <- precede(tssmp.gr,tssp.gr)
  # Delete TssM that are at the end of a chromosome
  toDel <- which(is.na(idxTssp))
  tssmp.gr <- tssmp.gr[-toDel]
  idxTsspnoNA <- idxTssp[-toDel]
  
  tssp_following_tssm.gr <- tssp.gr[idxTsspnoNA]
  
  dst <- start(tssp_following_tssm.gr) - start(tssmp.gr) 
  
  tssmDIV.gr <- tssmp.gr[which(dst<=maxgap & dst > 0)]
  tsspDIV.gr <- tssp_following_tssm.gr[which(dst<=maxgap & dst > 0)]
  if(geneID) divgenes_gap.gr <- GRanges(seqnames=seqnames(tssmDIV.gr),IRanges(start=start(tssmDIV.gr),end=end(tsspDIV.gr)),gIDm=tssmDIV.gr$gene_id,gIDp=tsspDIV.gr$gene_id)
  else divgenes_gap.gr <- GRanges(seqnames=seqnames(tssmDIV.gr),IRanges(start=start(tssmDIV.gr),end=end(tsspDIV.gr)))
  folm <- findOverlaps(tss.gr,tssmDIV.gr)
  folp <- findOverlaps(tss.gr,tsspDIV.gr)
  toDel <- unique(c(queryHits(folm),queryHits(folp)))
  nodivgenes.gr <- tss.gr[-toDel]
  return(list(divgenes=divgenes_gap.gr,nodivgenes=nodivgenes.gr))
}




fill_uncovered_domain_GR <- function(mydom.gr){
  mappable.gr <- import("resources/ext_data/Database/Drosophila/Genome_info/genome/DM3/mappable_regions/iHMM.M1K16.fly_L3.bed")
  mappable.gr <- addSeqinfo(mappable.gr)
  plotlist <- NULL
  for(nme in unique(mappable.gr$name)){
    assign(paste0(nme,".gr"),mappable.gr[mappable.gr$name==nme])
    plotlist <- c(plotlist,paste0(nme,".gr"))
  }
  xtendUnmp.gr <- GenomicRanges::reduce(`17_Unmap.gr`+ 1000)
  gap.gr <- setdiff(dm3_mychr.gr,mydom.gr)
  mygr_filled.gr <- GenomicRanges::reduce(c(mydom.gr,subsetByOverlaps(gap.gr,xtendUnmp.gr,type="within")))
  return(mygr_filled.gr)
}


gap_gr <- function(my.gr,sqinf) {
  mygap.gr <- gaps(my.gr)
  mygap.gr <- subset(mygap.gr,width %ni% seqlengths(sqinf))
  return(mygap.gr)
}



bw_signal <- function(trt.bwp,ctl.bwp,anc.gr) {
  trt.cov <- importData(trt.bwp,"BigWig",anc.gr)[[1]]
  mcols(anc.gr)["trtSIGNAL"] <- rowSums(as.matrix(trt.cov[anc.gr]),na.rm = T)/width(anc.gr)
  ctl.cov <- importData(ctl.bwp,"BigWig",anc.gr)[[1]]
  mcols(anc.gr)["ctlSIGNAL"] <- rowSums(as.matrix(ctl.cov[anc.gr]),na.rm = T)/width(anc.gr)
  return(anc.gr)
}

bw_signal_simple <- function(ctl.bwp,anc.gr,name="") {
  ctl.cov <- importData(ctl.bwp,"BigWig",anc.gr)[[1]]
  mcols(anc.gr)[paste0(name)] <- rowSums(as.matrix(ctl.cov[anc.gr]),na.rm = T)/width(anc.gr)
  return(anc.gr)
}

bw_matrix_simple <- function(sig.bwp,anc.gr,upstream=0, downstream=0) {
  anc.gr <- extendGR(anc.gr,upstream = upstream,downstream = downstream)
  sig.cov <- importData(sig.bwp,"BigWig",anc.gr)[[1]]
  return(as.matrix(sig.cov[anc.gr]))
}

lfc_zs_dif_decile <- function(anc.gr,NTILE,ASYM=T,nm="") {
  anc.gr$LFC <- log2((anc.gr$trtSIGNAL + 1e-5 ) / (anc.gr$ctlSIGNAL + 1e-5))
  anc.gr$ZS <- (anc.gr$trtSIGNAL -  anc.gr$ctlSIGNAL)/ sqrt(mean(c(anc.gr$trtSIGNAL,anc.gr$ctlSIGNAL)))
  anc.gr$DIF <- (anc.gr$trtSIGNAL -  anc.gr$ctlSIGNAL)
  print("ok")
  if(ASYM){
    anc_sup.gr <- subset(anc.gr,LFC >= 0)
    anc_min.gr <- subset(anc.gr,LFC < 0)
    anc_sup.gr <- GRanges(as.data.frame(anc_sup.gr) %>% mutate(decLFC=(NTILE+1)-ntile(LFC,NTILE)))
    anc_min.gr <- GRanges(as.data.frame(anc_min.gr) %>% mutate(decLFC=(NTILE*2+1)-ntile(LFC,NTILE)))
    anc.gr <- c(anc_sup.gr,anc_min.gr)
    
    anc_sup.gr <- subset(anc.gr,ZS >= 0)
    anc_min.gr <- subset(anc.gr,ZS < 0)
    anc_sup.gr <- GRanges(as.data.frame(anc_sup.gr) %>% mutate(decZS=(NTILE+1)-ntile(ZS,NTILE)))
    anc_min.gr <- GRanges(as.data.frame(anc_min.gr) %>% mutate(decZS=(NTILE*2+1)-ntile(ZS,NTILE)))
    anc.gr <- c(anc_sup.gr,anc_min.gr)
    
    anc_sup.gr <- subset(anc.gr,DIF >= 0)
    anc_min.gr <- subset(anc.gr,DIF < 0)
    anc_sup.gr <- GRanges(as.data.frame(anc_sup.gr) %>% mutate(decDIF=(NTILE+1)-ntile(DIF,NTILE)))
    anc_min.gr <- GRanges(as.data.frame(anc_min.gr) %>% mutate(decDIF=(NTILE*2+1)-ntile(DIF,NTILE)))
    anc.gr <- c(anc_sup.gr,anc_min.gr)
    
  }else{
    anc.gr <- GRanges(as.data.frame(anc.gr) %>% mutate(decLFC=(NTILE+1)-ntile(LFC,NTILE)))
    anc.gr <- GRanges(as.data.frame(anc.gr) %>% mutate(decZS=(NTILE+1)-ntile(ZS,NTILE)))
    anc.gr <- GRanges(as.data.frame(anc.gr) %>% mutate(decDIF=(NTILE+1)-ntile(DIF,NTILE)))
  }
  clnm <- colnames(mcols(anc.gr))
  lfc_c <- which(clnm == "LFC");clnm[lfc_c] <- paste0("LFC",nm)
  lfc_c <- which(clnm == "decLFC");clnm[lfc_c] <- paste0("decLFC",nm)
  lfc_c <- which(clnm == "ZS");clnm[lfc_c] <- paste0("ZS",nm)
  lfc_c <- which(clnm == "decZS");clnm[lfc_c] <- paste0("decZS",nm)
  lfc_c <- which(clnm == "DIF");clnm[lfc_c] <- paste0("DIF",nm)
  lfc_c <- which(clnm == "decDIF");clnm[lfc_c] <- paste0("decDIF",nm)
  
  colnames(mcols(anc.gr)) <- clnm
  return(anc.gr)
}


# NA aware
Range2 <- function(vec){
  return(max(vec,na.rm = T)-min(vec,na.rm=T))
}
# Range of a vector
Range <- function(vec){
  return(max(vec)-min(vec))
}

log10(5/0.2)*sqrt(5.2/2)
log10(2500/100)*sqrt(2600/2)





# GRAPHICS ----------------------------------------------------------------
theme_geometry <- function(xvals, yvals, xgeo = 0, ygeo = 0, 
                           color = "black", size = 1, 
                           xlab = "x", ylab = "y",
                           ticks = 10,
                           textsize = 3,
                           xlimit = max(abs(xvals),abs(yvals)),
                           ylimit = max(abs(yvals),abs(xvals)),
                           epsilon = max(xlimit,ylimit)/50){
  
  #INPUT:
  #xvals .- Values of x that will be plotted
  #yvals .- Values of y that will be plotted
  #xgeo  .- x intercept value for y axis
  #ygeo  .- y intercept value for x axis
  #color .- Default color for axis
  #size  .- Line size for axis
  #xlab  .- Label for x axis
  #ylab  .- Label for y axis
  #ticks .- Number of ticks to add to plot in each axis
  #textsize .- Size of text for ticks
  #xlimit .- Limit value for x axis 
  #ylimit .- Limit value for y axis
  #epsilon .- Parameter for small space
  
  
  #Create axis 
  xaxis <- data.frame(x_ax = c(-xlimit, xlimit), y_ax = rep(ygeo,2))
  yaxis <- data.frame(x_ax = rep(xgeo, 2), y_ax = c(-ylimit, ylimit))
  
  #Add axis
  theme.list <- 
    list(
      theme_void(), #Empty the current theme
      geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = xaxis),
      geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = yaxis),
      annotate("text", x = xlimit + 2*epsilon, y = ygeo, label = xlab, size = 2*textsize),
      annotate("text", x = xgeo, y = ylimit + 4*epsilon, label = ylab, size = 2*textsize),
      xlim(-xlimit - 7*epsilon, xlimit + 7*epsilon), #Add limits to make it square
      ylim(-ylimit - 7*epsilon, ylimit + 7*epsilon)  #Add limits to make it square
    )
  
  #Add ticks programatically
  ticks_x <- round(seq(-xlimit, xlimit, length.out = ticks),2)
  ticks_y <- round(seq(-ylimit, ylimit, length.out = ticks),2)
  
  #Add ticks of x axis
  nlist <- length(theme.list)
  for (k in 1:ticks){
    
    #Create data frame for ticks in x axis
    xtick <- data.frame(xt = rep(ticks_x[k], 2), 
                        yt = c(xgeo + epsilon, xgeo - epsilon))
    
    #Create data frame for ticks in y axis
    ytick <- data.frame(xt = c(ygeo + epsilon, ygeo - epsilon), 
                        yt = rep(ticks_y[k], 2))
    
    #Add ticks to geom line for x axis
    theme.list[[nlist + 4*k-3]] <- geom_line(aes(x = xt, y = yt), 
                                             data = xtick, size = size, 
                                             color = color)
    
    #Add labels to the x-ticks
    theme.list[[nlist + 4*k-2]] <- annotate("text", 
                                            x = ticks_x[k], 
                                            y = ygeo - 2.5*epsilon,
                                            size = textsize,
                                            label = paste(ticks_x[k]))
    
    
    #Add ticks to geom line for y axis
    theme.list[[nlist + 4*k-1]] <- geom_line(aes(x = xt, y = yt), 
                                             data = ytick, size = size, 
                                             color = color)
    
    #Add labels to the y-ticks
    theme.list[[nlist + 4*k]] <- annotate("text", 
                                          x = xgeo - 2.5*epsilon, 
                                          y = ticks_y[k],
                                          size = textsize,
                                          label = paste(ticks_y[k]))
  }
  
  #Add theme
  #theme.list[[3]] <- 
  return(theme.list)
}

#   -----------------------------------------------------------------------



ntile_na <- function(x,n)
{
  notna <- !is.na(x)
  out <- rep(NA_real_,length(x))
  out[notna] <- ntile(x[notna],n)
  return(out)
}

# OVERLAPS ---
# input : gr1, gr2
# output : list(gr1xgr2, gr1nogr2)
myOverlaps <- function(n1="n1",n2="n2",gr1.GR,gr2.GR,gap=1,ign.str=T,save=F,outdir=""){
  
  if(ign.str){
    strand(gr1.GR) <- "*"
    strand(gr2.GR) <- "*"
  }
  gr1xgr2.GR <- unique(subsetByOverlaps(gr1.GR, gr2.GR,maxgap=gap))
  gr2xgr1.GR <- unique(subsetByOverlaps( gr2.GR,gr1.GR,maxgap=gap))
  gr1nogr2.GR <- gr1.GR[gr1.GR %ni% gr1xgr2.GR]
  gr2nogr1.GR <- gr2.GR[gr2.GR %ni% gr2xgr1.GR]
  if(save){
    saveRDS(gr1xgr2.GR,paste0(outdir,"/",n1,"_",n2,"_g",gap,"_Kc167_dm3.gr.RData"))
    saveRDS(gr2xgr1.GR,paste0(outdir,"/",n2,"_",n1,"_g",gap,"_Kc167_dm3.gr.RData"))
    saveRDS(gr1nogr2.GR,paste0(outdir,"/",n1,"_no",n2,"_g",gap,"_Kc167_dm3.gr.RData"))
    saveRDS(gr2nogr1.GR,paste0(outdir,"/",n2,"_no",n1,"_g",gap,"_Kc167_dm3.gr.RData"))
    return()
  }
  return(list(gr1xgr2=gr1xgr2.GR,gr2xgr1=gr2xgr1.GR,gr1nogr2=gr1nogr2.GR,gr2nogr1=gr2nogr1.GR))
}


removeOOB <- function(f.gr){
  toDel <- GenomicRanges:::get_out_of_bound_index(f.gr)
  if(length(toDel)) f.gr <- f.gr[-GenomicRanges:::get_out_of_bound_index(f.gr)]
  return(f.gr)
}


BinToDec <- function(x) 
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))

bcast <- function(... , agg.name =NA){ res <- cast(...)
if(!is.na(agg.name)){ names(res)[length(res)] <- agg.name } 
res}

mutate_cond <- function(.data, condition, ..., new_init = NA, envir = parent.frame()) {
  # Initialize any new variables as new_init
  new_vars <- substitute(list(...))[-1]
  new_vars %<>% sapply(deparse) %>% names %>% setdiff(names(.data))
  .data[, new_vars] <- new_init
  
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data %>% filter(condition) %>% mutate(...)
  .data
}

mutate_cond_ <- function(.data, condition, ..., new_init = NA, envir = parent.frame()) {
  # Initialize any new variables as new_init
  new_vars <- substitute(list(...))[-1]
  new_vars %<>% sapply(deparse) %>% names %>% setdiff(names(.data))
  .data[, new_vars] <- new_init
  
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}


domainsBorder <- function(tad.gr) {
  bordL.gr <- resize(tad.gr,1,"start")
  bordR.gr <- resize(tad.gr,1,"end")
  bord.gr <- reduce(c(bordL.gr,bordR.gr))
  
}


mrc_fisherTest <- function(gr.cl,xtend,cellT,xtBool) {
  grlf <- list()
  for(gr in gr.cl[xtBool]){
    c.gr <- get(gr)
    mcols(c.gr) <- NULL
    gn <- strsplit(gr,"\\.")[[1]][2]
    nme <- strsplit(gr,"\\.")[[1]][1]
    c.gr <- addSeqinfo(resize(c.gr,1,"center")+ xtend,gn,cellT)
    strand(c.gr) <- "*"
    grlf[[nme]] <- c.gr
  }
  for(gr in gr.cl[!xtBool]){
    c.gr <- get(gr)
    mcols(c.gr) <- NULL
    gn <- strsplit(gr,"\\.")[[1]][2]
    nme <- strsplit(gr,"\\.")[[1]][1]
    c.gr <- addSeqinfo(c.gr,gn,cellT)
    strand(c.gr) <- "*"
    grlf[[nme]] <- c.gr
  }
  set.gr <- NULL
  for(nm in names(grlf)){
    set.gr <- c(set.gr,grlf[[nm]])
  }
  
  set.gr <- GenomicRanges::trim(Reduce("c",set.gr))
  mrc.df <- data.frame(idx=1:length(set.gr))
  print(nrow(mrc.df))
  print("COOL")
  shl <- NULL
  for(nm in names(grlf)){
    asgn <- paste0("sh",nm)
    # assign(paste0("sh",nm),subjectHits(findOverlaps(grlf[[nm]],set.gr,type="within")),envir = .GlobalEnv)
    assign(paste0("sh",nm),findOverlaps(set.gr,grlf[[nm]])@from,envir = .GlobalEnv)
    shl[[nm]] <- get(asgn)
    print(nm)
    print(length(unique(get(asgn))))
  }
  mrc.df$all <- 1
  for(nm in names(grlf)){
    mrc.df[,nm] <- 0
    mrc.df[shl[[nm]],nm] <- 1
  }
  
  nmB <- NULL
  for(nm in colnames(mrc.df)[-1]){
    assn <- toupper(nm)
    assign(assn,mrc.df[,nm]==1,envir = .GlobalEnv)
    nmB[[assn]] <- get(assn)
  }
  return(list(mrc=mrc.df,BOOL=nmB))
}





# EXAMPLE
# variable must have the following shape NAME.GENOME.WHATEVER...
# If some features MUST NOT be resized set the boolean vector xtBool to F at appropriate places
# myL <- c("z1.hg19.gr","mtr4.hg19.gr","ctcf.hg19.gr","k4me1.hg19.gr","prompts.hg19.500.gr")
# tst <- mrc_fisherTest(myL,500,"hela",c(T,T,T,T,F))

generateVennTHREE <- function(myWhole.l,xtVec,xtend=500,myVenn.vec,f_dir="/home/",fisher=T,GP=NULL,celltype="hela") {
  F1 <- toupper(myVenn.vec[1])
  noF1 <- paste0("no",F1)
  F2 <- toupper(myVenn.vec[2])
  F3 <- toupper(myVenn.vec[3])
  
  toMRC <- mrc_fisherTest(myWhole.l,xtend,celltype,xtVec)
  mrc <- toMRC$mrc
  # browser()
  print("tot")
  
  # Compute Fisher test for the given SET
  for(nm in names(toMRC$BOOL)){
    assign(nm,toMRC[["BOOL"]][[nm]],envir = .GlobalEnv)
    assign(paste0("no",nm),!get(nm),envir = .GlobalEnv)
  }
  
  if(fisher){
    myObjP <- toupper(myVenn.vec[-1])
    myObjN <- paste0("no",toupper(myVenn.vec[-1]))
    expG <- expand.grid(0:1,  0:1)
    expG <- expG[-1,]
    myLst <- t(apply(expG==1,1,function(x)c(myObjP[x],myObjN[!x])))
    fd_dir <- create(paste0(f_dir,"/details/SET_",glue::glue_collapse(c(F1,F2,F3),"_"),"/"))
    myDF.l <- apply(myLst,1,
                    function(x){
                      print(x)
                      ok <- isEnrichAny(x,"ALL",F1,fd_dir,mrc)
                    })
    myDF.l <- apply(myLst,1,
                    function(x){
                      print(x)
                      ok <-  isEnrichAny(x,"ALL",noF1,fd_dir,mrc)
                    })
    
    F2_noF1 <- c(F2,noF1)
    F3_noF1 <- c(F3,noF1)
    print(F2_noF1)
    print(F3_noF1)
    ok <- isEnrichAny(F2_noF1,"ALL",F3_noF1,fd_dir,mrc)
  }
  # browser()
  
  
  # COMPUTE PROPORTIONNAL VENN DIAGRAM
  colnames(mrc) <- toupper(colnames(mrc))
  dt_mrc <- data.table::as.data.table(mrc)
  n111 <- dt_mrc[get(F1)==1L & get(F2)==1L & get(F3)==1L, .N]
  n110 <- dt_mrc[get(F1)==1L & get(F2)==1L & get(F3)==0L, .N]
  n101 <- dt_mrc[get(F1)==1L & get(F2)==0L & get(F3)==1L, .N]
  n100 <- dt_mrc[get(F1)==1L & get(F2)==0L & get(F3)==0L, .N]
  n011 <- dt_mrc[get(F1)==0L & get(F2)==1L & get(F3)==1L, .N]
  n010 <- dt_mrc[get(F1)==0L & get(F2)==1L & get(F3)==0L, .N]
  n001 <- dt_mrc[get(F1)==0L & get(F2)==0L & get(F3)==1L, .N]
  myVenn <- Venn(SetNames=c(F1,F2,F3),
                 Weight=c(`111`=n111,`110`=n110,`101`=n101,`100`=n100,`011`=n011,`010`=n010,`001`=n001))
  
  v_dir <- create(paste0(f_dir,"/plot_venn/SET_",F1,"_",F2,"_",F3,"/"))
  tiff(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2,F3),"_"),"_Figure.tiff"),width=2000, height=2000, res=2000)
  plot(myVenn,doWeights=T,show=list(FaceText="",SetLabels=F),gp=GP)
  dev.off()
  pdf(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2,F3),"_"),"_View.pdf"),width=40, height=40)
  plot(myVenn,doWeights=T,gp=GP)
  dev.off()
  pdf(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2,F3),"_"),"_Figure.pdf"),width=40, height=40)
  plot(myVenn,doWeights=T,show=list(FaceText="",SetLabels=F),gp=GP)
  dev.off()
  jpeg(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2,F3),"_"),"_View.jpeg"),width=2000, height=2000)
  plot(myVenn,doWeights=T,gp=GP)
  dev.off()
  # rm(list=ls(),envir=test.env)
  return(myVenn)
}


generateVennTHREE_NEW <- function(myWhole.l,xtVec,xtend=500,myVenn.vec,f_dir="/home/",fisher=F,genome="hg19",celltype="hela") {
  F1 <- toupper(myVenn.vec[1])
  noF1 <- paste0("no",F1)
  F2 <- toupper(myVenn.vec[2])
  F3 <- toupper(myVenn.vec[3])
  
  # toMRC <- mrc_fisherTest(myWhole.l,xtend,celltype,xtVec)
  # mrc <- toMRC$mrc
  
  all.gr <- trim(GenomicRanges::reduce(do.call("c",sapply(1:length(myWhole.l),function(X){
    gr <- get(myWhole.l[X])
    mcols(gr) <- NULL
    if(xtVec[X]){print("REDUCE"); ret_res <- addSeqinfo(resize(gr,1,"center")+xtend,genome,celltype);strand(ret_res) <- "*"}
    else{ret_res <- addSeqinfo(gr,genome,celltype);strand(ret_res) <- "*"}
    return(ret_res)}))))
  print("ok")
  
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
  
  mrc <- apa_DF
  
  # Compute Fisher test for the given SET
  for(nm in colnames(mrc)[-1]){
    assign(nm,mrc[,nm]==1,envir = .GlobalEnv)
    assign(paste0("no",nm),!get(nm),envir = .GlobalEnv)
  }
  
  if(fisher){
    myObjP <- toupper(myVenn.vec[-1])
    myObjN <- paste0("no",toupper(myVenn.vec[-1]))
    expG <- expand.grid(0:1,  0:1)
    expG <- expG[-1,]
    myLst <- t(apply(expG==1,1,function(x)c(myObjP[x],myObjN[!x])))
    fd_dir <- create(paste0(f_dir,"/details/SET_",glue::glue_collapse(c(F1,F2,F3),"_"),"/"))
    
    myDF.l <- apply(myLst,1,
                    function(x){
                      print(x)
                      ok <- isEnrichAny(x,"ALL",F1,fd_dir,mrc)
                    })
    myDF.l <- apply(myLst,1,
                    function(x){
                      print(x)
                      ok <-  isEnrichAny(x,"ALL",noF1,fd_dir,mrc)
                    })
    
    F2_noF1 <- c(F2,noF1)
    F3_noF1 <- c(F3,noF1)
    print(F2_noF1)
    print(F3_noF1)
    ok <- isEnrichAny(F2_noF1,"ALL",F3_noF1,fd_dir,mrc)
  }
  
  # COMPUTE PROPORTIONNAL VENN DIAGRAM
  colnames(mrc) <- toupper(colnames(mrc))
  dt_mrc <- data.table::as.data.table(mrc)
  n111 <- dt_mrc[get(F1)==1L & get(F2)==1L & get(F3)==1L, .N]
  n110 <- dt_mrc[get(F1)==1L & get(F2)==1L & get(F3)==0L, .N]
  n101 <- dt_mrc[get(F1)==1L & get(F2)==0L & get(F3)==1L, .N]
  n100 <- dt_mrc[get(F1)==1L & get(F2)==0L & get(F3)==0L, .N]
  n011 <- dt_mrc[get(F1)==0L & get(F2)==1L & get(F3)==1L, .N]
  n010 <- dt_mrc[get(F1)==0L & get(F2)==1L & get(F3)==0L, .N]
  n001 <- dt_mrc[get(F1)==0L & get(F2)==0L & get(F3)==1L, .N]
  myVenn <- Venn(SetNames=c(F1,F2,F3),
                 Weight=c(`111`=n111,`110`=n110,`101`=n101,`100`=n100,`011`=n011,`010`=n010,`001`=n001))
  
  v_dir <- create(paste0(f_dir,"/plot_venn/SET_",F1,"_",F2,"_",F3,"/"))
  print(v_dir)
  pdf(paste0(v_dir,"Proportionnal_Venn_setO_",glue::glue_collapse(c(F1,F2,F3),"_"),"_View.pdf"),width=40, height=40)
  print(plot(myVenn,doWeights=T))
  dev.off()
  pdf(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2,F3),"_"),"_Figure.pdf"),width=40, height=40)
  print(plot(myVenn,doWeights=T,show=list(FaceText="",SetLabels=F)))
  dev.off()
  return(myVenn)
}

generateVennTHREE_NEW_G2I <- function(myWhole.l,xtVec,xtend=500,myVenn.vec,f_dir="/home/",fisher=F,genome="hg19",celltype="hela") {
  F1 <- toupper(myVenn.vec[1])
  noF1 <- paste0("no",F1)
  F2 <- toupper(myVenn.vec[2])
  F3 <- toupper(myVenn.vec[3])
  
  # toMRC <- mrc_fisherTest(myWhole.l,xtend,celltype,xtVec)
  # mrc <- toMRC$mrc
  
  all.gr <- trim(GenomicRanges::reduce(do.call("c",sapply(1:length(myWhole.l),function(X){
    gr <- get(myWhole.l[X])
    mcols(gr) <- NULL
    if(xtVec[X]){print("REDUCE"); ret_res <- addSeqinfo(resize(gr,1,"center")+xtend,genome,celltype);strand(ret_res) <- "*"}
    else{ret_res <- addSeqinfo(gr,genome,celltype);strand(ret_res) <- "*"}
    return(ret_res)}))))
  print("ok")
  
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
  
  mrc <- apa_DF
  
  # Compute Fisher test for the given SET
  for(nm in colnames(mrc)[-1]){
    assign(nm,mrc[,nm]==1,envir = .GlobalEnv)
    assign(paste0("no",nm),!get(nm),envir = .GlobalEnv)
  }
  
  if(fisher){
    myObjP <- toupper(myVenn.vec[-1])
    myObjN <- paste0("no",toupper(myVenn.vec[-1]))
    expG <- expand.grid(0:1,  0:1)
    expG <- expG[-1,]
    myLst <- t(apply(expG==1,1,function(x)c(myObjP[x],myObjN[!x])))
    fd_dir <- create(paste0(f_dir,"/details/SET_",glue::glue_collapse(c(F1,F2,F3),"_"),"/"))
    
    myDF.l <- apply(myLst,1,
                    function(x){
                      print(x)
                      ok <- isEnrichAny(x,"ALL",F1,fd_dir,mrc)
                    })
    myDF.l <- apply(myLst,1,
                    function(x){
                      print(x)
                      ok <-  isEnrichAny(x,"ALL",noF1,fd_dir,mrc)
                    })
    
    F2_noF1 <- c(F2,noF1)
    F3_noF1 <- c(F3,noF1)
    print(F2_noF1)
    print(F3_noF1)
    ok <- isEnrichAny(F2_noF1,"ALL",F3_noF1,fd_dir,mrc)
  }
  
  # COMPUTE PROPORTIONNAL VENN DIAGRAM
  colnames(mrc) <- toupper(colnames(mrc))
  dt_mrc <- data.table::as.data.table(mrc)
  n111 <- dt_mrc[get(F1)==1L & get(F2)==1L & get(F3)==1L, .N]
  n110 <- dt_mrc[get(F1)==1L & get(F2)==1L & get(F3)==0L, .N]
  n101 <- dt_mrc[get(F1)==1L & get(F2)==0L & get(F3)==1L, .N]
  n100 <- dt_mrc[get(F1)==1L & get(F2)==0L & get(F3)==0L, .N]
  n011 <- dt_mrc[get(F1)==0L & get(F2)==1L & get(F3)==1L, .N]
  n010 <- dt_mrc[get(F1)==0L & get(F2)==1L & get(F3)==0L, .N]
  n001 <- dt_mrc[get(F1)==0L & get(F2)==0L & get(F3)==1L, .N]
  myVenn <- Venn(SetNames=c(F1,F2,F3),
                 Weight=c(`111`=n111,`110`=n110,`101`=n101,`100`=n100,`011`=n011,`010`=n010,`001`=n001))
  
  v_dir <- paste0(f_dir,".VENN.SET_",F1,"_",F2,"_",F3)
  png(paste0(v_dir,".View.png"),width=1800, height=1200)
  plot(myVenn,doWeights=T)
  dev.off()
  png(paste0(v_dir,".Blank.png"),width=1800, height=1200)
  plot(myVenn,doWeights=T,show=list(FaceText="",SetLabels=F))
  dev.off()
  return(myVenn)
}


VennFisherThreeGrp <- function(baseSET.chgrv,tiles.gr,xtVec.v,myVenn.vec,f_dir="/home/",fisher=F,genome="hg19",celltype="hela") {
  F1 <- toupper(myVenn.vec[1])
  noF1 <- paste0("no",F1)
  F2 <- toupper(myVenn.vec[2])
  F3 <- toupper(myVenn.vec[3])
  
  base.df <- data.frame(IDX=1:length(tiles.gr),ALL=1)
  for(X in 1:length(baseSET.chgrv)){
    feat.nm <- baseSET.chgrv[X]
    feat.gr <- get(feat.nm)
    mcols(feat.gr) <- NULL;strand(feat.gr) <- "*"
    nm <- toupper(gsub(paste0("(.*).",genome,".gr"),"\\1",feat.nm))
    # feat.gr <- addSeqinfo(resize(feat.gr,1,"center"),genome,celltype)
    feat.gr <- addSeqinfo(feat.gr,genome,celltype)
    featTiles.IDX <- unique(findOverlaps(tiles.gr,feat.gr,ignore.strand=T)@from)
    base.df[,nm] <- 0
    base.df[featTiles.IDX,nm] <- 1
  }
  
  mrc <- base.df  
  
  # Compute Fisher test for the given SET
  for(nm in colnames(mrc)[-1]){
    assign(nm,mrc[,nm]==1,envir = .GlobalEnv)
    assign(paste0("no",nm),!get(nm),envir = .GlobalEnv)
  }
  
  
  if(fisher){
    mrcF <- mrc[,-c(1,2)] # remove ALL column
    toKEEP <- which(mrcF==1)
    mrcF <- mrcF[toKEEP,]
    myObjP <- toupper(myVenn.vec[-1])
    myObjN <- paste0("no",toupper(myVenn.vec[-1]))
    expG <- expand.grid(0:1,  0:1)
    expG <- expG[-1,]
    myLst <- t(apply(expG==1,1,function(x)c(myObjP[x],myObjN[!x])))
    fd_dir <- create(paste0(f_dir,"/details/SET_",glue::glue_collapse(c(F1,F2,F3),"_"),"/"))
    myDF.l <- apply(myLst,1,
                    function(x){
                      print(x)
                      ok <- isEnrichAny(x,"ALL",F1,fd_dir,mrcF)
                    })
    myDF.l <- apply(myLst,1,
                    function(x){
                      print(x)
                      ok <-  isEnrichAny(x,"ALL",noF1,fd_dir,mrcF)
                    })
    
    F2_noF1 <- c(F2,noF1)
    F3_noF1 <- c(F3,noF1)
    print(F2_noF1)
    print(F3_noF1)
    ok <- isEnrichAny(F2_noF1,"ALL",F3_noF1,fd_dir,mrcF)
  }
  
  
  
  # COMPUTE PROPORTIONNAL VENN DIAGRAM
  colnames(mrc) <- toupper(colnames(mrc))
  dt_mrc <- data.table::as.data.table(mrc)
  n111 <- dt_mrc[get(F1)==1L & get(F2)==1L & get(F3)==1L, .N]
  n110 <- dt_mrc[get(F1)==1L & get(F2)==1L & get(F3)==0L, .N]
  n101 <- dt_mrc[get(F1)==1L & get(F2)==0L & get(F3)==1L, .N]
  n100 <- dt_mrc[get(F1)==1L & get(F2)==0L & get(F3)==0L, .N]
  n011 <- dt_mrc[get(F1)==0L & get(F2)==1L & get(F3)==1L, .N]
  n010 <- dt_mrc[get(F1)==0L & get(F2)==1L & get(F3)==0L, .N]
  n001 <- dt_mrc[get(F1)==0L & get(F2)==0L & get(F3)==1L, .N]
  myVenn <- Venn(SetNames=c(F1,F2,F3),
                 Weight=c(`111`=n111,`110`=n110,`101`=n101,`100`=n100,`011`=n011,`010`=n010,`001`=n001))
  
  
  v_dir <- create(paste0(f_dir,"/plot_venn/SET_",F1,"_",F2,"_",F3,"/"))
  print(v_dir)
  pdf(paste0(v_dir,"Proportionnal_Venn_setO_",glue::glue_collapse(c(F1,F2,F3),"_"),"_View.pdf"),width=40, height=40)
  print(plot(myVenn,doWeights=T,gp=color_Venn))
  dev.off()
  pdf(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2,F3),"_"),"_Figure.pdf"),width=40, height=40)
  print(plot(myVenn,doWeights=T,show=list(FaceText="",SetLabels=F)))
  dev.off()
  
  
} # END 

generateVennTWO <- function(myWhole.l,xtVec,xtend=500,myVenn.vec,f_dir="/home/",GP=NULL,celltype="dm3") {
  F1 <- toupper(myVenn.vec[1])
  F2 <- toupper(myVenn.vec[2])
  toMRC <- mrc_fisherTest(myWhole.l,xtend,celltype,xtVec)
  mrc <- toMRC$mrc
  
  
  # COMPUTE PROPORTIONNAL VENN DIAGRAM
  colnames(mrc) <- toupper(colnames(mrc))
  dt_mrc <- data.table::as.data.table(mrc)
  n11 <- dt_mrc[get(F1)==1L & get(F2)==1L, .N]
  n10 <- dt_mrc[get(F1)==1L & get(F2)==0L, .N]
  n01 <- dt_mrc[get(F1)==0L & get(F2)==1L, .N]
  
  myVenn <- Venn(SetNames=c(F1,F2),
                 Weight=c(`11`=n11,`10`=n10,`01`=n01))
  
  
  print("toto")
  v_dir <- create(paste0(f_dir,"/plot_venn/SET_",F1,"_",F2,"/"))
  tiff(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2),"_"),"_Figure.tiff"),width=2000, height=2000, res=2000)
  plot(myVenn,doWeights=T,show=list(FaceText="",SetLabels=F),gp=GP)
  dev.off()
  pdf(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2),"_"),"_View.pdf"),width=40, height=40)
  plot(myVenn,doWeights=T,gp=GP)
  dev.off()
  pdf(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2),"_"),"_Figure.pdf"),width=40, height=40)
  plot(myVenn,doWeights=T,show=list(FaceText="",SetLabels=F),gp=GP)
  dev.off()
  jpeg(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2),"_"),"_View.jpeg"),width=2000, height=2000)
  plot(myVenn,doWeights=T,gp=GP)
  dev.off()
  return(myVenn)
}

generateVennTWO_NEW_G2I <- function(myWhole.l,xtVec,xtend=500,myVenn.vec,f_dir="/home/",GP=NULL,fisher = F,genome="hg19",celltype="hela") {
  F1 <- toupper(myVenn.vec[1])
  F2 <- toupper(myVenn.vec[2])
  
  all.gr <- trim(GenomicRanges::reduce(do.call("c",sapply(1:length(myWhole.l),function(X){
    gr <- get(myWhole.l[X])
    mcols(gr) <- NULL
    print(gr)
    if(xtVec[X]){print("REDUCE"); ret_res <- addSeqinfo(resize(gr,1,"center")+xtend,genome,celltype);strand(ret_res) <- "*"}
    else{ret_res <- addSeqinfo(gr,genome,celltype);strand(ret_res) <- "*"}
    return(ret_res)}))))
  
  print(all.gr)
  
  apa_DF <- data.frame(IDX=1:length(all.gr),ALL=1)
  for(feat in myWhole.l){
    nm <- toupper(gsub(paste0("(.*).",genome,".gr"),"\\1",feat))
    print(nm)
    apa_DF[,nm] <- 0
    gr <- get(feat)
    mcols(gr) <- NULL
    gr <- addSeqinfo(gr,genome,celltype)
    strand(gr) <- "*"
    IDX <- unique(findOverlaps(gr,all.gr)@to)
    apa_DF[IDX,nm] <- 1
  }
  
  mrc <- apa_DF
  
  
  for(nm in colnames(mrc)[-1]){
    assign(nm,mrc[,nm]==1,envir = .GlobalEnv)
    assign(paste0("no",nm),!get(nm),envir = .GlobalEnv)
  }
  
  if(fisher){
    fish_dir <- paste0(f_dir,".FISHER.SET_",F1,"_",F2)
    isEnrichAny_G2I(F2,"ALL",F1,out_d = fish_dir,mrc)
  }
  
  
  # COMPUTE PROPORTIONNAL VENN DIAGRAM
  colnames(mrc) <- toupper(colnames(mrc))
  print(colnames(mrc))
  dt_mrc <- data.table::as.data.table(mrc)
  n11 <- dt_mrc[get(F1)==1L & get(F2)==1L, .N]
  n10 <- dt_mrc[get(F1)==1L & get(F2)==0L, .N]
  n01 <- dt_mrc[get(F1)==0L & get(F2)==1L, .N]
  
  myVenn <- Venn(SetNames=c(F1,F2),
                 Weight=c(`11`=n11,`10`=n10,`01`=n01))
  
  
  print("toto")
  v_dir <- paste0(f_dir,".VENN.SET_",F1,"_",F2)
  png(paste0(v_dir,".View.png"),width=1800, height=1200)
  plot(myVenn,doWeights=T,gp=GP)
  dev.off()
  png(paste0(v_dir,".Blank.png"),width=1800, height=1200)
  plot(myVenn,doWeights=T,show=list(FaceText="",SetLabels=F),gp=GP)
  dev.off()
  return(myVenn)
}


generateVennTWO_NEW <- function(myWhole.l,xtVec,xtend=500,myVenn.vec,f_dir="/home/",GP=NULL,fisher = F,genome="hg19",celltype="hela") {
  F1 <- toupper(myVenn.vec[1])
  F2 <- toupper(myVenn.vec[2])
  
  
  all.gr <- trim(GenomicRanges::reduce(do.call("c",sapply(1:length(myWhole.l),function(X){
    gr <- get(myWhole.l[X])
    mcols(gr) <- NULL
    print(gr)
    if(xtVec[X]){print("REDUCE"); ret_res <- addSeqinfo(resize(gr,1,"center")+xtend,genome,celltype);strand(ret_res) <- "*"}
    else{ret_res <- addSeqinfo(gr,genome,celltype);strand(ret_res) <- "*"}
    return(ret_res)}))))
  
  print(all.gr)
  apa_DF <- data.frame(IDX=1:length(all.gr))
  for(feat in myWhole.l){
    nm <- toupper(gsub(paste0("(.*).",genome,".gr"),"\\1",feat))
    print(nm)
    apa_DF[,nm] <- 0
    gr <- get(feat)
    mcols(gr) <- NULL
    gr <- addSeqinfo(gr,genome,celltype)
    strand(gr) <- "*"
    IDX <- unique(findOverlaps(gr,all.gr)@to)
    apa_DF[IDX,nm] <- 1
  }
  
  mrc <- apa_DF
  
  
  
  for(nm in colnames(mrc)[-1]){
    assign(nm,mrc[,nm]==1,envir = .GlobalEnv)
    assign(paste0("no",nm),!get(nm),envir = .GlobalEnv)
  }
  
  if(fisher){
    fd_dir <- create(paste0(f_dir,"/details/SET_",glue::glue_collapse(c(F1,F2),"_"),"/"))
    isEnrichAny(F2,"ALL",F1,fd_dir,mrc)
  }
  
  
  # COMPUTE PROPORTIONNAL VENN DIAGRAM
  colnames(mrc) <- toupper(colnames(mrc))
  print(colnames(mrc))
  dt_mrc <- data.table::as.data.table(mrc)
  n11 <- dt_mrc[get(F1)==1L & get(F2)==1L, .N]
  n10 <- dt_mrc[get(F1)==1L & get(F2)==0L, .N]
  n01 <- dt_mrc[get(F1)==0L & get(F2)==1L, .N]
  
  myVenn <- Venn(SetNames=c(F1,F2),
                 Weight=c(`11`=n11,`10`=n10,`01`=n01))
  
  
  print("toto")
  v_dir <- create(paste0(f_dir,"/plot_venn/SET_",F1,"_",F2,"/"))
  tiff(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2),"_"),"_Figure.tiff"),width=2000, height=2000, res=2000)
  plot(myVenn,doWeights=T,show=list(FaceText="",SetLabels=F),gp=GP)
  dev.off()
  pdf(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2),"_"),"_View.pdf"),width=40, height=40)
  plot(myVenn,doWeights=T,gp=GP)
  dev.off()
  pdf(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2),"_"),"_Figure.pdf"),width=40, height=40)
  plot(myVenn,doWeights=T,show=list(FaceText="",SetLabels=F),gp=GP)
  dev.off()
  jpeg(paste0(v_dir,"Proportionnal_Venn_set_",glue::glue_collapse(c(F1,F2),"_"),"_View.jpeg"),width=2000, height=2000)
  plot(myVenn,doWeights=T,gp=GP)
  dev.off()
  return(myVenn)
}





fisher_ntiles <- function(DF=NULL,COLN=NULL,ROWN=NULL,NTILE=4,NAME=NULL) {
  fish.df <- NULL
  
  for(quart in 1:NTILE){
    myDF <- as.data.frame(DF %>%
                            group_by(UQ(as.name(COLN)),UQ(as.name(ROWN)))  %>% 
                            summarize(tot=n()))
    
    pq <- as.numeric(myDF %>% filter(UQ(as.name(COLN))==quart,UQ(as.name(ROWN))==1) %>% 
                       select_("tot") %>% summarize_all(sum))
    pnq <- as.numeric(myDF %>% filter(!(UQ(as.name(COLN))==quart ),UQ(as.name(ROWN))==1) %>% 
                        select_("tot") %>% summarize_all(sum))
    npq <- as.numeric(myDF %>% filter(UQ(as.name(COLN))==quart ,UQ(as.name(ROWN))==0) %>% 
                        select_("tot") %>% summarize_all(sum))
    npnq <- as.numeric(myDF %>% filter(!(UQ(as.name(COLN))==quart),UQ(as.name(ROWN))==0) %>% 
                         select_("tot") %>% summarize_all(sum))
    
    f_mat <- matrix(c(L_C=pq,L_nC=pnq,C_nL=npq,nC_nL=npnq),nrow=2,ncol=2)
    ft <- fisher.test(f_mat)
    
    mypv <- formatC(as.numeric(ft$p.value, format = "e", digits = 2))
    lfc <- round(log2(ft$estimate),2)
    sign <- ifelse(ft$p.value < 0.05,ifelse(ft$p.value < 0.01,ifelse(ft$p.value < 0.001,"***","**"),"*"),"NS")
    fish.df <- rbind.data.frame(fish.df,cbind.data.frame(XNAME=as.character(quart),pv=mypv,lfcb=lfc,sign=sign,YNAME=NAME))
  }
  fish.df$intSign <- sign(fish.df$lfcb)
  fish.df <- fish.df %>% mutate(lfc=ifelse(abs(lfcb) > 2,intSign*2,intSign*abs(lfcb)))
  return(fish.df)
  
}



fisher_classic <- function(DF=NULL,COLN=NULL,ROWN=NULL,XNAME=NULL,YNAME=NULL) {
  fish.df <- NULL
  
  myDF <- as.data.frame(DF %>%
                          group_by(UQ(as.name(COLN)),UQ(as.name(ROWN)))  %>% 
                          summarize(tot=n()))
  
  pq <- as.numeric(myDF %>% filter(UQ(as.name(COLN))==1,UQ(as.name(ROWN))==1) %>% 
                     select_("tot") %>% summarize_all(sum))
  pnq <- as.numeric(myDF %>% filter(UQ(as.name(COLN))==0 ,UQ(as.name(ROWN))==1) %>% 
                      select_("tot") %>% summarize_all(sum))
  npq <- as.numeric(myDF %>% filter(UQ(as.name(COLN))==1 ,UQ(as.name(ROWN))==0) %>% 
                      select_("tot") %>% summarize_all(sum))
  npnq <- as.numeric(myDF %>% filter(UQ(as.name(COLN))==0,UQ(as.name(ROWN))==0) %>% 
                       select_("tot") %>% summarize_all(sum))
  
  f_mat <- matrix(c(L_C=pq,L_nC=pnq,C_nL=npq,nC_nL=npnq),nrow=2,ncol=2)
  ft <- fisher.test(f_mat)
  
  mypv <- formatC(as.numeric(ft$p.value, format = "e", digits = 2))
  lfc <- round(log2(ft$estimate),2)
  sign <- ifelse(ft$p.value < 0.05,ifelse(ft$p.value < 0.01,ifelse(ft$p.value < 0.001,"***","**"),"*"),"NS")
  fish.df <- rbind.data.frame(fish.df,cbind.data.frame(YNAME=YNAME,pv=mypv,lfcb=lfc,sign=sign,XNAME=XNAME))
  fish.df$intSign <- sign(fish.df$lfcb)
  fish.df <- fish.df %>% mutate(lfc=ifelse(abs(lfcb) > 2,intSign*2,intSign*abs(lfcb)))
  return(fish.df)
}



theme_vertical <- theme(axis.text.x=element_text(angle = 90))
# 
# 
# fisher_ntiles_plot <- function(F.df) {
#   colfunc <- colorRampPalette(c("firebrick1", "white","dodgerblue"))
#   
#   p <- ggplot(F.df, aes(x = XNAME, y = YNAME,fill=lfc)) +      geom_tile(col="black") +
#     theme_minimal()+
#     coord_equal(ratio=.15) +
#     scale_y_discrete(expand=c(0,0))+
#     scale_x_discrete(position = "top") +
#     theme(axis.text.x=element_text(size=25, vjust=1, hjust=.5,face = "bold",
#                                    margin=margin(1,0,0,0)),
#           axis.text.y=element_text(size=20, margin=margin(0,-3,0,0),face = "bold"),
#           axis.title.x = element_blank(),
#           legend.title = element_text(size=12, vjust=-2, hjust=0.1,face = "bold"),
#           legend.background = element_rect(color="grey",size=0.8, linetype="dashed"),
#           axis.title.y = element_blank(),
#           legend.title.align=0.1,
#           plot.margin = unit(c(0,1,0,1), "cm"),
#           panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank()) +
#     
#     geom_text(aes(label =sign),size=4) + 
#     scale_fill_gradientn(colours = colfunc(3) ,
#                          guide = guide_colorbar(barwidth = 0.8,
#                                                 title = 'Log2 \nOddsRatio',
#                                                 label.theme = element_text(size=10, vjust=1, hjust=.5,face = "bold",angle = 0),
#                                                 barheight = 10,
#                                                 nbin = 10,
#                                                 draw.ulim = FALSE, 
#                                                 draw.llim = FALSE,
#                                                 ticks = FALSE),
#                          breaks=c(-2,0.2,2),
#                          labels=c("-1.5","0","1.5"),
#                          limits=c(-2,2)) 
#   return(p)
# }



fisher_ntiles_plot <- function(F.df,X,Y,ratio=.75) {
  colfunc <- colorRampPalette(c("firebrick1", "white","dodgerblue"))
  
  p <- ggplot(F.df, aes(x = X, y = Y,fill=lfc)) +      geom_tile(col="black") +
    theme_minimal()+
    coord_equal(ratio=ratio) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top") +
    theme(axis.text.x=element_text(angle = 90,size=25, vjust=1, hjust=.5,face = "bold",
                                   margin=margin(1,0,0,0)),
          axis.text.y=element_text(size=20, margin=margin(0,-3,0,0),face = "bold"),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12, vjust=-2, hjust=0.1,face = "bold"),
          legend.background = element_rect(color="grey",size=0.8, linetype="dashed"),
          axis.title.y = element_blank(),
          legend.title.align=0.1,
          plot.margin = unit(c(0,1,0,1), "cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    
    geom_text(aes(label =sign),size=4) + 
    scale_fill_gradientn(colours = colfunc(3) ,
                         guide = guide_colorbar(barwidth = 0.8,
                                                title = 'Log2 \nOddsRatio',
                                                label.theme = element_text(size=10, vjust=1, hjust=.5,face = "bold",angle = 0),
                                                barheight = 10,
                                                nbin = 10,
                                                draw.ulim = FALSE, 
                                                draw.llim = FALSE,
                                                ticks = FALSE),
                         breaks=c(-2,0.2,2),
                         labels=c("-1.5","0","1.5"),
                         limits=c(-2,2)) 
  return(p)
}




boxplot_ntiles <- function(DF=NULL,COLN=NULL,ROWN=NULL,NTILE=4) {
  bp.df <- NULL
  COLQN <- paste0("QUART_",COLN)
  
  for(quart in 1:NTILE){
    p.vec <- (DF %>% filter(UQ(as.name(COLQN))==quart & UQ(as.name(ROWN))==1))[,COLN]
    np.vec <- (DF %>% filter(UQ(as.name(COLQN))==quart & UQ(as.name(ROWN))==0))[,COLN]
    bp.df <- rbind.data.frame(bp.df,
                              cbind.data.frame(VAL=p.vec,NAME=ROWN,QUART=quart),
                              cbind.data.frame(VAL=np.vec,NAME=paste0("NO_",ROWN),QUART=quart)
    )
  }
  
  return(bp.df)
  
}




tile_plot <-function(f_df){
  p <- ggplot(f_df, aes(x = col, y = line)) +
    geom_tile(fill="black",col="white") +
    # ggtitle(title) +
    theme_minimal() +
    # coord_equal(ratio=0.15) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top") +
    labs(x="",y="") +
    theme(axis.text.x=element_text(size=6, vjust=1, hjust=.5,
                                   margin=margin(1,0,0,0)),
          axis.text.y=element_text(size=6, margin=margin(0,-3,0,0)),
          plot.margin=unit(c(0,0,2,1),"lines"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    geom_text(aes(label = perc),col="white",size=5); p }

flip_plot <- function(df){
  
  hp <- ggplot(df,aes(x=type,y=perc)) + geom_bar(position=position_dodge(width=0.2),stat="identity",width=0.1,fill="firebrick1") + theme_classic()
  t <- hp + coord_flip() + scale_y_reverse(limits=c(100,0)) + theme(axis.line.y = element_blank(),
                                                                    axis.title.y=element_blank(),
                                                                    axis.text.y=element_blank(),
                                                                    axis.ticks.y=element_blank(),
                                                                    axis.title.x=element_blank(),
                                                                    plot.margin=unit(c(1.5,0,2,1),"lines")) +
    geom_text(aes(label=paste0(round(perc,1)," %")), position=position_identity(),hjust = 5, size=2)
  t
}

top_plot <- function(df){
  ggplot(df,aes(x=type,y=perc)) + geom_bar(position=position_dodge(width=0.2),stat="identity",width=0.2,fill="firebrick1") + 
    ylim(0,100) + theme_classic() + theme(axis.line.x = element_blank(),
                                          axis.title.x=element_blank(),
                                          axis.text.x=element_blank(),
                                          axis.ticks.x=element_blank(),
                                          plot.margin=unit(c(1,1,1,3),"lines"))+
    geom_text(aes(label=paste0(round(perc,1)," %")), position=position_identity(),vjust = -5, size=2)
}

isEnrichAny <- function(ft,base,col,m_dir,matrice){
  dfb <- NULL
  featN <- glue::glue_collapse(ft,"_")
  for(base_set in base){
    t_set <- get(base_set) 
    for(test_set in col){
      title <- paste0("Line=",featN,"_Col=",test_set,"_Where=",base_set)
      out_d <- create(paste0(m_dir,"/BASE=",base_set,"/",featN,"_",test_set,"/"))
      t_set <- t_set & get(test_set)
      inter_set <- get(base_set)
      for(i in 1:length(ft)){
        inter_set <- inter_set & get(ft[i])
      }
      f_df <- NULL
      a_mat1 <- nrow(matrice[which(get(base_set) & inter_set & get(test_set)),])
      b_mat1 <- nrow(matrice[which(get(base_set) & !inter_set & get(test_set)),])
      a_mat2 <- nrow(matrice[which(get(base_set) & inter_set & !get(test_set)),])
      b_mat2 <- nrow(matrice[which(get(base_set) & !(inter_set) & !get(test_set)),])
      myFish.df <- myFisher(a_mat1,b_mat1,a_mat2,b_mat2,"two.sided")
      mypv <- formatC(as.numeric(myFish.df$ft.pv, format = "e", digits = 2))
      odds <- round(myFish.df$ft.fc,2)
      mycol <- ifelse(as.numeric(mypv) < 0.05,"darkgreen","red")
      myAnot <- paste0("Pvalue = ",mypv,"\n")
      myAnot <- paste0(myAnot,"Odds ratio = ",odds,"\n")
      myAnot <- paste0(myAnot,"Base set = ",base_set," genes")
      
      f_df <- data.frame(line=c(paste0("1_",featN),
                                paste0("1_",featN),
                                rep(paste0("2_no",featN),2)),
                         col=rep(c(paste0("1_",test_set),
                                   paste0("2_no",test_set)),2),
                         perc=with(myFish.df,
                                   c(a_mat1,a_mat2,b_mat1,b_mat2)))
      
      f_df$line <- factor(f_df$line)
      f_df$line <- factor(f_df$line,levels(f_df$line)[c(2,1)])
      dfl <- data.frame(perc=c(myFish.df$line2_perc,myFish.df$line1_perc),type=c("line1","line2"))
      dfc <- data.frame(perc=c(myFish.df$col1_perc,myFish.df$col2_perc),type=c("col1","col2"))
      
      mid <- tile_plot(f_df)
      lft <- flip_plot(dfl)
      top <- top_plot(dfc)
      anot <- blankPlot + annotate("text",label=myAnot,x=1,y=1,parse=F,col=mycol)
      ga <- arrangeGrob(anot ,top, blankPlot, 
                        lft, mid,blankPlot,
                        blankPlot,blankPlot,blankPlot,
                        blankPlot,blankPlot,blankPlot,
                        ncol=3, nrow=4, widths=c(2,2,1), heights=c(3,3,1,1))
      print("ok")
      ggsave(ga,file = paste0(out_d,"BaseSet=",base_set,"_genes_TestSet=",test_set,"_InterSet=",featN,".pdf"),width = 10,height = 10)
      # cat(paste0("BaseSet=",base_set,"_genes_TestSet=",test_set,"_InterSet=",featN,"\n \n"))
      # cat("Contingency matrix \n")
      # cat("----------------------");cat("\n")
      # print(matrix(c(a_mat1,a_mat2,b_mat1,b_mat2),ncol=2,nrow=2,byrow=T,dimnames = list(c(featN,paste0("no",featN)),c(test_set,paste0("no",test_set)))))
      # cat("\n\n")
      # cat("Fisher statistic \n")
      # cat("----------------------");cat("\n")
      # print(myFish.df)
      # cat(rep("\n",4))
      
      l_perc <- myFish.df$line1_perc
      lgt <- a_mat1
      mypv <- as.numeric(mypv)
      sign <- ifelse(mypv < 0.05,ifelse(mypv < 0.01,ifelse(mypv < 0.001,"***","**"),"*"),"NS")
      dfb <- rbind.data.frame(dfb,cbind.data.frame(test_set=test_set,inter_set=featN,base_set=base_set,
                                                   pv=mypv,odds=odds,lperc=l_perc,sign=sign))
      
      # sink()
    }
  }
  return(dfb)
}
# EXEMPLE
# myL <- c("z1.hg19.gr","mtr4.hg19.gr","ctcf.hg19.gr","rad21.hg19.gr",
#          "prompts.hg19.gr","prom2Kb.hg19.gr","enhCAGE.hg19.gr","ctcf.hg19.gr",
#          "genesBody.hg19.gr","h3k4me1.hg19.gr")
# xtVec <- c(T,T,T,T,F,F,T,T,F,T)   <- Granges of the associated vector you want to center and extend from xtend value
# xtend <- 500
# 
# plot(generateVennTHREE(myL,xtVec,xtend,c("ENHCAGE","Z1","MTR4")),doWeights=T)

isEnrichAny_G2I <- function(ft,base,col,out_d,matrice){
  dfb <- NULL
  featN <- glue::glue_collapse(ft,"_")
  for(base_set in base){
    t_set <- get(base_set) 
    for(test_set in col){
      title <- paste0("Line=",featN,"_Col=",test_set,"_Where=",base_set)
      t_set <- t_set & get(test_set)
      inter_set <- get(base_set)
      for(i in 1:length(ft)){
        inter_set <- inter_set & get(ft[i])
      }
      f_df <- NULL
      a_mat1 <- nrow(matrice[which(get(base_set) & inter_set & get(test_set)),])
      b_mat1 <- nrow(matrice[which(get(base_set) & !inter_set & get(test_set)),])
      a_mat2 <- nrow(matrice[which(get(base_set) & inter_set & !get(test_set)),])
      b_mat2 <- nrow(matrice[which(get(base_set) & !(inter_set) & !get(test_set)),])
      myFish.df <- myFisher(a_mat1,b_mat1,a_mat2,b_mat2,"two.sided")
      mypv <- formatC(as.numeric(myFish.df$ft.pv, format = "e", digits = 2))
      odds <- round(myFish.df$ft.fc,2)
      mycol <- ifelse(as.numeric(mypv) < 0.05,"darkgreen","red")
      myAnot <- paste0("Pvalue = ",mypv,"\n")
      myAnot <- paste0(myAnot,"Odds ratio = ",odds,"\n")
      myAnot <- paste0(myAnot,"Base set = ",base_set," genes")
      
      f_df <- data.frame(line=c(paste0("1_",featN),
                                paste0("1_",featN),
                                rep(paste0("2_no",featN),2)),
                         col=rep(c(paste0("1_",test_set),
                                   paste0("2_no",test_set)),2),
                         perc=with(myFish.df,
                                   c(a_mat1,a_mat2,b_mat1,b_mat2)))
      
      f_df$line <- factor(f_df$line)
      f_df$line <- factor(f_df$line,levels(f_df$line)[c(2,1)])
      dfl <- data.frame(perc=c(myFish.df$line2_perc,myFish.df$line1_perc),type=c("line1","line2"))
      dfc <- data.frame(perc=c(myFish.df$col1_perc,myFish.df$col2_perc),type=c("col1","col2"))
      
      mid <- tile_plot(f_df)
      lft <- flip_plot(dfl)
      top <- top_plot(dfc)
      anot <- blankPlot + annotate("text",label=myAnot,x=1,y=1,parse=F,col=mycol)
      ga <- arrangeGrob(anot ,top, blankPlot, 
                        lft, mid,blankPlot,
                        blankPlot,blankPlot,blankPlot,
                        blankPlot,blankPlot,blankPlot,
                        ncol=3, nrow=4, widths=c(2,2,1), heights=c(3,3,1,1))
      ggsave(ga,file = paste0(out_d,".png"),width = 10,height = 10)
      l_perc <- myFish.df$line1_perc
      lgt <- a_mat1
      mypv <- as.numeric(mypv)
      sign <- ifelse(mypv < 0.05,ifelse(mypv < 0.01,ifelse(mypv < 0.001,"***","**"),"*"),"NS")
      dfb <- rbind.data.frame(dfb,cbind.data.frame(test_set=test_set,inter_set=featN,base_set=base_set,
                                                   pv=mypv,odds=odds,lperc=l_perc,sign=sign))
      
    }
  }
  return(dfb)
}



#   -----------------------------------------------------------------------


# MY SEQPLOT FUN ----------------------------------------------------------
# Do a profile plot with multiple usefull option
# type = pf (start), mf (middle), ef (endpoint), af (normalize all GRanges to size 1, usefull for profile over genes where size may varies)
# bin = size of smoothing bin
# err = If you want to plot the standard error
# sd = number of Standard Deviation to allow in the plot (SD=3 -> 99.5 % of your values)
# smooth = If you want a cute plot shape (modulate spar options to change the smoothness)
# gnme = Genome of the BigWig to plot

seqPlotSDoutliers <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin,
                              sd=3,err=F,type="mf",smooth=F,spar=0.35,gnme="hg38",
                              ignore.strand=F,leg_pos="topright",
                              COLVEC=brewer.pal(n = length(bw.l)*length(gr.v)+1 , name = "Paired"),
                              main=""){
  mySMPLseq <- seq(0,1,1e-3)
  bw.n <- NULL
  o.tmp <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  for(mygr in gr.v){
    sze <- length(get(mygr))
    print(mygr)
    o.tmp <- c(o.tmp,toString(export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed"))))
  }
  print("ok")
  gpsa <- getPlotSetArray(bw.l,o.tmp,gnme,bin = bin,ignore_strand = ignore.strand,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type)
  gpsa.data <- gpsa$data
  print("totp")
  for(mygr in gr.v){
    for(my.bw in bw.n){
      sze <- length(get(mygr))
      # To avoid column entirely covered with 0 (issue with scale and smooth)
      gpsa.mtx <- data.frame(apply(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]]
                                   ,c(1,2),function(X) X + sample(mySMPLseq,1)))
      
      gpsa.scl.mtx <- gpsa.mtx %>% mutate_all(scale) # scale the data (center reduce)
      gpsa.scl.mtx[abs(gpsa.scl.mtx) > sd] <- NA # Remove value X SD away (sd = 3 by default ~ 98% of the data)
      gpsa.scl.mtx[apply(gpsa.scl.mtx,c(1,2),is.nan)] <- NA
      means <- colMeans(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      if(smooth){
        means = smooth.spline(1:(length(means)), means, spar=spar)$y
      }
      
      
      # OLD
      
      # gpsa.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]])
      # gpsa.scl.mtx <- gpsa.mtx %>% mutate_all(scale) # scale the data (center reduce)
      # gpsa.scl.mtx[abs(gpsa.scl.mtx) > sd] <- NA # Remove value X SD away (sd = 3 by default ~ 98% of the data)
      # means <- colMeans(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      # if(smooth){
      #   means = smooth.spline(1:(length(means)), means, spar=spar)$y
      # }
      
      
      # TEMPORARY START
      # if(grepl("_WT_*",my.bw)){ means <-  means +1}
      # TEMPORARY END
      stderror <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,2,function(n){
        sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
      })
      conint <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx, 2, function(n) {
        qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
      })
      stderror[is.na(stderror)] <- 0
      conint[is.na(conint)] <- 0
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means # change the means vector from getPlotSetArray object
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["stderror"]] <- stderror # change the means vector from getPlotSetArray object
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["conint"]] <- conint  # change the means vector from getPlotSetArray object
    }
    
  }
  file.remove(paste0(o.tmp))
  plotAverage(gpsa,xlab='Relative position [bp]', 
              ylim=ylim, ylab='Signal',
              main = main, 
              keepratio = F,error.estimates = err,
              frame.plot=F,
              cex.legend = 8,
              legend_pos = leg_pos,
              pointsize = 16,
              colvec = COLVEC)
  
}






seqPlotFastSD <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin,pscale="linear",
                          sd=3,err=F,type="mf",smooth=F,spar=0.35,gnme="hg38",
                          ignore.strand=F,leg_pos="topright",
                          COLVEC=brewer.pal(n = length(bw.l)*length(gr.v)+1 , name = "Paired"),
                          main="",sum_type="mean",shift=NULL){
  bw.n <- NULL
  o.tmp <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  for(mygr in gr.v){
    sze <- length(get(mygr))
    print(mygr)
    o.tmp <- c(o.tmp,toString(export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed"))))
  }
  print("ok")
  gpsa <- getPlotSetArray(bw.l,o.tmp,gnme,bin = bin,ignore_strand = ignore.strand,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type)
  gpsa.data <- gpsa$data
  print("totp")
  for(mygr in gr.v){
    for(my.bw in bw.n){
      sze <- length(get(mygr))
      gpsa.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]])
      
      gpsa.scl.mtx <- gpsa.mtx %>% mutate_all(scale) # scale the data (center reduce)
      gpsa.scl.mtx[abs(gpsa.scl.mtx) > sd] <- NA # Remove value X SD away (sd = 3 by default ~ 98% of the data)
      if(sum_type=="mean"){
        means <- colMeans(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      }else{
        mat_c <- gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx
        means <- apply(mat_c,2,function(X){
          median(X,na.rm=T)
        })
      }
      
      if(smooth){
        means = smooth.spline(1:(length(means)), means, spar=spar)$y
      }
      stderror <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,2,function(n){
        sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
      })
      conint <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx, 2, function(n) {
        qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
      })
      stderror[is.na(stderror)] <- 0
      conint[is.na(conint)] <- 0
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means # change the means vector from getPlotSetArray object
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["stderror"]] <- stderror # change the means vector from getPlotSetArray object
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["conint"]] <- conint  # change the means vector from getPlotSetArray object
    }
    
  }
  file.remove(paste0(o.tmp))
  plotAverage(gpsa,xlab='Relative position [bp]', 
              ylim=ylim, ylab='Signal',
              main = main, 
              plotScale = pscale,
              keepratio = F,error.estimates = err,
              frame.plot=F,
              cex.legend = 8,
              legend_pos = leg_pos,
              pointsize = 16,
              colvec = COLVEC) 
  
}


seqPlotFastPERC <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin,pscale="linear",centile=1:100,
                            normFAC=1,
                            sd=3,err=F,type="mf",smooth=F,spar=0.35,gnme="hg38",
                            ignore.strand=F,leg_pos="topright",
                            COLVEC=brewer.pal(n = length(bw.l)*length(gr.v)+1 , name = "Paired"),
                            main="",sum_type="mean"){
  bw.n <- NULL
  o.tmp <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  for(mygr in gr.v){
    sze <- length(get(mygr))
    print(mygr)
    o.tmp <- c(o.tmp,toString(export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed"))))
  }
  print("ok")
  gpsa <- getPlotSetArray(bw.l,o.tmp,gnme,bin = bin,ignore_strand = ignore.strand,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type)
  gpsa.data <- gpsa$data
  print("totp")
  
  for(mygr in gr.v){
    for(IDX in 1:length(bw.n)){
      my.bw <- bw.n[IDX]
      normF <- normFAC[IDX]
      sze <- length(get(mygr))
      gpsa.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]])
      
      RS <- data.frame(RS=rowSums(gpsa.mtx,na.rm=T))
      RS <- RS %>% mutate(CENT=ntile(RS,100))
      ROW_KEEP <- RS[,"CENT"] %in% centile[IDX]
      
      gpsa.mtx <- gpsa.mtx[ROW_KEEP,] 
      if(sum_type=="mean"){
        means <- colMeans(gpsa.mtx,na.rm=T)*normF # Now you can do the mean on original data without 3 SD away outliers
      }else{
        mat_c <- gpsa.mtx 
        means <- apply(mat_c,2,function(X){
          median(X,na.rm=T)
        })*normF
      }
      
      if(smooth){
        means = smooth.spline(1:(length(means)), means, spar=spar)$y
      }
      stderror <- apply(gpsa.mtx,2,function(n){
        sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
      })
      conint <- apply(gpsa.mtx , 2, function(n) {
        qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
      })
      stderror[is.na(stderror)] <- 0
      conint[is.na(conint)] <- 0
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means # change the means vector from getPlotSetArray object
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["stderror"]] <- stderror # change the means vector from getPlotSetArray object
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["conint"]] <- conint  # change the means vector from getPlotSetArray object
      
    }
    
  }
  file.remove(paste0(o.tmp))
  plotAverage(gpsa,xlab='Relative position [bp]', 
              ylim=ylim, ylab='Signal',
              main = main, 
              plotScale = pscale,
              keepratio = F,error.estimates = err,
              frame.plot=F,
              cex.legend = 8,
              legend_pos = leg_pos,
              pointsize = 16,
              colvec = COLVEC) 
  
}


seqPlotFastPERC_LST <- function(struct.ll,tmp,ylim,xlim=xlim,bin=bin,pscale="linear",
                                err=F,type="mf",smooth=F,spar=0.35,gnme="hg38",
                                ignore.strand=F,leg_pos="topright",
                                COLVEC=brewer.pal(n = length(bw.l)*length(gr.v)+1 , name = "Paired"),
                                main="",sum_type="mean", lwd=2){
  
  bw.n <- NULL
  o.tmp <- NULL
  BW.l <- NULL
  trash <- sapply(names(struct.ll),function(X){
    subS.l <- struct.ll[[X]]
    BW <- subS.l[["BW"]]
    GR <- subS.l[["GR"]]
    bw.n[[X]] <<- gsub("(.*).bw","\\1",basename(BW))
    
    sze <- length(get(GR))
    o.tmp <<- unique(c(o.tmp,toString(export.bed(get(GR),paste0(tmp,"/",GR,"_#",sze,"peaks.bed")))))
    BW.l <<- unique(c(BW.l,BW))
  })
  
  gpsa <- getPlotSetArray(BW.l,o.tmp,gnme,bin = bin,ignore_strand = ignore.strand,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type)
  gpsa.data <- gpsa$data
  
  trash <- sapply(names(struct.ll),function(X){
    subS.l <- struct.ll[[X]]
    BW <- subS.l[["BW"]]
    GR <- subS.l[["GR"]]
    bw.n[[X]] <<- gsub("(.*).bw","\\1",basename(BW))
    normF <- subS.l[["NORM"]]
    CENT <- subS.l[["CENTILE"]]
    
    sze <- length(get(GR))
    gpsa.mtx <- data.frame(gpsa.data[[paste0(GR,"_#",sze,"peaks")]][[bw.n[[X]]]][["heatmap"]])
    RS <- data.frame(RS=rowSums(gpsa.mtx,na.rm=T))
    RS <- RS %>% mutate(CENT=ntile(RS,100))
    ROW_KEEP <- RS[,"CENT"] %in% CENT
    
    gpsa.mtx <- gpsa.mtx[ROW_KEEP,] 
    if(sum_type=="mean"){
      if(normF[["BOOL"]]==T){
        if(normF[["TYPE"]]=="exp"){
          means <- exp(colMeans(gpsa.mtx,na.rm=T))*normF[["LINFAC"]] # Now you can do the mean on original data without 3 SD away outliers
        }else if(normF[["TYPE"]]=="pow"){
          means <- `^`(colMeans(gpsa.mtx,na.rm=T),normF[["FAC"]])*normF[["LINFAC"]]
        }else if(normF[["TYPE"]]=="log"){
          means <- log(colMeans(gpsa.mtx,na.rm=T),base = normF[["FAC"]])*normF[["LINFAC"]]
        }else if(normF[["TYPE"]]=="MIN"){
          means <- colMeans(gpsa.mtx,na.rm=T)-normF[["LINFAC"]]
        }else{
          means <- colMeans(gpsa.mtx,na.rm=T)*normF[["LINFAC"]]
        }
      }else{
        means <- colMeans(gpsa.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      }
    }else{
      mat_c <- gpsa.mtx 
      
      if(normF[["BOOL"]]==T){
        if(normF[["TYPE"]]=="exp"){
          means <- exp(apply(mat_c,2,function(X){
            median(X,na.rm=T)
          }))*normF[["LINFAC"]] # Now you can do the mean on original data without 3 SD away outliers
        }else if(normF[["TYPE"]]=="pow"){
          means <- `^`(apply(mat_c,2,function(X){
            median(X,na.rm=T)
          }),normF[["FAC"]])*normF[["LINFAC"]]
        }else if(normF[["TYPE"]]=="log"){
          means <- log(apply(mat_c,2,function(X){
            median(X,na.rm=T)
          }),base = normF[["FAC"]])*normF[["LINFAC"]]
        }else if(normF[["TYPE"]]=="log"){
          means <- apply(mat_c,2,function(X){
            median(X,na.rm=T)
          })-normF[["LINFAC"]]
        }else{
          means <- apply(mat_c,2,function(X){
            median(X,na.rm=T)
          })*normF[["LINFAC"]]
        }
      }else{
        means <- colMeans(gpsa.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      }
    }
    print(means)
    if(smooth){
      means = smooth.spline(1:(length(means)), means, spar=spar)$y
    }
    stderror <- apply(gpsa.mtx,2,function(n){
      sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
    })
    conint <- apply(gpsa.mtx , 2, function(n) {
      qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
    })
    stderror[is.na(stderror)] <- 0
    conint[is.na(conint)] <- 0
    gpsa$data[[paste0(GR,"_#",sze,"peaks")]][[bw.n[[X]]]][["means"]] <<- means # change the means vector from getPlotSetArray object
    gpsa$data[[paste0(GR,"_#",sze,"peaks")]][[bw.n[[X]]]][["stderror"]] <<- stderror # change the means vector from getPlotSetArray object
    gpsa$data[[paste0(GR,"_#",sze,"peaks")]][[bw.n[[X]]]][["conint"]] <<- conint  # change the means vector from getPlotSetArray object
    
  })
  file.remove(paste0(o.tmp))
  plotAverage(gpsa,xlab='Relative position [bp]', 
              ylim=ylim, ylab='Signal',
              main = main, 
              plotScale = pscale,
              keepratio = F,error.estimates = err,
              frame.plot=F,
              cex.legend = 8,
              legend_pos = leg_pos,
              pointsize = 16,
              ln.v=F,
              colvec = COLVEC) 
}




seqPlotSmooth <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin,
                          err=F,type="mf",smooth=F,spar=0.35,gnme="hg38",
                          ignore.strand=F,leg_pos="topright",
                          COLVEC=brewer.pal(n = length(bw.l)*length(gr.v)+1 , name = "Paired"),
                          main=""){
  bw.n <- NULL
  o.tmp <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  for(mygr in gr.v){
    sze <- length(get(mygr))
    print(mygr)
    o.tmp <- c(o.tmp,toString(export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed"))))
  }
  print("ok")
  gpsa <- getPlotSetArray(bw.l,o.tmp,gnme,bin = bin,ignore_strand = ignore.strand,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type)
  gpsa.data <- gpsa$data
  print("totp")
  for(mygr in gr.v){
    for(my.bw in bw.n){
      sze <- length(get(mygr))
      # To avoid column entirely covered with 0 (issue with smooth)
      gpsa.mtx <- gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]]
      means <- colMeans(gpsa.mtx, na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      if(smooth){
        means = smooth.spline(1:(length(means)), means, spar=spar)$y
      }
      gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means # change the means vector from getPlotSetArray object
    }
    
  }
  file.remove(paste0(o.tmp))
  plotAverage(gpsa,xlab='Relative position [bp]', 
              ylim=ylim, ylab='Signal',
              main = main, 
              keepratio = F,error.estimates = err,
              frame.plot=F,
              cex.legend = 8,
              legend_pos = leg_pos,
              pointsize = 16,
              colvec = COLVEC)
  
}





seqPlotSDoutliersDIV <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin,
                                 sd=3,err=F,type="mf",smooth=F,spar=0.35,gnme="hg38",
                                 ignore.strand=F,
                                 COLVEC=brewer.pal(n = length(bw.l)*length(gr.v)+1 , name = "Paired"),
                                 main=""){
  bw.n <- NULL
  o.tmp <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  for(mygr in gr.v){
    sze <- length(get(mygr))
    print(mygr)
    o.tmp <- c(o.tmp,toString(export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed"))))
  }
  print("ok")
  gpsa <- getPlotSetArray(bw.l,o.tmp,gnme,bin = bin,ignore_strand = ignore.strand,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type)
  gpsa.data <- gpsa$data
  for(mygr in gr.v){
    for(my.bw in bw.n){
      sze <- length(get(mygr))
      if(grepl("wtxh3",my.bw)){
        print("WT")
        print(my.bw)
        
        gpsa.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]])
        gpsa.scl.mtx <- gpsa.mtx %>% mutate_all(scale) # scale the data (center reduce)
        gpsa.scl.mtx[abs(gpsa.scl.mtx) > sd] <- NA # Remove value X SD away (sd = 3 by default ~ 98% of the data)
        means <- colMeans(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
        # if(smooth){
        #   means = smooth.spline(1:(length(means)), means, spar=spar)$y
        # }
        # TEMPORARY START
        # if(grepl("_WT_*",my.bw)){ means <-  means +1}
        # TEMPORARY END
        stderror <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,2,function(n){
          sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
        })
        conint <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx, 2, function(n) {
          qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
        })
        stderror[is.na(stderror)] <- 0
        conint[is.na(conint)] <- 0
        BWT_NAME <- my.bw
        WT_MEANS <- means # change the means vector from getPlotSetArray object
        print(WT_MEANS)
        WT_ERR <- stderror # change the means vector from getPlotSetArray object
        WT_CONINT <- conint  # change the means vector from getPlotSetArray object
      }else if(grepl("kdxh3",my.bw)){
        print("kd")
        print(my.bw)
        gpsa.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]])
        gpsa.scl.mtx <- gpsa.mtx %>% mutate_all(scale) # scale the data (center reduce)
        gpsa.scl.mtx[abs(gpsa.scl.mtx) > sd] <- NA # Remove value X SD away (sd = 3 by default ~ 98% of the data)
        means <- colMeans(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
        # if(smooth){
        #   means = smooth.spline(1:(length(means)), means, spar=spar)$y
        # }
        # TEMPORARY START
        # if(grepl("_WT_*",my.bw)){ means <-  means +1}
        # TEMPORARY END
        stderror <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,2,function(n){
          sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
        })
        conint <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx, 2, function(n) {
          qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
        })
        stderror[is.na(stderror)] <- 0
        conint[is.na(conint)] <- 0
        BKD_NAME <- my.bw
        
        KD_MEANS <- means # change the means vector from getPlotSetArray object
        print(KD_MEANS)
        
        KD_ERR <- stderror # change the means vector from getPlotSetArray object
        KD_CONINT <- conint  # change the means vector from getPlotSetArray object
      }else{
        gpsa.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]])
        gpsa.scl.mtx <- gpsa.mtx %>% mutate_all(scale) # scale the data (center reduce)
        gpsa.scl.mtx[abs(gpsa.scl.mtx) > sd] <- NA # Remove value X SD away (sd = 3 by default ~ 98% of the data)
        means <- colMeans(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
        if(smooth){
          means = smooth.spline(1:(length(means)), means, spar=spar)$y
        }
        # TEMPORARY START
        # if(grepl("_WT_*",my.bw)){ means <-  means +1}
        # TEMPORARY END
        stderror <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,2,function(n){
          sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
        })
        conint <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx, 2, function(n) {
          qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
        })
        stderror[is.na(stderror)] <- 0
        conint[is.na(conint)] <- 0
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means # change the means vector from getPlotSetArray object
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["stderror"]] <- stderror # change the means vector from getPlotSetArray object
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["conint"]] <- conint  # change the means vector from getPlotSetArray object
      }
    }
  }
  
  
  myMEAN <- apply(cbind(KD_MEANS,WT_MEANS),1,mean)
  
  myVEC <- 2*(KD_MEANS - WT_MEANS)/sqrt(abs(myMEAN))
  myVEC = smooth.spline(1:(length(myVEC)), myVEC, spar=spar)$y + 0.2
  
  
  gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[BKD_NAME]][["means"]] <- myVEC
  gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[BKD_NAME]][["stderror"]] <- WT_ERR # change the means vector from getPlotSetArray object
  gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[BKD_NAME]][["conint"]] <- WT_CONINT  # change the means vector from getPlotSetArray object
  
  gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[BWT_NAME]][["means"]] <-myVEC
  gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[BWT_NAME]][["stderror"]] <- WT_ERR # change the means vector from getPlotSetArray object
  gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[BWT_NAME]][["conint"]] <- WT_CONINT  # change the means vector from getPlotSetArray object
  
  file.remove(paste0(o.tmp))
  plotAverage(gpsa,xlab='Relative position [bp]', 
              ylim=ylim, ylab='Signal',
              main = main, 
              keepratio = F,error.estimates = err,
              frame.plot=F,
              cex.legend = 8,
              legend_pos = "topright",
              pointsize = 16,lwd=12,
              colvec = COLVEC)
  
}



diff_HeatmapPlot <- function(bw1,bw2,tmp,mygr,Xlim=c(5000,5000),bin=100,type="mf",smooth=F,spar=0.8,genome="hg19",ignore.strand=F){
  
  bw1.n <- gsub("(.*).bw","\\1",basename(bw1))
  bw2.n <- gsub("(.*).bw","\\1",basename(bw2))
  sze <- length(get(mygr))
  print(mygr)
  o.tmp <- toString(export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed")))
  gpsa1 <- getPlotSetArray(bw1,o.tmp,genome,bin = bin,ignore_strand = ignore.strand,xmin = Xlim[1],xmax=Xlim[2],rm0 = F,type=type)
  gpsa2 <- getPlotSetArray(bw2,o.tmp,genome,bin = bin,ignore_strand = ignore.strand,xmin = Xlim[1],xmax=Xlim[2],rm0 = F,type=type)
  
  
  
  gpsa.data1 <- gpsa1$data
  gpsa.data2 <- gpsa2$data
  
  gpsaONE.mtx <- as.matrix(gpsa.data1[[paste0(mygr,"_#",sze,"peaks")]][[bw1.n]][["heatmap"]]) + 1 
  gpsaTWO.mtx <- as.matrix(gpsa.data2[[paste0(mygr,"_#",sze,"peaks")]][[bw2.n]][["heatmap"]]) + 1
  
  X <- list(gpsaONE.mtx, gpsaTWO.mtx)
  Y <- do.call(cbind, X)
  Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
  
  meanMat <- abs(apply(Y, c(1, 2), mean, na.rm = TRUE))
  ZS.HTMP <- (gpsaONE.mtx - gpsaTWO.mtx)/sqrt(meanMat)
  OK.HTMP <- ZS.HTMP
  OK.HTMP[is.na(OK.HTMP)] <- 0
  if(smooth){
    OK.HTMP <- t(apply(OK.HTMP,1,function(X)smooth.spline(1:(length(X)), X, spar=0.8)$y))
  }
  gpsa1$data[[paste0(mygr,"_#",sze,"peaks")]][[bw1.n]][["heatmap"]] <- OK.HTMP
  
  file.remove(paste0(o.tmp))
  return(gpsa1)
}


heatmapStrandAware <- function(bw.l,tmp,mygr,xlim=xlim,bin=bin,type="mf",genome="hg19",invertscorerev=F,invertscorefwd=F){
  bw.n <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  sze <- length(get(mygr))
  print(mygr)
  o.tmp <- toString(export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed")))
  gpsa <- getPlotSetArray(bw.l,o.tmp,genome,bin = bin,ignore_strand = T,xmin = xlim[1],xmax=xlim[2],rm0 = T,type=type)
  gpsa.data <- gpsa$data
  for(myREV.bw in bw.n[grep("rev",bw.n)]){
    
    
    myFWD.bw <- gsub("(rev)","fwd",myREV.bw)
    gpsaREV.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[myREV.bw]][["heatmap"]])
    gpsaFWD.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[myFWD.bw]][["heatmap"]])
    if(invertscorerev) gpsaREV.mtx <- gpsaREV.mtx * -1   
    if(invertscorefwd) gpsaFWD.mtx <- gpsaFWD.mtx * -1   
    
    selNeg <- !(get(mygr)$str)
    tmporary <- gpsaREV.mtx[selNeg,ncol(gpsaREV.mtx):1]
    gpsaREV.mtx[selNeg,] <- gpsaFWD.mtx[selNeg,ncol(gpsaFWD.mtx):1]
    gpsaFWD.mtx[selNeg,] <- tmporary
    
    gpsaREV.mtx <- as.matrix(gpsaREV.mtx)
    gpsaFWD.mtx <- as.matrix(gpsaFWD.mtx)
    gpsaREV.mtx[which(is.na(gpsaREV.mtx))] <- 0
    gpsaFWD.mtx[which(is.na(gpsaFWD.mtx))] <- 0
    gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[myREV.bw]][["heatmap"]] <- gpsaREV.mtx
    gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[myFWD.bw]][["heatmap"]] <- gpsaFWD.mtx
  }
  file.remove(paste0(o.tmp))
  return(gpsa)
}


heatmapStrandAwareMerged <- function(bw.l,tmp,mygr,xlim=xlim,bin=bin,type="mf",genome="hg19",ign.str =F,invertscorerev=F,invertscorefwd=F){
  bw.n <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  sze <- length(get(mygr))
  print(mygr)
  o.tmp <- toString(export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed")))
  gpsa <- getPlotSetArray(bw.l,o.tmp,genome,bin = bin,ignore_strand = T,xmin = xlim[1],xmax=xlim[2],rm0 = T,type=type)
  gpsa.data <- gpsa$data
  for(myREV.bw in bw.n[grep("rev",bw.n)]){
    
    
    myFWD.bw <- gsub("(rev)","fwd",myREV.bw)
    gpsaREV.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[myREV.bw]][["heatmap"]])
    gpsaFWD.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[myFWD.bw]][["heatmap"]])
    if(invertscorerev) gpsaREV.mtx <- gpsaREV.mtx * -1   
    if(invertscorefwd) gpsaFWD.mtx <- gpsaFWD.mtx * -1   
    
    if(!ign.str){
      selNeg <- !(get(mygr)$str)
      tmporary <- gpsaREV.mtx[selNeg,ncol(gpsaREV.mtx):1]
      gpsaREV.mtx[selNeg,] <- gpsaFWD.mtx[selNeg,ncol(gpsaFWD.mtx):1]
      gpsaFWD.mtx[selNeg,] <- tmporary
    }
    gpsaREV.mtx <- as.matrix(gpsaREV.mtx)
    gpsaFWD.mtx <- as.matrix(gpsaFWD.mtx)
    gpsaREV.mtx[which(is.na(gpsaREV.mtx))] <- 0
    gpsaFWD.mtx[which(is.na(gpsaFWD.mtx))] <- 0
    gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[myREV.bw]][["heatmap"]] <- gpsaREV.mtx + gpsaFWD.mtx
    gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[myFWD.bw]][["heatmap"]] <- gpsaFWD.mtx + gpsaREV.mtx
  }
  file.remove(paste0(o.tmp))
  return(gpsa)
}





heatmapClassic <- function(bw.l,tmp,mygr,xlim=xlim,bin=bin,type="pf",ign.str=F,genome="dm3"){
  bw.n <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  sze <- length(get(mygr))
  print(mygr)
  o.tmp <- toString(export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed")))
  gpsa <- getPlotSetArray(bw.l,o.tmp,genome,bin = bin,ignore_strand = ign.str,xmin = xlim[1],xmax=xlim[2],rm0 = T,type=type)
  return(gpsa)
}


boxplotStrandAware <- function(bw.l,tmp,mygr,xlim=xlim,bin=bin,type="mf"){
  bw.n <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  sze <- length(get(mygr))
  print(mygr)
  o.tmp <- toString(export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed")))
  gpsa <- getPlotSetArray(bw.l,o.tmp,"hg38",bin = bin,ignore_strand = T,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type)
  gpsa.data <- gpsa$data
  for(myREV.bw in bw.n[grep("rev",bw.n)]){
    myFWD.bw <- gsub("(rev)","fwd",myREV.bw)
    gpsaREV.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[myREV.bw]][["heatmap"]])
    gpsaFWD.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[myFWD.bw]][["heatmap"]])
    selNeg <- !(get(mygr)$str)
    gpsaREV.mtx[selNeg,] <- gpsaFWD.mtx[selNeg,ncol(gpsaFWD.mtx):1]
    gpsaFWD.mtx[selNeg,] <- gpsaREV.mtx[selNeg,ncol(gpsaREV.mtx):1]
    
    gpsaREV.mtx <- as.matrix(gpsaREV.mtx)
    gpsaFWD.mtx <- as.matrix(gpsaFWD.mtx)
    
    gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[myREV.bw]][["heatmap"]] <- gpsaREV.mtx
    gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[myFWD.bw]][["heatmap"]] <- gpsaFWD.mtx
  }
  file.remove(paste0(o.tmp))
  return(gpsa)
}



myplotfun <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin,err=F,type="mf"){
  for(mygr in gr.v){
    sze <- length(get(mygr))
    print(mygr)
    o.tmp <- toString(export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed")))
    plotTMP <- getPlotSetArray(bw.l,o.tmp,"hg38",bin = bin,ignore_strand = F,xmin = xlim[1],xmax=xlim[2],rm0 = F,type = type)
    file.remove(paste0(o.tmp))
    print(plotAverage(plotTMP,xlab='Relative position [bp]', ylim=ylim, ylab='Signal',main = "ALL SIGNAL",keepratio = F,error.estimates = err,colvec =c("dodgerblue","firebrick2") ))
  }
}
#   -----------------------------------------------------------------------




# GROUPING FUNCTION -------------------------------------------------------
# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')
# exemple g(c, d) %=% c(1,2)

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}


# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}


#   -----------------------------------------------------------------------


#  RANGE FUNCTION ---------------------------------------------------------


# WHICH DOMAINS ---
# Input : ChIPseq GenomicRanges, TADs Genomic Ranges
# Return ChIPseq GenomicRanges with the indice and size of the corresponding domain it falls in
whichDomain <- function(peak.GR,dom.GR){
  nPeak.GR <- resize(peak.GR,1,fix="center")
  peak_dom.OL <- findOverlaps(nPeak.GR,dom.GR,type ="within")
  peak_to_keep.GR <- peak.GR[queryHits(peak_dom.OL)]
  dom.IDX <- subjectHits(peak_dom.OL)
  peak_to_keep.GR$dom <- dom.IDX
  peak_to_keep.GR$domSize <- width(dom.GR[dom.IDX])/2
  return(peak_to_keep.GR)
}


#   -----------------------------------------------------------------------

#  ROC CURVE ---------------------------------------------------------------
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

#   -----------------------------------------------------------------------

#  FISHER ------------------------------------------------------------------

myFisher <- function(a_mat1,b_mat1,a_mat2,b_mat2,meth){
  col1_perc <- (a_mat1*100)/(a_mat1+b_mat1) # 58%
  line1_perc <- (a_mat1*100)/(a_mat1+a_mat2) # 58%
  col2_perc <- (a_mat2*100)/(a_mat2+b_mat2) # 58%
  line2_perc <- (b_mat1*100)/(b_mat1+b_mat2) # 58%
  
  
  
  f_mat <- matrix(c(line_col=a_mat1,line_nocol=a_mat2,noline_col=b_mat1,noline_nocol=b_mat2),nrow=2,ncol=2,byrow=T)
  ft <- fisher.test(f_mat,alternative = meth)
  ft.pv <- ft$p.value # 0.93
  ft.fc <- ft$estimate # 0.97
  return(list(ft.pv=ft.pv,ft.fc=ft.fc,
              line1_perc=line1_perc,
              col1_perc=col1_perc,
              line2_perc=line2_perc,
              col2_perc=col2_perc,
              a_num=a_mat1,b_num=b_mat1,matF=f_mat))
}

fisherDistPlotBothSided <- function(anchor,bait,inter,maxrange,step,window,out){
  # Output setup
  main <- paste0("Anchor=",anchor,"_Bait=",bait,"_Intersect=",inter,"_maxrange=",maxrange,"KB_step=",step,"_Win=",window)
  out <- create(paste0(out,"/Anchor=",anchor,"/Filter=",bait,"/Range=",maxrange,"/Window=",window,"/"))
  wdth <- 25000/(maxrange/step)
  o_pdf <- paste0(out,"/",main,".pdf")
  #Get the data : from character to vector
  anchor <- get(anchor);bait <- get(bait);inter <- get(inter)
  # Turn TSS anchor point to Granges vector
  anchor.gr <- mrc_tss.gr[anchor]
  # Create a reference set of matrecap GRanges index that fall into the range of analysis
  refset_idx <- unique(subjectHits(findOverlaps(trim(anchor.gr+(maxrange+window)),mrc_tss.gr)))
  # get filter set
  bait_idx <- which(bait==1)
  # Apply refset and filterset
  tokeep <- intersect(refset_idx,bait_idx)
  mrc_filt_refeset.df <- matrecap[tokeep,]
  mrc_filt_refset.gr <- mrc_tss.gr[tokeep]
  refset_size <- length(mrc_filt_refset.gr)
  inter_filt_refset.gid <- mrc_filt_refset.gr[which(inter[tokeep]==1)]$gene_id
  inter_size <- length(inter_filt_refset.gid)
  shft <- seq(0,maxrange,step)
  window.gr <- resize(anchor.gr,window,"start",ignore.strand=T)
  myvec <- sapply(shft,
                  function(sshft){
                    sbstp.gr <- subsetByOverlaps(mrc_filt_refset.gr,trim(shift(window.gr,sshft)))
                    sbstn.gr <- subsetByOverlaps(mrc_filt_refset.gr,trim(shift(window.gr,-(sshft+window)))) # Added
                    sbst.gr <- unique(c(sbstp.gr,sbstn.gr)) # Added
                    if(length(sbst.gr)){
                      win_set <- sbst.gr$gene_id
                      a1 <- length(which(win_set %in% inter_filt_refset.gid))
                      a2 <- length(win_set) - a1
                      b1 <- inter_size - a1
                      b2 <- refset_size - (inter_size + length(win_set) - a1)
                      # print(matrix(c(a1,a2,b1,b2),2,2))
                      pv <- fisher.test(matrix(c(a1,a2,b1,b2),2,2),alternative = "greater")$p.value
                      perc <- 100*(a1/(a1+a2))
                      return(cbind(pv,perc))
                    }
                    else
                      return(cbind(1,NA))
                  })
  perc.df <- NULL
  myvec <- t(myvec)
  
  trsh <- sapply(1:(nrow(myvec)-3),
                 function(idx){
                   # print(idx)
                   perc.df <<- rbind.data.frame(perc.df,cbind.data.frame(mean(myvec[c(idx:(idx+3)),2]),sd(myvec[c(idx:(idx+3)),2])))
                 })
  perc.df <- cbind(perc.df,shft[-((length(shft)-2):length(shft))]+window/2)
  colnames(perc.df) <- c("mean","sd","x")
  toPlot <- data.frame(x=(shft+(window/2)),y=myvec[,1])
  toPlot$padj <- p.adjust(toPlot$y,method = "BH",n=length(toPlot$y))
  toPlot$logpv <- -log10(toPlot$padj)
  toPlot$logpv <- ifelse(toPlot$logpv > 8,8,toPlot$logpv)
  toPlot$thr <- ifelse(toPlot$logpv > -log10(0.05),toPlot$logpv,0)
  toPlot$thr <- ifelse(toPlot$thr > 8,8,toPlot$thr)
  
  
  a <- ggplot(toPlot,aes(x=x,y=padj)) + ggtitle(main) + ylim(0,1) + geom_line() + ylab("Padj") + theme(plot.title = element_text(size = 10, face = "bold"))
  b <- ggplot(toPlot,aes(x=x,y=logpv)) + ggtitle(main) + geom_line() + ylim(0,8) + ylab("logPVadj")+ theme(plot.title = element_text(size = 10, face = "bold"))
  c <- ggplot(toPlot,aes(x=x,y=thr)) + ggtitle(main) + geom_line() + ylim(0,8) + ylab("logPVadjThresholded")+ theme(plot.title = element_text(size = 10, face = "bold"))
  d <- ggplot(toPlot,aes(x=x,y=y)) + ggtitle(main) + geom_point(size=0.05) + ylim(0,1) + ylab("raw_pv")+ theme(plot.title = element_text(size = 10, face = "bold"))
  e <- ggplot(perc.df, aes(x=x , y=mean) ) + geom_bar(position=position_dodge(.2),stat="identity",width = wdth)   +
    scale_y_continuous(name="% DR/UR", limits=c(0, 100)) + theme_bw() + rotate_x_text(45) + 
    geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),size=0.5,width=.15,position=position_dodge(.2))
  
  pdf(o_pdf)
  print(a);print(b);print(c);print(d);print(e)
  dev.off()
  return(perc.df)
}


fisherDPorientated <- function(anchor,focus,refForw,refRev,bait,inter,maxrange,step,window,out){
  main <- paste0("Anchor=",anchor,"_Focus=",focus,"_Bait=",bait,"_Intersect=",inter,"_maxrange=",maxrange,"KB_step=",step,"_Win=",window)
  out <- create(paste0(out,"/Focus=",focus,"/Anchor=",anchor,"/Bait=",bait,"/Range=",maxrange,"/Window=",window,"/"))
  wdth <- 25000/(maxrange/step)
  o_pdf <- paste0(out,"/",main,".pdf")
  refForw <- get(refForw);refRev <- get(refRev);bait <- get(bait);inter <- get(inter)
  # get the reference set
  refF.gr <- mrc_tss.gr[refForw]
  refR.gr <- mrc_tss.gr[refRev]
  refsetF_idx <- unique(subjectHits(findOverlaps(trim(resize(refF.gr,fix = "start",width=(maxrange+window),ignore.strand=T)),mrc_tss.gr)))
  refsetR_idx <- unique(subjectHits(findOverlaps(trim(resize(refR.gr,fix="start",width=(maxrange+window),ignore.strand=T)),mrc_tss.gr)))
  refset_idx <- unique(c(refsetF_idx,refsetR_idx))
  # get filter set
  bait_idx <- which(bait==1)
  # Apply refset and filterset
  tokeep <- intersect(refset_idx,bait_idx)
  mrc_bait_refeset.df <- matrecap[tokeep,]
  mrc_bait_refset.gr <- mrc_tss.gr[tokeep]
  refset_size <- length(mrc_bait_refset.gr)
  inter_bait_refset.gid <- mrc_bait_refset.gr[inter[tokeep]]$gene_id
  inter_size <- length(inter_bait_refset.gid)
  shft <- seq(0,maxrange,step)
  refF.gr <- resize(refF.gr,window,"start",ignore.strand=T)
  refR.gr <- resize(shift(refF.gr,-window),window,"start",ignore.strand=T)
  myvec <- sapply(shft,
                  function(sshft){
                    sbstf.gr <- subsetByOverlaps(mrc_bait_refset.gr,trim(shift(refF.gr,sshft)))
                    sbstr.gr <- subsetByOverlaps(mrc_bait_refset.gr,trim(shift(refR.gr,-(sshft)))) # Added
                    sbst.gr <- unique(c(sbstf.gr,sbstr.gr)) # Added
                    if(length(sbst.gr)){
                      win_set <- sbst.gr$gene_id
                      a1 <- length(which(win_set %in% inter_bait_refset.gid))
                      a2 <- length(win_set) - a1
                      b1 <- inter_size - a1
                      b2 <- refset_size - (inter_size + length(win_set) - a1)
                      # print(matrix(c(a1,a2,b1,b2),2,2))
                      pv <- fisher.test(matrix(c(a1,a2,b1,b2),2,2),alternative = "greater")$p.value
                      perc <- 100*(a1/(a1+a2))
                      return(cbind(pv,perc))
                    }
                    else
                      return(cbind(1,NA))
                  })
  perc.df <- NULL
  myvec <- t(myvec)
  trsh <- sapply(1:(nrow(myvec)-3),
                 function(idx){
                   perc.df <<- rbind.data.frame(perc.df,cbind.data.frame(mean(myvec[c(idx:(idx+3)),2]),sd(myvec[c(idx:(idx+3)),2])))
                 })
  perc.df <- cbind(perc.df,shft[-((length(shft)-2):length(shft))]+window/2)
  colnames(perc.df) <- c("mean","sd","x")
  toPlot <- data.frame(x=(shft+window/2),y=myvec[,1])
  toPlot$padj <- p.adjust(toPlot$y,method = "BH",n=length(toPlot$y))
  toPlot$logpv <- -log10(toPlot$padj)
  toPlot$logpv <- ifelse(toPlot$logpv > 8,8,toPlot$logpv)
  toPlot$thr <- ifelse(toPlot$logpv > -log10(0.05),toPlot$logpv,0)
  toPlot$thr <- ifelse(toPlot$thr > 8,8,toPlot$thr)
  
  
  a <- ggplot(toPlot,aes(x=x,y=padj)) + ggtitle(main) + ylim(0,1) + geom_line() + ylab("Padj") + theme(plot.title = element_text(size = 10, face = "bold"))
  b <- ggplot(toPlot,aes(x=x,y=logpv)) + ggtitle(main) + geom_line() + ylim(0,8) + ylab("logPVadj")+ theme(plot.title = element_text(size = 10, face = "bold"))
  c <- ggplot(toPlot,aes(x=x,y=thr)) + ggtitle(main) + geom_line() + ylim(0,8) + ylab("logPVadjThresholded")+ theme(plot.title = element_text(size = 10, face = "bold"))
  d <- ggplot(toPlot,aes(x=x,y=y)) + ggtitle(main) + geom_point(size=0.05) + ylim(0,1) + ylab("raw_pv")+ theme(plot.title = element_text(size = 10, face = "bold"))
  e <- ggplot(perc.df, aes(x=x , y=mean) ) + geom_bar(position=position_dodge(.2),stat="identity",width = wdth)   +
    scale_y_continuous(name="% DR/UR", limits=c(0, 100)) + theme_bw() + rotate_x_text(45) + 
    geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),size=0.5,width=.15,position=position_dodge(.2))
  
  pdf(o_pdf)
  print(a);print(b);print(c);print(d);print(e)
  dev.off()
}
fisherRangePlotBothSided <- function(anchor,bait,inter,set,maxrange,step,window,out){
  # Output setup
  main <- paste0("Anchor=",anchor,"_Bait=",bait,"_Intersect=",inter,"_set=",set,"_maxrange=",maxrange,"KB_step=",step,"_Win=",window)
  out <- create(paste0(out,"/Anchor=",anchor,"/Bait=",bait,"/Range=",maxrange,"/Window=",window,"/"))
  wdth <- 25000/(maxrange/step)
  o_pdf <- paste0(out,"/",main,".pdf")
  # Set if for extra filter like "NoBeaf" at all on the global set variables
  anchor <- get(anchor);bait <- get(bait) & get(set);inter <- get(inter) & get(set);set <- get(set)
  # get the TSS.gr of the reference set
  anchor.gr <- mrc_tss.gr[anchor]
  # get filter set
  # Apply refset and filterset
  bait.gr <- mrc_tss.gr[bait]
  inter.gr <- mrc_tss.gr[inter]
  set.gr <- mrc_tss.gr[set]
  
  # Do the intersections
  bait_inter.gid <- intersect(bait.gr$gene_id,inter.gr$gene_id)
  nobait_inter.gid <- setdiff(inter.gr$gene_id,bait.gr$gene_id)
  bait_nointer.gid <- setdiff(bait.gr$gene_id,inter.gr$gene_id)
  nobait_nointer.gid <- setdiff(set.gr$gene_id,unique(c(bait.gr$gene_id,inter.gr$gene_id)))
  
  shft <- seq(0,maxrange,step)
  anchor.gr <- resize(anchor.gr,window,"start",ignore.strand=T)
  myvec <- sapply(shft,
                  function(sshft){
                    sbstp.gr <- subsetByOverlaps(set.gr,trim(shift(anchor.gr,sshft)))
                    sbstn.gr <- subsetByOverlaps(set.gr,trim(shift(anchor.gr,-(sshft+window)))) # Added
                    sbst.gr <- unique(c(sbstp.gr,sbstn.gr)) # Added
                    if(length(sbst.gr)){
                      win_set <- sbst.gr$gene_id
                      a1 <- length(which(bait_inter.gid %in% win_set))
                      a2 <- length(which(nobait_inter.gid %in% win_set))
                      b1 <- length(which(bait_nointer.gid %in% win_set))
                      b2 <- length(which(nobait_nointer.gid %in% win_set))
                      # print(matrix(c(a1,a2,b1,b2),2,2))
                      pv <- fisher.test(matrix(c(a1,a2,b1,b2),2,2),alternative = "greater")$p.value
                      perc <- 100*(a1/(a1+a2))
                      return(cbind(pv,perc))
                    }
                    else
                      return(cbind(1,NA))
                  })
  perc.df <- NULL
  myvec <- t(myvec)
  trsh <- sapply(1:(nrow(myvec)-3),
                 function(idx){
                   # print(idx)
                   perc.df <<- rbind.data.frame(perc.df,cbind.data.frame(mean(myvec[c(idx:(idx+3)),2]),sd(myvec[c(idx:(idx+3)),2])))
                 })
  perc.df <- cbind(perc.df,shft[-((length(shft)-2):length(shft))]+window/2)
  colnames(perc.df) <- c("mean","sd","x")
  toPlot <- data.frame(x=(shft+(window/2)),y=myvec[,1])
  toPlot$padj <- p.adjust(toPlot$y,method = "BH",n=length(toPlot$y))
  toPlot$logpv <- -log10(toPlot$padj)
  toPlot$logpv <- ifelse(toPlot$logpv > 8,8,toPlot$logpv)
  toPlot$thr <- ifelse(toPlot$logpv > -log10(0.05),toPlot$logpv,0)
  toPlot$thr <- ifelse(toPlot$thr > 8,8,toPlot$thr)
  
  
  a <- ggplot(toPlot,aes(x=x,y=padj)) + ggtitle(main) + ylim(0,1) + geom_line() + ylab("Padj") + theme(plot.title = element_text(size = 10, face = "bold"))
  b <- ggplot(toPlot,aes(x=x,y=logpv)) + ggtitle(main) + geom_line() + ylim(0,8) + ylab("logPVadj")+ theme(plot.title = element_text(size = 10, face = "bold"))
  c <- ggplot(toPlot,aes(x=x,y=thr)) + ggtitle(main) + geom_line() + ylim(0,8) + ylab("logPVadjThresholded")+ theme(plot.title = element_text(size = 10, face = "bold"))
  d <- ggplot(toPlot,aes(x=x,y=y)) + ggtitle(main) + geom_point(size=0.05) + ylim(0,1) + ylab("raw_pv")+ theme(plot.title = element_text(size = 10, face = "bold"))
  e <- ggplot(perc.df, aes(x=x , y=mean) ) + geom_bar(position=position_dodge(.2),stat="identity",width = wdth)   +
    scale_y_continuous(name="% DR/UR", limits=c(0, 100)) + theme_bw() + rotate_x_text(45) + 
    geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),size=0.5,width=.15,position=position_dodge(.2))
  
  pdf(o_pdf)
  print(a);print(b);print(c);print(d);print(e)
  dev.off()
}
fisherRPorientated <- function(anchor,focus,refForw,refRev,bait,inter,set,maxrange,step,window,out){
  # Output setup
  main <- paste0("Anchor=",anchor,"Focus=",focus,"_Bait=",bait,"_Intersect=",inter,"_set=",set,"_maxrange=",maxrange,"KB_step=",step,"_Win=",window)
  out <- create(paste0(out,"/Focus=",focus,"/Anchor=",anchor,"/Bait=",bait,"/Range=",maxrange,"/Window=",window,"/"))
  wdth <- 25000/(maxrange/step)
  o_pdf <- paste0(out,"/",main,".pdf")
  # Set if for extra filter like "NoBeaf" at all on the global set variables
  refForw <- get(refForw);refRev <- get(refRev);bait <- get(bait) & get(set);inter <- get(inter) & get(set); set <- get(set)
  # get the reference set
  refF.gr <- mrc_tss.gr[refForw]
  refR.gr <- mrc_tss.gr[refRev]
  refsetF.gr <- trim(resize(refF.gr,fix = "start",width=window,ignore.strand=T))
  refsetR.gr <- trim(resize(refR.gr,fix="start",width=window,ignore.strand=T))
  # get filter set
  # Apply refset and filterset
  bait.gr <- mrc_tss.gr[bait]
  inter.gr <- mrc_tss.gr[inter]
  set.gr <- mrc_tss.gr[set]
  
  # Do the intersections
  bait_inter.gid <- intersect(bait.gr$gene_id,inter.gr$gene_id)
  nobait_inter.gid <- setdiff(inter.gr$gene_id,bait.gr$gene_id)
  bait_nointer.gid <- setdiff(bait.gr$gene_id,inter.gr$gene_id)
  nobait_nointer.gid <- setdiff(set.gr$gene_id,unique(c(bait.gr$gene_id,inter.gr$gene_id)))
  
  shft <- seq(0,maxrange,step)
  myvec <- sapply(shft,
                  function(sshft){
                    sbstp.gr <- subsetByOverlaps(set.gr,trim(shift(refsetF.gr,sshft)))
                    sbstn.gr <- subsetByOverlaps(set.gr,trim(shift(refsetR.gr,-(sshft+window)))) # Added
                    sbst.gr <- unique(c(sbstp.gr,sbstn.gr)) # Added
                    if(length(sbst.gr)){
                      win_set <- sbst.gr$gene_id
                      a1 <- length(which(bait_inter.gid %in% win_set))
                      a2 <- length(which(nobait_inter.gid %in% win_set))
                      b1 <- length(which(bait_nointer.gid %in% win_set))
                      b2 <- length(which(nobait_nointer.gid %in% win_set))
                      # print(matrix(c(a1,a2,b1,b2),2,2))
                      pv <- fisher.test(matrix(c(a1,a2,b1,b2),2,2),alternative = "greater")$p.value
                      perc <- 100*(a1/(a1+a2))
                      return(cbind(pv,perc))
                    }
                    else
                      return(cbind(1,NA))
                  })
  perc.df <- NULL
  myvec <- t(myvec)
  trsh <- sapply(1:(nrow(myvec)-3),
                 function(idx){
                   # print(idx)
                   perc.df <<- rbind.data.frame(perc.df,cbind.data.frame(mean(myvec[c(idx:(idx+3)),2]),sd(myvec[c(idx:(idx+3)),2])))
                 })
  perc.df <- cbind(perc.df,shft[-((length(shft)-2):length(shft))]+window/2)
  colnames(perc.df) <- c("mean","sd","x")
  toPlot <- data.frame(x=(shft+(window/2)),y=myvec[,1])
  toPlot$padj <- p.adjust(toPlot$y,method = "BH",n=length(toPlot$y))
  toPlot$logpv <- -log10(toPlot$padj)
  toPlot$logpv <- ifelse(toPlot$logpv > 8,8,toPlot$logpv)
  toPlot$thr <- ifelse(toPlot$logpv > -log10(0.05),toPlot$logpv,0)
  toPlot$thr <- ifelse(toPlot$thr > 8,8,toPlot$thr)
  
  
  a <- ggplot(toPlot,aes(x=x,y=padj)) + ggtitle(main) + ylim(0,1) + geom_line() + ylab("Padj") + theme(plot.title = element_text(size = 10, face = "bold"))
  b <- ggplot(toPlot,aes(x=x,y=logpv)) + ggtitle(main) + geom_line() + ylim(0,8) + ylab("logPVadj")+ theme(plot.title = element_text(size = 10, face = "bold"))
  c <- ggplot(toPlot,aes(x=x,y=thr)) + ggtitle(main) + geom_line() + ylim(0,8) + ylab("logPVadjThresholded")+ theme(plot.title = element_text(size = 10, face = "bold"))
  d <- ggplot(toPlot,aes(x=x,y=y)) + ggtitle(main) + geom_point(size=0.05) + ylim(0,1) + ylab("raw_pv")+ theme(plot.title = element_text(size = 10, face = "bold"))
  e <- ggplot(perc.df, aes(x=x , y=mean) ) + geom_bar(position=position_dodge(.2),stat="identity",width = wdth)   +
    scale_y_continuous(name="% DR/UR", limits=c(0, 100)) + theme_bw() + rotate_x_text(45) + 
    geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),size=0.5,width=.15,position=position_dodge(.2))
  
  pdf(o_pdf)
  print(a);print(b);print(c);print(d);print(e)
  dev.off()
}

fisherTADb_XtoYkb <- function(refForw,refRev,bait,inter,set,X,Y){
  # Set if for extra filter like "NoBeaf" at all on the global set variables
  refForw <- get(refForw);refRev <- get(refRev);bait <- get(bait) & get(set);inter <- get(inter) & get(set); set <- get(set)
  refF.gr <- mrc_tss.gr[refForw]
  refR.gr <- mrc_tss.gr[refRev]
  
  refsetF.gr <- trim(resize(trim(shift(refF.gr,shift = X)),Y-X,"start",ignore.strand=T))
  refsetR.gr <- trim(resize(trim(shift(refR.gr,shift = -Y)),Y-X,"start",ignore.strand=T))
  
  # get filter set
  # Apply refset and filterset
  
  bait.gr <- mrc_tss.gr[bait]
  inter.gr <- mrc_tss.gr[inter]
  set.gr <- mrc_tss.gr[set]
  
  bait_inter.gid <- intersect(bait.gr$gene_id,inter.gr$gene_id)
  nobait_inter.gid <- setdiff(inter.gr$gene_id,bait.gr$gene_id)
  bait_nointer.gid <- setdiff(bait.gr$gene_id,inter.gr$gene_id)
  nobait_nointer.gid <- setdiff(set.gr$gene_id,unique(c(bait.gr$gene_id,inter.gr$gene_id)))
  
  sbstp.gr <- subsetByOverlaps(set.gr,refsetF.gr)
  sbstn.gr <- subsetByOverlaps(set.gr,refsetR.gr) # Added
  sbst.gr <- unique(c(sbstp.gr,sbstn.gr)) # Added
  sbst.gr <- unique(c(sbstp.gr,sbstn.gr)) # Added
  win_set <- sbst.gr$gene_id
  
  a1 <- length(which(bait_inter.gid %in% win_set))
  a2 <- length(which(nobait_inter.gid %in% win_set))
  b1 <- length(which(bait_nointer.gid %in% win_set))
  b2 <- length(which(nobait_nointer.gid %in% win_set))
  
  ft <- fisher.test(matrix(c(a1,a2,b1,b2),2,2),alternative = "greater")
  a_perc <- 100*(a1/(a1+a2))
  ft.pv <- ft$p.value # 0.93
  ft.fc <- ft$estimate # 0.97
  a_perc <- (a1*100)/(a1+a2) # 58%
  b_perc <- (b1*100)/(b1+b2) # 58%
  return(list(ft.pv=ft.pv,ft.fc=ft.fc,a_perc=a_perc,b_perc=b_perc,a_num=a1,b_num=b1))
  
}

fisherTADb_XtoYkb_BothSided <- function(anchor,bait,inter,set,X,Y){
  # Set if for extra filter like "NoBeaf" at all on the global set variables
  anchor <- get(anchor);bait <- get(bait) & get(set);inter <- get(inter) & get(set); set <- get(set)
  anchor.gr <- mrc_tss.gr[anchor]
  
  refsetF.gr <- trim(resize(trim(shift(anchor.gr,shift = X)),Y-X,"start",ignore.strand=T))
  refsetR.gr <- trim(resize(trim(shift(anchor.gr,shift = -Y)),Y-X,"start",ignore.strand=T))
  
  # get filter set
  # Apply refset and filterset
  
  bait.gr <- mrc_tss.gr[bait]
  inter.gr <- mrc_tss.gr[inter]
  set.gr <- mrc_tss.gr[set]
  
  bait_inter.gid <- intersect(bait.gr$gene_id,inter.gr$gene_id)
  nobait_inter.gid <- setdiff(inter.gr$gene_id,bait.gr$gene_id)
  bait_nointer.gid <- setdiff(bait.gr$gene_id,inter.gr$gene_id)
  nobait_nointer.gid <- setdiff(set.gr$gene_id,unique(c(bait.gr$gene_id,inter.gr$gene_id)))
  
  sbstp.gr <- subsetByOverlaps(set.gr,refsetF.gr)
  sbstn.gr <- subsetByOverlaps(set.gr,refsetR.gr) # Added
  sbst.gr <- unique(c(sbstp.gr,sbstn.gr)) # Added
  sbst.gr <- unique(c(sbstp.gr,sbstn.gr)) # Added
  win_set <- sbst.gr$gene_id
  
  a1 <- length(which(bait_inter.gid %in% win_set))
  a2 <- length(which(nobait_inter.gid %in% win_set))
  b1 <- length(which(bait_nointer.gid %in% win_set))
  b2 <- length(which(nobait_nointer.gid %in% win_set))
  
  print(matrix(c(a1,a2,b1,b2)))
  ft <- fisher.test(matrix(c(a1,a2,b1,b2),2,2),alternative = "greater")
  a_perc <- 100*(a1/(a1+a2))
  ft.pv <- ft$p.value # 0.93
  ft.fc <- ft$estimate # 0.97
  a_perc <- (a1*100)/(a1+a2) # 58%
  b_perc <- (b1*100)/(b1+b2) # 58%
  return(list(ft.pv=ft.pv,ft.fc=ft.fc,a_perc=a_perc,b_perc=b_perc,a_num=a1,b_num=b1))
  
}

nb_bait_on_inter_XtoYkb_BothSided <- function(anchor,bait,inter,set,X,Y){
  # Set if for extra filter like "NoBeaf" at all on the global set variables
  anchor <- get(anchor);bait <- get(bait) & get(set);inter <- get(inter) & get(set); set <- get(set)
  anchor.gr <- mrc_tss.gr[anchor]
  
  refsetF.gr <- trim(resize(trim(shift(anchor.gr,shift = X)),Y-X,"start",ignore.strand=T))
  refsetR.gr <- trim(resize(trim(shift(anchor.gr,shift = -Y)),Y-X,"start",ignore.strand=T))
  
  # get filter set
  # Apply refset and filterset
  
  bait.gr <- mrc_tss.gr[bait]
  inter.gr <- mrc_tss.gr[inter]
  set.gr <- mrc_tss.gr[set]
  
  bait_inter.gid <- intersect(bait.gr$gene_id,inter.gr$gene_id)
  
  sbstp.gr <- subsetByOverlaps(set.gr,refsetF.gr)
  sbstn.gr <- subsetByOverlaps(set.gr,refsetR.gr) # Added
  sbst.gr <- unique(c(sbstp.gr,sbstn.gr)) # Added
  sbst.gr <- unique(c(sbstp.gr,sbstn.gr)) # Added
  win_set <- sbst.gr$gene_id
  
  a1 <- length(which(bait_inter.gid %in% win_set))
  
  return(a1)
  
}
#   -----------------------------------------------------------------------


#  3D GENOMICS -------------------------------------------------------------

#                                               ---     APA     ---
# --------------------------------------------------------------------------------------------------------------#

# RETRIEVE THE Z-SCORE ---
# From an apa matrix, comppute the z-score and associated pvalue
# Input : apa_matrix (21x21)
# Output : zscore and associated pvalue
apa_zscore_comp <- function(apa_mat){
  central_bin <- mean(apa_mat[11,11])
  ctrl_bin <- apa_mat[16:21,1:6]
  zscore <- (central_bin - mean(ctrl_bin))/sd(ctrl_bin)
  pv <- (1-pnorm(zscore))
  return(c(`z-score`=zscore,pval=pv))
}

# THE FISHER TEST COMPARING APA MATRIX ---
# From 2 apa matrix, comppute the fisher statistic and associated pvalue between 2 features
# For instance between ORCacRad21 vs ORCacNoRad21
# the used value is the central bin
# Input : apa_matrix (21x21), apa_matrix (21x21)
# Output : odds ratio & associated pvalue
apa_fisher_comp <- function(apa_mat1,apa_mat2){
  # WITHH ALL THE MATRIX
  a_mat1 <- apa_mat1[11,11]
  b_mat1 <- sum(diag(apa_mat1)[-11])
  a_mat2 <- apa_mat2[11,11]
  b_mat2 <- sum(diag(apa_mat2)[-11])
  f_mat <- matrix(c(tst_feat=a_mat1,tst_noFeat=a_mat2,ctrl_feat=b_mat1,ctrl_noFeat=b_mat2),nrow=2,ncol=2)
  
  ft <- fisher.test(f_mat)
  ft.pv <- ft$p.value
  ft.fc <- ft$estimate
  return(list(`odds ratio`=ft.fc,pval=ft.pv,fmat=f_mat))
}

# THE POISSON TEST ---
# Poisson test between central tst bin and lower left corner controls
apa_poisson_test <- function(apa_mat){
  a <- apa_mat[11,11]
  b <- round(mean(apa_mat[c(16:21),c(1,6)]))
  c <- sum(apa_mat)
  return(c(pval=poisson.test(c(a,b),c(c,c),alternative = "greater")$p.value))
}

# THE BINOMIAL TEST ---
# Same as above but with a binomial test between central tst bin and lower left corner controls
apa_binom_test <- function(apa_path){
  a <- apa_mat[11,11]
  b <- round(mean(apa_mat[c(16:21),c(1,6)]))
  c <- sum(apa_mat)
  return(c(pval=binom.test(a,c,b/c,alternative = "greater")))
}

# EXTRACT APA MATRIX FROM APA.txt ---
# Input :  path to APA.txt
# Output : matrix of agregated contact (21x21)
apa_exctraction <- function(apa_path){
  require(gsubfn)
  mat_unparsed <- read.table(apa_path,h=F,sep=",")
  mat <- t(apply(mat_unparsed,1,function(x) as.numeric(gsubfn(".",list("["="","]"=""),x))))
  return(mat)
}



#                                      --- HIC TREATMENT ---
# --------------------------------------------------------------------------------------------------------------#

# CREATE COUPLE SEPARATED BY MAX DISTANCE OVER AN ENTIRE GENOME ---
# From a genomicRAnges object of ChIPseq data, create all possible couple
# within a certain min and max distance
# Input : GenomicRanges, maxDist, min Dist, Resolution of HiC
# Output : 2D bed (chr posL1 posR1 chr posL2 posR2)
createPksCpleByDist <- function(obj.gr,dstMax,dstMin,res=500){
  obj.gr$enzIDX <- ceiling((start(obj.gr))/res)
  obj.grl <- split(obj.gr,seqnames(obj.gr))
  dstMaxEnz <- ceiling(dstMax/500)
  dstMinEnz <- ceiling(dstMin/500)
  contact.intll <- lapply(obj.grl,
                          function(pksChrX.gr){
                            
                            contact.intl <- combn(pksChrX.gr,2,
                                                  function(cpl.grl){
                                                    cur_chr.str <- unique(seqnames(cpl.grl[1]))
                                                    if(((cpl.grl[2]$enzIDX - cpl.grl[1]$enzIDX) < dstMaxEnz) & ((cpl.grl[2]$enzIDX - cpl.grl[1]$enzIDX) > dstMinEnz)){
                                                      print("ok")
                                                      print(cbind(paste0(as.character(cur_chr.str)),
                                                                  start(cpl.grl[1]),
                                                                  end( cpl.grl[1]),
                                                                  paste0(as.character(cur_chr.str)),
                                                                  start( cpl.grl[2]),
                                                                  end( cpl.grl[2])))
                                                      return(cbind(paste0(as.character(cur_chr.str)),
                                                                   start(cpl.grl[1]),
                                                                   end( cpl.grl[1]),
                                                                   paste0(as.character(cur_chr.str)),
                                                                   start( cpl.grl[2]),
                                                                   end( cpl.grl[2])))
                                                    }
                                                    else return(F)
                                                  },simplify=F)
                            
                          })
  return(contact.intll)
}



# Create a 2dbed file that correspond to all potential interaction between one set of postion eg protein binding sites
# (int_fdir, path to GRanges object) in a given context (dom, path to a GRange object) at a given resolution (bin)
find_int_int_bydom <- function(int_fdir,bin,dom,o.dir){
  bin <- get(bin)
  dom.gr <- readRDS(get(dom))
  featName.str <- gsub("(.+)_(K|k)c167_dm3.(GR|gr).RData","\\1",basename(int_fdir))
  feat.gr <- resize(readRDS(int_fdir),1,fix="center") #resize the feature to do the overlap on domains without multi hits
  feat.gr$enzIDX <- floor((start(feat.gr))/500) # Bin are 0-based ie bp 500 for a prot will give bin 501-999 in the hic matrix
  df <- data.frame(findOverlapPairs(dom.gr,feat.gr))
  df$filt <- paste0(df$first.X.seqnames,"_",df$first.X.start)
  print(featName.str)
  feats_by_tad.lst <- by(df[,c("second.X.start","second.enzIDX","second.X.seqnames")],df$filt,function(x)x)
  # print(feats_by_tad.lst[1])
  lapply(feats_by_tad.lst,
         function(elm.vec){
           if(!(nrow(elm.vec)-1)){return()}
           combn(1:nrow(elm.vec),2,
                 function(idx.vec){
                   # print(elm.vec)
                   chr <- elm.vec[idx.vec[2],"second.X.seqnames"]
                   # We want to have 3x3 bins matrix size (bin=500bp)
                   if(elm.vec[idx.vec[2],"second.enzIDX"]-elm.vec[idx.vec[1],"second.enzIDX"] > 30){
                     s1 <- (floor(as.integer(elm.vec[idx.vec[1],"second.X.start"])/500)-((bin-1)/2))*500
                     e1 <- (floor(as.integer(elm.vec[idx.vec[1],"second.X.start"])/500)+((bin-1)/2))*500
                     s2 <- (floor(as.integer(elm.vec[idx.vec[2],"second.X.start"])/500)-((bin-1)/2))*500
                     e2 <- (floor(as.integer(elm.vec[idx.vec[2],"second.X.start"])/500)+((bin-1)/2))*500
                     
                     if(s1 < s2)
                       resu <- cbind(as.character(chr),s1,e1,as.character(chr),s2,e2)
                     else
                       resu <- cbind(as.character(chr),s2,e2,as.character(chr),s1,e1)
                     
                     write.table(resu,paste0(o.dir,"int=",featName.str,"_anc=",featName.str,"_dom=",dom,".bed"), quote=FALSE, sep="\t",
                                 row.names = FALSE,
                                 col.names = FALSE,
                                 append = T)
                   }
                 },simplify=F)
         })
}

# Create a 2dbed file that correspond to all potential interaction between two different set of postion eg 2 protein binding sites
# (int_fdir, anc_fdir path to GRanges object) in a given context (dom, path to a GRange object) at a given resolution (bin)
find_int_anc_bydom <- function(int_fdir,anc_fdir,bin,dom,o.dir,binS,smpl=NULL){
  # Limitating domain 
  dom.gr <- readRDS(get(dom))
  # Bin size 
  bin <- get(bin)
  #Interacting feature
  featNamei.str <- gsub("(.+).(ints|ancs).gr.hg19.rds","\\1",basename(int_fdir))
  feati.gr <- resize(readRDS(int_fdir),1,fix="center") #resize the feature to do the overlap on domains without multi hits
  feati.gr$enzIDX <- floor((start(feati.gr))/binS)
  dfi <- data.frame(findOverlapPairs(dom.gr,feati.gr))
  dfi$filt <- paste0(dfi$first.seqnames,"_",dfi$first.start)
  dfi <- dfi[,c("second.enzIDX","filt")]
  print(featNamei.str)
  #Anchoring feature
  featNamea.str <- gsub("(.+).(ints|ancs).gr.hg19.rds","\\1",basename(anc_fdir))
  feata.gr <- resize(readRDS(anc_fdir),1,fix="center") #resize the feature to do the overlap on domains without multi hits
  feata.gr$enzIDX <- floor((start(feata.gr))/binS)
  print(featNamea.str)
  dfa <- data.frame(findOverlapPairs(dom.gr,feata.gr))
  dfa$filt <- paste0(dfa$first.seqnames,"_",dfa$first.start)
  dfa <- dfa[,c("second.enzIDX","filt")]
  mia <- merge(dfi,dfa,by = "filt")
  colnames(mia) <- c("tadpos","enzi","enza")
  # Remove couples that are < 15kb/binS fragment away
  mia <- mia[-which(abs(mia$enzi-mia$enza)<(15000/binS)),]
  # Create couples
  print(nrow(mia))
  if(nrow(mia) < 2) next()
  if(smpl & (nrow(mia) > smpl))
    mia <- mia[sample(nrow(mia),smpl),]
  trash <- apply(mia, 1, 
                 function(lne){
                   chr <- strsplit(as.character(lne[1]),"_")[[1]][1]
                   s1 <- (as.integer(lne[2])-((bin-1)/2))*binS 
                   e1 <- (as.integer(lne[2])+((bin-1)/2))*binS
                   s2 <- (as.integer(lne[3])-((bin-1)/2))*binS
                   e2 <- (as.integer(lne[3])+((bin-1)/2))*binS
                   
                   
                   if(s1 < s2){
                     resu <- cbind(as.character(chr),s1,e1,as.character(chr),s2,e2)
                     write.table(resu,paste0(o.dir,"int=",featNamei.str,"_anc=",featNamea.str,"_dom=",dom,".bed"), quote=FALSE, sep="\t",
                                 row.names = FALSE,
                                 col.names = FALSE,
                                 append = T)
                   }
                   
                   else{
                     resu <- cbind(as.character(chr),s2,e2,as.character(chr),s1,e1)
                     write.table(resu,paste0(o.dir,"int=",featNamei.str,"_anc=",featNamea.str,"_dom=",dom,"_toRev.bed"), quote=FALSE, sep="\t",
                                 row.names = FALSE,
                                 col.names = FALSE,
                                 append = T)
                   }
                   
                   
                 })
  print("next")
}
#   -----------------------------------------------------------------------

#  GRAPHICS ---------------------------------------------------------------
blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

#                                         --- AGREGATION PLOTS ---
# --------------------------------------------------------------------------------------------------------------#

# UNNORMALIZED AGREGGATION PLOT ---
unorm_agreg_plot <- function(apa_mat,o_name){
  l1 <- apa_mat[11,]/mean(apa_mat[11,])
  lctrl <- NULL
  for(i in c(1:4,18:21)){
    lctrl <- rbind(lctrl,apa_mat[i,]/mean(apa_mat[i,]))
  }
  jpeg(paste0(o_name,"_unnorm_agrPlot.jpeg"),width=5, height=5, units="in", res=500)
  plot(x=-10:10,y=l1,col="red",type="l",lwd=1.7)
  for(i in 1:nrow(lctrl)){
    lines(x=-10:10,y=lctrl[i,],col="grey",type="l")
  }
  dev.off()
}
unorm_agreg_plot_paper <- function(apa_mat,o_name){
  l1 <- apa_mat[11,]/mean(apa_mat[11,])
  lctrl <- NULL
  for(i in c(1:4,18:21)){
    lctrl <- rbind(lctrl,apa_mat[i,]/mean(apa_mat[i,]))
  }
  jpeg(paste0(o_name,"_unnorm_agrPlot_paper.jpeg"),width=5, height=5, units="in", res=500)
  par(mar=c(0, 0, 0, 0),oma=c(0,0,0,0))
  
  plot(x=-10:10,y=l1,col="red",type="l",lwd=2,ylab='',xlab='',cex.axis=1.6,font=2,axes=FALSE,bty="n")
  for(i in 1:nrow(lctrl)){
    lines(x=-10:10,y=lctrl[i,],col="grey",type="l",lwd=1.4)
  }
  box("outer",lwd=3)
  
  dev.off()
}

# NORMALIZED AGREGGATION PLOT ---
norm_agreg_plot <- function(apa_mat,o_name){
  l1 <- apa_mat[11,]/mean(apa_mat[11,])
  lctrl <- NULL
  for(i in c(1:4,18:21)){
    lctrl <- rbind(lctrl,apa_mat[i,]/mean(apa_mat[i,]))
  }
  m_tot <- rbind(lctrl,l1)
  fac_norm <- apply(m_tot,2,mean)
  l1norm <- l1 / fac_norm
  lctrlNorm <- t(apply(lctrl,1,function(x) x/fac_norm))
  jpeg(paste0(o_name,"_norm_agrPlot.jpeg"))
  par(mfrow=c(2,1))
  plot(x=-10:10,y=l1norm,col="red",type="l",lwd=1.7,ylim=c(0.8,1.2))
  for(i in 1:nrow(lctrlNorm)){
    lines(x=-10:10,y=lctrlNorm[i,],col="grey",type="l")
  }
  
  mlctrl <- apply(lctrlNorm,2,mean)
  plot(x=-10:10,y=l1norm,col="red",type="l",lwd=1.7,ylim=c(0.8,1.2))
  lines(x=-10:10,y=mlctrl,col="grey",type="l")
  dev.off()
}

norm_agreg_plot_paper <- function(apa_mat,o_name){
  l1 <- apa_mat[11,]/mean(apa_mat[11,])
  lctrl <- NULL
  for(i in c(1:4,18:21)){
    lctrl <- rbind(lctrl,apa_mat[i,]/mean(apa_mat[i,]))
  }
  m_tot <- rbind(l1,lctrl)
  meanT <- (apply(m_tot,2,mean))
  sumT <- sum(meanT)
  l1norm <- l1/(sum(l1)/sumT)/meanT
  lctrlNorm <- t(apply(lctrl,1,function(x) (x/(sum(x)/sumT))/meanT))
  
  # x=-10:10
  # plot(x=x,y=l1norm,col="black",type="l",lwd=1.7,ylim=c(0.5,1.5),ylab='',xlab='',cex.axis=1.6,font=2,axes=FALSE,bty="n")
  # polygon(c(x[1],x,x[length(x)]),c(1,l1norm,1),col=adjustcolor("black",alpha.f=0.5),panel.first=abline(h=1,lty=3),border=NA)
  #
  
  jpeg(paste0(o_name,"_norm_agrPlot.jpeg"),width=5, height=5, units="in", res=500)
  par(mar=c(0, 0, 0, 0),oma=c(0,0,0,0))
  x=-10:10
  plot(x=x,y=l1norm,col="black",type="l",lwd=1.7,ylim=c(0.5,1.5),ylab='',xlab='',cex.axis=1.6,font=2,axes=FALSE,bty="n")
  polygon(c(x[1],x,x[length(x)]),c(1,l1norm,1),col=adjustcolor("black",alpha.f=0.5),panel.first=abline(h=1,lty=3),border=NA)
  for(i in 1:nrow(lctrlNorm)){
    lines(x=-10:10,y=lctrlNorm[i,],col="black",type="l",lty=6)
    polygon(c(x,rev(x)),c(lctrlNorm[i,],rep(1,21)),col=adjustcolor("grey",alpha.f=0.5),border=NA)
  }
  box("outer",lwd=3)
  dev.off()
  
  jpeg(paste0(o_name,"_norm_agrPlot_mean_paper.jpeg"),width=5, height=5, units="in", res=500)
  par(mar=c(0, 0, 0, 0),oma=c(0,0,0,0))
  mlctrl <- apply(lctrlNorm,2,mean)
  plot(x=-10:10,y=l1norm,col="red",type="l",lwd=2,ylim=c(0.8,1.2),ylab='',xlab='',cex.axis=1.6,font=2,axes=FALSE,bty="n")
  lines(x=-10:10,y=mlctrl,col="grey",type="l",lwd=1.4)
  abline(v=0,lty=3,col="darkgrey")
  box("outer",lwd=3)
  dev.off()
}

norm_agreg_plot_w20 <- function(apa_mat,o_name){
  l1 <- diag(apa_mat)[-c(1:10,31:40)]
  lctrl <- NULL
  for(i in seq(2,18,by=2)){
    dg_upper <- diag(apa_mat[-c(41:(42-i)),-c(1:i)])
    dg_lower <- diag(apa_mat[-c(1:i),-c(41:(42-i))])
    print(dg_lower[-c(1:(10-i/2),(length(dg_lower)-length(c(1:(10-i/2)))):length(dg_lower))])
    lctrl <- rbind(lctrl,dg_lower[-c(1:(10-i/2),((length(dg_lower)-length(c(1:(10-i/2))))+1):length(dg_lower))])
    lctrl <- rbind(lctrl,dg_upper[-c(1:(10-i/2),((length(dg_upper)-length(c(1:(10-i/2))))+1):length(dg_upper))])
  }
  m_tot <- rbind(lctrl,l1)
  fac_norm <- apply(m_tot,2,mean)
  l1norm <- l1 / fac_norm
  lctrlNorm <- t(apply(lctrl,1,function(x) x/fac_norm))
  jpeg(paste0(o_name,"_norm_agrPlot_w20.jpeg"))
  par(mfrow=c(2,1))
  x=-10:10
  plot(x=x,y=l1norm,col="black",type="l",lwd=1.7,ylim=c(0.5,1.5))
  polygon(c(x[1],x,x[length(x)]),c(1,l1norm,1),col=adjustcolor("orange",alpha.f=0.5),panel.first=abline(h=1,lty=3),border=NA)
  
  for(i in 1:nrow(lctrlNorm)){
    lines(x=-10:10,y=lctrlNorm[i,],col="black",type="l",lty=6)
    polygon(c(x,rev(x)),c(lctrlNorm[i,],rep(1,21)),col=adjustcolor("grey",alpha.f=0.5),border=NA)
  }
  
  mlctrl <- apply(lctrlNorm,2,mean)
  plot(x=-10:10,y=l1norm,col="red",type="l",lwd=1.7,ylim=c(0.5,2))
  lines(x=-10:10,y=mlctrl,col="grey",type="l")
  dev.off()
}

norm_agreg_plot_w20_paper <- function(apa_mat,o_name){
  l1 <- diag(apa_mat)[-c(1:10,31:40)]
  lctrl <- NULL
  for(i in seq(6,18,by=2)){
    dg_upper <- diag(apa_mat[-c(41:(42-i)),-c(1:i)])
    dg_lower <- diag(apa_mat[-c(1:i),-c(41:(42-i))])
    print(dg_lower[-c(1:(10-i/2),(length(dg_lower)-length(c(1:(10-i/2)))):length(dg_lower))])
    lctrl <- rbind(lctrl,dg_lower[-c(1:(10-i/2),((length(dg_lower)-length(c(1:(10-i/2))))+1):length(dg_lower))])
    lctrl <- rbind(lctrl,dg_upper[-c(1:(10-i/2),((length(dg_upper)-length(c(1:(10-i/2))))+1):length(dg_upper))])
  }
  
  m_tot <- rbind(l1,lctrl)
  meanT <- (apply(m_tot,2,mean))
  sumT <- sum(meanT)
  l1norm <- l1/(sum(l1)/sumT)/meanT
  lctrlNorm <- t(apply(lctrl,1,function(x) (x/(sum(x)/sumT))/meanT))
  
  # x=-10:10
  # plot(x=x,y=l1norm,col="black",type="l",lwd=1.7,ylim=c(0.5,1.5),ylab='',xlab='',cex.axis=1.6,font=2,axes=FALSE,bty="n")
  # polygon(c(x[1],x,x[length(x)]),c(1,l1norm,1),col=adjustcolor("black",alpha.f=0.5),panel.first=abline(h=1,lty=3),border=NA)
  
  
  jpeg(paste0(o_name,"_norm_agrPlot_w20.jpeg"),width=5, height=5, units="in", res=500)
  par(mar=c(0, 0, 0, 0),oma=c(0,0,0,0))
  x=-10:10
  plot(x=x,y=l1norm,col="black",type="l",lwd=1.7,ylim=c(0.5,1.5),ylab='',xlab='',cex.axis=1.6,font=2,axes=FALSE,bty="n")
  polygon(c(x[1],x,x[length(x)]),c(1,l1norm,1),col=adjustcolor("black",alpha.f=0.5),panel.first=abline(h=1,lty=3),border=NA)
  max_y <- apply(lctrlNorm,2,max)
  min_y <- apply(lctrlNorm,2,min)
  lines(x=-10:10,y=max_y,col=adjustcolor("red",alpha.f=0.3),type="l",lty=3)
  lines(x=-10:10,y=min_y,col=adjustcolor("red",alpha.f=0.3),type="l",lty=3)
  polygon(c(x,rev(x)),c(max_y,rev(min_y)),col=adjustcolor("grey",alpha.f=0.6),border=NA)
  for(i in 1:nrow(lctrlNorm)){
    lines(x=-10:10,y=lctrlNorm[i,],col=adjustcolor("black",alpha.f=0.6),type="l",lty=6)
    # polygon(c(x,rev(x)),c(lctrlNorm[i,],rep(1,21)),col=adjustcolor("grey",alpha.f=0.3),border=NA)
  }
  # for(i in 1:nrow(lctrlNorm)){
  # lines(x=-10:10,y=lctrlNorm[i,],col="black",type="l",lty=6)
  # polygon(c(x,rev(x)),c(lctrlNorm[i,],rep(1,21)),col=adjustcolor("grey",alpha.f=0.1),border=NA)
  # }
  box("outer",lwd=3)
  dev.off()
}

norm_agreg_plot_paper <- function(apa_mat,o_name){
  l1 <- apa_mat[11,]/mean(apa_mat[11,])
  lctrl <- NULL
  for(i in c(1:4,18:21)){
    lctrl <- rbind(lctrl,apa_mat[i,]/mean(apa_mat[i,]))
  }
  m_tot <- rbind(lctrl,l1)
  fac_norm <- apply(m_tot,2,mean)
  l1norm <- l1 / fac_norm
  lctrlNorm <- t(apply(lctrl,1,function(x) x/fac_norm))
  jpeg(paste0(o_name,"_norm_agrPlot_paper.jpeg"),width=5, height=5, units="in", res=500)
  par(mar=c(0, 0, 0, 0),oma=c(0,0,0,0))
  plot(x=-10:10,y=l1norm,col="red",type="l",lwd=2,ylim=c(0.8,1.2),ylab='',xlab='',cex.axis=1.6,font=2,axes=FALSE,bty="n")
  for(i in 1:nrow(lctrlNorm)){
    lines(x=-10:10,y=lctrlNorm[i,],col="grey",type="l",lwd=1.4)
  }
  abline(v=0,lty=3,col="darkgrey")
  box("outer",lwd=3)
  
  dev.off()
  
  jpeg(paste0(o_name,"_norm_agrPlot_mean_paper.jpeg"),width=5, height=5, units="in", res=500)
  par(mar=c(0, 0, 0, 0),oma=c(0,0,0,0))
  mlctrl <- apply(lctrlNorm,2,mean)
  plot(x=-10:10,y=l1norm,col="red",type="l",lwd=2,ylim=c(0.8,1.2),ylab='',xlab='',cex.axis=1.6,font=2,axes=FALSE,bty="n")
  lines(x=-10:10,y=mlctrl,col="grey",type="l",lwd=1.4)
  abline(v=0,lty=3,col="darkgrey")
  box("outer",lwd=3)
  dev.off()
}

#DIAGONAL PLOT ---
diag_agreg_plot_paper <- function(apa_mat,o_name){
  toPl <- diag(apa_mat)
  toPl <- toPl/max(toPl)
  jpeg(paste0(o_name,"_diag_agrPlot_paper.jpeg"),width=5, height=5, units="in", res=500)
  par(mar=c(0, 0, 0, 0),oma=c(0,0,0,0))
  plot(x=-10:10,y=toPl,col="red",type="l",lwd=2,ylim=c(0,1.2),ylab='',xlab='',cex.axis=1.6,font=2,axes=FALSE,bty="n")
  abline(v=0,lty=3,col="darkgrey")
  box("outer",lwd=3)
  dev.off()
}




#                                         --- BARPLOTS ---
# --------------------------------------------------------------------------------------------------------------#

simplePlot <- function(o_name,dat){
  p <- ggplot(dat,aes(x=Row,y=value)) +
    geom_bar(stat = "identity",position = "dodge")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(filename=paste0(o_name,".pdf"),width = 10,plot=p)
  ggsave(filename=paste0(o_name,".jpg"),width = 10,plot=p)
  
}

# TOGPLOT BAR PLOT COMPARISON
toGplotPaper <- function(o_name,i_vec,dat){
  
  colnames(dat) <- c("TAD","Var","COEF","upper","lower")
  dat <- as.data.frame(dat)
  dat$TAD <- as.factor(dat$TAD)
  #dat$Var must be the same as fill=Var,levels = ...
  dat$Var <- i_vec
  
  # print(i_vec)
  p <- ggplot(dat, aes(x=TAD ,y=COEF,fill=factor(Var,levels = dat$Var <- i_vec)),environment = environment()) +
    geom_bar(colour="black",width=.5, size=0.9,position=position_dodge(.7), stat="identity") + ylim(0,0.15) +
    geom_errorbar(aes(ymin=lower, ymax=upper),width=.05, size=0.7, position=position_dodge(.7))+ # ymin=as.numeric(COEF)-as.numeric(sd)
    scale_fill_manual(values=c("#e10000","darkgrey","#e10000","darkgrey","#e10000","darkgrey","#e10000","darkgrey","#e10000","darkgrey","#e10000","darkgrey"))  +
    
    theme(panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.margin = unit(c(0, 0, 0, 0), "cm"),
          plot.margin = rep(unit(0,"null"),4),
          axis.ticks=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.length = unit(0,"null"),
          axis.ticks.margin = unit(0,"null"),
          plot.background = element_rect(fill = "transparent",colour="black",size=1),
          # plot.margin = unit(c(-1, -1.2, -1.2, -1.5), "cm"),  # Edited code
          legend.position = 'none') +
    labs(x=NULL, y=NULL)
  
  ggsave(filename=paste0(o_name,".pdf"),plot=p)
  ggsave(filename=paste0(o_name,".jpg"), plot=p)
  
}

#TOGPLOT LOGFC ODDRATIO PLOT COMPARISON
toGplotEnrichment <- function(o_name,dat,paper=F,wd=10){
  
  p <- ggplot(dat,aes(x = Input,y = value,fill=Var)) +
    geom_bar(stat = "identity",position = "dodge") +
    geom_errorbar(aes(ymin=sdlo, ymax=sdup),width=.05, size=0.7, position=position_dodge(.9))+
    scale_fill_manual(values = c("#f0adad","#f08181","#ec5f5f","#e52a2a"))
  if(paper){
    p <- p +    theme(panel.background=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      panel.margin = unit(c(0, 0, 0, 0), "cm"),
                      plot.margin = rep(unit(0,"null"),4),
                      axis.ticks=element_blank(),
                      axis.text.x=element_blank(),
                      axis.text.y=element_blank(),
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank(),
                      axis.ticks.length = unit(0,"null"),
                      axis.ticks.margin = unit(0,"null"),
                      plot.background = element_rect(fill = "transparent",colour="black",size=1),
                      # plot.margin = unit(c(-1, -1.2, -1.2, -1.5), "cm"),  # Edited code
                      legend.position = 'none') +
      labs(x=NULL, y=NULL)
  }
  
  
  ggsave(filename=paste0(o_name,".pdf"),width = wd,plot=p)
  ggsave(filename=paste0(o_name,".jpg"), width = wd,plot=p)
  
}

# FISHER BARPLOT RESUME (PERCENTAGE OF THE INTERSECTION)
# df columns => x = col y=row fill=logpv, pv, odds, perc
fisherBplot <- function(df,title,out_f){
  o_tiff <- paste0(out_f,"_BP_PAPER.tiff")
  o_perc <- paste0(out_f,"_BP_PERC.pdf")
  p <- ggbarplot(df, x="col" , y="perc" , position = position_dodge(0), fill="type", palette=c(alpha("dodgerblue",0.3),alpha("firebrick1",0.3)),width = 0.3)   +
    scale_y_continuous(name="% DR/UR", limits=c(0, 100)) + theme_bw() + rotate_x_text(45)
  p <- facet(p,facet.by="row",ncol=length(unique(df$col)),axis.text.x=element_text(angle = 90, hjust = 1))
  ggexport(p,filename = o_perc,height = 2,width=10)
  
  tiff(o_tiff,units="in", width=5, height=5, res=300)
  p <- ggbarplot(df, x="row" , y="perc" , position = position_dodge(0), fill="type", palette=c(alpha("dodgerblue",0.3),alpha("firebrick1",0.3)),width = 0.3)   +
    scale_y_continuous(name="% DR/UR", limits=c(0, 100))+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank(),legend.position = "none")
  
  p <- facet(p,facet.by="col",ncol=length(unique(df$col)))
  print(p)
  dev.off()
  
}
# FISHER BARPLOT RESUME WITH MEAN AND SD OF 1/2,5/6, and 9/10 decile (PERCENTAGE OF THE INTERSECTION)
# df columns => x = col y=row fill=logpv, pv, odds, perc
fisherBSDplot <- function(df,title,out_f){
  o_tiff <- paste0(out_f,"_BP_SD_PAPER.tiff")
  o_perc <- paste0(out_f,"_BP_SD_PERC.pdf")
  p <- ggplot(df, aes(x=type , y=mean,fill=type) ) + geom_bar(position = "dodge",stat="identity",width = 0.3)   +
    scale_y_continuous(name="% DR/UR", limits=c(0, 100)) +facet_grid(~row)+ theme_bw() + rotate_x_text(45)
  p <- p + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),size=0.5,
                         width=.15,position=position_dodge(.9)) +  facet_grid(~row)
  ggexport(p,filename = o_perc,height = 2,width=10)
  
  
  tiff(o_tiff,units="in", width=5, height=5, res=300)
  p <-  ggplot(df, aes(x=type , y=mean,fill=type) ) + geom_bar(position = "dodge",stat="identity",width = 0.3)   +
    scale_y_continuous(name="% DR/UR", limits=c(0, 100))+ theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_blank(),
          panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          
          legend.position = "none") + facet_grid(~row)
  
  p <- p + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),size=0.5,
                         width=.15,position=position_dodge(.9)) +  facet_grid(~row)
  print(p)
  dev.off()
  
}


#                                         --- HEATMAP ---
# --------------------------------------------------------------------------------------------------------------#

# APA DIFFERENTIAL HEATMAP ---
diff_heatmap_plot <- function(apa1.mtx,apa2.mtx,qInf=0.025,qSup=0.975,o_name){
  apa1_norm.mtx <- apa1.mtx * (mean(apa2.mtx)/mean(apa1.mtx))
  diff_apa.mtx <- apa1_norm.mtx / apa2.mtx
  
  qInf.val <- quantile(diff_apa.mtx,qInf)
  qSup.val <- quantile(diff_apa.mtx,qSup)
  
  diff_apa.mtx[diff_apa.mtx<qInf.val] <- qInf.val
  diff_apa.mtx[diff_apa.mtx>qSup.val] <- qSup.val
  
  my_palette <- colorRampPalette(c("white","firebrick1"))(n = 256)
  jpeg(paste0(o_name,"_diff_heatmap_qInf=",qInf,"_qSup=",qSup,".jpeg"))
  heatmap(diff_apa.mtx,Rowv=NA,Colv=NA,keep.dendro=F,scale="none",revC=T,col=my_palette)
  dev.off()
}

# APA CLASSICAL HEATMAP ---
classic_heatmap_plot <- function(apa1.mtx,qInf=0.025,qSup=0.975,o_name){
  
  qInf.val <- quantile(apa1.mtx,qInf)
  qSup.val <- quantile(apa1.mtx,qSup)
  
  apa1.mtx[apa1.mtx<qInf.val] <- qInf.val
  apa1.mtx[apa1.mtx>qSup.val] <- qSup.val
  
  my_palette <- colorRampPalette(c("white","firebrick1"))(n = 256)
  jpeg(paste0(o_name,"_classic_heatmap.jpeg"))
  heatmap(apa1.mtx,Rowv=NA,Colv=NA,keep.dendro=F,scale="none",revC=T,col=my_palette)
  dev.off()
  
}

# FISHER GGPLOT
# df columns => x = col y=row fill=logpv, pv, odds, perc
fisherGplot <- function(df,title,out_f,ratio=.15){
  p <- ggplot(df, aes(x = col, y = row)) +
    geom_tile(aes(fill = logpv),col="black") +
    # ggtitle(title) +
    theme_minimal() +
    coord_equal(ratio=ratio) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top") +
    labs(x="",y="") +
    theme(axis.text.x=element_text(size=4, vjust=1, hjust=.5,
                                   margin=margin(1,0,0,0)),
          axis.text.y=element_text(size=10, margin=margin(0,-3,0,0)),
          # panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    geom_text(aes(label = paste0(pv,"\n",odds)),col="black",size=2) +
    scale_fill_gradient2(low = "firebrick",mid="white", high = "dodgerblue") +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10,col="black",ticks = FALSE))
  ggexport(p,filename=paste0(out_f,"_MAT_PVAL.pdf"))
  # Percentage % colored scale
  p <- ggplot(df, aes(x = col, y = row)) +
    geom_tile(aes(fill = perc),col="black") +
    # ggtitle(title) +
    theme_minimal() +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top") +
    coord_equal(ratio=ratio) +
    labs(x="",y="") +
    theme(axis.text.x=element_text(size=4, vjust=1, hjust=.5,
                                   margin=margin(1,0,0,0)),
          axis.text.y=element_text(size=10, margin=margin(0,-3,0,0)),
          # panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    geom_text(aes( label = paste0(round(perc)," % \n",lgt)),col="black",size=2) +
    scale_fill_gradient(low = "white", high = "orange") +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10,col="black",ticks = FALSE))
  ggexport(p,filename=paste0(out_f,"_MAT_PERC.pdf"))
  
  tiff(paste0(out_f,"_MAT_PVAL_PAPER.tiff"),units="in", width=5, height=5, res=300)
  p <- ggplot(df, aes(x = col, y = row)) +
    geom_tile(aes(fill = logpv),col="black") +
    # ggtitle(title) +
    theme_minimal() +
    coord_equal(ratio=ratio) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top") +
    labs(x="",y="") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_blank(),
          legend.title = element_blank()) +
    scale_fill_gradient2(low = "firebrick",mid="white", high = "dodgerblue") +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10,col="black",ticks = FALSE))
  print(p)
  dev.off()
  
  tiff(paste0(out_f,"_MAT_PERC_PAPER.tiff"),units="in", width=5, height=5, res=300)
  p <- ggplot(df, aes(x = col, y = row)) +
    geom_tile(aes(fill = perc),col="black") +
    # ggtitle(title) +
    theme_minimal() +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top") +
    coord_equal(ratio=ratio) +
    labs(x="",y="") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_blank(),
          legend.title = element_blank()) +
    scale_fill_gradient(low = "white", high = "orange") +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10,col="black",ticks = FALSE))
  print(p)
  dev.off()
  
  
}


# fisherGBplot <- function(df,bpdf,title,out_f,ratio=.15){
fisherGBplot <- function(df,title,out_f,ratio=.15){
  l_mat <- matrix(c(1,2),nrow=2,byrow=T)
  mylay <- layout(mat=l_mat,
                  heights = c(2,10))
  
  # MAT PVAL PLOT ---
  p <- ggplot(df, aes(x = col, y = row)) +
    geom_tile(aes(fill = logpv),col="black") +
    theme_minimal() +
    coord_equal(ratio=ratio) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top") +
    labs(x="",y="") +
    theme(axis.text.x=element_text(size=10, vjust=1, hjust=.5,
                                   margin=margin(1,0,0,0)),
          axis.text.y=element_text(size=10, margin=margin(0,-3,0,0)),
          plot.margin = unit(c(0,1,0,1), "cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    geom_text(aes(label = paste0(pv,"\n",odds)),col="black",size=2) +
    
    # scale_fill_gradient2(low = "firebrick",mid="white", high = "dodgerblue") +
    scale_fill_gradientn(colours=c("firebrick","white","dodgerblue"),
                         breaks=c(-2,0,2),labels=c("-2","0","+2"),
                         limits=c(-2,2))
  # guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10,col="black",ticks = FALSE))
  # bp <- ggplot(na.omit(bpdf), aes(x=x, y=y)) +
  #   geom_boxplot()+theme_minimal()+
  #   theme(panel.grid.minor = element_blank(),
  #         plot.margin = unit(c(0,1,0,1), "cm"),
  #         panel.grid.major = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.text.x=element_blank(),
  #         axis.ticks = element_blank(),
  #         axis.title.x=element_blank(),
  #         axis.line.x = element_blank())
  fp <- ggarrange(p,nrow=1,heights = c(0.3,0.7))
  # fp <- ggarrange(bp,p,nrow=2,heights = c(0.3,0.7))
  ggexport(fp,filename=paste0(out_f,"_MAT_PVAL_BOXPL.pdf"))
  # 
  # #  MAT PERCENTAGE
  # p <- ggplot(df, aes(x = col, y = row)) +
  #   geom_tile(aes(fill = perc),col="black") +
  #   theme_minimal() +
  #   scale_y_discrete(expand=c(0,0))+
  #   scale_x_discrete(position = "top") +
  #   coord_equal(ratio=ratio) +
  #   labs(x="",y="") +
  #   theme(axis.text.x=element_text(size=10, vjust=1, hjust=.5,
  #                                  margin=margin(1,0,0,0)),
  #         axis.text.y=element_text(size=10, margin=margin(0,-3,0,0)),
  #         panel.grid.minor = element_blank(),
  #         panel.grid.major = element_blank()) +
  #   geom_text(aes( label = paste0(round(perc)," % \n",lgt)),col="black",size=2) +
  #   scale_fill_gradient(low = "white", high = "orange") +
  #   guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10,col="black",ticks = FALSE))
  # bp <- ggplot(na.omit(bpdf), aes(x=x, y=y)) +
  #   geom_boxplot()+theme_minimal()+
  #   theme(panel.grid.minor = element_blank(),
  #         plot.margin = unit(c(0,1,0,1), "cm"),
  #         panel.grid.major = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.text.x=element_blank(),
  #         axis.ticks = element_blank(),
  #         axis.title.x=element_blank(),
  #         axis.line.x = element_blank())
  # fp <- ggarrange(bp,p,nrow=2,heights = c(0.3,0.7))
  # ggexport(fp,filename=paste0(out_f,"_MAT_PERC_BOXPL.pdf"))
  # 
  
  # TIFF PART
  tiff(paste0(out_f,"_MAT_PVAL_PAPER.tiff"),units="in", width=5, height=5, res=300)
  p <- ggplot(df, aes(x = col, y = row)) +
    geom_tile(aes(fill = logpv),col="black") +
    # ggtitle(title) +
    theme_minimal() +
    coord_equal(ratio=ratio) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top") +
    labs(x="",y="") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.text = element_blank(),
          legend.title = element_blank()) +
    scale_fill_gradient2(low = "firebrick",mid="white", high = "dodgerblue") +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10,col="black",ticks = FALSE))
  print(p)
  dev.off()
  
  # tiff(paste0(out_f,"_MAT_PERC_PAPER.tiff"),units="in", width=5, height=5, res=300)
  # p <- ggplot(df, aes(x = col, y = row)) +
  #   geom_tile(aes(fill = perc),col="black") +
  #   # ggtitle(title) +
  #   theme_minimal() +
  #   scale_y_discrete(expand=c(0,0))+
  #   scale_x_discrete(position = "top") +
  #   coord_equal(ratio=ratio) +
  #   labs(x="",y="") +
  #   theme(axis.title.x=element_blank(),
  #         axis.text.x=element_blank(),
  #         axis.title.y=element_blank(),
  #         axis.text.y = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.grid.major = element_blank(),
  #         legend.text = element_blank(),
  #         legend.title = element_blank()) +
  #   scale_fill_gradient(low = "white", high = "orange") +
  #   guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10,col="black",ticks = FALSE))
  # print(p)
  # dev.off()
  
  
}



fisher_plot_pdf <- function(DF,main=""){
  colfunc <- colorRampPalette(c("firebrick1", "white","dodgerblue"))
  p <- ggplot(DF, aes(x = COL, y = ROW,fill=lfc)) +      geom_tile(col="black") +
    theme_minimal()+
    coord_equal(ratio=.75) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top") +
    theme(axis.text.x=element_text(angle = 90,size=25, vjust=1, hjust=.5,face = "bold",
                                   margin=margin(1,0,0,0)),
          axis.text.y=element_text(size=20, margin=margin(0,-3,0,0),face = "bold"),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12, vjust=-2, hjust=0.1,face = "bold"),
          legend.background = element_rect(color="grey",size=0.8, linetype="dashed"),
          axis.title.y = element_blank(),
          legend.title.align=0.1,
          plot.margin = unit(c(0,1,0,1), "cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    
    geom_text(aes(label =sign),size=4) + 
    ggtitle(main) +
    scale_fill_gradientn(colours = colfunc(3) ,
                         guide = guide_colorbar(barwidth = 0.8,
                                                title = 'Log2 \nOddsRatio',
                                                label.theme = element_text(size=10, vjust=1, hjust=.5,face = "bold",angle = 0),
                                                barheight = 10,
                                                nbin = 10,
                                                draw.ulim = FALSE, 
                                                draw.llim = FALSE,
                                                ticks = FALSE),
                         breaks=c(-2,0.2,2),
                         labels=c("-1.5","0","1.5"),
                         limits=c(-2,2)) 
  return(p)
}


fisher_plot_tiff <- function(DF){
  colfunc <- colorRampPalette(c("firebrick1", "white","dodgerblue"))
  p <- ggplot(fish.df, aes(x = COL, y = ROW,fill=lfc)) +      geom_tile(col="black") +
    theme_minimal()+
    coord_equal(ratio=.75) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top") +
    theme(axis.text.x=element_text(angle = 90,size=25, vjust=1, hjust=.5,face = "bold",
                                   margin=margin(1,0,0,0)),
          axis.text.y=element_text(size=20, margin=margin(0,-3,0,0),face = "bold"),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12, vjust=-2, hjust=0.1,face = "bold"),
          legend.background = element_rect(color="grey",size=0.8, linetype="dashed"),
          axis.title.y = element_blank(),
          legend.title.align=0.1,
          plot.margin = unit(c(0,1,0,1), "cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    theme_tiff+
    scale_fill_gradientn(colours = colfunc(3) ,
                         guide = guide_colorbar(barwidth = 0.8,
                                                title = 'Log2 \nOddsRatio',
                                                label.theme = element_text(size=10, vjust=1, hjust=.5,face = "bold",angle = 0),
                                                barheight = 10,
                                                nbin = 10,
                                                draw.ulim = FALSE, 
                                                draw.llim = FALSE,
                                                ticks = FALSE),
                         breaks=c(-2,0.2,2),
                         labels=c("-1.5","0","1.5"),
                         limits=c(-2,2)) 
  return(p)
}



plotEnrichmentGSEA <- function (pathway, stats, gseaParam = 1) 
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) +
    #  geom_point(color = "green",size = 0.1) + 
    geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = min(bottoms),colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = 0, colour = "black") + 
    geom_line(color = "firebrick",size=0.6) + 
    # theme_bw() + 
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = -diff/8, xend = x, yend = diff/8),width=0.01, size = 0.05) +
    # theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + theme_bw()+
    labs(x = "rank", y = "enrichment score")
  g
}


#TOGPLOT PERCENTAGE NO SD
toGplotEnrichmentLoop <- function(o_name,dat,paper=F,wd=10){
  
  p <- ggplot(dat,aes(x = Input,y = value,fill=Var)) +
    geom_bar(stat = "identity",position = "dodge") +
    scale_fill_manual(values = c("#f0adad","#f08181","#42f47d","#17ef5f"))
  if(paper){
    p <- p +    theme(panel.background=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 1),
                      panel.margin = unit(c(0, 0, 0, 0), "cm"),
                      plot.margin = rep(unit(0,"null"),4),
                      axis.ticks=element_blank(),
                      axis.text.y=element_blank(),
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank(),
                      axis.ticks.length = unit(0,"null"),
                      axis.ticks.margin = unit(0,"null"),
                      plot.background = element_rect(fill = "transparent",colour="black",size=1),
                      # plot.margin = unit(c(-1, -1.2, -1.2, -1.5), "cm"),  # Edited code
                      legend.position = 'none') +
      labs(x=NULL, y=NULL)
  }else{
    p <- p +    theme(panel.background=element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 1),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      panel.margin = unit(c(0, 0, 0, 0), "cm"),
                      plot.margin = rep(unit(0,"null"),4),
                      axis.ticks=element_blank(),
                      axis.ticks.length = unit(0,"null"),
                      axis.ticks.margin = unit(0,"null"),
                      plot.background = element_rect(fill = "transparent",colour="black",size=1))
  }
  
  
  ggsave(filename=paste0(o_name,".pdf"),width = wd,plot=p)
  ggsave(filename=paste0(o_name,".jpg"), width = wd,plot=p)
  
}

toGplotCorrelation <- function(o_name,dat,paper=F,wd=10){
  p <- ggplot(data =  dat, aes(x = row, y = col)) +
    geom_tile(aes(fill = ratio), colour = "white") +
    geom_text(aes(fill = ratio, label = round(-log10(pv)))) +
    scale_fill_gradient(low = "green", high = "red") +
    coord_fixed(ratio=0.7)
  p <- p +    theme(panel.background=element_blank(),
                    axis.text.x = element_text(angle = 90, hjust = 1),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.margin = unit(c(0, 0, 0, 0), "cm"),
                    plot.margin = rep(unit(0,"null"),4),
                    axis.ticks=element_blank(),
                    axis.ticks.length = unit(0,"null"),
                    axis.ticks.margin = unit(0,"null"),
                    plot.background = element_rect(fill = "transparent",colour="black",size=1))
  ggsave(filename=paste0(o_name,".pdf"),width = 10,plot=p)
  ggsave(filename=paste0(o_name,".jpg"), width = 10,plot=p)
  
}

#                                         --- SCATTER PLOTS ---
# --------------------------------------------------------------------------------------------------------------#

toGplotScatter <- function(o_name,dat,paper=F){
  p <- ggplot(df,aes(x=X,y=Y))+ geom_point(aes(size = size)) +
    geom_text(aes(label=.id,col="red"), size=3)+facet_grid(~ split_f)+ theme_bw()
  ggsave(filename=paste0(o_name,".pdf"),width = 10,plot=p)
  ggsave(filename=paste0(o_name,".jpg"),width = 10,plot=p)
}

# PROFILES
# Input : Granges of feature of interest, GRanges of Domains, name of output, name of domains, outputdirectory
# Output : Graphic of Relative ditribution (1 Domain = 100 bins) of the feature inside the domain borders
feat_at_domain_profile <- function(feat.GR, tad.GR,name,tadName,out_dir){
  
  feat.GR <- whichDomain(feat.GR,tad.GR)
  dist.intv <- start(resize(feat.GR,1,fix="center"))-start(tad.GR[feat.GR$dom])
  class.intv <- ceiling((dist.intv/width(tad.GR[feat.GR$dom]))*100)
  
  feat.LB.DIST <-  data.frame(dist=class.intv[class.intv %in% 1:50]-1,type='left')
  feat.RB.DIST <- data.frame(dist= class.intv[class.intv %in% 51:100]-1,type='right')
  
  distFeat.Tot.DF <- rbind(feat.LB.DIST,feat.RB.DIST)
  # To use for fills, add
  pdf(paste0(out_dir,"/Distribution_of_",name,"_peaks_at_",tadName,"_domain_borders.pdf"))
  # pdf(paste0(out_dir,"/Distribution of toto.pdf"))
  p <- ggplot(distFeat.Tot.DF, aes(dist, fill = type)) + geom_histogram(aes(y=..density..),binwidth=1) + geom_density(alpha=.2) +ylim(c(0,0.2))
  p <- p + theme_classic()
  p <- p + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + scale_x_continuous(waiver(),labels=c("Left Border","","Middle","","Right Border"))
  print(p)
  dev.off()
}

# K27 profiles
# Input :vector of features that will be used as a get(filter[1]) etc.
# Input : matrecap, Granges of tss, output directory, directory of bigwig K27 file that will be plotted




seqplot_prof <- function(filt_l,mrc,tss.gr,bw.l,out_f,tmp){
  pdf(paste0(out_f,".pdf"),width=20)
  for(filt in filt_l){
    g_fbgn <- mrc[get(filt),"Fbgnid"]
    sze <- length(g_fbgn)
    print(sze)
    if(!(sze)) next()
    o.tmp <- toString(export.bed(tss.gr[names(tss.gr) %in% g_fbgn],paste0(tmp,"prof_",filt,"_#",sze,".bed")))
    plotTMP <- getPlotSetArray(bw.l,o.tmp,"dm3",bin = 20L,xmin = 2000,xmax=2000,rm0 = T)
    file.remove(paste0(o.tmp))
    print(plotAverage(plotTMP,xlab='Relative position [bp]', ylim=c(0,5), ylab='Signal',error.estimates = T))
  }
  dev.off()
}


theme_tiff <- theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.title.y=element_blank(),
                    legend.position = "none"
)


#                                         --- 3D PLOTS ---
# --------------------------------------------------------------------------------------------------------------#
surf.colors <- function(x, col = terrain.colors(20)) {
  
  # First we drop the 'borders' and average the facet corners
  # we need (nx - 1)(ny - 1) facet colours!
  x.avg <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] +
              x[-(nrow(x) -1), -1] + x[-(nrow(x) -1), -(ncol(x) - 1)]) / 4
  
  # Now we construct the actual colours matrix
  colors = col[cut(x.avg, breaks = length(col), include.lowest = T)]
  
  return(colors)
}


apa3dPlot <- function(feat_apa.mtx,o_plot){
  # z <- feat_apa.mtx
  z <- feat_apa.mtx
  jpeg(paste0(o_plot,"_3d_plot.jpeg"),width=5, height=5, units="in", res=500)
  persp(z,col = surf.colors(z,col = rev(heat.colors(80))),theta = 145,phi = 30, border = NA, shade = 2.5)
  dev.off()
  
  jpeg(paste0(o_plot,"_plot3D_plot.jpeg"),width=5, height=5, units="in", res=500)
  par(mar=c(0, 0, 0, 0),oma=c(0,0,0,0))
  myMat <- feat_apa.mtx
  hist3D(z = t(myMat), x=1:nrow(myMat),y=1:ncol(myMat),scale = T, expand = 1, bty = "g", phi = 35,theta=10,
         border = "black", ltheta = 5,
         space = 0.01, d = 1,axes=FALSE,plot=F,colkey = FALSE)
  
  # arrows3D(11,11,myMat[11,11]+max(myMat)/3,
  #          z1=myMat[11,11]+max(myMat)/8,
  #          col="black",
  #          border="white",
  #          lwd = 5,
  #          add=T,
  #          plot=T,bty = "10",type = "cone")
  dev.off()
  jpeg(paste0(o_plot,"_plot3D_legend_plot.jpeg"),width=5, height=5, units="in", res=500)
  myMat <- feat_apa.mtx
  hist3D(z = t(myMat), x=1:nrow(myMat),y=1:ncol(myMat),scale = T, expand = 1, bty = "g", phi = 35,theta=10,
         border = "black", ltheta = 5,
         space = 0.01, d = 1,axes=FALSE,plot=F,colkey = list(side = 2, length = 0.7))
  
  arrows3D(3,3,myMat[3,3]+max(myMat)/3,
           z1=myMat[3,3]+max(myMat)/8,
           col="black",
           border="white",
           lwd = 5,
           add=T,
           plot=T,bty = "10",type = "cone")
  dev.off()
}

simple3Dplot <-  function(feat_apa.mtx,o_plot){
  # z <- feat_apa.mtx
  z <- feat_apa.mtx
  jpeg(paste0(o_plot,"_3d_plot.jpeg"),width=5, height=5, units="in", res=500)
  persp(z,col = surf.colors(z,col = rev(heat.colors(80))),theta = 145,phi = 30, border = NA, shade = 2.5)
  dev.off()
  
  jpeg(paste0(o_plot,"_plot3D_plot.jpeg"),width=5, height=5, units="in", res=500)
  par(mar=c(0, 0, 0, 0),oma=c(0,0,0,0))
  myMat <- feat_apa.mtx
  hist3D(z = t(myMat), x=1:nrow(myMat),y=1:ncol(myMat),scale = T, expand = 1, bty = "g", phi = 35,theta=10,
         border = "black", ltheta = 5,
         space = 0.01, d = 1,axes=FALSE,plot=F,colkey = FALSE)
  
  # arrows3D(11,11,myMat[11,11]+max(myMat)/3,
  #          z1=myMat[11,11]+max(myMat)/8,
  #          col="black",
  #          border="white",
  #          lwd = 5,
  #          add=T,
  #          plot=T,bty = "10",type = "cone")
  dev.off()
  jpeg(paste0(o_plot,"_plot3D_legend_plot.jpeg"),width=5, height=5, units="in", res=500)
  myMat <- feat_apa.mtx
  hist3D(z = t(myMat), x=1:nrow(myMat),y=1:ncol(myMat),scale = T, expand = 1, bty = "g", phi = 35,theta=10,
         border = "black", ltheta = 5,
         space = 0.01, d = 1,axes=FALSE,plot=F,colkey = list(side = 2, length = 0.7))
  
  arrows3D(3,3,myMat[3,3]+max(myMat)/3,
           z1=myMat[3,3]+max(myMat)/8,
           col="black",
           border="white",
           lwd = 5,
           add=T,
           plot=T,bty = "10",type = "cone")
  dev.off()
}

diff3dPlot <- function(apa_with,apa_without,o_plot){
  norm <- sum(apa_with)/sum(apa_without)
  apa_without <- apa_without * norm
  z <- apa_with - apa_without
  jpeg(paste0(o_plot,"_3d_diff_plot.jpeg"),width=5, height=5, units="in", res=500)
  persp(z,col = surf.colors(z,col = rev(heat.colors(80))),theta = 145,phi = 30, border = NA, shade = 2.5)
  dev.off()
  jpeg(paste0(o_plot,"_plot3D_diff_plot.jpeg"),width=5, height=5, units="in", res=500)
  par(mar=c(0, 0, 0, 0),oma=c(0,0,0,0))
  myMat <- z[11:31,11:31]
  hist3D(z = t(myMat), x=1:nrow(myMat),y=1:ncol(myMat),scale = T, expand = 1, bty = "g", phi = 35,theta=10,
         border = "black", ltheta = 5,
         space = 0.01, d = 1,axes=FALSE,plot=F,colkey = FALSE)
  
  arrows3D(11,11,myMat[11,11]+max(myMat)/3,
           z1=myMat[11,11]+max(myMat)/8,
           col="black",
           border="white",
           lwd = 5,
           add=T,
           plot=T,bty = "10",type = "cone")
  dev.off()
  jpeg(paste0(o_plot,"_plot3D_diff_legend_plot.jpeg"),width=5, height=5, units="in", res=500)
  myMat <- z[11:31,11:31]
  hist3D(z = t(myMat), x=1:nrow(myMat),y=1:ncol(myMat),scale = T, expand = 1, bty = "g", phi = 35,theta=10,
         border = "black", ltheta = 5,
         space = 0.01, d = 1,axes=FALSE,plot=F,colkey = list(side = 2, length = 0.7))
  
  arrows3D(11,11,myMat[11,11]+max(myMat)/3,
           z1=myMat[11,11]+max(myMat)/8,
           col="black",
           border="white",
           lwd = 5,
           add=T,
           plot=T,bty = "10",type = "cone")
  dev.off()
}


#   -----------------------------------------------------------------------






# UNCLASSED ---------------------------------------------------------------
# CORREALTION MATRIX CONSTRUCTION
cor_matrix_fcomp <- function(f_df,uv_vec,out_d){
  cor.l <- lapply(f_df,function(x){cor <- cor.test(x,uv_vec)$estimate})
  cor_df <- ldply(cor.l,data.frame)
  colnames(cor_df) <- c(".id","correlation")
  cor_df$.id <- as.factor(cor_df$.id)
  
  n <- length(uv_vec)
  col.df <- NULL
  row.df <- NULL
  cval.df <- NULL
  pval.df <- NULL
  trash <- combn(1:nrow(cor_df),2,
                 function(x){
                   row.df <<- rbind(row.df, data.frame(row=cor_df[x[1],]$.id))
                   col.df <<- rbind(col.df,data.frame(col=cor_df[x[2],]$.id))
                   c1 <- cor_df[x[1],]$correlation
                   c2 <- cor_df[x[2],]$correlation
                   z <- (c1 - c2) / sqrt(2*(1 / (n - 3)))
                   pval.df <<- c(pval.df,(1-pnorm(abs(z)))*2 )
                   cval.df <<- c(cval.df,c1/c2)
                 },simplify = F)
  
  mydf <- cbind(row=row.df,col=col.df,ratio=cval.df,pv=pval.df)
  o_name <- paste0(out_d,"_matrix")
  toGplotCorrelation(o_name,mydf)
  
  recap_mat <- NULL
  for(i in unique(c(unique(as.character(mydf[,1])),unique(as.character(mydf[,2]))))){
    recap_mat[[i]] <- 0
  }
  for(x in 1:nrow(mydf)){
    add1 <- ifelse((mydf[x,3] > 1) & (mydf[x,4] < 0.05),1,ifelse((mydf[x,4] < 0.05),-1,0))
    recap_mat[[as.character(mydf[x,1])]] <- recap_mat[[as.character(mydf[x,1])]]+add1
    add2 <-  ifelse((mydf[x,3] > 1) & (mydf[x,4] < 0.05),-1,ifelse((mydf[x,4] < 0.05),1,0))
    recap_mat[[as.character(mydf[x,2])]] <- recap_mat[[as.character(mydf[x,2])]] + add2
  }
  recap.df <- data.frame(Row=names(recap_mat),value=recap_mat)
  recap.df$Row <- factor(recap.df$Row, levels = recap.df$Row[order(recap.df$value)])
  o_name=paste0(out_d,"_ranking")
  simplePlot(o_name,recap.df)
  
}


#  IMPORT BEDGRAPH --------------------------------------------------------

read.bedgraph <- function(file) {
  dat <- scan(file=file,
              what=list(character(),integer(),integer(),numeric()), sep="\t", skip=1)
  dat <- data.frame(chr=dat[[1]], start=dat[[2]], end=dat[[3]],
                    val=dat[[4]])
  return(dat)
}


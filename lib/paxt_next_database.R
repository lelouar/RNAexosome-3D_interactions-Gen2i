#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-


# TSS TXDB ANNOTATION
genes.hg19.gr <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genesB.hg19.gr <- subset(genes.hg19.gr,width > 1200)
genesBody.hg19.gr <- reduceGR(genesB.hg19.gr,upstream = 1000,downstream = 0)
genesBody.hg19.gr <- addSeqinfo(genesBody.hg19.gr,g_name = "hg19",celltype = "hela",chrTYPE = "NCBI")

tss.hg19.gr <- addSeqinfo(resize(genes.hg19.gr,1,"start"),g_name = "hg19",celltype = "hela",chrTYPE = "NCBI")
genes.hg19.sglTSS.pm2kb.gr <- readRDS("resources/ext_data/Projects/PROMPTS/lncRNA_Z1/Parametrization/tss_sense_aware_generator/gene.hg19.sglTSS.pm2kb.gr")
tss.hg19.sglTSS.pm2kb.gr <- resize(genes.hg19.sglTSS.pm2kb.gr,1,"start")
TSSnoDiv.hg19.sglTSS.2kb0.gr <- g2i:::setSeqinfo(readRDS("resources/ext_data/Projects/PROMPTS/lncRNA_Z1/Parametrization/tss_sense_aware_generator/gene.hg19.TSSnodivgenes.2kb0.gr"),"hg19","NCBI")
# prom2Kb.hg19.gr <- extendGR(TSSnoDiv.hg19.sglTSS.2kb0.gr,2000,0)
# prom2Kb.hg19.gr <- extendGR(tss.hg19.sglTSS.pm2kb.gr,2000,0)
tss_noSTR.hg19.gr <- tss.hg19.gr
strand(tss_noSTR.hg19.gr) <- "*"

# DESEQ2 ANALYSIS
# /!\ UR & DR are correct here (lfc > 0 = UpReg in KD condition)
# BUt i Inverted it in analysis
tx_dir <- "resources/ext_data/Projects/PAXT_NEXT/PROPER/LOAD_DATA/"
(tx.gr <- resize(readRDS(paste0(tx_dir,"matrecap_genes.rds")),1,"start"))


# RNASeQ
MTR4_ctl1_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/MTR4_RUN1_SE/BW/BW_REVERSED/Scr2_library1_rep2_201408_RNAseqSE.Aligned.sortedByCoord.out.dedup.fwd_RPGC.bw"
MTR4_ctl1_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/MTR4_RUN1_SE/BW/BW_REVERSED/Scr2_library1_rep2_201408_RNAseqSE.Aligned.sortedByCoord.out.dedup.rev_RPGC.bw"
MTR4_ctl2_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/MTR4_RUN1_SE/BW/BW_REVERSED/Scr1_library2_rep1_201408_RNAseqSE.Aligned.sortedByCoord.out.dedup.fwd_RPGC.bw"
MTR4_ctl2_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/MTR4_RUN1_SE/BW/BW_REVERSED/Scr1_library2_rep1_201408_RNAseqSE.Aligned.sortedByCoord.out.dedup.rev_RPGC.bw"

MTR4_kd1_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/MTR4_RUN1_SE/BW/BW_REVERSED/Mtr4-2_library14_rep2_201408_RNAseqSE.Aligned.sortedByCoord.out.dedup.fwd_RPGC.bw"
MTR4_kd1_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/MTR4_RUN1_SE/BW/BW_REVERSED/Mtr4-2_library14_rep2_201408_RNAseqSE.Aligned.sortedByCoord.out.dedup.rev_RPGC.bw"
MTR4_kd2_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/MTR4_RUN1_SE/BW/BW_REVERSED/Mtr4-1_library13_rep1_201408_RNAseqSE.Aligned.sortedByCoord.out.dedup.fwd_RPGC.bw"
MTR4_kd2_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/MTR4_RUN1_SE/BW/BW_REVERSED/Mtr4-1_library13_rep1_201408_RNAseqSE.Aligned.sortedByCoord.out.dedup.rev_RPGC.bw"

# Z1 DEPLETED + CONTROL
Z1_ctl1_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/WT_JENSEN/BW/SRR3757038_39_merge.Aligned.sortedByCoord.out.dedup.fwd_RPGC.bw"
Z1_ctl1_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/WT_JENSEN/BW/SRR3757038_39_merge.Aligned.sortedByCoord.out.dedup.rev_RPGC.bw"
Z1_ctl2_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/WT_JENSEN/BW/SRR3757040_41_merge.Aligned.sortedByCoord.out.dedup.fwd_RPGC.bw"
Z1_ctl2_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/WT_JENSEN/BW/SRR3757040_41_merge.Aligned.sortedByCoord.out.dedup.rev_RPGC.bw"

Z1_kd1_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/Z1_JENSEN/BW/SRR3757088_89_merge.Aligned.sortedByCoord.out.fwd_RPGC.bw"
Z1_kd1_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/Z1_JENSEN/BW/SRR3757088_89_merge.Aligned.sortedByCoord.out.rev_RPGC.bw"
Z1_kd2_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/Z1_JENSEN/BW/SRR3757090_91_merge.Aligned.sortedByCoord.out.fwd_RPGC.bw"
Z1_kd2_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/Z1_JENSEN/BW/SRR3757090_91_merge.Aligned.sortedByCoord.out.rev_RPGC.bw"

# Z8 DEPLETED
Z8_kd1_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/Z8/BW/SRR3757077_1.fq.gz_Q10.Aligned.sortedByCoord.out.dedup.fwd_RPGC_scale_1.0.bw"
Z8_kd1_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/Z8/BW/SRR3757077_1.fq.gz_Q10.Aligned.sortedByCoord.out.dedup.rev_RPGC_scale_1.0.bw"
Z8_kd2_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/Z8/BW/SRR375707879_1.fq.gz_Q10.Aligned.sortedByCoord.out.dedup.fwd_RPGC_scale_1.0.bw"
Z8_kd2_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/Z8/BW/SRR375707879_1.fq.gz_Q10.Aligned.sortedByCoord.out.dedup.rev_RPGC_scale_1.0.bw"

ZWT_1_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/WT_JENSEN/BW/SRR3757038_39_merge.Aligned.sortedByCoord.out.dedup.fwd_RPGC.bw"
ZWT_1_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/WT_JENSEN/BW/SRR3757038_39_merge.Aligned.sortedByCoord.out.dedup.rev_RPGC.bw"
ZWT_2_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/WT_JENSEN/BW/SRR3757040_41_merge.Aligned.sortedByCoord.out.dedup.fwd_RPGC.bw"
ZWT_2_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/WT_JENSEN/BW/SRR3757040_41_merge.Aligned.sortedByCoord.out.dedup.rev_RPGC.bw"


MWT_1_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/MTR4_RUN1_SE/BW/BW_REVERSED/Scr2_library1_rep2_201408_RNAseqSE.Aligned.sortedByCoord.out.dedup.fwd_RPGC.bw"
MWT_1_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/MTR4_RUN1_SE/BW/BW_REVERSED/Scr2_library1_rep2_201408_RNAseqSE.Aligned.sortedByCoord.out.dedup.rev_RPGC.bw"
MWT_2_fwd.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/MTR4_RUN1_SE/BW/BW_REVERSED/Scr1_library2_rep1_201408_RNAseqSE.Aligned.sortedByCoord.out.dedup.fwd_RPGC.bw"
MWT_2_rev.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/rna-seq/MTR4_RUN1_SE/BW/BW_REVERSED/Scr1_library2_rep1_201408_RNAseqSE.Aligned.sortedByCoord.out.dedup.rev_RPGC.bw"


# IPSC RNAseq
IPSC_ctl1_fwd.bwp <- "resources/ext_data//Database/Human/iPSC/hg19/rnaseq/CORNELIADELANGE/WT/WT_REP1/BW/WT_REP1_R1.fastq.gz_Q10.Aligned.sortedByCoord.out.dedup.fwd_RPGC_scale_1.0_no_input.bw"
IPSC_ctl1_rev.bwp <- "resources/ext_data//Database/Human/iPSC/hg19/rnaseq/CORNELIADELANGE/WT/WT_REP1/BW/WT_REP1_R1.fastq.gz_Q10.Aligned.sortedByCoord.out.dedup.rev_RPGC_scale_1.0_no_input.bw"
IPSC_ctl2_fwd.bwp <- "resources/ext_data//Database/Human/iPSC/hg19/rnaseq/CORNELIADELANGE/WT/WT_REP2/BW/WT_REP2_R1.fastq.gz_Q10.Aligned.sortedByCoord.out.dedup.fwd_RPGC_scale_1.0_no_input.bw"
IPSC_ctl2_rev.bwp <- "resources/ext_data//Database/Human/iPSC/hg19/rnaseq/CORNELIADELANGE/WT/WT_REP2/BW/WT_REP2_R1.fastq.gz_Q10.Aligned.sortedByCoord.out.dedup.rev_RPGC_scale_1.0_no_input.bw"

IPSC_kd1_fwd.bwp <- "resources/ext_data//Database/Human/iPSC/hg19/rnaseq/CORNELIADELANGE/CDLS/CDLS_REP1/BW/CDLS_REP1_R1.fastq.gz_Q10.Aligned.sortedByCoord.out.dedup.fwd_RPGC_scale_1.0_no_input.bw"
IPSC_kd1_rev.bwp <- "resources/ext_data//Database/Human/iPSC/hg19/rnaseq/CORNELIADELANGE/CDLS/CDLS_REP1/BW/CDLS_REP1_R1.fastq.gz_Q10.Aligned.sortedByCoord.out.dedup.rev_RPGC_scale_1.0_no_input.bw"
IPSC_kd2_fwd.bwp <- "resources/ext_data//Database/Human/iPSC/hg19/rnaseq/CORNELIADELANGE/CDLS/CDLS_REP3/BW/CDLS_REP3_R1.fastq.gz_Q10.Aligned.sortedByCoord.out.dedup.fwd_RPGC_scale_1.0_no_input.bw"
IPSC_kd2_rev.bwp <- "resources/ext_data//Database/Human/iPSC/hg19/rnaseq/CORNELIADELANGE/CDLS/CDLS_REP3/BW/CDLS_REP3_R1.fastq.gz_Q10.Aligned.sortedByCoord.out.dedup.rev_RPGC_scale_1.0_no_input.bw"



# CHIPSEQ
# Z1
z1.hg19.gr <- addSeqinfo(import("resources/ext_data//Database/Human/helaS3/hg19/chip-seq/Z1/Kiernan_Raoul_July19/L_Z1/IP_L_Z1_peaks.narrowPeak.cutted.bed"),"hg19","hela") 
z1.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/Z1/Kiernan_Raoul_July19/L_Z1/IP_L_Z1.input_normalised.bw"
z1NI.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/Z1/Kiernan_Raoul_July19/L_Z1/IP_L_Z1.bw"
z1.shMTR4.hg19.gr <- addSeqinfo(import("resources/ext_data//Database/Human/helaS3/hg19/chip-seq/Z1/Kiernan_Raoul_July19/IP_shMTR4_Z1_normalized/peak_calling/hg19/narrow_nomodel/Z_shMTR4_Z1_summits.bed"),"hg19","hela")
z1.shMTR4.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/Z1_Z8_MTR4_Raoul_pipeline/Z_Rad21/Analyse/Z_shMTR4_Z1.input_normalised.bw"


# MTR4
# mtr4.hg19.gr <- addSeqinfo(import("resources/ext_data/Database/Human/helaS3/hg19/chip-seq/mtr4/Kiernan_IGH/BED/Kiernan_F44-MTR4_S1_hg19_peaks.narrowPeak.cutted.bed"),"hg19","hela")
mtr4.hg19.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/chip-seq/mtr4/Kiernan_IGH/BW/Kiernan_F44-MTR4_S1_hg19.sorted_RPGC.bw"
mtr4N.hg19.gr <- addSeqinfo(import("resources/ext_data//Database/Human/helaS3/hg19/chip-seq/mtr4/Kiernan_RaoulPIPELINE_July19/F44_MTR4/IP_F44_peaks.narrowPeak.cutted.bed"),"hg19","hela")
mtr4N.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/mtr4/Kiernan_RaoulPIPELINE_July19/F44_MTR4/IP_F44.input_normalised.bw"
mtr4NNI.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/mtr4/Kiernan_RaoulPIPELINE_July19/F44_MTR4/IP_F44.bw"


# Z8
# z8_RK.hg38.gr <- import("resources/ext_data/Database/Human/helaS3//hg38/chipseq/Z8/Rosemary_Pipeline/IP_F5_Z8/peak_calling/hg38/narrow_nomodel/Kiernan_F5-Z8_S6_L004_R1_001_summits.bed")
# z8_RK.hg38.gr <- g2i:::setSeqinfo(z8_RK.hg38.gr,"hg38","NCBI") # Work better with this replicate of Z8, the other (E2) is nonsense, peak are not in TSS nor in enhancer nor CTCF/RAD21 
# z8_RK.hg19.gr <- liftover(z8_RK.hg38.gr,hg38TOhg19,sqLvlstyle = "NCBI")
# z8.hg19.gr <- trim(addSeqinfo(z8_RK.hg19.gr,"hg19","hela")) # Work better with this replicate of Z8, the other (E2) is nonsense, peak are not in TSS nor in enhancer nor CTCF/RAD21 
z8.hg19.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/chip-seq/Z8_PAPG/BOWTIE/BW/Kiernan_12_F5-Z8-s2_submitted_S6_L004_R1_001_Q10_sorted_dedup_RPGC_ratio_inputReadCountNorm.bw"
# Work better with this replicate of Z8, the other (E2) is nonsense, peak are not in TSS nor in enhancer nor CTCF/RAD21 
z8N.hg19.gr <- addSeqinfo(import("resources/ext_data//Database/Human/helaS3/hg19/chip-seq/Z8/Kiernan_Raoul_July19/F5_Z8/IP_Z8_F5_peaks.narrowPeak.cutted.bed"),"hg19","hela")
z8N.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/Z8/Kiernan_Raoul_July19/F5_Z8/IP_Z8_F5.input_normalised.bw"
z8NNI.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/Z8/Kiernan_Raoul_July19/F5_Z8/IP_Z8_F5.bw"



# PROMPTS
ncRNA_Z1kd.hg19.gr <- readRDS(paste0("resources/ext_data/Projects/PROMPTS/lncRNA_Z1/find_prompts/normR_dedup_Jensen/prompts_all12_filtered_bs500_fdr05.rds"))
ncRNA_MTR4kd.hg19.gr <- readRDS(paste0("resources/ext_data/Projects/PROMPTS/lncRNA_hg19/find_prompts/normR/prompts_all12_filtered_bs500_fdr05.rds"))
ncRNA_Z8kd.hg19.gr <- readRDS(paste0("resources/ext_data/Projects/PROMPTS/INTEGRATION_PAPOLG_Z8/PROMPTS_DETECTION/prompts_all12_filtered_bs500_fdr05.rds"))



# promptsZ1.hg19.gr <- myOverlaps("ncRNAZ1","PROM",ncRNA_Z1kd.hg19.gr,prom2Kb.hg19.gr)$gr1xgr2
# promptsMTR4.hg19.gr <- myOverlaps("ncRNAMTR4","PROM",ncRNA_MTR4kd.hg19.gr,prom2Kb.hg19.gr)$gr1xgr2
# promptsZ8.hg19.gr <- myOverlaps("ncRNAZ8","PROM",ncRNA_Z8kd.hg19.gr,prom2Kb.hg19.gr)$gr1xgr2

# myLIST <- c("resources/ext_data//Database/Human/helaS3/hg19/chip-seq/mtr4/Kiernan_RaoulPIPELINE_July19/F44_MTR4/IP_F44.bw",
#             "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/Z1/Kiernan_Raoul_July19/L_Z1/IP_L_Z1.input_normalised.bw",
#             "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/Z1/Kiernan_Raoul_July19/L_Z1/IP_L_Z1.bw",
#             "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/Z8/Kiernan_Raoul_July19/F5_Z8/IP_Z8_F5.input_normalised.bw",
#             "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/Z8/Kiernan_Raoul_July19/F5_Z8/IP_Z8_F5.bw",
#             "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Raoul_July19/IP_shC_Rad21_normalized/bigwig/hg19/Z_shC_Rad21.bw",
#             "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Raoul_July19/IP_shC_Rad21_normalized/bigwig/hg19/Z_shC_Rad21.input_normalised.bw",
#             )
#
# for(item in myLIST){
#   tst <- import(item)
#   NM <- gsub("(.*)")
#   seqlevelsStyle(tst) <- "UCSC"
#   export.bw(tst,"resources/ext_data//Database/Human/helaS3/hg19/chip-seq/mtr4/Kiernan_RaoulPIPELINE_July19/F44_MTR4/IP_F44.input_normalised_chr.bw")
#
# }


# RAD21
rad21.shC.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Raoul_July19/IP_shC_Rad21_normalized/bigwig/hg19/Z_shC_Rad21.input_normalised.bw"
rad21.shC_NOINP.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Raoul_July19/IP_shC_Rad21_normalized/bigwig/hg19/Z_shC_Rad21.bw"
rad21.shC.hg19.gr <- addSeqinfo(import("resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Raoul_July19/IP_shC_Rad21_normalized/peak_calling/hg19/narrow_nomodel/Z_shC_Rad21_summits.bed"),"hg19","hela")
rad21.shMTR4.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Raoul_July19/IP_shMTR4_RAD21/bigwig/Z_shMTR4_Rad21.input_normalised.bw"
rad21.shMTR4.hg19.gr <- addSeqinfo(import("resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Raoul_July19/IP_shMTR4_RAD21/peakcalling/Z_shMTR4_Rad21_peaks.narrowPeak.cutted.bed"),"hg19","hela")
rad21.shZ1.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Xavier_2020/Zj4_shZ1_Rad21.input_normalised.bw"
rad21.shZ1noINP.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Xavier_2020/Zj4_shZ1_Rad21.bw"
# rad21.shZ1.hg19.gr <- addSeqinfo(import("resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Xavier_2020/BED/Zj4_shZ1_Rad21_peaks.merged_nomodel.bed"),"hg19","hela")


rad212.shC.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Xavier_2020/Zj10_shC_Rad21.input_normalised.bw"
rad212.shC_NOINP.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Xavier_2020/Zj10_shC_Rad21.bw"
# rad212.shC.hg19.gr <- addSeqinfo(import("resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Xavier_2020/BED/Zj10_shC_Rad21_peaks.merged_nomodel.bed"),"hg19","hela")
rad212.shMTR4.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Xavier_2020/Zj16_shMTR4_Rad21.input_normalised.bw"
# rad212.shMTR4.hg19.gr <- addSeqinfo(import("resources/ext_data//Database/Human/helaS3/hg19/chip-seq/RAD21/Kiernan_Xavier_2020/BED/Zj16_shMTR4_Rad21_peaks.merged_nomodel.bed"),"hg19","hela")




# We have checked that the overlap with existing CTCF in litterature is reproducible so we will take the intersection of both
# Rad21 - Modencode
# rad21.hg19.gr <- addSeqinfo(GenomicRanges::reduce(import("resources/ext_data/Database/Human/helaS3//hg19/chip-seq/apbs/rad21/GSM935571/GSM935571_hg19_wgEncodeSydhTfbsHelas3Rad21IggrabPk.narrowPeak.cutted.bed")),"hg19","hela")

# CTCF - Kiernan 78148
# ctcfK.hg19.bwp <- "resources/ext_data//Database/Human/helaS3/hg19/chip-seq/CTCF/Kiernan_Raoul_July19/L_CTCF/IP_CTCF_L.input_normalised.bw"
ctcfK.hg19.gr <- addSeqinfo(GenomicRanges::reduce(import("resources/ext_data//Database/Human/helaS3/hg19/chip-seq/CTCF/Kiernan_Raoul_July19/L_CTCF/IP_CTCF_L_peaks.narrowPeak.cutted.bed")),"hg19","hela")
ctcfG.hg19.gr <- addSeqinfo(GenomicRanges::reduce(import("resources/ext_data/Database/Human/helaS3//hg19/chip-seq/apbs/ctcf/GSM749729/Hela-DS11552.peaks.fdr0.01.hg19_cutted.bed")),"hg19","hela")
ctcf.hg19.gr <- uniqueOvlp(ctcfK.hg19.gr,ctcfG.hg19.gr,1)


#HISTONE - GEO
h3k4me1.hg19.gr <- addSeqinfo(import("resources/ext_data/Database/Human/helaS3//hg19/chip-seq/histone/H3K4me1/GSM798322/GSM798322_hg19_wgEncodeBroadHistoneHelas3H3k04me1StdPk.broadPeak.bed"),"hg19","hela")
k4me1.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/chip-seq/histone/H3K4me1/GSM798322/BW/GSM798322_hg19_wgEncodeBroadHistoneHelas3H3k04me1StdSig.bigWig"
k4me3.hg19.gr <- addSeqinfo(import("resources/ext_data/Database/Human/helaS3//hg19/chip-seq/histone/H3K4me3/GSM733682/GSM733682_hg19_wgEncodeBroadHistoneHelas3H3k4me3StdPk.broadPeak.bed"),"hg19","hela")
k27ac.hg19.gr <- addSeqinfo(import("resources/ext_data/Database/Human/helaS3//hg19/chip-seq/histone/H3K27ac/GSM733684/GSM733684_hg19_wgEncodeBroadHistoneHelas3H3k27acStdPk.broadPeak.bed"),"hg19","hela")
k27me3.hg19.gr <- addSeqinfo(import("resources/ext_data/Database/Human/helaS3//hg19/chip-seq/histone/H3K27me3/wgEncodeBroadHistoneHelas3H3k27me3StdPk.bed"),"hg19","hela")
k9me3.hg19.gr <- addSeqinfo(import("resources/ext_data/Database/Human/helaS3//hg19/chip-seq/histone/H3K9me3/wgEncodeBroadHistoneHelas3H3k09me3Pk.bed"),"hg19","hela")
k9ac.hg19.gr <- addSeqinfo(import("resources/ext_data/Database/Human/helaS3//hg19/chip-seq/histone/H3K9ac/wgEncodeBroadHistoneHelas3H3k9acStdPk.bed"),"hg19","hela")
k36me3.hg19.gr <- addSeqinfo(import("resources/ext_data/Database/Human/helaS3//hg19/chip-seq/histone/H3K36me3/wgEncodeBroadHistoneHelas3H3k36me3StdPk.bed"),"hg19","hela")


# ENNHANCER
enh.hg19.gr <-  addSeqinfo(import("resources/ext_data/Database/Human/helaS3/hg19/starr-seq/GSE100432/GSE100432_peaks_inhibitor_correctedEnrichment4_supp.table3.bed"),"hg19","hela")
enhCAGE.hg19.gr <-  addSeqinfo(readRDS("resources/ext_data/Database/Human/helaS3/hg19/cage/HACER_DB/CAGE-HelaS3_ENCODE.rds"),"hg19","hela")


# BW CTL IGG
ctrl.bwp <- "resources/ext_data/Database/Human/helaS3/hg19/chip-seq/Igg/GSM935560_hg19_wgEncodeSydhTfbsHelas3Rfx5200401194IggrabSig.bigWig"

# BIGWIG LIST
Z1_RNA.bwl <- list(REP1KD=list(FWD=Z1_kd1_fwd.bwp,REV=Z1_kd1_rev.bwp),
                   REP1WT=list(FWD=ZWT_1_fwd.bwp,REV=ZWT_1_rev.bwp),
                   REP2KD=list(FWD=Z1_kd2_fwd.bwp,REV=Z1_kd2_rev.bwp),
                   REP2WT=list(FWD=ZWT_2_fwd.bwp,REV=ZWT_2_rev.bwp))
Z8_RNA.bwl <- list(REP1KD=list(FWD=Z8_kd1_fwd.bwp,REV=Z8_kd1_rev.bwp),
                   REP1WT=list(FWD=ZWT_1_fwd.bwp,REV=ZWT_1_rev.bwp),
                   REP2KD=list(FWD=Z8_kd2_fwd.bwp,REV=Z8_kd2_rev.bwp),
                   REP2WT=list(FWD=ZWT_2_fwd.bwp,REV=ZWT_2_rev.bwp))
MTR4_RNA.bwl <- list(REP1KD=list(FWD=MTR4_kd1_fwd.bwp,REV=MTR4_kd1_rev.bwp),
                     REP1WT=list(FWD=MWT_1_fwd.bwp,REV=MWT_1_rev.bwp),
                     REP2KD=list(FWD=MTR4_kd2_fwd.bwp,REV=MTR4_kd2_rev.bwp),
                     REP2WT=list(FWD=MWT_2_fwd.bwp,REV=MWT_2_rev.bwp))
CDLS_RNA.bwl <- list(REP1KD=list(FWD=IPSC_kd1_fwd.bwp,REV=IPSC_kd1_rev.bwp),
                     REP1WT=list(FWD=IPSC_ctl1_fwd.bwp,REV=IPSC_ctl1_rev.bwp),
                     REP2KD=list(FWD=IPSC_kd2_fwd.bwp,REV=IPSC_kd2_rev.bwp),
                     REP2WT=list(FWD=IPSC_ctl2_fwd.bwp,REV=IPSC_ctl2_rev.bwp))

# PARAM
BS <- 500
xtend <- 500



# SUB functions
fisher.grad.grad <- function(df,grad_row,grad_col,verbose=T) {
  r_df <- NULL
  for(i in na.omit(unique(df[[grad_row]]))){
    for(j in na.omit(unique(df[[grad_col]]))){
      cr <- as.numeric(nrow(as.data.frame(df) %>% filter(!!sym(grad_row)==i & !!sym(grad_col)==j)))
      cnr <- as.numeric(nrow(as.data.frame(df) %>% filter(!!sym(grad_row)!=i & !!sym(grad_col)==j)))
      ncr <- as.numeric(nrow(as.data.frame(df) %>% filter(!!sym(grad_row)==i & !!sym(grad_col)!=j)))
      ncnr <- as.numeric(nrow(as.data.frame(df) %>% filter(!!sym(grad_row)!=i & !!sym(grad_col)!=j)))
      f_mat <- matrix(c(L_C=cr,L_nC=cnr,C_nL=ncr,nC_nL=ncnr),nrow=2,ncol=2)
      if(verbose) print(f_mat)
      ft <- fisher.test(f_mat)
      mypv <- formatC(as.numeric(ft$p.value, format = "e", digits = 2))
      lfc <- round(log2(ft$estimate),2)
      # deal with inf cases
      lfc <- ifelse(is.finite(lfc),lfc,sign(lfc)*10)
      sign <- ifelse(ft$p.value < 0.05,ifelse(ft$p.value < 0.01,ifelse(ft$p.value < 0.001,"***","**"),"*"),"NS")
      r_df <- rbind.data.frame(r_df,
                               cbind.data.frame(ROW=as.integer(i),COL=as.integer(j),
                                                pv=mypv,lfc=lfc,sign=sign))
    }
  }
  r_df$ROW <- factor(r_df$ROW,levels=sort(unique(as.numeric(as.character(r_df$ROW))),decreasing = T))
  r_df$COL <- factor(r_df$COL,levels=sort(unique(as.numeric(as.character(r_df$COL)))))
  r_df
}

fisher.bin.grad <- function(df,bin_row,grad_col,verbose=T) {
  r_df <- NULL
  for(i in na.omit(unique(df[[bin_row]]))){
    for(j in na.omit(unique(df[[grad_col]]))){
      cr <- as.numeric(nrow(as.data.frame(df) %>% filter(!!sym(bin_row)==i & !!sym(grad_col)==j)))
      cnr <- as.numeric(nrow(as.data.frame(df) %>% filter(!!sym(bin_row)!=i & !!sym(grad_col)==j)))
      ncr <- as.numeric(nrow(as.data.frame(df) %>% filter(!!sym(bin_row)==i & !!sym(grad_col)!=j)))
      ncnr <- as.numeric(nrow(as.data.frame(df) %>% filter(!!sym(bin_row)!=i & !!sym(grad_col)!=j)))
      f_mat <- matrix(c(L_C=cr,L_nC=cnr,C_nL=ncr,nC_nL=ncnr),nrow=2,ncol=2)
      if(verbose) print(f_mat)
      ft <- fisher.test(f_mat)
      mypv <- formatC(as.numeric(ft$p.value, format = "e", digits = 2))
      lfc <- round(log2(ft$estimate),2)
      # deal with inf cases
      lfc <- ifelse(is.finite(lfc),lfc,sign(lfc)*10)
      sign <- ifelse(ft$p.value < 0.05,ifelse(ft$p.value < 0.01,ifelse(ft$p.value < 0.001,"***","**"),"*"),"NS")
      r_df <- rbind.data.frame(r_df,
                               cbind.data.frame(ROW=i,COL=as.integer(j),
                                                pv=mypv,lfc=lfc,sign=sign))
    }
  }
  r_df$COL <- factor(r_df$COL,levels=sort(unique(as.numeric(as.character(r_df$COL)))))
  r_df
}



fisher_plot <- function(DF,main="",range=c(-2,2), invert = F, title.size=22, xl="Undefined Quantiles", yl="Undefined Quantiles"){
  if (!is.null(range)){
    min_lim <- min(range) # New, before it was set automtcally to -2
    max_lim <- max(range) # New, before it was set automtcally to 2
    DF <- DF %>% mutate_cond(lfc < 0, lfc=rangeMinMax(x = lfc,min=min_lim,max=0)) %>% 
      mutate_cond(lfc >= 0, lfc=rangeMinMax(lfc,0,max_lim))
  } else {
    min_lim <- min(DF$lfc)
    max_lim <- max(DF$lfc)
  }
  if (invert){
    colnames(DF) <- c("COL","ROW","pv","lfc","sign")
  }
  # colfunc <- colorRampPalette(c("mediumblue", "white","firebrick"))
  colfunc <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))
  p <- ggplot(DF, aes(x = COL, y = ROW,fill=lfc)) +      geom_tile(col=NA) +
    theme_minimal()+
    coord_equal(ratio=1) +
    xlab(xl)+ylab(yl)+
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top",expand=c(-1,0)) +
    geom_text(aes(label = sign),size=7) + 
    ggtitle(main) +
    scale_fill_gradientn(colours = colfunc(10) ,
                         guide = guide_colorbar(barwidth = 0.8,
                                                title = 'Normalized Log2 \nOddsRatio',
                                                label.theme = element_text(size=10, vjust=1, hjust=.5,face = "bold",angle = 0),
                                                barheight = 10,
                                                nbin = 20,
                                                draw.ulim = FALSE, 
                                                draw.llim = FALSE,
                                                ticks = FALSE),
                         # breaks=c(min(range),0,max(range)),
                         breaks=c(min_lim,0,max_lim),
                         # labels=c("-2","0","2"),
                         labels=c(as.character(min_lim),"0",as.character(max_lim)),
                         # limits=c(-2,2)) 
                         limits=c(min_lim,max_lim)) 
  
  return(p)
}

rangeMinMax <- function(x,min=-1,max=1)
{
  return(min+((x- min(x,na.rm=T))*(max-min)) /(max(x,na.rm=T)-min(x,na.rm=T)))
}# TODO example : rangeMinMax(-5:5,-1,1)
#-----#

rangeMinMaxAsym <- function(x,nmin=-1,nmax=0,pmin=0,pmax=1)
{
  neg <- x[x<=0]
  pos <- x[x>=0]
  if(length(neg)){
    if(length(unique(neg))-1)
      ret_neg <- nmin+((neg- min(neg,na.rm=T))*(nmax-nmin)) /(max(neg,na.rm=T)-min(neg,na.rm=T))
    else
      ret_neg <- rep(nmin,length(neg))
  }else{
    ret_neg <- NULL
  }
  if(length(pos)){
    if(length(unique(pos))-1)
      ret_pos <- pmin+((pos- min(pos,na.rm=T))*(pmax-pmin)) /(max(pos,na.rm=T)-min(pos,na.rm=T))
    else
      ret_pos <- rep(pmax,length(pos))
  }else{
    ret_pos <- NULL
  }
  x[x<=0] <- ret_neg
  x[x>=0] <- ret_pos
  x
}# TODO example : rangeMinMax(-5:5,-1,1)


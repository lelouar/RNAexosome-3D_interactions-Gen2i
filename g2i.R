# STAR_METHODS/g2i.R
# Minimal standalone replacements for g2i / g2ismk package functions.
# Sourced by g2i_apa_lib.R, fig4.AB.R, fig5.AB.R.

# hg19 chromosome lengths (non-circular, NCBI seqlevel style)
.hg19_seqlengths <- c(
  "1"  = 249250621L, "2"  = 243199373L, "3"  = 198022430L, "4"  = 191154276L,
  "5"  = 180915260L, "6"  = 171115067L, "7"  = 159138663L, "8"  = 146364022L,
  "9"  = 141213431L, "10" = 135534747L, "11" = 135006516L, "12" = 133851895L,
  "13" = 115169878L, "14" = 107349540L, "15" = 102531392L, "16" =  90354753L,
  "17" =  81195210L, "18" =  78077248L, "19" =  59128983L, "20" =  63025520L,
  "21" =  48129895L, "22" =  51304566L, "X"  = 155270560L, "Y"  =  59373566L
)

# Build hg19 10kb tile GRanges with chr_bin annotation (required by APA engine).
# chr_bin = per-chromosome 1-based bin index used as row/col index into the H5 matrix.
makeTile10kb <- function() {
  si <- GenomeInfoDb::Seqinfo(
    seqnames   = names(.hg19_seqlengths),
    seqlengths = .hg19_seqlengths,
    genome     = "hg19"
  )
  gr <- GenomicRanges::tileGenome(
    GenomeInfoDb::seqlengths(si),
    tilewidth              = 10000L,
    cut.last.tile.in.chrom = TRUE
  )
  gr <- GenomeInfoDb::keepSeqlevels(gr, names(.hg19_seqlengths), pruning.mode = "coarse")
  gr$chr_bin <- unlist(lapply(split(gr, GenomeInfoDb::seqnames(gr)), seq_along))
  gr$name    <- paste(as.character(GenomeInfoDb::seqnames(gr)),
                      GenomicRanges::start(gr),
                      GenomicRanges::end(gr), sep = "_")
  names(gr) <- gr$name
  gr
}

# Load a BED file (3-6 columns) as GRanges with NCBI seqlevel style.
# Handles both UCSC ("chr1") and NCBI ("1") input — seqlevelsStyle converts automatically.
loadBED <- function(path) {
  dt <- data.table::fread(path, header = FALSE)
  gr <- GenomicRanges::GRanges(
    seqnames = as.character(dt[[1]]),
    ranges   = IRanges::IRanges(as.integer(dt[[2]]) + 1L, as.integer(dt[[3]]))
  )
  suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI")
  gr
}

# Load a named TAD / domain constraint GRanges.
# Used by pairwiseInteractionDT() when cons_nm is neither "intra" nor "inter".
loadConstraint <- function(cons_nm) {
  paths <- list(
    tad_rao_hela = "resources/ext_data/hg19/NGS/HiC/rk_tad_data/GSE63525_HeLa_Arrowhead_domainlist.bed"
  )
  if (!cons_nm %in% names(paths))
    stop(paste0("Unknown constraint: '", cons_nm,
                "'. Known: ", paste(names(paths), collapse = ", ")), call. = FALSE)
  loadBED(paths[[cons_nm]])
}

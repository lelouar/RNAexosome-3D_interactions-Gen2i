# STAR METHODS

## METHOD DETAILS

### ChIP-seq Peak Analysis and Genomic Distribution
ChIP-seq peaks for MTR4, Z1, Z8, CTCF, and RAD21 were mapped to the hg19 human genome assembly. For genomic distribution analysis, peaks were categorized based on their overlap with annotated genomic features including Transcription Start Sites (TSS), Transcription End Sites (TES), Gene Bodies, and Intergenic regions. 

The TSS and TES regions were defined as 1 bp sites at the start and end of genes, respectively, extended by ±500 bp for overlap analysis. Gene bodies were defined as the region between the TSS and TES. "True Intergenic" regions were defined as genomic areas located outside of annotated gene bodies. 

For ChIP-seq peak-centric analysis, the proportion of peaks overlapping each feature was calculated. Additionally, intergenic peaks were further characterized by their overlap with H3K4me1-marked regions and CAGE-defined enhancers to identify putative distal regulatory elements.

### ncRNA Characterization and Overlap with ChIP-seq
ncRNA transcripts (PROMPTS and eRNAs) were identified and their relationship with ChIP-seq peaks was analyzed. PROMPTS were defined as ncRNA transcripts located within 2 kb upstream of a TSS. eRNAs were defined as ncRNA transcripts overlapping with CAGE-defined enhancer regions. The percentage of TSS regions associated with at least one overlapping ncRNA event was also quantified.

### Bin-based Overlap and Statistical Significance
To ensure an unbiased comparison of overlaps between genomic features of varying sizes and proximity, a bin-based approach was employed. A unified set of genomic bins was created by merging all features of interest (ChIP-seq peaks, ncRNA regions, and regulatory annotations) and applying `GenomicRanges::reduce` with a 500 bp extension around feature centers. This process generated a discrete set of genomic locations where each feature's presence or absence was recorded as a binary value.

Venn diagrams were generated to visualize the intersection of these binary features. Statistical significance of the overlaps was assessed using Fisher’s exact tests performed on the binned data, ensuring that the results were not skewed by peak length or density.

### Hi-C Interaction Analysis and Aggregate Peak Analysis (APA)
Hi-C interaction data for control (SCR) and MTR4-depleted HeLa cells was processed to analyze changes in 3D genomic organization. Interaction matrices were binned at 10 kb resolution. A consolidated metrics database was generated (containing `METRICS.df` and `METRICS_QUANT.df`) to store interaction strengths between all anchor pairs. To evaluate the relationship between transcriptional activity and chromatin interactions, bait regions were defined based on TSS coordinates and annotated with the presence of associated ncRNA (PROMPTS). 

Aggregate Peak Analysis (APA) was performed using a list of 2D sub-matrices (`myMAT`) centered on sets of genomic interacting pairs (e.g., enhancers, TSS, or ChIP-seq peaks). These matrices were averaged to visualize the global enrichment of interaction signal. Statistical significance of the changes in interaction frequencies across different quintiles was assessed using Fisher's exact tests performed on the binned interaction metrics.

## QUANTIFICATION AND STATISTICAL ANALYSIS

All bioinformatic analyses were performed using R (v3.4.2). Genomic range operations, including overlaps and reductions, were conducted using the `GenomicRanges` and `rtracklayer` packages. Hi-C data analysis and APA were performed using custom scripts and functions integrated with `data.table` for efficient processing of large interaction matrices. Data manipulation and statistical testing were performed using `data.table`, `dplyr`, and base R functions. Fisher’s exact tests were used to calculate the significance of overlaps, enrichment across genomic quantiles, and changes in interaction frequencies. Figures were generated using `ggplot2`, `ggpubr`, and `pheatmap`.

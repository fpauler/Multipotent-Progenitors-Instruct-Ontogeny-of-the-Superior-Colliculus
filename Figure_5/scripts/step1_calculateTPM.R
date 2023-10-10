# 2 cell clone analysis
# first step: calculate TPM and prepare a sample list with alignment stats
# checked 10.10.2023

library (openxlsx)
library (GenomicFeatures)
library (GenomicAlignments)
library (scater)

# NOTE: this analysis is partially based on the initial PatchSeq analysis, therefore give initial output folder here
# based on the initial meta data file: initial_meta_data_v1.xlsx

base_initial <- "../Figure_4"
analysis_base <- "./"

# this is the folder with the alignment
gene_count_base <- paste(analysis_base, "TwoCellClonePatchSeq/", sep="" )

#read previously created GRanges objects
mm.introns <- readRDS( file=paste(base_initial, "annotation/gencode.vM27.annotation.GRangesListIntrons.RDS", sep="/" ) )
mm.exons <- readRDS( file=paste(base_initial, "annotation/gencode.vM27.annotation.GRangesListExons.RDS", sep="/" ) )

#read all BAM files
files <- c(list.files(path =gene_count_base, pattern=glob2rx("*.bam"), recursive = T, full.names = T))

#counts reads overlapping exons
mm.exons.counts <- summarizeOverlaps(features = mm.exons, 
                                     reads = files,
                                     singleEnd=TRUE, mode = "IntersectionNotEmpty", ignore.strand = T, inter.feature = T)

#count reads overlapping introns 
mm.intron.counts <- summarizeOverlaps(features = mm.introns, 
                                      reads = files,
                                      singleEnd=TRUE, mode = "IntersectionNotEmpty", ignore.strand = T, inter.feature = T)

# set identical order of rows and columns in exon and intron count matrix
mm.intron.counts <- mm.intron.counts[rownames(mm.exons.counts), colnames(mm.exons.counts)]

#obtain cDNA length - in kb
mm.exon.len <- sum(width(mm.exons)) / 1000

#calculate TPM
exon.matrix.tpm <- calculateTPM(mm.exons.counts, lengths=mm.exon.len)

#obtain total intron length
mm.intron.len <- sum(width(mm.introns)) / 1000
#zero intron length is not possible for TPM calculation add 1bp
mm.intron.len[which(mm.intron.len == 0)] <- 1

#several single exon genes are missing from intron.matrix.tpm
intron.matrix.tpm <- calculateTPM(mm.intron.counts, lengths=mm.intron.len)

#save TPM matrix - sum of itron and exon TPM - like in Cadwell eLife 2020 paper
comb_matrix.tpm <- intron.matrix.tpm + exon.matrix.tpm

col_names <- colnames(comb_matrix.tpm)
col_names <- gsub(pattern = "STAR\\.", replacement = "", x = col_names)
col_names <- gsub(pattern = "\\.Aligned\\.sortedByCoord\\.out\\.bam", replacement = "", x = col_names)

colnames(comb_matrix.tpm) <- col_names 

# save to csv file - this is the file uploaded to GEO
write.csv(x = comb_matrix.tpm, file = paste(analysis_base, "RDS_files/comb_matrix_TPM_exon_plus_intron.csv", sep="/"), row.names = T)

# prepare a sample list
#analyse STAR Logs
files.STAR <- c(list.files(path=gene_count_base, pattern=glob2rx("*.Log.final.out"), recursive = T))

STAR_df2plot <- data.frame()
for (file_name in files.STAR) {
  STAR.Log <- read.table(file = paste(gene_count_base, file_name, sep=""), header = F, fill = T, sep="\t", stringsAsFactors = F)
  
  #isolate alignment folder name
  folders <- strsplit(x = file_name, split = "\\/")[[1]] 
  
  #if the folder name contains a dot this indicates that the sample ID needs to be extracted
  #otherwise the sample_id is in the filename
  sample_name <- strsplit(x = folders[[length(folders)-1]], split = "\\.")[[1]][2]
  
  STAR_df2plot <- rbind(STAR_df2plot, data.frame(sample_name=sample_name, STAR.unique.aligned=STAR.Log$V2[8], 
                                                 STAR.tot=STAR.Log$V2[5], STAR.perc = STAR.Log$V2[9], stringsAsFactors = F))
}

STAR_df2plot$STAR.unique.aligned <- as.numeric( STAR_df2plot$STAR.unique.aligned )
STAR_df2plot$STAR.tot <- as.numeric( STAR_df2plot$STAR.tot )
STAR_df2plot$STAR.perc <- STAR_df2plot$STAR.unique.aligned / STAR_df2plot$STAR.tot

sample_list <- read.xlsx( xlsxFile = paste(analysis_base, "initial_meta_data_v1.xlsx", sep=""))
sample_list <- merge( sample_list, STAR_df2plot, by="sample_name" )

wb <- createWorkbook()
sheetName <- "meta_data"
addWorksheet(wb, sheetName)

writeData(wb, sheetName, sample_list)
addFilter(wb, sheetName, row = 1, cols = 1:ncol(sample_list))
setColWidths(wb, sheetName, cols = 1:ncol(sample_list), widths="auto")

saveWorkbook(wb, file = paste(analysis_base, "/Supplement/initial_meta_data_for_reporting.xlsx", sep=""), overwrite = T) 


sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.6 LTS
# 
# Matrix products: default
# BLAS:   /opt/R/4.0.3/lib/R/lib/libRblas.so
# LAPACK: /opt/R/4.0.3/lib/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_IE.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_IE.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_IE.UTF-8      
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_IE.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] openxlsx_4.2.4              scater_1.18.6               ggplot2_3.3.5               SingleCellExperiment_1.12.0 GenomicAlignments_1.26.0    Rsamtools_2.6.0             Biostrings_2.58.0          
# [8] XVector_0.30.0              SummarizedExperiment_1.20.0 MatrixGenerics_1.2.1        matrixStats_0.60.1          GenomicFeatures_1.42.3      AnnotationDbi_1.52.0        Biobase_2.50.0             
# [15] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7         IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.1        
# 
# loaded via a namespace (and not attached):
#   [1] viridis_0.6.1             httr_1.4.2                BiocSingular_1.6.0        viridisLite_0.4.0         bit64_4.0.5               DelayedMatrixStats_1.12.3 scuttle_1.0.4            
# [8] assertthat_0.2.1          askpass_1.1               BiocFileCache_1.14.0      blob_1.2.2                vipor_0.4.5               GenomeInfoDbData_1.2.4    progress_1.2.2           
# [15] pillar_1.6.2              RSQLite_2.2.8             lattice_0.20-41           beachmat_2.6.4            glue_1.4.2                colorspace_2.0-2          Matrix_1.3-4             
# [22] plyr_1.8.6                XML_3.99-0.8              pkgconfig_2.0.3           biomaRt_2.46.3            zlibbioc_1.36.0           purrr_0.3.4               scales_1.1.1             
# [29] BiocParallel_1.24.1       tibble_3.1.4              openssl_1.4.5             generics_0.1.0            ellipsis_0.3.2            cachem_1.0.6              withr_2.4.2              
# [36] magrittr_2.0.1            crayon_1.4.1              memoise_2.0.0             fansi_0.5.0               xml2_1.3.2                beeswarm_0.4.0            tools_4.0.3              
# [43] prettyunits_1.1.1         hms_1.1.1                 lifecycle_1.0.0           stringr_1.4.0             munsell_0.5.0             zip_2.2.0                 irlba_2.3.3              
# [50] DelayedArray_0.16.3       packrat_0.7.0             compiler_4.0.3            rsvd_1.0.5                rlang_0.4.11              grid_4.0.3                RCurl_1.98-1.4           
# [57] BiocNeighbors_1.8.2       rstudioapi_0.13           rappdirs_0.3.3            bitops_1.0-7              gtable_0.3.0              DBI_1.1.1                 curl_4.3.2               
# [64] reshape2_1.4.4            R6_2.5.1                  gridExtra_2.3             dplyr_1.0.7               rtracklayer_1.50.0        fastmap_1.1.0             bit_4.0.4                
# [71] utf8_1.2.2                ggbeeswarm_0.6.0          stringi_1.7.4             Rcpp_1.0.7                vctrs_0.3.8               dbplyr_2.1.1              tidyselect_1.1.1         
# [78] sparseMatrixStats_1.2.1  

# prepare a GRanges object with exons and introns and performs read counting
# this is the first step in the PatchSeq analysis

# NOTE: to save space all alignment details, including BAM files are stored in this file: Patch_Seq_M27.tar.gz
#       extract prior to running this script

# NOTE: TPM values can be use to do data integration as suggested here:
# https://github.com/satijalab/seurat/issues/2854

#checked 6.10.2023

library (GenomicFeatures)
library(GenomicAlignments)
library(scater)

# define base folder
base <- "./"

get_introns <- function (x) {
  
  tmp <- gaps(x)
  # if one intron is present we have 2 gaps
  if (length(tmp) > 1) {
    # Note: first gap is from chromosome start to first exon - discard
    tmp <- tmp[2:length(tmp)]
    return (tmp)
  } else {
   # no intron: return empty granges object
    GRanges(seqnames= NULL, ranges = NULL, strand = NULL)
  }
  
}

if (!file.exists(paste(base, "annotation/gencode.vM27.annotation.GRangesListIntrons.RDS", sep="/"))) {
  
  # read transcript info from gtf file
  txdb <- makeTxDbFromGFF(file = paste(base, "annotation/gencode.vM27.annotation.gtf", sep="/"))
  # create GRanges list with exon info for each gene
  mm.exons <- GenomicFeatures::exonsBy(txdb, by="gene")
  # merge overlapping exons - can happen if many transcript versions are present for one gene
  mm.exons <- reduce(mm.exons)
  
  #create an intronic GRanges List - long wait!
  mm.introns <- endoapply(mm.exons, get_introns)
  
  saveRDS( object =  mm.introns, file=paste(base, "annotation/gencode.vM27.annotation.GRangesListIntrons.RDS", sep="/" ))
  saveRDS( object =  mm.exons, file=paste(base, "annotation/gencode.vM27.annotation.GRangesListExons.RDS", sep="/" ))
  
}
# stop and restart to free memory - only necessary if you have limited RAM

#read previously created GRanges objects
mm.introns <- readRDS( file=paste(base, "annotation/gencode.vM27.annotation.GRangesListIntrons.RDS", sep="/" ) )
mm.exons <- readRDS( file=paste(base, "annotation/gencode.vM27.annotation.GRangesListExons.RDS", sep="/" ) )

#read all BAM files
files <- c(list.files(path = paste(base, '/Patch_Seq_M27', sep=""), pattern=glob2rx("*.bam"), recursive = T, full.names = T))

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

#save TPM matrix - sum of intron and exon TPM - like in Cadwell eLife 2020 paper
comb_matrix.tpm <- intron.matrix.tpm + exon.matrix.tpm

col_names <- colnames(comb_matrix.tpm)
col_names <- gsub(pattern = "STAR\\.", replacement = "", x = col_names)
col_names <- gsub(pattern = "\\.Aligned\\.sortedByCoord\\.out\\.bam", replacement = "", x = col_names)

colnames(comb_matrix.tpm) <- col_names 

# save to csv file - this is the file uploaded to GEO
write.csv(x = comb_matrix.tpm, file = paste(base, "RDS_files/comb_matrix_TPM_exon_plus_intron.csv", sep="/"), row.names = T)

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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_IE.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_IE.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_IE.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_IE.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scater_1.18.6               ggplot2_3.3.5               SingleCellExperiment_1.12.0 GenomicAlignments_1.26.0    Rsamtools_2.6.0            
# [6] Biostrings_2.58.0           XVector_0.30.0              SummarizedExperiment_1.20.0 MatrixGenerics_1.2.1        matrixStats_0.60.1         
# [11] GenomicFeatures_1.42.3      AnnotationDbi_1.52.0        Biobase_2.50.0              GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
# [16] IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.1        
# 
# loaded via a namespace (and not attached):
#   [1] viridis_0.6.1             httr_1.4.2                BiocSingular_1.6.0        viridisLite_0.4.0         bit64_4.0.5              
# [6] DelayedMatrixStats_1.12.3 scuttle_1.0.4             assertthat_0.2.1          askpass_1.1               BiocFileCache_1.14.0     
# [11] blob_1.2.2                vipor_0.4.5               GenomeInfoDbData_1.2.4    progress_1.2.2            pillar_1.6.2             
# [16] RSQLite_2.2.8             lattice_0.20-41           beachmat_2.6.4            glue_1.4.2                colorspace_2.0-2         
# [21] Matrix_1.3-4              XML_3.99-0.8              pkgconfig_2.0.3           biomaRt_2.46.3            zlibbioc_1.36.0          
# [26] purrr_0.3.4               scales_1.1.1              BiocParallel_1.24.1       tibble_3.1.4              openssl_1.4.5            
# [31] generics_0.1.0            ellipsis_0.3.2            cachem_1.0.6              withr_2.4.2               magrittr_2.0.1           
# [36] crayon_1.4.1              memoise_2.0.0             fansi_0.5.0               xml2_1.3.2                beeswarm_0.4.0           
# [41] tools_4.0.3               prettyunits_1.1.1         hms_1.1.1                 lifecycle_1.0.0           stringr_1.4.0            
# [46] munsell_0.5.0             irlba_2.3.3               DelayedArray_0.16.3       packrat_0.7.0             compiler_4.0.3           
# [51] rsvd_1.0.5                rlang_0.4.11              grid_4.0.3                RCurl_1.98-1.4            BiocNeighbors_1.8.2      
# [56] rstudioapi_0.13           rappdirs_0.3.3            bitops_1.0-7              gtable_0.3.0              DBI_1.1.1                
# [61] curl_4.3.2                R6_2.5.1                  gridExtra_2.3             dplyr_1.0.7               rtracklayer_1.50.0       
# [66] fastmap_1.1.0             bit_4.0.4                 utf8_1.2.2                ggbeeswarm_0.6.0          stringi_1.7.4            
# [71] Rcpp_1.0.7                vctrs_0.3.8               dbplyr_2.1.1              tidyselect_1.1.1          sparseMatrixStats_1.2.1 

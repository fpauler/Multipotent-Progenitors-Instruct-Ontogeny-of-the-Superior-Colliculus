# 2 cell clone analysis
# 2nd step: follow initial analysis to calculate NMS scores
# checked 10.10.2023

library (Seurat)
library (openxlsx)
library (dplyr)
library (ggplot2)
library(ggrepel)

set.seed(2401)

base_initial <- "../Figure_4/"
analysis_base <- "./"

# prepare the PatchSeq dataset
if (!file.exists(paste(analysis_base, "RDS_files/comb_matrix_TPM_exon_plus_intron_Seurat.RDS", sep="/"))) {
  
  #read gene ID conversion table 
  ensmusg_symbol <- readRDS(file = paste(base_initial, "RDS_files/ensmusg_symbol.RDS", sep="/"))
  
  counts.matrix <- read.csv(file = paste(analysis_base, "RDS_files/comb_matrix_TPM_exon_plus_intron.csv", sep="/"), header = T)
  
  # process count matrix for further analyses
  rownames(counts.matrix) <- counts.matrix$X
  counts.matrix <- counts.matrix[,2:ncol(counts.matrix)]
  colnames(counts.matrix) <- gsub(pattern = "X", replacement = "", x = colnames(counts.matrix))
  
 #give unique symbol IDs to TPM count matrix genes
  gn <- sapply(rownames(counts.matrix), function (x) ensmusg_symbol[[x]])
  
  length(gn)
  # 55359
  
  #remove duplicated symbols
  gn <- gn[!(duplicated(gn) | duplicated(gn, fromLast = TRUE))]
  
  length(gn)
  # 52266
  
  counts.matrix <- counts.matrix[names(gn),]
  rownames(counts.matrix) <- as.character(gn[rownames(counts.matrix)])
  
  # NOTE: I omit normalisation as TPM values are already normalized!
  
  # log2 transform count matrix
  counts.matrix <- log2(counts.matrix + 1)
  
  sample_list <- read.xlsx( xlsxFile = paste(analysis_base, "Supplement/initial_meta_data_for_reporting.xlsx", sep="/"), sheet=1)
  
  nrow(sample_list)
  # [1] 36
  
  # all cells meet the criteria for read coverage and alignment rate
  length( which(sample_list$STAR.perc > 0.05 & sample_list$STAR.unique.aligned > 26000) )
  # [1] 36
  
  # create the Seurat object
  PatchSeq_Seurat <- CreateSeuratObject(counts = counts.matrix[,as.character(sample_list$sample_name)], min.cells = 5, min.features = 500)
  
  # normalization already done by TPM!
  
  # do scaling for all genes
  all.genes <- rownames(PatchSeq_Seurat)
  
  PatchSeq_Seurat <- ScaleData(PatchSeq_Seurat, features = all.genes)
  
  rownames(sample_list) <- sample_list$sample_name
  PatchSeq_Seurat <- AddMetaData(object = PatchSeq_Seurat, metadata = sample_list)
  
  # An object of class Seurat 
  # 19272 features across 36 samples within 1 assay 
  # Active assay: RNA (19272 features, 0 variable features)
  
  saveRDS(object = PatchSeq_Seurat, file = paste(analysis_base, "RDS_files/comb_matrix_TPM_exon_plus_intron_Seurat.RDS", sep="/"))
}

if ( !file.exists( paste(analysis_base, "RDS_files/comb_matrix_TPM_exon_plus_intron_Seurat_NMS.RDS", sep="/") ) ) {
  
  # read in data from first PatchSeq analysis
  MBd_Seurat <- readRDS(file = paste(base_initial, "RDS_files/MBd_Seurat_wholeMBD_final.RDS", sep="/"))
  all.markers <-  readRDS( file = paste(base_initial, "RDS_files/MBd_Seurat_wholeMBD_RAW_markers.RDS", sep="/") )
  
  # read Seurat object prepared above
  PatchSeq_Seurat <- readRDS(file = paste(analysis_base, "RDS_files/comb_matrix_TPM_exon_plus_intron_Seurat.RDS", sep="/"))
  
  all.markers %>%
    group_by(cluster) %>%
    top_n(n = 200, wt = avg_log2FC) -> top_markers
  
  cell_types <- sapply(unique(MBd_Seurat@meta.data$tax_3), function (cell_type) {
    length(rownames(MBd_Seurat@meta.data[which(MBd_Seurat@meta.data$tax_3 == cell_type),]))
  })
  
  #remove small cell types 
  cell_types <- cell_types[which(cell_types > 100)]
  
  # Oligodendrocytes Di- and mesencephalon neurons          Astroependymal cells                Vascular cells 
  #             4366                          4253                          3682                          1759 
  
  # Immune cells 
  #          939 
  
  # calculate the mean expression for each gene in each cell type
  medians <- sapply( names(cell_types), function (cell_type) {
    apply(MBd_Seurat[["RNA"]]@data[unique(top_markers$gene),
                                   rownames(MBd_Seurat@meta.data[which(MBd_Seurat@meta.data$tax_3 == cell_type),])], 1, mean)
  })
  
  length(unique(top_markers$gene))
  # 986 genes
  
  # remove genes where the median is 0
  med2plot <- medians
  med2plot <- med2plot[which(apply(med2plot, 1, max) > 0),]
  
  #check min number of markers per tissue
  min_n_gl <- min(table(apply(med2plot, 1, function (x) which(x == max(x)))))
  # [1] 188
  
  #determine cell type with highest expression
  marker_list <- apply(med2plot, 1, function (x) which(x == max(x)))
  marker_list_df <- data.frame(gene = names(marker_list), idx = as.numeric(marker_list))
  marker_list_df$cell_type <- colnames(med2plot)[marker_list_df$idx]
  
  #cut the gene list by taking most expressed genes
  cut_gl <- lapply(unique(marker_list_df$cell_type), function (cell_type) {
    sort(med2plot[marker_list_df[which(marker_list_df$cell_type == cell_type), "gene"], cell_type], decreasing = T)[1:min_n_gl]
  })
  
  # prepare a character vector
  cut_gl <- names(unlist(cut_gl))
  
  # filter marker list 
  marker_list_df <- marker_list_df[which(marker_list_df$gene %in% cut_gl),]
  
  #determine mean marker expression in reference cell type
  ref_med_expr <- sapply( unique(marker_list_df$cell_type), function (cell_types) {
    tmp <- apply(MBd_Seurat[["RNA"]]@data[ marker_list_df$gene[which(marker_list_df$cell_type == cell_types)],
                                           rownames(MBd_Seurat@meta.data[which(MBd_Seurat@meta.data$tax_3 == cell_types),]) ], 1, mean)
    return(mean(tmp))
  })
  
  # clean up marker list to only contain genes also informative in PatchSeq
  marker_list_df <- marker_list_df[which(marker_list_df$gene %in% rownames(PatchSeq_Seurat[["RNA"]]@data)),]
  
  # compare over all expression levels in ref and PatchSeq
  ref_expr <- apply(MBd_Seurat[["RNA"]]@data, 1, mean)
  PatchSeq_expr <- apply(PatchSeq_Seurat[["RNA"]]@data, 1, mean)
  
  rel_expr <- PatchSeq_expr[marker_list_df$gene] / ref_expr[marker_list_df$gene]
  
  mean(rel_expr)
  #[1] 13.15191
  
  # expression matrix from PatchSeq with marker genes
  lPS <- PatchSeq_Seurat[["RNA"]]@data[marker_list_df$gene,]
  
  # calculate mean marker expression for each cell
  mean_marker <- sapply(colnames(lPS), function (cell){
    
    sapply ( unique(marker_list_df$cell_type), function (x) {
      tmp_counts <- lPS[marker_list_df[which(marker_list_df$cell_type == x), "gene"], cell]
      mean(tmp_counts)
    })
  })
  
  #now the NMS calculation
  for (cell_type in unique(marker_list_df$cell_type) ) {
    mean_marker[cell_type,] <- mean_marker[cell_type,] / ref_med_expr[cell_type]
  }
  
  # renaming vector
  ren_vec <- c("NMS_Oligo", "NMS_Vascular", "NMS_Immune", "NMS_Neurons", "NMS_Astro")
  names(ren_vec) <- c("Oligodendrocytes", "Vascular cells", "Immune cells", 
                      "Di- and mesencephalon neurons", "Astroependymal cells")
  
  #prepare a df with all NMS scores for each cell type
  NMS_df <- data.frame(sample_id = colnames( mean_marker ), NMS =  mean_marker[1,] )
  colnames(NMS_df) <- c("sample_id", ren_vec[rownames(mean_marker)[1]])
  
  for (i in 2: nrow(mean_marker)) {
    tmp <- data.frame(sample_id = colnames( mean_marker ), NMS =  mean_marker[i,] )
    colnames(tmp) <- c("sample_id", ren_vec[rownames(mean_marker)[i]])
    NMS_df <- merge(NMS_df, tmp, by="sample_id")
  }
  
  rownames(NMS_df) <- NMS_df$sample_id
  PatchSeq_Seurat <- AddMetaData(object = PatchSeq_Seurat, metadata = NMS_df[2:ncol(NMS_df)])
  
  saveRDS(object = PatchSeq_Seurat, file = paste(analysis_base, "RDS_files/comb_matrix_TPM_exon_plus_intron_Seurat_NMS.RDS", sep="/"))
  
} else {
  PatchSeq_Seurat <- readRDS( file = paste(analysis_base, "RDS_files/comb_matrix_TPM_exon_plus_intron_Seurat_NMS.RDS", sep="/") )
}

# identify the maximum NMS score from non-neuronal (aka contaminating) t-types
NMS_cont <- apply( PatchSeq_Seurat@meta.data[, c("NMS_Oligo", "NMS_Vascular", "NMS_Immune", "NMS_Astro")], 1, max)

# plot contaminating vs neuron NMS scores
df2plot <- data.frame ( NMS_neuron = PatchSeq_Seurat@meta.data$NMS_Neurons, NMS_cont = as.numeric(NMS_cont), sample_id = names(NMS_cont))
ggplot( df2plot, aes ( x= NMS_neuron, y = NMS_cont, label=sample_id)) + geom_point() + 
  xlim(0,6) + ylim(0,6) + geom_text_repel()
ggsave(filename = paste(analysis_base, "QC/NMS_score_plot.pdf", sep="/"))

# I do no filtering here - looks reasonably OK

# do SCTransform to correct for alignment rate, as done before
PatchSeq_Seurat_SCT <- SCTransform(object = PatchSeq_Seurat, vars.to.regress = c("STAR.perc"), 
                                   residual.features = rownames(PatchSeq_Seurat) )

# save the file
saveRDS( object = PatchSeq_Seurat_SCT, file = paste(analysis_base, "RDS_files/PatchSeq_Seurat_TPM_SCT.RDS", sep="/"))


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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggrepel_0.9.1      ggplot2_3.3.5      dplyr_1.0.7        openxlsx_4.2.4     SeuratObject_4.0.2 Seurat_4.0.4.9000 
# 
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-149          matrixStats_0.60.1    spatstat.sparse_2.0-0 RcppAnnoy_0.0.19      RColorBrewer_1.1-2    httr_1.4.2            sctransform_0.3.2     tools_4.0.3           utf8_1.2.2           
# [10] R6_2.5.1              irlba_2.3.3           rpart_4.1-15          KernSmooth_2.23-17    uwot_0.1.10           mgcv_1.8-33           DBI_1.1.1             lazyeval_0.2.2        colorspace_2.0-2     
# [19] withr_2.4.2           tidyselect_1.1.1      gridExtra_2.3         compiler_4.0.3        plotly_4.9.4.1        labeling_0.4.2        scales_1.1.1          lmtest_0.9-38         spatstat.data_2.1-0  
# [28] ggridges_0.5.3        pbapply_1.5-0         goftest_1.2-2         stringr_1.4.0         digest_0.6.27         spatstat.utils_2.2-0  pkgconfig_2.0.3       htmltools_0.5.2       parallelly_1.28.1    
# [37] fastmap_1.1.0         htmlwidgets_1.5.4     rlang_0.4.11          shiny_1.6.0           farver_2.1.0          generics_0.1.0        zoo_1.8-9             jsonlite_1.7.2        ica_1.0-2            
# [46] zip_2.2.0             magrittr_2.0.1        patchwork_1.1.1       Matrix_1.3-4          Rcpp_1.0.7            munsell_0.5.0         fansi_0.5.0           abind_1.4-5           reticulate_1.20      
# [55] lifecycle_1.0.0       stringi_1.7.4         MASS_7.3-53           Rtsne_0.15            plyr_1.8.6            grid_4.0.3            parallel_4.0.3        listenv_0.8.0         promises_1.2.0.1     
# [64] crayon_1.4.1          miniUI_0.1.1.1        deldir_0.2-10         lattice_0.20-41       cowplot_1.1.1         splines_4.0.3         tensor_1.5            pillar_1.6.2          igraph_1.2.6         
# [73] spatstat.geom_2.2-2   future.apply_1.8.1    reshape2_1.4.4        codetools_0.2-16      leiden_0.3.9          glue_1.4.2            packrat_0.7.0         data.table_1.14.0     png_0.1-7            
# [82] vctrs_0.3.8           httpuv_1.6.3          gtable_0.3.0          RANN_2.6.1            purrr_0.3.4           spatstat.core_2.3-0   polyclip_1.10-0       tidyr_1.1.3           scattermore_0.7      
# [91] future_1.22.1         assertthat_0.2.1      mime_0.11             xtable_1.8-4          later_1.3.0           survival_3.2-7        viridisLite_0.4.0     tibble_3.1.4          cluster_2.1.0        
# [100] globals_0.14.0        fitdistrplus_1.1-5    ellipsis_0.3.2        ROCR_1.0-11  

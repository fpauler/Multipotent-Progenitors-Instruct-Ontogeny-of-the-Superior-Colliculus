# second step of the Pten analysis
# produces meta data linking the Pten P0 data to embryo trajectory data

# checked 4.10.2023

library (Seurat)
library (ggplot2)

set.seed (2401)

# define base folder
base <- "/."

  # classify based on Figure 1 data - combined Midbrain specific trajectory
  embryonic_ref <- readRDS( file = paste(base, "/RDS_files/LaManno_embryonic.combined.CCscored.RDS", sep=""))
  # first do a low-res clustering
  embryonic_ref <- FindClusters(embryonic_ref, resolution = 0.2)
  
  # since we focus only on the endpoints of neuronal production, remove RGPs and immature neurons
  embryonic_ref <- subset(embryonic_ref, idents = c("0", "6", "11"), invert = T)
    
  # read the raw Pten data
  Pten_seurat_list <- readRDS ( paste(base, "/RDS_files/P0_Pten_ctrl_RAW.RDS", sep="" ) )
  
  # read the neuron only analysis
  Pten.neuron.seurat <- readRDS( file = paste(base, "/RDS_files/Pten.neuron.seurat.RDS", sep="" ) )
  
  ctrl_cellID <- rownames( Pten.neuron.seurat@meta.data[which(Pten.neuron.seurat@meta.data$genotype == "ctrl"),] )
  ctrl_cellID <- gsub( pattern = "_1", replacement = "", x = ctrl_cellID)
  sparse_cellID <- rownames( Pten.neuron.seurat@meta.data[which(Pten.neuron.seurat@meta.data$genotype == "sparse"),] )
  sparse_cellID <- gsub( pattern = "_2", replacement = "", x = sparse_cellID)
  
  if (!file.exists(paste(base, "/RDS_files/Pten_ctrl_anchors.RDS", sep=""))) {
    embryo.anchors <- FindTransferAnchors(reference = embryonic_ref, query = subset(Pten_seurat_list[["ctrl"]], cells = ctrl_cellID), dims = 1:30)
    saveRDS( embryo.anchors, paste(base, "/RDS_files/Pten_ctrl_anchors.RDS", sep="") )

    # check that clustering stays the same
    umap_plot <- DimPlot( embryonic_ref, label = T) + NoLegend()
    ggsave ( paste(  base, "/QC/embryo_umap_Pten_integration_1.pdf", sep="" ))
 
    stop( "first step done, rstart")
  }
  
  if (!file.exists(paste(base, "/RDS_files/Pten_sparse_anchors.RDS", sep=""))) {
    embryo.anchors <- FindTransferAnchors(reference = embryonic_ref, query = subset(Pten_seurat_list[["sparse"]], cells = sparse_cellID), dims = 1:30)
    saveRDS( embryo.anchors, paste(base, "/RDS_files/Pten_sparse_anchors.RDS", sep="") )

    # check that clustering stays the same
    umap_plot <- DimPlot( embryonic_ref, label = T) + NoLegend()
    ggsave ( paste(  base, "/QC/embryo_umap_Pten_integration_2.pdf", sep="" ))
 
    stop( "second step done, rstart")
  }

# check that clustering stays the same
umap_plot <- DimPlot( embryonic_ref, label = T) + NoLegend()
ggsave ( paste(  base, "/QC/embryo_umap_Pten_integration_3.pdf", sep="" ))

ctrl.anchors <- readRDS( paste(base, "/RDS_files/Pten_ctrl_anchors.RDS", sep="") )
pred_ctrl <- TransferData(anchorset = ctrl.anchors, refdata = embryonic_ref$seurat_clusters, dims = 1:30)
pred_ctrl <- pred_ctrl[1]
colnames(pred_ctrl) <- "embryo_traj"
rownames(pred_ctrl) <- paste( rownames(pred_ctrl), "1", sep="_" )

sparse.anchors <- readRDS( paste(base, "/RDS_files/Pten_sparse_anchors.RDS", sep="") )
pred_sparse <- TransferData(anchorset = sparse.anchors, refdata = embryonic_ref$seurat_clusters, dims = 1:30)
pred_sparse <- pred_sparse[1]
colnames(pred_sparse) <- "embryo_traj"
rownames(pred_sparse) <- paste( rownames(pred_sparse), "2", sep="_" )

saveRDS( rbind( pred_ctrl, pred_sparse), paste(base, "/RDS_files/Pten_embryonic_link.RDS", sep="") )

sessionInfo()


# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.6 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_IE.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_IE.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_IE.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_IE.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.3.6      sp_1.5-0           SeuratObject_4.1.2 Seurat_4.2.0      
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.16            colorspace_2.0-3      deldir_1.0-6          ellipsis_0.3.2        ggridges_0.5.4        rstudioapi_0.14       spatstat.data_3.0-0  
# [8] leiden_0.4.3          listenv_0.8.0         farver_2.1.1          ggrepel_0.9.1         fansi_1.0.3           codetools_0.2-18      splines_4.2.1        
# [15] polyclip_1.10-4       jsonlite_1.8.3        packrat_0.8.1         ica_1.0-3             cluster_2.1.3         png_0.1-7             rgeos_0.5-9          
# [22] pheatmap_1.0.12       uwot_0.1.14           shiny_1.7.3           sctransform_0.3.5     spatstat.sparse_3.0-0 compiler_4.2.1        httr_1.4.4           
# [29] Matrix_1.5-1          fastmap_1.1.0         lazyeval_0.2.2        cli_3.4.1             later_1.3.0           htmltools_0.5.3       tools_4.2.1          
# [36] igraph_1.3.5          gtable_0.3.1          glue_1.6.2            RANN_2.6.1            reshape2_1.4.4        dplyr_1.0.10          Rcpp_1.0.9           
# [43] scattermore_0.8       vctrs_0.5.0           nlme_3.1-157          progressr_0.11.0      lmtest_0.9-40         spatstat.random_2.2-0 stringr_1.4.1        
# [50] globals_0.16.1        mime_0.12             miniUI_0.1.1.1        lifecycle_1.0.3       irlba_2.3.5.1         goftest_1.2-3         future_1.28.0        
# [57] MASS_7.3-57           zoo_1.8-11            scales_1.2.1          spatstat.core_2.4-4   ragg_1.2.4            promises_1.2.0.1      spatstat.utils_3.0-1 
# [64] parallel_4.2.1        RColorBrewer_1.1-3    reticulate_1.26       pbapply_1.5-0         gridExtra_2.3         rpart_4.1.16          stringi_1.7.8        
# [71] rlang_1.0.6           pkgconfig_2.0.3       systemfonts_1.0.4     matrixStats_0.62.0    lattice_0.20-45       ROCR_1.0-11           purrr_0.3.5          
# [78] tensor_1.5            patchwork_1.1.2       htmlwidgets_1.5.4     labeling_0.4.2        cowplot_1.1.1         tidyselect_1.2.0      parallelly_1.32.1    
# [85] RcppAnnoy_0.0.20      plyr_1.8.7            magrittr_2.0.3        R6_2.5.1              generics_0.1.3        pillar_1.8.1          withr_2.5.0          
# [92] mgcv_1.8-40           fitdistrplus_1.1-8    survival_3.3-1        abind_1.4-5           tibble_3.1.8          future.apply_1.9.1    KernSmooth_2.23-20   
# [99] utf8_1.2.2            spatstat.geom_3.0-3   plotly_4.10.0         grid_4.2.1            data.table_1.14.4     digest_0.6.30         xtable_1.8-4         
# [106] tidyr_1.2.1           httpuv_1.6.6          textshaping_0.3.6     munsell_0.5.0         viridisLite_0.4.1   

# various plots and analyses for revision
# checked 4.10.2023

library (Seurat)
library (ggplot2)

# define base folder
base <- "./"

################################################
##
# plot Pten expression during RGP development
##
################################################

# read the data from Figure 1
neuron_analysis_base <- "./"
integrated_data <- readRDS(file = paste(neuron_analysis_base, "/RDS_files/LaManno_embryonic.combined.CCscored.RDS", sep=""))

Idents(integrated_data) <- "cell_type"
RGP_seurat <- subset(integrated_data, idents = "Radial glia")

# An object of class Seurat 
# 24755 features across 4527 samples within 2 assays 
# Active assay: integrated (2000 features, 2000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

# extract raw data for RGPs to redo analysis and plotting
rgp_raw <- GetAssayData( object = RGP_seurat, slot = "count", assay = "RNA")

RGP_Seurat <- CreateSeuratObject(counts = rgp_raw, min.cells = 5, min.features = 500)
RGP_Seurat <- NormalizeData( RGP_Seurat)
RGP_Seurat <- FindVariableFeatures( RGP_Seurat )
all_genes <- rownames(RGP_Seurat)
RGP_Seurat <- ScaleData(object = RGP_Seurat, features = all_genes)
RGP_Seurat <- AddMetaData( object = RGP_Seurat, metadata = RGP_seurat@meta.data)

# An object of class Seurat 
# 16524 features across 4527 samples within 1 assay 
# Active assay: RNA (16524 features, 2000 variable features)

RGP_Seurat@meta.data$Age <- gsub( pattern = "e9.0", replacement = "e09.0", x = RGP_Seurat@meta.data$Age )

RGP_Seurat$num_age <- as.numeric ( gsub(pattern = "e", replacement = "", x = RGP_Seurat$Age) )
RGP_Seurat$cor_age <- ceiling( as.numeric(RGP_Seurat$num_age) )

# plot a corrected age tag
Idents( RGP_Seurat ) <- "cor_age"
VlnPlot( object = subset( x = RGP_Seurat, idents= c(10:14)), features = "Pten", group.by = "cor_age", assay = "RNA", slot = "data") + NoLegend()
ggsave(filename = paste(base, "Supplement/FigureS8A.pdf", sep=""))


################################################
##
# plot GABAergic markers in adult inhibitory ttypes
##
################################################

#PatchSeq base
base <- "../Figure_4"
figure_base <- "./"

# read raw data from Zeisel Midbrain Neurons
cti_Seurat <- readRDS(file = paste(base, "RDS_files/cti_Seurat_RAW.RDS", sep=""))
Idents(cti_Seurat) <- "tax_5"

# convert Zeisel nomenclature to the one used in the paper
Zeisel_to_Cheung <- sapply( c(1:10), function (x) paste("SCINH", x, sep=""))
names( Zeisel_to_Cheung ) <- c("MEINH2", "MEINH3", "MEINH5", "MEINH6", "MEINH7", "MEINH8", "MEINH9", "MEINH10", "MEINH11", "MEINH12")

Inh_Seurat <- subset( x = cti_Seurat, idents = names(Zeisel_to_Cheung))
# An object of class Seurat 
# 15081 features across 2396 samples within 1 assay 
# Active assay: RNA (15081 features, 0 variable features)

Inh_Seurat <- NormalizeData( Inh_Seurat)
Inh_Seurat <- FindVariableFeatures( Inh_Seurat )
all_genes <- rownames(Inh_Seurat)
Inh_Seurat <- ScaleData(object = Inh_Seurat, features = all_genes)
Inh_Seurat <- RunPCA(Inh_Seurat, npcs = 30, verbose = FALSE)

ElbowPlot( Inh_Seurat, ndims = 30)

Inh_Seurat <- RunUMAP(Inh_Seurat, reduction = "pca", dims = 1:15)

Inh_Seurat@meta.data$CheungID <- Zeisel_to_Cheung[ as.character(Inh_Seurat@meta.data$tax_5) ]

all_inh_umap <- DimPlot( object = Inh_Seurat, reduction = "umap", group.by = "CheungID", label=T, label.size = 7 ) + NoLegend()
all_inh_markers <- FeaturePlot( object = Inh_Seurat, features = c("Sst", "Pvalb", "Calb2", "Gad2"), 
             order = T, min.cutoff = "q10", max.cutoff = "q95", slot = "scale.data")

all_inh_umap / all_inh_markers
ggsave(filename = paste(figure_base, "Figures/FigureS8IJ.pdf", sep=""), height = 10)

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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_IE.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_IE.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_IE.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_IE.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.3.6      sp_1.5-0           SeuratObject_4.1.2 Seurat_4.2.0      
# 
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-157          matrixStats_0.62.0    spatstat.sparse_3.0-0 RcppAnnoy_0.0.20      RColorBrewer_1.1-3    httr_1.4.4           
# [7] sctransform_0.3.5     tools_4.2.1           utf8_1.2.2            R6_2.5.1              irlba_2.3.5.1         rpart_4.1.16         
# [13] KernSmooth_2.23-20    uwot_0.1.14           mgcv_1.8-40           rgeos_0.5-9           lazyeval_0.2.2        colorspace_2.0-3     
# [19] withr_2.5.0           tidyselect_1.2.0      gridExtra_2.3         compiler_4.2.1        progressr_0.11.0      textshaping_0.3.6    
# [25] cli_3.4.1             plotly_4.10.0         labeling_0.4.2        scales_1.2.1          lmtest_0.9-40         spatstat.data_3.0-0  
# [31] ggridges_0.5.4        pbapply_1.5-0         systemfonts_1.0.4     goftest_1.2-3         stringr_1.4.1         digest_0.6.30        
# [37] spatstat.utils_3.0-1  pkgconfig_2.0.3       htmltools_0.5.3       parallelly_1.32.1     fastmap_1.1.0         htmlwidgets_1.5.4    
# [43] rlang_1.0.6           rstudioapi_0.14       shiny_1.7.3           farver_2.1.1          generics_0.1.3        zoo_1.8-11           
# [49] jsonlite_1.8.3        spatstat.random_2.2-0 ica_1.0-3             dplyr_1.0.10          magrittr_2.0.3        patchwork_1.1.2      
# [55] Matrix_1.5-1          Rcpp_1.0.9            munsell_0.5.0         fansi_1.0.3           abind_1.4-5           reticulate_1.26      
# [61] lifecycle_1.0.3       stringi_1.7.8         MASS_7.3-57           Rtsne_0.16            plyr_1.8.7            grid_4.2.1           
# [67] parallel_4.2.1        listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1         deldir_1.0-6          miniUI_0.1.1.1       
# [73] lattice_0.20-45       cowplot_1.1.1         splines_4.2.1         tensor_1.5            pillar_1.8.1          igraph_1.3.5         
# [79] spatstat.geom_3.0-3   future.apply_1.9.1    reshape2_1.4.4        codetools_0.2-18      leiden_0.4.3          glue_1.6.2           
# [85] packrat_0.8.1         data.table_1.14.4     png_0.1-7             vctrs_0.5.0           httpuv_1.6.6          polyclip_1.10-4      
# [91] gtable_0.3.1          RANN_2.6.1            purrr_0.3.5           spatstat.core_2.4-4   tidyr_1.2.1           scattermore_0.8      
# [97] future_1.28.0         mime_0.12             xtable_1.8-4          later_1.3.0           ragg_1.2.4            survival_3.3-1       
# [103] viridisLite_0.4.1     tibble_3.1.8          cluster_2.1.3         globals_0.16.1        fitdistrplus_1.1-8    ellipsis_0.3.2       
# [109] ROCR_1.0-11 

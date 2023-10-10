# analyse glia celly in Midbrain 
# creates plots for Figure 3
# checked 4.10.2023

library (Seurat)
library (GenBinomApps)
library (ggplot2)
library (openxlsx)

set.seed(2401)

# define an output folder
base <- "./"

if (!file.exists(paste (base, "/RDS_files/Frz10_seurat.RDS", sep=""))) {
  
  neuron_analysis_base <- "../Figure_1/"
  
  # this part is to isolate RGPs and immature glia from the Frz10 analysis
  # initially all Frz10 cells were integrated with RGP/Neurons from LaManno
  # here I test whether also immature glia are present in the data and clusterd with RGP in the initial analysis
  integrated_data <- readRDS(file = paste(neuron_analysis_base, "/RDS_files/LaManno_embryonic.combined.CCscored.RDS", sep=""))
  
  # check whether glia progenitors are present in the data
  # note that the Frz data was never annotated/filtered for glia progenitors
  DimPlot(object = integrated_data, reduction = "umap", group.by = "seurat_clusters", label=T) + NoLegend()
  ggsave( filename = paste(base, "QC/initial_clustering.pdf", sep=""))
  FeaturePlot(object = integrated_data, features = c("Aldh1l1", "Sox9", "Olig2", "Olig1"), order = T, min.cutoff = "q10")
  ggsave( filename = paste(base, "QC/initial_marker_expression.pdf", sep=""))
  
  # this analysis shows that glial markers are mainly in 4 clusters
  # NOTE: since I filtered glia cell out of the LaManno data, these can only be Frz10 cells
  rgp_clusters <- subset(x = integrated_data, idents = c("4", "21", "22", "24"))
  
  # extract Frz10 cells
  Idents(rgp_clusters) <- "Age"
  Frz_cells <- subset(rgp_clusters, idents = "embryonic_pool")
  Frz_cells_count_matrix <- GetAssayData( object = Frz_cells, slot = "counts", assay = "RNA")
  
  # create a new Seurat object
  Frz_cells_Seurat <- CreateSeuratObject(counts = Frz_cells_count_matrix, min.cells = 5, min.features = 500)
  Frz_cells_Seurat <- NormalizeData( Frz_cells_Seurat)
  Frz_cells_Seurat <- FindVariableFeatures( Frz_cells_Seurat )
  
  # save to free memory
  saveRDS( object = Frz_cells_Seurat, file = paste (base, "/RDS_files/Frz10_seurat.RDS", sep=""))
}

Frz_cells_Seurat <-readRDS( file = paste (base, "/RDS_files/Frz10_seurat.RDS", sep="") )
# An object of class Seurat 
# 15592 features across 1072 samples within 1 assay 
# Active assay: RNA (15592 features, 2000 variable features)

# read the glia data prepared before data
LaManno <- readRDS( paste (base, "RDS_files/seurat.LaManno.glia.RDS", sep="") )
# An object of class Seurat 
# 17358 features across 6374 samples within 1 assay 
# Active assay: RNA (17358 features, 0 variable features)

# Run the standard workflow for visualization and clustering
LaManno <- NormalizeData(LaManno)
LaManno <- FindVariableFeatures(LaManno)

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list(LaManno, Frz_cells_Seurat))
anchors <- FindIntegrationAnchors(object.list = list(LaManno, Frz_cells_Seurat), anchor.features = features)
# this command creates an 'integrated' data assay
seurat.combined <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seurat.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat.combined <- ScaleData(seurat.combined, verbose = FALSE)
seurat.combined <- RunPCA(seurat.combined, npcs = 30, verbose = FALSE)

ElbowPlot( seurat.combined, ndims = 30)

seurat.combined <- RunUMAP(seurat.combined, reduction = "pca", dims = 1:15)
seurat.combined <- FindNeighbors(seurat.combined, reduction = "pca", dims = 1:30)
seurat.combined <- FindClusters(seurat.combined, resolution = 1)

seurat.combined$origin <- "reference"
seurat.combined$origin[which( is.na(seurat.combined$Age))] <- "Frz10"

# create numeric age for grouping
seurat.combined$num_age <- as.numeric ( gsub(pattern = "e", replacement = "", x = seurat.combined$Age) )

# add a broad age - early/late for simpler plotting
seurat.combined$age_group <- "early"
seurat.combined$age_group [ which( seurat.combined$num_age > 12.5 ) ] <- "late"

# extract UMAP coordinated and prepare plots
umap_df <- as.data.frame( Embeddings(seurat.combined, reduction = "umap") )
umap_df$cell_type <- seurat.combined$cell_type
umap_df$age_group <- seurat.combined$age_group
umap_df$origin <- seurat.combined$origin

cell_type_plot <- ggplot( ) + 
  geom_point( data = umap_df[which( !is.na(umap_df$cell_type) & umap_df$age_group == "late"),], aes(x=UMAP_1, y=UMAP_2, color=cell_type ), size = 0.5 ) +
  geom_point( data = umap_df[which( !is.na(umap_df$cell_type) & umap_df$age_group == "early"),], aes(x=UMAP_1, y=UMAP_2, color=cell_type ), size = 0.5 ) +
  scale_color_manual( values = c( "Glioblast" = "#2D52A4", "Radial glia" = "#FF8080")) +
  theme_classic()

age_plot <- ggplot( ) + 
  geom_point( data = umap_df[which( !is.na(umap_df$cell_type) & umap_df$age_group == "late"),], aes(x=UMAP_1, y=UMAP_2, color=age_group ), size = 0.5 ) +
  geom_point( data = umap_df[which( !is.na(umap_df$cell_type) & umap_df$age_group == "early"),], aes(x=UMAP_1, y=UMAP_2, color=age_group ), size = 0.5 ) +
  scale_color_manual( values = c( "early" = "darkviolet", "late" = "goldenrod")) +
  theme_classic()

overlap_plot <- ggplot( umap_df, aes(x=UMAP_1, y=UMAP_2, color=origin )) + 
  geom_point( size = 0.5 ) + 
  scale_color_manual( values = c( "reference" = "grey80", "Frz10" = "black")) +
  theme_classic()

expression_plot <- FeaturePlot(object = seurat.combined, features = c("Mki67", "Sox9", "Aldh1l1", "Olig2"), order = T, min.cutoff = "q10", max.cutoff = "q95") & NoLegend()

# save plots separately for processing in the paper
cell_type_plot <- cell_type_plot + theme(legend.position = "none") + coord_fixed()
ggsave(filename = paste(base, "/Figures/Figure3B_NoLegend.png", sep=""), width = 6, plot = cell_type_plot)
ggsave(filename = paste(base, "/Figures/Figure3B_NoLegend.pdf", sep=""), width = 6, plot = cell_type_plot)

age_plot <- age_plot + theme(legend.position = "none") + coord_fixed()
ggsave(filename = paste(base, "/Figures/Figure3C_NoLegend.png", sep=""), width = 6, plot = age_plot)
ggsave(filename = paste(base, "/Figures/Figure3C_NoLegend.pdf", sep=""), width = 6, plot = age_plot)

ggsave(filename = paste(base, "/Figures/Figure3D_NoLegend.png", sep=""), width = 6, plot = expression_plot)
ggsave(filename = paste(base, "/Figures/Figure3D_NoLegend.pdf", sep=""), width = 6, plot = expression_plot)

overlap_plot <- overlap_plot + theme(legend.position = "none") + coord_fixed()
ggsave(filename = paste(base, "/Figures/Figure3E_NoLegend.png", sep=""), width = 6, plot = overlap_plot)
ggsave(filename = paste(base, "/Figures/Figure3E_NoLegend.pdf", sep=""), width = 6, plot = overlap_plot)

# redundant to first script - to determine total number of cells in LaManno data
library (loomR)

# do the same filtering as before to focus on Dorsal Midbrain

# read location information - defines which cell ś to keep
location2keep <- read.table(paste(base, "other_data/LaManno.Location-Giselle.csv", sep=""), sep="\t", header = T)
location2keep <- location2keep[which(location2keep$use == "yes"),]

# process loom file
lfile <- connect(filename = paste(base, "other_data/dev_all.loom", sep=""), mode = "r+", skip.validate = T)

# get samples indices from Midbrain or Dorsal Midbrain
mbdIDX <- which( lfile$col.attrs$Tissue[] == "Midbrain" | lfile$col.attrs$Tissue[] == "MidbrainDorsal" )
mbdIDX <- mbdIDX[which(!(duplicated(lfile$col.attrs$CellID[mbdIDX]) | duplicated(lfile$col.attrs$CellID[mbdIDX], fromLast = T)))]

# additional filtering based on cluster annotation
sclassMBD <- which( grepl( x=lfile$col.attrs$Subclass[], pattern="midbrain", ignore.case = T) | grepl( x=lfile$col.attrs$Subclass[], pattern="Mixed" ) )
mbdIDX <- intersect( mbdIDX, sclassMBD)

# identify which cells from e9-11 should be removed
locIDX <-  which( !lfile$col.attrs$Location[] %in% c(location2keep$location) )
locIDX <- intersect( locIDX, which( lfile$col.attrs$Age[] %in% c("e9.0", "e10.0", "e11.0")))

# remove cells that are not from the location of interest
cti_MBd <- setdiff(mbdIDX, locIDX)

# count all cells per age
cell_count_df <- reshape2::melt( table( lfile$col.attrs$Age[cti_MBd] ) )

# merge with glioblast count per age
cell_count_df <- merge( cell_count_df, reshape2::melt( table(LaManno$Age[which(LaManno$cell_type == "Glioblast")]) ), by="Var1", all=T )

# make simpler age groups
colnames(cell_count_df) <- c("age", "total", "glia")
cell_count_df$age <- gsub( pattern = "e", replacement = "", x = cell_count_df$age)
cell_count_df$cor_age <- ceiling( as.numeric(cell_count_df$age) )
cell_count_df$glia[ which(is.na(cell_count_df$glia))] <- 0

# do the plotting
df2plot <- data.frame()
for (age in unique(cell_count_df$cor_age)) {
  tmp <- cell_count_df[ which( cell_count_df$cor_age == age ), ]
  df2plot <- rbind( df2plot, data.frame( age = unique(tmp$cor_age), total = sum(tmp$total), glia = sum(tmp$glia) ) )
}

df2plot$rel <- ( df2plot$glia / df2plot$total ) * 100

# add confidence interval as in Cadwell et al., Elife 2020
conf_int <- t(sapply (1:nrow(df2plot), function (x) {
  tmp <- clopper.pearson.ci(df2plot$glia[x], df2plot$total[x], alpha = 0.05, CI = "two.sided")
}))

df2plot <- cbind(df2plot, conf_int)

df2plot$Lower.limit <- as.numeric(df2plot$Lower.limit) * 100
df2plot$Upper.limit <- as.numeric(df2plot$Upper.limit) * 100

all_age_rel_plot <- ggplot(df2plot, aes(x=age, y=rel)) + geom_bar(stat="identity") + 
                      theme_classic() + ggtitle( "%glia in data" ) + ylab("") + xlab("age rounded to next full day") +
                      geom_errorbar(aes(ymin=Lower.limit, ymax=Upper.limit))

ggsave(plot = all_age_rel_plot, filename = paste(base, "/Figures/Figure3A.pdf", sep=""), width = 8)

# save for external plotting
wb <- createWorkbook()
sheetName <- "relative abdundance"
addWorksheet(wb, sheetName)

writeData(wb, sheetName, df2plot )
addFilter(wb, sheetName, row = 1, cols = 1:ncol(df2plot))
setColWidths(wb, sheetName, cols = 1:ncol(df2plot), widths="auto")

saveWorkbook(wb, file = paste(base, "/Supplement/Figure3A.xlsx", sep=""), overwrite = T) 

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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] patchwork_1.1.2             shiny_1.7.3                 SeuratWrappers_0.3.1        monocle3_1.2.9             
# [5] SingleCellExperiment_1.18.1 SummarizedExperiment_1.26.1 GenomicRanges_1.48.0        GenomeInfoDb_1.32.4        
# [9] MatrixGenerics_1.8.1        matrixStats_0.62.0          openxlsx_4.2.5.1            ggplot2_3.3.6              
# [13] GenBinomApps_1.2            org.Mm.eg.db_3.15.0         AnnotationDbi_1.58.0        IRanges_2.30.1             
# [17] S4Vectors_0.34.0            Biobase_2.56.0              BiocGenerics_0.42.0         clusterProfiler_4.4.4      
# [21] sp_1.5-0                    SeuratObject_4.1.2          Seurat_4.2.0               
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.2             R.utils_2.12.0         reticulate_1.26        lme4_1.1-30            tidyselect_1.2.0      
# [6] RSQLite_2.2.18         htmlwidgets_1.5.4      grid_4.2.1             BiocParallel_1.30.4    Rtsne_0.16            
# [11] scatterpie_0.1.8       munsell_0.5.0          ragg_1.2.4             codetools_0.2-18       ica_1.0-3             
# [16] future_1.28.0          miniUI_0.1.1.1         withr_2.5.0            spatstat.random_2.2-0  colorspace_2.0-3      
# [21] GOSemSim_2.22.0        progressr_0.11.0       rstudioapi_0.14        ROCR_1.0-11            tensor_1.5            
# [26] DOSE_3.22.1            listenv_0.8.0          labeling_0.4.2         GenomeInfoDbData_1.2.8 polyclip_1.10-4       
# [31] bit64_4.0.5            farver_2.1.1           downloader_0.4         parallelly_1.32.1      vctrs_0.5.0           
# [36] treeio_1.20.2          generics_0.1.3         R6_2.5.1               graphlayouts_0.8.3     rsvd_1.0.5            
# [41] bitops_1.0-7           spatstat.utils_3.0-1   cachem_1.0.6           fgsea_1.22.0           gridGraphics_0.5-1    
# [46] DelayedArray_0.22.0    assertthat_0.2.1       promises_1.2.0.1       scales_1.2.1           ggraph_2.1.0          
# [51] enrichplot_1.16.2      rgeos_0.5-9            gtable_0.3.1           globals_0.16.1         goftest_1.2-3         
# [56] tidygraph_1.2.2        rlang_1.0.6            systemfonts_1.0.4      splines_4.2.1          lazyeval_0.2.2        
# [61] spatstat.geom_3.0-3    BiocManager_1.30.19    reshape2_1.4.4         abind_1.4-5            httpuv_1.6.6          
# [66] qvalue_2.28.0          tools_4.2.1            ggplotify_0.1.0        ellipsis_0.3.2         spatstat.core_2.4-4   
# [71] jquerylib_0.1.4        RColorBrewer_1.1-3     proxy_0.4-27           ggridges_0.5.4         Rcpp_1.0.9            
# [76] plyr_1.8.7             zlibbioc_1.42.0        purrr_0.3.5            RCurl_1.98-1.9         rpart_4.1.16          
# [81] deldir_1.0-6           pbapply_1.5-0          viridis_0.6.2          cowplot_1.1.1          zoo_1.8-11            
# [86] ggrepel_0.9.1          cluster_2.1.3          magrittr_2.0.3         data.table_1.14.4      scattermore_0.8       
# [91] DO.db_2.9              lmtest_0.9-40          RANN_2.6.1             packrat_0.8.1          fitdistrplus_1.1-8    
# [96] mime_0.12              xtable_1.8-4           gridExtra_2.3          compiler_4.2.1         tibble_3.1.8          
# [101] KernSmooth_2.23-20     crayon_1.5.2           shadowtext_0.1.2       R.oo_1.25.0            minqa_1.2.5           
# [106] htmltools_0.5.3        ggfun_0.0.7            mgcv_1.8-40            later_1.3.0            tidyr_1.2.1           
# [111] aplot_0.1.8            DBI_1.1.3              tweenr_2.0.2           MASS_7.3-57            boot_1.3-28           
# [116] leidenbase_0.1.12      Matrix_1.5-1           cli_3.4.1              R.methodsS3_1.8.2      parallel_4.2.1        
# [121] igraph_1.3.5           pkgconfig_2.0.3        terra_1.6-17           plotly_4.10.0          spatstat.sparse_3.0-0 
# [126] ggtree_3.4.4           bslib_0.4.0            XVector_0.36.0         yulab.utils_0.0.5      stringr_1.4.1         
# [131] digest_0.6.30          sctransform_0.3.5      RcppAnnoy_0.0.20       spatstat.data_3.0-0    Biostrings_2.64.1     
# [136] leiden_0.4.3           fastmatch_1.1-3        tidytree_0.4.1         uwot_0.1.14            nloptr_2.0.3          
# [141] lifecycle_1.0.3        nlme_3.1-157           jsonlite_1.8.3         limma_3.52.4           viridisLite_0.4.1     
# [146] fansi_1.0.3            pillar_1.8.1           lattice_0.20-45        KEGGREST_1.36.3        fastmap_1.1.0         
# [151] httr_1.4.4             survival_3.3-1         GO.db_3.15.0           remotes_2.4.2          glue_1.6.2            
# [156] zip_2.2.2              png_0.1-7              bit_4.0.4              sass_0.4.2             ggforce_0.4.1         
# [161] stringi_1.7.8          blob_1.2.3             textshaping_0.3.6      memoise_2.0.1          dplyr_1.0.10          
# [166] irlba_2.3.5.1          future.apply_1.9.1     ape_5.6-2  

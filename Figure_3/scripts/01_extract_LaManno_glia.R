# first step of analysis
# extracts a set of glia and glioblast from LaManno for further analyses
# checked 4.10.2023

library (loomR)
library (Seurat)

# define a base folder
base <- "./"

# read location information - defines which cells to keep
# this part is identical to Figure 1 analysis
location2keep <- read.table(paste(base, "other_data/LaManno.Location-Giselle.csv", sep=""), sep="\t", header = T)
location2keep <- location2keep[which(location2keep$use == "yes"),]

# process loom file
lfile <- connect(filename = paste(base, "other_data/dev_all.loom", sep=""), mode = "r+", skip.validate = T)

# get samples indices from Midbrain or Dorsal Midbrain
mbdIDX <- which( lfile$col.attrs$Tissue[] == "Midbrain" | lfile$col.attrs$Tissue[] == "MidbrainDorsal" )
#mbdIDX <- which( grepl( x = lfile$col.attrs$Tissue[], pattern = "midbrain", ignore.case = T ))
mbdIDX <- mbdIDX[which(!(duplicated(lfile$col.attrs$CellID[mbdIDX]) | duplicated(lfile$col.attrs$CellID[mbdIDX], fromLast = T)))]

# additional filtering based on cluster annotation
sclassMBD <- which( grepl( x=lfile$col.attrs$Subclass[], pattern="midbrain", ignore.case = T) | grepl( x=lfile$col.attrs$Subclass[], pattern="Mixed" ) )
mbdIDX <- intersect( mbdIDX, sclassMBD)

# identify which cells from e9-11 should be removed
locIDX <-  which( !lfile$col.attrs$Location[] %in% c(location2keep$location) )
locIDX <- intersect( locIDX, which( lfile$col.attrs$Age[] %in% c("e9.0", "e10.0", "e11.0")))

# focus only on relevant cell types
cti_MBd <- intersect( mbdIDX, which( lfile$col.attrs$Class[] %in% c("Glioblast", "Radial glia")))

# remove cells that are not from the location of interest
cti_MBd <- setdiff(cti_MBd, locIDX)

message("retained ", length(cti_MBd), " cells")
  
#extract the reads from these cells - NOTE row/col are changed according to loom convention!
data.subset <- lfile[["matrix"]][cti_MBd, ]
  
#give cellID as row
rownames(data.subset) <- lfile$col.attrs$CellID[cti_MBd]
#give Gene name as col
colnames(data.subset) <- lfile$row.attrs$Gene[]
  
#create the Seurat object
MBd_Seurat <- CreateSeuratObject(counts = as.matrix(t(data.subset)), min.cells = 5, min.features = 500)
  
#check reads in mitocondrial genome
#really just a check - no QC on this data as already pre-analysed
MBd_Seurat[["percent.mt"]] <- PercentageFeatureSet(MBd_Seurat, pattern = "^mt-")
  
#add meta data
meta_data <- list()
meta_data[["cell_type"]] <- lfile$col.attrs$Class[cti_MBd]
names(meta_data[["cell_type"]]) <- lfile$col.attrs$CellID[cti_MBd]
  
meta_data[["Age"]] <- lfile$col.attrs$Age[cti_MBd]
names(meta_data[["Age"]]) <- lfile$col.attrs$CellID[cti_MBd]

meta_data[["Tissue"]] <- lfile$col.attrs$Tissue[cti_MBd]
names(meta_data[["Tissue"]]) <- lfile$col.attrs$CellID[cti_MBd]  

meta_data[["Location"]] <- lfile$col.attrs$Location_E9_E11[cti_MBd]
names(meta_data[["Location"]]) <- lfile$col.attrs$CellID[cti_MBd]

meta_data[["clusterID"]] <- lfile$col.attrs$ClusterName[cti_MBd]
names(meta_data[["clusterID"]]) <- lfile$col.attrs$CellID[cti_MBd]

meta_data[["subclass"]] <- lfile$col.attrs$Subclass[ cti_MBd ]
names(meta_data[["subclass"]]) <- lfile$col.attrs$CellID[cti_MBd]

#add these labels to the Seurat object
for (i in names(meta_data)){
  MBd_Seurat <- AddMetaData(
    object = MBd_Seurat,
    metadata = meta_data[[i]],
    col.name = i
  )  
}

# write out before integration
saveRDS(object = MBd_Seurat, file = paste(base, "/RDS_files/seurat.LaManno.glia.RDS", sep="") )
 
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
#   [1] sp_1.5-0           SeuratObject_4.1.2 Seurat_4.2.0       loomR_0.2.1.9000   hdf5r_1.3.7        R6_2.5.1          
# 
# loaded via a namespace (and not attached):
# [1] nlme_3.1-157          matrixStats_0.62.0    spatstat.sparse_3.0-0 bit64_4.0.5           RcppAnnoy_0.0.20      RColorBrewer_1.1-3   
# [7] httr_1.4.4            sctransform_0.3.5     tools_4.2.1           utf8_1.2.2            irlba_2.3.5.1         rpart_4.1.16         
# [13] KernSmooth_2.23-20    uwot_0.1.14           mgcv_1.8-40           rgeos_0.5-9           lazyeval_0.2.2        colorspace_2.0-3     
# [19] tidyselect_1.2.0      gridExtra_2.3         bit_4.0.4             compiler_4.2.1        progressr_0.11.0      cli_3.4.1            
# [25] plotly_4.10.0         scales_1.2.1          lmtest_0.9-40         spatstat.data_3.0-0   ggridges_0.5.4        pbapply_1.5-0        
# [31] goftest_1.2-3         stringr_1.4.1         digest_0.6.30         spatstat.utils_3.0-1  pkgconfig_2.0.3       htmltools_0.5.3      
# [37] parallelly_1.32.1     fastmap_1.1.0         htmlwidgets_1.5.4     rlang_1.0.6           rstudioapi_0.14       shiny_1.7.3          
# [43] generics_0.1.3        zoo_1.8-11            jsonlite_1.8.3        spatstat.random_2.2-0 ica_1.0-3             dplyr_1.0.10         
# [49] magrittr_2.0.3        patchwork_1.1.2       Matrix_1.5-1          Rcpp_1.0.9            munsell_0.5.0         fansi_1.0.3          
# [55] abind_1.4-5           reticulate_1.26       lifecycle_1.0.3       stringi_1.7.8         MASS_7.3-57           Rtsne_0.16           
# [61] plyr_1.8.7            grid_4.2.1            parallel_4.2.1        listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1        
# [67] deldir_1.0-6          miniUI_0.1.1.1        lattice_0.20-45       cowplot_1.1.1         splines_4.2.1         tensor_1.5           
# [73] pillar_1.8.1          igraph_1.3.5          spatstat.geom_3.0-3   future.apply_1.9.1    reshape2_1.4.4        codetools_0.2-18     
# [79] leiden_0.4.3          glue_1.6.2            packrat_0.8.1         data.table_1.14.4     png_0.1-7             vctrs_0.5.0          
# [85] httpuv_1.6.6          polyclip_1.10-4       gtable_0.3.1          RANN_2.6.1            purrr_0.3.5           spatstat.core_2.4-4  
# [91] tidyr_1.2.1           scattermore_0.8       future_1.28.0         ggplot2_3.3.6         mime_0.12             xtable_1.8-4         
# [97] later_1.3.0           survival_3.3-1        viridisLite_0.4.1     tibble_3.1.8          cluster_2.1.3         globals_0.16.1       
# [103] fitdistrplus_1.1-8    ellipsis_0.3.2        ROCR_1.0-11      

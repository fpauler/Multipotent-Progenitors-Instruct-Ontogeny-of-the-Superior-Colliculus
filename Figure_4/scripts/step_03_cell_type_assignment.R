# third step - assigning cell types
# checked 6.10.2023

library (Seurat)
library (dplyr)
library (batchelor)
library (ggplot2)
library (openxlsx)
library (umap)
library(scales)

# fastMNN analysis Seurat style - modified code from fast_mnn.R from Seurat

RunFastMNN_FMP <- function(
  object.list,
  assay = NULL,
  features = features,
  reduction.name = "mnn",
  reduction.key = "mnn_",
  reconstructed.assay = "mnnreconstructed",
  verbose = TRUE,
  ...
) {
  
  objects.sce <- lapply(
    X = object.list,
    FUN = function(x, f) {
      # this is my modification - use scale.data
      return(SingleCellExperiment(list(logcounts = x[["SCT"]]@scale.data[f,]), 
                                  colData=DataFrame(cellID=colnames(x[["SCT"]]@scale.data[f,])),
                                  rowData=DataFrame(rownames(x[["SCT"]]@scale.data[f,]))))
    },
    f = features
  )
  
  integrated <- merge(
    x = object.list[[1]],
    y = object.list[2:length(x = object.list)]
  )
  
  out <- do.call(
    what = batchelor::fastMNN,
    args = c(
      objects.sce,
      list(...)
    )
  )
  
  rownames(x = SingleCellExperiment::reducedDim(x = out)) <- colnames(x = integrated)
  colnames(x = SingleCellExperiment::reducedDim(x = out)) <- paste0(reduction.key, 1:ncol(x = SingleCellExperiment::reducedDim(x = out)))
  integrated[[reduction.name]] <- CreateDimReducObject(
    embeddings = SingleCellExperiment::reducedDim(x = out),
    loadings = as.matrix(SingleCellExperiment::rowData(x = out)),
    assay = DefaultAssay(object = integrated),
    key = reduction.key
  )
  
  # Add reconstructed matrix (gene x cell)
  integrated[[reconstructed.assay]] <- CreateAssayObject(
    data = as(object = SummarizedExperiment::assay(x = out), Class = "sparseMatrix"),
  )
  # Add variable features
  VariableFeatures(object = integrated[[reconstructed.assay]]) <- features
  Tool(object = integrated) <- S4Vectors::metadata(x = out)
  integrated <- LogSeuratCommand(object = integrated)
  return(integrated)
}

# closest centroid for bootstrapping 
cor_analysis <- function ( gl_boot_list, ref, patchSeq, cells2clusters  ) {
  
  cl <- makeCluster(5)
  clusterExport(cl, c("gl_boot_list", "ref", "patchSeq", "cells2clusters"), envir=environment())
  
  out <- pbapply::pblapply(cl = cl, X = 1:length(gl_boot_list), FUN = function(x) {
    
    gl_boot <- gl_boot_list[[x]]
    
    sapply ( unique(cells2clusters$cluster), function(cluster) {
      
      cells <- cells2clusters$sample_id[which(cells2clusters$cluster == cluster)]
      
      av_ref <- apply(ref[,cells], 1, mean)
      
      # calculate Pearson correlation between each PatchSeq cell and the average expression of the Zeisel clusters
      
      sapply ( colnames(patchSeq), function (y) {
        tmp <- cor.test(av_ref[gl_boot], patchSeq[gl_boot, y], alternative = "greater")
        return(tmp$estimate)
      })
    })
  }) 
  
  stopCluster(cl)
  
  return (out)
  
}

# define a base folder
base <- "./"

set.seed(2401)

# define Cell Clusters / t-types of Interest - based on Zeisel annotation (http://mousebrain.org/celltypes/)
# all those Tell types predicted to be located in our regions of interest 

cgi <- c("MEGLU1", "MEGLU2", "MEGLU3", "MEGLU4", "MEGLU5", "MEGLU6",
         "MEINH2", "MEINH3", "MEINH5", "MEINH6", "MEINH7", "MEINH8", "MEINH9", "MEINH10", "MEINH11", "MEINH12")

# prepare a reference dataset with only t-types of interest
if (!file.exists(paste(base, "RDS_files/cti_Seurat_RAW.RDS", sep="/"))) {
  
  library (loomR)
  
  lfile <- connect(filename = paste(base, "other_data/l5_all.loom", sep="/"), mode = "r+")
  
  ctiIDX <- which(lfile$col.attrs$ClusterName[] %in% cgi)
  length(ctiIDX)
  # [1] 3698
  
  #remove duplicated CellIDs - only 1
  ctiIDX <- ctiIDX[which(!(duplicated(lfile$col.attrs$CellID[ctiIDX]) | duplicated(lfile$col.attrs$CellID[ctiIDX], fromLast = T)))]
  length(ctiIDX)
  # [1] 3696
  
  #extract the reads from these cells - NOTE row/col are changed according to convention!
  data.subset <- lfile[["matrix"]][ctiIDX, ]
  
  #give cellID as row
  rownames(data.subset) <- lfile$col.attrs$CellID[ctiIDX]
  #give Gene name as col
  colnames(data.subset) <- lfile$row.attrs$Gene[]
  
  #create the Seurat object
  cti_Seurat <- CreateSeuratObject(counts = as.matrix(t(data.subset)), min.cells = 5, min.features = 500)
  
  #add meta data - essentially cluster labels
  taxonomy <- list()
  taxonomy[[1]] <- lfile$col.attrs$TaxonomyRank1[ctiIDX]
  names(taxonomy[[1]]) <- lfile$col.attrs$CellID[ctiIDX]
  
  taxonomy[[2]] <- lfile$col.attrs$TaxonomyRank2[ctiIDX]
  names(taxonomy[[2]]) <- lfile$col.attrs$CellID[ctiIDX]
  
  taxonomy[[3]] <- lfile$col.attrs$TaxonomyRank3[ctiIDX]
  names(taxonomy[[3]]) <- lfile$col.attrs$CellID[ctiIDX]
  
  taxonomy[[4]] <- lfile$col.attrs$TaxonomyRank4[ctiIDX]
  names(taxonomy[[4]]) <- lfile$col.attrs$CellID[ctiIDX]
  
  taxonomy[[5]] <- lfile$col.attrs$ClusterName[ctiIDX]
  names(taxonomy[[5]]) <- lfile$col.attrs$CellID[ctiIDX]
  
  index <- lfile$col.attrs$SampleID[ctiIDX]
  names(index) <- lfile$col.attrs$CellID[ctiIDX]
  
  #add these labels to the Seurat object
  for (i in 1:5){
    cti_Seurat <- AddMetaData(
      object = cti_Seurat,
      metadata = taxonomy[[i]],
      col.name = paste('tax', i, sep='_')
    )  
  }
  
  #add index to object
  cti_Seurat <- AddMetaData(
    object = cti_Seurat,
    metadata = index,
    col.name = "idx"
  )  
  
  saveRDS(object = cti_Seurat, file = paste(base, "RDS_files/cti_Seurat_RAW.RDS", sep="/")  )
  
}
# stop and restart point - if you have limited memory

if (!file.exists(paste(base, "RDS_files/cti_Seurat_SCT.RDS", sep="/"))) {
  
  # read in saved object
  cti_Seurat <- readRDS( file = paste(base, "RDS_files/cti_Seurat_RAW.RDS", sep="/") )
  
  # SCT workflow
  cti_Seurat_standard_sct <- SCTransform(cti_Seurat)
  cti_Seurat_standard_sct <- RunPCA(cti_Seurat_standard_sct, assay = "SCT")
  cti_Seurat_standard_sct <- RunUMAP(cti_Seurat_standard_sct, dims = 1:20)
  cti_Seurat_standard_sct_plot <- DimPlot(cti_Seurat_standard_sct, group.by = "tax_5", label = T)
  
  # save this plot for later
  saveRDS(object = cti_Seurat_standard_sct_plot, file = paste(base, "RDS_files/cti_Seurat_standard_sct_plot.RDS", sep="/")  )
  
  # do analysis with all genes for inegration later
  all.genes <- rownames(cti_Seurat)

  cti_Seurat <- SCTransform(cti_Seurat, residual.features = all.genes)
  
  cti_Seurat <- RunPCA(cti_Seurat, assay = "SCT")
  cti_Seurat <- RunUMAP(cti_Seurat, dims = 1:25)
  
  saveRDS(object = cti_Seurat, file = paste(base, "RDS_files/cti_Seurat_SCT.RDS", sep="/")  )
  
  # prepare marker genes for the different cell types
  Idents(cti_Seurat) <- "tax_5"
  all.markers <- FindAllMarkers(cti_Seurat, logfc.threshold = 0.2)
  
  saveRDS(object = all.markers, file = paste(base, "RDS_files/all.markers_tax5.RDS", sep="/")  )
  
} else {
  cti_Seurat<- readRDS( file = paste(base, "RDS_files/cti_Seurat_SCT.RDS", sep="/")  )
  all.markers <- readRDS( file = paste(base, "RDS_files/all.markers_tax5.RDS", sep="/")  )
  cti_Seurat_standard_sct_plot <- readRDS( file = paste(base, "RDS_files/cti_Seurat_standard_sct_plot.RDS", sep="/")  )  
}

PatchSeq_Seurat_SCT <- readRDS(file = paste(base, "RDS_files/PatchSeq_Seurat_TPM_SCT.RDS", sep="/") )
orig.varFeatures.PatchSeq <- VariableFeatures(PatchSeq_Seurat_SCT)  

all.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top_markers

varFeatures <- unique(top_markers$gene)

# filter varFeatures for genes present in both Seurat objects
varFeatures <- varFeatures[which(varFeatures %in% rownames(cti_Seurat))]
varFeatures <- varFeatures[which(varFeatures %in% rownames(PatchSeq_Seurat_SCT))]

VariableFeatures(cti_Seurat) <- varFeatures
VariableFeatures(PatchSeq_Seurat_SCT) <- varFeatures

################################################
##
# closest/nearest centroid with fastMNN data integration 
##
################################################

# fastMNN batch correction
# get the overlap between the variabe genes from the 2 datasets

# NOTE: more genes performed better in my hands
all.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top_markers

varFeatures <- unique(top_markers$gene)
varFeatures <- varFeatures[which(varFeatures %in% rownames(cti_Seurat[["SCT"]]@scale.data))]
varFeatures <- varFeatures[which(varFeatures %in% rownames(PatchSeq_Seurat_SCT[["SCT"]]@scale.data))]

# get the number of overlapping genes
# these genes will be used in the folowing analyses
length (varFeatures)
# 996

# here the fastMNN integration is done
seurat_integrated <- RunFastMNN_FMP( object.list = c(cti_Seurat, PatchSeq_Seurat_SCT), features = varFeatures )

# add the PatchSeq cells onto the clean UMAP from Zeisel data - inspired by code found here:
# https://github.com/AllenInstitute/patchseq_human_L23/blob/master/scripts/Code4_mouse_patchseq_visualization_L23exc.Rmd
data_patch <- FetchData(seurat_integrated, vars=paste0("mnn_",1:dims), cells=colnames(PatchSeq_Seurat_SCT), slot="scale.data") 
data_FACS  <- FetchData(seurat_integrated, vars=paste0("mnn_",1:dims), cells=colnames(cti_Seurat), slot="scale.data")

umap_FACs  <- umap(data_FACS)

umap_FACs$layout <- as.matrix(cti_Seurat_standard_sct_plot$data[,c("UMAP_1", "UMAP_2")])

# Predict Patch-seq cell locations using "predict" function
umap_patch <- predict(umap_FACs,data_patch)

seurat_integrated <- RunUMAP( seurat_integrated, reduction = "mnn", dims = 1:20 )

# Update the umap coordinates from above.
seurat_integrated@reductions$umap@cell.embeddings <- rbind(umap_FACs$layout,umap_patch)

seurat_integrated$Genotype[which(is.na(seurat_integrated$Genotype))] <- "reference"

# show that Genotype does not majorly influence distribution
UMAP_integrated <- DimPlot(seurat_integrated, group.by = "Genotype")

ggplot(UMAP_integrated$data, aes(x=UMAP_1, y=UMAP_2, color=Genotype, size = Genotype)) + geom_point() +
  scale_color_manual(values = c("reference" = "grey60", "M11GTTG-Fz10-CreER" = "red", "M11GTTG-Sox2-CreER" = "green")) +
  scale_size_manual(values = c("reference" = 1, "M11GTTG-Fz10-CreER" = 2.5, "M11GTTG-Sox2-CreER" = 2.5)) +
  theme_classic()
ggsave(filename = paste(base, "QC/UMAP_integrated_genotypes.pdf", sep="/"), width=10, height=8)

# prepare a dataframe to plot some meta-data from the reference
plot_tax4 <- DimPlot(seurat_integrated, group.by = "tax_4")
plot_tax5 <- DimPlot(seurat_integrated, group.by = "tax_5")
plot_gen <- DimPlot(seurat_integrated, group.by = "Genotype")
# cell order is identical in both plots
identical(rownames(plot_tax4$data), rownames(plot_tax5$data))

df2plot <- cbind(plot_tax4$data, plot_tax5$data[,"tax_5"], plot_gen$data[,"Genotype"])
colnames(df2plot) <- c("UMAP_1", "UMAP_2", "tax4", "tax5", "Genotype")
# color palette for UMAP
hex_codes2 <- hue_pal()(length(cgi)) 
names(hex_codes2) <- cgi
hex_codes2 <- c(hex_codes2, c("PatchSeq" = "black"))

df2plot$tax5 <- as.character(df2plot$tax5)
df2plot$tax5[which(is.na(df2plot$tax5))] <- "PatchSeq"

# plot a UMAP with Patch-Seq cells - not final!
ggplot(df2plot, aes(x=UMAP_1, y=UMAP_2, color=tax5, size = Genotype)) + geom_point() +
  scale_color_manual(values = hex_codes2) +
  scale_size_manual(values = c("reference" = 1, "M11GTTG-Fz10-CreER" = 2.5, "M11GTTG-Sox2-CreER" = 2.5)) +
  theme_classic()

ggsave(filename = paste(base, "QC/UMAP_integrated.pdf", sep="/"), width=10, height=8)

# save the meta data for the plot
saveRDS(object = df2plot, file = paste(base, "RDS_files/integrated_umap_df.RDS", sep="/"))

#extract batch corrected expression values
ref <- as.data.frame( GetAssayData(seurat_integrated, assay = "mnnreconstructed")[varFeatures, rownames(cti_Seurat@meta.data)] )
patchSeq <-  as.data.frame( GetAssayData(seurat_integrated, assay = "mnnreconstructed")[varFeatures, rownames(PatchSeq_Seurat_SCT@meta.data)] )

cells2cluster <- cti_Seurat@meta.data
cells2cluster$sample_id <- rownames(cells2cluster)
cells2cluster <- cells2cluster[, c("sample_id", "tax_5")]

colnames(cells2cluster) <- c("sample_id", "cluster")

out <- sapply( unique(cells2cluster$cluster), FUN = function( cluster ) {
  
  cells <- cells2cluster$sample_id[which(cells2cluster$cluster == cluster)]
  
  av_ref <- apply(ref[,cells], 1, mean)
  
  # calculate Pearson correlation between each PatchSeq cell and the average expression of the Zeisel clusters
  
  sapply ( colnames(patchSeq), function (x) {
    tmp <- cor.test(av_ref[varFeatures], patchSeq[varFeatures,x], alternative = "two.sided", method = "pearson")
    return(paste(tmp$estimate, tmp$p.value, sep="_"))
  })
})


# create a list of additional infos about corellations
assign_cc_cor_list <- list()

# do the cell type assignment
assign_cc_cor_list[[length(assign_cc_cor_list)+1]] <- apply ( out, 1, function(x) {
  correlation <- sapply(x, function (y) as.numeric(strsplit(y, "_")[[1]][1]))
  names(which(correlation == max(correlation)))
})

# record the top correlation
assign_cc_cor_list[[length(assign_cc_cor_list)+1]] <- apply ( out, 1, function(x) {
  correlation <- sapply(x, function (y) as.numeric(strsplit(y, "_")[[1]][1]))
  tmp <- sort(correlation, decreasing = T)
  tmp[1]
})

# significance of correlation
assign_cc_cor_list[[length(assign_cc_cor_list)+1]] <- apply ( out, 1, function(x) {
  correlation <- sapply(x, function (y) as.numeric(strsplit(y, "_")[[1]][1]))
  significance <- sapply(x, function (y) as.numeric(strsplit(y, "_")[[1]][2]))
  tmp <- names(sort(correlation, decreasing = T))
  significance[tmp[1]]
})

# calculate the difference of the correlation of the best and the next-best matching cell type
assign_cc_cor_list[[length(assign_cc_cor_list)+1]] <- apply ( out, 1, function(x) {
  correlation <- sapply(x, function (y) as.numeric(strsplit(y, "_")[[1]][1]))
  tmp <- sort(correlation, decreasing = T)
  tmp[1] - tmp[2]
})


# finally botstrapping to test the influence of the gene list on the final assignment
# save the gene list for reproducibility
if (!file.exists( paste(base, "RDS_files/gl_boot_list.RDS", sep="/") )) {
  gl_boot_list <- lapply(1:100, function (x) return(sample(x = varFeatures, size = length(varFeatures), replace = T)))
  saveRDS(object = gl_boot_list, file = paste(base, "RDS_files/gl_boot_list.RDS", sep="/") )
} else {
  gl_boot_list <- readRDS(file = paste(base, "RDS_files/gl_boot_list.RDS", sep="/")  )
}

cor_boot <- cor_analysis ( gl_boot_list = gl_boot_list, ref = ref, patchSeq = patchSeq, cells2clusters = cells2cluster )

# assign cell types from bootstrapping

# first get the cell type with highest correlation from each bootstrap
assign_mat_boot <- sapply (1:100, function (i){
  apply (cor_boot[[i]], 1, function(x) {
    names(which(x == max(x)))
  })
})

# final cell type assignment - most often assigned label
# NOTE: this can include 2 labels per cell as it is only the frequency of occurence
cell_type_assign_boot <- sapply ( rownames(assign_mat_boot), function (x) {
  tmp <- table( assign_mat_boot[x, ] ) 
  names(tmp)[which(tmp == max(tmp))]
})

cell_type_assign_boot <- sapply(cell_type_assign_boot, function (x) paste(x, collapse=","))
names(cell_type_assign_boot) <- gsub(pattern = ".cor", replacement = "", x = names(cell_type_assign_boot))

# calculate the overlap between the bootstrap and the actual cell type assignment
length(which(cell_type_assign_boot == assign_cc_cor_list[[1]][names(cell_type_assign_boot)])) / length(assign_cc_cor_list[[1]])
# [1] 0.9615385

# score shows how often the final ttype was identified using bootstraped gene lists
assign_cc_cor_list[[length(assign_cc_cor_list)+1]] <- sapply ( rownames(assign_mat_boot), function (x) {
  tmp <- table( assign_mat_boot[x, ] ) 
  final_ttype <- as.character(assign_cc_cor_list[[1]][gsub(pattern = ".cor", replacement = "", x = x)])
  return(tmp[final_ttype] / sum(tmp))
})

# make a named list
names(assign_cc_cor_list) <- c("predicted_ttype", "top_cor", "top_sig", "delta_cor", "boot_score")

# prepare a meta data table to add to Seurat object
assign_cc_cor_meta <- data.frame()
for (i in 1:length(assign_cc_cor_list)) {
  names(assign_cc_cor_list[[i]]) <- sapply(names(assign_cc_cor_list[[i]]), function (x) strsplit(x = x, split = "\\.")[[1]][1])
  tmp_df <- data.frame( cell_id = names(assign_cc_cor_list[[i]]), value = assign_cc_cor_list[[i]])
  colnames(tmp_df) <- c("cell_id", names(assign_cc_cor_list)[i])
  if (i == 1) {
    assign_cc_cor_meta <- tmp_df
  } else {
    assign_cc_cor_meta <- merge(assign_cc_cor_meta, tmp_df, by="cell_id")
  }
}

# add the score to the Seurat object
rownames(assign_cc_cor_meta) <- assign_cc_cor_meta$cell_id
assign_cc_cor_meta <- assign_cc_cor_meta[,2:ncol(assign_cc_cor_meta)]

PatchSeq_Seurat_SCT <- AddMetaData( object = PatchSeq_Seurat_SCT,
                                    metadata = assign_cc_cor_meta )

saveRDS(object = PatchSeq_Seurat_SCT, file = paste(base, "RDS_files/PatchSeq_Seurat_SCT_cellTypeAssignment.RDS", sep="/") )

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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scales_1.1.1                umap_0.2.7.0                openxlsx_4.2.4              ggplot2_3.3.5               batchelor_1.6.3            
# [6] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0 Biobase_2.50.0              GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
# [11] IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.1         MatrixGenerics_1.2.1        matrixStats_0.60.1         
# [16] dplyr_1.0.7                 SeuratObject_4.0.2          Seurat_4.0.4.9000          
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15                colorspace_2.0-2          deldir_0.2-10             ellipsis_0.3.2            ggridges_0.5.3           
# [6] scuttle_1.0.4             XVector_0.30.0            BiocNeighbors_1.8.2       rstudioapi_0.13           spatstat.data_2.1-0      
# [11] farver_2.1.0              leiden_0.3.9              listenv_0.8.0             ggrepel_0.9.1             RSpectra_0.16-0          
# [16] fansi_0.5.0               sparseMatrixStats_1.2.1   codetools_0.2-16          splines_4.0.3             polyclip_1.10-0          
# [21] jsonlite_1.7.2            ResidualMatrix_1.0.0      packrat_0.7.0             ica_1.0-2                 cluster_2.1.0            
# [26] png_0.1-7                 uwot_0.1.10               shiny_1.6.0               sctransform_0.3.2         spatstat.sparse_2.0-0    
# [31] compiler_4.0.3            httr_1.4.2                assertthat_0.2.1          Matrix_1.3-4              fastmap_1.1.0            
# [36] lazyeval_0.2.2            cli_3.0.1                 BiocSingular_1.6.0        later_1.3.0               htmltools_0.5.2          
# [41] tools_4.0.3               rsvd_1.0.5                igraph_1.2.6              gtable_0.3.0              glue_1.4.2               
# [46] GenomeInfoDbData_1.2.4    RANN_2.6.1                reshape2_1.4.4            Rcpp_1.0.7                scattermore_0.7          
# [51] vctrs_0.3.8               nlme_3.1-149              DelayedMatrixStats_1.12.3 lmtest_0.9-38             stringr_1.4.0            
# [56] globals_0.14.0            beachmat_2.6.4            mime_0.11                 miniUI_0.1.1.1            lifecycle_1.0.0          
# [61] irlba_2.3.3               goftest_1.2-2             future_1.22.1             MASS_7.3-53               zlibbioc_1.36.0          
# [66] zoo_1.8-9                 spatstat.core_2.3-0       promises_1.2.0.1          spatstat.utils_2.2-0      RColorBrewer_1.1-2       
# [71] reticulate_1.20           pbapply_1.5-0             gridExtra_2.3             rpart_4.1-15              stringi_1.7.4            
# [76] zip_2.2.0                 BiocParallel_1.24.1       rlang_0.4.11              pkgconfig_2.0.3           bitops_1.0-7             
# [81] lattice_0.20-41           ROCR_1.0-11               purrr_0.3.4               tensor_1.5                labeling_0.4.2           
# [86] patchwork_1.1.1           htmlwidgets_1.5.4         cowplot_1.1.1             tidyselect_1.1.1          parallelly_1.28.1        
# [91] RcppAnnoy_0.0.19          plyr_1.8.6                magrittr_2.0.1            R6_2.5.1                  generics_0.1.0           
# [96] DelayedArray_0.16.3       DBI_1.1.1                 withr_2.4.2               pillar_1.6.2              mgcv_1.8-33              
# [101] fitdistrplus_1.1-5        survival_3.2-7            abind_1.4-5               RCurl_1.98-1.4            tibble_3.1.4             
# [106] future.apply_1.8.1        crayon_1.4.1              KernSmooth_2.23-17        utf8_1.2.2                spatstat.geom_2.2-2      
# [111] plotly_4.9.4.1            grid_4.0.3                data.table_1.14.0         digest_0.6.27             xtable_1.8-4             
# [116] tidyr_1.1.3               httpuv_1.6.3              openssl_1.4.5             munsell_0.5.0             viridisLite_0.4.0        
# [121] askpass_1.1  

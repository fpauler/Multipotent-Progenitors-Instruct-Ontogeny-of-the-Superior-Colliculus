# second step of analysis
# prepares reference and PatchSeq cells

# checked 6.10.2023

library (Seurat)
library (openxlsx)
library (dplyr)
library (ggplot2)

# define a base folder
base <- "./"

set.seed(2401)

#prepare environment hash table as lookup for ENSMUSG - Symbol

if (!file.exists(paste(base, "RDS_files/ensmusg_symbol.RDS", sep="/"))) {
  
  ensmusg_symbol <- new.env()
  
  #read conversion table
  ensmusg_symbol_chr <- read.table(paste(base, "annotation/exchange_table.tsv", sep="/"), fill = T, header = F, stringsAsFactors = F)
  colnames(ensmusg_symbol_chr) <- c("ENSMUSG", "SYMBOL", "chr")
  
  # sanity check for unique ENSMUSG IDs
  table(table(ensmusg_symbol_chr$ENSMUSG))
  
  # 1 
  # 55359 
  
  #fill conversion table - handled just like a list but apparently much faster
  for (x in 1:nrow(ensmusg_symbol_chr)) {
    #deal here with multi mappers - one SYMBOL mapped to many ENSMUSG
    #a comma in the SYMBOL name indicates these multi-mappers
    ensmusg_symbol[[ensmusg_symbol_chr$ENSMUSG[x]]] <- ensmusg_symbol_chr$SYMBOL[x]
  }
  
  saveRDS(ensmusg_symbol, paste(base, "RDS_files/ensmusg_symbol.RDS", sep="/")) 
} 

# NOTE: not run in the same environment as other parts
# this part was done on the IST cluster - needs extensive memory!

if ( !file.exists(  paste(base, "RDS_files/MBd_Seurat_wholeMBD_RAW.RDS", sep="/")) ) {
  library (loomR)
  

  lfile <- connect(filename = paste(base, "other_data/l5_all.loom", sep="/"), mode = "r+")

  # get samples indices from cortex
  mbdIDX <- which(grepl(pattern = 'MBd', x = lfile$col.attrs$Tissue[]))
  mbdIDX <- mbdIDX[which(!(duplicated(lfile$col.attrs$CellID[mbdIDX]) | duplicated(lfile$col.attrs$CellID[mbdIDX], fromLast = T)))]
  
  age <- unique(lfile$col.attrs$Age[mbdIDX])
  # [1] "p20"    "p25-27" "p22"
  
  #extract the reads from these cells - NOTE row/col are changed according to convention!
  data.subset <- lfile[["matrix"]][mbdIDX, ]
  
  #give cellID as row
  rownames(data.subset) <- lfile$col.attrs$CellID[mbdIDX]
  #give Gene name as col
  colnames(data.subset) <- lfile$row.attrs$Gene[]
  
  #create the Seurat object
  MBd_Seurat <- CreateSeuratObject(counts = as.matrix(t(data.subset)), min.cells = 5, min.features = 500)
  
  #check reads in mitocondrial genome
  #really just a check - no QC on this data as already pre-analysed
  MBd_Seurat[["percent.mt"]] <- PercentageFeatureSet(MBd_Seurat, pattern = "^mt-")
  
  #add meta data - essentially cluster labels
  taxonomy <- list()
  taxonomy[[1]] <- lfile$col.attrs$TaxonomyRank1[mbdIDX]
  names(taxonomy[[1]]) <- lfile$col.attrs$CellID[mbdIDX]
  
  taxonomy[[2]] <- lfile$col.attrs$TaxonomyRank2[mbdIDX]
  names(taxonomy[[2]]) <- lfile$col.attrs$CellID[mbdIDX]
  
  taxonomy[[3]] <- lfile$col.attrs$TaxonomyRank3[mbdIDX]
  names(taxonomy[[3]]) <- lfile$col.attrs$CellID[mbdIDX]
  
  taxonomy[[4]] <- lfile$col.attrs$TaxonomyRank4[mbdIDX]
  names(taxonomy[[4]]) <- lfile$col.attrs$CellID[mbdIDX]
  
  taxonomy[[5]] <- lfile$col.attrs$ClusterName[mbdIDX]
  names(taxonomy[[5]]) <- lfile$col.attrs$CellID[mbdIDX]
  
  index <- lfile$col.attrs$SampleID[mbdIDX]
  names(index) <- lfile$col.attrs$CellID[mbdIDX]
  
  #add these labels to the Seurat object
  for (i in 1:5){
    MBd_Seurat <- AddMetaData(
      object = MBd_Seurat,
      metadata = taxonomy[[i]],
      col.name = paste('tax', i, sep='_')
    )  
  }
  
  #add index to object
  MBd_Seurat <- AddMetaData(
    object = MBd_Seurat,
    metadata = index,
    col.name = "idx"
  )  

  saveRDS(object = MBd_Seurat, file = paste(base, "RDS_files/MBd_Seurat_wholeMBD_RAW.RDS", sep="/" ))
}

if(!file.exists(paste(base, "RDS_files/MBd_Seurat_wholeMBD_final.RDS", sep="/"))) {
  
  MBd_Seurat <- readRDS( file = paste(base, "RDS_files/MBd_Seurat_wholeMBD_RAW.RDS", sep="/"))

  MBd_Seurat <- NormalizeData(MBd_Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  all.genes <- rownames(MBd_Seurat)
  MBd_Seurat <- ScaleData(MBd_Seurat, features = all.genes)
  
  # quick sanity check
  MBd_Seurat <- FindVariableFeatures(MBd_Seurat)
  MBd_Seurat <- RunPCA(MBd_Seurat, npcs = 30, verbose = FALSE)
  MBd_Seurat <- RunUMAP(MBd_Seurat, reduction = "pca", dims = 1:30)
  
  # UMAP  
  DimPlot( object = MBd_Seurat, group.by = "tax_4", label = T) + NoLegend()
  ggsave(paste(base, "QC/Zeisel_MBd_UMAP.pdf", sep="/"), width = 10, height=10)
  
  # plot to show exclusivity of excitatory and inhibitory markers
  FeaturePlot(object = MBd_Seurat, features = c("Gad1", "Slc17a6"), order = T, blend = T, 
              min.cutoff = "q10", cols = c("lightgrey", "#ff0000", "blue"), blend.threshold = 0.3)
  ggsave(paste(base, "QC/Zeisel_MBd_UMAP_featurePlot.pdf", sep="/"), width = 25, height=8)
  
  # some infos on this Seurat object:
  # An object of class Seurat 
  # 17672 features across 15071 samples within 1 assay 
  # Active assay: RNA (17672 features, 2000 variable features)
  # 2 dimensional reductions calculated: pca, umap
  
  # save the object for further analysis
  saveRDS(object = MBd_Seurat, file = paste(base, "RDS_files/MBd_Seurat_wholeMBD_final.RDS", sep="/"))
}

# prepare the PatchSeq dataset
if (!file.exists(paste(base, "RDS_files/comb_matrix_TPM_exon_plus_intron_Seurat.RDS", sep="/"))) {
  
    #read gene ID conversion table 
    ensmusg_symbol <- readRDS(file = paste(base, "RDS_files/ensmusg_symbol.RDS", sep="/"))
    
    # to do: read from csv file
    counts.matrix <- read.csv(file = paste(base, "RDS_files/comb_matrix_TPM_exon_plus_intron.csv", sep="/"), header = T)
    
    # process count matrix for further analyses
    rownames(counts.matrix) <- counts.matrix$X
    counts.matrix <- counts.matrix[,2:ncol(counts.matrix)]
    colnames(counts.matrix) <- gsub(pattern = "X", replacement = "", x = colnames(counts.matrix))
    
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
    
    # read the sample list with added STAR alignment stats
    # this was an initial meta data list 
    sample_list <- read.xlsx( xlsxFile = paste(base, "Supplement/PatchSeq_sample_list_final.xlsx", sep="/"), sheet=1)
 
    # mget some stats from presumed negative controls - gives an idea how to set cutoffs
    median(sample_list[which(sample_list$Clone %in% c("Negative control")), "STAR.unique.aligned"])
    # [1] 26299.5
    median(sample_list[which(sample_list$Clone %in% c("Negative control")), "STAR.perc"])
    # [1] 0.02323241
   
    nrow(sample_list)
    # [1] 429
    
    # remove negative controls
    sample_list <- sample_list[which(grepl(pattern = "PS", x = sample_list$Clone)),]
    nrow(sample_list)
    # [1] 399
    
    # remove cells with really low read coverage
    sample_list <- sample_list[ which(sample_list$STAR.perc > 0.05 & sample_list$STAR.unique.aligned > 26000),]
    nrow(sample_list)
    # [1] 352
    
    PatchSeq_Seurat <- CreateSeuratObject(counts = counts.matrix[,as.character(sample_list$sample_id)], min.cells = 5, min.features = 500)
    
    # normalization already done by TPM!
    
    # do scaling for all genes
    all.genes <- rownames(PatchSeq_Seurat)
    
    PatchSeq_Seurat <- ScaleData(PatchSeq_Seurat, features = all.genes)
    
    rownames(sample_list) <- sample_list$sample_id
    PatchSeq_Seurat <- AddMetaData(object = PatchSeq_Seurat, metadata = sample_list)

    # An object of class Seurat 
    # 28817 features across 352 samples within 1 assay 
    # Active assay: RNA (28817 features, 0 variable features)
    
    saveRDS(object = PatchSeq_Seurat, file = paste(base, "RDS_files/comb_matrix_TPM_exon_plus_intron_Seurat.RDS", sep="/"))
} 
# stop and restart point - if you have limited memory

# identify marker genes for broad classes
if (!file.exists( paste(base, "RDS_files/MBd_Seurat_wholeMBD_RAW_markers.RDS", sep="/") )) {
  
  MBd_Seurat <- readRDS(file = paste(base, "RDS_files/MBd_Seurat_wholeMBD_final.RDS", sep="/"))
  Idents(MBd_Seurat) <- "tax_3"
  table( MBd_Seurat@meta.data$tax_3)
  
  # Astroependymal cells            Cerebellum neurons Di- and mesencephalon neurons                  Immune cells 
  #                 3682                             1                          4253                           939 
  #
  # Neural crest-like glia              Oligodendrocytes                Vascular cells 
  #                     71                          4366                          1759 
  
  #remove some minority cell types
  MBd_Seurat <- subset ( MBd_Seurat, idents = c("Neural crest-like glia", "Cerebellum neurons"), invert = T)
  
  # prepare marker genes for the different cell types
  all.markers <- FindAllMarkers(MBd_Seurat, only.pos = T)
  saveRDS(object = all.markers, file = paste(base, "RDS_files/MBd_Seurat_wholeMBD_RAW_markers.RDS", sep="/"))
} else {
  all.markers <-  readRDS( file = paste(base, "RDS_files/MBd_Seurat_wholeMBD_RAW_markers.RDS", sep="/") )
}
# stop and restart point - if you have limited memory

# determine Normalized Marker Sum (NMS) as used for example here: https://elifesciences.org/articles/65482

if ( !file.exists( paste(base, "RDS_files/comb_matrix_TPM_exon_plus_intron_Seurat_NMS.RDS", sep="/") ) ) {
  
  # log normalized data
  MBd_Seurat <- readRDS(file = paste(base, "RDS_files/MBd_Seurat_wholeMBD_final.RDS", sep="/"))
  PatchSeq_Seurat <- readRDS(file = paste(base, "RDS_files/comb_matrix_TPM_exon_plus_intron_Seurat.RDS", sep="/"))

  all.markers %>%
    group_by(cluster) %>%
    top_n(n = 200, wt = avg_log2FC) -> top_markers
  
  cell_types <- sapply(unique(MBd_Seurat@meta.data$tax_3), function (cell_type) {
    length(rownames(MBd_Seurat@meta.data[which(MBd_Seurat@meta.data$tax_3 == cell_type),]))
  })
  
  #remove minor cell types 
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
  #[1] 8.099442
  
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
  
  # add NMS info to Seurat object
  # add more meta data to the Seurat object
  for (col_name in colnames(NMS_df)[2:ncol(NMS_df)]) {
    meta_to_add <- NMS_df[,col_name]
    names(meta_to_add) <- NMS_df$sample_id
    
    PatchSeq_Seurat <- AddMetaData(
      object = PatchSeq_Seurat,
      metadata = meta_to_add,
      col.name = col_name
    )  
  }
  
  saveRDS(object = PatchSeq_Seurat, file = paste(base, "RDS_files/comb_matrix_TPM_exon_plus_intron_Seurat_NMS.RDS", sep="/"))
  
}

# once the preparative work is done
# stop and restart point - if you have limited memory

PatchSeq_Seurat <- readRDS(file = paste(base, "RDS_files/comb_matrix_TPM_exon_plus_intron_Seurat_NMS.RDS", sep="/"))
nrow(PatchSeq_Seurat@meta.data)
# 352
    
# identify the maximum NMS score from non-neuronal (aka contaminating) t-types
NMS_cont <- apply( PatchSeq_Seurat@meta.data[, c("NMS_Oligo", "NMS_Vascular", "NMS_Immune", "NMS_Astro")], 1, max)
    
# plot contaminating vs neuron NMS scores - there are only few outliers
df2plot <- data.frame ( NMS_neuron = PatchSeq_Seurat@meta.data$NMS_Neurons, NMS_cont = as.numeric(NMS_cont))
ggplot( df2plot, aes ( x= NMS_neuron, y = NMS_cont)) + geom_point()
ggsave(filename = paste(base, "QC/NMS_score_plot.pdf", sep="/"))
    
# set contaminating NMS score cutoff
cont_cutoff <- 4
    
# typically the cutoff is set to 0.4, which assumes equal expression levels in ref and PatchSeq
# as we have 8x more expression in PatchSeq, the cutoff would be 3.2 - set to 3 
neuro_cutoff <- 3
    
PatchSeq_Seurat <- subset(PatchSeq_Seurat, NMS_Neurons > neuro_cutoff)
nrow(PatchSeq_Seurat@meta.data)
# 320
    
PatchSeq_Seurat <- subset(PatchSeq_Seurat, NMS_Oligo < cont_cutoff & NMS_Vascular < cont_cutoff & NMS_Immune < cont_cutoff & NMS_Astro < cont_cutoff)
nrow(PatchSeq_Seurat@meta.data)
# 312
    
# do batch correction - some notes on this:
# for the gene expression calculation of PatchSeq I use exon TPM + intron TPM as suggested
# I correct the PatchSeq values based on the % aligned reads as hinted at here: 
# Front. Mol. Neurosci., 08 October 2018 | https://doi.org/10.3389/fnmol.2018.00363
    
# using TPM for SCTransform seems to be OK as discussed here: https://github.com/satijalab/seurat/issues/2854
# caveats of using integrated gene expression data seems to be mainly concerning DEG calculations
    
PatchSeq_Seurat_SCT <- SCTransform(object = PatchSeq_Seurat, vars.to.regress = c("Plate", "STAR.perc"), 
                                   residual.features = rownames(PatchSeq_Seurat) )
    
PatchSeq_Seurat_SCT <- RunPCA(PatchSeq_Seurat_SCT)
PatchSeq_Seurat_SCT <- RunUMAP(PatchSeq_Seurat_SCT, dims = 1:12)
    
PatchSeq_Seurat_SCT <- FindNeighbors(PatchSeq_Seurat_SCT, dims = 1:9, verbose = FALSE)
PatchSeq_Seurat_SCT <- FindClusters(PatchSeq_Seurat_SCT, resolution = 0.6, verbose = FALSE)
    
# check if clustering PatchSeq cells on their own separates excitatory and inhibitory neurons - not really
FeaturePlot(PatchSeq_Seurat_SCT, features = c("Gad1", "Slc17a6"), order = T, pt.size = 1.5, blend = T)
ggsave( file = paste(base, "QC/PatchSeq_Seurat_TPM_SCT_UMAP_expr.pdf", sep="/"), width = 15 )
DimPlot(PatchSeq_Seurat_SCT, group.by = "Genotype", pt.size = 1.5)
ggsave( file = paste(base, "QC/PatchSeq_Seurat_TPM_SCT_UMAP_genotype.pdf", sep="/") )
    
# save for further analyses
saveRDS( object = PatchSeq_Seurat_SCT, file = paste(base, "RDS_files/PatchSeq_Seurat_TPM_SCT.RDS", sep="/"))

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
    #   [1] stats     graphics  grDevices utils     datasets  methods   base     
    # 
    # other attached packages:
    #   [1] ggplot2_3.3.5      dplyr_1.0.7        openxlsx_4.2.4     SeuratObject_4.0.2 Seurat_4.0.4.9000 
    # 
    # loaded via a namespace (and not attached):
    #   [1] nlme_3.1-149          matrixStats_0.60.1    spatstat.sparse_2.0-0 RcppAnnoy_0.0.19      RColorBrewer_1.1-2    httr_1.4.2           
    # [7] sctransform_0.3.2     tools_4.0.3           utf8_1.2.2            R6_2.5.1              irlba_2.3.3           rpart_4.1-15         
    # [13] KernSmooth_2.23-17    uwot_0.1.10           mgcv_1.8-33           DBI_1.1.1             lazyeval_0.2.2        colorspace_2.0-2     
    # [19] withr_2.4.2           tidyselect_1.1.1      gridExtra_2.3         compiler_4.0.3        plotly_4.9.4.1        labeling_0.4.2       
    # [25] scales_1.1.1          lmtest_0.9-38         spatstat.data_2.1-0   ggridges_0.5.3        pbapply_1.5-0         goftest_1.2-2        
    # [31] stringr_1.4.0         digest_0.6.27         spatstat.utils_2.2-0  pkgconfig_2.0.3       htmltools_0.5.2       parallelly_1.28.1    
    # [37] fastmap_1.1.0         htmlwidgets_1.5.4     rlang_0.4.11          rstudioapi_0.13       shiny_1.6.0           farver_2.1.0         
    # [43] generics_0.1.0        zoo_1.8-9             jsonlite_1.7.2        ica_1.0-2             zip_2.2.0             magrittr_2.0.1       
    # [49] patchwork_1.1.1       Matrix_1.3-4          Rcpp_1.0.7            munsell_0.5.0         fansi_0.5.0           abind_1.4-5          
    # [55] reticulate_1.20       lifecycle_1.0.0       stringi_1.7.4         MASS_7.3-53           Rtsne_0.15            plyr_1.8.6           
    # [61] grid_4.0.3            parallel_4.0.3        listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1         crayon_1.4.1         
    # [67] miniUI_0.1.1.1        deldir_0.2-10         lattice_0.20-41       cowplot_1.1.1         splines_4.0.3         tensor_1.5           
    # [73] pillar_1.6.2          igraph_1.2.6          spatstat.geom_2.2-2   future.apply_1.8.1    reshape2_1.4.4        codetools_0.2-16     
    # [79] leiden_0.3.9          glue_1.4.2            packrat_0.7.0         data.table_1.14.0     png_0.1-7             vctrs_0.3.8          
    # [85] httpuv_1.6.3          gtable_0.3.0          RANN_2.6.1            purrr_0.3.4           spatstat.core_2.3-0   polyclip_1.10-0      
    # [91] tidyr_1.1.3           scattermore_0.7       future_1.22.1         assertthat_0.2.1      mime_0.11             xtable_1.8-4         
    # [97] RSpectra_0.16-0       later_1.3.0           survival_3.2-7        viridisLite_0.4.0     tibble_3.1.4          cluster_2.1.0        
    # [103] globals_0.14.0        fitdistrplus_1.1-5    ellipsis_0.3.2        ROCR_1.0-11 

# Analysis of Pten 10X data
# first step

# checked 4.10.2023

library (dplyr)
library (Seurat)
library (patchwork)
library (ggplot2)
library (GenBinomApps)
library (clusterProfiler)
library (org.Mm.eg.db)
library (openxlsx)

set.seed (2401)

# function to calculate relative abundances of sub-groups in different samples
# samples are given as split
# typically I use the seurat meta data df for this analysis
# grouping variable subdivides the split - typically I use the genotype as split and the group for clusters/cell types etc.
# the reference is the sample/split that should be used in significance calculations, typically the control

prep_rel_df <- function ( meta_data = meta_data, group = group, split = split, reference = reference) {
   
  # to maximise flexibility identify column indexes here
  group_idx <- which(colnames(meta_data) == group)
  split_idx <- which(colnames(meta_data) == split)
  
  # calculate relative abundance of cell broad types
  all_rel_group <- data.frame()
  for ( split in unique(meta_data[, split_idx]) ) {
    tmp <- reshape2::melt( table( meta_data[, group_idx][which(  meta_data[, split_idx] == split )] ) )
    tmp$rel <- tmp$value / sum(tmp$value)
    tmp$split <- split
    tmp$total <- sum(tmp$value)
    all_rel_group <- rbind(all_rel_group, tmp)
  }
  
  # confidence intervals are calculated here
  conf_int <- t(sapply (1:nrow(all_rel_group), function (x) {
    tmp <- clopper.pearson.ci(all_rel_group$value[x], all_rel_group$total[x], alpha = 0.05, CI = "two.sided")
  }))
  
  conf_int <- as.data.frame( conf_int )
  
  all_rel_group <- cbind( all_rel_group, conf_int )
  all_rel_group$Lower.limit <- as.numeric( all_rel_group$Lower.limit )
  all_rel_group$Upper.limit <- as.numeric( all_rel_group$Upper.limit )
  
  # do all pairwise comparisons
  pvalues <- data.frame()
  
  for (cluster in unique(all_rel_group$Var1)) {
    tmp <- all_rel_group[which(all_rel_group$Var1 == cluster),]
    
    # make sure that there both splits are >0
    if (nrow(tmp) == length(unique(all_rel_group$split))) {
      ctrl_idx <- which( tmp$split == reference )
      for ( split in unique(tmp$split)[ which( !unique(tmp$split) == reference ) ] ) {
        idx <- which(tmp$split == split)
        tmp_table <- as.table( rbind( c(tmp[ctrl_idx, "value"], tmp[ctrl_idx, "total"] - tmp[ctrl_idx, "value"]),
                                      c(tmp[idx, "value"], tmp[idx, "total"] - tmp[idx, "value"]) )
        )
        chisq_out <- chisq.test( tmp_table, correct = T, rescale.p = F, simulate.p.value = F )
        pvalues <- rbind( pvalues, data.frame(cluster = cluster, split = split, pvalue = chisq_out$p.value, stringsAsFactors = F) )
      }
    }
  }
  
  pvalues$p.adjust <- p.adjust( p = pvalues$pvalue, method = "BH" )
  
  pvalues$sig <- sapply(pvalues$p.adjust, function (x) {
    if (x < 0.001) {
      return ("***")
    } else if ( x < 0.01) {
      return ("**")
    } else if (x < 0.05) {
      return ("*")
    } else {
      return ("")
    }
  })
  
  pvalues$new_id <- paste(pvalues$cluster, pvalues$sig)
  
  all_rel_group <- merge( all_rel_group, pvalues[, c("cluster", "pvalue", "p.adjust", "sig", "new_id")], by.x="Var1", by.y="cluster")
    
  return (all_rel_group)
}

# define a base folder
base <- "./"

if (! file.exists ( paste(base, "/RDS_files/seurat_P0_ctrl_Pten_combined.RDS", sep=""))) {
  
  # count files
  # give full path to the output files from cellRanger: available on GEO
  count_files <- c("Pten_MBd_Ctrl_filtered_feature_bc_matrix.h5 ", "Pten_MBd_KO_filtered_feature_bc_matrix.h5")
  genotype_list <- c("ctrl", "sparse")
  seurat.object.list <- list()
  
  for ( i in c(1,2)) {
    
    # Initialize the Seurat object with the raw (non-normalized data).
    seurat.object <- CreateSeuratObject(counts = Read10X_h5( filename = count_files[i]), min.cells = 3, min.features = 200)
 
    # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
    seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "^mt-")
    
    # Visualize QC metrics as a violin plot
    VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
    # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
    
    plot1 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    plot1 + plot2
    
    seurat.object <- subset(seurat.object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
 
    plot1 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    plot1 + plot2
    
    seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 2000)
    
    seurat.object@meta.data$genotype <- genotype_list[i]
    
    seurat.object.list[[ genotype_list[i] ]] <- seurat.object
  }
  
  # $ctrl
  # An object of class Seurat 
  # 21603 features across 4212 samples within 1 assay 
  # Active assay: RNA (21603 features, 2000 variable features)
  # 
  # $sparse
  # An object of class Seurat 
  # 21748 features across 4713 samples within 1 assay 
  # Active assay: RNA (21748 features, 2000 variable features)
  # 
  saveRDS( object = seurat.object.list, file = paste(base, "/RDS_files/P0_Pten_ctrl_RAW.RDS", sep=""))
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = seurat.object.list)
  
  anchors <- FindIntegrationAnchors(object.list = seurat.object.list, anchor.features = features)
  
  # this command creates an 'integrated' data assay
  seurat.combined <- IntegrateData( anchorset = anchors )
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(seurat.combined) <- "integrated"
  
  saveRDS(object = seurat.combined, file = paste(base, "/RDS_files/seurat_P0_ctrl_Pten_combined.RDS", sep="") )
  stop ("raw data processing done. proceed with integration now")
}

seurat.combined <- readRDS( file = paste(base, "/RDS_files/seurat_P0_ctrl_Pten_combined.RDS", sep="") )

Idents( seurat.combined ) <- "Age"
qc_stats_df <- seurat.combined@meta.data
qc_stats_df <- qc_stats_df[,c("nFeature_RNA", "nCount_RNA", "percent.mt", "genotype")]

colnames(  seurat.combined@meta.data )
# [1] "nFeature_RNA" "nCount_RNA"   "percent.mt" 

# first classify cells based on e18 LaManno data
# note this was done on a different machine - mainly to remove cells from tissues other than Midbrain
LaManno_e18_metadata <- readRDS( file = paste(base, "integrate_Pten_e18/LaManno_e18_metadata.RDS", sep=""))

# details on production of these anchors in the integrate_Pten_e18 folder

# first classify tissue of origin
# first for control
embryo.anchors.ctrl <- readRDS( file = paste(base, "integrate_Pten_e18/Pten_ctrl_anchors.RDS", sep=""))
predictions_ctrl <- TransferData(anchorset = embryo.anchors.ctrl, refdata = LaManno_e18_metadata$Tissue, dims = 1:30)
# modify rownames to fit to merged dataset 
rownames(predictions_ctrl) <- paste( rownames(predictions_ctrl), "1", sep="_" )

# do the same for the sparse PtenKO
embryo.anchors.sparse <- readRDS( file = paste(base, "integrate_Pten_e18/Pten_sparse_anchors.RDS", sep=""))
predictions_sparse <- TransferData(anchorset = embryo.anchors.sparse, refdata = LaManno_e18_metadata$Tissue, dims = 1:30)
rownames(predictions_sparse) <- paste( rownames(predictions_sparse), "2", sep="_" )

predictions <- rbind( predictions_ctrl, predictions_sparse )

column2add <- predictions[1]
colnames(column2add) <- "Tissue"
seurat.combined <- AddMetaData( seurat.combined, metadata = column2add)

# second classify by cell type
# first for control
predictions_ctrl <- TransferData(anchorset = embryo.anchors.ctrl, refdata = LaManno_e18_metadata$cell_type, dims = 1:30)
# modify rownames to fit the merged dataset
rownames(predictions_ctrl) <- paste( rownames(predictions_ctrl), "1", sep="_" )
# same for the sparse PtenKO
predictions_sparse <- TransferData(anchorset = embryo.anchors.sparse, refdata = LaManno_e18_metadata$cell_type, dims = 1:30)
rownames(predictions_sparse) <- paste( rownames(predictions_sparse), "2", sep="_" )
predictions <- rbind( predictions_ctrl, predictions_sparse )

# create one dataframe to add to the Seurat object
column2add <- predictions[1]
colnames(column2add) <- "cell_type"
seurat.combined <- AddMetaData( seurat.combined, metadata = column2add)

# write out meta data including classification

wb <- createWorkbook()

sheetName <- "meta_data"
  addWorksheet(wb, sheetName)
  tmp <- seurat.combined@meta.data
  tmp$cellID <- rownames(tmp)
  writeData(wb, sheetName, tmp )
  addFilter(wb, sheetName, row = 1, cols = 1:ncol( tmp ))
  setColWidths(wb, sheetName, cols = 1:ncol( tmp ), widths="auto")

saveWorkbook(wb, file = paste(base, "/Supplement/Pten_analysis_meta_data.xlsx", sep=""), overwrite = T) 

table(seurat.combined$genotype)
# cell number before filtering
# ctrl sparse 
# 4212   4713 

# now do the filtering based on the classifications
Idents(seurat.combined) <- "Tissue"
seurat.combined <- subset( seurat.combined, idents = "Midbrain")

table(seurat.combined$genotype)
# cell number after tissue of origin filtering
# ctrl sparse 
# 3049   1982 

# remove contaminating cell types
Idents(seurat.combined) <- "cell_type"
seurat.combined <- subset( seurat.combined, idents = c("Immune", "Fibroblast", "Vascular", "Subcommissural organ", "Blood"), invert=T)

# we know that at e18 there are no RGPs any more
# the Ependymal as well as RGP cells cluster close to glioblasts
# so it is reasonable to assume these are glia progenitors
seurat.combined@meta.data$cell_type[ which(  seurat.combined@meta.data$cell_type %in% c("Radial glia", "Ependymal"))] <- "glia progenitor"

table(seurat.combined$genotype)
# cell number after removing minority populations
# this is the final cell number
# ctrl sparse 
# 2985   1889 

# quick calculation of minority/contaminating cell types
(1 - 2985 / 3049) *100
# 2.099049
(1 - 1889 / 1982) * 100
# 4.69223

# Run the standard workflow for visualization and clustering
seurat.combined <- ScaleData(seurat.combined, verbose = FALSE)
seurat.combined <- RunPCA(seurat.combined, npcs = 40, verbose = FALSE)

ElbowPlot( seurat.combined, ndims = 40)
ggsave( paste(base, "QC/Pten.first.Elbow.pdf", sep=""))

seurat.combined <- RunUMAP(seurat.combined, reduction = "pca", dims = 1:35)
seurat.combined <- FindNeighbors(seurat.combined, reduction = "pca", dims = 1:35)
seurat.combined <- FindClusters(seurat.combined, resolution = 0.5)

first_overview_umap <- DimPlot( seurat.combined, group.by = "seurat_clusters", label = T)
Pten_first_overview_umap_split <- DimPlot( seurat.combined, group.by = "cell_type", label = F, split.by = "genotype")

# print UMAP with different point sized to find optimal plot for the paper
for (pt.size in seq(0.5, 2, 0.5)) {
  Pten_first_overview_cellTypes <- DimPlot( seurat.combined, group.by = "cell_type", label = F, pt.size = pt.size)
  ggsave (filename = paste(base, "/Figures/Figure6O.", pt.size, ".jpg", sep="" ), plot = Pten_first_overview_cellTypes, width = 10, height = 10 )
}


seurat.combined@meta.data$super_group <- "glia"
seurat.combined@meta.data$super_group[ which( seurat.combined@meta.data$cell_type %in% c("Neuron", "Neuroblast"))] <- "neuron"

Idents(seurat.combined) <- "cell_type"
# calculate relative abundance changes and significances
all_rel_ttype <- prep_rel_df( meta_data = seurat.combined@meta.data, group = "cell_type", split = "genotype", reference="ctrl")

#         Var1 value         rel  split total Confidence.Interval Lower.limit Upper.limit alpha       pvalue     p.adjust sig              new_id
# 1  glia progenitor    12 0.004020101   ctrl  2985           two.sided 0.002078917 0.007011766  0.05 5.704317e-03 0.0095071944  **  glia progenitor **
# 2  glia progenitor    21 0.011116993 sparse  1889           two.sided 0.006894399 0.016943541  0.05 5.704317e-03 0.0095071944  **  glia progenitor **
# 3        Glioblast   827 0.437797777 sparse  1889           two.sided 0.415271885 0.460517105  0.05 8.806188e-01 0.8806187934              Glioblast 
# 4        Glioblast  1299 0.435175879   ctrl  2985           two.sided 0.417293833 0.453185137  0.05 8.806188e-01 0.8806187934              Glioblast 
# 5       Neuroblast   131 0.069348862 sparse  1889           two.sided 0.058303577 0.081752693  0.05 3.035968e-01 0.3794959930             Neuroblast 
# 6       Neuroblast   232 0.077721943   ctrl  2985           two.sided 0.068367346 0.087913793  0.05 3.035968e-01 0.3794959930             Neuroblast 
# 7           Neuron   436 0.230809952 sparse  1889           two.sided 0.211973882 0.250485413  0.05 4.495469e-05 0.0001123867 ***          Neuron ***
# 8           Neuron   848 0.284087102   ctrl  2985           two.sided 0.267959365 0.300639061  0.05 4.495469e-05 0.0001123867 ***          Neuron ***
# 9  Oligodendrocyte   594 0.198994975   ctrl  2985           two.sided 0.184805527 0.213776857  0.05 2.289651e-05 0.0001123867 *** Oligodendrocyte ***
# 10 Oligodendrocyte   474 0.250926416 sparse  1889           two.sided 0.231509640 0.271119465  0.05 2.289651e-05 0.0001123867 *** Oligodendrocyte ***

Idents(seurat.combined) <- "cell_type"
# check for skew in Neuroblast/Neuron abundance between genotypes
all_rel_NeuronNeuroblast <- prep_rel_df( meta_data = subset( seurat.combined, idents = c("Neuron", "Neuroblast"))@meta.data, group = "cell_type", split = "genotype", reference="ctrl")

#         Var1 value       rel  split total Confidence.Interval Lower.limit Upper.limit alpha    pvalue  p.adjust sig      new_id
# 1 Neuroblast   232 0.2148148   ctrl  1080           two.sided   0.1906619   0.2405321  0.05 0.4887982 0.4887982     Neuroblast 
# 2 Neuroblast   131 0.2310406 sparse   567           two.sided   0.1969422   0.2679688  0.05 0.4887982 0.4887982     Neuroblast 
# 3     Neuron   848 0.7851852   ctrl  1080           two.sided   0.7594679   0.8093381  0.05 0.4887982 0.4887982         Neuron 
# 4     Neuron   436 0.7689594 sparse   567           two.sided   0.7320312   0.8030578  0.05 0.4887982 0.4887982         Neuron 

ggplot( all_rel_NeuronNeuroblast, aes(x=split, y=rel)) + 
  geom_bar( stat="identity") + 
  geom_errorbar(aes(ymin=Lower.limit, ymax=Upper.limit)) +
  facet_grid(~Var1) +
  theme_classic()
ggsave( filename = paste(base, "/Figures/Figure6P.pdf", sep="" ) )

# save the raw data
all_rel_NeuronNeuroblast$Confidence.Interval <- unlist( all_rel_NeuronNeuroblast$Confidence.Interval )
all_rel_NeuronNeuroblast$alpha <- unlist( all_rel_NeuronNeuroblast$alpha )
write.csv( x = all_rel_NeuronNeuroblast, file = paste( base, "/Supplement/Figure6P.csv", sep=""))

# extract Neuroblasts and analyse
Pten.neuroblast.seurat <- subset( seurat.combined, idents = c("Neuroblast"))

table(Pten.neuroblast.seurat@meta.data$genotype)
# ctrl sparse 
# 232    131 

test_diff <- list()

Idents(Pten.neuroblast.seurat) <- "genotype"
test_diff[["NB"]] <- FindMarkers( object = Pten.neuroblast.seurat, ident.1 = "sparse", ident.2 = "ctrl", test.use = "MAST", logfc.threshold = 0.1)

# extract neurons only and perform DEG analysis
Idents(seurat.combined) <- "cell_type"
Pten.neuron.seurat <- subset( seurat.combined, idents = c("Neuron"))

Idents(Pten.neuron.seurat) <- "genotype"
test_diff[["Neuron"]] <- FindMarkers( object = Pten.neuron.seurat, ident.1 = "sparse", ident.2 = "ctrl", test.use = "MAST", logfc.threshold = 0.1)

# write out DEG analysis
wb <- createWorkbook()

for ( name in names(test_diff) ){
  sheetName <- name
  addWorksheet(wb, sheetName)
  tmp <- test_diff[[name]]
  tmp$gene <- rownames(tmp)
  writeData(wb, sheetName, tmp )
  addFilter(wb, sheetName, row = 1, cols = 1:ncol( tmp ))
  setColWidths(wb, sheetName, cols = 1:ncol( tmp ), widths="auto")
}

saveWorkbook(wb, file = paste(base, "/Supplement/Pten_DEG_list.xlsx", sep=""), overwrite = T) 


deg_list <- list()
for (cell_type in names(test_diff)) {
  deg_list[[cell_type]] <- list()
  for (dir in c("up", "down")) {
    if (dir == "up") {
      deg_list[[cell_type]][[dir]] <- rownames(test_diff[[cell_type]][which(test_diff[[cell_type]]$avg_log2FC > 0 & test_diff[[cell_type]]$p_val_adj < 0.01),])
    } else {
      deg_list[[cell_type]][[dir]] <- rownames(test_diff[[cell_type]][which(test_diff[[cell_type]]$avg_log2FC < 0 & test_diff[[cell_type]]$p_val_adj < 0.01),])
    }
  }
}

# plot number of DEGs
nDEG_plot_df <- data.frame()
for (cell_type in names(deg_list)) {
  for (dir in names(deg_list[[cell_type]])) {
    nDEG_plot_df <- rbind( nDEG_plot_df, data.frame( cell_type = cell_type, dir = dir, n = length(deg_list[[cell_type]][[dir]]) ) )
  }
}

ggplot( nDEG_plot_df, aes(x=cell_type, y=n, fill=dir)) + 
  geom_bar(stat="identity", position="dodge") + 
  theme_classic() + xlab("") + ylab("") + 
  ggtitle("#DEGs p.adj < 0.01 sparseKO / control")
ggsave( filename = paste(base, "/Figures/Figure6Q.pdf", sep="" ) )

write.csv( x = nDEG_plot_df, file = paste(base, "/Figures/Figure6Q.csv", sep=""))

# perform GO enrichment
eGO <- list()
for (cell_type in names(deg_list)) {
  for (dir in names(deg_list[[cell_type]])) {
    
    gene_list <- bitr( geneID = deg_list[[cell_type]][[dir]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db )
    eGO[[paste(cell_type, dir, sep="_")]] <- enrichGO( gene = gene_list$ENTREZID, OrgDb = org.Mm.eg.db, ont = "ALL", readable = T)
    
  }
}

# save for supplement
wb <- createWorkbook()

for (name in names(eGO)) {
  
  sheetName <- name
  addWorksheet(wb, sheetName)
  
  writeData(wb, sheetName, eGO[[name]] )
  addFilter(wb, sheetName, row = 1, cols = 1:ncol(eGO[[name]]))
  setColWidths(wb, sheetName, cols = 1:ncol(eGO[[name]]), widths="auto")
  
}

saveWorkbook(wb, file = paste(base, "/Supplement/Pten_GO_analysis.xlsx", sep=""), overwrite = T) 

#########################
###
# analyse GO enrichment
###
#########################

pattern_list <- list()
pattern_list[["INH"]] <- "GABA|interneuron|inhibitory"
pattern_list[["CC"]] <- "cell cycle"
pattern_list[["DEV"]] <- c("differentiation", "neuron")

go2plot <- data.frame()
go2genes <- list()
for (name in names(eGO)) {
  for (pattern in names(pattern_list)) {
    # filter for significance
    tmp <- eGO[[name]]@result[which(eGO[[name]]@result$p.adjust < 0.05),]
    # filter for Description tags
    idx <- c()
    for (i in 1:length(pattern_list[[pattern]])) {
      idx <- c( idx, which(grepl(pattern = pattern_list[[pattern]][i], x = tmp$Description)) )
    }
    idx <- as.numeric(names( which( table(idx) == length(pattern_list[[pattern]]) ) ) )
    
    tmp <- tmp[idx,]
    nGOterms <- nrow(tmp)
    go2genes[[paste(name, pattern, sep="_")]] <- unique( unlist(lapply(tmp$geneID, function (x) strsplit(x = x, split = "/")[[1]])) )
    
    go2plot <- rbind( go2plot, data.frame(group = name, pattern=pattern, nGOterms=nGOterms ) )
  }
}

# plot number of sig. enriched GO terms for main figure
ggplot( go2plot, aes(y=pattern, x=nGOterms)) + geom_bar(stat="identity") + 
  facet_wrap(~group, ncol = 1) + theme_classic() + ggtitle("# GO term in category") +
  xlab("") + ylab("")
ggsave( filename = paste(base, "/Figures/Figure6R.pdf", sep="" ), width = 2 )

# save the raw data
# save for supplement
wb <- createWorkbook()

sheetName <- "GO group count"
addWorksheet(wb, sheetName)
  
writeData(wb, sheetName, go2plot )
addFilter(wb, sheetName, row = 1, cols = 1:ncol( go2plot ))
setColWidths(wb, sheetName, cols = 1:ncol( go2plot ), widths="auto")

saveWorkbook(wb, file = paste(base, "/Supplement/Figure6R.xlsx", sep=""), overwrite = T) 

deg_plot_df <- data.frame()
# plot actual DEGs in the groups
for (name in names(go2genes)) {
  name_parts <- strsplit(x = name, split = "_")[[1]]
  tmp <- test_diff[[name_parts[1]]][go2genes[[name]],]
  tmp$score <- log10(tmp$p_val_adj)
  tmp$score[which(tmp$avg_log2FC > 0)] <- tmp$score[which(tmp$avg_log2FC > 0)] * -1
  tmp$cell_type <- name_parts[1]
  tmp$dir <- name_parts[2]
  tmp$group <- name_parts[3]
  tmp$gene_name <- rownames(tmp)
  deg_plot_df <- rbind( deg_plot_df, tmp )
}

tmp_df2plot <- deg_plot_df[which(deg_plot_df$cell_type == "NB"),]
tmp_df2plot$gene_name <- factor( tmp_df2plot$gene_name, levels = unique(tmp_df2plot$gene_name[ order( tmp_df2plot$score ) ]))

ggplot(tmp_df2plot, aes(x=gene_name, y=score, fill=group)) + 
  geom_bar(stat="identity", position="dodge") + 
  theme_classic() + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("DEG score of Neuroblast DEGs in GO groups")
ggsave( filename = paste(base, "/Supplement/FigureS8D.pdf", sep="" ), width=10 )

tmp_df2plot <- deg_plot_df[which(deg_plot_df$cell_type == "Neuron"),]
#tmp_df2plot$gene_name <- rownames(tmp_df2plot)
tmp_df2plot$gene_name <- factor( tmp_df2plot$gene_name, levels = unique(tmp_df2plot$gene_name[ order( tmp_df2plot$score ) ]))
ggplot(tmp_df2plot, aes(x=gene_name, y=score, fill=group)) + 
  geom_bar(stat="identity", position="dodge") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) + ggtitle("DEG score of Neurons DEGs in GO groups")
ggsave( filename = paste(base, "/Supplement/FigureS8E.pdf", sep="" ), width = 6 )

############################
###
# re-cluster neurons only
###
############################

Pten.neuron.seurat <- ScaleData( Pten.neuron.seurat, verbose = FALSE)
Pten.neuron.seurat <- RunPCA( Pten.neuron.seurat, npcs = 40, verbose = FALSE)

ElbowPlot( Pten.neuron.seurat, ndims = 40 )
ggsave( paste(base, "QC/Pten.Neuron.Elbow.pdf", sep=""))

Pten.neuron.seurat <- RunUMAP( Pten.neuron.seurat, reduction = "pca", dims = 1:30)

Pten.neuron.seurat <- FindNeighbors(Pten.neuron.seurat, reduction = "pca", dims = 1:30)
Pten.neuron.seurat <- FindClusters(Pten.neuron.seurat, resolution = 0.5)

Pten_neurons_overview_umap <- DimPlot( Pten.neuron.seurat, group.by = "seurat_clusters", label = T) + NoLegend()
Pten.neuron.seurat_overview_expression <- VlnPlot( object = Pten.neuron.seurat, features = c("Gad2", "Sox14" , "Slc17a6" , "Pou4f1"), ncol = 2 )

ggsave (filename = paste(base, "/Supplement/FigureS8F.pdf", sep="" ), plot = Pten_neurons_overview_umap, width = 5, height = 5 )
ggsave (filename = paste(base, "/Supplement/FigureS8G.pdf", sep="" ), plot = Pten.neuron.seurat_overview_expression, width = 5, height = 5 )

# prepare broad group
neuron_to_group <- c( "exc", "exc/inh", "exc", "inh", "inh", "inh", "inh")
names( neuron_to_group ) <- as.character(c(0:6))

Pten.neuron.seurat@meta.data$broad_group <- neuron_to_group [ as.character ( Pten.neuron.seurat@meta.data$seurat_clusters ) ]

# change to RGB Inh: 109/144/202, Exc:216/164/40
FeaturePlot( Pten.neuron.seurat, features = c("Gad2", "Slc17a6"), blend = T, 
             order=T, min.cutoff = "q10", max.cutoff = "q95", cols = c( "#6d90ca", "#d8a428"))
ggsave( filename = paste(base, "/Figures/Figure6S_S8H.jpg", sep="" ), width = 15, height = 5)

# calculate relative abundance changes and significances
all_rel_ttype <- prep_rel_df( meta_data = Pten.neuron.seurat@meta.data, group = "broad_group", split = "genotype", reference="ctrl")

#      Var1 value       rel  split total Confidence.Interval Lower.limit Upper.limit alpha       pvalue     p.adjust sig    new_id
# 1     exc   460 0.5424528   ctrl   848           two.sided   0.5082337   0.5763760  0.05 3.298177e-06 9.894530e-06 ***   exc ***
# 2     exc   176 0.4036697 sparse   436           two.sided   0.3572603   0.4513954  0.05 3.298177e-06 9.894530e-06 ***   exc ***
# 3 exc/inh   108 0.1273585   ctrl   848           two.sided   0.1056569   0.1516867  0.05 4.807460e-02 4.807460e-02   * exc/inh *
# 4 exc/inh    74 0.1697248 sparse   436           two.sided   0.1356844   0.2083205  0.05 4.807460e-02 4.807460e-02   * exc/inh *
# 5     inh   280 0.3301887   ctrl   848           two.sided   0.2985833   0.3629795  0.05 8.338996e-04 1.250849e-03  **    inh **
# 6     inh   186 0.4266055 sparse   436           two.sided   0.3796705   0.4745431  0.05 8.338996e-04 1.250849e-03  **    inh **
  
  
rel_plot <- ggplot( all_rel_ttype, aes(x=split, y=rel, fill=split)) + 
  geom_bar(stat="identity", position="dodge") + 
  geom_errorbar( aes( ymin = Lower.limit, ymax = Upper.limit) ) + 
  facet_wrap(~new_id, ncol = 4) + theme_classic() + theme(legend.position = "none") 

ggsave( filename = paste(base, "/Figures/Figure6T.pdf", sep="" ), plot = rel_plot, width = 5, height = 5)

all_rel_ttype$Confidence.Interval <- unlist(all_rel_ttype$Confidence.Interval)
all_rel_ttype$alpha <- unlist( all_rel_ttype$alpha )

write.csv( x = all_rel_ttype, file = paste( base, "/Supplement/Figure6T.csv", sep=""))

# save the neuron analysis for the next step
saveRDS( Pten.neuron.seurat, file = paste(base, "/RDS_files/Pten.neuron.seurat.RDS", sep="" ) )

### end of first step ###

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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_IE.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_IE.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_IE.UTF-8      
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_IE.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] openxlsx_4.2.5.1      org.Mm.eg.db_3.15.0   AnnotationDbi_1.58.0  IRanges_2.30.1        S4Vectors_0.34.0      Biobase_2.56.0        BiocGenerics_0.42.0   clusterProfiler_4.4.4
# [9] GenBinomApps_1.2      ggplot2_3.3.6         patchwork_1.1.2       sp_1.5-0              SeuratObject_4.1.2    Seurat_4.2.0          dplyr_1.0.10         
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.2                  reticulate_1.26             tidyselect_1.2.0            RSQLite_2.2.18              htmlwidgets_1.5.4           grid_4.2.1                 
# [7] BiocParallel_1.30.4         Rtsne_0.16                  scatterpie_0.1.8            munsell_0.5.0               codetools_0.2-18            ragg_1.2.4                 
# [13] ica_1.0-3                   future_1.28.0               miniUI_0.1.1.1              withr_2.5.0                 spatstat.random_2.2-0       colorspace_2.0-3           
# [19] GOSemSim_2.22.0             progressr_0.11.0            rstudioapi_0.14             SingleCellExperiment_1.18.1 ROCR_1.0-11                 tensor_1.5                 
# [25] DOSE_3.22.1                 listenv_0.8.0               MatrixGenerics_1.8.1        labeling_0.4.2              GenomeInfoDbData_1.2.8      polyclip_1.10-4            
# [31] bit64_4.0.5                 farver_2.1.1                downloader_0.4              parallelly_1.32.1           vctrs_0.5.0                 treeio_1.20.2              
# [37] generics_0.1.3              R6_2.5.1                    GenomeInfoDb_1.32.4         ggbeeswarm_0.6.0            graphlayouts_0.8.3          DelayedArray_0.22.0        
# [43] bitops_1.0-7                spatstat.utils_3.0-1        cachem_1.0.6                fgsea_1.22.0                gridGraphics_0.5-1          promises_1.2.0.1           
# [49] scales_1.2.1                ggraph_2.1.0                enrichplot_1.16.2           beeswarm_0.4.0              rgeos_0.5-9                 gtable_0.3.1               
# [55] globals_0.16.1              goftest_1.2-3               tidygraph_1.2.2             rlang_1.0.6                 systemfonts_1.0.4           splines_4.2.1              
# [61] lazyeval_0.2.2              spatstat.geom_3.0-3         reshape2_1.4.4              abind_1.4-5                 httpuv_1.6.6                qvalue_2.28.0              
# [67] tools_4.2.1                 ggplotify_0.1.0             ellipsis_0.3.2              spatstat.core_2.4-4         RColorBrewer_1.1-3          ggridges_0.5.4             
# [73] Rcpp_1.0.9                  plyr_1.8.7                  progress_1.2.2              zlibbioc_1.42.0             purrr_0.3.5                 RCurl_1.98-1.9             
# [79] prettyunits_1.1.1           rpart_4.1.16                deldir_1.0-6                pbapply_1.5-0               viridis_0.6.2               cowplot_1.1.1              
# [85] zoo_1.8-11                  SummarizedExperiment_1.26.1 ggrepel_0.9.1               cluster_2.1.3               magrittr_2.0.3              data.table_1.14.4          
# [91] scattermore_0.8             DO.db_2.9                   lmtest_0.9-40               RANN_2.6.1                  packrat_0.8.1               fitdistrplus_1.1-8         
# [97] matrixStats_0.62.0          hms_1.1.2                   mime_0.12                   xtable_1.8-4                gridExtra_2.3               compiler_4.2.1             
# [103] tibble_3.1.8                KernSmooth_2.23-20          crayon_1.5.2                shadowtext_0.1.2            htmltools_0.5.3             ggfun_0.0.7                
# [109] mgcv_1.8-40                 later_1.3.0                 tidyr_1.2.1                 aplot_0.1.8                 DBI_1.1.3                   tweenr_2.0.2               
# [115] MASS_7.3-57                 MAST_1.22.0                 Matrix_1.5-1                cli_3.4.1                   parallel_4.2.1              igraph_1.3.5               
# [121] GenomicRanges_1.48.0        pkgconfig_2.0.3             plotly_4.10.0               spatstat.sparse_3.0-0       ggtree_3.4.4                vipor_0.4.5                
# [127] XVector_0.36.0              yulab.utils_0.0.5           stringr_1.4.1               digest_0.6.30               sctransform_0.3.5           RcppAnnoy_0.0.20           
# [133] spatstat.data_3.0-0         Biostrings_2.64.1           leiden_0.4.3                fastmatch_1.1-3             tidytree_0.4.1              uwot_0.1.14                
# [139] shiny_1.7.3                 lifecycle_1.0.3             nlme_3.1-157                jsonlite_1.8.3              viridisLite_0.4.1           fansi_1.0.3                
# [145] pillar_1.8.1                lattice_0.20-45             ggrastr_1.0.1               KEGGREST_1.36.3             fastmap_1.1.0               httr_1.4.4                 
# [151] survival_3.3-1              GO.db_3.15.0                glue_1.6.2                  zip_2.2.2                   png_0.1-7                   bit_4.0.4                  
# [157] ggforce_0.4.1               stringi_1.7.8               blob_1.2.3                  textshaping_0.3.6           memoise_2.0.1               irlba_2.3.5.1              
# [163] future.apply_1.9.1          ape_5.6-2   

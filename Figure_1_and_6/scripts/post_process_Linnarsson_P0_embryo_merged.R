# produces all plots for Figure 1
# checked 4.10.2023

library (Seurat)
library (ggplot2)
library (clusterProfiler)
library (org.Mm.eg.db)
library (cowplot)
library (patchwork)
library (openxlsx)
library (GenBinomApps)
library (plyr)
library (ggrepel)
library (dplyr)
library (monocle3)
library (SeuratWrappers)

# simulate ggplot default colors to match UMAP coloring
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

set.seed(2401)

# define a base folder
base <- "./"

if (!file.exists(paste(base, "/RDS_files/LaManno_embryonic.combined.CCscored.RDS", sep=""))) {
  
  # details for preparation of this file in LaManno_integration folder
  seurat.combined <- readRDS( paste(base, "/LaManno_integration/seurat.combined.RDS", sep="") )
  
  seurat.combined <- ScaleData(seurat.combined, verbose = FALSE)
  seurat.combined <- RunPCA(seurat.combined, npcs = 30, verbose = FALSE)
  
  # determine dimensionality of the data
  ElbowPlot(seurat.combined, ndims = 30)
  ggsave(filename = paste(base, "/QC/initial_Elbow_plot.pdf", sep=""))
  
  # initial clustering
  seurat.combined <- RunUMAP(seurat.combined, reduction = "pca", dims = 1:26)
  seurat.combined <- FindNeighbors(seurat.combined, reduction = "pca", dims = 1:23)
  seurat.combined <- FindClusters(seurat.combined, resolution = 1.8)
  
  seurat.combined@meta.data$cell_type[which(is.na(seurat.combined@meta.data$cell_type))] <- "nan"
  
  DimPlot( object = seurat.combined, group.by = "seurat_clusters", label = T) + NoLegend()
  ggsave(filename = paste(base, "/QC/LaManno_EmbryoPool_integrated_initial_UMAP.png", sep=""))
  
  # determine cell cycle phase - as in Seurat vignette
  cell_cycle_phase_gene <- read.table(file = paste(base, "/other_data/Mus_musculus.csv", sep=""), sep = ",", header = T)
  gene_conversion <- bitr(geneID = cell_cycle_phase_gene$geneID, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
  cell_cycle_phase_gene <- merge(cell_cycle_phase_gene, gene_conversion, by.x="geneID", by.y="ENSEMBL")
  
  #score cell cycle
  seurat.combined <- CellCycleScoring(
    object = seurat.combined,
    g2m.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "G2/M"), "SYMBOL"],
    s.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "S"), "SYMBOL"]
  )
  
  # plot number of genes detected in each cluster
  cdata <- ddply(seurat.combined@meta.data, c("seurat_clusters"), summarise,
                 N    = length(nFeature_RNA),
                 mean_gene = mean(nFeature_RNA),
                 sd_gene   = sd(nFeature_RNA),
                 se_gene   = sd_gene / sqrt(N)
  )
  
  cdata$seurat_clusters <- factor( cdata$seurat_clusters, levels = cdata$seurat_clusters[order(cdata$mean_gene, decreasing = T)])
  
  ggplot(cdata, aes(x=seurat_clusters, y=mean_gene)) + 
    geom_bar(stat="identity") + geom_errorbar(aes(ymin=mean_gene-se_gene, ymax=mean_gene+se_gene)) + theme_classic()
  
  ggsave(filename = paste(base, "/QC/nFeature_dist.pdf", sep=""))
  
  # remove two clusters with the least genes detected - indicates a technical issue as the reason for clustering
  Idents (seurat.combined) <- "seurat_clusters"
  seurat.combined <- subset( seurat.combined, idents=c("46", "30"), invert=T)
 
  # re-do PCA, and clustering
  seurat.combined <- RunPCA( seurat.combined, npcs = 50, verbose = FALSE )
  seurat.combined <- RunUMAP( seurat.combined, reduction = "pca", dims = 1:30 )
  seurat.combined <- FindNeighbors(seurat.combined, reduction = "pca", dims = 1:30)
  seurat.combined <- FindClusters(seurat.combined, resolution = 1.3)
  DimPlot( seurat.combined, label = T ) + NoLegend()
  
  # the Seurat object:
  # An object of class Seurat 
  # 24755 features across 32084 samples within 2 assays 
  # Active assay: integrated (2000 features, 2000 variable features)
  # 1 other assay present: RNA
  # 2 dimensional reductions calculated: pca, umap
  
  saveRDS(object = seurat.combined, file = paste(base, "/RDS_files/LaManno_embryonic.combined.CCscored.RDS", sep=""))
  ggsave(filename = paste(base, "/QC/LaManno_EmbryoPool_integrated_final_UMAP.png", sep=""))

  ##################################################
  ###
  #   write out files needed for velocity analysis 
  ###
  #################################################
  
  Idents( seurat.combined ) <- "Age"
  embryo_subset <- subset(seurat.combined, idents = "embryonic_pool")
  # save components of the Seurat object for analysis with Velocyto - following https://github.com/basilkhuder/Seurat-to-RNA-Velocity
  write.csv(Cells( embryo_subset ), file =  paste(base, "velocity_analysis_files/cellID_obs.csv", sep=''))
  write.csv(Embeddings(embryo_subset, reduction = "umap"), file =paste(base, "velocity_analysis_files/cell_embeddings.csv", sep=''))
  #link cellID to some meta data
  meta_data <- data.frame( seurat_clusters = embryo_subset@meta.data$seurat_clusters )
  meta_data$cellID <- rownames(embryo_subset@meta.data)
  write.csv(meta_data, file = paste(base, "velocity_analysis_files/meta_data.csv", sep=''))
  
  # I restart R at this point to free memory
  stop("first step finished - restart")
}

# read the filtered Seurat object
seurat.combined <- readRDS( file = paste(base, "/RDS_files/LaManno_embryonic.combined.CCscored.RDS", sep=""))

# Note: In the paper we give the number of Reference cells only!
# seurat.combined also includes the Fzd10 cells!
# here the calculation of the Reference only cells
table(seurat.combined@meta.data$Age)["embryonic_pool"]

# embryonic_pool 
#           5552

32084 - 5552
# [1] 26532

# save numbers of cells per Age - NOTE: this includes the embryo pool data, which is not shown in the paper
age_cellNumber_dist <- reshape2::melt( table ( seurat.combined@meta.data$Age) )

wb <- createWorkbook()
sheetName <- "overlap analysis"
addWorksheet(wb, sheetName)

writeData(wb, sheetName, age_cellNumber_dist )
addFilter(wb, sheetName, row = 1, cols = 1:ncol(age_cellNumber_dist))
setColWidths(wb, sheetName, cols = 1:ncol(age_cellNumber_dist), widths="auto")

saveWorkbook(wb, file = paste(base, "/Supplement/FigureS1A.xlsx", sep=""), overwrite = T) 

# QC Plots for embryo pool data
Idents( seurat.combined ) <- "Age"
qc_stats_df <- subset( seurat.combined, idents = "embryonic_pool" )@meta.data
qc_stats_df <- qc_stats_df[,c("nFeature_RNA", "nCount_RNA", "percent.mt")]

colnames( qc_stats_df )
qc_plot1 <- ggplot( qc_stats_df, aes(x="pool", y=nFeature_RNA)) + geom_violin() + geom_jitter( size = 0.1 ) + ggtitle("genes") + theme_classic()
qc_plot2 <- ggplot( qc_stats_df, aes(x="pool", y=nCount_RNA)) + geom_violin() + geom_jitter( size = 0.1 ) + ggtitle("transcripts") + theme_classic()
qc_plot3 <- ggplot( qc_stats_df, aes(x="pool", y=percent.mt)) + geom_violin() + geom_jitter( size = 0.1 ) + ggtitle("% mt") + theme_classic()

comb_qc_plot <- qc_plot1 + qc_plot2 + qc_plot3

ggsave(filename = paste(base, "/Supplement/FigureS2MNO.png", sep=""), plot = comb_qc_plot, width = 7)
ggsave(filename = paste(base, "/Supplement/FigureS2MNO.pdf", sep=""), plot = comb_qc_plot, width = 7)

# add an inh/exc tag where possible
# NOTE: this is ONLY done for cells with the respective tag in "subclass"
seurat.combined@meta.data$broad_group <- "unknown"
seurat.combined@meta.data$broad_group[ which( grepl(pattern = "glutamatergic", x = seurat.combined@meta.data$subclass))] <- "exc"
seurat.combined@meta.data$broad_group[ which( grepl(pattern = "GABA", x = seurat.combined@meta.data$subclass))] <- "inh"

# combine cell type and neuro transmitter type
seurat.combined@meta.data$broad_group1 <- paste(seurat.combined@meta.data$cell_type, seurat.combined@meta.data$broad_group, sep="_")
seurat.combined@meta.data$broad_group1[which(seurat.combined@meta.data$broad_group == "unknown")] <- "unknown"
seurat.combined@meta.data$broad_group1[ which( seurat.combined@meta.data$cell_type == "Radial glia") ] <- "RGP"

# plot distribution of cell type and Neurons
all_age_rel <- data.frame()
all_age_Neuron_type_rel <-data.frame()

for ( age in unique( seurat.combined@meta.data$Age ) ) {

  # only for LaManno data!
  if (!age == "embryonic_pool") {
    
    # all values are relative to all cells from the respective age
    tmp <- seurat.combined@meta.data [ which( seurat.combined@meta.data$Age %in%  age ),]
    
    # use broad_group1 here - Neurons with neurotransmitter tag
    tmp_abs <- reshape2::melt( table(tmp$broad_group1) )
    tmp_abs$total <- sum( tmp_abs$value )
    tmp_abs$rel <- tmp_abs$value / sum( tmp_abs$value )
    tmp_abs$age <- age
    
    # NOTE: this number is relative to to the total number of cells at this age!
    all_age_Neuron_type_rel <- rbind( all_age_Neuron_type_rel, tmp_abs[which(tmp_abs$Var1 %in% c("Neuron_exc", "Neuron_inh")),])
    
    # use cell type here - broadest category
    tmp <- reshape2::melt( table(tmp$cell_type) )
    
    tmp$rel <- tmp$value / sum(tmp$value)
    tmp$age <- age
    
    tmp$total <- sum( tmp$value )
    
    all_age_rel <- rbind(all_age_rel, tmp)
    
  }
  
} 

# make the age an integer - for correct x-axis scale
all_age_rel$age <- gsub( pattern = "e", replacement = "", x = all_age_rel$age )
all_age_rel$age <- as.numeric(all_age_rel$age)

all_age_Neuron_type_rel$age <- gsub( pattern = "e", replacement = "", x = all_age_Neuron_type_rel$age )
all_age_Neuron_type_rel$age <- as.numeric(all_age_Neuron_type_rel$age)

# add confidence interval as in Cadwell et al., Elife 2020
conf_int <- t(sapply (1:nrow(all_age_rel), function (x) {
  
  tmp <- clopper.pearson.ci(all_age_rel$value[x], all_age_rel$total[x], alpha = 0.05, CI = "two.sided")
  
}))

all_age_rel <- cbind(all_age_rel, conf_int)

all_age_rel$Lower.limit <- as.numeric(all_age_rel$Lower.limit)
all_age_rel$Upper.limit <- as.numeric(all_age_rel$Upper.limit)

ggplot( all_age_rel, aes(x=age, y=rel, color=Var1, group = Var1)) + 
  geom_line() + geom_point() + theme_classic() + geom_errorbar(aes(ymin=Lower.limit, ymax=Upper.limit))

ggsave(filename = paste(base, "/Figures/Figure1C.pdf", sep=""), width = 6)


conf_int <- t(sapply (1:nrow(all_age_Neuron_type_rel), function (x) {
  
  tmp <- clopper.pearson.ci(all_age_Neuron_type_rel$value[x], all_age_Neuron_type_rel$total[x], alpha = 0.05, CI = "two.sided")
  
}))

all_age_Neuron_type_rel <- cbind(all_age_Neuron_type_rel, conf_int)

all_age_Neuron_type_rel$Lower.limit <- as.numeric(all_age_Neuron_type_rel$Lower.limit)
all_age_Neuron_type_rel$Upper.limit <- as.numeric(all_age_Neuron_type_rel$Upper.limit)

ggplot( all_age_Neuron_type_rel, aes(x=age, y=rel, color=Var1, group = Var1)) + 
  geom_line() + geom_point() + theme_classic() + geom_errorbar(aes(ymin=Lower.limit, ymax=Upper.limit))

ggsave(filename = paste(base, "/Figures/Figure1E.pdf", sep=""), width = 6)

# write table out for Giselle to plot herself
wb <- createWorkbook()
sheetName <- "LaManno all"
addWorksheet(wb, sheetName)

writeData(wb, sheetName,all_age_rel )
addFilter(wb, sheetName, row = 1, cols = 1:ncol(all_age_rel))
setColWidths(wb, sheetName, cols = 1:ncol(all_age_rel), widths="auto")

sheetName <- "LaManno neurons"
addWorksheet(wb, sheetName)

writeData(wb, sheetName, all_age_Neuron_type_rel )
addFilter(wb, sheetName, row = 1, cols = 1:ncol(all_age_Neuron_type_rel))
setColWidths(wb, sheetName, cols = 1:ncol(all_age_Neuron_type_rel), widths="auto")

saveWorkbook(wb, file = paste(base, "/Figures/Figure1CE.xlsx", sep=""), overwrite = T) 

#######################################
###
# plot overlap of our data and LaManno
###
#######################################

# save the initial clustering
seurat.combined@meta.data$seurat_clusters_high <- seurat.combined@meta.data$seurat_clusters

# do a low-res clustering
seurat.combined <- FindClusters(seurat.combined, resolution = 0.2)

umap_df <- as.data.frame( Embeddings(seurat.combined, reduction = "umap") )
umap_df$cellID <- rownames( umap_df )
tmp_df <- seurat.combined@meta.data[,c("Age", "seurat_clusters", "seurat_clusters_high", "broad_group", "broad_group1", "cell_type")]
tmp_df$cellID <- rownames(tmp_df)
umap_df <- merge ( tmp_df, umap_df , by="cellID" )
umap_df$origin <- "LaManno"
umap_df$origin[ which( umap_df$Age == "embryonic_pool" ) ] <- "embryo_pool"

ggplot( umap_df, aes(x=UMAP_1, y=UMAP_2, color=origin )) + 
  geom_point( size = 0.5 ) + 
  scale_color_manual( values = c( "LaManno" = "grey80", "embryo_pool" = "black")) +
  scale_size_manual( values = c( "LaManno" = 1.0, "embryo_pool" = 1.5 )) +
  theme_classic()

ggsave(filename = paste(base, "/Supplement/FigureS2P.png", sep=""), width = 6)

# plot without legend for preparing UMAPs with similar size
ggplot( umap_df, aes(x=UMAP_1, y=UMAP_2, color=origin )) + 
  geom_point( size = 0.5 ) + 
  scale_color_manual( values = c( "LaManno" = "grey80", "embryo_pool" = "black")) +
  scale_size_manual( values = c( "LaManno" = 1.0, "embryo_pool" = 1.5 )) +
  theme_classic() + theme(legend.position = "none")
ggsave(filename = paste(base, "/Supplement/FigureS2P_noLegend.png", sep=""), width = 6)

# plot UMAP based on larger clusters
cluster2color <- gg_color_hue( length(unique(seurat.combined@meta.data$seurat_clusters)))
names( cluster2color) <- as.character( sort( unique(seurat.combined@meta.data$seurat_clusters) ) )

# calculate positions of labels
tmp <- umap_df[which(!umap_df$Age == "embryonic_pool"),]
centroid_umap1 <- sapply( unique(tmp$seurat_clusters), function (x) median( umap_df$UMAP_1[which(umap_df$seurat_clusters == x)] ) )
centroid_umap2 <- sapply( unique(tmp$seurat_clusters), function (x) median( umap_df$UMAP_2[which(umap_df$seurat_clusters == x)] ) )

cluster_label_df <- data.frame( UMAP_1 = centroid_umap1, UMAP_2 = centroid_umap2, label = as.character(unique(tmp$seurat_clusters)) )

# NOTE: This is the only plot with Frz cells
lumped_clusters_plot <- ggplot( ) + 
  geom_point( data = umap_df[which(umap_df$Age == "embryonic_pool"),], aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters), size = 0.5 ) + 
  # geom_text(data = cluster_label_df, aes(x=UMAP_1, y=UMAP_2, label=label)) +
  scale_color_manual( values = cluster2color) +
  theme_classic() + theme( legend.position = "none" )

rel_cell_type_df <- data.frame()
for ( cluster in unique(seurat.combined@meta.data$seurat_clusters) ) {

  reference <- length( which( seurat.combined@meta.data$seurat_clusters %in% cluster & !seurat.combined@meta.data$Age == "embryonic_pool") )
  embryo_pool <- length( which( seurat.combined@meta.data$seurat_clusters %in% cluster & seurat.combined@meta.data$Age == "embryonic_pool") )
  
  rel_cell_type_df <- rbind( rel_cell_type_df, data.frame(reference = reference, Fzd = embryo_pool, clusterID = cluster) )
  
}

rel_cell_type_df$ref_rel <- rel_cell_type_df$reference / sum( rel_cell_type_df$reference )
rel_cell_type_df$rel_Fzd <- rel_cell_type_df$Fzd / sum( rel_cell_type_df$Fzd )

cor_fraction <- ggplot( rel_cell_type_df, aes(x=ref_rel, y=rel_Fzd, label=clusterID)) + 
  geom_point() + coord_fixed() + geom_text_repel() + theme_classic() + 
  xlab("fraction reference cells") + ylab("fraction Frz+ cells")

lumped_clusters_plot + cor_fraction
ggsave(filename = paste(base, "/Supplement/FigureS2QR.pdf", sep=""), width = 8)

cor_test_df <- cor.test( rel_cell_type_df$ref_rel, rel_cell_type_df$rel_Fzd)

# Pearson's product-moment correlation
# 
# data:  rel_cell_type_df$ref_rel and rel_cell_type_df$rel_Fzd
# t = 4.4265, df = 10, p-value = 0.001281
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4498668 0.9458845
# sample estimates:
#       cor 
# 0.8136904 

# save as xlsx for Giselle to analyse further
# write table out for Giselle to plot herself
wb <- createWorkbook()
sheetName <- "overlap analysis"
addWorksheet(wb, sheetName)

writeData(wb, sheetName, rel_cell_type_df )
addFilter(wb, sheetName, row = 1, cols = 1:ncol(rel_cell_type_df))
setColWidths(wb, sheetName, cols = 1:ncol(rel_cell_type_df), widths="auto")

sheetName <- "correlation test"
addWorksheet(wb, sheetName)

tmp <- as.data.frame( unlist(cor_test_df) )
tmp$description <- rownames( tmp )
writeData(wb, sheetName, tmp )
addFilter(wb, sheetName, row = 1, cols = 1:ncol(tmp))
setColWidths(wb, sheetName, cols = 1:ncol(tmp), widths="auto")

saveWorkbook(wb, file = paste(base, "/Supplement/FigureS2R.xlsx", sep=""), overwrite = T) 

#########
##
# plot broad groups
##
#########

# classify - give each Seurat cluster a tag
cluster2broadGroup <- sapply ( unique(umap_df$seurat_clusters), function (cluster) {
  tmp <- umap_df[which(umap_df$seurat_clusters == cluster),]
  tmp <- reshape2::melt( table(tmp$broad_group1) )
  # remove unknown cells
  tmp <- tmp[which(!tmp$Var1 == "unknown"),]
  # majority vote for broad group
  as.character( tmp[which( tmp$value == max(tmp$value)), "Var1"] )
})

names( cluster2broadGroup ) <- unique(umap_df$seurat_clusters)
umap_df$broad_group <- cluster2broadGroup[ as.character(umap_df$seurat_clusters)]

# for this set of figures use only LaManno cells for plotting
ggplot( umap_df[which(umap_df$origin == "LaManno"),], aes(x=UMAP_1, y=UMAP_2, color=broad_group)) + 
  geom_point( size = 0.5) + 
  scale_color_manual( values = c( "Neuron_exc" = "#d8a428", "Neuron_inh" = "#6d90ca", 
                                  "Neuroblast_exc" = "grey80", "Neuroblast_inh" = "grey80",
                                  "RGP" = "grey80")) + 
  theme_classic()

ggsave(filename = paste(base, "/Figures/Figure1D.png", sep=""))

# plot without legend
ggplot( umap_df[which(umap_df$origin == "LaManno"),], aes(x=UMAP_1, y=UMAP_2, color=broad_group)) + 
  geom_point( size = 0.5) + 
  scale_color_manual( values = c( "Neuron_exc" = "#d8a428", "Neuron_inh" = "#6d90ca", 
                                  "Neuroblast_exc" = "grey80", "Neuroblast_inh" = "grey80",
                                  "RGP" = "grey80")) +
  theme_classic() + theme(legend.position = "none") + coord_fixed()

ggsave(filename = paste(base, "/Figures/Figure1D_NoLegend.png", sep=""), width = 6)

# color broad cell types
ggplot( umap_df[which(umap_df$origin == "LaManno"),], aes(x=UMAP_1, y=UMAP_2, color=cell_type)) + 
  geom_point( size = 0.5 ) + 
  scale_color_manual( values = c( "Radial glia" = "#FF8080", "Neuroblast" = "#0F99B2", "Neuron" = "violet")) +
  # scale_size_manual( values = c( "unknown" = 0.5, "exc" = 1.5, "inh" = 1.5 )) +
  theme_classic() + coord_fixed()

ggsave(filename = paste(base, "/Figures/Figure1B.png", sep=""), width = 6)

# color broad cell types
ggplot( umap_df[which(umap_df$origin == "LaManno"),], aes(x=UMAP_1, y=UMAP_2, color=cell_type)) + 
  geom_point( size = 0.5 ) + 
  scale_color_manual( values = c( "Radial glia" = "#FF8080", "Neuroblast" = "#0F99B2", "Neuron" = "violet")) +
  # scale_size_manual( values = c( "unknown" = 0.5, "exc" = 1.5, "inh" = 1.5 )) +
  theme_classic() + theme(legend.position = "none") + coord_fixed()

ggsave(filename = paste(base, "/Figures/Figure1B_NoLegend.png", sep=""), width = 6)

# plot some marker genes for overview in supplement - Note: not all used for final figure
Idents(seurat.combined) <- "Age"

# use only LaManno cells for plotting
LaManno.seurat <- subset( seurat.combined, idents="embryonic_pool", invert=T)

# An object of class Seurat 
# 24755 features across 26532 samples within 2 assays 
# Active assay: integrated (2000 features, 2000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

FeaturePlot( object = LaManno.seurat, 
             features = c("Gad2", "Sox14", "Lhx1", "Pou4f1", "Slc17a6", "Nrn1", "Mki67", "Sox2", "Hes5"), 
             order = T, min.cutoff = "q25", slot = "scale.data", ncol = 3) & NoLegend()
ggsave( paste(base, "/Supplement/FigureS1B.png", sep=""), width=10, height=10 )

DimPlot( object = LaManno.seurat, group.by = "seurat_clusters", label=T) + NoLegend()
ggsave( paste(base, "/Supplement/FigureS1C.png", sep=""), width=6 )

DimPlot( object = LaManno.seurat, group.by = "seurat_clusters", label=F) + NoLegend()
ggsave( paste(base, "/Supplement/FigureS1C_noLabel.png", sep=""), width=6 )


##################################################
###
#   monocle3 analysis of the whole embryo dataset
###
#################################################

# this chunk can be used to test multiple parameters for Monocle3 analysis
do_test = F

if (do_test) {
  
  Idents( seurat.combined ) <- "seurat_clusters"
   
  cds <- as.cell_data_set(seurat.combined )
  cds <- cluster_cells(cds)
  plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
  
  i <- 0
  for (gdr in seq(0.1, 0.3, 0.05)) {
    for (edr in seq(1, 2, 0.25)) {
  
      message(i, " ", gdr, " ", edr)
      
      lgc <- list(ncenter=1000, nn.k=25, geodesic_distance_ratio=gdr, euclidean_distance_ratio=edr)
      
      cds <- as.cell_data_set(seurat.combined )
      cds <- cluster_cells(cds)
      cds <- learn_graph( cds, learn_graph_control = lgc, verbose = F, close_loop = T, use_partition = T )
      
      traj_plot <- plot_cells( cds )
      ggsave(filename = paste(base, "/QC/traj_test.noloop.", i, ".jpg", sep=""), plot = traj_plot)
  
      i<-i+1
     }
  }
  
  ?learn_graph
  ?set_nn_control
}

# do the actual Monocle3 analysis
lgc <- list(ncenter=1000, nn.k=25, geodesic_distance_ratio=0.25, euclidean_distance_ratio=1.5)

cds <- as.cell_data_set(seurat.combined )
cds <- cluster_cells(cds)
cds <- learn_graph( cds, learn_graph_control = lgc, verbose = F, close_loop = T, use_partition = T )

# Note: start points chosen arbitrarily, but can be seen on figures
cds <- order_cells( cds, reduction_method = "UMAP" )

# use only LaManno cells for plotting
traj_plot <- plot_cells( cds[,rownames(seurat.combined@meta.data[which(!seurat.combined$Age == "embryonic_pool"),])], color_cells_by = "seurat_clusters", graph_label_size = 6, label_branch_points = F, 
                         show_trajectory_graph = T, label_cell_groups = F )
ggsave(filename = paste(base, "/Supplement/FigureS1D.pdf", sep=""), width = 8, height = 6, plot = traj_plot)

traj_plot <- plot_cells( cds[,rownames(seurat.combined@meta.data[which(!seurat.combined$Age == "embryonic_pool"),])], color_cells_by = "pseudotime", graph_label_size = 3, show_trajectory_graph = T, label_branch_points = F )

ggsave(filename = paste(base, "/Figures/Figure1F.pdf", sep=""), width = 8, height = 6, plot = traj_plot)

#############################
###
# link to adult cell types
###
#############################

# only use e18 cells with a Neuron tag - only LaManno data
# extract cell names
cell_ids <- rownames( seurat.combined@meta.data[ ( which( seurat.combined@meta.data$Age %in% c("e18.0") & seurat.combined@meta.data$cell_type == "Neuron") ),] )
length(cell_ids)
# [1] 3741

# extract raw data and start from scratch for the integration analysis
e18_Neurons <- CreateSeuratObject(counts = GetAssayData( object = seurat.combined, slot = "count", assay = "RNA")[,cell_ids], project = "e18_Neurons", min.cells = 3, min.features = 200)
e18_Neurons <- NormalizeData( e18_Neurons, verbose = F)
e18_Neurons <- FindVariableFeatures( e18_Neurons, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# An object of class Seurat 
# 17149 features across 3741 samples within 1 assay 
# Active assay: RNA (17149 features, 2000 variable features)

# this RDS file comes from PatchSeq analysis
Zeisel_reference <- readRDS( "~/Documents/PatchSeq_R403_v1/final_1st_submission/RDS_files/cti_Seurat_RAW.RDS")

# An object of class Seurat 
# 15081 features across 3685 samples within 1 assay 
# Active assay: RNA (15081 features, 0 variable features)
# NOTE: this cell number is after creating a Seurat object, which is anadditional filter!

Zeisel_reference <- NormalizeData( Zeisel_reference, verbose = F)
Zeisel_reference <- FindVariableFeatures( Zeisel_reference, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
Zeisel_reference <- ScaleData(Zeisel_reference, verbose = FALSE)
Zeisel_reference <- RunPCA(Zeisel_reference, npcs = 30, verbose = FALSE)
Zeisel_reference <- RunUMAP(Zeisel_reference, reduction = "pca", dims = 1:20, verbose = FALSE)

Zeisel_MBd_UMAP <- DimPlot(Zeisel_reference, group.by = "tax_5", label = T) + NoLegend()
Zeisel_MBd_UMAP <- Zeisel_MBd_UMAP + ggtitle("Zeisel et al. MBd cell types")

ggsave(filename = paste(base, "/Supplement/Zeisel_MBd_UMAP.pdf", sep=""), width = 8, height = 6, plot = Zeisel_MBd_UMAP)

# do label transfer
transfer.anchors <- FindTransferAnchors(reference = Zeisel_reference, query = e18_Neurons,
                                        dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = transfer.anchors, refdata = Zeisel_reference$tax_5,
                            dims = 1:20)

# add predictions to Seurat object
LaManno.seurat <- AddMetaData(LaManno.seurat, metadata = predictions)

# calculate the average score for each cell type in each cluster
# we are interested in end-points, so do not analyse RGPs and immature neuron clusters 0, 6, 11
mat <- sapply ( as.character( c(1,2,3,4,5,7,8,9,10) ), function (cluster) {
  sapply( as.character( colnames(LaManno.seurat@meta.data)[22:37] ), function (prediction){
    mean( LaManno.seurat@meta.data[ which( LaManno.seurat@meta.data$seurat_clusters  == cluster), prediction], na.rm = T )
  })
})

rownames(mat) <- gsub(pattern = "prediction.score.", replacement = "", x = rownames(mat))

# determine cluster with max average score for this cell type
cellType_cluster <- reshape2::melt ( apply(mat,1, function (x) which(x==max(x)) ) )
# not that the above analysis only gives the index of the cluster, not the real name!
cellType_cluster$real_cluster_name <- c(1,2,3,4,5,7,8,9,10)[cellType_cluster$value]
cellType_cluster$id <- gsub(pattern = "prediction.score.", replacement = "", x = rownames(cellType_cluster))

saveRDS( object = cellType_cluster, file = paste(base, "/RDS_files/linktoadult.RDS", sep="") )

adult2cluster <- sapply( as.character(unique(cellType_cluster$real_cluster_name)), function (x) paste(cellType_cluster[which(cellType_cluster$real_cluster_name == x), "id"], collapse=","))
LaManno.seurat@meta.data$adult_link <- adult2cluster[ as.character(LaManno.seurat@meta.data$seurat_clusters) ]

adult_cellType_extrapolation <- pheatmap::pheatmap(mat, scale="row", silent = T)
pheatmap::pheatmap(mat, scale="row", silent = T, filename =  paste(base, "/Figures/Figure1H.pdf", sep=""))

LaManno_broad_seurat_clusters <- DimPlot(object = LaManno.seurat, group.by = "seurat_clusters", label = T) + NoLegend() 

#match the link names to default colors
manual_colors <- gg_color_hue(12)
names(manual_colors) <- as.character(c(0:11))
manual_colors <- manual_colors[ names(adult2cluster) ]
names(manual_colors) <- as.character( adult2cluster )

# do the plot
LaManno_adult_link_UMAP <- DimPlot(object = LaManno.seurat, group.by = "adult_link", label = T) + NoLegend() 
LaManno_adult_link_UMAP <- LaManno_adult_link_UMAP + scale_color_manual( values = manual_colors)

# save separately!
ggsave(filename = paste(base, "/Supplement/FigureS1F.pdf", sep=""), width = 6, height = 5, plot = Zeisel_MBd_UMAP)
ggsave(filename = paste(base, "/Figures/Figure1I.pdf", sep=""), width = 6,  height = 5, plot = LaManno_adult_link_UMAP)

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
#   [1] shiny_1.7.3                 SeuratWrappers_0.3.1        monocle3_1.2.9              SingleCellExperiment_1.18.1 SummarizedExperiment_1.26.1 GenomicRanges_1.48.0       
# [7] GenomeInfoDb_1.32.4         MatrixGenerics_1.8.1        matrixStats_0.62.0          dplyr_1.0.10                ggrepel_0.9.1               plyr_1.8.7                 
# [13] GenBinomApps_1.2            openxlsx_4.2.5.1            patchwork_1.1.2             cowplot_1.1.1               org.Mm.eg.db_3.15.0         AnnotationDbi_1.58.0       
# [19] IRanges_2.30.1              S4Vectors_0.34.0            Biobase_2.56.0              BiocGenerics_0.42.0         clusterProfiler_4.4.4       ggplot2_3.3.6              
# [25] sp_1.5-0                    SeuratObject_4.1.2          Seurat_4.2.0               
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.2             R.utils_2.12.0         reticulate_1.26        lme4_1.1-30            tidyselect_1.2.0       RSQLite_2.2.18         htmlwidgets_1.5.4      grid_4.2.1            
# [9] BiocParallel_1.30.4    Rtsne_0.16             scatterpie_0.1.8       munsell_0.5.0          ragg_1.2.4             codetools_0.2-18       ica_1.0-3              future_1.28.0         
# [17] miniUI_0.1.1.1         withr_2.5.0            spatstat.random_2.2-0  colorspace_2.0-3       GOSemSim_2.22.0        progressr_0.11.0       rstudioapi_0.14        ROCR_1.0-11           
# [25] tensor_1.5             DOSE_3.22.1            listenv_0.8.0          labeling_0.4.2         GenomeInfoDbData_1.2.8 polyclip_1.10-4        bit64_4.0.5            farver_2.1.1          
# [33] pheatmap_1.0.12        downloader_0.4         parallelly_1.32.1      vctrs_0.5.0            treeio_1.20.2          generics_0.1.3         R6_2.5.1               graphlayouts_0.8.3    
# [41] rsvd_1.0.5             DelayedArray_0.22.0    bitops_1.0-7           spatstat.utils_3.0-1   cachem_1.0.6           fgsea_1.22.0           gridGraphics_0.5-1     assertthat_0.2.1      
# [49] promises_1.2.0.1       scales_1.2.1           ggraph_2.1.0           enrichplot_1.16.2      rgeos_0.5-9            gtable_0.3.1           globals_0.16.1         goftest_1.2-3         
# [57] tidygraph_1.2.2        rlang_1.0.6            systemfonts_1.0.4      splines_4.2.1          lazyeval_0.2.2         spatstat.geom_3.0-3    BiocManager_1.30.19    reshape2_1.4.4        
# [65] abind_1.4-5            httpuv_1.6.6           qvalue_2.28.0          tools_4.2.1            ggplotify_0.1.0        ellipsis_0.3.2         spatstat.core_2.4-4    jquerylib_0.1.4       
# [73] RColorBrewer_1.1-3     proxy_0.4-27           ggridges_0.5.4         Rcpp_1.0.9             zlibbioc_1.42.0        purrr_0.3.5            RCurl_1.98-1.9         rpart_4.1.16          
# [81] deldir_1.0-6           pbapply_1.5-0          viridis_0.6.2          zoo_1.8-11             cluster_2.1.3          magrittr_2.0.3         data.table_1.14.4      scattermore_0.8       
# [89] DO.db_2.9              lmtest_0.9-40          RANN_2.6.1             packrat_0.8.1          fitdistrplus_1.1-8     mime_0.12              xtable_1.8-4           gridExtra_2.3         
# [97] compiler_4.2.1         tibble_3.1.8           KernSmooth_2.23-20     crayon_1.5.2           shadowtext_0.1.2       R.oo_1.25.0            minqa_1.2.5            htmltools_0.5.3       
# [105] ggfun_0.0.7            mgcv_1.8-40            later_1.3.0            tidyr_1.2.1            aplot_0.1.8            DBI_1.1.3              tweenr_2.0.2           MASS_7.3-57           
# [113] boot_1.3-28            leidenbase_0.1.12      Matrix_1.5-1           cli_3.4.1              R.methodsS3_1.8.2      parallel_4.2.1         igraph_1.3.5           pkgconfig_2.0.3       
# [121] terra_1.6-17           plotly_4.10.0          spatstat.sparse_3.0-0  ggtree_3.4.4           bslib_0.4.0            XVector_0.36.0         yulab.utils_0.0.5      stringr_1.4.1         
# [129] digest_0.6.30          sctransform_0.3.5      RcppAnnoy_0.0.20       spatstat.data_3.0-0    Biostrings_2.64.1      leiden_0.4.3           fastmatch_1.1-3        tidytree_0.4.1        
# [137] uwot_0.1.14            nloptr_2.0.3           lifecycle_1.0.3        nlme_3.1-157           jsonlite_1.8.3         viridisLite_0.4.1      fansi_1.0.3            pillar_1.8.1          
# [145] lattice_0.20-45        KEGGREST_1.36.3        fastmap_1.1.0          httr_1.4.4             survival_3.3-1         GO.db_3.15.0           remotes_2.4.2          glue_1.6.2            
# [153] zip_2.2.2              png_0.1-7              bit_4.0.4              sass_0.4.2             ggforce_0.4.1          stringi_1.7.8          blob_1.2.3             textshaping_0.3.6     
# [161] memoise_2.0.1          irlba_2.3.5.1          future.apply_1.9.1     ape_5.6-2   

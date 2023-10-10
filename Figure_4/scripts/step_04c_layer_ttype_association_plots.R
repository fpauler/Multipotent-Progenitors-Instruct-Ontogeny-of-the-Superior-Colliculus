# cell pairing analysis

library (dendextend)
library (openxlsx)
library (pheatmap)
library (pvclust)
library (dplyr)
library (ggplot2)
library (ggupset)
library (ggbeeswarm)

set.seed(2401)

# checked 6.10.2023

# define a base folder
base <- "./"

source(paste(base, "other_data/analysis_functions.R", sep=""))

total_cells_clone <- get_total_clone_size(base, "PatchSeq_sample_list_final.xlsx")

#####################################################################################
##
# read PatchSeq data from previous analysis and create a grouping variable for excitatory/inhibitory neurons
# add ID for subclones - Clone + color
##
#####################################################################################

PatchSeq_Seurat_SCT <- readRDS(file = paste(base, "RDS_files/PatchSeq_Seurat_SCT_cellTypeAssignment.RDS", sep="/") )
PatchSeq_Seurat_SCT@meta.data$Color <- toupper(PatchSeq_Seurat_SCT@meta.data$Color)
PatchSeq_Seurat_SCT@meta.data$subclone <- paste( PatchSeq_Seurat_SCT@meta.data$Clone, PatchSeq_Seurat_SCT@meta.data$Color, sep=".")

# broad grouping
PatchSeq_Seurat_SCT@meta.data$broad_class <- ""
PatchSeq_Seurat_SCT@meta.data$broad_class[ which( grepl(pattern = "INH", x = PatchSeq_Seurat_SCT@meta.data$predicted_ttype )) ] <- "INH"
PatchSeq_Seurat_SCT@meta.data$broad_class[ which( grepl(pattern = "GLU", x = PatchSeq_Seurat_SCT@meta.data$predicted_ttype )) ] <- "EXC"

clone_type <- "subclone"
cc_cor <- 0.1 
cc_score <- 0.5

meta_data <- extract_cells_from_meta_data( seurat_meta = PatchSeq_Seurat_SCT@meta.data, total_clone_size = total_cells_clone, 
                                           cc_cor = cc_cor, cc_score = cc_score, clone_coverage = 0.3, clone_type = clone_type)

# identify the column of the current clone type
clone_type_idx <- which(colnames(meta_data) == clone_type)

###############################################################################
##  
# layer spanning t-type linkage 
##
###############################################################################

meta_data$predicted_ttype.layer <- paste(meta_data$predicted_ttype, meta_data$Subregion, sep=".")

# use the z-Score to identify layer specific t-types (done in main script, just use a vector here)
ttype_layer_association <- c("MEGLU4.sSC", "MEGLU5.sSC", "MEGLU6.sSC", "MEINH10.sSC", "MEINH12.sSC", "MEINH5.sSC", "MEINH9.sSC",
                             "MEGLU1.dSC", "MEINH2.dSC", "MEINH7.dSC", "MEGLU2.PAG", "MEGLU3.PAG", "MEINH3.PAG")

# create a binary matrix for clone <-> t-type association
ttypeLayer_binary_matrix <- sapply(unique(meta_data[,clone_type_idx]), function (clone_id) {
  sapply(ttype_layer_association, function (ttype) {
    any( meta_data[,clone_type_idx] == clone_id & meta_data$predicted_ttype.layer == ttype)
  })
})

# convert the binary matrix to numeric - necessary for pheatmap
ttypeLayer_binary_matrix[which(ttypeLayer_binary_matrix)] <- 1
# remove clones with no cell with the correct layer association
ttypeLayer_binary_matrix <- ttypeLayer_binary_matrix[, which(!apply(ttypeLayer_binary_matrix, 2, max) == 0) ]


#prepare and save heatmap
pheatmap( mat = ttypeLayer_binary_matrix, cluster_rows = F, 
          clustering_distance_cols = "binary", clustering_method = "ward.D", silent = T, legend = F,
          color = colorRampPalette(colors = c("white", "black"))(10), 
          filename = paste(base, "/Supplement/FigureS7N.", clone_type, ".pdf", sep=""), width=22.6, height = 16.2 )

pvclust_output <- pvclust( t(ttypeLayer_binary_matrix), method.dist = "binary", method.hclust="ward.D", nboot = 10000 )

pdf( file = paste(base, "/Supplement/FigureS7O.", clone_type, ".pdf", sep="") )
plot(pvclust_output)
dev.off()

# check layer distribution of patched cells

## missed cells per layer
missed_cells <- read.xlsx(xlsxFile = paste(base, "PatchSeq_sample_list_final.xlsx", sep="/"), sheet = 2)

missed_cells$layer.clone <- paste(missed_cells$Clone, missed_cells$Subregion, sep=".")
missed_cells_df <- reshape2::melt(table(missed_cells$layer.clone))

meta_data$layer.subclone <-  paste(meta_data$subclone, meta_data$Subregion, sep=".")

total_cells_clone_layers <- get_total_clone_size_layers(base, "PatchSeq_sample_list_final.xlsx")

layer_dist_tot <- reshape2::melt( total_cells_clone_layers[["subclone"]][unique(meta_data$layer.subclone)] )
layer_dist_tot$Var1 <- as.character(layer_dist_tot$Var1)
layer_dist_tot$layer <- sapply( layer_dist_tot$Var1, function (x) strsplit(x = x, split = "\\.")[[1]][3] )

tot_layer_dist <- plyr::ddply(layer_dist_tot, "layer", summarise,
                                total = sum(value)
                    )

tot_layer_dist$rel <- tot_layer_dist$total / sum(tot_layer_dist$total)

layer_dist_patch <- reshape2::melt( table(meta_data$Subregion) )
layer_dist_patch$rel <- layer_dist_patch$value / sum(layer_dist_patch$value)
colnames(layer_dist_patch) <- c("layer", "total", "rel")


layer_dist_plot <- rbind(tot_layer_dist, layer_dist_patch)

layer_dist_plot$group <- c("total", "total", "total", "patched", "patched", "patched")
layer_dist_plot$layer <- factor(layer_dist_plot$layer, levels=c("sSC", "dSC", "PAG"))

# test if the layer distribution of patched cells is different than the expected proportion in each layer based on total cells
test.diff <- chisq.test(x = layer_dist_patch$total, p = tot_layer_dist$rel ) 

# Chi-squared test for given probabilities
# 
# data:  layer_dist_patch$total
# X-squared = 2.2664, df = 2, p-value = 0.322

ggplot(layer_dist_plot, aes(x=layer, y=rel, fill=group)) + 
  geom_bar(stat="identity", position="dodge") + ylab("") + xlab("") +
  theme_classic() + scale_fill_manual( values = c("patched"= "grey60", "total" = "grey20")) +
  ggtitle( paste("chi-square", test.diff$p.value))
ggsave( paste(base, "/Supplement/FigureS6H.pdf", sep=""))

write.csv(x = layer_dist_plot, file = paste(base, "/Supplement/FigureS6H.csv", sep=""))

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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_IE.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_IE.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_IE.UTF-8      
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_IE.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] SeuratObject_4.0.2 ggbeeswarm_0.6.0   ggupset_0.3.0      ggplot2_3.3.5      dplyr_1.0.7        pvclust_2.2-0      pheatmap_1.0.12    openxlsx_4.2.4     dendextend_1.15.1 
# 
# loaded via a namespace (and not attached):
#   [1] zip_2.2.0          Rcpp_1.0.7         vipor_0.4.5        RColorBrewer_1.1-2 pillar_1.6.2       compiler_4.0.3     plyr_1.8.6         viridis_0.6.1      tools_4.0.3        digest_0.6.27     
# [11] packrat_0.7.0      lattice_0.20-41    lifecycle_1.0.0    tibble_3.1.4       gtable_0.3.0       viridisLite_0.4.0  pkgconfig_2.0.3    rlang_0.4.11       Matrix_1.3-4       DBI_1.1.1         
# [21] beeswarm_0.4.0     gridExtra_2.3      withr_2.4.2        stringr_1.4.0      generics_0.1.0     vctrs_0.3.8        grid_4.0.3         tidyselect_1.1.1   glue_1.4.2         R6_2.5.1          
# [31] fansi_0.5.0        farver_2.1.0       reshape2_1.4.4     purrr_0.3.4        magrittr_2.0.1     scales_1.1.1       ellipsis_0.3.2     assertthat_0.2.1   colorspace_2.0-2   labeling_0.4.2    
# [41] utf8_1.2.2         stringi_1.7.4      munsell_0.5.0      crayon_1.4.1  

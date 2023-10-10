# third step Pten Neuron analysis

# checked 4.10.2023

library (Seurat)
library (GenBinomApps)
library (plyr)
library (ggplot2)

set.seed (2401)

#define base folder
base <- "./"

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
  
  # add one dummy cell to avoid missing clusters
  for (cluster in unique(all_rel_group$Var1)) {
    tmp <-   all_rel_group[which(all_rel_group$Var1 == cluster),]
    if (nrow(tmp) == 1) {
      if (any(tmp$split == "sparse")) {
        # ctrl is missing
        total <- unique(all_rel_group[which(all_rel_group$split == "ctrl"), "total"])
        all_rel_group <- rbind(all_rel_group, data.frame( Var1 = cluster, value = 1, rel = 1/total, split = "ctrl", total = total))
      } else {
        total <- unique(all_rel_group[which(all_rel_group$split == "sparse"), "total"])
        all_rel_group <- rbind(all_rel_group, data.frame( Var1 = cluster, value = 1, rel = 1/total, split = "sparse", total = total))
        
      }
    }
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
    } else if (x < 0.1) {
      return ("*")
    } else {
      return ("")
    }
  })
  
  all_rel_group <- merge( all_rel_group, pvalues[, c("cluster", "pvalue", "p.adjust", "sig")], by.x="Var1", by.y="cluster")
  
  return (all_rel_group)
}

#read the data from previous steps and combine
embryo_link_df <- readRDS( paste(base, "/RDS_files/Pten_embryonic_link.RDS", sep="") )
Pten.neuron.seurat <- readRDS( file = paste(base, "/RDS_files/Pten.neuron.seurat.RDS", sep="" ) )
Pten.neuron.seurat <- AddMetaData(object = Pten.neuron.seurat, metadata = embryo_link_df)

# read the embryonic data for final plotting
embryonic_ref <- readRDS( file = paste(base, "/RDS_files/LaManno_embryonic.combined.CCscored.RDS", sep=""))
embryonic_ref <- FindClusters(embryonic_ref, resolution = 0.2)

# sanity check that custering is identical to before
DimPlot( embryonic_ref, label = T)
ggsave ( paste(  base, "/QC/embryo_umap_Pten_integration_4.pdf", sep="" ))

rel_diff <- prep_rel_df( meta_data = Pten.neuron.seurat@meta.data, group = "embryo_traj", split = "genotype", reference="ctrl")
cluster2adult <- readRDS( file = paste(base, "/RDS_files/linktoadult.RDS", sep="") )
adult2cluster <- sapply( as.character(unique(cluster2adult$real_cluster_name)), function (x) paste(cluster2adult[which(cluster2adult$real_cluster_name == x), "id"], collapse=","))
rel_diff$clusterID <- adult2cluster[ as.character(rel_diff$Var1) ]
  
ggplot( rel_diff[which(!rel_diff$sig == ""),], aes( x=split, y=rel, fill = split)) + 
  geom_bar(stat="identity") + 
  geom_text(aes(y=0.5, label=sig)) +
  geom_errorbar( aes(ymax=Upper.limit, ymin=Lower.limit)) +
  scale_fill_manual( values = c(ctrl="grey60", sparse="grey30")) +
  facet_wrap(~Var1, nrow=1) + theme_classic()
ggsave ( paste(  base, "/Figures/Figure6V.pdf", sep="" ), width = 6 )

# write out raw analysis
rel_diff$Confidence.Interval <- unlist( rel_diff$Confidence.Interval )
rel_diff$alpha <- unlist( rel_diff$alpha )
write.csv( x = rel_diff, file = paste(base, "/Supplement/Figure6V.csv", sep=""))

rel_diff_df <- data.frame()
for (cluster in unique(rel_diff$Var1)) {
  tmp <- rel_diff[which(rel_diff$Var1 == cluster),]
  sparse_rel <- tmp[which(tmp$split == "sparse"), "rel"]
  ctrl_rel <- tmp[which(tmp$split == "ctrl"), "rel"]
  ratio <- log2( sparse_rel / ctrl_rel )
  rel_diff_df <- rbind( rel_diff_df,data.frame( cluster = cluster, ratio = ratio, sig = tmp$sig[1]))
}

rownames( rel_diff_df ) <- as.character( rel_diff_df$cluster)


rel_diff_df$sig_ratio <- rel_diff_df$ratio
rel_diff_df$sig_ratio[ which(!grepl(pattern =  "\\*", x = rel_diff_df$sig)) ] <- 0
rel_diff_df$sig_cluster <- paste( rel_diff_df$cluster, rel_diff_df$sig, sep="")

embryonic_ref@meta.data$Pten_ratio <- rel_diff_df[ as.character( embryonic_ref@meta.data$seurat_clusters ), "sig_ratio" ]
embryonic_ref@meta.data$Pten_ratio[which(is.na(embryonic_ref@meta.data$Pten_ratio))] <- 0

embryonic_ref@meta.data$Pten_ratio[which(embryonic_ref@meta.data$Pten_ratio < -1)] <- -1
umap_df <- as.data.frame( Embeddings(embryonic_ref, reduction = "umap") )
umap_df$cellID <- rownames( umap_df )
tmp_df <- embryonic_ref@meta.data[,c("Age", "seurat_clusters", "Pten_ratio")]
tmp_df$cellID <- rownames(tmp_df)
umap_df <- merge ( tmp_df, umap_df , by="cellID" )


umap_df$bin_ratio <- sapply( umap_df$Pten_ratio, function (x) {
  if (x > 0) {
    return (1)
  } else if (x < 0) {
    return (-1)
  }else {
    return(0)
  }
})
umap_df$bin_ratio[which(umap_df$seurat_clusters %in% c(0,6,11))] <- 2

umap_df$bin_ratio <- as.character(umap_df$bin_ratio)

# calculate centroids of each cluster for labeling
label_df <- ddply(umap_df, c("seurat_clusters"), summarise,
               UMAP_1 = median(UMAP_1),
               UMAP_2 = median(UMAP_2)
)
label_df$sig_cluster <- rel_diff_df[ as.character(label_df$seurat_clusters), "sig_cluster"]
label_df$sig_cluster[ which(is.na(label_df$sig_cluster)) ] <- as.character( label_df$seurat_clusters[ which(is.na(label_df$sig_cluster)) ] )

# up-regulation / up in Pten: 127/63/151 #7f3f97
# down-regulation /up in ctrl: 0/166/156 #00a69c
# do it binary
ggplot( ) + 
  geom_point( data = umap_df, aes(x=UMAP_1, y=UMAP_2, color=bin_ratio ), size = 0.5 )+ 
  scale_color_manual( values = c("-1" = "#00a69c", "0" = "grey60", "1" = "#7f3f97", "2" = "grey80")) +
  geom_text( data = label_df, aes( x=UMAP_1, y=UMAP_2, label=sig_cluster )) +
  theme_classic()

ggsave ( paste(  base, "/Figures/Figure6U.pdf", sep="" ), width=6)

ggplot( ) + 
  geom_point( data = umap_df, aes(x=UMAP_1, y=UMAP_2, color=bin_ratio ), size = 0.5 )+ 
  scale_color_manual( values = c("-1" = "#00a69c", "0" = "grey60", "1" = "#7f3f97", "2" = "grey80")) +
  theme_classic()

ggsave ( paste(  base, "/Figures/Figure6U_noLabel.png", sep="" ), width=6)

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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.3.6      plyr_1.8.7         GenBinomApps_1.2   sp_1.5-0           SeuratObject_4.1.2 Seurat_4.2.0      
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.16            colorspace_2.0-3      deldir_1.0-6          ellipsis_0.3.2        ggridges_0.5.4        rstudioapi_0.14       spatstat.data_3.0-0   leiden_0.4.3         
# [9] listenv_0.8.0         farver_2.1.1          ggrepel_0.9.1         fansi_1.0.3           codetools_0.2-18      splines_4.2.1         polyclip_1.10-4       jsonlite_1.8.3       
# [17] packrat_0.8.1         ica_1.0-3             cluster_2.1.3         png_0.1-7             rgeos_0.5-9           uwot_0.1.14           shiny_1.7.3           sctransform_0.3.5    
# [25] spatstat.sparse_3.0-0 compiler_4.2.1        httr_1.4.4            Matrix_1.5-1          fastmap_1.1.0         lazyeval_0.2.2        cli_3.4.1             later_1.3.0          
# [33] htmltools_0.5.3       tools_4.2.1           igraph_1.3.5          gtable_0.3.1          glue_1.6.2            RANN_2.6.1            reshape2_1.4.4        dplyr_1.0.10         
# [41] Rcpp_1.0.9            scattermore_0.8       vctrs_0.5.0           nlme_3.1-157          progressr_0.11.0      lmtest_0.9-40         spatstat.random_2.2-0 stringr_1.4.1        
# [49] globals_0.16.1        mime_0.12             miniUI_0.1.1.1        lifecycle_1.0.3       irlba_2.3.5.1         goftest_1.2-3         future_1.28.0         MASS_7.3-57          
# [57] zoo_1.8-11            scales_1.2.1          spatstat.core_2.4-4   ragg_1.2.4            promises_1.2.0.1      spatstat.utils_3.0-1  parallel_4.2.1        RColorBrewer_1.1-3   
# [65] reticulate_1.26       pbapply_1.5-0         gridExtra_2.3         rpart_4.1.16          stringi_1.7.8         rlang_1.0.6           pkgconfig_2.0.3       systemfonts_1.0.4    
# [73] matrixStats_0.62.0    lattice_0.20-45       ROCR_1.0-11           purrr_0.3.5           tensor_1.5            patchwork_1.1.2       htmlwidgets_1.5.4     labeling_0.4.2       
# [81] cowplot_1.1.1         tidyselect_1.2.0      parallelly_1.32.1     RcppAnnoy_0.0.20      magrittr_2.0.3        R6_2.5.1              generics_0.1.3        pillar_1.8.1         
# [89] withr_2.5.0           mgcv_1.8-40           fitdistrplus_1.1-8    survival_3.3-1        abind_1.4-5           tibble_3.1.8          future.apply_1.9.1    KernSmooth_2.23-20   
# [97] utf8_1.2.2            spatstat.geom_3.0-3   plotly_4.10.0         grid_4.2.1            data.table_1.14.4     digest_0.6.30         xtable_1.8-4          tidyr_1.2.1          
# [105] httpuv_1.6.6          textshaping_0.3.6     munsell_0.5.0         viridisLite_0.4.1

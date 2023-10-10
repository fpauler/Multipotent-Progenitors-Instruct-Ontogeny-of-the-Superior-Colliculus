# step4 - produces Figures - focuses on bulk analysis

# checked 6.10.2023

library (Seurat)
library (ggplot2)
library (ggbeeswarm)
library (openxlsx)
library (dplyr)
library (patchwork)
library (GenBinomApps)
library (cowplot)
library (pheatmap)
library (scales)
library (data.table)

# for parallelization of Seurat - speeds up marker detection
library (future)
# upsetplot
library (ggupset)
library (pvclust)

# define a base folder
base <- "./"
source(paste(base, "other_data/analysis_functions.R", sep=""))

set.seed(2401)

############
##
# ANALYSIS from here
##
############

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

#####################################################

# first prepare all possible combinations - identical to length(ncol(combn(16, 2)))
cell_type_combinations <- combn( unique(PatchSeq_Seurat_SCT@meta.data$predicted_ttype), 2)
cell_type_combinations <- sapply(1:ncol(cell_type_combinations), function (x) paste(sort(cell_type_combinations[,x]), collapse=".")) 
cell_type_combinations <- length(cell_type_combinations)

# enable parallelization of Seurat
plan("multiprocess", workers = 4)

#######################################################################################
##
# read reference Seurat Object and do DEG calculations of Reference cells
##
#######################################################################################

cti_Seurat<- readRDS( file = paste(base, "RDS_files/cti_Seurat_SCT.RDS", sep="/")  )

# DEG analysis between exc and inh neurons in reference 
Idents(cti_Seurat) <- "tax_4"
exc.inh.markers <- FindMarkers(cti_Seurat, ident.1 = "Di- and mesencephalon excitatory neurons", 
                               ident.2 = "Di- and mesencephalon inhibitory neurons", logfc.threshold = 0.2)
colnames(exc.inh.markers) <- sapply(X = colnames(exc.inh.markers), FUN = function (x) paste("ref", x, sep="."))
exc.inh.markers$gene <- rownames(exc.inh.markers)

# DEG analysis between all t-types in Reference
Idents(cti_Seurat) <- "tax_5"
tax5.cti_Seurat <- FindAllMarkers( cti_Seurat, logfc.threshold = 0.2)

# extract top 20 marker genes
top_markers_tax5_Reference <- tax5.cti_Seurat %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

# write out the DEG stats
write.csv(x = exc.inh.markers, file = paste(base, "/Supplement/Figure4C_Reference_DEG_stats.", clone_type, ".csv", sep="") )
write.csv(x = tax5.cti_Seurat, file = paste(base, "/Supplement/FigureS6IJ_Reference_DEG_stats.", clone_type, ".csv", sep="") )

# Initially we did all analyses both at the clonal and subclone level
# NOTE: this makes a difference to all analyses steps due to the clone_coverage cutoff!
# results are very similar, just there are much less clones than subclones
# also at E10 mainly symmetric division happen, so subclones also originate from RGPs

# for final version, to avoid confusion only use subclone!

  clone_type <-  "subclone"
  cc_cor <- 0.1 
  cc_score <- 0.5
  
  meta_data <- extract_cells_from_meta_data( seurat_meta = PatchSeq_Seurat_SCT@meta.data, total_clone_size = total_cells_clone, 
                                        cc_cor = cc_cor, cc_score = cc_score, clone_coverage = 0.3, clone_type = clone_type)
  
  ######################################################################
  ##
  # prepare some basic statistics of Patch-Seq cells - Supplemental Figure
  ##
  ######################################################################
  
  basic_stat_1 <- ggplot( meta_data, aes(x="number of genes", y = nFeature_RNA / 1000)) + 
    geom_boxplot(outlier.shape = NA, color="grey40") + geom_beeswarm(size=0.25, cex = 1.5) + 
    ylim(0, max(meta_data$nFeature_RNA / 1000)) +
    theme_classic() + ggtitle(paste("number of genes x 1000\nmedian:\n", median(meta_data$nFeature_RNA))) + 
    ylab("") + xlab("")
  
  basic_stat_2 <- ggplot( meta_data, aes(x="number of reads", y = STAR.unique.aligned / 1000000)) + 
    geom_boxplot(outlier.shape = NA, color="grey40") + geom_beeswarm(size=0.25, cex = 1.5) + 
    ylim(0, max(meta_data$STAR.unique.aligned / 1000000)) +
    scale_y_continuous(breaks = seq(0,3,1)) +
    theme_classic() + ggtitle(paste("unique M reads \nmedian:\n", median(meta_data$STAR.unique.aligned))) + 
    ylab("") + xlab("")
  
  basic_stat_3 <- ggplot( meta_data, aes(x="alignment rate", y = STAR.perc)) + 
    geom_boxplot(outlier.shape = NA, color="grey40") + geom_beeswarm(size=0.25, cex = 1.5) + ylim(0,1) + 
    theme_classic() + scale_y_continuous(breaks = seq(0,1,0.2)) +
    ggtitle(paste("alignment rate (unique reads)\nmedian:\n", round(median(meta_data$STAR.perc), 2))) + 
    ylab("") + xlab("")
  
  hist_data <- hist( meta_data$boot_score, plot = F, breaks = seq(0,1,0.05))
  hist_datadf<- data.frame(x = hist_data$breaks[2:length(hist_data$breaks)], y = hist_data$counts)
  
  boot_score_plot <- ggplot(hist_datadf, aes(x=x, y=y)) + 
    geom_bar(stat="identity", position = "dodge") + 
    scale_y_continuous(breaks = seq(0,140,50)) + scale_x_continuous(breaks = seq(0,1,0.2)) + 
    theme_classic() + ggtitle("bootstrap score\n#cells") + ylab("") + xlab("")
  
  hist_data <- hist(meta_data$top_cor, plot = F, breaks = seq(0,1,0.05))
  hist_datadf<- data.frame(x = hist_data$breaks[2:length(hist_data$breaks)], y = hist_data$counts)
  
  pearson_cor_plot <- ggplot(hist_datadf, aes(x=x, y=y)) + geom_bar(stat="identity", position = "dodge") + 
    scale_x_continuous(breaks = seq(0,1,0.2)) + 
    theme_classic() + ggtitle("Pearson correlation\n# cells") + ylab("") + xlab("")
  
  hist_data <- hist(meta_data$NMS_Neurons, plot = F, breaks = seq(0,6,0.1))
  hist_datadf<- data.frame(x = hist_data$breaks[2:length(hist_data$breaks)], y = hist_data$counts, group = "neuron")
  
  hist_data <- hist(apply(meta_data[,c("NMS_Oligo", "NMS_Vascular", "NMS_Immune", "NMS_Astro")], 1, max), plot = F, breaks = seq(0,6,0.1))
  hist_datadf<- rbind(hist_datadf, data.frame(x = hist_data$breaks[2:length(hist_data$breaks)], y = hist_data$counts, group = "cont"))
  
  NMS_plot <- ggplot(hist_datadf, aes(x=x, y=y, fill=group)) + geom_bar(stat="identity", position = "dodge", width = 0.1) + 
    theme_classic() + ggtitle("NMS scores\n# cells") + ylab("") + xlab("")
  
  plot_basic_stats <- (basic_stat_1 + basic_stat_2 + basic_stat_3) / (NMS_plot + boot_score_plot + pearson_cor_plot)
  ggsave(paste(base, "/Supplement/FigureS6BtoG.", clone_type, ".pdf", sep=""), plot = plot_basic_stats, width = 6, height=6)
  
  # save some stats for the text
  cor_stats_df <- data.frame()
  for ( cor_stats in c("nFeature_RNA", "STAR.unique.aligned", "STAR.perc", "boot_score", "top_cor")) {
    mean_data <- mean(meta_data[,cor_stats])
    sd_data <- sd(meta_data[,cor_stats])
    sem <- sd_data/sqrt(nrow(meta_data))
    tmp_df <- data.frame(stat = c("mean", "sd", "sem"), value = c(mean_data, sd_data, sem))
    colnames(tmp_df) <- c("stat", cor_stats)
    if (nrow(cor_stats_df) == 0) {
      cor_stats_df <- tmp_df
    } else {
      cor_stats_df <- merge(cor_stats_df, tmp_df, by="stat")
    }
  }
  write.csv(x = cor_stats_df, file = paste(base, "/Supplement/basic_stats.", clone_type, ".csv", sep=""))
  
  
  # identify the column of the current clone type
  clone_type_idx <- which(colnames(meta_data) == clone_type)
  
  ############################################
  ##
  # UMAP of PatchSeq cells on the Zeisel data
  ##
  ############################################
  
  # this comes from the third analysis step and is mainly a re-plot of the already integrated UMAP  
  UMAP_df <- readRDS(file = paste(base, "/RDS_files/integrated_umap_df.RDS", sep=""))
  
  # only show cells that are actually in the analysis
  UMAP_df <- UMAP_df[which(UMAP_df$Genotype == "reference" | rownames(UMAP_df) %in% rownames(meta_data)),]
  
  # color palette
  # color palette for UMAP
  hex_codes2 <- hue_pal()(length(unique(UMAP_df$tax5))) 
  #names(hex_codes2) <- sort(unique(UMAP_df$tax5))
  names(hex_codes2) <- c("MEGLU1", "MEGLU2", "MEGLU3", "MEGLU4", "MEGLU5", "MEGLU6", 
                          "MEINH2", "MEINH3", "MEINH5", "MEINH6", "MEINH7",
                         "MEINH8", "MEINH9", "MEINH10", "MEINH11", "MEINH12", "PatchSeq")
  
  hex_codes2["PatchSeq"] <- "black"
  
  # prepare a label for the clusters
  centroid_label <- data.frame()
  for (ttype in unique(UMAP_df$tax5) ) {
    tmp_df <- UMAP_df[which(UMAP_df$tax5 == ttype & UMAP_df$Genotype == "reference"),]
 
    #extract only cells under investigation - minus the outliers as defined above
    core_tmp_df <- tmp_df[which(!tmp_df$tax5 == 'out'),]

    # label df
    centroid_label <- rbind(centroid_label, 
                            data.frame(x = mean(core_tmp_df$UMAP_1), y = mean(core_tmp_df$UMAP_2), label = ttype, stringsAsFactors = F))
  }
  
  UMAP_main <- ggplot() + 
    geom_point( data = UMAP_df, aes(x=UMAP_1, y=UMAP_2, color=tax5, size = Genotype) ) +
    scale_color_manual(values = hex_codes2) +
    geom_text(data = centroid_label, aes(x=x, y=y, label=label)) +
    scale_size_manual(values = c("reference" = 1, "M11GTTG-Fz10-CreER" = 2.5, "M11GTTG-Sox2-CreER" = 2.5)) +
    theme_classic() + theme(legend.position = "none")
  
  ggsave(paste(base, "/Figures/Figure4B_", clone_type, "_tax5.pdf", sep=""), width=12, height=12, plot = UMAP_main)
  ggsave(paste(base, "/Figures/Figure4B_", clone_type, "_tax5.eps", sep=""), width=12, height=12, plot = UMAP_main)
  
  # Warning message:
  # Removed 1 rows containing missing values (geom_text). 
  # this comes from PatchSeq cells that are NA for the label
  
  #####################################################################
  ##
  # DEG analysis PatchSeq
  # calculated here as cells might be excluded based on the filtering
  ##
  #####################################################################
  
  Idents(PatchSeq_Seurat_SCT) <- "broad_class"
  exc.inh.PatchSeq <- FindMarkers( subset(PatchSeq_Seurat_SCT, cells = rownames(meta_data)), ident.1 = "EXC", ident.2 = "INH", logfc.threshold = 0.2)
  
  Idents(PatchSeq_Seurat_SCT) <- "predicted_ttype"
  tax5.PatchSeq <- FindAllMarkers( subset(PatchSeq_Seurat_SCT, cells = rownames(meta_data)), logfc.threshold = 0.2)
  
  # write out the DEG stats
  write.csv(x = exc.inh.PatchSeq, file = paste(base, "/Supplement/Figure4C_PatchSeq_DEG_stats.", clone_type, ".csv", sep="") )
  write.csv(x = tax5.PatchSeq, file = paste(base, "/Supplement/FigureS6IJ_PatchSeq_DEG_stats.", clone_type, ".csv", sep="") )
  
  # prepare a nice heatmap - do a DEG score based on the log10 of the adjusted p-value
  colnames(exc.inh.PatchSeq) <- sapply(X = colnames(exc.inh.PatchSeq), FUN = function (x) paste("PatchSeq", x, sep="."))
  exc.inh.PatchSeq$gene <- rownames(exc.inh.PatchSeq)
  
  # join the DEG analysis of Reference and PatchSeq - only genes in both data sets will be reported
  joined_deg <- merge(exc.inh.markers, exc.inh.PatchSeq, by="gene")
  joined_deg <- joined_deg[which(joined_deg$PatchSeq.p_val_adj < 0.01),]
  
  joined_deg$PatchSeq_score <- log10(joined_deg$PatchSeq.p_val_adj)
  # cut the score at -5 for nicer heatmap
  joined_deg$PatchSeq_score[which(joined_deg$PatchSeq_score < -5)] <- -5
  joined_deg$PatchSeq_score[which(joined_deg$PatchSeq.avg_log2FC > 0)] <- joined_deg$PatchSeq_score[which(joined_deg$PatchSeq.avg_log2FC > 0)] * -1
  
  joined_deg$ref_score <- log10(joined_deg$ref.p_val_adj)
  joined_deg$ref_score[which(is.infinite(joined_deg$ref_score))] <- -5
  joined_deg$ref_score[which(joined_deg$ref_score < -5)] <- -5
  joined_deg$ref_score[which(joined_deg$ref.avg_log2FC > 0)] <- joined_deg$ref_score[which(joined_deg$ref.avg_log2FC > 0)] * -1
  
  deg_heatmap_ggplot <- reshape2::melt(joined_deg[,c("gene", "ref_score", "PatchSeq_score" )])
  deg_heatmap_ggplot <- deg_heatmap_ggplot[order(deg_heatmap_ggplot$value, decreasing = T),]
  
  gene_order <- deg_heatmap_ggplot[which(deg_heatmap_ggplot$variable == "PatchSeq_score"), ]
  
  deg_heatmap_ggplot$gene <- factor(deg_heatmap_ggplot$gene, levels = gene_order$gene )
  
  # heatmap of exc/inh marker genes in reference and real data
  main_a <- ggplot(deg_heatmap_ggplot, aes(x=variable, y=gene, fill=value)) + geom_tile() + 
    scale_fill_gradient2(low = "cornflowerblue", mid = "white", high = "goldenrod") + xlab("") + ylab("") +
    ggtitle("DEG score") + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(paste(base, "/Figures/Figure4C_", clone_type, ".pdf", sep=""), plot = main_a, width=6, height=6)
  
  ###################################
  ##
  # analyse EXC and INH cell types
  ##
  ###################################
  
  # count exc and inh cells per genotype
  exc_inh_genotype_df <- plyr::ddply(meta_data, c("Genotype"), summarise,
              exc = length( which( broad_class == "EXC" ) ),
              inh = length( which( broad_class == "INH" ) ),
              genotype = unique(Genotype)
  )
  #             Genotype exc inh           genotype
  # 1 M11GTTG-Fz10-CreER  50  65 M11GTTG-Fz10-CreER
  # 2 M11GTTG-Sox2-CreER  67  71 M11GTTG-Sox2-CreER
  
  # plot to show abundance of exc and inh cells - independent of genotype
  main_b <- ggplot( reshape2::melt(exc_inh_genotype_df), aes(x = variable, y=value, fill = variable)) + 
    geom_bar(stat = "identity") +
    scale_fill_manual( values = c( "inh" = "cornflowerblue", "exc" = "goldenrod") ) +
    theme_classic() + ylab("") + xlab("") + ggtitle("# neurons")
  
  # Supplement: plot to show abundance of exc and inh cells separated by cre-driver
  main_b_supplement <- ggplot( reshape2::melt(exc_inh_genotype_df), aes(x = variable, y=value, fill = variable)) + 
    geom_bar(stat = "identity") +
    scale_fill_manual( values = c( "inh" = "cornflowerblue", "exc" = "goldenrod") ) +
    theme_classic() + ylab("") + xlab("") + ggtitle("# neurons") + facet_grid(~Genotype)
  
  # write out raw data
  write.csv( x = exc_inh_genotype_df, file = paste( base, "/Supplement/Figure4D_S6K.csv", sep=""))
  
  # count exc and inh cells per layer
  exc_inh_subregion_df <- plyr::ddply(meta_data, c( "Subregion" ), summarise,
                                     exc = length( which( broad_class == "EXC" ) ),
                                     inh = length( which( broad_class == "INH" ) )
  )
  
  
  
  # calculate Chi-square goodness-of-fit test to test for equal distribution among layers
  expected <- c(1/3, 1/3, 1/3) #must add up to 1
  
  # EXC:
  chisq.test(x=exc_inh_subregion_df$exc, p=expected, simulate.p.value = T, )
  
  # Chi-squared test for given probabilities with simulated p-value (based on 2000 replicates)
  # 
  # data:  exc_inh_subregion_df$exc
  # X-squared = 3.2308, df = NA, p-value = 0.2049  
  
  #INH:
  chisq.test(x=exc_inh_subregion_df$inh, p=expected, simulate.p.value = T)
  
  # Chi-squared test for given probabilities with simulated p-value (based on 2000 replicates)
  # 
  # data:  exc_inh_subregion_df$inh
  # X-squared = 43.868, df = NA, p-value = 0.0004998
  
  exc_inh_subregion_df$exc_perc <- exc_inh_subregion_df$exc / sum( exc_inh_subregion_df$exc )
  exc_inh_subregion_df$inh_perc <- exc_inh_subregion_df$inh / sum( exc_inh_subregion_df$inh )
  
  exc_inh_subregion_df$exc_perc <- exc_inh_subregion_df$exc_perc * -1
  
  exc_inh_subregion_df <- reshape2::melt( exc_inh_subregion_df[,c("Subregion", "exc_perc", "inh_perc")] )
  
  exc_inh_subregion_df$Subregion <- factor( exc_inh_subregion_df$Subregion, levels = c("PAG", "dSC", "sSC"))
  
  # write out the raw plotting data
  write.csv( x = exc_inh_subregion_df, file = paste( base, "/Supplement/Figure4E.csv", sep=""))
  
  # main figure - show all cells
  main_c <- ggplot(exc_inh_subregion_df, aes(y = Subregion, x=value, fill = variable)) + 
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual( values = c( "inh_perc" = "cornflowerblue", "exc_perc" = "goldenrod") ) +
    theme_classic() + ylab("") + xlab("") + xlim(-0.65, 0.65) +
    ggtitle("relative abundance in layers")
  
  # Supplement - separate data by genotype
  exc_inh_subregion_df <- plyr::ddply(meta_data, c( "Subregion", "Genotype"), summarise,
                                      exc = length( which( broad_class == "EXC" ) ),
                                      inh = length( which( broad_class == "INH" ) )
  )
  
  # Chi-square goodness-of-fit test
  sapply( c("M11GTTG-Fz10-CreER", "M11GTTG-Sox2-CreER"), function (cre){
    sapply ( c("exc", "inh"), function (type) {
      tmp <- exc_inh_subregion_df[which(exc_inh_subregion_df$Genotype == cre),]
      p_calc <- chisq.test(x=tmp[,type], p=expected, simulate.p.value = T)
      return(p_calc$p.value)
    })
  })
  
  # M11GTTG-Fz10-CreER M11GTTG-Sox2-CreER
  # exc       0.1249375312       0.3238380810
  # inh       0.0004997501       0.0004997501  
  
  df2plot <- data.frame()
  for (genotype in unique (exc_inh_subregion_df$Genotype)) {
    tmp <- exc_inh_subregion_df[which(exc_inh_subregion_df$Genotype == genotype),]
    perc_exc <- tmp$exc / sum(tmp$exc)
    perc_inh <- tmp$inh / sum(tmp$inh)
    df2plot <- rbind(df2plot, data.frame(genotype = genotype, subregion = tmp$Subregion, type=c("exc"), perc=perc_exc*-1))
    df2plot <- rbind(df2plot, data.frame(genotype = genotype, subregion = tmp$Subregion, type=c("inh"), perc=perc_inh))
  }
  
  df2plot$subregion <- factor( df2plot$subregion, levels = c("PAG", "dSC", "sSC"))
  
  main_c_supplement <- ggplot(df2plot, aes(y = subregion, x=perc, fill = type)) + 
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual( values = c( "inh" = "cornflowerblue", "exc" = "goldenrod") ) +
    theme_classic() + ylab("") + xlab("") + xlim(-0.65, 0.65) +
    ggtitle("relative abundance in layers") + facet_grid(~genotype)

  #write out the raw plotting data
  write.csv( x = df2plot, file = paste( base, "/Supplement/FigureS6L.csv", sep=""))
  
  ###########################
  ##
  # DEG heatmap for t-types
  ##
  ###########################
  
  # create a convenient t-type order
  ttype_list <- c( "MEGLU6", "MEGLU2", "MEGLU3", "MEGLU1", "MEGLU5", "MEGLU4",
                   "MEINH8", "MEINH3", "MEINH2", "MEINH7", "MEINH5", "MEINH12", "MEINH11", "MEINH9", "MEINH10", "MEINH6")
  
  # initially S6J was a main figure and contained PatchSeq and Reference data
  # eventually this was changed - I left the code and only adapted figure numbers
 
  for (figure_type in c("Supplement", "main")) {
    
    # create a convenient gene order - NOTE: this is a hand picked list of genes created by investigating DEG lists in Patch-Seq and Reference
    if (figure_type == "main") {
      gene_list <- c("Sln", "Ngb", "Cartpt", "Nrn1", "Calb1", "Nrgn", "Sst", "Vgf", "Syt2", "Pax7", "Tnnt1", "Pcp4", "C1ql3", "Cbln4", "Shisa8", "Gm2694")
    } else {
      gene_list <- unique(top_markers_tax5_Reference$gene)[ which( unique(top_markers_tax5_Reference$gene) %in% tax5.PatchSeq$gene ) ]
    }
    
    # match the gene list to the smaller PatchSeq DEG analysis   
    gene_list <- gene_list[which(gene_list %in% tax5.PatchSeq$gene)]
    
    df2plot <- tax5.PatchSeq[which(tax5.PatchSeq$gene %in% gene_list),]
    df2plot$score <- log2( df2plot$p_val )
    df2plot$score[which(df2plot$avg_log2FC > 0)] <- df2plot$score[which(df2plot$avg_log2FC > 0)] * -1
    df2plot$score[which(df2plot$score > 30)] <- 30
    
    df2plot$gene <- factor( df2plot$gene, levels = gene_list)
    df2plot$cluster <- factor( df2plot$cluster, levels = ttype_list)
    df2plot <- df2plot[which(df2plot$score > 0),]
    
    main_d <- ggplot(df2plot, aes(x = cluster, y = gene, fill = score)) + 
      geom_tile() + 
      scale_fill_gradient( low = "grey80", high = "black") +
      theme_classic() + ylab("") + xlab("") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.text.y = element_text(size=5)) +
      ggtitle("t-type DEG score Patch-Seq")
    
    # do a comparative plot for reference  -  only positive scores!
    df2plot <- tax5.cti_Seurat[which(tax5.cti_Seurat$gene %in% gene_list),]
    df2plot$score <- log2( df2plot$p_val )
    df2plot$score[which(df2plot$avg_log2FC > 0)] <- df2plot$score[which(df2plot$avg_log2FC > 0)] * -1
    df2plot$score[which(df2plot$score > 200)] <- 200
    df2plot <- df2plot[which(df2plot$score > 0),]
    
    df2plot$gene <- factor( df2plot$gene, levels = gene_list)
    df2plot$cluster <- factor( df2plot$cluster, levels = ttype_list)
    
    ttype_heat_plot_Reference <- ggplot(df2plot, aes(x = cluster, y = gene, fill = score)) + 
      geom_tile() + 
      scale_fill_gradient( low = "grey80", high = "black") +
      theme_classic() + ylab("") + xlab("") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.text.y=element_blank()) +
      ggtitle("t-type DEG score reference")
    
    comb_heatmap <- main_d + ttype_heat_plot_Reference
    
    if (figure_type == "main") {
      ggsave(paste(base, "/Supplement/FigureS6J_", clone_type, ".eps", sep=""), width=12, height=8)  
    } else {
      ggsave(paste(base, "/Supplement/FigureS6I_", clone_type, ".eps", sep=""), width=12, height=8)  
    }
    
  }

  ##############################################
  ##
  # relative abundance of t-types in the data
  ##
  ##############################################
  
  cdata_PatchSeq <- plyr::ddply(meta_data, c( "predicted_ttype" ), summarise,
                                N    = length(predicted_ttype)
  )
  
  cdata_PatchSeq$relN <- cdata_PatchSeq$N / sum(cdata_PatchSeq$N)
  
  # add conficende interval as in Cadwell et al., Elife 2020 - not used at the moment
  conf_int <- t(sapply (1:nrow(cdata_PatchSeq), function (x) {
    
    tmp <- clopper.pearson.ci(cdata_PatchSeq$N[x], sum(cdata_PatchSeq$N), alpha = 0.05, CI = "two.sided")
    
  }))
  
  cdata_PatchSeq <- cbind(cdata_PatchSeq, conf_int)
  
  cdata_PatchSeq$Lower.limit <- as.numeric(cdata_PatchSeq$Lower.limit)
  cdata_PatchSeq$Upper.limit <- as.numeric(cdata_PatchSeq$Upper.limit)
  
  tax5_order <- cdata_PatchSeq$predicted_ttype[order(cdata_PatchSeq$relN)]
  
  cdata_PatchSeq$predicted_ttype <- factor(cdata_PatchSeq$predicted_ttype, levels = tax5_order)
  
  # broad class tag - not used at the moment
  cdata_PatchSeq$broad_class <- ""
  cdata_PatchSeq$broad_class[which(grepl(pattern = "INH", x = cdata_PatchSeq$predicted_ttype))] <- "INH"
  cdata_PatchSeq$broad_class[which(grepl(pattern = "GLU", x = cdata_PatchSeq$predicted_ttype))] <- "EXC"

  cdata_PatchSeq$Confidence.Interval <- unlist( cdata_PatchSeq$Confidence.Interval)
  cdata_PatchSeq$alpha <- unlist( cdata_PatchSeq$alpha)
  
  #write out the raw plotting data
  write.csv( x = cdata_PatchSeq, file = paste( base, "/Supplement/FigureS6M.csv", sep=""))
  
  # plot showing relative abundance of t-types
  main_e <- ggplot( cdata_PatchSeq, aes(y=predicted_ttype, x=relN )) + geom_bar(stat="identity", position = "dodge") +
    theme_classic() + xlab("") + ylab("") +
    ggtitle ("relative\nt-type abundance")
  
  # Supplement compares the t-type abundance between genotypes
  cdata_PatchSeq <- plyr::ddply(meta_data, c("Genotype", "predicted_ttype"), summarise,
                                N    = length(Genotype)
  )
  
  genotype_length <- table(meta_data$Genotype)
  
  # add conficende interval as in Cadwell et al., Elife 2020
  conf_int <- t(sapply (1:nrow(cdata_PatchSeq), function (x) {
    
    tmp <- clopper.pearson.ci(cdata_PatchSeq$N[x], genotype_length[cdata_PatchSeq$Genotype[x]], alpha = 0.05, CI = "two.sided")
    
  }))
  
  cdata_PatchSeq <- cbind(cdata_PatchSeq, conf_int)
  
  cdata_PatchSeq$Lower.limit <- as.numeric(cdata_PatchSeq$Lower.limit)
  cdata_PatchSeq$Upper.limit <- as.numeric(cdata_PatchSeq$Upper.limit)
  
  cdata_PatchSeq$predicted_ttype <- factor(cdata_PatchSeq$predicted_ttype, levels = tax5_order)

  cdata_PatchSeq$relN <- cdata_PatchSeq$N / genotype_length[cdata_PatchSeq$Genotype]
  cdata_PatchSeq$Genotype <- gsub(pattern = "M11GTTG-", replacement = "", x = cdata_PatchSeq$Genotype)
  cdata_PatchSeq$Genotype <- gsub(pattern = "-CreER", replacement = "", x = cdata_PatchSeq$Genotype)
  
  # test significance of differences
  gen1 <- "Fz10"
  gen2 <- "Sox2"
  
  data1 <- cdata_PatchSeq[ which(cdata_PatchSeq$Genotype == gen1), ]
  rownames(data1) <- data1$predicted_ttype
  data2 <- cdata_PatchSeq[ which(cdata_PatchSeq$Genotype == gen2), ]
  rownames(data2) <- data2$predicted_ttype
  
  # determine total number of cells per genotype
  genotype_length <- c( gen1 = sum(data1$N), gen2 = sum(data2$N) )
  names(genotype_length) <- c(gen1, gen2)
  
  # merge the dataframes and add zeros for undetected cell types
  merged_data <- merge(data1[,c("predicted_ttype", "N")], data2[,c("predicted_ttype", "N")], by="predicted_ttype", all=T)
  merged_data[is.na(merged_data)] <- 0
  
  # restructure data and do the statistics
  count_table <- as.table(rbind( merged_data$N.x , merged_data$N.y))
  rownames(count_table) <- c(gen1, gen2)
  colnames(count_table) <- merged_data$predicted_ttype
  
  cnames <- colnames(count_table)
  
  tmp_pvalues <- sapply( cnames, function (x) {
    tmp_table <- as.table( rbind( c(count_table[gen1, x], genotype_length[gen1] - count_table[gen1, x]),
                                  c(count_table[gen2, x], genotype_length[gen2] - count_table[gen2, x]))
    )
    
    chisq_out <- chisq.test( tmp_table, correct = F, rescale.p = F, simulate.p.value = T )
    return(chisq_out$p.value)
  })
  
  # write out p-values
  tmp_p_out <- p.adjust( p = tmp_pvalues, method = "bonferroni")
  tmp_p_out <- cbind( as.data.frame(tmp_pvalues), as.data.frame(tmp_p_out) )
  tmp_p_out$ttype <- rownames(tmp_p_out)
  colnames(tmp_p_out) <- c("raw_p", "padj", "ttype")
  write.csv(x = tmp_p_out, file = paste(base, "/Supplement/FigureS6N_pvalue.", clone_type,".csv", sep=""))
  
  # write out raw plot data
  cdata_PatchSeq$Confidence.Interval <- unlist ( cdata_PatchSeq$Confidence.Interval )
  cdata_PatchSeq$alpha <- unlist ( cdata_PatchSeq$alpha )
  
  write.csv(x = cdata_PatchSeq, file = paste(base, "/Supplement/FigureS6N.", clone_type,".csv", sep=""))
  
  
  # plot showing relative abundance of t-types
  main_e_supplement <- ggplot( cdata_PatchSeq, aes(x=Genotype, y=relN )) + geom_point() +
    theme_classic() + xlab("") + ylab("") + geom_errorbar(aes(ymin=Upper.limit, ymax=Lower.limit)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), strip.text.x.top = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
    ggtitle ("relative\nt-type abundance") + facet_grid(~predicted_ttype)
  
  #####################################
  ##
  # cell type distribution per layer
  ##
  #####################################
  
  ttype_Subregion <- plyr::ddply(meta_data, c("predicted_ttype"), summarise,
                                     sSC    = length(which(Subregion == "sSC")),
                                     dSC    = length(which(Subregion == "dSC")),
                                     PAG    = length(which(Subregion == "PAG"))
  )
  
  # randomize the layer information 1000 times to identify significance of layer association
  ttype_rand_Subregion <- lapply (1:1000, function (x) {
    meta_rand <- meta_data
    
    # randomize the subregion - random association of ttype with subregion
    meta_rand$Subregion <- sample(meta_rand$Subregion)
    
    tmp_ttype_count <- plyr::ddply(meta_rand, c("predicted_ttype"), summarise,
                                            sSC    = length(which(Subregion == "sSC")),
                                            dSC    = length(which(Subregion == "dSC")),
                                            PAG    = length(which(Subregion == "PAG"))
    )
    
    rownames(tmp_ttype_count) <- tmp_ttype_count$predicted_ttype
    tmp_ttype_count <- tmp_ttype_count[,2:4]
    return(tmp_ttype_count)
  })
  
  # calculate mean and SD from these randomizations
  Subregion_dist_SD <- sd_from_mat_list(ttype_rand_Subregion)
  colnames(Subregion_dist_SD) <- colnames(ttype_rand_Subregion[[1]])
  rownames(Subregion_dist_SD) <- rownames(ttype_rand_Subregion[[1]])
  
  Subregion_dist_mean <- Reduce('+', ttype_rand_Subregion) / length(ttype_rand_Subregion)

  rownames(ttype_Subregion) <- ttype_Subregion$predicted_ttype
  ttype_Subregion <- ttype_Subregion[,c(2:4)]
  
  # sanity checks
  identical(rownames(ttype_Subregion), rownames(Subregion_dist_SD))
  identical(rownames(ttype_Subregion), rownames(Subregion_dist_mean))
  
  identical(colnames(ttype_Subregion), colnames(Subregion_dist_SD))
  identical(colnames(ttype_Subregion), colnames(Subregion_dist_mean))
  
  # z-score calculation
  Subregion_dist_zScore <- ttype_Subregion - Subregion_dist_mean
  Subregion_dist_zScore <- Subregion_dist_zScore / Subregion_dist_SD
  
  # get max z-score for each ttype
  max_z_score <- sapply( rownames(Subregion_dist_zScore), function (x){
    tmp <- reshape2::melt(as.matrix(Subregion_dist_zScore))
    tmp <- tmp[which(tmp$Var1 == x),]
    return(tmp$value[which(tmp$value == max(tmp$value))])
  })
  
  # give a symbol for the signifiance based on z-score
  # probably not so useful
  sig_symbol <- sapply(max_z_score, function (x) {
    if (x > 1.3) {
      if (x > 1.96) {
        if (x > 3.58) {
          return ("***")
        } else {
          return ("**") 
        } 
      } else {
        return ("*")
      } 
    } else {
      return ("ns")
    }
  })
  
  # calculate a p-value from the z-score
  p_value_from_Z <- sapply(max_z_score, function (x) { 
    pvalue <- pnorm(q = x, lower.tail=FALSE) 
    pvalue <- formatC(pvalue, format = "e", digits = 2)
  })
  
  # cut the z-Score for better visualisation
  Subregion_dist_zScore[Subregion_dist_zScore < -2] <- -2
  Subregion_dist_zScore[Subregion_dist_zScore > 2] <- 2
  
  # create a heatmap - mainly to ge the nice clustering and a quick look
  subregion_zScore_plot <- pheatmap(Subregion_dist_zScore, scale = "none", silent = T)
  row_order <- subregion_zScore_plot$tree_row$labels[subregion_zScore_plot$tree_row$order]
  
  # arbitrary row order for nice display
  row_order_arb <- c( "MEINH6", "MEINH11", "MEINH8", "MEGLU6", "MEINH12", "MEINH5", "MEGLU4", "MEGLU5", "MEINH10", "MEINH9",
                      "MEINH2", "MEINH7", "MEGLU1", "MEINH3", "MEGLU2", "MEGLU3") 
  Subregion_dist_zScore <- reshape2::melt( as.matrix(Subregion_dist_zScore) )
  Subregion_dist_zScore$Var1 <- factor(Subregion_dist_zScore$Var1, levels = row_order_arb)
  Subregion_dist_zScore$label <- ""
  # add significance symbol
  Subregion_dist_zScore <- rbind( Subregion_dist_zScore, data.frame(Var1 = row_order_arb, Var2="sig", value = -2, label=sig_symbol[row_order_arb]) )
  Subregion_dist_zScore <- rbind( Subregion_dist_zScore, data.frame(Var1 = row_order_arb, Var2="p", value = -2, label=p_value_from_Z[row_order_arb]) )
  
  Subregion_dist_zScore$Var2 <- factor(Subregion_dist_zScore$Var2, levels = c("PAG", "dSC", "sSC", "sig", "p"))
    
  # write out raw plotting data
  write.csv(x = Subregion_dist_zScore, file = paste(base, "/Supplement/FigureS6O.", clone_type,".csv", sep=""))
  # write out raw z-score
  write.csv(x = as.data.frame( max_z_score ), file = paste(base, "/Supplement/FigureS6O_zScore.", clone_type,".csv", sep=""), row.names = T)
  
  main_f <- ggplot( Subregion_dist_zScore , aes(x=Var2, y=Var1, fill = value, label=label)) + geom_tile() + geom_text() +
    scale_fill_gradient2( low = "white", mid = "grey80", high = "black") + theme_classic() +
    xlab("") + ylab("") + ggtitle("t-type layer distribution")
  
  ggsave(paste(base, "/Supplement/FigureS6O_", clone_type, ".pdf", sep=""), width=6, height=6, plot = main_f)  
  
  # do a combined main plot
  layout <- 
  "ABCDEF"
  
  plot_main_first <- main_a + main_b + main_c + main_d + main_e + main_f + plot_layout( design = layout )
  plot_main_first <- plot_main_first & theme(legend.position= "none")
  ggsave(filename = paste(base, "/Figures/Figure4CDE_S6JMO.", clone_type,".pdf", sep=""), plot = plot_main_first, width = 17)
  
  # do a Supplemental plot
  layout <- "AB\nCC\nDD"
  main_b_supplement <- main_b_supplement + NoLegend()
  main_c_supplement <- main_c_supplement + NoLegend()
  
  plot_sup_main_first <- main_b_supplement + main_c_supplement + main_e_supplement + plot_layout( design = layout, guides = "collect" )
  ggsave(filename = paste(base, "/Supplement/FigureS6KLN.", clone_type,".pdf", sep=""), plot = plot_sup_main_first, width = 6, height = 8)



wb <- createWorkbook()
sheetName <- "final meta_data"
addWorksheet(wb, sheetName)

writeData(wb, sheetName, meta_data)
addFilter(wb, sheetName, row = 1, cols = 1:ncol(meta_data))
setColWidths(wb, sheetName, cols = 1:ncol(meta_data), widths="auto")

saveWorkbook(wb, file = paste(base, "/Supplement/SupTable_PatchSeq_meta_data_253cells.xlsx", sep=""), overwrite = T) 

### end of analysis ###  

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
#   [1] BiRewire_3.22.0    Matrix_1.3-4       tsne_0.1-3         slam_0.1-48        igraph_1.2.6       pvclust_2.2-0      ggupset_0.3.0      future_1.22.1      data.table_1.14.0  scales_1.1.1      
# [11] pheatmap_1.0.12    cowplot_1.1.1      GenBinomApps_1.1   patchwork_1.1.1    dplyr_1.0.7        openxlsx_4.2.4     ggbeeswarm_0.6.0   ggplot2_3.3.5      SeuratObject_4.0.2 Seurat_4.0.4.9000 
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_2.0-2      deldir_0.2-10         ellipsis_0.3.2        ggridges_0.5.3        rstudioapi_0.13       spatstat.data_2.1-0   farver_2.1.0         
# [9] leiden_0.3.9          listenv_0.8.0         ggrepel_0.9.1         fansi_0.5.0           codetools_0.2-16      splines_4.0.3         polyclip_1.10-0       jsonlite_1.7.2       
# [17] packrat_0.7.0         ica_1.0-2             cluster_2.1.0         png_0.1-7             uwot_0.1.10           shiny_1.6.0           sctransform_0.3.2     spatstat.sparse_2.0-0
# [25] compiler_4.0.3        httr_1.4.2            assertthat_0.2.1      fastmap_1.1.0         lazyeval_0.2.2        cli_3.0.1             limma_3.46.0          later_1.3.0          
# [33] htmltools_0.5.2       tools_4.0.3           gtable_0.3.0          glue_1.4.2            RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.7            scattermore_0.7      
# [41] vctrs_0.3.8           nlme_3.1-149          lmtest_0.9-38         stringr_1.4.0         globals_0.14.0        mime_0.11             miniUI_0.1.1.1        lifecycle_1.0.0      
# [49] irlba_2.3.3           goftest_1.2-2         MASS_7.3-53           zoo_1.8-9             spatstat.core_2.3-0   promises_1.2.0.1      spatstat.utils_2.2-0  parallel_4.0.3       
# [57] RColorBrewer_1.1-2    reticulate_1.20       pbapply_1.5-0         gridExtra_2.3         rpart_4.1-15          stringi_1.7.4         zip_2.2.0             rlang_0.4.11         
# [65] pkgconfig_2.0.3       matrixStats_0.60.1    lattice_0.20-41       ROCR_1.0-11           purrr_0.3.4           tensor_1.5            labeling_0.4.2        htmlwidgets_1.5.4    
# [73] tidyselect_1.1.1      parallelly_1.28.1     RcppAnnoy_0.0.19      plyr_1.8.6            magrittr_2.0.1        R6_2.5.1              generics_0.1.0        DBI_1.1.1            
# [81] pillar_1.6.2          withr_2.4.2           mgcv_1.8-33           fitdistrplus_1.1-5    survival_3.2-7        abind_1.4-5           tibble_3.1.4          future.apply_1.8.1   
# [89] crayon_1.4.1          KernSmooth_2.23-17    utf8_1.2.2            spatstat.geom_2.2-2   plotly_4.9.4.1        grid_4.0.3            digest_0.6.27         xtable_1.8-4         
# [97] tidyr_1.1.3           httpuv_1.6.3          munsell_0.5.0         beeswarm_0.4.0        viridisLite_0.4.0     vipor_0.4.5          

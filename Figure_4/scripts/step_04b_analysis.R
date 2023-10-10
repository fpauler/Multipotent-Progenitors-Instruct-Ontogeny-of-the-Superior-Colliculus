# produces Figures - focuses on clone specific analyses
#checked 6.10.2023

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

# upsetplot
library (ggupset)
library (pvclust)

library (dendextend)

# prepare randomisation based on Ratz et al., bioRvix 2022, which cites Iorio et al., BMC Bioinformatics 2016
# We randomized the clone-cell type associations, while preserving the number of cell types related
# to each clone, and the number of clones related to each cell type, to create 1000 randomized data sets
# this is done with this package:
library (BiRewire)

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

# Initially we did all analyses both at the clonal and subclone level
# NOTE: this makes a difference to all analyses steps due to the clone_coverage cutoff!
# results are very similar, just there are much less clones than subclones
# also at E10 mainly symmetric division happen, so subclones also originate from RGPs
  
# to simplify for final submission only analyse subclones here

  clone_type <- "subclone"
  cc_cor <- 0.1 
  cc_score <- 0.5
  
  meta_data <- extract_cells_from_meta_data( seurat_meta = PatchSeq_Seurat_SCT@meta.data, total_clone_size = total_cells_clone, 
                                             cc_cor = cc_cor, cc_score = cc_score, clone_coverage = 0.3, clone_type = clone_type)
  
  # identify the column of the current clone type
  clone_type_idx <- which(colnames(meta_data) == clone_type)
  
  ########################################
  ##
  # clone specific analysis
  ##
  ########################################
  
  # simple statistics about the clones
  clone_count <- plyr::ddply(meta_data, clone_type_idx, summarise,
                             N    = length(Clone),
                             inh  = length(which(grepl(pattern = "INH", predicted_ttype ))),
                             exc  = length(which(grepl(pattern = "GLU", predicted_ttype ))),
                             total_celltype = length(unique(predicted_ttype)),
                             genotype = unique(Genotype),
                             cell_types = length( unique ( broad_class )),
                             layers = paste(sort(unique(Subregion)), collapse=",")
  )
  
  clone_count$coverage <- clone_count$N / total_cells_clone[[clone_type]][as.character(clone_count[,1])]
  clone_count$total_clone_size <- total_cells_clone[[clone_type]][as.character(clone_count[,1])]
  
  # some clone details for the text:
  max (clone_count$N) # 12
  min (clone_count$N) # 2
  mean(clone_count$N) # 4.362069
  
  max (clone_count$coverage) # 1
  min (clone_count$coverage) # 0.31
  mean(clone_count$coverage) # 0.6
  
  cells <- nrow(meta_data)
  clones <- length(unique(meta_data[,clone_type_idx]))
  
  # main plot - basic clone stats - not used as such for now
  main_cells_per_clone <- ggplot() + 
    geom_boxplot(  data = clone_count, aes(x="a", y=N) ) +
    geom_beeswarm( data = clone_count, aes(x="a", y=N), cex = 5 ) +
    geom_text(  aes(x=1, y=max(clone_count$N) + 1, label = nrow(meta_data))) +
    geom_text(  aes(x=1, y=max(clone_count$N) + 2, label = length(unique(meta_data[,clone_type_idx])))) + 
    ylim(0, (max(clone_count$N) + 2)) +
    theme_classic() + ggtitle("# informative cells\nper clone") + xlab("") + ylab("")
  
  main_cov_per_clone <- ggplot() + 
    geom_boxplot(  data = clone_count, aes(x="a", y=coverage) ) +
    geom_beeswarm( data = clone_count, aes(x="a", y=coverage), cex = 5 ) + ylim(0,1.2) +
    theme_classic() + ggtitle("# fraction informative cells\nper clone") + xlab("") + ylab("")
  
  # Supplement plots - basic stats
  # duplicate the data to also plot all clones together
  stats2plot <- clone_count[,c("N", "coverage", "genotype")]
  tmp <- stats2plot
  tmp$genotype <- "all"
  stats2plot <- rbind(stats2plot, tmp)
  stats2plot$genotype <- gsub(pattern = "M11GTTG-", replacement = "", x = stats2plot$genotype)
  stats2plot$genotype <- gsub(pattern = "-CreER", replacement = "", x = stats2plot$genotype)
  
  # test for significant difference cells per clone
  t.test( x = stats2plot[which(stats2plot$genotype == "Fz10"), "N"], y = stats2plot[which(stats2plot$genotype == "Sox2"), "N"], paired = F)
  
  # Welch Two Sample t-test
  # 
  # data:  stats2plot[which(stats2plot$genotype == "Fz10"), "N"] and stats2plot[which(stats2plot$genotype == "Sox2"), "N"]
  # t = -1.6328, df = 55.389, p-value = 0.1082
  # alternative hypothesis: true difference in means is not equal to 0
  # 95 percent confidence interval:
  #   -2.4392886  0.2488125
  # sample estimates:
  #   mean of x mean of y 
  # 3.833333  4.928571 
  
  #coverage per clone
  # test for significant difference
  t.test( x = asin(sqrt(stats2plot[which(stats2plot$genotype == "Fz10"), "coverage"])), 
          y = asin(sqrt(stats2plot[which(stats2plot$genotype == "Sox2"), "coverage"])), paired = F)
  
  # Welch Two Sample t-test
  # 
  # data:  asin(sqrt(stats2plot[which(stats2plot$genotype == "Fz10"), "coverage"])) and asin(sqrt(stats2plot[which(stats2plot$genotype == "Sox2"), "coverage"]))
  # t = -1.7521, df = 55.616, p-value = 0.08527
  # alternative hypothesis: true difference in means is not equal to 0
  # 95 percent confidence interval:
  #   -0.22260152  0.01490327
  # sample estimates:
  #   mean of x mean of y 
  # 0.8532729 0.9571220 
  
  Supplement_cells_per_clone <- ggplot() + 
    geom_boxplot(  data = stats2plot, aes(x="a", y=N) ) +
    geom_beeswarm( data = stats2plot, aes(x="a", y=N), cex = 5 ) +
    ylim(0, (max(stats2plot$N) + 2)) + theme_classic() +
    facet_grid(~genotype) + ggtitle("# informative cells\nper clone") + xlab("") + ylab("")
  
  Supplement_cov_per_clone <- ggplot() + 
    geom_boxplot(  data = stats2plot, aes(x="a", y=coverage) ) +
    geom_beeswarm( data = stats2plot, aes(x="a", y=coverage), cex = 5 ) + ylim(0,1.2) +
    theme_classic() + facet_grid(~genotype) + 
    ggtitle("# fraction informative cells\nper clone") + xlab("") + ylab("")
  
  write.csv(x = stats2plot, file = paste(base, "/Supplement/FigureS7AB.csv", sep=""))
  
  # some stats for adding to figure/text
  tmp_stats <- plyr::ddply(stats2plot, "genotype", summarise,
                           median_cov = median(coverage),
                           median_cells_per_clone = median(N)
  )
  
  write.csv(x = tmp_stats, file = paste(base, "/Supplement/stats.", clone_type,".csv", sep=""))
  
  ##########################################
  ##
  # EXC, INH cell types per clone analysis
  ##
  ##########################################
  
  # order the data separately by genotype
  df2plot <- data.frame()
  for (gen in unique(clone_count$genotype)) {
    tmp <- clone_count[which(clone_count$genotype == gen),]
    tmp$inh <- tmp$inh / tmp$N
    tmp$exc <- tmp$exc / tmp$N
    tmp <- tmp[order(tmp$inh),]
    tmp$idx <- 1:nrow(tmp)
    df2plot <- rbind(df2plot, tmp)
  }
  
  clone_order_all <- df2plot[,1][order(df2plot$inh)]
  
  df2plot <- reshape2::melt(df2plot[, which(colnames(df2plot) %in% c(clone_type, "genotype", "inh", "exc", "idx")) ], 
                            id.vars=c(clone_type, "genotype", "idx"))
  df2plot[,1] <- factor(df2plot[,1], levels = clone_order_all)
  df2plot <- df2plot[order(df2plot[,1]),]
  
  # for plotting all clones I need to prepare a new idx variable
  # this idx has to be identical for the exc and inh fraction
  clone2idx <- 1:length(df2plot[which(df2plot$variable == "inh"), 1])
  names(clone2idx) <- as.character(df2plot[which(df2plot$variable == "inh"), 1])
  df2plot$all_idx <- clone2idx[df2plot[,1]]
  
  #write out raw plotting data
  write.csv(x = df2plot, file = paste(base, "/Supplement/Figure4H_S7C.csv", sep=""))
  
  # plot the fraction of exc/inh cells in each clone - combined for main figure
  Clone_composition_main <- ggplot ( df2plot, aes(x=all_idx, y = value, fill = variable)) + 
    geom_bar(stat="identity") + 
    scale_fill_manual(values = c("inh" = "goldenrod", "exc" = "cornflowerblue")) + theme_classic() +
    xlab("") + ylab("") + theme(legend.position = "none")
  
  ggsave(plot = Clone_composition_main, filename = paste(base, "/Figures/Figure4H.", clone_type, ".pdf", sep=""), width=2260, height = 1620, units = "px")
  
  # plot the fraction of exc/inh cells in each clone - for supplement separated by genotype
  Supplement_clone_composition <- ggplot ( df2plot, aes(x=idx, y = value, fill = variable)) + 
    geom_bar(stat="identity") + facet_grid(~ genotype) +
    scale_fill_manual(values = c("inh" = "goldenrod", "exc" = "cornflowerblue")) + theme_classic() +
    xlab("") + ylab("") + theme(legend.position = "bottom")
  
  # focus on clones with exc and inh cells - not used at the moment
  clone_count_ExcInh <- clone_count[which(clone_count$cell_types == 2),]
  cells <- sum(clone_count_ExcInh$N)
  clones <- nrow(clone_count_ExcInh)
  
  # overview analysis of clones with exc/inh t-types
  clone_count_plot_mat <- reshape2::melt(table(clone_count$cell_types))
  clone_count_plot_mat$Var1 <- as.factor(clone_count_plot_mat$Var1)
  clone_count_plot_mat$perc <- clone_count_plot_mat$value / sum(clone_count_plot_mat$value)
  
  # Var1 value      perc
  # 1    1    13 0.2241379
  # 2    2    45 0.7758621
  
  write.csv(x = clone_count_plot_mat, file = paste(base, "/Supplement/Figure4I.csv", sep=""))
  
  exc_inh_plot <- ggplot( clone_count_plot_mat , aes( x = "# cell types", y=perc, fill = Var1) ) +
    geom_bar(stat="identity", color="black") + 
    scale_fill_manual(values = c("1" = "white", "2" = "black")) +
    xlab("") + ylab("") + theme_classic() +
    ggtitle( "# clones with\nexh and inh neurons")
  
  ggsave(plot = exc_inh_plot, filename = paste(base, "/Figures/Figure4I.", clone_type, ".pdf", sep=""), width=760, height = 1000, units = "px")
  
  # Supplemental plot to show the difference in the clone composition between clones with
  # 1 or 2 exc/inh t-types
  df2plot <- reshape2::melt(clone_count[,c("cell_types", "coverage", "total_clone_size", "N")], id.vars = c("cell_types"))
  df2plot$variable <- as.character(df2plot$variable)
  df2plot$variable[which(df2plot$variable == "N")] <- "# informative cells"
  
  # test for significance
  median(clone_count[which(clone_count$cell_types == 1), "N"])
  # 2
  median(clone_count[which(clone_count$cell_types == 2), "N"])
  # 4
  
  median(clone_count[which(clone_count$cell_types == 1), "total_clone_size"])
  # 4
  median(clone_count[which(clone_count$cell_types == 2), "total_clone_size"])
  # 6
  
  t.test(x = clone_count[which(clone_count$cell_types == 2), "total_clone_size"], y = clone_count[which(clone_count$cell_types == 1), "total_clone_size"], paired = F)
  
  # Welch Two Sample t-test
  # 
  # data:  clone_count[which(clone_count$cell_types == 2), "total_clone_size"] and clone_count[which(clone_count$cell_types == 1), "total_clone_size"]
  # t = 2.3155, df = 30.053, p-value = 0.02759
  # alternative hypothesis: true difference in means is not equal to 0
  # 95 percent confidence interval:
  #   0.3392891 5.4077195
  # sample estimates:
  #   mean of x mean of y 
  # 8.488889  5.615385 
  
  t.test(x = clone_count[which(clone_count$cell_types == 2), "N"], y = clone_count[which(clone_count$cell_types == 1), "N"], paired = F)
                   
  # Welch Two Sample t-test
  # 
  # data:  clone_count[which(clone_count$cell_types == 2), "N"] and clone_count[which(clone_count$cell_types == 1), "N"]
  # t = 3.2671, df = 38.083, p-value = 0.002305
  # alternative hypothesis: true difference in means is not equal to 0
  # 95 percent confidence interval:
  #   0.7055451 3.0038566
  # sample estimates:
  #   mean of x mean of y 
  # 4.777778  2.923077 
  
  clone_count$cov_asin <- asin(sqrt(clone_count$coverage))
  t.test(x = clone_count[which(clone_count$cell_types == 2), "cov_asin"], y = clone_count[which(clone_count$cell_types == 1), "cov_asin"], paired = F)
  
  # Welch Two Sample t-test
  # 
  # data:  clone_count[which(clone_count$cell_types == 2), "cov_asin"] and clone_count[which(clone_count$cell_types == 1), "cov_asin"]
  # t = 0.31974, df = 18.672, p-value = 0.7527
  # alternative hypothesis: true difference in means is not equal to 0
  # 95 percent confidence interval:
  #   -0.134157  0.182469
  # sample estimates:
  #   mean of x mean of y 
  # 0.9088212 0.8846652 
  
  # write out raw plotting data
  write.csv(x = df2plot, file = paste(base, "/Supplement/FigureS7DEF.csv", sep=""))
  
  Sup_exc_inh_plot <- ggplot( df2plot, aes(x=as.factor(cell_types), y=value, group = cell_types)) + 
    geom_boxplot() + geom_beeswarm(cex=3) +
    theme_classic() + xlab("# exh/inh types in clone") + ylab("") +
    facet_wrap(~variable, scales = "free_y")
  
  ###############################################################
  ##
  # test for number of different cell types found in each clone with a random background
  ##
  ###############################################################
  
  cells_per_clone <- table(meta_data[, clone_type_idx ]) 
  
  # prepare a set of random clones to test for unique cell types in each clone as a function of clone size
  ct_rand_clone <-  sapply(1:1000, function (x) {
    sapply(sort(unique(cells_per_clone)), function (clone) {
      cell_types <- sample( x = unique(meta_data$predicted_ttype), size = clone, replace = T )
      length(unique(cell_types))
    })
  })
  
  rand_cell_type_clone <- data.frame(cells_per_clone = sort(unique(cells_per_clone)), 
                                     mean_ct_per_clone = apply(ct_rand_clone, 1, mean),
                                     sd_ct_per_clone = apply(ct_rand_clone, 1, sd))
  
  # calculate the real number of different cell types in each clone
  cell_type_clone <- plyr::ddply(meta_data, clone_type_idx, summarise,
                                 string    = paste(unique(predicted_ttype), collapse=","),
                                 N_uniq    = length(unique(predicted_ttype)),
                                 N_tot     = length(predicted_ttype))
  
  # compress the data for better visualisation - not used at the moment
  cell_type_clone <- plyr::ddply(cell_type_clone, c("N_tot", "N_uniq"), summarise,
                                 N_uniq    = unique(N_uniq),
                                 N_tot_sum     = length(N_tot))
  
  ct_clone_plot <- ggplot() + 
    geom_ribbon( data = rand_cell_type_clone, aes(x=cells_per_clone, y=mean_ct_per_clone, 
                                                  ymin=mean_ct_per_clone-sd_ct_per_clone, 
                                                  ymax=mean_ct_per_clone+sd_ct_per_clone), fill="grey80") +
    geom_line(data = rand_cell_type_clone, aes(x=cells_per_clone, y=mean_ct_per_clone), color="grey60") + 
    geom_hline(yintercept = c(2,3,4), color="grey80") +
    geom_point( data = cell_type_clone, aes(y=N_uniq, x=N_tot, size=3)) +
    theme_classic() + ylab("") + ggtitle(paste("# cell types in", length(cells_per_clone), "clones")) +
    theme(legend.position="none") 
  
  ggsave(plot = ct_clone_plot, filename = paste(base, "/Figures/Figure4L.", clone_type, ".pdf", sep=""), width=2260, height = 1620, units = "px")
  
  # write out raw plotting data
  write.csv(x = rand_cell_type_clone, file = paste(base, "/Supplement/Figure4L_rand.csv", sep=""), row.names = F)
  write.csv(x = cell_type_clone, file = paste(base, "/Supplement/Figure4L_real.csv", sep=""), row.names = F)
  
  # calculate z-score and associated p-value
  # re-calculate the real number of different cell types in each clone
  cell_type_clone <- plyr::ddply(meta_data, clone_type_idx, summarise,
                                 string    = paste(unique(predicted_ttype), collapse=","),
                                 N_uniq    = length(unique(predicted_ttype)),
                                 N_tot     = length(predicted_ttype))
  
  #calculate the mean of t-types in clones in cell_per_clone bin 
  av_cell_types <- plyr::ddply(cell_type_clone, c("N_tot"), summarise,
                               real_average    = mean(N_uniq),
                               real_sd    = sd(N_uniq))
  
  # add the mean and stdev from randomisation
  av_cell_types <- merge(av_cell_types, rand_cell_type_clone, by.x="N_tot", by.y="cells_per_clone")
  av_cell_types$z_score <- (av_cell_types$real_average - av_cell_types$mean_ct_per_clone) / av_cell_types$sd_ct_per_clone
  
  # calculate a p-value based on z-score, two-sided
  av_cell_types$p <- sapply(av_cell_types$z_score, function (x) { 
    if (x > 0) {
      pvalue <- 2*pnorm(q = x, lower.tail=F) 
    } else {
      pvalue <- 2*pnorm(q = x, lower.tail=T) 
    }
    pvalue <- formatC(pvalue, format = "e", digits = 2)
    return (pvalue)
  })
  
  av_cell_types$p.adjust <- p.adjust(p = av_cell_types$p, method = "bonferroni")
  
  write.csv(x = av_cell_types, file = paste(base, "/Supplement/Figure4L_pvalues.", clone_type, ".csv", sep=""), row.names = F)
  
  #############################################################################
  ##
  # analyse the total number of t-type pairs found in the data and at random
  ##
  #############################################################################
  
  # create a binary matrix for clone <-> t-type association
  ttype_binary_matrix <- sapply(unique(meta_data[,clone_type_idx]), function (clone_id) {
    sapply(unique(meta_data$predicted_ttype), function (ttype) {
      any( meta_data[,clone_type_idx] == clone_id & meta_data$predicted_ttype == ttype)
    })
  })
  
  # analyse the type and frequency of t-type pairings in real data
  real_tt_comb_count <- ttype_comb_count( t(ttype_binary_matrix) )
  perc_real_tt_comb_count <- length(real_tt_comb_count) / ncol(combn(16,2))
  # 118/120
  
  #define output path for randomisations
  rand_path <- paste(base, '/biwire', sep='')
  
  # do 1000 randomisations
  birewire.sampler.bipartite( t(ttype_binary_matrix), K=1000, path=rand_path, verbose=F, write.sparse = F)
  
  # obtain a data.table with combined randomizations
  comb_rand_df <- read_randomizations ( rand_path )
  
  comb_rand_df[is.na(comb_rand_df)] <- 0
  
  # each row is one t-type combination, each column one randomisation
  # non-zero cases are t-type combinations detected in the data
  perc_tt_comb_rand <- apply( comb_rand_df[,2:ncol(comb_rand_df)], 2, function (x) length(which(x > 0)) / ncol(combn(16,2)) )
  
  # calculate mean and sd for randomisations
  mean_perc_rand <- mean ( perc_tt_comb_rand )
  sd_perc_rand <- sd ( perc_tt_comb_rand )
  
  perc_ct_pairs_df <- data.frame(x = c("real", "random"), 
                                 y = c(perc_real_tt_comb_count, mean_perc_rand), 
                                 sd = c(0, sd_perc_rand))
  
  perc_ct_pairs_plot <- ggplot(perc_ct_pairs_df, aes(x=x, y=y, fill = x)) + 
    geom_bar( stat="identity", color = "black") +
    geom_errorbar( aes(ymin=y-sd, ymax=y+sd), width = 0.5) +
    scale_fill_manual(values = c("random" = "white", "real" = "grey60")) +
    ylab("") + xlab("") +
    theme_classic() + theme(legend.position = "none") + ggtitle("% possible cell type\npairings in clones")
  
  ggsave(plot = perc_ct_pairs_plot, filename = paste(base, "/Figures/Figure4M.", clone_type, ".pdf", sep=""), width=760, height = 1000, units = "px")
  
  ###############################################################################
  ##  
  # layer spanning t-type linkage 
  ##
  ###############################################################################
  
  meta_data$predicted_ttype.layer <- paste(meta_data$predicted_ttype, meta_data$Subregion, sep=".")
  meta_data$broad_class.layer <- paste(meta_data$broad_class, meta_data$Subregion, sep=".")
  
  # potentially loop to produce different type of plots - not used in final analysis
  # only one type of plot is produced
  
    ttype_class <- "broad_class.layer"

    ttype_idx <- which( colnames(meta_data) == ttype_class)
    
    # create a binary matrix for clone <-> t-type association
    ttypeLayer_binary_matrix <- sapply(unique(meta_data[,clone_type_idx]), function (clone_id) {
      sapply( unique( meta_data[,ttype_idx] ), function (ttype) {
         any( meta_data[,clone_type_idx] == clone_id & meta_data[,ttype_idx] == ttype)
      })
    })
    
    # convert the binary matrix to numeric - necessary for pheatmap
    ttypeLayer_binary_matrix[which(ttypeLayer_binary_matrix)] <- 1
    
    # prepare a row dendrogram with defined order
    dst <- dist(x = ttypeLayer_binary_matrix, method = "binary")
    hc_row <- hclust(d = dst, method = "ward.D")
    ol_dend <- rotate(as.dendrogram(hc_row), order = c("EXC.sSC", "INH.sSC", "EXC.dSC", "INH.dSC", "EXC.PAG", "INH.PAG"))
      
    pheatmap( mat = ttypeLayer_binary_matrix, cluster_rows=as.hclust(ol_dend), 
                clustering_distance_cols = "binary", clustering_method = "ward.D", silent = F, legend = F,
                color = colorRampPalette(colors = c("white", "black"))(10),
                filename = paste(base, "/Supplement/FigureS7K.",ttype_class, ".", clone_type, ".pdf", sep=""), width=10, height = 8 )

    # supplement to show significance of t-type clustering
    pvclust_output <- pvclust(t(ttypeLayer_binary_matrix), method.dist = "binary", method.hclust="ward.D", nboot = 10000)
    pdf(file = paste(base, "/Supplement/FigureS7K_stats.", clone_type, ".pdf", sep=""), width = 12, height = 9)
    plot(pvclust_output, hang = -1, cex = 0.5)
    pvrect(pvclust_output, alpha=0.95)
    dev.off()
  
  #######################################################
  ##
  # upset plot showing how clones link ttypes in layers
  ##
  #######################################################
  
  clone_layer_link <- plyr::ddply(meta_data, clone_type_idx, summarise,
                                  layers    = paste(sort(unique(Subregion)), collapse="."),
                                  n_layers  = length(unique(Subregion)),
                                  ttypes = paste(unique(predicted_ttype.layer), collapse=".")
  )
  
  write.csv(x = clone_layer_link, file = paste(base, "/Supplement/FigureS7J.", clone_type, ".csv", sep=""), row.names = F)
    
  table (clone_layer_link$n_layers)
  #  1  2  3 
  # 17 27 14   
  
  # prepare a tibble with layer associations as a list - necessary for the upset plot
  upset_layer_link <- tibble()
  for (i in 1:nrow(clone_layer_link)) {
    tmp <- tibble(layer_span = list(strsplit(clone_layer_link$layers[i], split = "\\.")[[1]]))
    upset_layer_link <- rbind(upset_layer_link, tmp)
  }
  
  upsetplot_layers <- ggplot(upset_layer_link, aes( layer_span )) + geom_bar() + 
    scale_x_upset(sets = c( "sSC", "dSC", "PAG" )) + theme_classic() + 
    scale_y_continuous( breaks = seq(2, 12, 2) ) +
    ylab("") + xlab("") + ggtitle("# clones with\nlayer specific t-types")
  ggsave(plot = upsetplot_layers, filename = paste(base, "/Supplement/FigureS7J.", clone_type, ".pdf", sep=""), width=760, height = 1000, units = "px")
  
  ################
  ##
  #  Supplement
  ##
  ################
  
  layout <- 
    "AB\nCC\nDD"
  
  sup_plot_second <- Supplement_cells_per_clone + Supplement_cov_per_clone + Supplement_clone_composition + 
    Sup_exc_inh_plot + plot_layout( design = layout )
  
  ggsave(filename = paste(base, "/Supplement/FigureS7ABCDEF.", clone_type, ".pdf", sep=""), plot = sup_plot_second, width=10, height=15)


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
#   [1] BiRewire_3.22.0    Matrix_1.3-4       tsne_0.1-3         slam_0.1-48        igraph_1.2.6       dendextend_1.15.1  pvclust_2.2-0      ggupset_0.3.0      scales_1.1.1       pheatmap_1.0.12   
# [11] cowplot_1.1.1      GenBinomApps_1.1   patchwork_1.1.1    dplyr_1.0.7        openxlsx_4.2.4     ggbeeswarm_0.6.0   ggplot2_3.3.5      SeuratObject_4.0.2 Seurat_4.0.4.9000 
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15            colorspace_2.0-2      deldir_0.2-10         ellipsis_0.3.2        ggridges_0.5.3        rstudioapi_0.13       spatstat.data_2.1-0   farver_2.1.0         
# [9] leiden_0.3.9          listenv_0.8.0         ggrepel_0.9.1         fansi_0.5.0           codetools_0.2-16      splines_4.0.3         polyclip_1.10-0       jsonlite_1.7.2       
# [17] packrat_0.7.0         ica_1.0-2             cluster_2.1.0         png_0.1-7             uwot_0.1.10           shiny_1.6.0           sctransform_0.3.2     spatstat.sparse_2.0-0
# [25] compiler_4.0.3        httr_1.4.2            assertthat_0.2.1      fastmap_1.1.0         lazyeval_0.2.2        cli_3.0.1             later_1.3.0           htmltools_0.5.2      
# [33] tools_4.0.3           gtable_0.3.0          glue_1.4.2            RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.7            scattermore_0.7       vctrs_0.3.8          
# [41] nlme_3.1-149          lmtest_0.9-38         stringr_1.4.0         globals_0.14.0        mime_0.11             miniUI_0.1.1.1        lifecycle_1.0.0       irlba_2.3.3          
# [49] goftest_1.2-2         future_1.22.1         MASS_7.3-53           zoo_1.8-9             spatstat.core_2.3-0   promises_1.2.0.1      spatstat.utils_2.2-0  parallel_4.0.3       
# [57] RColorBrewer_1.1-2    reticulate_1.20       pbapply_1.5-0         gridExtra_2.3         rpart_4.1-15          stringi_1.7.4         zip_2.2.0             rlang_0.4.11         
# [65] pkgconfig_2.0.3       matrixStats_0.60.1    lattice_0.20-41       ROCR_1.0-11           purrr_0.3.4           tensor_1.5            labeling_0.4.2        htmlwidgets_1.5.4    
# [73] tidyselect_1.1.1      parallelly_1.28.1     RcppAnnoy_0.0.19      plyr_1.8.6            magrittr_2.0.1        R6_2.5.1              generics_0.1.0        DBI_1.1.1            
# [81] pillar_1.6.2          withr_2.4.2           mgcv_1.8-33           fitdistrplus_1.1-5    survival_3.2-7        abind_1.4-5           tibble_3.1.4          future.apply_1.8.1   
# [89] crayon_1.4.1          KernSmooth_2.23-17    utf8_1.2.2            spatstat.geom_2.2-2   plotly_4.9.4.1        viridis_0.6.1         grid_4.0.3            data.table_1.14.0    
# [97] digest_0.6.27         xtable_1.8-4          tidyr_1.1.3           httpuv_1.6.3          munsell_0.5.0         beeswarm_0.4.0        viridisLite_0.4.0     vipor_0.4.5 

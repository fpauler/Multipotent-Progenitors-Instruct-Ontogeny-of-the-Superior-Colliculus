# 2 cell clone analysis
# step 4: prepare plots
# checked 10.10.2023

library(ggbeeswarm)
library (Seurat)
library (dplyr)
library (openxlsx)

set.seed(2401)

analysis_base <- "./"

PatchSeq_Seurat_SCT <- readRDS( file = paste(analysis_base, "RDS_files/PatchSeq_Seurat_SCT_cellTypeAssignment.RDS", sep="/") )

meta_data <- PatchSeq_Seurat_SCT@meta.data

wb <- createWorkbook()
sheetName <- "final meta_data"
addWorksheet(wb, sheetName)

writeData(wb, sheetName, meta_data)
addFilter(wb, sheetName, row = 1, cols = 1:ncol(meta_data))
setColWidths(wb, sheetName, cols = 1:ncol(meta_data), widths="auto")

saveWorkbook(wb, file = paste(analysis_base, "Supplement/final_metadata.xlsx", sep=""), overwrite = T) 

###########################
##
#  plot some basic stats
##
###########################

basic_stat_1 <- ggplot( meta_data, aes(x="number of genes", y = nFeature_RNA / 1000)) + 
  geom_boxplot(outlier.shape = NA, color="grey40") + geom_beeswarm(size=0.25, cex = 1.5) + 
  ylim(0, max(meta_data$nFeature_RNA / 1000)) +
  theme_classic() + ggtitle(paste("number of genes x 1000\nmedian:\n", median(meta_data$nFeature_RNA))) + 
  ylab("") + xlab("")

basic_stat_2 <- ggplot( meta_data, aes(x="number of reads", y = STAR.unique.aligned / 1000000)) + 
  geom_boxplot(outlier.shape = NA, color="grey40") + geom_beeswarm(size=0.25, cex = 1.5) + 
  ylim(0, max(meta_data$STAR.unique.aligned / 1000000)) + ylim(0,4) +
  theme_classic() + ggtitle(paste("unique M reads \nmedian:\n", median(meta_data$STAR.unique.aligned))) + 
  ylab("") + xlab("")

basic_stat_3 <- ggplot( meta_data, aes(x="alignment rate", y = STAR.perc)) + 
  geom_boxplot(outlier.shape = NA, color="grey40") + geom_beeswarm(size=0.25, cex = 1.5) + ylim(0,1) + 
  theme_classic() +
  ggtitle(paste("alignment rate (unique reads)\nmedian:\n", round(median(meta_data$STAR.perc), 2))) + 
  ylab("") + xlab("")

hist_data <- hist( meta_data$boot_score, plot = F, breaks = seq(0,1,0.05))
hist_datadf<- data.frame(x = hist_data$breaks[2:length(hist_data$breaks)], y = hist_data$counts)

boot_score_plot <- ggplot(hist_datadf, aes(x=x, y=y)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_y_continuous(breaks = seq(0,18,6)) + scale_x_continuous(breaks = seq(0,1,0.2)) + 
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

# Note: This plot did not make it into the final paper due to limited supplemental space
plot_basic_stats <- (basic_stat_1 + basic_stat_2 + basic_stat_3) / (NMS_plot + boot_score_plot + pearson_cor_plot)
ggsave(paste(analysis_base, "/Supplement/FigureSX_cell_type_assignment_stats.pdf", sep=""), plot = plot_basic_stats, width = 6, height=6)

# analyse exc/inh neurons first
meta_data$broad_class <- ""
meta_data$broad_class[ which( grepl(pattern = "INH", x = meta_data$predicted_ttype )) ] <- "INH"
meta_data$broad_class[ which( grepl(pattern = "GLU", x = meta_data$predicted_ttype )) ] <- "EXC"

# plot the broad class
df2plot <- as.data.frame( table(meta_data[, "broad_class"]) )
ggplot( df2plot, aes(x=Var1, y=Freq, fill=Var1)) + geom_bar(stat = "identity") + 
  scale_fill_manual( values = c("EXC" = "goldenrod", "INH" = "cornflowerblue")) + 
  theme_classic() + ylab( "number of cells") + xlab("")
ggsave(filename = paste(analysis_base, "/Figures/Figure_5C.pdf", sep="") )

broad_type_count <- plyr::ddply(meta_data, c("Clone.ID"), summarise,
                                   broad_type = length(unique(broad_class) ),
                                exc = length( which( broad_class == "EXC" ) ),
                                inh = length( which( broad_class == "INH" ) ),
                                gen = unique(Genotype)
)

broad_type_count$frac_exc <- broad_type_count$exc / (broad_type_count$exc + broad_type_count$inh)
broad_type_count$frac_inh <- 1-broad_type_count$frac_exc

tmp <- broad_type_count[, c("Clone.ID", "frac_exc", "gen")]
colnames(tmp) <- c("Clone", "frac", "genotype")
tmp$type <- "exc"

df2plot <- broad_type_count[, c("Clone.ID", "frac_inh", "gen")]
colnames(df2plot) <- c("Clone", "frac", "genotype")
df2plot$type <- "inh"

df2plot <- rbind(df2plot, tmp)

broad_type_count$class <- "both"
broad_type_count$class[ which(broad_type_count$exc == 0)] <- "inh"
broad_type_count$class[ which(broad_type_count$inh == 0)] <- "exc"

broad_type_count$class <- factor(broad_type_count$class, levels=c("exc", "inh", "both"))
broad_type_count <- broad_type_count[order(broad_type_count$class),]
broad_type_count$Clone.ID <- factor(broad_type_count$Clone.ID, levels=broad_type_count$Clone.ID)

df2plot$Clone <- factor(df2plot$Clone, levels=broad_type_count$Clone.ID)

broad_type_plot <- ggplot(df2plot, aes(x=Clone, y=frac, fill=type)) + geom_bar(stat="identity")+
  scale_fill_manual( values = c("exc" = "goldenrod", "inh" = "cornflowerblue")) + theme_classic()

ttype_count <- plyr::ddply(meta_data, c("Clone.ID"), summarise,
                                n_ttype = length(unique(predicted_ttype))
)

ttype_count$Clone.ID <- factor( ttype_count$Clone.ID, levels = broad_type_count$Clone.ID) 
nttype_plot <- ggplot(ttype_count, aes(x=Clone.ID, y=n_ttype)) + geom_bar(stat="identity") + theme_classic() +
  ylab("# distinct cell types in clone")

plot_dist_cellTypes <- broad_type_plot / nttype_plot

number_clones <- length(unique(meta_data$Clone.ID))

# do randomization of a data set with the same number of clones as in the real data
# consisting of 2 cells each with a given probability that one cell type is chosen
rand_out <-   sapply(1:10000, function (y) {
   rand_out <- sapply( 1:number_clones, function (x) {
      tmp <- sample(x = c("E", "I"), size = 2,  replace = T)
      if (length(unique(tmp)) == 1) {
        return(unique(tmp))
      } else {
        return (2)
      }
    })
   
   tmp <- table(rand_out) / sum(table(rand_out))
   
   # deal with exceptions where one category is missing
   if (!length(tmp) == 3) {
     missing <- c("2", "E", "I")[ which(! c("2", "E", "I") %in% names(tmp)) ]
     old_names <- names(tmp)
     
     tmp <- c(tmp, rep(0, length(missing)))
     
     names(tmp) <- c(old_names, missing)
     tmp <- tmp[c("2", "E", "I")]
   }
   
   return (tmp)
   
  })

df2plot <- reshape2::melt(table( broad_type_count$class ) / number_clones)
df2plot$class <- "real"

tmp <- reshape2::melt( apply(rand_out, 1, mean) )
tmp$Var1 <- c("both", "exc", "inh")
tmp$class <- "rand"

df2plot <- rbind(df2plot, tmp)

sd_df <- apply(rand_out, 1, sd)
sd_df <-  reshape2::melt( sd_df )
sd_df$Var1 <- c("both", "exc", "inh")
sd_df$class <- "rand"

# add zero sd for the real data - only necessary for plotting
sd_df <- rbind( sd_df, data.frame( value = c(0,0,0), Var1 = c("both", "exc", "inh"), class="real"))
colnames(sd_df) <- c("sd", "Var1_sd", "class_sd")

df2plot$Var1 <- as.character(df2plot$Var1)

sd_df <- sd_df[order(sd_df$Var1, sd_df$class),]
df2plot <- df2plot[order(df2plot$Var1, df2plot$class),]
df2plot <- cbind(df2plot, sd_df)

exc_inh_plot <- ggplot( df2plot, aes(x=class, y=value, fill=class)) + geom_bar(stat="identity", position="dodge", color="black") + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd)) + facet_grid(~Var1) +
  theme_classic() + scale_fill_manual( values = c(rand="white", real = "grey30")) + theme(legend.position = "none")

# prepare a table with observations
tmp <- table( broad_type_count$class )
missing <- c("exc", "inh", "both")[ which(! c("exc", "inh", "both") %in% names(tmp)) ]

# this part catches exceptions where not all categories are identified
# also makes sure the order of observations is correct
old_names <- names(tmp)

tmp <- c(tmp, rep(0, length(missing)))

names(tmp) <- c(old_names, missing)
tmp <- tmp[c("exc", "inh", "both")]

observed <- tmp
expected <- c(0.25, 0.25, 0.5)

# do the test
chisq.test(x=observed, p=expected, simulate.p.value = F )
# 
# Chi-squared test for given probabilities
# 
# data:  observed
# X-squared = 1, df = 2, p-value = 0.6065

# random assignment of 16 cell types
rand_out <-   sapply(1:10000, function (y) {
  rand_out <- sapply( 1:number_clones, function (x) {
    sample(x = c(1:16), size = 2, replace = T)
  })
})

count_rand <- mean( apply(rand_out, 2, function (x) length(unique(x))))
sd_rand <- sd( apply(rand_out, 2, function (x) length(unique(x))))

count_real <- length(unique(meta_data$predicted_ttype))

# prepare a dataframe for plotting
df2plot <- data.frame( count = c(count_rand, count_real), sd = c(sd_rand, 0), group = c("rand", "real"))

#   count       sd group
# 14.4457 1.030222  rand
# 15.0000 0.000000  real

cellType_plot <- ggplot(df2plot, aes(x=group, y=count, fill=group)) + geom_bar(stat="identity", position="dodge", color="black") + 
  geom_errorbar( aes(ymin=count-sd, ymax=count+sd)) + ylab("number of celltypes in all cells") + theme_classic() +
  scale_fill_manual( values = c(rand="white", real = "grey30")) + theme(legend.position = "none")

plot_dist_cellTypes / (exc_inh_plot + cellType_plot)

# Note that here one more plot is given, which did not make it into the paper as such
# In the paper the fact that all clones, except for one, contain 2 distict cell types is given as a %
# 15/16 = 0.9375

ggsave(paste(analysis_base, "/Figures/Figure_5D-F.pdf", sep=""), width = 6, height=8)

# calculate a z-score to test for significance of difference of observed to random
zscore <- (df2plot$count[2] - df2plot$count[1]) / df2plot$sd[1]

if (zscore > 0) {
  pvalue <- 2*pnorm(q = zscore, lower.tail=F) 
} else {
  pvalue <- 2*pnorm(q = zscore, lower.tail=T) 
}
pvalue <- formatC(pvalue, format = "e", digits = 2)

# not significant
# [1] "5.91e-01"

# which cell types are missing? - also not in the paper
cgi <- c("MEGLU1", "MEGLU2", "MEGLU3", "MEGLU4", "MEGLU5", "MEGLU6",
         "MEINH2", "MEINH3", "MEINH5", "MEINH6", "MEINH7", "MEINH8", 
         "MEINH9", "MEINH10", "MEINH11", "MEINH12")

cgi[ which(!cgi %in% unique( meta_data$predicted_ttype )) ]

# MEINH9

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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggbeeswarm_0.6.0            scales_1.1.1                umap_0.2.7.0                batchelor_1.6.3             SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0 Biobase_2.50.0             
# [8] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7         IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.1         MatrixGenerics_1.2.1        matrixStats_0.60.1         
# [15] ggrepel_0.9.1               ggplot2_3.3.5               dplyr_1.0.7                 openxlsx_4.2.4              SeuratObject_4.0.2          Seurat_4.0.4.9000          
# 
# loaded via a namespace (and not attached):
#   [1] plyr_1.8.6                igraph_1.2.6              lazyeval_0.2.2            splines_4.0.3             BiocParallel_1.24.1       listenv_0.8.0             scattermore_0.7          
# [8] digest_0.6.27             htmltools_0.5.2           fansi_0.5.0               magrittr_2.0.1            tensor_1.5                cluster_2.1.0             ROCR_1.0-11              
# [15] globals_0.14.0            askpass_1.1               spatstat.sparse_2.0-0     colorspace_2.0-2          crayon_1.4.1              RCurl_1.98-1.4            jsonlite_1.7.2           
# [22] spatstat.data_2.1-0       survival_3.2-7            zoo_1.8-9                 glue_1.4.2                polyclip_1.10-0           gtable_0.3.0              zlibbioc_1.36.0          
# [29] XVector_0.30.0            leiden_0.3.9              DelayedArray_0.16.3       BiocSingular_1.6.0        future.apply_1.8.1        abind_1.4-5               DBI_1.1.1                
# [36] miniUI_0.1.1.1            Rcpp_1.0.7                viridisLite_0.4.0         xtable_1.8-4              reticulate_1.20           spatstat.core_2.3-0       rsvd_1.0.5               
# [43] ResidualMatrix_1.0.0      htmlwidgets_1.5.4         httr_1.4.2                RColorBrewer_1.1-2        ellipsis_0.3.2            ica_1.0-2                 pkgconfig_2.0.3          
# [50] farver_2.1.0              scuttle_1.0.4             uwot_0.1.10               deldir_0.2-10             utf8_1.2.2                tidyselect_1.1.1          labeling_0.4.2           
# [57] rlang_0.4.11              reshape2_1.4.4            later_1.3.0               munsell_0.5.0             tools_4.0.3               cli_3.0.1                 generics_0.1.0           
# [64] ggridges_0.5.3            stringr_1.4.0             fastmap_1.1.0             goftest_1.2-2             fitdistrplus_1.1-5        zip_2.2.0                 purrr_0.3.4              
# [71] RANN_2.6.1                packrat_0.7.0             pbapply_1.5-0             future_1.22.1             nlme_3.1-149              sparseMatrixStats_1.2.1   mime_0.11                
# [78] compiler_4.0.3            rstudioapi_0.13           beeswarm_0.4.0            plotly_4.9.4.1            png_0.1-7                 spatstat.utils_2.2-0      tibble_3.1.4             
# [85] stringi_1.7.4             RSpectra_0.16-0           lattice_0.20-41           Matrix_1.3-4              vctrs_0.3.8               pillar_1.6.2              lifecycle_1.0.0          
# [92] spatstat.geom_2.2-2       lmtest_0.9-38             RcppAnnoy_0.0.19          BiocNeighbors_1.8.2       data.table_1.14.0         cowplot_1.1.1             bitops_1.0-7             
# [99] irlba_2.3.3               httpuv_1.6.3              patchwork_1.1.1           R6_2.5.1                  promises_1.2.0.1          KernSmooth_2.23-17        gridExtra_2.3            
# [106] vipor_0.4.5               parallelly_1.28.1         codetools_0.2-16          MASS_7.3-53               assertthat_0.2.1          openssl_1.4.5             withr_2.4.2              
# [113] sctransform_0.3.2         GenomeInfoDbData_1.2.4    mgcv_1.8-33               grid_4.0.3                rpart_4.1-15              beachmat_2.6.4            tidyr_1.1.3              
# [120] DelayedMatrixStats_1.12.3 Rtsne_0.15                shiny_1.6.0 

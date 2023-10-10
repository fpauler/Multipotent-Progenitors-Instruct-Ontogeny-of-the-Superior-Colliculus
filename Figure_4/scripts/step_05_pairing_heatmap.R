# thi script was added during revision

library (openxlsx)
library (pheatmap)
library (ggplot2)
library (parallel)
library (pbapply)

# define base folder
base <- "./"

set.seed( 2401 )

# convert Zeisel nomenclature to the one used in the paper
Zeisel_to_Cheung <- sapply( c(1:10), function (x) paste("SCINH", x, sep=""))
Zeisel_to_Cheung <- c( Zeisel_to_Cheung, sapply( c(1:6), function (x) paste("SCEXC", x, sep="")) )
names( Zeisel_to_Cheung ) <- c("MEINH2", "MEINH3", "MEINH5", "MEINH6", "MEINH7", "MEINH8", "MEINH9", "MEINH10", "MEINH11", "MEINH12",
                               "MEGLU1", "MEGLU2", "MEGLU3", "MEGLU4", "MEGLU5", "MEGLU6")

# from here, 2-tailed: https://www.statology.org/p-value-of-z-score-r/
p_from_z <- function( zscore ) {
  if (is.na(zscore)) {
    return(NA)
  } else {
    if (zscore > 0) {
      pvalue <- 2*pnorm(q = zscore, lower.tail=F) 
    } else {
      pvalue <- 2*pnorm(q = zscore, lower.tail=T) 
    }
    return (pvalue)
  }
}

# function to prepare a matrix of pairwise frequencies
get_real_freq_mat <- function ( ttype_order = ttype_order, meta_data = meta_data, cloneID = cloneID, ttype = ttype, comb_message = F) {
  
  cloneIDX <- which( colnames(meta_data) == cloneID )
  if (length(cloneIDX) == 0) {
    stop ("cloneID not found in meta data column names")
  }
  
  ttypeIDX <- which( colnames(meta_data) == ttype )
  if (length(ttypeIDX) == 0) {
    stop ("ttype not found in meta data column names")
  }
  
  # calculate frequencies of pairings
  real_comb_clone <- sapply( unique(meta_data[,cloneIDX]), function (clone) {
    all_cells <- meta_data[,ttypeIDX][which(meta_data[,cloneIDX] == clone)]
    all_combinations <- apply(combn( all_cells, 2),2, function(comb) paste(sort(comb), collapse=","))
    unique(all_combinations)
  })
  
  # report the number of unique combinations found - if requested
  n_unique_combinations <- length(unique(unlist(real_comb_clone)))
  if (comb_message) {
    message( "found ", n_unique_combinations, " pairwise ttype combinations")
  }
  
  # prepare a frequency table
  comb_freq <- table( unlist(real_comb_clone) )
  
  # create a matrix for the pairs
  real_freq_mat <- sapply( ttype_order, function (ttype1) {
    sapply( ttype_order, function (ttype2) {
      comb_name <- paste( sort( c(ttype1, ttype2) ), collapse = ",")
      as.numeric(comb_freq[comb_name])
    })
  })
  
  real_freq_mat[is.na(real_freq_mat)] <- 0
  
  return( real_freq_mat )
}

# function to prepare a matrix of pairwise frequencies with randomised data
# currently only the column that is used to determine frequencies is randomised

get_random_freq_mat <- function ( ttype_order = ttype_order, meta_data = meta_data, cloneID = cloneID, ttype = ttype, iter = iter) {
  
  cloneIDX <- which( colnames(meta_data) == cloneID )
  if (length(cloneIDX) == 0) {
    stop ("cloneID not found in meta data column names")
  }
  
  message("cloneIDX ", cloneIDX)
  
  ttypeIDX <- which( colnames(meta_data) == ttype )
  if (length(ttypeIDX) == 0) {
    stop ("ttype not found in meta data column names")
  }
  
  cl <- parallel::makeCluster( 5 )
  parallel::clusterExport( cl=cl, varlist=c("meta_data", "ttype_order", "cloneIDX", "ttypeIDX", "iter"), envir=environment())
  
  # randomization
  rand_list <- pblapply(cl = cl, X=1:iter, FUN=function (x) {
    
    tmp <- meta_data
    tmp[,ttypeIDX] <- sample(tmp[,ttypeIDX])
    
    rand_comb_clone <- sapply( unique(tmp[,cloneIDX]), function (clone) {
      all_cells <- tmp[,ttypeIDX][which(tmp[,cloneIDX] == clone)]
      all_combinations <- apply(combn( all_cells, 2), 2, function(comb) paste(sort(comb), collapse=","))
      unique(all_combinations)
    })
    
    rand_comb_freq <- table( unlist(rand_comb_clone) )
    
    freq_mat <- sapply( ttype_order, function (ttype1) {
      sapply( ttype_order, function (ttype2) {
        comb_name <- paste( sort( c(ttype1, ttype2) ), collapse = ",")
        as.numeric(rand_comb_freq[comb_name])
      })
    })
    
    freq_mat[ is.na(freq_mat) ] <- 0
    
    return( freq_mat )
    
  })
  
  stopCluster(cl)
  
  
  return( rand_list )
}

# read MADM Clone-Seq data
meta_data_first <- read.xlsx(xlsxFile = paste(base_e10_PatchSeq, "Supplement/SupTable_PatchSeq_meta_data_253cells.xlsx", sep="") )
data2analyse <- meta_data_first[, c("predicted_ttype", "subclone")]
colnames(data2analyse) <- c("predicted_ttype", "Clone.ID")

ttype_order <- c("MEGLU1", "MEGLU2", "MEGLU3", "MEGLU4", "MEGLU5", "MEGLU6",
                 "MEINH2", "MEINH3", "MEINH5" ,"MEINH6", "MEINH7", "MEINH8", "MEINH9", "MEINH10", "MEINH11", "MEINH12"  )

# calculate real and randomized matrices for frequency of pairs
real_freq_mat <- get_real_freq_mat (ttype_order = ttype_order, meta_data = data2analyse, cloneID = "Clone.ID", ttype = "predicted_ttype", comb_message = T )
rand_freq_mat <- get_random_freq_mat( ttype_order = ttype_order, meta_data = data2analyse, cloneID = "Clone.ID", ttype = "predicted_ttype", iter = 10000)

# calculate mean and sd of the randomized matrices
mean_mat <- apply(simplify2array(rand_freq_mat), 1:2, function (x) mean( x, na.rm = T))
sd_mat <- apply(simplify2array(rand_freq_mat), 1:2, function (x) sd( x, na.rm = T))

# calculate the z-score
real_freq_mat[is.na(real_freq_mat)] <- 0
real_zscore_mat <- (real_freq_mat - mean_mat ) / sd_mat 

# calculate the p-value from the z-score
mat2plot <- reshape2::melt( real_zscore_mat )
mat2plot$pvalue <- sapply( mat2plot$value, function (x) p_from_z(x) )
mat2plot$pvalue_cor <- p.adjust( mat2plot$pvalue, method = "fdr")

# add a star to the significantly enriched pairs
mat2plot$star <- ""
mat2plot$star[which( mat2plot$pvalue_cor < 0.05 )] <- "*"
mat2plot$value2plot <- mat2plot$value
mat2plot$value2plot[which(!mat2plot$star == "*")] <- NA

# convert the ttype IDs and prepare for plotting
mat2plot$conv_ttype_1 <- Zeisel_to_Cheung[ as.character(mat2plot$Var1) ]
mat2plot$conv_ttype_2 <- Zeisel_to_Cheung[ as.character(mat2plot$Var2) ]

mat2plot$conv_ttype_1 <- factor( mat2plot$conv_ttype_1, levels = Zeisel_to_Cheung)
mat2plot$conv_ttype_2 <- factor( mat2plot$conv_ttype_2, levels = Zeisel_to_Cheung)


ggplot(mat2plot, aes(x=conv_ttype_1, y=conv_ttype_2, fill=value, label=star)) + 
  geom_tile() + geom_text() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_gradient2(low = "blue", mid="white", high = "yellow") + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("")

ggsave( filename = paste(base, "Supplement/FigureS7G.pdf", sep=""))

write.csv(x = mat2plot, file = paste(base, "Supplement/FigureS7G.csv", sep=""))

min(mat2plot$pvalue_cor)
# [1] 0.2824221

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
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pbapply_1.5-0   ggplot2_3.3.5   pheatmap_1.0.12 openxlsx_4.2.4 
# 
# loaded via a namespace (and not attached):
#   [1] zip_2.2.0          Rcpp_1.0.7         pillar_1.6.2       compiler_4.0.3     RColorBrewer_1.1-2 plyr_1.8.6         tools_4.0.3       
# [8] digest_0.6.27      packrat_0.7.0      lifecycle_1.0.0    tibble_3.1.4       gtable_0.3.0       pkgconfig_2.0.3    rlang_0.4.11      
# [15] rstudioapi_0.13    cli_3.0.1          DBI_1.1.1          withr_2.4.2        stringr_1.4.0      dplyr_1.0.7        generics_0.1.0    
# [22] vctrs_0.3.8        grid_4.0.3         tidyselect_1.1.1   glue_1.4.2         R6_2.5.1           fansi_0.5.0        reshape2_1.4.4    
# [29] purrr_0.3.4        farver_2.1.0       magrittr_2.0.1     scales_1.1.1       ellipsis_0.3.2     assertthat_0.2.1   colorspace_2.0-2  
# [36] labeling_0.4.2     utf8_1.2.2         stringi_1.7.4      munsell_0.5.0      crayon_1.4.1 

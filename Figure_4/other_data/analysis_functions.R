# source with a couple of useful functions

############
##
# Functions
##
############

# from here: https://stackoverflow.com/questions/39351013/standard-deviation-over-a-list-of-matrices-in-r
sd_from_mat_list <- function(lst) {
  n <- length(lst)
  rc <- dim(lst[[1]])
  ar1 <- array(unlist(lst), c(rc, n))
  apply(ar1, c(1, 2), sd)
}

# function to extract meta data from Seurat object based on mapping cut offs
# seurat_meta: meta data table from Seurat object
# total_clone_size: named list with clone sizes
# clone_type: has to be one column of seurat_meta and one of the names of total_clone_size

extract_cells_from_meta_data <- function ( seurat_meta = seurat_meta, total_clone_size = total_clone_size, 
                                           cc_cor = 0.1 , cc_score = 0.5, clone_coverage = 0.3, clone_type = "Clone") {
  
  if (!any(colnames(seurat_meta) == clone_type)) stop ("clone_type must match one column in seurat_meta ")
  if (!any(names(total_clone_size) == clone_type)) stop ("clone_type must match one name in total_clone_size ")
  
  message("extracting clones with: clone type: ", clone_type, " cor: ", cc_cor, " score: ",cc_score, " coverage: ", clone_coverage, " cutoff")
  
  # identify the column with the clone_ids information
  idx <- which(colnames(seurat_meta) == clone_type)
  
  # filter cells and meta data based on current cutoff
  filtered_meta <- seurat_meta[which(seurat_meta$cor_ttype > cc_cor & seurat_meta$boot_score_ttype > cc_score),]
  
  # filter based on coverage cutoff
  clone_cov <- table(filtered_meta[,idx])
  clone_cov <- clone_cov / total_clone_size[[clone_type]][names(clone_cov)]
  clones2keep <- names(which(clone_cov > clone_coverage))
  
  filtered_meta <- filtered_meta[which(filtered_meta[,idx] %in% clones2keep),]
  
  # filter out clones with only a single cell
  clones2keep <- names(which(table(filtered_meta[,idx]) > 1))
  filtered_meta <- filtered_meta[which(filtered_meta[,idx] %in% clones2keep),]
  
  return(filtered_meta)
  
}

# function to obtain t-type pairings from binary matrix
ttype_comb_count <- function (ct_matrix = ct_matrix) {
  ttype_combs <- apply( ct_matrix, 1, function (x) {
    ttypes <- names(which(x))
    if (length(ttypes) > 1) {
      ttype_combn <- combn(ttypes, 2)
      sapply(1:ncol(ttype_combn), function (col){
        paste(sort(ttype_combn[,col]), collapse=".")
      })
    } else {
      return(NA)
    }
  })
  
  table(unlist(ttype_combs))
}

# read and combine randomized data ccreated with BiRewire
read_randomizations <- function ( rand_path = rand_path) {
  
  # read all files with randomisations and put into a list of data tables for efficient handling
  files = list.files(path = rand_path, pattern = "network_*", recursive = T)
  
  rand_list <- list()
  i <- 1
  for (file in files) {
    rand_table <- read.table(file = paste(rand_path, "/", file, sep=""), header = T)
    tmp <- ttype_comb_count( rand_table )
    tmp <- data.table::as.data.table (tmp)
    data.table::setnames(tmp, c("ct_comb", i))
    rand_list[[i]] <- tmp
    
    # sanity check that all constraints are met
    if (  any(!rowSums(rand_table) == rowSums( t(ttype_binary_matrix) )) |   any(!colSums(rand_table) == colSums( t(ttype_binary_matrix) )) ) {
      warning(i, " had a problem")
    }
    
    i <- i + 1
  }
  
  # merge all data.tables and modify for further analysis
  message ("merging data table")
  comb_rand_df <- Reduce(function(dtf1, dtf2) data.table::merge.data.table(dtf1, dtf2, all = TRUE), rand_list)
  
  return (comb_rand_df)
}

# prepare a vector with total clone size - important to plot coverage!
# do this for Clones and Subclones
# NOTE: this is tailored for the specific sample_list unsed here!
get_total_clone_size <- function (base_dir=base_dir, sample_list_file=sample_list_file) {
  message("using: ", paste(base_dir, sample_list_file, sep="/"), " as sample list")
  tot_cells_clone <- list()
  for (column in c("Clone", "subclone")) {
    # get details of total number of neurons in each clone
    missed_clones <- read.xlsx(xlsxFile = paste(base_dir, sample_list_file, sep="/"), sheet = 2)
    missed_clones$Color <- toupper(missed_clones$Color)
    missed_clones <- missed_clones[,1:2]
    missed_clones$subclone <- paste(missed_clones$Clone, missed_clones$Color, sep=".")
    
    idx <- which(colnames(missed_clones) == column)
    
    missed_clones <- table(missed_clones[,idx])
    
    patched_clones <- read.xlsx(xlsxFile =paste(base_dir, sample_list_file, sep="/"), sheet = 1)
    patched_clones <- patched_clones[which(grepl(pattern = "PS", x = patched_clones$Clone)),]
    patched_clones$Color <- toupper(patched_clones$Color)
    
    patched_clones$subclone <- paste(patched_clones$Clone, patched_clones$Color, sep=".")
    idx <- which(colnames(patched_clones) == column)
    patched_clones <- table(patched_clones[,idx])
    
    # add zeros to patched clones
    clones2add <- rep(0, length(names(patched_clones)[which(!names(patched_clones) %in% names(missed_clones))]))
    names(clones2add) <- names(patched_clones)[which(!names(patched_clones) %in% names(missed_clones))]
    missed_clones <- c(missed_clones, clones2add)
    
    # these clones miss one subclone completely - none of the cells of this subclone were patched!
    names(missed_clones)[which(!names(missed_clones) %in% names(patched_clones))]
    # output for subclone, Clone is empty
    # [1] "PS01.GREEN" "PS05.GREEN" "PS06.RED"   "PS08.RED"   "PS21.GREEN" "PS27.RED"   "PS44.GREEN"
    
    # cleanup this discrepancy
    missed_clones <- missed_clones[names(patched_clones)]
    tot_cells_clone[[column]] <- missed_clones + patched_clones[names(missed_clones)]
  }
  
  return(tot_cells_clone)
} 

# function for counting layer specific total number of cells
get_total_clone_size_layers <- function (base_dir=base_dir, sample_list_file=sample_list_file) {
  message("using: ", paste(base_dir, sample_list_file, sep="/"), " as sample list")
  tot_cells_clone <- list()
  for (column in c("Clone", "subclone")) {
    
    # get details of total number of neurons in each clone
    missed_clones <- read.xlsx(xlsxFile = paste(base_dir, sample_list_file, sep="/"), sheet = 2)
    missed_clones$Color <- toupper(missed_clones$Color)
   
    missed_clones$subclone <- paste(missed_clones$Clone, missed_clones$Color, sep=".")
    
    idx <- which(colnames(missed_clones) == column)
    
    missed_clones$clone_group <- paste(missed_clones[,idx], missed_clones$Subregion, sep=".")
    
    missed_clones <- table(missed_clones$clone_group)
    
    patched_clones <- read.xlsx(xlsxFile =paste(base_dir, sample_list_file, sep="/"), sheet = 1)
    patched_clones <- patched_clones[which(grepl(pattern = "PS", x = patched_clones$Clone)),]
    patched_clones$Color <- toupper(patched_clones$Color)
    
    patched_clones$subclone <- paste(patched_clones$Clone, patched_clones$Color, sep=".")
    idx <- which(colnames(patched_clones) == column)
    
    patched_clones$clone_group <- paste(patched_clones[,idx], patched_clones$Subregion, sep=".")
    
    patched_clones <- table(patched_clones$clone_group)
    
    # add zeros to patched clones
    clones2add <- rep(0, length(names(patched_clones)[which(!names(patched_clones) %in% names(missed_clones))]))
    names(clones2add) <- names(patched_clones)[which(!names(patched_clones) %in% names(missed_clones))]
    missed_clones <- c(missed_clones, clones2add)
    
    # these clones miss one subclone completely - none of the cells of this subclone were patched!
    names(missed_clones)[which(!names(missed_clones) %in% names(patched_clones))]
    # output for subclone, Clone is empty
    # [1] "PS01.GREEN" "PS05.GREEN" "PS06.RED"   "PS08.RED"   "PS21.GREEN" "PS27.RED"   "PS44.GREEN"
    
    # cleanup this discrepancy
    missed_clones <- missed_clones[names(patched_clones)]
    tot_cells_clone[[column]] <- missed_clones + patched_clones[names(missed_clones)]
  }
  
  return(tot_cells_clone)
} 


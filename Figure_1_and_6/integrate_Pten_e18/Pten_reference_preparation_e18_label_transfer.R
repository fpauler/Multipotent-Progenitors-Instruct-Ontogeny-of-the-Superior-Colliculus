# label transfer from LaManno e18 to Pten data

library ( Seurat )
library( loomR )

# define base folder
base <- "./"

################################  
###
# prepare a LaManno reference
###
################################

lfile <- connect(filename = paste(base, "other_data/dev_all.loom", sep=""), mode = "r+", skip.validate = T)
  
# get samples indices from Midbrain or Hindbrain
mbdIDX <- which( lfile$col.attrs$Tissue[] == "Midbrain" | lfile$col.attrs$Tissue[] == "Hindbrain")
mbdIDX <- mbdIDX[which(!(duplicated(lfile$col.attrs$CellID[mbdIDX]) | duplicated(lfile$col.attrs$CellID[mbdIDX], fromLast = T)))]
  
# focus on e18
mbdIDX <- intersect( mbdIDX, which( lfile$col.attrs$Age[] %in% c("e18.0")))
  
# this allows to recycle older code chunks       
cti_MBd <- mbdIDX
  
#extract the reads from these cells - NOTE row/col are changed according to convention!
message ("extracting count data")
data.subset <- lfile[["matrix"]][cti_MBd, ]
  
#give cellID as row
rownames(data.subset) <- lfile$col.attrs$CellID[cti_MBd]
#give Gene name as col
colnames(data.subset) <- lfile$row.attrs$Gene[]
  
#create the Seurat object
MBd_Seurat <- CreateSeuratObject(counts = as.matrix(t(data.subset)), min.cells = 5, min.features = 500)
  
#check reads in mitocondrial genome
#really just a check - no QC on this data as already pre-analysed
MBd_Seurat[["percent.mt"]] <- PercentageFeatureSet(MBd_Seurat, pattern = "^mt-")
  
#add meta data
meta_data <- list()
meta_data[["cell_type"]] <- lfile$col.attrs$Class[cti_MBd]
names(meta_data[["cell_type"]]) <- lfile$col.attrs$CellID[cti_MBd]
 
meta_data[["Age"]] <- lfile$col.attrs$Age[cti_MBd]
names(meta_data[["Age"]]) <- lfile$col.attrs$CellID[cti_MBd]
  
meta_data[["Tissue"]] <- lfile$col.attrs$Tissue[cti_MBd]
names(meta_data[["Tissue"]]) <- lfile$col.attrs$CellID[cti_MBd]  
 
meta_data[["Location"]] <- lfile$col.attrs$Location_E9_E11[cti_MBd]
names(meta_data[["Location"]]) <- lfile$col.attrs$CellID[cti_MBd]
  
meta_data[["clusterID"]] <- lfile$col.attrs$ClusterName[cti_MBd]
names(meta_data[["clusterID"]]) <- lfile$col.attrs$CellID[cti_MBd]
  
meta_data[["subclass"]] <- lfile$col.attrs$Subclass[ cti_MBd ]
names(meta_data[["subclass"]]) <- lfile$col.attrs$CellID[cti_MBd]
  
  #add these labels to the Seurat object
  for (i in names(meta_data)){
    MBd_Seurat <- AddMetaData(
      object = MBd_Seurat,
      metadata = meta_data[[i]],
      col.name = i
    )  
  }
  
MBd_Seurat <- NormalizeData( MBd_Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
MBd_Seurat <- FindVariableFeatures( MBd_Seurat, selection.method = "vst", nfeatures = 2000)
  
# read pre-filtered data from Pten experiment
# this is a list
Pten_seurat_list <- readRDS ( paste(base, "RDS_files/P0_Pten_ctrl_RAW.RDS", sep="") )
  
# add e18 to the list
Pten_seurat_list[["e18"]] <- MBd_Seurat
# add genotype tag for LaManno data
Pten_seurat_list[["e18"]]@meta.data$genotype <- "e18"

# label transfer
message("integrating control")
embryo.anchors <- FindTransferAnchors(reference = MBd_Seurat, Pten_seurat_list[["ctrl"]], dims = 1:30)
saveRDS( embryo.anchors, paste(base, "/RDS_files/Pten_ctrl_anchors.RDS", sep="") )
 
message("integrating KO")
embryo.anchors <- FindTransferAnchors(reference = MBd_Seurat, query = Pten_seurat_list[["sparse"]], dims = 1:30)
saveRDS( embryo.anchors, paste(base, "/RDS_files/Pten_sparse_anchors.RDS", sep="") )

# save the meta data
saveRDS( MBd_Seurat@meta.data, paste(base, "/RDS_files/LaManno_e18_metadata.RDS", sep="") )

stop ("first step done - restart")


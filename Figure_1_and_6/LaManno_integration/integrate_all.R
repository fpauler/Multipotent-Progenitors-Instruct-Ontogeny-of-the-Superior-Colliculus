# script to integrate LaManno et al., Nature 2021 and Fzd embryo pool data

library (loomR)
library (Seurat)

# read location information - defines which cells to keep
# manually defined 
location2keep <- read.table("LaManno.Location-Giselle.csv", sep="\t", header = T)
location2keep <- location2keep[which(location2keep$use == "yes"),]

# process loom file
lfile <- connect(filename = "dev_all.loom", mode = "r+", skip.validate = T)

# get samples indices from Midbrain or Dorsal Midbrain
mbdIDX <- which( lfile$col.attrs$Tissue[] == "Midbrain" | lfile$col.attrs$Tissue[] == "MidbrainDorsal" )
mbdIDX <- mbdIDX[which(!(duplicated(lfile$col.attrs$CellID[mbdIDX]) | duplicated(lfile$col.attrs$CellID[mbdIDX], fromLast = T)))]

# additional filtering based on cluster annotation
sclassMBD <- which( grepl( x=lfile$col.attrs$Subclass[], pattern="midbrain", ignore.case = T) | grepl( x=lfile$col.attrs$Subclass[], pattern="Mixed" ) )
mbdIDX <- intersect( mbdIDX, sclassMBD)
  
# identify which cells from e9-11 should be removed
locIDX <-  which( !lfile$col.attrs$Location[] %in% c(location2keep$location) )
locIDX <- intersect( locIDX, which( lfile$col.attrs$Age[] %in% c("e9.0", "e10.0", "e11.0")))

# focus only on relevant cell types
cti_MBd <- intersect( mbdIDX, which( lfile$col.attrs$Class[] %in% c("Neuroblast", "Neuron", "Radial glia")))
  
# remove cells that are not from the location of interest
cti_MBd <- setdiff(cti_MBd, locIDX)

message("retained ", length(cti_MBd), " cells")
  
#extract the reads from these cells - NOTE row/col are changed according to convention!
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

# proces our own data - simple QC filtering
# Load the dataset: note that this file is the result of the cellRangre pipeline and available via GEO:

embryo.data <- Read10X_h5( filename = "Fzd10_MBd_filtered_feature_bc_matrix.h5")
    
# Initialize the Seurat object with the raw (non-normalized data).
embryo.seurat <- CreateSeuratObject(counts = embryo.data, min.cells = 3, min.features = 200)

embryo.seurat[["percent.mt"]] <- PercentageFeatureSet(embryo.seurat, pattern = "^mt-")

# remove low quality cells    
embryo.seurat <- subset(embryo.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

embryo.seurat <- NormalizeData(embryo.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
    
embryo.seurat <- FindVariableFeatures(embryo.seurat, selection.method = "vst", nfeatures = 2000)

embryo.seurat@meta.data$Age <- "embryonic_pool"
 
# split LaManno into separate Ages
Idents(MBd_Seurat) <- "Age"

# split into distinct ages and do normalisation
seurat.obj.list <- list()
  for (age in unique(MBd_Seurat@meta.data$Age)) {
    seurat.obj.list[[age]] <- subset(MBd_Seurat, idents = age, invert=F)
    seurat.obj.list[[age]] <- NormalizeData( seurat.obj.list[[age]], normalization.method = "LogNormalize", scale.factor = 10000)
    seurat.obj.list[[age]] <- FindVariableFeatures( seurat.obj.list[[age]], selection.method = "vst", nfeatures = 2000)    
  }

# add our own data
seurat.obj.list[["embryo"]] <- embryo.seurat

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seurat.obj.list)
  
anchors <- FindIntegrationAnchors(object.list = seurat.obj.list, anchor.features = features)
  
# this command creates an 'integrated' data assay
seurat.combined <- IntegrateData( anchorset = anchors )
  
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seurat.combined) <- "integrated"
  
# this file will be the base of downstream analyses
saveRDS(object = seurat.combined, file = "seurat.combined.RDS")

sessionInfo()
  

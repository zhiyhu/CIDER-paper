# preprocessing dendritic dataset
# Zhiyuan Hu
# 25 Sep 2020

library(data.table)
library(DoubletFinder)

dirsave <- "dendritic"
setwd("~/OneDrive - Nexus365/Project ClinCluster/evaluate-integration/evaluate-integration")

dir.create(file.path(paste0("plots/",dirsave)))
dir.create(file.path(paste0("plots/",dirsave), "preprocessing"))
dir.create(file.path(paste0("rdata"), dirsave))

## list of files
list_files <- list.files(paste0("raw_data/", dirsave))
list_files
# [1] "dataset1_sm_uc3.txt.gz" "sample_sm_uc3.txt.gz" 

## load data
mat <- read.table(paste0("raw_data/", dirsave,"/",list_files[1]), header = TRUE, row.names = 1)
df_cell <- read.table(paste0("raw_data/", dirsave,"/",list_files[2]), header = TRUE, row.names = 1, sep = "\t")
all(rownames(df_cell) == colnames(mat))
# [1] TRUE

## count total features
df_cell$total_features <- colSums(mat > 0)

## filtering by a minimum 500 genes detected
mat <- mat[,df_cell$total_features >= 500]
df_cell <- df_cell[df_cell$total_features >= 500,]
head(df_cell)

## create seurat obj
seu <- CreateSeuratObject(counts = mat, project = "dataset", min.cells = 0, min.features = 500)
seu$Batch <-  df_cell$batch              
seu$Group <-  df_cell$celltype

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "Group")
ggsave(paste0("plots/",dirsave,"/preprocessing/qc_vln.pdf"))


table(seu$Batch)
# Batch1 Batch2 
# 281    283 

## - - - - - - - - - - - - - - - #
## Seurat preprocessing pipeline #
## - - - - - - - - - - - - - - - #
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 20, verbose = FALSE)
ElbowPlot(seu)
seu <- RunTSNE(seu, reduction = "pca", dims = 1:10)

table(seu$Batch)
# Batch1 Batch2 
# 281    283 

## save rds
saveRDS(seu, file = paste0("rdata/",dirsave,"/seurat_object_preprocessed.rds"), compress = TRUE)



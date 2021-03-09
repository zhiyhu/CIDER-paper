## Preprocess human and mouse pancreas data
## Zhiyuan Hu
## 15 Sep 2020

library(data.table)
dirsave <- "pancreas_twospecies"
setwd("/Users/zhiyuan/OneDrive - Nexus365/Project ClinCluster/evaluate-integration/evaluate-integration")

## list of files
list_files <- list.files("raw_data/pancreas_twospecies/", pattern = "*.csv.gz")

## human data
df_human <- c()
counts_human <- c()
for(i in 1:4) {
  mat <- read.table(paste0("raw_data/pancreas_twospecies/",list_files[i]), sep = ",",
                    header = TRUE, row.names = 1)
  df_cell <- mat[,1:2]
  df_human <- rbind(df_human, df_cell)
  mat <- mat[,-1:-2]
  mat <- as.matrix(t(mat))
  counts_human <- cbind(counts_human, mat)
}

## mouse data
df_mouse <- c()
counts_mouse <- c()
for(i in 5:6) {
  mat <- read.table(paste0("raw_data/pancreas_twospecies/",list_files[i]), sep = ",",
                    header = TRUE, row.names = 1)
  df_cell <- mat[,1:2]
  df_mouse <- rbind(df_mouse, df_cell)
  mat <- mat[,-1:-2]
  mat <- as.matrix(t(mat))
  counts_mouse <- cbind(counts_mouse, mat)
}

rm(mat, df_cell)
rm(list_files)

## count total features
df_mouse$total_features <- colSums(counts_mouse > 0)
df_human$total_features <- colSums(counts_human > 0)

## filtering by a minimum 500 genes detected
counts_mouse <- counts_mouse[,df_mouse$total_features >= 500]
df_mouse <- df_mouse[df_mouse$total_features >= 500,]

counts_human <- counts_human[,df_human$total_features >= 500]
df_human <- df_human[df_human$total_features >= 500,]

###########
## Human ##
###########

## create seurat obj
seu <- CreateSeuratObject(counts = counts_human, project = "dataset", min.cells = 0, min.features = 500)
seu$Batch <-   substr(colnames(counts_human), 1, 6)                               
seu$Group <-   df_human$assigned_cluster    
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "Group")

## Prepare for doubletfinder
seu_list <- SplitObject(seu, split.by = "Batch")
for(i in 1:length(seu_list)){
  seu_list[[i]] <- NormalizeData(seu_list[[i]], verbose = FALSE)
  seu_list[[i]] <- FindVariableFeatures(seu_list[[i]], selection.method = "vst",nfeatures = 2000, verbose = FALSE)
  seu_list[[i]] <- ScaleData(seu_list[[i]], verbose = FALSE)
  seu_list[[i]] <- RunPCA(seu_list[[i]], npcs = 20, verbose = FALSE)
  seu_list[[i]] <- RunTSNE(seu_list[[i]], reduction = "pca", dims = 1:10)
}

## Double Detection
plist <- list()
for(i in 1:length(seu_list)){
  sweep.res.list <- paramSweep_v3(seu_list[[i]], PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  head(bcmvn)
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  annotations <- seu_list[[i]]@meta.data$Group
  homotypic.prop <- modelHomotypic(annotations)       
  nExp_poi <- round(0.05*ncol(seu_list[[i]]))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seu_list[[i]] <- doubletFinder_v3(seu_list[[i]], PCs = 1:10, pN = 0.25, 
                                    pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  seu_list[[i]] <- doubletFinder_v3(seu_list[[i]], PCs = 1:10, pN = 0.25, 
                                    pK = pK, nExp = nExp_poi.adj,
                                    reuse.pANN = colnames(seu_list[[i]]@meta.data)[grep("pANN_0.25",colnames(seu_list[[i]]@meta.data))], sct = FALSE)
  
  colname <- colnames(seu_list[[i]]@meta.data)[grep("DF.classifications",colnames(seu_list[[i]]@meta.data))][2]
  seu_list[[i]]$doublet <- seu_list[[i]]@meta.data[,colname]
  plist[[i]] <- DimPlot(seu_list[[i]], group.by = "doublet")
}

cowplot::plot_grid(plotlist = plist, ncol = 4)
ggsave("plots/pancreas_twospecies/preprocessing/human_doublets.png", width = 14, height = 3)

## Doublet filtering
seu$doublet <- NA
for(i in 1:length(seu_list)){
  seu$doublet[match(colnames(seu_list[[i]]), colnames(seu))] <- seu_list[[i]]$doublet
}
table(seu$doublet)

seu <- seu[,which(seu$doublet == "Singlet")]

## Seurat preprocessing pipeline
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 20, verbose = FALSE)
seu <- RunTSNE(seu, reduction = "pca", dims = 1:10)
table(seu$Batch)

## save rds
saveRDS(seu, file = paste0("rdata/",dirsave,"/seurat_object_preprocessed_human.rds"), compress = TRUE)
seuh <- seu

###########
## Mouse ##
###########

## create seurat obj
seu <- CreateSeuratObject(counts = counts_mouse, project = "dataset", min.cells = 0, min.features = 500)
seu$Batch <- substr(colnames(counts_mouse), 1, 6)                               
seu$Group <- df_mouse$assigned_cluster    
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "Group")

## Prepare for doubletfinder
# seu_list <- SplitObject(seu, split.by = "Batch")
# for(i in 1:length(seu_list)){
#   seu_list[[i]] <- NormalizeData(seu_list[[i]], verbose = FALSE)
#   seu_list[[i]] <- FindVariableFeatures(seu_list[[i]], selection.method = "vst",nfeatures = 2000, verbose = FALSE)
#   seu_list[[i]] <- ScaleData(seu_list[[i]], verbose = FALSE)
#   seu_list[[i]] <- RunPCA(seu_list[[i]], npcs = 20, verbose = FALSE)
#   seu_list[[i]] <- RunTSNE(seu_list[[i]], reduction = "pca", dims = 1:10)
# }

## Double Detection
# plist <- list()
# for(i in 1:length(seu_list)){
#   sweep.res.list <- paramSweep_v3(seu_list[[i]], PCs = 1:10, sct = FALSE)
#   sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
#   bcmvn <- find.pK(sweep.stats)
#   head(bcmvn)
#   pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
#   
#   annotations <- seu_list[[i]]@meta.data$Group
#   homotypic.prop <- modelHomotypic(annotations)       
#   nExp_poi <- round(0.05*ncol(seu_list[[i]]))
#   nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#   seu_list[[i]] <- doubletFinder_v3(seu_list[[i]], PCs = 1:10, pN = 0.25, 
#                                     pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#   
#   seu_list[[i]] <- doubletFinder_v3(seu_list[[i]], PCs = 1:10, pN = 0.25, 
#                                     pK = pK, nExp = nExp_poi.adj,
#                                     reuse.pANN = colnames(seu_list[[i]]@meta.data)[grep("pANN_0.25",colnames(seu_list[[i]]@meta.data))], sct = FALSE)
#   
#   colname <- colnames(seu_list[[i]]@meta.data)[grep("DF.classifications",colnames(seu_list[[i]]@meta.data))][2]
#   seu_list[[i]]$doublet <- seu_list[[i]]@meta.data[,colname]
#   plist[[i]] <- DimPlot(seu_list[[i]], group.by = "doublet")
# }
# 
# cowplot::plot_grid(plotlist = plist, ncol = 4)
# ggsave("plots/pancreas_twospecies/preprocessing/mouse_doublets.png", width = 14, height = 3)


## Doublet filtering
# seu$doublet <- NA
# for(i in 1:length(seu_list)){
#   seu$doublet[match(colnames(seu_list[[i]]), colnames(seu))] <- seu_list[[i]]$doublet
# }
# table(seu$doublet)
# seu <- seu[,which(seu$doublet == "Singlet")]

## Seurat preprocessing pipeline
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 20, verbose = FALSE)
seu <- RunTSNE(seu, reduction = "pca", dims = 1:10)
table(seu$Batch)


## save rds
saveRDS(seu, file = paste0("rdata/",dirsave,"/seurat_object_preprocessed_mouse.rds"), compress = TRUE)
seum <- seu
rm(seu)

###################
## merge two seu obj
hu_genes <- rownames(seuh)
ms_genes <- rownames(seum)
sum(toupper(ms_genes) %in% hu_genes)
sum(hu_genes %in% toupper(ms_genes))


## adjust the genenames
shared_genes <- toupper(ms_genes)[toupper(ms_genes) %in% hu_genes]
shared_genes <- unique(shared_genes)
shared_genes1 <- c(shared_genes, "INS1")
shared_genes2 <- c(shared_genes, "INS")

all(shared_genes1 %in% toupper(rownames(seum)))

mat_ms <- seum@assays$RNA@counts[match(shared_genes1, toupper(rownames(seum))),]
mat_hu <- seuh@assays$RNA@counts[match(shared_genes2, toupper(rownames(seuh))),]

rownames(mat_ms) <- toupper(rownames(mat_ms))
rownames(mat_hu)[which(rownames(mat_hu) != rownames(mat_ms))] <- "INS1"

## add species infomation
seuh$sample <- seuh$Batch
seum$sample <- seum$Batch
seuh$Batch <- "human"
seum$Batch <- "mouse"

info_hu <- seuh@meta.data
info_ms <- seum@meta.data
info_ms$doublet <- "Singlet"

seu <- CreateSeuratObject(counts = cbind(mat_ms, mat_hu), meta.data = rbind(info_ms, info_hu))
seu$Group <- tolower(seu$Group)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 20, verbose = FALSE)
seu <- RunTSNE(seu, reduction = "pca", dims = 1:15)

## save merged
saveRDS(seu, paste0("rdata/",dirsave,"/seurat_object_preprocessed.rds"))




# Seurat Louvain on COVID data
# 25 Jan 2021
# Zhiyuan Hu
# Last modified 3 Aug 2021

# Set up ----
library(Seurat)
library(mclust)

verbose <- FALSE
dirsave <- "covid19"
source("/home/z/zhu/cider/functions_pipelines.R")

# read data ----
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 

# Seurat -----
runtime <- system.time({
	seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, features = VariableFeatures(object = seu),verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:15)
  seu <- FindClusters(seu, resolution = 0.4)
})[3]
print(runtime)

seu <- RunTSNE(seu, dims = 1:15)

# Save results ----
results <- data.frame(method = "seurat",
                     runtime  = runtime,
                     ARI = mclust::adjustedRandIndex(seu$seurat_clusters, seu$Group),
                     batch_ARI = mclust::adjustedRandIndex(seu$seurat_clusters, seu$Batch))
clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu$seurat_clusters,
                             dim1 = seu@reductions$tsne@cell.embeddings[,1],
                             dim1 = seu@reductions$tsne@cell.embeddings[,2])

saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/seurat_ARI_res.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/seurat_clustering_res.rds"))



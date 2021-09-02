# Seurat Louvain on BC data
# 6 Aug 2021
# Zhiyuan Hu
# Last modified 6 Aug 2021

# Set up ----
library(Seurat)
library(mclust)

verbose <- FALSE
dirsave <- "pancan"
source("/home/z/zhu/cider/functions_pipelines.R")

# read data ----
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 

# Seurat -----
runtime <- system.time({
	seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 1000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, features = VariableFeatures(object = seu),verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:10)
  seu <- FindClusters(seu, resolution = 0.3)
})[3]
print(runtime)

seu <- RunTSNE(seu, dims = 1:20, check_duplicates = FALSE)

# Save results ----
results <- data.frame(method = "seurat",
                     runtime  = runtime,
                     ARI = mclust::adjustedRandIndex(seu$seurat_clusters, seu$Group),
                     batch_ARI = mclust::adjustedRandIndex(seu$seurat_clusters, seu$Batch))
clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu$seurat_clusters,
                             group = seu$Group,
                             batch = seu$Batch,
                             dim1 = seu@reductions$tsne@cell.embeddings[,1],
                             dim2 = seu@reductions$tsne@cell.embeddings[,2])
results

saveRDS(as.data.frame(seu@meta.data), paste0("/home/z/zhu/cider/rdata/", dirsave, "/seurat_object_metadata.rds"))
saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/seurat_ARI_res.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/seurat_clustering_res.rds"))



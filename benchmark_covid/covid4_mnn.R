# MNN on COVID data
# 25 Jan 2021
# Zhiyuan Hu
# Last modified 25 Jan 2021

# Set up ----
library(Seurat)
library(batchelor)
library(scater)
library(scran)
source("functions_simulate.R")
verbose <- FALSE
dirsave <- "covid19"

# read data ----
seu <- readRDS(paste0("rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
gc()

## Run MNN ----
runtime <- system.time({
  f.out <- mnn_integration(seu)
  seu@reductions$pca@cell.embeddings <- reducedDim(f.out, "corrected")
  seu <- FindNeighbors(seu, dims = 1:15)
  seu <- FindClusters(seu, resolution = 0.4)
})
f.out$cluster <- seu$cluster <- as.factor(seu$seurat_clusters)

# Save results ----
results <- data.frame(method = "MNN",
                      runtime  = as.numeric(runtime[3]),
                      ARI = mclust::adjustedRandIndex(seu$Group, seu$cluster),
                      batch_ARI = mclust::adjustedRandIndex(seu$Batch, seu$cluster))

clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu$cluster)

saveRDS(f.out, paste0("rdata/", dirsave, "/MNN_intergrated_SCEObject.rds"))
saveRDS(results, paste0("rdata/", dirsave, "/MNN_ARI_res.rds"))
saveRDS(clustering_res, paste0("rdata/", dirsave, "/MNN_clustering_res.rds"))


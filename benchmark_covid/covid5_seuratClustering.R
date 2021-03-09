# Seurat Louvain on COVID data
# 25 Jan 2021
# Zhiyuan Hu
# Last modified 25 Jan 2021

# Set up ----
library(Seurat)
library(mclust)

verbose <- FALSE
dirsave <- "covid19"
source("functions_pipelines.R")

# read data ----
seu <- readRDS(paste0("rdata/",dirsave,"/seurat_object_preprocessed.rds")) 

# Seurat -----
runtime <- system.time({
  seu <- seurat_clustering(seu_obj = seu, res = 0.4)
})[3]


# Save results ----
results <- data.frame(method = "seurat",
                     runtime  = runtime,
                     ARI = mclust::adjustedRandIndex(seu$seurat_cluster, seu$Group),
                     batch_ARI = mclust::adjustedRandIndex(seu$seurat_cluster, seu$Batch))
clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu$seurat_cluster)

saveRDS(results, paste0("rdata/", dirsave, "/seurat_ARI_res.rds"))
saveRDS(clustering_res, paste0("rdata/", dirsave, "/seurat_clustering_res.rds"))



# Scanorama on pbmc data
# 25 Jan 2021
# Zhiyuan Hu
# Last modified 2 Aug 2021

# Set up ----
library(Seurat)
library(reticulate)

verbose <- FALSE
dirsave <- "pbmc"

# read data ----
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
gc()

## Run scanorama ----
datasets <- list(t(as.matrix(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[1]])),
                 t(as.matrix(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[2]])))

# List of gene lists:
genes_list <- list(rownames(seu),rownames(seu))
scanorama <- import('scanorama')
# Integration
runtime <- system.time({
  integrated.data <- scanorama$integrate(datasets, genes_list)
  
  tmp <- t(integrated.data[[1]][[1]])
  tmp <- cbind(tmp,  t(integrated.data[[1]][[2]]))
  
  colnames(tmp) <- c(colnames(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[1]]), 
                     colnames(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[2]]))
  tmp <- t(tmp)
  tmp <- tmp[match(colnames(seu), rownames(tmp)),]
  seu@reductions$pca@cell.embeddings <- tmp
  seu <- FindNeighbors(seu, dims = 1:10)
  seu <- FindClusters(seu, resolution = 0.4)
})[3]

print(runtime)

# Save results ----
results <- data.frame(method = "scanorama",
                      runtime  = as.numeric(runtime),
                      ARI = mclust::adjustedRandIndex(seu$Group, seu$seurat_clusters),
                      batch_ARI = mclust::adjustedRandIndex(seu$Batch, seu$seurat_clusters))

clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu$seurat_clusters)
clustering_res <- cbind(clustering_res, seu@reductions$pca@cell.embeddings)

# saveRDS(tmp, paste0("rdata/", dirsave, "/scanorama_corrected_PCA_downsampled.rds"))
saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/scanorama_ARI_res.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/scanorama_clustering_res.rds"))


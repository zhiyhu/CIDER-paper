# Scanorama on COVID data
# 25 Jan 2021
# Zhiyuan Hu
# Last modified 19 July 2021

# Set up ----
library(Seurat)
library(reticulate)

verbose <- FALSE
dirsave <- "covid19"

# read data ----
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
gc()

## Run scanorama ----
# set.seed(12345)
# select <- sample(1:ncol(seu), ncol(seu)/10, replace = FALSE)
# n_groups <- length(unique(seu[,select]$Group))
# seu <- seu[, select] # ncol 1191

datasets <- list(t(as.matrix(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[1]])),
                 t(as.matrix(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[2]])),
                 t(as.matrix(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[3]])),
                 t(as.matrix(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[4]])),
                 t(as.matrix(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[5]]))
                 )

# List of gene lists:
genes_list <- list(rownames(seu),rownames(seu),rownames(seu),rownames(seu),rownames(seu))
scanorama <- import('scanorama')
# Integration
runtime <- system.time({
  integrated.data <- scanorama$integrate(datasets, genes_list)
  
  tmp <- t(integrated.data[[1]][[1]])
  tmp <- cbind(tmp,  t(integrated.data[[1]][[2]]))
  tmp <- cbind(tmp,  t(integrated.data[[1]][[3]]))
  tmp <- cbind(tmp,  t(integrated.data[[1]][[4]]))
  tmp <- cbind(tmp,  t(integrated.data[[1]][[5]]))
  
  colnames(tmp) <- c(colnames(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[1]]), 
                     colnames(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[2]]),
                     colnames(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[3]]),
                     colnames(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[4]]),
                     colnames(seu@assays$RNA@counts[,seu$Batch == unique(seu$Batch)[5]]))
  tmp <- t(tmp)
  tmp <- tmp[match(colnames(seu), rownames(tmp)),]
  seu@reductions$pca@cell.embeddings <- tmp
  seu <- FindNeighbors(seu, dims = 1:15)
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


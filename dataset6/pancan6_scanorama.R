# Scanorama on cancer data
# 8 Aug 2021
# Zhiyuan Hu
# Last modified 8 Aug 2021

# Set up ----
library(Seurat)
library(reticulate)

verbose <- FALSE
dirsave <- "pancan"
scanorama <- import('scanorama')

# read data ----
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
gc()

## Run scanorama ----
# set.seed(12345)
# select <- sample(1:ncol(seu), ncol(seu)/10, replace = FALSE)
# n_groups <- length(unique(seu[,select]$Group))
# seu <- seu[, select] # ncol 1191

# create PCA 
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", 
                            nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 20, verbose = FALSE)

# prepare data for Scanorama
seu.list <- SplitObject(seu, split.by = "Batch")
datasets <- list()
genes_list <- list()
for(i in 1:length(seu.list)){
  datasets[[i]] <- t(as.matrix(seu.list[[i]]@assays$RNA@counts))
  genes_list[[i]] <- rownames(seu)# List of gene lists
}

# Integration
runtime <- system.time({
  integrated.data <- scanorama$integrate(datasets, genes_list)

  tmp <- c()
  tmp_colnames <- c()
  for(i in 1:length(seu.list)){
    tmp <- cbind(tmp, t(integrated.data[[1]][[i]]))
    tmp_colnames <- c(tmp_colnames, colnames(seu.list[[i]]@assays$RNA@counts))
  }

  colnames(tmp) <- tmp_colnames
  tmp <- t(tmp)
  tmp <- tmp[match(colnames(seu), rownames(tmp)),]
  seu@reductions$pca@cell.embeddings <- tmp
  seu <- FindNeighbors(seu, dims = 1:20)
  seu <- FindClusters(seu, resolution = 0.3)
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


# Combat on pancreas_twospecies data
# 2 Aug 2021
# Zhiyuan Hu
# Last modified 2 Aug 2021

# Set up ----
library(Seurat)
library(sva)

verbose <- FALSE
dirsave <- "pancreas_twospecies"

# read data ----
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
gc()

# Integration
runtime <- system.time({
  mat <- log2(as.matrix(seu@assays$RNA@counts) + 1)
  mode(mat)
  dim(mat)
  range(rowSums(mat))
  
  mat <- mat[rowSums(mat) >= 5,] # have to filter out no-expressing genes
  
  corrected <- ComBat(
    dat = as.matrix(mat), 
    batch = seu$Batch, 
    mod = NULL,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  seu.integrated <- CreateSeuratObject(counts = corrected)
  seu.integrated@assays$RNA@scale.data <- corrected
  seu.integrated <- FindVariableFeatures(seu.integrated, selection.method = "vst",
                                         nfeatures = 2000, verbose = FALSE)
  seu.integrated <- RunPCA(object = seu.integrated, pcs.compute = 20)
  
  seu.integrated <- FindNeighbors(seu.integrated, dims = 1:15)
  seu.integrated <- FindClusters(seu.integrated, resolution = 0.4)
})[3]

print(runtime)

# Save results ----
results <- data.frame(method = "combat",
                      runtime  = as.numeric(runtime),
                      ARI = mclust::adjustedRandIndex(seu$Group, seu.integrated$seurat_clusters),
                      batch_ARI = mclust::adjustedRandIndex(seu$Batch, seu.integrated$seurat_clusters))

clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu.integrated$seurat_clusters)
clustering_res <- cbind(clustering_res, seu.integrated@reductions$pca@cell.embeddings)

# saveRDS(tmp, paste0("rdata/", dirsave, "/scanorama_corrected_PCA_downsampled.rds"))
saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/combat_ARI_res.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/combat_clustering_res.rds"))


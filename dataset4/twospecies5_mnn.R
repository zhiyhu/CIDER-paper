# MNN on cross species data
# 25 Jan 2021
# Zhiyuan Hu
# Last modified 3 Aug 2021

# Set up ----
library(Seurat)
library(batchelor)
library(scater)
library(scran)

verbose <- FALSE
dirsave <- "pancreas_twospecies"

# read data ----
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
gc()

## Run MNN ----
runtime <- system.time({
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst",
                              nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  sce <- SingleCellExperiment(assays = list(counts = seu@assays$RNA@counts, 
                                            logcounts = seu@assays$RNA@data))
  # logcounts(sce) <- log2(calculateCPM(sce) + 1)
  sce$Batch  <- seu$Batch
  dec <- modelGeneVar(sce)
  chosen.hvgs <- dec$bio > 0
  
  set.seed(1000)
  f.out <- fastMNN(sce, batch = sce$Batch, subset.row=chosen.hvgs)
  # str(reducedDim(f.out, "corrected"))
  f.out$Group <- seu$Group
  f.out$Batch <- seu$Batch
  seu@reductions$pca@cell.embeddings <- reducedDim(f.out, "corrected")
  seu <- FindNeighbors(seu, dims = 1:15)
  seu <- FindClusters(seu, resolution = 0.4)
})
f.out$cluster <- seu$cluster <- as.factor(seu$seurat_clusters)

print("runtime: ")
print(runtime[3])

# Save results ----
results <- data.frame(method = "MNN",
                      runtime  = as.numeric(runtime[3]),
                      ARI = mclust::adjustedRandIndex(seu$Group, seu$cluster),
                      batch_ARI = mclust::adjustedRandIndex(seu$Batch, seu$cluster))

clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu$cluster)
clustering_res <- cbind(clustering_res, seu@reductions$pca@cell.embeddings)

# saveRDS(f.out, paste0("/home/z/zhu/cider/rdata/", dirsave, "/MNN_intergrated_SCEObject.rds"))
saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/MNN_ARI_res.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/MNN_clustering_res.rds"))


# run rpca on cancer data
# 13 Aug 2021
# Zhiyuan Hu
# last modified 13 Aug 2021

# Set up ----
library(Seurat)

verbose <- FALSE
dirsave <- "pancan"

# read data ----
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
seu <- Seurat::CreateSeuratObject(counts = seu@assays$RNA@counts, 
                                  meta.data = data.frame(Batch = seu$Batch))
gc()

## Run cca ----
runtime <- system.time({
  seu.list <- SplitObject(seu, split.by = "Batch")
  
  for (i in 1:length(seu.list)) {
    seu.list[[i]] <- NormalizeData(seu.list[[i]], verbose = FALSE)
    seu.list[[i]] <- FindVariableFeatures(seu.list[[i]], selection.method = "vst", 
                                          nfeatures = 2000, verbose = FALSE)
    
  }
  # normalize and identify variable features for each dataset independently
  seu.list <- lapply(X = seu.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  
  features <- SelectIntegrationFeatures(object.list = seu.list)
  seu.list <- lapply(X = seu.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  seu.anchors <- FindIntegrationAnchors(object.list = seu.list, anchor.features = features, reduction = "rpca")
  rm(seu.list, seu)
  gc()
  seu.integrated <- IntegrateData(anchorset = seu.anchors, k.weight=50)
  
  DefaultAssay(seu.integrated) <- "integrated"
  # Run the standard workflow for visualization and clustering
  seu.integrated <- ScaleData(seu.integrated, verbose = FALSE)
  seu.integrated <- RunPCA(seu.integrated, npcs = 20, verbose = FALSE)
  seu.integrated <- FindNeighbors(seu.integrated, dims = 1:20)
  seu.integrated <- FindClusters(seu.integrated, resolution = 0.3)
  
})
print("runtime:")
print(runtime)

seu.integrated <- RunTSNE(seu.integrated, reduction = "pca", dims = 1:20)
seu.integrated$cluster <- as.factor(seu.integrated$seurat_clusters)

rm(seu.anchors)
gc()

seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
seu.integrated$Group <- seu$Group

# Save results ----
results <- data.frame(method = "rpca",
                      runtime  = as.numeric(runtime[3]),
                      ARI = mclust::adjustedRandIndex(seu.integrated$Group, seu.integrated$cluster),
                      batch_ARI = mclust::adjustedRandIndex(seu.integrated$Batch, seu.integrated$cluster))

clustering_res <- data.frame(sample = colnames(seu.integrated),
                             clusters = seu.integrated$cluster,
                             tSNE1 = seu.integrated@reductions$tsne@cell.embeddings[,1],
                             tSNE2 = seu.integrated@reductions$tsne@cell.embeddings[,2])

saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/rpca_ARI_res.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/rpca_clustering_res.rds"))

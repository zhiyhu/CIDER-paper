# Harmony on simulated data
# 10 Aug 2021
# Zhiyuan Hu
# Last modified 10 Aug 2021

# Set up ----
library(Seurat)
library(harmony)

verbose <- FALSE
dirsave <- "simulation"

for(itor in 1:20){
  # read data----
  seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed",itor,".rds"))
  gc()
  
  ## Run Harmony integration ----
  runtime <- system.time({
    seu <- NormalizeData(seu, verbose = FALSE)
    seu <- FindVariableFeatures(seu, selection.method = "vst",
                                nfeatures = 2000, verbose = FALSE)
    seu <- ScaleData(seu, verbose = FALSE)
    seu <- RunPCA(seu, npcs = 20, verbose = FALSE)
    seu.integrated <- RunHarmony(seu, group.by.vars = "Batch", assay.use = "RNA")
    seu.integrated <- FindNeighbors(seu.integrated, dims = 1:10, reduction = "harmony")
    seu.integrated <- FindClusters(seu.integrated, resolution = 0.3)
  })
  seu.integrated <- RunTSNE(seu.integrated, reduction = "harmony", dims = 1:10)
  seu.integrated$cluster <- as.factor(seu.integrated$seurat_clusters)
  gc()
  
  print("Harmony runtime:")
  print(runtime)
  
  # Save results ----
  results <- data.frame(method = "harmony",
                        runtime  = runtime[3],
                        ARI = mclust::adjustedRandIndex(seu.integrated$Group, seu.integrated$cluster),
                        batch_ARI = mclust::adjustedRandIndex(seu.integrated$Batch, seu.integrated$cluster))
  clustering_res <- data.frame(sample = colnames(seu.integrated),
                               clusters = seu.integrated$cluster,
                               dim1 = seu.integrated@reductions$tsne@cell.embeddings[,1],
                               dim2 = seu.integrated@reductions$tsne@cell.embeddings[,2])
  
  saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/harmony_ARI_res",itor,".rds"))
  saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/harmony_clustering_res",itor,".rds"))
}
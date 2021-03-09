# Seurat CCA on COVID data
# 24 Jan 2021
# Zhiyuan Hu
# Last modified 25 Jan 2021
options(future.globals.maxSize = 4000 * 1024^5) # https://www.biostars.org/p/432389/
# Set up ----
library(Seurat)

verbose <- FALSE
dirsave <- "covid19"

# source("functions_simulate.R")

# read data ----
seu <- readRDS(paste0("rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
seu <- Seurat::CreateSeuratObject(counts = seu@assays$RNA@counts, 
                                  meta.data = data.frame(Batch = seu$Batch))
gc()
## Run cca ----

start_time <- Sys.time()
seu.list <- SplitObject(seu, split.by = "Batch")
for (i in 1:length(seu.list)) {
  seu.list[[i]] <- NormalizeData(seu.list[[i]], verbose = FALSE)
  seu.list[[i]] <- FindVariableFeatures(seu.list[[i]], selection.method = "vst", 
                                        nfeatures = 1000, verbose = FALSE)
}
seu.anchors <- FindIntegrationAnchors(object.list = seu.list, dims = 1:15)
rm(seu.list, seu)
gc()
seu.integrated <- IntegrateData(anchorset = seu.anchors, dims = 1:15)

DefaultAssay(seu.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
seu.integrated <- ScaleData(seu.integrated, verbose = FALSE)
seu.integrated <- RunPCA(seu.integrated, npcs = 15, verbose = FALSE)
seu.integrated <- RunTSNE(seu.integrated, reduction = "pca", dims = 1:15)
seu.integrated <- FindNeighbors(seu.integrated, dims = 1:15)
seu.integrated <- FindClusters(seu.integrated, resolution = 0.4)
end_time <- Sys.time()
runtime <- end_time - start_time
seu.integrated$cluster <- as.factor(seu.integrated$seurat_clusters)

rm(seu.anchors)
gc()

seu <- readRDS(paste0("rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
seu.integrated$Group <- seu$Group

# Save results ----
results <- data.frame(method = "CCA",
                      runtime  = as.numeric(runtime) * 60 * 60,
                      ARI = mclust::adjustedRandIndex(seu.integrated$Group, seu.integrated$cluster),
                      batch_ARI = mclust::adjustedRandIndex(seu.integrated$Batch, seu.integrated$cluster))

clustering_res <- data.frame(sample = colnames(seu.integrated),
                             clusters = seu.integrated$cluster)

saveRDS(seu.integrated, paste0("rdata/", dirsave, "/CCA_intergrated_SeuratObject.rds"))
saveRDS(results, paste0("rdata/", dirsave, "/CCA_ARI_res.rds"))
saveRDS(clustering_res, paste0("rdata/", dirsave, "/CCA_clustering_res.rds"))

# Harmony on COVID data
# 6 Feb 2021
# Zhiyuan Hu
# Last modified 6 Feb 2021
# Set up ----
library(Seurat)
library(harmony)

verbose <- FALSE
dirsave <- "covid19"

# source("functions_simulate.R")

# read data ----
seu <- readRDS(paste0("rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
gc()

## Run Harmony integration ----
start_time <- Sys.time()
seu.integrated <- RunHarmony(seu, group.by.vars = "Batch", assay.use = "RNA")
seu.integrated <- FindNeighbors(seu.integrated, dims = 1:15, reduction = "harmony")
seu.integrated <- FindClusters(seu.integrated, resolution = 0.4)
end_time <- Sys.time()
runtime <- end_time - start_time
# seu.integrated <- RunTSNE(seu.integrated, reduction = "harmony", dims = 1:15)
seu.integrated$cluster <- as.factor(seu.integrated$seurat_clusters)
gc()


# Save results ----
results <- data.frame(method = "harmony",
                      runtime  = as.numeric(runtime) * 60,
                      ARI = mclust::adjustedRandIndex(seu.integrated$Group, seu.integrated$cluster),
                      batch_ARI = mclust::adjustedRandIndex(seu.integrated$Batch, seu.integrated$cluster))

clustering_res <- data.frame(sample = colnames(seu.integrated),
                             clusters = seu.integrated$cluster)

saveRDS(seu.integrated, paste0("rdata/", dirsave, "/harmony_intergrated_SeuratObject.rds"))
saveRDS(results, paste0("rdata/", dirsave, "/harmony_ARI_res.rds"))
saveRDS(clustering_res, paste0("rdata/", dirsave, "/harmony_clustering_res.rds"))

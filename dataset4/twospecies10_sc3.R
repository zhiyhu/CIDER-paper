# SC3 on cross species data 
# 25 Jan 2021
# Zhiyuan Hu
# Last modified 19 July 2021

# Set up ----
library(SingleCellExperiment)
library(SC3)
library(scater)

verbose <- FALSE
dirsave <- "pancreas_twospecies"
source("/home/z/zhu/cider/functions_pipelines.R")

# read data ----
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
seu <- Seurat::CreateSeuratObject(counts = seu@assays$RNA@counts, 
                                  meta.data = data.frame(Batch = seu$Batch,
                                                         Group = seu$Group))
gc()

# SC3 ------
n_groups <- length(unique(seu$Group))
runtime <- system.time({
  seu <- sc3_clustering(seu_obj = seu, g = n_groups, n.cores = 8) 
})[3]

# Save results ----
results <- data.frame(method = "SC3",
                      runtime  = runtime,
                      ARI = mclust::adjustedRandIndex(seu$sc3_cluster, seu$Group),
                      batch_ARI = mclust::adjustedRandIndex(seu$sc3_cluster, seu$Batch))
clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu$sc3_cluster)

saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/sc3_ARI_res.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/sc3_clustering_res.rds"))


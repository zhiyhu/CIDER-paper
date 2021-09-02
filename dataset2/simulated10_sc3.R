# SC3 on simulated data
# 10 Aug 2021
# Zhiyuan Hu
# Last modified 10 Aug 2021

# Set up ----
library(SingleCellExperiment)
library(SC3)
library(scater)

verbose <- FALSE
dirsave <- "simulation"
source("/home/z/zhu/cider/functions_pipelines.R")

for(itor in 1:20){
  # read data ----
  seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed",itor,".rds")) 
  seu <- Seurat::CreateSeuratObject(counts = seu@assays$RNA@counts, 
                                        meta.data = data.frame(Batch = seu$Batch,
                                                               Group = seu$Group))
  gc()
  
  # SC3 ------
  n_groups <- length(unique(seu$Group))
  runtime <- system.time({
    seu <- sc3_clustering(seu_obj = seu, g = n_groups) 
  })[3]
  
  # Save results ----
  results <- data.frame(method = "SC3",
                        runtime  = runtime,
                        ARI = mclust::adjustedRandIndex(seu$sc3_cluster, seu$Group),
                        batch_ARI = mclust::adjustedRandIndex(seu$sc3_cluster, seu$Batch))
  clustering_res <- data.frame(sample = colnames(seu),
                               clusters = seu$sc3_cluster)
  
  saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/sc3_ARI_res",itor,".rds"))
  saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/sc3_clustering_res",itor,".rds"))
}

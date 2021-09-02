# Monocle on pbmc data
# 2 Aug 2021
# Zhiyuan Hu
# Last modified 2 Aug 2021

# Set up ----
library(Seurat)
library(monocle3)
library(dplyr)

verbose <- FALSE
dirsave <- "pbmc"

# read data ----
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
seu <- Seurat::CreateSeuratObject(counts = seu@assays$RNA@counts, 
                                  meta.data = data.frame(Batch = seu$Batch,
                                                         Group = seu$Group))
gc()

n_groups <- length(unique(seu$Group))

# clstering-----
runtime <- system.time({
  cds <- new_cell_data_set(seu@assays$RNA@counts,
                           cell_metadata = seu@meta.data
                           )
  cds <- preprocess_cds(cds, num_dim = 100)
  cds <- reduce_dimension(cds)
  cds = align_cds(cds, num_dim = 100, alignment_group = "Batch")
  cds = reduce_dimension(cds)
  cds = cluster_cells(cds, resolution=1e-5)
  
})[3]

seu$cluster <- clusters(cds)

# Save results ----
results <- data.frame(method = "raceid",
                      runtime  = runtime,
                      ARI = mclust::adjustedRandIndex(seu$cluster, seu$Group),
                      batch_ARI = mclust::adjustedRandIndex(seu$cluster, seu$Batch))
clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu$cluster)

saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/monocle_ARI_res.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/monocle_clustering_res.rds"))


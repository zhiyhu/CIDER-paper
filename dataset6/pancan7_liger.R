# LIGER on cancer data
# 8 Aug 2021
# Zhiyuan Hu
# Last modified 8 Aug 2021

# Set up ----
library(Seurat)
library(rliger)

verbose <- FALSE
dirsave <- "pancan"

# read data ----
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
gc()


# Integration
runtime <- system.time({
  lgr <- seuratToLiger(objects = seu, combined.seurat = TRUE, meta.var = "Batch")
  lgr <- normalize(lgr)
  lgr <- selectGenes(lgr)
  lgr <- scaleNotCenter(lgr)
  lgr <- optimizeALS(lgr, k = 20)
  lgr <- quantile_norm(lgr)
  lgr <- louvainCluster(lgr, resolution = 0.3)
})[3]

print(runtime)

# seu.integrated <- ligerToSeurat(lgr)



# Save results ----
results <- data.frame(method = "liger",
                      runtime  = as.numeric(runtime),
                      ARI = mclust::adjustedRandIndex(seu$Group, lgr@clusters),
                      batch_ARI = mclust::adjustedRandIndex(seu$Batch, lgr@clusters))

clustering_res <- data.frame(sample = colnames(seu),
                             clusters = lgr@clusters)
# clustering_res <- cbind(clustering_res, seu.integrated@reductions$inmf@cell.embeddings)

saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/liger_ARI_res.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/liger_clustering_res.rds"))


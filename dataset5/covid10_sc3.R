# SC3 on COVID data
# 25 Jan 2021
# Zhiyuan Hu
# Last modified 19 July 2021

# Set up ----
library(SingleCellExperiment)
library(SC3)
library(scater)

verbose <- FALSE
dirsave <- "covid19"
source("/home/z/zhu/cider/functions_pipelines.R")

# read data ----
seu <- readRDS(paste0("/home/z/zhu/cider/rdata/",dirsave,"/seurat_object_preprocessed.rds")) 
seu <- Seurat::CreateSeuratObject(counts = seu@assays$RNA@counts, 
                                  meta.data = data.frame(Batch = seu$Batch,
                                                         Group = seu$Group))
gc()

idx <- sample(1:ncol(seu), replace = FALSE, size = ncol(seu)/5)
seu <- seu[,idx]
# SC3 ------
n_groups <- length(unique(seu$Group)) # number of groups
runtime <- system.time({
  g <- n_groups
  sim.groups <- SingleCellExperiment(assays = list(counts = as.matrix(seu@assays$RNA@counts)))
  logcounts(sim.groups) <- log2(calculateCPM(sim.groups) + 1)
  rowData(sim.groups)$feature_symbol <- rownames(sim.groups)
  sim.groups <- SC3::sc3(sim.groups, ks = g, biology = F)
  seu$sc3_cluster <- colData(sim.groups)[,ncol(colData(sim.groups))]
  # seu <- sc3_clustering(seu_obj = seu, g = n_groups) 
})[3]

head(colData(sim.groups))
head(colData(sim.groups)[,ncol(colData(sim.groups))])
head(seu$sc3_cluster)

# Save results ----
results <- data.frame(method = "SC3",
                      runtime  = runtime,
                      ARI = mclust::adjustedRandIndex(seu$sc3_cluster, seu$Group),
                      batch_ARI = mclust::adjustedRandIndex(seu$sc3_cluster, seu$Batch))
clustering_res <- data.frame(sample = colnames(seu),
                             clusters = seu$sc3_cluster)

saveRDS(results, paste0("/home/z/zhu/cider/rdata/", dirsave, "/sc3_ARI_res.rds"))
saveRDS(clustering_res, paste0("/home/z/zhu/cider/rdata/", dirsave, "/sc3_clustering_res.rds"))


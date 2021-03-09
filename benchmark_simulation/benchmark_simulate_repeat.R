library(splatter)
library(Seurat)
# BiocManager::install("batchelor")
library(batchelor)
library(scater)
library(scran)
library(SingleCellExperiment)
library(wesanderson)

library(limma)
library(harmony)
library(reticulate)
scanorama <- import('scanorama')

# devtools::install_github("sctyner/geomnet")
library(igraph)
library(geomnet)

library(mclust)

source("functions_simulate.R")


## set seed -----------
set.seed(1234)
seeds <- sample(1:10000, size = 20, replace = FALSE)
# [1] 7452 8016 7162 8086 7269 9196  623  934 2948 2146 2774 5722 2374 1103 4366 3454 2622 3972 1685 5332

## repeat 20 times -----
df_ari <- data.frame(rep = 1:20,
                     our = NA,
                     cca = NA,
                     mnn = NA,
                     scan = NA,
                     harmony = NA)

for(itor in 1:20){
  params <- newSplatParams()
  params <- setParams(params, seed = seeds[itor], 
                      group.prob = rep(1/5, 5), 
                      batchCells = c(2000,2000,2000)
  )
  
  sim.groups <- splatSimulate(params, method = "groups", verbose = FALSE)
  
  sim.groups <- sim.groups[,paste0(sim.groups$Batch, sim.groups$Group) %in% c("Batch1Group1","Batch1Group2", "Batch1Group3","Batch2Group2","Batch2Group3","Batch2Group4","Batch3Group3","Batch3Group4","Batch3Group5")]
  logcounts(sim.groups) <- log2(calculateCPM(sim.groups) + 1)
  seu <- seurat_preprocess(sim.groups)
  
  # our pipeline
  metadata <- data.frame(label = seu$Group,
                         batch = seu$Batch,
                         V1 = seu@reductions$tsne@cell.embeddings[,1],
                         V2 = seu@reductions$tsne@cell.embeddings[,2])
  
  n_size <- 50
  select <- downsampling(metadata = metadata, n_size = n_size)
  
  matrix <- as.matrix(seu@assays$RNA@counts[,select])
  list_dist <- calculate_ider(matrix = matrix, metadata = metadata)
  dist_p <- list_dist[[3]]
  
  hc <- hclust(as.dist(1-(dist_p + t(dist_p))))
  res <- cutree(hc, k = 5)
  df_res <- data.frame(res)
  df_res$initial_cluster <- rownames(df_res)
  metadata$final_cluster <-  df_res$res[match(paste(metadata$label, metadata$batch, sep = "_"), df_res$initial_cluster)]
  df_ari$our[itor] <- adjustedRandIndex(metadata$label, metadata$final_cluster)
  
  ## CCA 
  seu.integrated <- cca_integration(seu)
  res <- dbscan::hdbscan(seu.integrated@reductions$tsne@cell.embeddings[,1:2], minPts = 25)
  seu.integrated$cluster <- as.factor(res$cluster)
  df_ari$cca[itor] <- adjustedRandIndex(seu.integrated$Group, seu.integrated$cluster)
  
  ## MNN
  f.out <- mnn_integration(seu)
  res <- dbscan::hdbscan(reducedDim(f.out,"TSNE"), minPts = 75)
  f.out$cluster <- as.factor(res$cluster)
  df_ari$mnn[itor] <- adjustedRandIndex(f.out$Group, f.out$cluster)
  
  ## scanorama
  datasets <- list(t(as.matrix(seu@assays$RNA@counts[,seu$Batch == "Batch1"])),
                   t(as.matrix(seu@assays$RNA@counts[,seu$Batch == "Batch2"])),
                   t(as.matrix(seu@assays$RNA@counts[,seu$Batch == "Batch3"])))
  genes_list <- list(rownames(seu),rownames(seu),rownames(seu))
  integrated.data <- scanorama$integrate(datasets, genes_list)
  
  tmp <- t(integrated.data[[1]][[1]])
  tmp <- cbind(tmp,  t(integrated.data[[1]][[2]]))
  tmp <- cbind(tmp,  t(integrated.data[[1]][[3]]))
  df <- RunTSNE(t(tmp)[,1:50])
  res <- dbscan::hdbscan(df@cell.embeddings, minPts = 75)
  df_ari$scan[itor] <- adjustedRandIndex(seu.integrated$Group, res$cluster)
  
  ## Harmony
  seu <- RunHarmony(seu, "Batch")
  seu <- FindNeighbors(seu, dims = 1:15, reduction = "harmony")
  seu <- FindClusters(seu, resolution = 0.4)
  seu$cluster <- as.factor(seu$seurat_clusters)
  df_ari$harmony[itor] <- adjustedRandIndex(seu$Group, seu$cluster)
}

# saveRDS(df_ari, "plots/simulation_repeat/df_ari20200807.rds")
saveRDS(df_ari, "plots/simulation_repeat/df_ari20210206.rds")